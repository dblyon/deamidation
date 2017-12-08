from __future__ import print_function
import re, argparse, sys, os

import numpy as np
import pandas as pd



class MQcolnames(object):

    def __init__(self, df):
        self.columns = df.columns.tolist()
        self.new2old_colnames_dict = {"Potential contaminant": "Contaminant",
                                      "Modified sequence": "Modified Sequence",
                                      "Raw file": "Raw File",
                                      "Razor + unique peptides": "Razor + unique Peptides",
                                      "Leading razor protein": "Leading Razor Protein"}

    def check_if_colname_or_similar_exists(self, colname):
        """
        return the correct column name or False
        :param colname: String
        :return: String or Boolean
        """
        if colname in self.columns:
            return colname
        else:
            if colname in self.new2old_colnames_dict:
                return self.new2old_colnames_dict[colname]
            else:
                print("Column doesn't exist: {}".format(colname))
                import ipdb
                ipdb.set_trace()
                return False

def add_num_N_Q_cols(df):
    """
    add 2 columns with Integers to DataFrame, count the number of N and Q from Sequence column
    :param df: DataFrame
    :return: DataFrame
    """
    aaseq2count_N_Q_dict = {}
    aaseq_arr = df["Sequence"].unique()
    for aaseq in aaseq_arr:
        aaseq2count_N_Q_dict[aaseq] = count_N_Q(aaseq)
    df['num_N'] = df["Sequence"].apply(lambda aaseq: aaseq2count_N_Q_dict[aaseq][0], 1)
    df['num_Q'] = df["Sequence"].apply(lambda aaseq: aaseq2count_N_Q_dict[aaseq][1], 1)
    return df

def count_N_Q(aaseq):
    N_count = aaseq.count('N')
    Q_count = aaseq.count('Q')
    return N_count, Q_count

def add_num_N2D_Q2E_ratios(df, mqc=None):
    """
    add 4 columns to DataFrame: Int, Int, Float, Float
    number of N modified to D and analogous Q to E,
    the ratio of N / (N + N2D) and
    analogous the ratio of Q / (Q + Q2E)
    :param df: DataFrame
    :param mqc: MQcolumns instance
    :return: DataFrame
    """
    if not mqc:
        mqc = MQcolnames(df)
    colname = mqc.check_if_colname_or_similar_exists("Modified sequence")
    aaseq2count_N2D_Q2E_dict = {}
    aaseq_mod_arr = df[colname].unique()
    for aaseq in aaseq_mod_arr:
        aaseq2count_N2D_Q2E_dict[aaseq] = count_N2D_Q2E(aaseq)
    df['num_N2D'] = df[colname].apply(lambda aaseq: aaseq2count_N2D_Q2E_dict[aaseq][0], 1)
    df['num_Q2E'] = df[colname].apply(lambda aaseq: aaseq2count_N2D_Q2E_dict[aaseq][1], 1)
    # sanity check: "Deamidation (NQ)" == ( "num_N2D" + "num_Q2E" )
    assert (df["Deamidation (NQ)"] == ( df["num_N2D"] + df["num_Q2E"] )).all()

    df['ratio_N2D'] = df['num_N2D']*1.0 / df['num_N']
    df['ratio_Q2E'] = df['num_Q2E']*1.0 / df['num_Q']
    return df

def count_N2D_Q2E(aaseq, char2find=r"\(de\)"):
    """
    count the number of N modified to D and Q modified to E
    e.g.
    aaseq = r"_AAIAAFNAQ(de)N(de)N(de)GSNFQIEEISR_"
    count_N2D_Q2E(aaseq)
    :param aaseq: String
    :param char2find: String
    :return: Tuple(Int, Int)
    """
    N_Q_list = [aaseq[i.start() - 1] for i in re.finditer(char2find, aaseq)]
    N2D_count = N_Q_list.count("N")
    Q2E_count = N_Q_list.count("Q")
    return N2D_count, Q2E_count

def percent_deamidation(group, abundance_colname="Intensity"):
    """
    calculate the deamidation rate (percentage) for N to D and Q to E;
    D = number of Aspartic acids (mod. N)
    N = number of Asparagines
    abu = abundance
    each per charge state
    (D0*abu_0 + D1*abu_1 + D2*abu_2) / ((D0+N0)*abu_0 + (D1+N1)*abu_1 + (D2+N2)*abu_2) # long
    (ratio_0*abu_0 + ratio_1*abu_1 + ratio_2*abu_2) / sum([abu_0, abu_1, abu_2]) # short
    abundance_colname: "Intensity", "MS/MS Count"
    ca. 14 times faster :) than conditional indexing
    :param group: DataFrame(grouped)
    :param abundance_colname: String
    :return: DataFrame(grouped
    """
    group["perc_DeamN"] = 100.0 * sum(group["ratio_N2D"] * group[abundance_colname]) / sum(group[abundance_colname])
    group["perc_DeamQ"] = 100.0 * sum(group["ratio_Q2E"] * group[abundance_colname]) / sum(group[abundance_colname])
    return group

# def avg_PSM2Pep_Pep2Prot_Prot2Rawfile(df, mqc):
#     """
#     this is how we do it!
#     average PSMs to Peptides
#     average Peptides to Proteins
#     average Proteins to RawFile
#     :param df: DataFrame (with perc_DeamN col --> calc deamidation per RawFile, Sequence and Charge)
#     :return: DataFrame
#     """
#     colname_RawFile = mqc.check_if_colname_or_similar_exists("Raw file")
#     colname_LeadingRazorProtein = mqc.check_if_colname_or_similar_exists("Leading razor protein")
#     # subset columns needed
#     dfx = df[[colname_RawFile, colname_LeadingRazorProtein, "Sequence", "Charge", "perc_DeamN", "perc_DeamQ"]].drop_duplicates()
#     # average PSMs to Peptides (merge charge states of same sequence)
#     dfx = dfx.groupby([colname_RawFile, colname_LeadingRazorProtein, "Sequence"])[["perc_DeamN", "perc_DeamQ"]].agg("mean")
#     dfx = dfx.reset_index()
#     # averge Peptides to Proteins (merge sequences of same protein)
#     dfx = dfx.groupby([colname_RawFile, colname_LeadingRazorProtein])[["perc_DeamN", "perc_DeamQ"]].agg("mean")
#     dfx = dfx.reset_index()
#     # average Proteins to RawFiles
#     # first confidence interval then mean and std
#     dfx = dfx.groupby([colname_RawFile])[["perc_DeamN", "perc_DeamQ"]].agg(["mean", "std"])
#     dfx = dfx.reset_index()
#     dfx.columns = ['_'.join(col).strip() for col in dfx.columns.values]
#     dfx.columns = [colname_RawFile, "perc_DeamN", "perc_DeamN_std", "perc_DeamQ", "perc_DeamQ_std"]
#     return dfx

def avg_PSM2Pep(df, colname_RawFile, colname_LeadingRazorProtein, mean_or_median="mean"):
    """
    average PSMs to Peptides
    :param df: DataFrame (with perc_DeamN col --> calc deamidation per RawFile, Sequence and Charge)
    :param colname_RawFile: String
    :param colname_LeadingRazorProtein: String
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: DataFrame (with reduced shape)
    """
    # subset columns needed
    dfx = df[[colname_RawFile, colname_LeadingRazorProtein, "Sequence", "Charge", "perc_DeamN", "perc_DeamQ"]].drop_duplicates()
    # average PSMs to Peptides (merge charge states of same sequence)
    dfx = dfx.groupby([colname_RawFile, colname_LeadingRazorProtein, "Sequence"])[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)
    dfx = dfx.reset_index()
    return dfx

def avg_Pep2Prot(df, colname_RawFile, colname_LeadingRazorProtein, mean_or_median="mean"):
    """
    average Peptides to Proteins
    :param df: DataFrame (with perc_DeamN col --> calc deamidation per RawFile, Sequence and Charge)
    :param colname_RawFile: String
    :param colname_LeadingRazorProtein: String
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: DataFrame (with reduced shape)
    """
    # averge Peptides to Proteins (merge sequences of same protein)
    dfx = df.groupby([colname_RawFile, colname_LeadingRazorProtein])[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)
    dfx = dfx.reset_index()
    return dfx

def avg_Prot2RawFile(df, colname_RawFile, std=False, mean_or_median="mean"):
    """
    average Proteins to RawFile
    :param df: DataFrame (with perc_DeamN col --> calc deamidation per RawFile, Sequence and Charge)
    :param colname_RawFile: String
    :param std: Bool (flag for standard deviation)
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: DataFrame (with reduced shape)
    """
    # average Proteins to RawFiles
    # first confidence interval then mean and std
    if std:
        dfx = df.groupby([colname_RawFile])[["perc_DeamN", "perc_DeamQ"]].agg([mean_or_median, "std"])
        dfx = dfx.reset_index()
        dfx.columns = ['_'.join(col).strip() for col in dfx.columns.values]
        dfx.columns = [colname_RawFile, "perc_DeamN", "perc_DeamN_std", "perc_DeamQ", "perc_DeamQ_std"]
    else:
        dfx = df.groupby([colname_RawFile])[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)
        dfx = dfx.reset_index()
    return dfx

def deam_per_RawFile_with_CI(df, mqc, ci=95, sampling="peptides", num_bootstraps=1000, groupby_="Raw file", mean_or_median="mean", average_method="Peptides_2_Proteins_per_groupby_"):
    """
    first collapse PSMs to Peptides by averaging
    then either bootstrap sampling either of Peptides
    or collapse Pepties to Proteins and bootstrap Proteins
    :param df: DataFrame (evidence.txt with deam on PSM level)
    :param mqc: MQcolnames Instance (check old vs. new colnames in tools-module)
    :param ci: Float or Int (Confidence Interval)
    :param sampling: String ("peptides" or "proteins")
    :param num_bootstraps: Int
    :param groupby_: String ('Raw file' or 'Experiment') --> if Experiment names are set
                properly, then using Experiment will result in useful names for the plot
    :param mean_or_median: String(Flag to select mean/average or median)
    :param average_method: String("Peptides_2_Proteins_per_groupby_", "Peptides_per_groupby_". Flag to select average Peptides to Proteins and average Proteins per RawFile/Folder OR average Peptides per RawFile/Folder)
    :return: DataFrame (Raw File, N_Q, percDeam, CI_low, CI_up)
    """
    ci_low = (100 - ci) / 2.0
    ci_up = ci + ci_low

    groupby_ = mqc.check_if_colname_or_similar_exists(groupby_)
    colname_LeadingRazorProtein = mqc.check_if_colname_or_similar_exists("Leading razor protein")

    # average PSMs to Peptides
    dfx = avg_PSM2Pep(df, groupby_, colname_LeadingRazorProtein, mean_or_median=mean_or_median)
    if sampling == "peptides":
        # sampling Peptides for every RawFile
        pass
    elif sampling == "proteins":
        # sampling proteins for every RawFile
        # average Peptides to Proteins
        dfx = avg_Pep2Prot(dfx, groupby_, colname_LeadingRazorProtein, mean_or_median=mean_or_median)
    groupby_ = mqc.check_if_colname_or_similar_exists(groupby_)
    grouped = dfx.groupby(groupby_)
    name_list, n_lower_list, n_upper_list, q_lower_list, q_upper_list = [], [], [], [], []
    for name, group in grouped:
        print(name)
        group_index = group.index.tolist()
        group_n = []
        group_q = []
        for random_index in yield_random_combination_with_replacement(group_index, num_bootstraps):
            dft = dfx.loc[random_index, :]
            ### sanity check
            cond = dft.index == random_index
            assert cond.all() == True
            if average_method == "Peptides_2_Proteins_per_groupby_":
                ser = bootstrap_RFgrouped_avg_Pep2Prot_avg_Prot(dft, colname_LeadingRazorProtein, mean_or_median=mean_or_median)
            ### average Peptides per RawFile/Folder
            elif average_method == "Peptides_per_groupby_":
                ser = bootstrap_RFgrouped_avg_Pep(dft, mean_or_median=mean_or_median)
            group_n.append(ser["perc_DeamN"])
            group_q.append(ser["perc_DeamQ"])
        group_n = np.array(group_n)
        group_q = np.array(group_q)
        group_n = group_n[np.isfinite(group_n)]
        group_q = group_q[np.isfinite(group_q)]
        if group_n.size > 0:
            group_n_lower, group_n_upper = np.percentile(group_n, ci_low), np.percentile(group_n, ci_up)
        else:
            group_n_lower, group_n_upper = np.nan, np.nan
        if group_q.size > 0:
            group_q_lower, group_q_upper = np.percentile(group_q, ci_low), np.percentile(group_q, ci_up)
        else:
            group_q_lower, group_q_upper = np.nan, np.nan
        name_list.append(name)
        n_lower_list.append(group_n_lower)
        n_upper_list.append(group_n_upper)
        q_lower_list.append(group_q_lower)
        q_upper_list.append(group_q_upper)
    df2merge = pd.DataFrame({groupby_: name_list, "DeamN_low": n_lower_list,
                             "DeamN_up": n_upper_list, "DeamQ_low": q_lower_list, "DeamQ_up": q_upper_list})
    if average_method == "Peptides_2_Proteins_per_groupby_":
        dfx = avg_Pep2Prot(dfx, groupby_, colname_LeadingRazorProtein, mean_or_median=mean_or_median)
        dfx = avg_Prot2RawFile(dfx, groupby_, std=False, mean_or_median=mean_or_median)
    elif average_method == "Peptides_per_groupby_":
        dfx = dfx.groupby(groupby_)[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)
        dfx = dfx.reset_index()
    dfm = pd.merge(dfx, df2merge, how="outer")

    dfm1 = pd.melt(dfm, id_vars=[groupby_], value_vars=["perc_DeamN", "perc_DeamQ"], var_name="N_Q", value_name="percDeam")
    dfm1.loc[dfm1["N_Q"] == "perc_DeamN", "N_Q"] = "N"
    dfm1.loc[dfm1["N_Q"] == "perc_DeamQ", "N_Q"] = "Q"

    dfm2 = pd.melt(dfm, id_vars=[groupby_], value_vars=["DeamN_low", "DeamQ_low"], var_name="N_Q", value_name="CI_low")
    dfm3 = pd.melt(dfm, id_vars=[groupby_], value_vars=["DeamN_up", "DeamQ_up"], var_name="N_Q", value_name="CI_up")

    dfm2.loc[dfm2["N_Q"] == "DeamN_low", "N_Q"] = "N"
    dfm2.loc[dfm2["N_Q"] == "DeamQ_low", "N_Q"] = "Q"

    dfm3.loc[dfm3["N_Q"] == "DeamN_up", "N_Q"] = "N"
    dfm3.loc[dfm3["N_Q"] == "DeamQ_up", "N_Q"] = "Q"

    dfm_p1 = pd.merge(dfm1, dfm2, how='outer')
    dfmx = pd.merge(dfm_p1, dfm3, how='outer')
    return dfmx

def yield_random_combination_with_replacement(iterable, n_times):
    """
    should be even more memory efficient :)
    random selection with replacement from iterable of len(iterable) --> produces results
    then randomly select n_times from results (since results are sorted and
    single values are always recurring at beginning of results)
    :param iterable: Iterable
    :param n_times: Int (number of random combinations)
    :return: Generator (with yields iterable)
    """
    for _ in range(0, n_times):
        yield np.random.choice(iterable, len(iterable), replace=True)

# def bootstrap_RFgrouped_avg_Prot(df):
#     """
#     calculate deamidation for N and Q by averaging Proteins
#     expects a DataFrame of one RawFile
#     :param df: DataFrame(evidence of one RawFile)
#     :return: Series
#     """
#     return df[["perc_DeamN", "perc_DeamQ"]].mean()

def bootstrap_RFgrouped_avg_Pep2Prot_avg_Prot(df, colname_LeadingRazorProtein, mean_or_median="mean"):
    """
    calculate deamidation for N and Q by averaging Peptides to Proteins
    and averaging Proteins
    expects a DataFrame of one RawFile
    :param df: DataFrame(evidence of one RawFile)
    :param colname_LeadingRazorProtein: String
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: Series
    """
    dft = df.groupby(colname_LeadingRazorProtein)[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)
    return dft[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)

def bootstrap_RFgrouped_avg_Pep(df, mean_or_median="mean"):
    """
    calculate deamidation for N and Q by averaging Peptides to Proteins
    and averaging Proteins
    expects a DataFrame of one RawFile
    :param df: DataFrame(evidence of one RawFile)
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: Series
    """
    return df[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median)

def run(fn_evidence, abundance_colname, ci, sampling, num_bootstraps, output_dir, colname_proteins):
    ###### use a MaxQuant evidence.txt as input
    ## merging data per RawFile
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df = pd.read_csv(fn_evidence, sep='\t')
    df = add_num_N_Q_cols(df)
    df = add_num_N2D_Q2E_ratios(df)
    ## calc deamidation per RawFile, Sequence and Charge
    df = df.groupby(by=["Raw file", "Sequence", "Charge"], axis=0).apply(percent_deamidation, abundance_colname)
    mqc = MQcolnames(df)
    dfmx = deam_per_RawFile_with_CI(df, mqc, ci=ci, sampling=sampling, num_bootstraps=num_bootstraps)
    fn_out = os.path.join(output_dir, "Deamidation_RawFile.txt")
    dfmx.to_csv(fn_out, sep='\t', header=True, index=False)
    ## merging data per Protein
    colname_RawFile = mqc.check_if_colname_or_similar_exists("Raw file")
    colname_LeadingRazorProtein = mqc.check_if_colname_or_similar_exists(colname_proteins)
    dfx = avg_PSM2Pep(df, colname_RawFile, colname_LeadingRazorProtein)
    dfx = avg_Pep2Prot(dfx, colname_RawFile, colname_LeadingRazorProtein)
    fn_out = os.path.join(output_dir, "Deamidation_Protein.txt")
    dfx.to_csv(fn_out, sep='\t', header=True, index=False)
    ############################################################################################################################################


def error_(parser):
    sys.stderr.write("The arguments passed are invalid.\nPlease check the input parameters.\n\n")
    parser.print_help()
    sys.exit(2)

if __name__ == '__main__':
    cmd_line = True

    if not cmd_line:
        ############################################################################################################################################
        ### input options
        fn_evidence = r"/Users/dblyon/Downloads/evidence.txt"
        output_dir = r"/Users/dblyon/Downloads/Output"
        abundance_colname = "Intensity"
        sampling = "peptides"
        num_bootstraps = 1000
        ci = 95
        colname_proteins = "Leading razor protein"

    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("-e", "--evidence", help="MaxQuant evidence.txt file (absolute path)", type=str)
        parser.add_argument("-o", "--outputdirectory", help="Output directory for results (absolute path)", type=str)
        parser.add_argument("-a", "--abundancecolname", help="Abundance column name ", type=str, default="Intensity")
        parser.add_argument("-s", "--samplepeptidesorproteins", help="Sample 'peptides' or 'proteins' for Confidence Interval calculation (default='peptides')", type=str, default="peptides")
        parser.add_argument("-ci", "--confidenceinterval", help="Confidence Interval (default=95)", type=int, default=95)
        parser.add_argument("-b", "--numbootstraps", help="Number of bootstrap iterations (default=1000)", type=int, default=1000)
        parser.add_argument("-c", "--colnameproteins", help="Column name for Protein identifiers (default='Leading razor protein')", type=str, default="Leading razor protein")

        args = parser.parse_args()
        fn_evidence = args.evidence
        abundance_colname = args.abundancecolname
        ci = args.confidenceinterval
        sampling = args.samplepeptidesorproteins
        num_bootstraps = args.numbootstraps
        output_dir = args.outputdirectory
        colname_proteins = args.colnameproteins

        # for arg in sorted(vars(args)):
        #     if getattr(args, arg) is None:
        #         error_(parser)
        
        if fn_evidence is None:
            error_(parser)

        if output_dir is None:
            output_dir = os.path.dirname(fn_evidence)

        print("#" * 80)
        for arg in sorted(vars(args)):
            print(arg, ": ", getattr(args, arg))

    run(fn_evidence, abundance_colname, ci, sampling, num_bootstraps, output_dir, colname_proteins)

