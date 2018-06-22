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

def run(fn_evidence, abundance_colname, ci, sampling, num_bootstraps, output_dir, colname_proteins, protein_bootstraps=False):
    """
    :param fn_evidence: String (absolute path to 'evidence.txt' file)
    :param abundance_colname: String (column name of abundance data e.g. 'Intensity')
    :param ci: Integer (Confidence Interval e.g. 95)
    :param sampling: String ("peptides" or "proteins")
    :param num_bootstraps: Integer (Number of Bootstrap iterations e.g. 1000)
    :param output_dir: String (absolute path to output directory)
    :param colname_proteins: String (column name of Protein Accession Numbers/Identifierse.g. 'Leading razor protein')
    :param protein_bootstraps: Boolean(flag to get additional output --> bootstrap Peptides per Protein per RawFile
        (usually there are not enough data for this to be meaningful and thus this can lead to misleading results,
        this is NOT recommended as a default analysis)
    :return: None
    """
    ###### use a MaxQuant evidence.txt as input
    ## merging data per RawFile
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Reading file")
    df = pd.read_csv(fn_evidence, sep='\t')
    mqc = MQcolnames(df)
    colname_proteins = mqc.check_if_colname_or_similar_exists(colname_proteins)
    colname_RawFile = mqc.check_if_colname_or_similar_exists("Raw file")
    print("Calculating N2D and Q2E ratios")
    df = add_num_N_Q_cols(df)
    df = add_num_N2D_Q2E_ratios(df)
    ## calc deamidation per RawFile, Sequence and Charge
    df = df.groupby(by=[colname_RawFile, "Sequence", "Charge"], axis=0).apply(percent_deamidation, abundance_colname)

    fn_out_num_peptides = os.path.join(output_dir, "Number_of_Peptides_per_RawFile.txt")
    print("Writing {}".format(fn_out_num_peptides))
    _ = number_of_N_Q_peptides_with_DeamPerc(df, colname_RawFile, fn_out_num_peptides)

    ### Protein level deamidation --> without CI
    dfp = avg_PSM2Pep(df, colname_RawFile, colname_proteins)
    dfp = avg_Pep2Prot(dfp, colname_RawFile, colname_proteins)
    fn_out = os.path.join(output_dir, "Protein_deamidation.txt")
    print("Writing {}".format(fn_out))
    dfp.to_csv(fn_out, sep='\t', header=True, index=False)

    print("Bootstrapping")
    fn_out_bootstrapped = os.path.join(output_dir, "Bootstrapped_values.txt")
    df_bs = deam_per_RawFile_returnBootstrappedVals(df, mqc, colname_proteins, sampling=sampling, num_bootstraps=num_bootstraps, groupby_=colname_RawFile, mean_or_median="mean", average_method="Peptides_per_groupby_")

    df_bs = df_bs.reset_index(drop=True)
    print("Writing results to: {}".format(fn_out_bootstrapped))
    df_bs.to_csv(fn_out_bootstrapped, sep='\t', header=True, index=False)

    fn_out = os.path.join(output_dir, "Deamidation.txt")
    print("Writing {}".format(fn_out))
    calculate_mean_and_CIs(df_bs, ci, colname_RawFile, fn_out)

    if protein_bootstraps:
        print("Calculating protein level deamidation per RawFile using bootstraps")
        df_protein_boots = deam_per_RawFile_returnBootstrappedVals_bootstrapPeptidesPerProtein(df, mqc, colname_proteins=colname_proteins, num_bootstraps=num_bootstraps, groupby_=colname_RawFile)
        fn_out_bootstrapped_proteins = os.path.join(output_dir, "Bootstrapped_values_proteins.txt")
        print("Writing results to: {}".format(fn_out_bootstrapped_proteins))
        df_protein_boots.to_csv(fn_out_bootstrapped_proteins, sep='\t', header=True, index=False)

        fn_out_proteins = os.path.join(output_dir, "Deamidation_proteins.txt")
        print("Writing {}".format(fn_out_proteins))
        cond_not_NaN = df_protein_boots["percDeam"].notnull() # need to be removed for CI calculation, but left them in deamidation vals
        calculate_mean_and_CIs_proteins(df_protein_boots[cond_not_NaN], ci, colname_RawFile, colname_proteins, fn_out_proteins)

        fn_out_num_peptides_proteins = os.path.join(output_dir, "Number_of_Peptides_per_Protein_per_RawFile.txt")
        print("Writing {}".format(fn_out_num_peptides_proteins))
        _ = number_of_N_Q_peptides_with_DeamPerc_per_Protein(df, colname_RawFile, colname_proteins, fn_out_num_peptides_proteins)


def deam_per_RawFile_returnBootstrappedVals(df, mqc, colname_proteins="Leading razor protein", sampling="peptides", num_bootstraps=1000, groupby_="Raw file", mean_or_median="mean", average_method="Peptides_per_groupby_"):
    """
    calculate deamidation-MEAN, but save values for BoxPlot
    additionally calculate median and add to plot

    first collapse PSMs to Peptides by averaging
    then either bootstrap sampling either of Peptides
    or collapse Pepties to Proteins and bootstrap Proteins
    :param df: DataFrame (evidence.txt with deam on PSM level)
    :param mqc: MQcolnames Instance (check old vs. new colnames in tools-module)
    :param colname_proteins: String (column name e.g. 'Leading razor protein')
    :param sampling: String ("peptides" or "proteins")
    :param num_bootstraps: Int
    :param groupby_: String ('Raw file' or 'Experiment') --> if Experiment names are set
                properly, then using Experiment will result in useful names for the plot
    :param mean_or_median: String(Flag to select mean/average or median)
    :param average_method: String(Flag to select average Peptides to Proteins and average Proteins per RawFile/Folder OR average Peptides per RawFile/Folder)
    :return: DataFrame (Raw File, N_Q, percDeam, CI_low, CI_up)
    """
    colname_RawFile = mqc.check_if_colname_or_similar_exists(groupby_)
    colname_proteins = mqc.check_if_colname_or_similar_exists(colname_proteins)

    # average PSMs to Peptides
    dfx = avg_PSM2Pep(df, colname_RawFile, colname_proteins, mean_or_median=mean_or_median)
    if sampling == "peptides":
        # sampling Peptides for every RawFile
        pass
    elif sampling == "proteins":
        # sampling proteins for every RawFile
        # average Peptides to Proteins
        dfx = avg_Pep2Prot(dfx, colname_RawFile, colname_proteins, mean_or_median=mean_or_median)
    colname_RawFile = mqc.check_if_colname_or_similar_exists(groupby_)
    grouped = dfx.groupby(colname_RawFile)
    df_list = []
    for name, group in grouped:
        print(name)
        df_temp = pd.DataFrame()
        # df_temp[groupby_] = groupby_
        group_index = group.index.tolist()
        group_n = []
        group_q = []
        for random_index in yield_random_combination_with_replacement(group_index, num_bootstraps):
            dft = dfx.loc[random_index, :]
            ### sanity check
            cond = dft.index == random_index
            assert cond.all() == True
            ### average Peptides to Proteins and average Proteins per RawFile/Folder
            if average_method == "Peptides_2_Proteins_per_groupby_":
                ser = bootstrap_RFgrouped_avg_Pep2Prot_avg_Prot(dft, colname_proteins, mean_or_median=mean_or_median)
            ### average Peptides per RawFile/Folder
            elif average_method == "Peptides_per_groupby_":
                ser = bootstrap_RFgrouped_avg_Pep(dft, mean_or_median=mean_or_median)
            else:
                print("Method doesn't exist")
                raise StopIteration
            group_n.append(ser["perc_DeamN"])
            group_q.append(ser["perc_DeamQ"])
        group_n = np.array(group_n)
        group_q = np.array(group_q)
        group_n = group_n[np.isfinite(group_n)]
        group_q = group_q[np.isfinite(group_q)]
        if group_n.size == 0:
            group_n = np.nan
        if group_q.size == 0:
            group_q = np.nan
        df_temp["perc_DeamN"] = pd.Series(group_n)
        df_temp["perc_DeamQ"] = pd.Series(group_q)
        df_temp[groupby_] = name
        df_list.append(df_temp)
    dfm = pd.concat(df_list)
    dfm = pd.melt(dfm, id_vars=[groupby_], value_vars=["perc_DeamN", "perc_DeamQ"], var_name="N_Q", value_name="percDeam")
    dfm.loc[dfm["N_Q"] == "perc_DeamN", "N_Q"] = "N"
    dfm.loc[dfm["N_Q"] == "perc_DeamQ", "N_Q"] = "Q"
    return dfm

def deam_per_RawFile_returnBootstrappedVals_bootstrapPeptidesPerProtein(df, mqc, colname_proteins="Leading razor protein", num_bootstraps=1000, groupby_="Raw file", mean_or_median="mean"):
    """
    calculate deamidation-MEAN, but save values for BoxPlot
    additionally calculate median and add to plot

    first collapse PSMs to Peptides by averaging
    then bootstrap sampling of Peptides per Proteins --> CAREFUL THIS IS NOT RECOMMENDED AS A DEFAULT METHOD
    :param df: DataFrame (evidence.txt with deam on PSM level)
    :param mqc: MQcolnames Instance (check old vs. new colnames in tools-module)
    :param colname_proteins: String (column name e.g. 'Leading razor protein')
    :param num_bootstraps: Int
    :param groupby_: String ('Raw file' or 'Experiment') --> if Experiment names are set
                properly, then using Experiment will result in useful names for the plot
    :param mean_or_median: String(Flag to select mean/average or median)
    :return: DataFrame (Raw File, N_Q, percDeam, CI_low, CI_up)
    """
    colname_RawFile = mqc.check_if_colname_or_similar_exists(groupby_)
    colname_proteins = mqc.check_if_colname_or_similar_exists(colname_proteins)
    # average PSMs to Peptides
    dfx = avg_PSM2Pep(df, colname_RawFile, colname_proteins, mean_or_median=mean_or_median)
    colname_RawFile = mqc.check_if_colname_or_similar_exists(groupby_)
    grouped = dfx.groupby(colname_RawFile)
    df_list = []
    for name, group in grouped:
        print(name)
        group_index = group.index.tolist()
        for random_index in yield_random_combination_with_replacement(group_index, num_bootstraps):
            dft = dfx.loc[random_index, :]
            ### sanity check
            cond = dft.index == random_index
            assert cond.all() == True
            df_deam_per_protein = dft.groupby(colname_proteins)[["perc_DeamN", "perc_DeamQ"]].agg(mean_or_median).reset_index()
            df_deam_per_protein[groupby_] = name
            df_list.append(df_deam_per_protein)
    dfm = pd.concat(df_list)
    dfm = pd.melt(dfm, id_vars=[groupby_, colname_proteins], value_vars=["perc_DeamN", "perc_DeamQ"], var_name="N_Q", value_name="percDeam")
    dfm.loc[dfm["N_Q"] == "perc_DeamN", "N_Q"] = "N"
    dfm.loc[dfm["N_Q"] == "perc_DeamQ", "N_Q"] = "Q"
    return dfm

def number_of_N_Q_peptides_with_DeamPerc(df, colname_RawFile, fn_out):
    """
    number of Peptides for which a deamidation percentage exists, and that can therefore be used for bootstrapping.
    """
    ser_list = []

    # already calculated previously
    # # calc deamidation per RawFile, Sequence and Charge
    # df = df.groupby(by=[colname_RawFile, "Sequence", "Charge"], axis=0).apply(percent_deamidation, abundance_colname)

    grouped = df[[colname_RawFile, "perc_DeamN", "perc_DeamQ"]].groupby(colname_RawFile)
    for name, group in grouped:
        ser = group[["perc_DeamN", "perc_DeamQ"]].apply(lambda df: sum(df.notnull()))
        ### doesn't work on groupED object but on group: grouped[["perc_DeamN", "perc_DeamQ"]].apply(lambda df: sum(df.notnull()))
        ser.name = name
        ser_list.append(ser)
    df = pd.concat(ser_list, axis=1).T
    df = df.reset_index()
    df.columns = [colname_RawFile, "num_N_peptides", "num_Q_peptides"]
    dfm = df.melt(id_vars=[colname_RawFile], value_vars=["num_N_peptides", "num_Q_peptides"], var_name="N_Q", value_name="num_peptides")
    cond_n = dfm["N_Q"] == "num_N_peptides"
    dfm.loc[cond_n, "N_Q"] = "N"
    dfm.loc[-cond_n, "N_Q"] = "Q"
    print("Writing results to: {}".format(fn_out))
    dfm.to_csv(fn_out, sep='\t', header=True, index=False)
    return dfm

def number_of_N_Q_peptides_with_DeamPerc_per_Protein(df, colname_RawFile, colname_proteins, fn_out):
    """
    number of Peptides for which a deamidation percentage exists, and that can therefore be used for bootstrapping.
    """
    ser_list = []
    grouped = df[[colname_RawFile, colname_proteins, "perc_DeamN", "perc_DeamQ"]].groupby([colname_RawFile, colname_proteins])
    for name, group in grouped:
        ser = group[["perc_DeamN", "perc_DeamQ"]].apply(lambda df: sum(df.notnull()))
        ### doesn't work on groupED object but on group: grouped[["perc_DeamN", "perc_DeamQ"]].apply(lambda df: sum(df.notnull()))
        ser.name = name
        ser_list.append(ser)
    df = pd.concat(ser_list, axis=1).T
    df = df.reset_index()
    df.columns = [colname_RawFile, colname_proteins, "num_N_peptides", "num_Q_peptides"]
    dfm = df.melt(id_vars=[colname_RawFile, colname_proteins], value_vars=["num_N_peptides", "num_Q_peptides"], var_name="N_Q", value_name="num_peptides")
    cond_n = dfm["N_Q"] == "num_N_peptides"
    dfm.loc[cond_n, "N_Q"] = "N"
    dfm.loc[-cond_n, "N_Q"] = "Q"
    print("Writing results to: {}".format(fn_out))
    dfm.to_csv(fn_out, sep='\t', header=True, index=False)
    return dfm



def calculate_mean_and_CIs(df, ci, groupby_, fn_out):
    ci_low = int((100 - ci) / 2.0)
    ci_up = int(ci + ci_low)

    percentile_low = df.groupby([groupby_, "N_Q"])["percDeam"].apply(np.percentile, ci_low).reset_index()
    percentile_low = percentile_low.rename(columns={"percDeam": "CI_low"})
    percentile_up = df.groupby([groupby_, "N_Q"])["percDeam"].apply(np.percentile, ci_up).reset_index()
    percentile_up = percentile_up.rename(columns={"percDeam": "CI_up"})
    df_ci = pd.merge(percentile_low, percentile_up, how='outer')

    dfx = df.groupby([groupby_, "N_Q"])["percDeam"].agg(["mean", "std"]).reset_index()
    df = pd.merge(dfx, df_ci, how='outer')
    df.to_csv(fn_out, sep='\t', header=True, index=False)

def calculate_mean_and_CIs_proteins(df, ci, groupby_, colname_proteins, fn_out):
    ci_low = int((100 - ci) / 2.0)
    ci_up = int(ci + ci_low)
    percentile_low = df.groupby([groupby_, colname_proteins, "N_Q"])["percDeam"].apply(np.percentile, ci_low).reset_index()
    percentile_low = percentile_low.rename(columns={"percDeam": "CI_low"})
    percentile_up = df.groupby([groupby_, colname_proteins, "N_Q"])["percDeam"].apply(np.percentile, ci_up).reset_index()
    percentile_up = percentile_up.rename(columns={"percDeam": "CI_up"})
    df_ci = pd.merge(percentile_low, percentile_up, how='outer')

    dfx = df.groupby([groupby_, colname_proteins, "N_Q"])["percDeam"].agg(["mean", "std"]).reset_index()
    df = pd.merge(dfx, df_ci, how='outer')
    df.to_csv(fn_out, sep='\t', header=True, index=False)
    return df

def error_(parser):
    sys.stderr.write("The arguments passed are invalid.\nPlease check the input parameters.\n\n")
    parser.print_help()
    sys.exit(2)

if __name__ == '__main__':
    cmd_line = True

    if not cmd_line:
        ### input options
        fn_evidence = r"/Users/dblyon/Downloads/evidence.txt"
        output_dir = r"/Users/dblyon/Downloads/Output"
        abundance_colname = "Intensity"
        sampling = "peptides"
        num_bootstraps = 1000
        ci = 95
        colname_proteins = "Leading razor protein"
        protein_bootstraps = True


    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("-e", "--evidence", help="MaxQuant evidence.txt file (absolute path)", type=str)
        parser.add_argument("-o", "--outputdirectory", help="Output directory for results (absolute path)", type=str)
        parser.add_argument("-a", "--abundancecolname", help="Abundance column name ", type=str, default="Intensity")
        parser.add_argument("-s", "--samplepeptidesorproteins", help="Sample 'peptides' or 'proteins' for Confidence Interval calculation (default='peptides')", type=str, default="peptides")
        parser.add_argument("-ci", "--confidenceinterval", help="Confidence Interval (default=95)", type=int, default=95)
        parser.add_argument("-b", "--numbootstraps", help="Number of bootstrap iterations (default=1000)", type=int, default=1000)
        parser.add_argument("-c", "--colnameproteins", help="Column name for Protein identifiers (default='Leading razor protein')", type=str, default="Leading razor protein")
        parser.add_argument("-p", "--protein_bootstraps", help="Bootstrap Peptides per Protein per RawFile. Usually there are not enough data for this to be meaningful and thus, this can lead to misleading results. This is NOT recommended as a default analysis. (default=False)", type=bool, default=False)

        args = parser.parse_args()
        fn_evidence = args.evidence
        abundance_colname = args.abundancecolname
        ci = args.confidenceinterval
        sampling = args.samplepeptidesorproteins
        num_bootstraps = args.numbootstraps
        output_dir = args.outputdirectory
        colname_proteins = args.colnameproteins
        protein_bootstraps = args.protein_bootstraps

        if fn_evidence is None:
            error_(parser)

        if output_dir is None:
            output_dir = os.path.dirname(fn_evidence)

        print("#" * 80)
        print("SETTINGS:")
        for arg in sorted(vars(args)):
            print(arg, ": ", getattr(args, arg))
        print("#" * 80)
    run(fn_evidence, abundance_colname, ci, sampling, num_bootstraps, output_dir, colname_proteins, protein_bootstraps)

