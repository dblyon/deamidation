# Readme for deamidation
This is a command-line script to calculate deamidation rates of peptides to estimate the deamidation of Asparagine to Aspartic acid and Glutamine to Glutamic acid per RawFile.
The program expects an "evidence.txt" file from <a href="http://www.biochem.mpg.de/5111795/maxquant">MaxQuant</a> as input, and will generate four tab delimited text files as output.

- Deamidation.txt (lists the RawFiles and their respective deamidation for N and Q, as mean, standard deviation, 95% confidence lower and upper limit)
- Number_of_Peptides_per_RawFile.txt (number of peptides available for calculation of deamidation for N and for Q)
- Bootstrapped_values.txt (all the deamidation percentages calculated by e.g. 1000 bootstrap iterations, which are subsequently used to calculate the mean, std, and CI for shown in "Deamidation.txt")
- Protein_deamidation.txt (deamidation on the protein level, to be used with restraint since there usually are few data to acquire meaningful results, therefore no bootstrapping is applied)
- Deamidation on the protein level can also be calculated, although this is NOT recommended as a default analysis (since there often are too little data for this to be meaningful/robust). Use "protein_bootstraps" flag to additionally generate the following files which are analogous the files above but grouped by RawFile and Proteins instead of only by RawFile.
    * Bootstrapped_values_proteins.txt
    * Deamidation_proteins.txt
    * Number_of_Peptides_per_Protein_per_RawFile.txt

The "R_plot.R" script contains code that uses "Deamidation.txt" and "Number_of_Peptides_per_RawFile.txt" to create a BarPlot with STD as error bars, and includes the number of peptides used for the respective deamidation calculation.  

# Citation
Please cite the following publication when using this tool:
Mackie et al., Palaeoproteomic Profiling of Conservation Layers on a 14th Century Italian Wall Painting., Angew Chem Int Ed Engl., 2018
- https://doi.org/10.1002/anie.201713020
- https://www.ncbi.nlm.nih.gov/pubmed/29603563

# Prerequisites
    - Python 2.x
    - Python packages: numpy, pandas

# Usage
type "-h" or "--help" for help.
e.g. type the following command in the terminal to see options.
<code>python deamidation.py -h</code>

## example for OSX
<code>python deamidation.py -e /Directory/evidence.txt -o /Directory/Output</code>

## example for Windows
<code>python deamidation.py -e C:\Directory\evidence.txt -o C:\Directory\Output</code>

# Support & Questions
feel free to contact <david.lyon@uzh.ch> for help


## Deamidation calculation explanation
See 'Deamidation_calculation.pdf' 
Peptide level deamidation rates are calculated separately for asparagine and glutamine taking relative abundance values into account, resulting in a deamidation rate between 0 and 1 for each unique peptide sequence and charge state. The various charge states are collapsed to the peptide level by calculating the mean, resulting in a separate deamidation rate for asparagine as well as glutamine, respectively. The peptide level deamidation rates are sampled with replacement (bootstrapped) 1000 times and therefrom the mean and the standard deviation calculated (in order to gain an estimate of the error of the calculation). Generally, the higher the number of peptides detected with the potential to be deamidated, the more accurate the measurement and thus the lower the error will be. Nevertheless, the method is intended for comparative purposes and not as an absolute measure of deamidation. 

