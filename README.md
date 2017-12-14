# Readme for deamidation
This is a command-line script to calculate deamidation rates of peptides to estimate the deamidation of Asparagine to Aspartic acid and Glutamine to Glutamic acid per RawFile.
The program expects an "evidence.txt" file from <a href="http://www.biochem.mpg.de/5111795/maxquant">MaxQuant</a> as input, and will generate four tab delimited text files as output.
    - Deamidation.txt (lists the RawFiles and their respective deamidation for N and Q, as mean, standard deviation, 95% confidence lower and upper limit)
    - Number_of_Peptides_per_RawFile.txt (number of peptides available for calculation of deamidation for N and for Q)
    - Bootstrapped_values.txt (all the deamidation percentages calculated by e.g. 1000 bootstrap iterations, which are subsequently used to calculate the mean, std, and CI for shown in "Deamidation.txt")
    - Protein_deamidation.txt (deamidation on the protein level, to be used with restraint since there usually are few data to acquire meaningful results, therefore no bootstrapping is applied)
    
# Citation
Publication will be referenced here.

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

# Support
feel free to contact <david.lyon@cpr.ku.dk> for help