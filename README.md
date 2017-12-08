# Readme for deamidation
This is a command-line script to calculate deamidation rates of peptides to estimate the deamidation of Asparagine to Aspartic acid and Glutamine to Glutamic acid per RawFile.
The program expects an "evidence.txt" file from <a href="http://www.biochem.mpg.de/5111795/maxquant">MaxQuant</a> as input, and will generate two tab delimited text files as output. 

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
<code>python deamidation.py -evidence /Directory/evidence.txt -o /Directory/Output</code>

## example for Windows
<code>python deamidation.py -fasta C:\Directory\evidence.txt -o C:\Directory\Output</code>

# Support
feel free to contact <david.lyon@cpr.ku.dk> for help