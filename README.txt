In silico README

### Dependencies ###

python3
python-dev (sudo apt-get install python-dev), 
regex (https://launchpad.net/ubuntu/+archive/primary/+files/python-regex_0.1.20170117.orig.tar.gz)

### goal ###

MLVA_finder.py is a python script designed to do Mutli loci VNTR analysis
(VNTR stand for Variable Number of Tandem Repeats ).
The goal of insilico is to find from fasta file(s), sequences of tandem repeat 
and assign an MLVA value for each locus.
This script use a list of primers to get loci's sequences.
The number of mismatch can be set in the command line. 

### input ###

-i or --input take the path of the directory containing all fasta file which will be used by the script 
-o or --output take the path of the results directory where results will be saved
-p or --primer take the path of the primers list in a csv file (delimited by '\t', ";" or ",")
[options]
-m or --mismatch take the number of mismatch allowed for each primer (default = 2)
-c or --contig : necessary if fasta files contains contigs
-r or --round : round MLVA score, take float, default is 0.25
-b or --binning : correct MLVA score for primers present in binning_file.csv 
(default binning_file contains brucella corrections, you can add exception by adding primers with insert_size and number of pattern repetitions in binning_file.csv. Binning file must be in Insilico directory)
--mixte : fasta file with one sequence will be considered as chromosome and fasta with sequences as contigs
(insert size for chromosome is calulated for circular chromosome) 
--full-locus-name : header will be full locus name instead of reduced locus name
--predicted-PCR-size-table : output a supplementary table with all predicted PCR size

primers format : primers name must be written as show below in primers_file :
<locus_name>_<pattern_size>bp_<theorical_insert>bp_<number_of_pattern_repetitions>U	forward_primer	reverse_primer

example : Lp03_96bp_941bp_8U	CAACCAATGAAGCAAAAGCA	AGGGGTTGATGGTCTCAATG

(mind that insert size contains both primers) 

### output ### 

MLVA_analysis_xxx.csv : csv file containing all MLVA value for each fasta file and all loci from the MLVA analysis 
	    this csv file is designed to be uploaded on http://microbesgenotyping.i2bc.paris-saclay.fr
output file : csv file containing All information from the MLVA analysis such as primers positions of match(s), size of insert, number of mismatch etc
mismatch file : txt file containing all different mismatch for each locus (only locus with mismatchs) found during MLVA analysis on input fasta sequence
predicted-pcr-size-table.csv : optional csv file containing predicted pcd size

### examples ###

python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/legionella_pneumophila/ -o . -p data_test/Legionella_pneumophila_primers.txt
python3.6 [/path/to]/MLVA/MLVA_finder.py --input data_test/legionella_pneumophila/ --output [/path/to]In_silico/ --primer data_test/Legionella_pneumophila_primers.txt (equivalent to the previous one)

with different number of mismatch allowed :
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/legionella_pneumophila/ -o . -p data_test/Legionella_pneumophila_primers.txt -m 1
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/legionella_pneumophila/ -o . -p data_test/Legionella_pneumophila_primers.txt -m 4

If fasta files contains contigs (and only contigs) :
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/legionella_pneumophila/ -o . -p data_test/Legionella_pneumophila_primers.txt -c
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/legionella_pneumophila/ -o . -p data_test/Legionella_pneumophila_primers.txt -c -m 0

#written by D.Christiany
