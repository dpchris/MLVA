MLVA_finder.py README

### Dependencies ###

python3
python-dev (sudo apt-get install python-dev), 
regex (https://launchpad.net/ubuntu/+archive/primary/+files/python-regex_0.1.20170117.orig.tar.gz)

### goal ###

MLVA_finder.py is a python script designed to do Multi loci VNTR analysis
(VNTR stands for Variable Number of Tandem Repeats ).
MLVA_finder.py performs an in silico PCR to extract sequences of tandem repeat from submitted fasta file(s)  
and call VNTR alleles.
This script use a list of primers to recover sequences from the VNTRs.
The number of allowed mismatches can be set in the command line. 

### input ###

-i or --input takes the path of the directory containing all fasta file which will be used by the script 
-o or --output takes the path of the results directory where results will be saved
-p or --primer takes the path of the primers list in a csv file (delimited by '\t', ";" or ",")
[options]
-m or --mismatch takes the number of mismatches allowed for each primer (default = 2)
-c or --contig : necessary if fasta files contain contigs
-r or --round : rounds MLVA allele values, take float, default is 0.25 (meaning that an allele value of an integer +- 0.25 will be rounded to this integer; larger deviations from this value will be called as "intermediate size alleles" i.e. all alleles between 1.26 and 1.74 will be called 1.5)
-b or --binning : corrects MLVA allele calls for primers present in binning_file.csv 
(default binning_file contains Brucella MLVA assay corrections as indicated in published allele numbering system version 3.6 (http://mlva.u-psud.fr/brucella/spip.php?article93) 
Rules are defined by adding primers with insert_size and number of pattern repetitions in binning_file.csv. Binning file must be located in the same directory as MLVA_finder.py)
--mixte : a fasta file with a single sequence will be considered as "chromosome" and fasta files with multiple sequences as contigs
(insert size for "chromosome" is calculated for circular chromosome) 
--full-locus-name : header will be full locus name instead of reduced locus name
--predicted-PCR-size-table : output a supplementary table with all predicted PCR sizes
--flanking-seq <int>: add flanking column in <output.csv>, flanking are the sequences before and after the insert (primers inculded), you can chose the size of flanking sequences <int>

primers format : primers name must be written as shown below in primers_file :
<locus_name>_<pattern_size>bp_<insert_size_in_reference_genome>bp_<corresponding allele coding convention>U	forward_primer	reverse_primer

example : Lp03_96bp_941bp_8U	CAACCAATGAAGCAAAAGCA	AGGGGTTGATGGTCTCAATG

(mind that indicated insert size contains both primers) 

### output ### 

MLVA_analysis_xxx.csv : csv file containing all MLVA values for each fasta file and all loci from the MLVA analysis 
	    this csv file is designed to be uploaded on http://microbesgenotyping.i2bc.paris-saclay.fr
output file : csv file containing All informations from the MLVA analysis such as primers positions of match(s), size of insert, number of mismatches etc
mismatch file : txt file containing all different mismatches for each locus (only locus with mismatches) found during MLVA analysis on input fasta sequence
predicted-pcr-size-table.csv : optional csv file containing predicted PCR size

### command line examples ###

python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Legionella_pneumophila/ -o . -p data_test/primers/Legionella_pneumophila_primers.txt
python3.6 [/path/to]/MLVA/MLVA_finder.py --input data_test/assemblies/Legionella_pneumophila/ --output [/path/to]/results/ --primer data_test/primers/Legionella_pneumophila_primers.txt (equivalent to the previous one)

#with different number of mismatch allowed :
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Legionella_pneumophila/ -o . -p data_test/primers/Legionella_pneumophila_primers.txt -m 1	#1 mismatch allowed
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Legionella_pneumophila/ -o . -p data_test/primers/Legionella_pneumophila_primers.txt -m 4	#4 mismatches allowed

#If fasta files contains contigs (and only contigs) :
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Legionella_pneumophila/ -o . -p data_test/primers/Legionella_pneumophila_primers.txt -c
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Legionella_pneumophila/ -o . -p data_test/primers/Legionella_pneumophila_primers.txt -c -m 0	#0 mismatches allowed

#Using binning option :
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Brucella -p data_test/primers/Brucella_primers.txt -o . -b data_test/Brucella_binning_file.csv 
python3.6 [/path/to]/MLVA/MLVA_finder.py -i data_test/assemblies/Brucella -p data_test/primers/Brucella_primers.txt -o . -b data_test/Brucella_binning_file.csv --predicted-PCR-size-table --flanking-seq 200 


#written by D.Christiany
