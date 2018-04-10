#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex, math, sys, os.path, getopt, csv, itertools

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}
#dico_degenerate = {'AG' : 'R','CT' : 'Y','CG' : 'S','AT':'W','GT':'K','AC':"M",'CGT' : 'B', "AGT" : "D",'ACT':'H',"ACG":"V","ACGT" : "N",".":"" }
dico_degenerate = {'R' : 'AG', 'Y' : 'CT', 'S' : 'CG', 'W':'AT','K':'GT','M':'AC','B':'CGT',"D":'AGT','H':'ACT',"V":"ACG","N":"ACGT",".":""}

def build_dictionnary(bin_file) : #build dictionnary for binning, take file from --binning  
	bin_file = open(bin_file,"r").read()
	bin_file = bin_file.replace("\t",";").replace(",",";").replace(" ",";").replace("\r","").split("\n")[:-1]
	bin_file=[primer.split(";") for primer in bin_file]

	dico_bin={}
	for primer in bin_file :
		if "-" in primer[1] :
			tmp= range(int(primer[1].split("-")[0]),int(primer[1].split("-")[1])+1)
			tmp=[(str(e),primer[2]) for e in tmp]

			if primer[0] in dico_bin :
				dico_bin[primer[0]].extend(tmp)
			else : 
				dico_bin[primer[0]]=tmp

		else :
			if primer[0] in dico_bin :
				dico_bin[primer[0]].append((primer[1],primer[2]))
			else : 
				dico_bin[primer[0]]=[(primer[1],primer[2])]

	return dico_bin

def clean_primers (primers_list) :
	tmp=[]
	for primers in primers_list :
		primers = primers.replace(" ",";").replace("\t",";").replace(",",";").split(";")
		primers=[primers[0]]+[primer.upper() for primer in primers[1:]]
		tmp.append(primers)
	return tmp

def inverComp (seq) :						#return the inversed complementary sequence : TTCGA -> TCGAA
	seq = seq.upper()                  		#acgt -> ACGT
	seq_comp="".join([dico_comp[nuc] for nuc in seq[::-1]])
	return (seq_comp)

#return all possible primers for a degenerated primer
def degenerated_primers (primer):
	degenerated_nucs=[]
	positions=[]
	primers=set()

	for i,nuc in enumerate(primer) :
		if nuc in dico_degenerate :
			positions.append(i)
			degenerated_nucs.append(dico_degenerate[nuc])
	combinations=list(itertools.product(*degenerated_nucs))
	primer=list(primer)

	for comb in combinations : 
		for pos, nuc in zip(positions,comb) :
			primer[pos]=nuc
		primers.add("".join(primer))

	return list(primers)


def binning_correction(primer,size,sizeU) : # correct sizeU if primer is in dico_bin
	if primer in dico_bin :
		tmp_size, tmp_sizeU = map(list,zip(*dico_bin[primer]))
		closest_size = min([int(e.split(' ')[0]) for e in tmp_size], key=lambda x:abs(x-size))
		pattern = int(primer.split('_')[1].replace('bp',''))
		if abs(closest_size - size) <= pattern : 
			sizeU = float(tmp_sizeU[tmp_size.index(str(closest_size))].replace('u',''))

	return sizeU

def get_flanking(seq,primers,pos1,pos2,splitted):
	if splitted is True :
		if pos1 < pos2 : 
			f1= inverComp(seq[pos2-flanking_len+len(primers[2]):pos2+len(primers[2])])
			f2= inverComp(seq[pos1:pos1+flanking_len])
		else :
			f1= seq[pos1-flanking_len+len(primers[1]):pos1+len(primers[1])]
			f2= seq[pos2:pos2+flanking_len]
	else :
		if pos1 < pos2 :
			f1= seq[pos1-flanking_len+len(primers[1]):pos1+len(primers[1])]
			f2= seq[pos2:pos2+flanking_len]
		else :
			f1=inverComp(seq[pos1:pos1+flanking_len])
			f2=inverComp(seq[pos2-flanking_len+len(primers[2]):pos2+len(primers[2])])
	return [f1,f2]

def positionsOfMatches (result,seq) : 				#get the matches positions in the fasta sequence
	pos = []  
	for res in result :
		pos.append([seq.find(res),res])
	return (pos)

def search_matches(nbmismatch,primer, seq) : 		#return all match(es) in the fasta sequence 
	return (regex.findall("("+primer+"){e<="+str(nbmismatch)+"}",seq,overlapped=True))

def remove_redondant_matches(match_list,primer) :
	reduced_list = [tmp for tmp in match_list if len(tmp)==len(primer)]
	if len(match_list) >0 and len(reduced_list) >0 :  
		return reduced_list
	else :
		return match_list

def pretty_mismatch (primer, found, nb_mismatch) : 	#lower the mismatches nucleotides 
	diff=abs(len(primer)-len(found))
	if len(found) == len(primer) and nb_mismatch>=len([e for e,i in zip(primer,found) if e!=i]) : #only mismatch(s)
		found="".join([n.lower() if n!=r else n for r,n in zip(primer,found)])

	elif len(found) < len(primer) : 				#deletion
		if found in primer :  						#only deletion, no mismatch(s), at the start or at the end
			diff = primer.find(found)		
			if diff == 0 :
				found = found+(len(primer)-len(found))*"."

			elif diff != -1 :
				found = (len(primer)-len(found))*"."+found
		elif len(primer)-len(found) == nb_mismatch : 	#only deletion, no mismatch(s)
			for i,n in enumerate(primer) :
				if i<len(found) :
					if found[i]!=n : found=found[:i]+"."+found[i:]
				else :
					found=found+'.'
		else : 											#deletion(s) and mismatch(s)
			i=0
			res=''
			tmp=''
			while found[i]==primer[i] :
				tmp+=found[i]
				i+=1 

			rev_found=list(reversed(found))
			rev_primer=list(reversed(primer))
			tmp2=''
			i=0
			while rev_found[i]==rev_primer[i] :
				tmp2+=rev_found[i]
				i+=1
			tmp2=tmp2[::-1]

			if (len(primer)-(len(tmp2)+len(tmp))) > 0 :
				res=tmp+"."*(len(primer)-(len(tmp2)+len(tmp)))+tmp2
			else : 
				exceed=abs(len(primer)-(len(tmp2)+len(tmp))) + (len(primer)-len(found))
				tmp2=tmp2[exceed:]
				res=tmp+"."*(len(primer)-(len(tmp2)+len(tmp)))+tmp2
			if len(res)==len(found) : found=res
			 
	elif len(found) > len(primer) :
		if diff == nb_mismatch : 				#only insertion(s)
			i=0
			tmp=primer
			while i < len(primer)+diff : 
				if found[i]!=tmp[i] : 
					found=found[:i]+found[i].lower()+found[i+1:]
					tmp=tmp[:i]+"."+tmp[i:]
				i+=1

		else : 									#insertion(s) + mismatch(s)
			res=''
			tmp = False
			for diff in range(nb_mismatch+1)[1:] :
				for i in range(len(found)) :
					if primer.find(found[:i])!=-1 and  primer.find(found[i+diff:])!=-1 :
						res = found[:i]+found[i:i+diff].lower()+found[i+diff:]
						tmp = True
						break 

				if tmp is True : break
			if len(res)==len(found) : found=res

	return found

#return the mismatches with differences in lower character from findfirst and findsec result
def clean_mismatches (nbprimer,primer,sense_list,found_list,nb_mismatch) : 

	res_found = []
	for sense,found in zip(sense_list,found_list) :
 
		if sense == "norm" :
			if nbprimer==2 : found = inverComp(found)
		else :
			if nbprimer==1 : found = inverComp(found)

		
		if found == primer : return ['']*len(found_list)
		res_found.append(pretty_mismatch(primer,found,nb_mismatch))

	return (res_found)

#function to find match(es) with at least two mismatches 
def mismatches (nb,primer,seq,nbmismatch) :
	match = []
	sense = []
	len_mismatch = []
	invP=inverComp(primer)	
	if nb==1 :
		tmpfind = search_matches(nbmismatch,primer,seq) 		#get all matches (with mismatches)
		tmpfind=remove_redondant_matches(tmpfind,primer)
		positions_f = positionsOfMatches(tmpfind,seq)			#get positions of matches
		if tmpfind == [] :
			tmpfind = search_matches(nbmismatch,invP,seq)		#get the match(es) with mismatches (use of regex.findall())
			tmpfind=remove_redondant_matches(tmpfind,primer)
			positions_r = positionsOfMatches(tmpfind,seq)		#get the position(s) of match(es)
			if tmpfind != [] :
				for i,res in enumerate(positions_r) :
					match.append(res[0])
					sense.append("inv")
					len_mismatch.append(len(tmpfind[i]))
		else :
			for i,res in enumerate(positions_f) :
				match.append(res[0])
				sense.append("norm")
				len_mismatch.append(len(tmpfind[i]))

	elif nb==2 :
		tmpfind = search_matches(nbmismatch,primer,seq) 		#get all matches
		tmpfind=remove_redondant_matches(tmpfind,primer)
		positions_f = positionsOfMatches(tmpfind,seq)			#get positions of matches
		if tmpfind != [] :
			for i,res in enumerate(positions_f) :
				match.append(res[0])
				sense.append("norm")
				len_mismatch.append(len(tmpfind[i]))


	return (match,sense,tmpfind,len_mismatch)					#return results (position of match + sense of primer)

#first search of the primer on the sequence (use inverComp() and mismatch())
def findFirst (primer,seq,nbmismatch) :

	match = []
	sense = []                    					#to store the sense of search (normal or inversed) 
	mismatchs = []
	len_mismatch=[]
	if nbmismatch == 0 :
		result=seq.find(primer)
		while (result!=-1) :          				#while the search has not been made on the entire sequence               
			match.append(result)             
			position=result+1     					#next search will start one nucleotide after the position of the last match
			sense.append("norm")
			result=seq.find(primer,position)		#next search, if no more match : resul = -1 -> end of the loop

		if match == [] :                         	#if no perfect match found with the regular primer 
			primer_inv=inverComp(primer)			#get the inversed complementary primer with inverComp()
			result=seq.find(primer_inv)
			while(result!=-1):						#same search with the converted primer
				match.append(result)                 
				position=result+1   
				sense.append("inv")
				result=seq.find(primer_inv,position) 

	elif nbmismatch >= 1 and match == [] :
		match,sense,mismatchs,len_mismatch = mismatches(1,primer,seq,nbmismatch)				

	return (match,sense,mismatchs,len_mismatch)

#search a match for the second primer
def findSec(primer,seq,sense,nbmismatch) :
	match = []
	mismatchs = []
	len_mismatch=[]
	if sense == "norm" : primer=inverComp(primer)
	if nbmismatch == 0 :
		result=seq.find(primer)
		while(result!=-1) :                      		#while there's a result
			match.append(result)          
			position=result+1							#get the position of the following nucleotide for the next search
			result=seq.find(primer,position)

	elif nbmismatch >= 1 and match == [] :
		match,trash,mismatchs,len_mismatch = mismatches(2,primer,seq,nbmismatch)

	return match,mismatchs, len_mismatch

#return the result of the matches
def find(primers,fasta,round,nbmismatch) : 
	
	fasta = fasta.replace(" ","").replace("\t","")		#delete spaces and tabulations                   
	
	sequences = fasta.split('>')[1:]					#split the fasta files into a list of fasta file
	dico_res = {}

	for s,seq in enumerate(sequences) :					#for each chromosome in the fasta file
		title_seq = seq.split("\n")[0]
		seq = "".join(seq.split("\n")[1:]).upper()

		if primers :									#if primers had been entered
			for p,primer in enumerate(primers) :		#for each couple of primers 
				primer_info = primer[0].split('_')
				primers_1 = degenerated_primers(primer[1])
				tmp1=[]
				tmp2=[]
				mismatchs=[]
				len_mismatch=[]
				for primer1 in primers_1 :
					tmp1_handle, tmp2_handle, mismatchs_handle, len_mismatch_handle = findFirst(primer1,seq,nbmismatch)
					tmp1.extend(tmp1_handle)
					tmp2.extend(tmp2_handle)
					mismatchs.extend(mismatchs_handle)
					len_mismatch.extend(len_mismatch_handle)
				#print (tmp1)
				if nbmismatch > 0 and mismatchs != [] : mismatchs = clean_mismatches(1,primer[1],tmp2,mismatchs,nbmismatch)
				first_match = tmp1, tmp2, len_mismatch								#search match(es) for the first primer 
				result = []
				insert=""
				for i,pos_match in enumerate(first_match[0]) :						#for each match of the first primer
					primers_2=degenerated_primers(primer[2])
					tmp=[]
					mismatchs2=[]
					len_mismatch2=[]
					for primer2 in primers_2 :
						tmp_handle, mismatchs2_handle, len_mismatch2_handle = findSec(primer2,seq,first_match[1][i],nbmismatch)
						tmp.extend(tmp_handle)
						mismatchs2.extend(mismatchs2_handle)
						len_mismatch2.extend(len_mismatch2_handle)													#search match(es) for the second primer
					if nbmismatch > 0 and mismatchs2 != [] : mismatchs2 = clean_mismatches(2,primer[2],[first_match[1][i]]*len(tmp),mismatchs2,nbmismatch) 	#lower the mismatched nucleotides
					if tmp != [] :													#if there is a match with the second primer on the complementary DNA sequence
						for m,pos_match2 in enumerate(tmp) :						#for each match found for the second primer
							if nbmismatch >0 : 
								len_match=first_match[2][i]
								len_match2=len_mismatch2[m]
								mismatch=mismatchs[i]								#get the coresponding mismatch of the 1st primer
								mismatch2=mismatchs2[m]								#get the coresponding mismatch of the second primer 
							else :
								len_match=len(primer[1])
								len_match2=len(primer[2])
								mismatch=''
								mismatch2=''

							#size is calculated with the primers size given, and not the match of the primers which may contain indel
							splitted = False
							if first_match[1][i] == "inv" :
								size = pos_match+len(primer[1])-(pos_match2+len_match2-len(primer[2]))
								size2 = pos_match+len(primer[1])+(len(seq)-(pos_match2+len_match2-len(primer[2])))	#if primers are separated by the splitted area in the sequence (circular chromosome)
							else :
								size = pos_match2+len(primer[2])-(pos_match+len_match-len(primer[1])) 
								size2 =pos_match2+len(primer[2])+(len(seq)-(pos_match+len_match-len(primer[1])))

							if contig is False and size > 0 and size2 < size : 
								size = size2	
								splitted=True

							if size > 0 : 							

								#insert may contains indel in primers area
								if splitted is True :									#if insert is separated by the splitted area in the sequence (circular chromosome)
									if pos_match < pos_match2 : insert = inverComp(seq[pos_match2:]+seq[:pos_match+len(primer[1])]) #reversed comp insert
									else : insert = seq[pos_match:]+seq[:pos_match2+len(primer[2])] 
								else :		
									if pos_match < pos_match2 : insert = seq[pos_match:pos_match2+len(primer[2])]
									else : insert = inverComp(seq[pos_match2:pos_match+len(primer[1])]) #reversed comp insert		

								sizeU = abs(float(primer_info[3].upper().replace("U",""))-\
									((float(primer_info[2].lower().replace("bp",""))-size)\
									/float(primer_info[1].lower().replace("bp",""))))					#computation of sizeU

								if binning is True and primer[0]in dico_bin : 					#if option binning selected, correction with the dictionary
									sizeU = binning_correction(primer[0],size,sizeU)

								if sizeU < 100 : result.append([primer[0],pos_match,pos_match2,size,sizeU,sequence+str(s+1),nbmismatch,primer[1],mismatch,primer[2],mismatch2,insert])

				if len(result) == 0 and primer_info[0] not in dico_res :	#if no result
					dico_res[primer_info[0]]=[primer[0],"","","","","",nbmismatch,primer[1],"",primer[2],"",""]	

				elif len(result) > 0 :											#if result(s)
					best_res = result[0]
					for res in result :											#keep the result with the minimum sizeU value
						if res[4]<best_res[4] : best_res=res	

					if round !="" and round>0 :                    				#round of the sizeU value
						sizeU=best_res[4]
						if sizeU>=math.floor(sizeU) and sizeU<(math.floor(sizeU)+round) :
							sizeU = math.floor(sizeU)
						elif sizeU <= math.ceil(sizeU) and sizeU>(math.ceil(sizeU)-round) :
							sizeU=math.ceil(sizeU)
						else :
							sizeU=math.floor(sizeU)+0.5
						if str(sizeU)[-2:]=='.0' : sizeU=int(sizeU)
						best_res[4]=sizeU									#set of the rounded sizeU value
					if primer_info[0] in dico_res and dico_res[primer_info[0]][4] != "" :
						best_res[5] = best_res[5]+", "+sequence+str(s+1) 	#if there's already a result with perfect matches
					dico_res[primer_info[0]]=best_res						#set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
					if flanking is True : 
						dico_flanking[primer_info[0]]=get_flanking(seq,primer,pos_match,pos_match2,splitted)

	return dico_res

#return primers with no result 
def get_empty_locus (dico_result) :
	tmpprimers = []
	for locus in dico_result :
		if dico_result[locus][4] == '' :
			tmpprimers.append([dico_result[locus][0],dico_result[locus][7],dico_result[locus][9]])
	return tmpprimers

#search MLVA with perfect match, then allow one additional mismatch for locus with no result until the number max of mismatch allowed 
def run (Primers,fasta,round,nbmismatch) : 
	tmp = len(Primers)
	tmpPrimers = Primers
	result = {}
	for mismatch_allowed in range(int(nbmismatch)+1) :
		tmp_dico = find(tmpPrimers,fasta,round,mismatch_allowed) 	#search with no mismacth
		result.update(tmp_dico)										#add results to the dictionnary
		tmpPrimers = get_empty_locus(result)						#only keep locus with no result
		nb_match = tmp -len(tmpPrimers) 
		print ("results with",mismatch_allowed,"mismatch(s): ",nb_match,"/",len(Primers))

		tmp = len(tmpPrimers)

		if len(tmpPrimers) == 0 :
			break 

	if len(tmpPrimers) != 0 : 
		print ("no match : ", tmp)

	return result

def usage() : #example of command to use insilico.py 
	print ("./insilico.py -i <input_directory> -o <output_directory> -p <primers_file> \n \
	[option -c for contigs] \n\
	[option -m x for number of mismatch, default = 2] \n\
	[option -r x : round MLVA score, default = 0.25] \n\
	[option -b : binning file is used to correct MLVA value for primers in binning_file.csv] \n\
	[option --mixte : fasta file with one sequence will be considered as chromosome and fasta with sequences as contigs] \n\
	[option --full-locus-name : header will be full locus name instead of reduced locus name] \n\
	[option --predicted-PCR-size-table : output a supplementary table with all predicted PCR size ] \n\
	[option --flanking-seq <int>: add flanking column in <output.csv>, flanking are the sequences before and after the insert (primers inculded), you can chose the size of flanking sequences <int> ")
  
def main() : #run find() for each genome file in the directory with all primers in the primers file

	if len(sys.argv)<2 : 
		usage()
		sys.exit(2)
	
	try:		#check if correct args
		opts, args = getopt.getopt(sys.argv[1:], "hm:i:o:p:cr:b:f:", ["help", "mismatch=", "input=", "output=", "primer=", "contig","round="\
			,"binning=","mixte","full-locus-name","predicted-PCR-size-table","flanking-seq="])
	except getopt.GetoptError as err:
		usage()
		sys.exit(2)
	nb_mismatch = 2		#default value for the number of mismatch allowed 
	global sequence
	sequence = "seq"
	global contig 
	contig = False
	global binning 
	binning = False
	global mixte
	mixte = False
	global full_locus_name
	full_locus_name = False
	global predicted_PCR_size_table
	predicted_PCR_size_table = False
	global flanking
	flanking = False
	global dico_flanking
	dico_flanking={} 
	round = 0.25
	for opt, arg in opts: #get args given by user
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			if os.path.exists(arg) is False : 
				print ("input directory path is invalid") 
				usage()
				sys.exit(2)
			fasta_path = arg
			if fasta_path[-1]!="/" : fasta_path=fasta_path+"/"
			files = os.listdir(fasta_path)
			files=sorted(files)
		elif opt in ("-p", "--primer"): 
			if os.path.exists(arg) is False : 
				print ("primers list csv path is invalid") 
				usage()
				sys.exit(2)
			Primers = open(arg,"r").read()
		elif opt in ("-o", "--output"):
			if os.path.exists(arg) is False : os.makedirs(arg)
			output_path = arg
			if output_path[-1] != "/" : output_path=output_path+"/"
		elif opt in ("-m","--mismatch"):
			nb_mismatch = int(arg)
		elif opt in ("-c", "--contig"):
			contig = True
			sequence = "contig"
		elif opt in ("-r", "--round"):
			round = float(arg)
		elif opt in ("-b", "--binning"):
			binning = True
			global dico_bin
			dico_bin = build_dictionnary(arg)
		elif opt in "--mixte" :
			mixte = True
		elif opt in "--full-locus-name" :
			full_locus_name = True
		elif opt in "--predicted-PCR-size-table" :
			predicted_PCR_size_table = True
		elif opt in ("--flanking-seq" ) :
			global flanking_len
			flanking_len=int(arg)
			flanking=True
			#print (flanking,flanking_len)
		else:
			assert False, "unhandled option"

	Primers = Primers.replace("\r","").split("\n")
	if Primers[-1]=="" : del Primers[-1]
	Primers=clean_primers(Primers)
	Primers_short = [pri[0].split("_")[0] for pri in Primers]
	Primers_full = [pri[0] for pri in Primers]

	#header of output_file
	for i,file in enumerate(files) :
		print (file, "\t strain ",i+1,"/",len(files))
		fasta = open(fasta_path+file,"r").read()
		fasta_names = []
		for line in fasta.split("\n") :
			if ">" in line :
				if '|' in line and len(line.split("|")) >=5 : 	#for fasta names like : >gi|1032812322|ref|NZ_CP015344.1| Legionella pneumophila strain D-7630, complete genome
					tmp_names=line.replace("\n","").split("|")
					fasta_names.append([tmp_names[3].split(".")[0],tmp_names[4][1:].split(",")[0]]) 
					del tmp_names
				else : 											#for fasta names like : >NC_003317.1 Brucella melitensis bv. 1 str. 16M chromosome I, complete sequence
					tmp_names=line.replace("\n","").split(",")[0].split(" ")
					fasta_names.append([tmp_names[0].replace(">","")," ".join(tmp_names[1:])])
		if len(fasta_names)==0 : 
				print ('no fasta header found in',file) 
				sys.exit()

		if mixte is True  :	#check for each fasta if there is more than one sequence (contigs) or not with --mixte option
			if fasta.count(">") > 1 : 
				contig = True
				sequence = "contig"
			else : 
				contig = False 
				sequence = "chr"

		result = run(Primers,fasta,round,nb_mismatch)			#use find for each number of mismatch

		#output file 
		if i ==0 : 
			output_file=[]
			header = ["strain","primer","position1","position2","size","allele","sequence","nb_mismatch","primer1","mismatch1","primer2","mismatch2","predicted PCR target"]
			if flanking is True : header.extend(['flanking1',"flanking2"])
			out = csv.writer(open(output_path+fasta_path.split("/")[-2]+"_output.csv","w",encoding='utf-8'), delimiter=';',quoting=csv.QUOTE_NONE)
			for row in [header] :
				out.writerow(row)
			out = csv.writer(open(output_path+fasta_path.split("/")[-2]+"_output.csv","a",encoding='utf-8'), delimiter=';',quoting=csv.QUOTE_NONE)

		cr=[]
		for Primer in Primers_short :
			if result[Primer][4]=='' : result[Primer][6]='ND'
			if flanking is True : 
				if Primer in dico_flanking : cr.append([file]+result[Primer]+dico_flanking[Primer])
				else : cr.append([file]+result[Primer]+["",""])
			else : cr.append([file]+result[Primer])

		for row in cr :
			out.writerow(row)

		output_file.extend(cr)

		#MLVA_analysis_file
		locus, mlva_score, mlva_insert = ([] for i in range(3))
		for primer_short,primer_full in zip(Primers_short,Primers_full) :
			if full_locus_name is False : locus.append(primer_short)
			else : locus.append(primer_full)
			mlva_score.append(str(result[primer_short][4]))		#scores 
			if predicted_PCR_size_table is True : mlva_insert.append(str(result[primer_short][3]))

		if i==0 :
			pathfile = output_path+"MLVA_analysis_"+fasta_path.split("/")[-2]+".csv"
			output = open(pathfile,"w",encoding='utf-8') 									#output is a csv file (delimiter=";")
			output.write(";".join(["key","Access_number"]+locus)+"\n")  	#header
			output = open(pathfile,"a",encoding='utf-8')

		output.write(";".join([str(i+1).zfill(3),file.split(".")[0]]+mlva_score)+"\n")

		#predicted PCR table
		if predicted_PCR_size_table is True :

			if i==0 :
				pathfile = output_path+"predicted_PCR_size_table_"+fasta_path.split("/")[-2]+".csv"
				output_pcr_size = open(pathfile,"w",encoding='utf-8') 									#output is a csv file (delimiter=";")
				output_pcr_size.write(";".join(["key","Access_number"]+locus)+"\n")  	#header
				output_pcr_size = open(pathfile,"a",encoding='utf-8')

			output_pcr_size.write(";".join([str(i+1).zfill(3),file.split(".")[0]]+mlva_insert)+"\n")

	output.close()
	if predicted_PCR_size_table is True : output_pcr_size.close()
	
	print ("MLVA analysis finished for "+fasta_path.split("/")[-2])

	##### creation of mismatch summary txt file #####

	dico_mismatch = {}
	for primer in Primers :
		dico_mismatch[primer[0]+"_FOR"] = set([])
		dico_mismatch[primer[0]+"_REV"] = set([])

	for line in output_file[1:] : #whithout the header
		if line[9] != "" : dico_mismatch[line[1]+"_FOR"].add(line[9])
		if line[11] != "" : dico_mismatch[line[1]+"_REV"].add(line[11])

	dico_mismatch = {key: value for key, value in dico_mismatch.items() if value != set() } 	#delete keys without value(s)
	
	output_mismatch = open(output_path+fasta_path.split("/")[-2]+"_mismatchs.txt","w")
	tmp_file =""
	for primer in Primers :
		if primer[0]+"_FOR" in dico_mismatch :
			tmp_file += primer[0]+"_FOR\r\n"+primer[1]+"\r\n"+"\r\n".join(list(dico_mismatch[primer[0]+"_FOR"]))+"\r\n\r\n"
		if primer[0]+"_REV" in dico_mismatch :
			tmp_file += primer[0]+"_REV\r\n"+primer[2]+"\r\n"+"\r\n".join(list(dico_mismatch[primer[0]+"_REV"]))+"\r\n\r\n"
	output_mismatch.write(tmp_file)
	output_mismatch.close()

if __name__ == "__main__" :
	main()