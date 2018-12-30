import requests, sys
import os
from Bio import SeqIO
import re


'''
def chromQuery(query):
	"get the chrom position from ensembl for each zebrafish protein variant"
	server = "https://rest.ensembl.org"
	#ext = "/map/translation/ENSDARP00000123118/60..60?"
	ext = "/map/translation/"+query+"?"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	 
	decoded = r.json() #type dictionary

	data = decoded["mappings"]
	chrom_name = data[0]["seq_region_name"]
	chrom_start = data[0]["start"]
	chrom_end = data[0]["end"]
	position = str(chrom_name)+"\t"+str(chrom_start)+"\t"+str(chrom_end) #get fasta non inclusisve
	position_flank20 = str(chrom_name)+"\t"+str(chrom_start-20)+"\t"+str(chrom_end+20)
	position_flank150 = str(chrom_name)+"\t"+str(chrom_start-150)+"\t"+str(chrom_end+150)
	return (position,position_flank20,position_flank150)
	#get three positions for better analysis of the sequence information


gene_dict = {}
f1 = open("zebrafish_protein.txt","r")
for line in f1.readlines()[1:]:
	gene = line.strip().split("\t")[2]
	proteinID = line.strip().split("\t")[1]
	gene_dict[gene] = proteinID
f1.close()

#head human_variant_VEP_extract_info.txt
#16:70268356-70268356	AARS	ENSG00000090861	ENST00000261772	986	329	R/H	cGc/cAc
#16:70269756-70269756	AARS	ENSG00000090861	ENST00000261772	824	275	G/D	gGt/gAt
#2:214982225-214982225	ABCA12	ENSG00000144452	ENST00000272895	4541	1514	R/H	cGc/cAc
#1:94000866-94000866	ABCA4	ENSG00000198691	ENST00000370225	6449	2150	C/Y	tGt/tAt
#1:94000878-94000878	ABCA4	ENSG00000198691	ENST00000370225	6437	2146	G/D	gGc/gAc


#head human_zebrafish_conserved_variant.txt 
#ACTN1	849	actn1	859	R	R
#NEXMIF	443	nexmifa	392	S	S
#NEXMIF	385	nexmifa	347	G	G

fa = open("human_zebrafish_conserved_variant.txt","r")

fc = open("human_zebrafish_conserved_variant_modified.txt","w+")

for line1 in fa.readlines():
	line1 = line1.strip().split("\t")
	line1_query = line1[0]+line1[1]
	fb = open("human_variant_VEP_extract_info.txt","r")
	for line2 in fb.readlines():
		line2 = line2.strip().split("\t")
		line2_query = line2[1]+line2[5]
		if line2_query == line1_query:
			fc.write(line1[0]+"\t"+line1[1]+"\t"+line1[2]+"\t"+line1[3]+"\t"+line2[6]+"\n")
	fb.close()
fa.close()
fc.close()			

#head human_zebrafish_conserved_variant_modified.txt 
#ACTN1	849	actn1	859	R/H
#NEXMIF	443	nexmifa	392	S/N
#NEXMIF	385	nexmifa	347	G/D
#B4GALT7	214	b4galt7	215	C/Y


f2 = open("human_zebrafish_conserved_variant.txt","r")
f3 = open("zebrafish_variant_chrom_position.bed","w+")
#Note: getfasta didn't not get the first nulceotide for each codon, but it doesn't matter, since it's G in the middle and the codon codes for arginine
#>cyp26b1:396
#GC
#>cyp26b1:362
#GT

f4 = open("zebrafish_variant_chrom_flank20.bed","w+")
f5 = open("zebrafish_variant_chrom_flank150.bed","w+")
for line in f2.readlines():
	#print line
	gene = line.strip().split("\t")[2]
	protein_pos = line.strip().split("\t")[3]
	if gene in gene_dict.keys():
		query = gene_dict[gene]+"/"+str(protein_pos)+".."+str(protein_pos)
		position = chromQuery(query)
		f3.write(position[0]+"\t"+gene+":"+protein_pos+"\n")
		f4.write(position[1]+"\t"+gene+":"+protein_pos+"\n")
		f5.write(position[2]+"\t"+gene+":"+protein_pos+"\n")
f2.close()
f3.close()
f4.close()
f5.close()



#os.system('sed -i.bak '/CHR_ALT_CTG9/d' zebrafish_variant_chrom_position.bed')
#zebrafish_chrom3.fa without modified
#7	25346667	25346669	cyp26b1:396
#7	25344786	25344788	cyp26b1:362
#23	44478150	44478152	actl6b:248
#5	23126324	23126326	nexmifa:230
#13	36330591	36330593	actn1:762

#os.system(' awk '$2=$2-1 {print $1"\t"$2"\t"$3"\t"$4}' zebrafish_variant_chrom_position.bed >zebrafish_variant_chrom_position_3.bed
')




os.system('bedtools getfasta -fi zv11.fa -bed zebrafish_variant_chrom_position_3.bed -name -fo zebrafish_chrom3.fa')
'''
#head zebrafish_chrom3.fa 
#>actn1:859
#CCT
#>nexmifa:392
#AGC
#>nexmifa:347
#GGC
#>b4galt7:215
#GCA
#>b4galt7:273
#AGC

#head human_zebrafish_conserved_variant_modified.txt 
#ACTN1	849	actn1	859	R/H
#NEXMIF	443	nexmifa	392	S/N
#NEXMIF	385	nexmifa	347	G/D
#B4GALT7	214	b4galt7	215	C/Y
'''
GC_seq_left = []
GC_seq_right = []
i = 0
for seq_record in SeqIO.parse("zebrafish_chrom3.fa","fasta"):
	header = seq_record.id.strip()
	seq = str(seq_record.seq)
	fd = open("human_zebrafish_conserved_variant_modified.txt","r")
	for line in fd.readlines():
		line = line.strip().split("\t")
		if line[2]+":"+line[3] == header and "GC" in seq: 
			i=i+1
			if line[4]=="A/V":
				GC_seq_right.append(header)
			else: 
				GC_seq_left.append(header)
print len(GC_seq_left) 	#1363
print len(GC_seq_right)	#0
print i #1363
# no A/V variant in human_zebrafish_conserved_variant_modified.txt actually
'''

GC_seq = [] #G to A couple with CCN PAM, #C to T couple with NGG PAM
for seq_record in SeqIO.parse("zebrafish_chrom3.fa","fasta"):
	seq = str(seq_record.seq)
	if "GC" in seq:
		GC_seq.append(seq_record.id)

#os.system('bedtools getfasta -fi zv11.fa -bed zebrafish_variant_chrom_flank20.bed -name -fo zebrafish_chrom_20.fa')

f8 = open("zebrafish_variant_for_further_analysis1.fa","w+")
for seq_record in SeqIO.parse("zebrafish_chrom_20.fa","fasta"):
	header = seq_record.id
	if header in GC_seq:	
		pattern1 = re.compile("CC")  #CCN PAM
		seq = str(seq_record.seq)
		if pattern1.search(seq[:17]): #lowercase for PAM and mutation C to T site
			seq = re.sub('CC','cc',seq[:17])+seq[17:19]+re.sub('GC','gC',seq[19:22])+seq[22:]
			f8.write(">"+header+"\n"+seq+"\n")
		#pattern2 = re.compile("GG")
		#if pattern2.search(seq[24:]): #NGG PAM
			#seq = seq[:19]+seq[19:22].lower()+seq[22:24]+re.sub('GG','gg',seq[24:])
			#f8.write(">"+header+"\n"+seq+"\n")
f8.close()


#fish_human_final_candidate.txt
#wdr34	ENSDART00000080347	WDR34	ENST00000372715
#gria3a	ENSDART00000137120	GRIA3	ENST00000616590

#human_zebrafish_conserved_variant.txt 
#CYP26B1	397	cyp26b1	396	R	R
#CYP26B1	363	cyp26b1	362	R	R

#VEFhumanVariantGtoA.txt 
#16:70268356-70268356	AARS	ENSG00000090861	ENST00000261772	986	329	R/H	cGc/cAc
#9:104791982-104791982	ABCA1	ENSG00000165029	ENST00000374736	5774	1925	R/Q	cGg/cAg



f11 = open("all_variants_again1.txt","w+") #final dataset
for seq_record in SeqIO.parse("zebrafish_variant_for_further_analysis1.fa","fasta"):
	header = seq_record.id
	fishgene = header.split(":")[0]
	fishpos = header.split(":")[1]
	seq = str(seq_record.seq)
	f10 = open("fish_human_final_candidate.txt","r")
	for line in f10.readlines():
		line = line.strip().split("\t")
		if fishgene == line[0]:
			f11.write(line[2]+"\t"+line[3]+"\t"+line[0]+"\t"+line[1]+"\t"+fishpos+"\t"+seq+"\n")
			f10.close()
f11.close()

#head all_variants_again.txt 
#NEXMIF	ENST00000055682	nexmifa	ENSDART00000149893	392	TGGTGAATTTAGTGATGACagcTCTTGCACggGCTCCCCAGA
#NEXMIF	ENST00000055682	nexmifa	ENSDART00000149893	347	AccATGCAGTAGAAAGGAGggcCCAAAAGAGAAACCAGACCA


#head human_zebrafish_conserved_variant_modified.txt 
#ACTN1	849	actn1	859	R/H
#NEXMIF	443	nexmifa	392	S/N
#NEXMIF	385	nexmifa	347	G/D
#B4GALT7	214	b4galt7	215	C/Y

f12 = open("all_variants_again1.txt","r")
f14 = open("all_variants_again_modified.txt","w+")
for line in f12.readlines():
	line = line.strip().split("\t")
	f13 = open("human_zebrafish_conserved_variant_modified.txt","r")
	for line1 in f13.readlines():
		line1 = line1.strip().split("\t")
		if line[0]+line[4] == line1[0]+line1[3]:
			f14.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line1[4]+"\t"+line[5]+"\n")
	f13.close()
f12.close()
f14.close()

#head all_variants_again_modified.txt 
#NEXMIF	ENST00000055682	nexmifa	ENSDART00000149893	347	G/D	AccATGCAGTAGAAAGGAGGgCCCAAAAGAGAAACCAGACCA
#B4GALT7	ENST00000029410	b4galt7	ENSDART00000170614	215	C/Y	AGCGGTTTGACATccCGTTgCACTTTAGGCAGAAATCACAAG
#B4GALT7	ENST00000029410	b4galt7	ENSDART00000170614	273	A/T	ATATTAccTGTTTTTGTGCAgCGATCCGCTTCTGGTCTCTTT
os.system('mv all_variants_again_modified.txt all_variants_negative.txt')


#After processing
#sort all_variants_again.txt > all_variants_again_sorted.txt
#sort RtoH.txt > RtoH_sorted.txt
#cut -f 1-5 all_variants_again_sorted.txt |sort > all_intermediate.txt
#cut -f 1-5 RtoH_sorted.txt |sort > RtoH_intermediate.txt
#comm -13 RtoH_intermediate.txt all_intermediate.txt > all_intermediate2.txt (1597 records)	
#cat all_intermediate2.txt | xargs -I {} grep {} all_variants_again_sorted.txt | sort | uniq > all_variants_final.txt (1133 unique sequence, 1626 records)


