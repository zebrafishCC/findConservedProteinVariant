from Bio import SeqIO
from Bio.Blast import NCBIXML

def blast(human_fasta,fish_fasta,human_protein_pos):			
	"given one human fasta file and zebrafish fasta file, find conserved human transcript and try to record the relevant zebrafish information"
	'''	
	conserved_file = open("human_zebrafish_conserved_variant.txt","w+")
	for seq_record in SeqIO.parse(human_fasta,"fasta"):
        	human_header = seq_record.id	
		human_gene = human_header.strip().split("|")[2]
	for seq_record in SeqIO.parse(fish_fasta,"fasta"):
		fish_header = seq_record.id
                fish_gene = fish_header.strip().split("|")[2]

	'''	
	#blast results:
	from Bio.Blast.Applications import NcbiblastpCommandline

	blastp_cline = NcbiblastpCommandline(query=human_fasta,subject=fish_fasta,outfmt=5,out="test.xml")

	#print blastp_cline
	blastp_cline() # run blastp

	result_handle = open("test.xml")

	#parse blast results
	blast_record = NCBIXML.read(result_handle)

	#function dealing with protein positions

	human_protein_pos = int(human_protein_pos)

	#login: make sure the human mutation site is conserved in fish protein level
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			#print hsp.query
			#print hsp.query_start
			i = 0 
			j = 0 #space in query sequence
			
			if human_protein_pos >= hsp.query_start: #make sure the alignment start position is on the left of human_protein_pos
				space = human_protein_pos-hsp.query_start+1
				for character in hsp.query:					
					if i < space:
						if character == "-": #calculate the occurrence of non amino acids
							j = j+1 
						elif character.isalpha():
							i = i+1	
				index = i+j
				if index < int(hsp.align_length): 
					#print hsp.align_length
					if hsp.match[index-1].isalpha():
						ii = 0
						k = 0 #space in subject sequence
						for character in hsp.sbjct:
							if ii+k < index:
								if character == "-":
									k = k+1
								else:
									ii=ii+1
						subject_pos = hsp.sbjct_start+ii-1 #need to adjust
						return (subject_pos,hsp.match[index-1])
				else:
					return None


#blast("humanAAAS.fasta","zebrafishaaas.fasta",329)

###from Bio import SeqIO

#build a dictionary for two sequence alignment
f1 = open("fish_human_final_candidate_positive.txt","r")
gene_dict = {}
for line in f1.readlines():
        line = line.strip().split("\t")
        zebrafish_transcript = line[1]
        human_transcript = line[3]
        gene_dict[human_transcript] = zebrafish_transcript
#print len(gene_dict)
f1.close()

pos_dict = {}
#f2 = open("VEFhumanVariant.txt","r")
#for line in f1.readlines():
for seq_record in SeqIO.parse("humanProtein_positive.fa","fasta"):
        human_header = seq_record.id
        #human_transcript = line.strip().split("\t")[4]
        human_transcript = human_header.strip().split("|")[1]
        pos_dict[human_transcript] = []
        f2 = open("human_variant_VEP_extract_info_positive.txt","r")
        for line in f2.readlines():  #readlines method can only be used once?
                transcript = line.strip().split("\t")[3]
                #print transcript
                if transcript == human_transcript:
                        human_protein_pos = line.strip().split("\t")[5]
                        #print transcript
                        #print human_protein_pos
                        pos_dict[human_transcript].append(human_protein_pos)
                else:
                        continue
        f2.close()

'''
for key in pos_dict.keys():
        if len(pos_dict[key])==0:
                del pos_dict[key]
'''
#make sure the GC pattern in the human transcript, which leads to only R-->Q and R-->H change, V-->A change
humanGCtranscript = []
f5 = open("human_selected_transcripts_positive.txt","r")
for line in f5.readlines():
	transcript = line.strip()
	humanGCtranscript.append(transcript)
f5.close()
#print len(pos_dict)
#print pos_dict.values()  
#pos_dict results shows that some human transcript for some gene didn't have variant. This problem arise because only choose transcript for one gene for simplicity.
###


#f3 = open("human_intermediate.fasta","w")
#f4 = open("fish_intermediate.fasta","w")
#main blast function

conserved_file = open("human_zebrafish_conserved_variant_positive.txt","w+")
for seq_record in SeqIO.parse("humanProtein_positive.fa","fasta"):
	f3 = open("human_intermediate.fasta","w")
	human_header = seq_record.id
	#print human_header
	human_transcript = human_header.strip().split("|")[1]
	if human_transcript in humanGCtranscript:
		human_gene = human_header.strip().split("|")[2]
		human_seq = str(seq_record.seq)
		if human_seq.isupper():
			f3.write(">"+human_header+"\n"+human_seq)
			f3.close()
			for seq_record in SeqIO.parse("zebrafishProtein_positive.fa","fasta"):
				f4 = open("fish_intermediate.fasta","w")
				zebrafish_header = seq_record.id
				#print zebrafish_header
				fish_gene = zebrafish_header.strip().split("|")[2]
				zebrafish_transcript = zebrafish_header.strip().split("|")[1]
				zebrafish_seq = str(seq_record.seq)
				#if zebrafish_seq != "Sequenceunavailable":
				if zebrafish_seq.isupper():
					f4.write(">"+zebrafish_header+"\n"+zebrafish_seq)
					f4.close()
					if gene_dict[human_transcript] == zebrafish_transcript:
						#print human_transcript+"\t"+zebrafish_transcript
						#print human_seq
						#print zebrafish_seq
						if len(pos_dict[human_transcript]) > 0:
							for human_protein_pos in pos_dict[human_transcript]:
								blast_result = blast("human_intermediate.fasta","fish_intermediate.fasta",human_protein_pos)
								if blast_result:
									subject_pos = blast_result[0]
									match_amino_acid = blast_result[1]
									human_amino_acid = human_seq[int(human_protein_pos)-1]
									conserved_file.write(human_gene+"\t"+human_protein_pos+"\t"+fish_gene+"\t"+str(subject_pos)+"\t"+human_amino_acid+"\t"+match_amino_acid+"\n")
				
				

conserved_file.close()

