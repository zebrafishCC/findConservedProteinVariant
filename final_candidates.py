f1 = open("zebrafish_selected_genes.txt","r")
f2 = open("zebrafish_human_common_set.txt","r")
f3 = open("zebrafish_id_conversion.txt","r")
f4 = open("human_id_conversion.txt","r")

fishgenes = []
for line in f1.readlines():
	fishgene = line.strip()
	fishgenes.append(fishgene)
print len(fishgenes)
print fishgenes[:10]

humangenes = []
genesPair = []
for line in f2.readlines():
	fishgene = line.strip().split("\t")[0]
	humangene = line.strip().split("\t")[2]
	if fishgene in fishgenes:
		humangenes.append(humangene)
		genesPair.append(fishgene+humangene)
		#print fishgene+"\t"+humangene
print len(humangenes)
print humangenes[:10]

fishTranscript = []
for line in f3.readlines():
	thisline = line.strip().split("\t")
	gene = thisline[2]
	transcript = thisline[1]
	if gene in fishgenes:
		fishgenes.remove(gene)
		fishTranscript.append(gene+"\t"+transcript)

print len(fishTranscript)

humanTranscript = []
for line in f4.readlines():
        thisline = line.strip().split("\t")
        gene = thisline[2]
        transcript = thisline[1]
        if gene in humangenes:
                humangenes.remove(gene)
                humanTranscript.append(gene+"\t"+transcript)


print len(humanTranscript)

f5 = open("fish_human_final_candidate.txt","w")
for fishGeneTranscript in fishTranscript:
	fishgene = fishGeneTranscript.split("\t")[0]
	for humanGeneTranscript in humanTranscript:
		humangene = humanGeneTranscript.split("\t")[0]
		fishHuman = fishgene+humangene
		if fishHuman in genesPair:
			f5.write(fishGeneTranscript+"\t"+humanGeneTranscript+"\n")

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
