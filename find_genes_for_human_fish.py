import sys

f1 = open(sys.argv[1], "r") #file 1 is the human genes with zebrafish homologues in variant file
f2 = open(sys.argv[2],"r") #file2 is zebrafish genes files that need to match the human genes
f3 = open(sys.argv[3],"w") #file3 write the common set genes into file


genes = []
for line in f1.readlines():
	gene = line.strip().split("\t")[1] #get human gene name from file1
	if gene not in genes:
		genes.append(gene)


previous = ""
#one human gene may coresspond to multiple zebrafish gene since human 3269 while fish 3261
fishgenes = []

for line in f2.readlines():
	gene = line.strip().split("\t")[2] #human gene
	fishgene = line.strip().split("\t")[0]
	if gene in genes:
		genes.remove(gene) #only use human gene name once
		if ":" not in fishgene: #remove zebrafish strange gene name
			if "." not in fishgene:
				if fishgene not in fishgenes:
					fishgenes.append(fishgene)
					f3.write(line)

f1.close()
f2.close()
f3.close()
		
