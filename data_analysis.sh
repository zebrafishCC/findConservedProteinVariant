#This file contains the major steps for the project. Some important statistics are also shown in the data set.


#screen for C T data
awk '$4=="C" && $5=="T" {print $0}' clinvar_20181202.vcf | wc -l
89225

awk '$4=="G" && $5=="A" {print $0}' clinvar_20181202.vcf | wc -l
88903

#only select variant file that is C to T
awk '$4=="C" && $5=="T" {print $0}' clinvar_20181202.vcf > humanCtoTvariant.vcf

#remove synonymous_variant 
grep -v "synonymous_variant" humanCtoTvariant.vcf > humanCtoTvariant_rm_synonymous.vcf
66080

#get human-fish orthologues from databases
awk -F "\t" '{print $4}' human_orthos_2018.12.06.txt | sed '1,2d' | sort | uniq| wc -l 
13187


#get human variant with zebrafish homologues
awk -F "\t" '{print $4}' human_orthos_2018.12.06.txt | sed '1,2d' | sort | uniq | awk '{print "GENEINFO="$0":"}'|xargs -I {}  grep {} humanCtoTvariant_rm_synonymous.vcf > humanCtoTvariant_with_homologues.vcf

61455

grep "missense" humanCtoTvariant_with_homologues.vcf | wc -l
38194

grep "nonsense" humanCtoTvariant_with_homologues.vcf | wc -l
4875

#get only coding region sequences
grep -E 'missense|nonsense' humanCtoTvariant_with_homologues.vcf | wc -l
43053

grep -E 'missense|nonsense' humanCtoTvariant_with_homologues.vcf > humanCtoTvariant_with_homologues_coding.vcf
43053
#use ensembl variant effect predictor to find the cds position for point mutation
mv EFaq8oyvRG1fl175.Consequence_is_missense_variant_or_Consequence_is_stop_gained.txt VEFhuman.txt

#extract useful information
cut -f 2,6,7,9,16,17,18,19 VEFhuman.txt >VEFhuman_extract_information.txt 

#only select the protein with the longest length for annotation
cat VEFhuman_extract_information.txt | awk '$1!=prev {print $0} {prev=$1}' | wc -l
42987

#store the final human variant information before comparing with zebrafish infomation
cat VEFhuman_extract_information.txt | awk '$1!=prev {print $0} {prev=$1}' > VEFhumanVariant.txt

cut -f 2,3,4,5 human_orthos_2018.12.06.txt | sort -k1 | uniq > zebrafish_human_orthologues.txt
16328

cut -f 2 VEFhumanVariant.txt | uniq | wc -l
4515


head zebrafish_human_orthologues.txt 
a1cf	apobec1 complementation factor	A1CF	APOBEC1 complementation factor
aaas	achalasia, adrenocortical insufficiency, alacrimia	AAAS	aladin WD repeat nucleoporin
aacs	acetoacetyl-CoA synthetase	AACS	acetoacetyl-CoA synthetase
aadac	arylacetamide deacetylase	AADAC	arylacetamide deacetylase
aadacl4	arylacetamide deacetylase-like 4	AADACL4	arylacetamide deacetylase like 4
aadat	aminoadipate aminotransferase	AADAT	aminoadipate aminotransferase
aagab	alpha and gamma adaptin binding protein	AAGAB	alpha and gamma adaptin binding protein
aak1a	AP2 associated kinase 1a	AAK1	AP2 associated kinase 1
aak1b	AP2 associated kinase 1b	AAK1	AP2 associated kinase 1
aamdc	adipogenesis associated, Mth938 domain containing	AAMDC	adipogenesis associated Mth938 domain containing



#find common genes from human zebrafish orthologues
ipython find_genes_for_human_fish.py VEFhumanVariant.txt zebrafish_human_orthologues.txt zebrafish_human_common_set.txt 

cat zebrafish_human_common_set.txt | wc -l
3200

#It should be note that ensembl gene ID transversion lags behind ZFIN annotation
cut -f 3 zebrafish_id_conversion.txt | sed '1d' | sort | uniq | wc -l
3091

#and for no reason bioMart convert some gene to UPPER case aftertransversion
ACSF3
ADAMTSL4
AGBL1
ARSB
ASS1
BFSP1
CCKAR
CLPB
COLEC10
DNAL1
DST
FAT4
GALNT3
GAN
GRIK2
GRXCR1
KCNJ6
LHX3
MYO9B
NPC1L1
PLCB4
RGS9BP
SEMA4D
SLC6A13
SNTA1
SPTBN4
TCF4
TMC1
TMEM216
TRAPPC9
TUBB4B
XRCC2
ZNF335
ZNF423
ZNF462
 
#final number
comm -12 zebrafish_test.txt zebrafish_test1.txt > zebrafish_test2.txt
3056
python final_candidates.py

fish_human_final_candidate.txt
wdr34	ENSDART00000080347	WDR34	ENST00000372715
gria3a	ENSDART00000137120	GRIA3	ENST00000616590
creb3l1	ENSDART00000024330	CREB3L1	ENST00000621158
cdc42	ENSDART00000155746	CDC42	ENST00000344548
syt2a	ENSDART00000135009	SYT2	ENST00000367267
adamts13	ENSDART00000146364	ADAMTS13	ENST00000628339
gnal	ENSDART00000140924	GNAL	ENST00000334049
tmem240a	ENSDART00000129910	TMEM240	ENST00000378733
osmr	ENSDART00000171701	OSMR	ENST00000502536
dnah3	ENSDART00000183613	DNAH3	ENST00000261383

#problems: uniq transcript retrieve from ensembl bioMart doens't show records in VEFhumanVariant.txt, which cause a great reduce of analysis variant
#only 1388 human transcripts with variant for blast analysis



#only search for codon with GC pair (C to T transition with PAM NGG) or CG pair (G to A transition with PAM CCN)

#GC pair amino acid change from V-->A
#CG pair amino acid change from A-->T,C-->Y,S-->N,R-->H,G-->D

no records of V-->A transition in VEFhumanVariant.txt
grep "A/T\|C/Y\|S/N\|R/H\|G/D" VEFhumanVariant.txt > VEFhumanVariantGtoA.txt
7837 records with 1330 Transcripts

cut -f 4 VEFhumanVariantGtoA.txt| sort | uniq > human_transcript_withGC.txt


python find_conserved_transcript.py

cat human_zebrafish_conserved_variant.txt | wc -l
1155
cut -f 1 human_zebrafish_conserved_variant.txt | sort | uniq | wc -l
352

python find_chromo_pos.py 

#

