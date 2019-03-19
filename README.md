# PyFuncover
PyFuncover : Full proteome search for a specific function using BLAST and PFAM. 

# Abstract :
Python Function uncover ( PyFuncover ) is a new bioinformatic tool able to search for protein with a specific function in a full proteome. The pipeline coded in python uses BLAST alignment and the sequences from a PFAM family as search seed. The methodology is based on Ada-BLAST which is no longer available. We tested PyFuncover using the FABP family Lipocalin_7 from PFAM (version 32, 2019) against the Homo sapiens NCBI proteome. After applying the scoring function in all the BLAST results, the data was classified and submitted to a GO-TERM analysis using bioDBnet. Analysis showed that all family of FABPs were ranked in the 900 and plus protein. Above this threshold were found families able to bind to hydrophobic molecules similar to fatty acid such as the retinol acid transporter and the cellular retinoic acid-binding protein.

# Workflow
![alt text](https://github.com/Tuisto59/PyFuncover/blob/master/image2.png "Workflow")

# Result

![alt text](https://github.com/Tuisto59/PyFuncover/blob/master/image1.png "Workflow")

# OS
  * Windows
  * Unix

# Python
  * Python 2.7+

# Library

  * Numpy
  * Pandas
  * BioPython
  * Matplotlib

# Software

  * NCBI-BLAST+

# Usage

### TO UPDATE THE DATABASES :

python PyFuncover.py --update

### === OBLIGATORY ARGUMENTS : ===

-pfam : List of PFAM familly ID : PF####
  each separated by a blank space

PyFuncover.py -pfam PF14651 PF#### ...


-taxid: The list of TaxID for each organism you want to download a proteome
Each separated by a space (for example Human and Yeast taxid)

  PyFuncover.py -taxid 9606 559492

Can be a Taxid that represent a node in the phylogenetic tree
(Eukaryotes : 2759 ; Insecta : 50557, ...)
He will retrieve all availlable assembly for them

## === OPTIONAL : ===

--update :
Download the last release of the NCBI Taxonomic Database
Download the last RefSeq, Prokaryote and Eukaryote Genome Assembly List

PyFuncover.py --update

--out : Filename output
Format are in CSV format (pandas.to_csv output)
Default : result.csv

--nb-blast : The number of parrallelized BLAST process (default : 10)
Be carefull, high number will use lot of memory and create a stck overflow !

--db : The list of choosen cross-ref number to retrieve data from bioDBnet database :
        default : 137 45 46 47 (UNIPROT ID, GO-TERMs Databases)

    WARNING ! : Too high number of requested cross-refs will occur a slow-mode request 1 by 1.
                If it get an error for 1 request with 1 protein in the mode 1 by 1,
                the program will ABORT with a too high number of choosen cross-ref exception !

--nb-prot : The number of protein per request to the bioDBnet Database:

    WARNING ! : Too high number of requested cross-refs will occur a slow-mode request 1 by 1.
                If it throw an error for 1 request in this mode, the program will ABORT with a
                too high number of cross ref choosen exception !

1. : Affy ID
2. : Agilent ID
3. : Allergome Code
4. : ApiDB_CryptoDB ID
5. : Biocarta Pathway Name
6. : BioCyc ID
7. : CCDS ID
8. : Chromosomal Location
9. : CleanEx ID
10. : CodeLink ID
11. : COSMIC ID
12. : CPDB Protein Interactor
13. : CTD Disease Info
14. : CTD Disease Name
15. : CYGD ID
16. : dbSNP ID
17. : dictyBase ID
18. : DIP ID
19. : DisProt ID
20. : DrugBank Drug ID
21. : DrugBank Drug Info
22. : DrugBank Drug Name
23. : EC Number
24. : EchoBASE ID
25. : EcoGene ID
26. : Ensembl Biotype
27. : Ensembl Gene ID
28. : Ensembl Gene Info
29. : Ensembl Protein ID
30. : Ensembl Transcript ID
31. : FlyBase Gene ID
32. : FlyBase Protein ID
33. : FlyBase Transcript ID
34. : GAD Disease Info
35. : GAD Disease Name
36. : GenBank Nucleotide Accession
37. : GenBank Nucleotide GI
38. : GenBank Protein Accession
39. : GenBank Protein GI
40. : Gene ID
41. : Gene Info
42. : Gene Symbol
43. : Gene Symbol and Synonyms
44. : Gene Symbol ORF
45. : Gene Synonyms
46. : GeneFarm ID
47. : GO - Biological Process
48. : GO - Cellular Component
49. : GO - Molecular Function
50. : GO ID
51. : GSEA Standard Name
52. : H-Inv Locus ID
53. : HAMAP ID
54. : HGNC ID
55. : HMDB Metabolite
56. : Homolog - All Ens Gene ID
57. : Homolog - All Ens Protein ID
58. : Homolog - All Gene ID
59. : Homolog - Human Ens Gene ID
60. : Homolog - Human Ens Protein ID
61. : Homolog - Human Gene ID
62. : Homolog - Mouse Ens Gene ID
63. : Homolog - Mouse Ens Protein ID
64. : Homolog - Mouse Gene ID
65. : Homolog - Rat Ens Gene ID
66. : Homolog - Rat Ens Protein ID
67. : Homolog - Rat Gene ID
68. : HomoloGene ID
69. : HPA ID
70. : HPRD ID
71. : HPRD Protein Complex
72. : HPRD Protein Interactor
73. : Illumina ID
74. : IMGT/GENE-DB ID
75. : InterPro ID
76. : IPI ID
77. : KEGG Disease ID
78. : KEGG Gene ID
79. : KEGG Orthology ID
80. : KEGG Pathway ID
81. : KEGG Pathway Info
82. : KEGG Pathway Title
83. : LegioList ID
84. : Leproma ID
85. : Locus Tag
86. : MaizeGDB ID
87. : MEROPS ID
88. : MGC(ZGC/XGC) ID
89. : MGC(ZGC/XGC) Image ID
90. : MGC(ZGC/XGC) Info
91. : MGI ID
92. : MIM ID
93. : MIM Info
94. : miRBase ID
95. : NCIPID Pathway Name
96. : NCIPID Protein Complex
97. : NCIPID Protein Interactor
98. : NCIPID PTM
99. : Orphanet ID
100. : PANTHER ID
101. : Paralog - Ens Gene ID
102. : PBR ID
103. : PDB ID
104. : PeroxiBase ID
105. : Pfam ID
106. : PharmGKB Drug Info
107. : PharmGKB Gene ID
108. : PIR ID
109. : PIRSF ID
110. : PptaseDB ID
111. : PRINTS ID
112. : ProDom ID
113. : PROSITE ID
114. : PseudoCAP ID
115. : PubMed ID
116. : Reactome ID
117. : Reactome Pathway Name
118. : REBASE ID
119. : RefSeq Genomic Accession
120. : RefSeq Genomic GI
121. : RefSeq mRNA Accession
122. : RefSeq ncRNA Accession
123. : RefSeq Nucleotide GI
124. : RefSeq Protein Accession
125. : RefSeq Protein GI
126. : Rfam ID
127. : RGD ID
128. : SGD ID
129. : SMART ID
130. : STRING Protein Interactor
131. : TAIR ID
132. : Taxon ID
133. : TCDB ID
134. : TIGRFAMs ID
135. : TubercuList ID
136. : UCSC ID
137. : UniGene ID
138. : UniProt Accession
139. : UniProt Entry Name
140. : UniProt Info
141. : UniProt Protein Name
142. : UniSTS ID
143. : VectorBase Gene ID
144. : VEGA Gene ID
145. : VEGA Protein ID
146. : VEGA Transcript ID
147. : WormBase Gene ID
148. : WormPep Protein ID
149. : XenBase Gene ID
150. : ZFIN ID
55. : HMDB Metabolite,
56. : Homolog - All Ens Gene ID,
57. : Homolog - All Ens Protein ID,
58. : Homolog - All Gene ID,
59. : Homolog - Human Ens Gene ID,
60. : Homolog - Human Ens Protein ID,
61. : Homolog - Human Gene ID,
62. : Homolog - Mouse Ens Gene ID,
63. : Homolog - Mouse Ens Protein ID,
64. : Homolog - Mouse Gene ID,
65. : Homolog - Rat Ens Gene ID,
66. : Homolog - Rat Ens Protein ID,
67. : Homolog - Rat Gene ID,
68. : HomoloGene ID,
69. : HPA ID,
70. : HPRD ID,
71. : HPRD Protein Complex,
72. : HPRD Protein Interactor,
73. : Illumina ID,
74. : IMGT/GENE-DB ID,
75. : InterPro ID,
76. : IPI ID,
77. : KEGG Disease ID,
78. : KEGG Gene ID,
79. : KEGG Orthology ID,
80. : KEGG Pathway ID,
81. : KEGG Pathway Info,
82. : KEGG Pathway Title,
83. : LegioList ID,
84. : Leproma ID,
85. : Locus Tag,
86. : MaizeGDB ID,
87. : MEROPS ID,
88. : MGC(ZGC/XGC) ID,
89. : MGC(ZGC/XGC) Image ID,
90. : MGC(ZGC/XGC) Info,
91. : MGI ID,
92. : MIM ID,
93. : MIM Info,
94. : miRBase ID,
95. : NCIPID Pathway Name,
96. : NCIPID Protein Complex,
97. : NCIPID Protein Interactor,
98. : NCIPID PTM,
99. : Orphanet ID,
100. : PANTHER ID,
101. : Paralog - Ens Gene ID,
102. : PBR ID,
103. : PDB ID,
104. : PeroxiBase ID,
105. : Pfam ID,
106. : PharmGKB Drug Info,
107. : PharmGKB Gene ID,
108. : PIR ID,
109. : PIRSF ID,
110. : PptaseDB ID,
111. : PRINTS ID,
112. : ProDom ID,
113. : PROSITE ID,
114. : PseudoCAP ID,
115. : PubMed ID,
116. : Reactome ID,
117. : Reactome Pathway Name,
118. : REBASE ID,
119. : RefSeq Genomic Accession,
120. : RefSeq Genomic GI,
121. : RefSeq mRNA Accession,
122. : RefSeq ncRNA Accession,
123. : RefSeq Nucleotide GI,
124. : RefSeq Protein Accession,
125. : RefSeq Protein GI,
126. : Rfam ID,
127. : RGD ID,
128. : SGD ID,
129. : SMART ID,
130. : STRING Protein Interactor,
131. : TAIR ID,
132. : Taxon ID,
133. : TCDB ID,
134. : TIGRFAMs ID,
135. : TubercuList ID,
136. : UCSC ID,
137. : UniGene ID,
138. : UniProt Accession,
139. : UniProt Entry Name,
140. : UniProt Info,
141. : UniProt Protein Name,
142. : UniSTS ID,
143. : VectorBase Gene ID,
144. : VEGA Gene ID,
145. : VEGA Protein ID,
146. : VEGA Transcript ID,
147. : WormBase Gene ID,
148. : WormPep Protein ID,
149. : XenBase Gene ID,
150. : ZFIN ID,
