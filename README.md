# PyFuncover
PyFuncover : Full proteome search for a specific function using BLAST and PFAM. 

# Abstract :
Python Function Uncover ( PyFuncover ) is a new bioinformatic tool able to search for protein with a specific function in a full proteome. The pipeline coded in python uses BLAST alignment and the sequences from a PFAM family as search seed. The methodology is based on Ada-BLAST which is no longer available. We tested PyFuncover using the FABP family Lipocalin_7 from PFAM (version 32, 2019) against the Homo sapiens NCBI proteome. After applying the scoring function in all the BLAST results, the data was classified and submitted to a GO-TERM analysis using bioDBnet. Analysis showed that all family of FABPs were ranked in the 900 and plus protein. Above this threshold were found families able to bind to hydrophobic molecules similar to fatty acid such as the retinol acid transporter and the cellular retinoic acid-binding protein.
