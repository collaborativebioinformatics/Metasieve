# Metasieve Report Dashboard
## Plan :
generate an interactive dashboard that will create an interactive report of the “De Novo seq” from metagenomics database analysis using the “MetaSieve” pipeline :

1. ui request for input files paths = kraken output file, seqscreen output file; esmfold output folder containing the .pdb files
2. if not data we have “Demonstration” option that uses data from ALFAFold database classified dark matter + generate simulation outputs for kraken2, seqscreen, esmfold // small dataset of 500 
3. It generates txt/pdf report  summary of the analysis , displayed in UI + option to download
4. Extract data of the “Unclassified” seqs form kraken output file, seqscreen output file; esmfold output + combine into a df + download as excel
5. Option Structural prediction using Chimera for selected seq

# Output - V0 version
1. Summarize the stats all the analyses  
   <img width="1025" height="546" alt="img1" src="https://github.com/user-attachments/assets/847fe7a1-b85c-4905-aacc-144d53d16a6a" />

2. Show the final unclassified sequences and their relevant data
   <img width="1859" height="442" alt="img2" src="https://github.com/user-attachments/assets/ccc0da42-1ec2-43bd-9e81-20b8532b2d48" />

3. Bonus : Uses Chimera to display predict structure 

<img width="1578" height="897" alt="img3" src="https://github.com/user-attachments/assets/c3711543-36f3-459b-8601-2a18d1b2e8aa" />

# How-to-Run

1. Download this folder

2. Download ChimeraX-1.12.exe from https://www.cgl.ucsf.edu/chimerax/download.html
   
3. Put ChimeraX-1.12.exe  in downloaded folder and set it as cd
   
4. Run main.py

Note : the current code uses test data for visualization purposes,  later it will use the main pipeline output files 
