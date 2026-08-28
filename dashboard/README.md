# Analysis Summary Dashboard

## Summary
MetaSieve comes with a dashboard for visualization of the Kraken2, SeqScreen, and ESMFold results. The dashboard requires output from all three programs:
1) Kraken2 report file (not output file)
2) SeqScreen CSV file
3) ESMFold Protein Structure PDB Files

From these files, the dashboard will provide multiple figures and statistics summarizing the combined information. Using [Chimera](https://www.cgl.ucsf.edu/chimerax/download.html), the dashboard also provides protein structure prediction and visualization. 


## Steps 

1. First, the user is required to import the output files of the main pipeline.


  <img height="300" alt="image" src="https://github.com/user-attachments/assets/342c05c3-d567-4e6d-8fba-3ba6f8e7244a" />

 
2. Demonstration Data could be used

 <img height="300" alt="image" src="https://github.com/user-attachments/assets/46a49b03-2fb2-42fb-a7c2-af771315927c" />


3.  The Kraken2 analysis summary features two main charts: a pie chart showing the dataset-wide taxonomic breakdown (including the 78.8\% unclassified reads count alongside bacteria, eukaryota, and archaea percentages) and a bar chart showing the read distribution per sample across all 21 analyzed samples.

   
   <img height="300" alt="image" src="https://github.com/user-attachments/assets/61168581-d8e3-48a3-b73e-e0b544d7faab" />
   
4. Next, the SeqScreen analysis section includes three charts and one table: the first chart displays the open reading frame count distribution across six translation frames, the second chart shows the strand orientation breakdown per frame, and the third chart presents the amino acid length distribution by frame. 

   <img width="1874" height="890" alt="image" src="https://github.com/user-attachments/assets/67e3494b-3dfe-4df7-83a6-33d84473d741" />

Followed by the SeqScreen output table detailing individual sequence records and annotations.
    <img width="1848" height="554" alt="image" src="https://github.com/user-attachments/assets/295ed376-0275-439f-864f-fe77b03f083c" />

5. Next, the ESMFold analysis summary section features a table displaying the final retained list of PDB models and their structural metrics, complete with an option to export the data as a CSV file.

   
   <img height="300" alt="image" src="https://github.com/user-attachments/assets/03c97a60-e842-4b10-b44a-cd326db41839" />

  Demo CSV file : https://github.com/collaborativebioinformatics/Metasieve/blob/main/dashboard/outputs/candidates.csv
  
6. The final section allows you to visualize each of these structures in 3D.


<img height="300" alt="image" src="https://github.com/user-attachments/assets/d1196c79-3540-49ea-9c6d-0187e3cf85d6" />

7. Finally, once everything is complete, you can download all the data displayed across the entire dashboard as a single PDF report.


   <img height="300" alt="image" src="https://github.com/user-attachments/assets/7561820d-d80e-4b45-9643-05bbf3f411b8" />

   Demo PDF file : https://github.com/collaborativebioinformatics/Metasieve/blob/main/dashboard/outputs/Metasieve_Analysis_Summary.pdf 

# How-to-Run

1. Download this directory ( https://github.com/collaborativebioinformatics/Metasieve/tree/main/dashboard) 

2. Download ChimeraX-1.12.exe (from https://www.cgl.ucsf.edu/chimerax/download.html) in the same folder folder and set it as cd
   
3. If demonstration data is needed, download the data and add the folder to cd : [https://drive.google.com/drive/folders/1jrNid351AQ11GuhHgEC8lHbqF3SifnLl?usp=sharing ](https://drive.google.com/drive/folders/1jrNid351AQ11GuhHgEC8lHbqF3SifnLl?usp=sharing)
4. Run main.py
