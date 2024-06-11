# BIO-3027-Python-Assignment
Assignment/Exam project for BIO-3027

This is the RNA-seq analysis pipeline carried out on frozen pancreatic/liver/skeletal tissue samples from Diversity Outbred mice.

The run code is in the Main.py:
Running it will open a basic GUI for choosing which dataset to analyse (pancreatic/liver/skeletal) before clicking the button to start the analysis

It should take 2-3 mins to complete due to delays introduced when requesting information from Entrez database to lessen traffic. 

Completion occurs when the results (3 text files & 2 png files are opened for viewing) 

The text files contain gene summary information while the two pngs are heatmap images of the genes. 

It is advised to let the program go to completion. 
All Files in the folder Python_Assign should be left as is lest errors arise 
Can refer to Assignment.pdf for more information


Created on Spyder ver 5.4.3


IMPORTANT: After re-running the code at 3 weeks post uploading, the PIL function for opening images seems to have hit an issue. Previously it would open the two png result files at the same time, now however, it only opens one, which must be exited for the second png file to open. I am unsure as to what changed in the 3 weeks and I truly cannot find a solution for this. There are no issues in the analysis, reading the raw data files and creating of the results both the pngs and textfiles, only opening of both the png files.  
