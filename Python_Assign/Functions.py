'''
Contains some of the functions called in the Main python file 

'''










import os

import pandas as pd 

import tkinter as tk
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
from PIL import Image






def data_conversion(data_frame, col_name):
    data_frame[col_name] = pd.to_numeric(data_frame[col_name], errors='coerce')
    data_frame[col_name] = data_frame[col_name].astype(pd.Int64Dtype())
    return data_frame
'''
 Converts specified column to numeric while removing any missing values then converting to Int64
 I am honestly unsure why I did this; got certain iteration errors and inaccurate entrez_ids if not done
 Which should have been solved if it was converted to string...
'''

def sorting_function(column_name): 
    if column_name.startswith("Sample"):
        return int(column_name[6:]) 
    return 0 
'''
This sorts the columns in the Islets vst dataset only.
Take in all column names starting with 'Sample', then taking the 6th character and after meaning
all the numerical parts of it, turning it into an interger for easier ordering
Any columns not starting with 'Sample' are designated as 0, ensuring they reamin at start of table
'''


def merging_function(df1, df2):
    merged_data = df1.merge(df2, how='outer', indicator=False)  
     
    return merged_data
# Merges two given df and returns missing values as NaN
# The merge column is not included (indicator=False)     
        
  
from Bio import Entrez
# Access the NCBI's Entrez database system 

Entrez.email = "nimantha.senavirathne@gmail.com" # provide email to contact if any issues arise

def get_gene_summary(entrez_id):
    
    handle = Entrez.esummary(db="gene", id=entrez_id)
    record = Entrez.read(handle)
    handle.close()
    return record["DocumentSummarySet"]['DocumentSummary']
'''
Function within the Entrez module to retrieve document summaries from gene database according to entrez_id
The handle contains the data retrieved (xml format) which is then read and stored in record before closing the 
handle 
The record is in key/value pairs, and returns specific sections of the data> documentSummarySet is one key which 
has the DocumentSummary that is required in this case 
It is in a dict form containing the summary themselves 
'''




def format_gene_summary(gene_summaries):
    formatted_summaries = []

    for summary in gene_summaries:
        formatted_summary = []
        for key, value in summary.items():
            if key == 'LocationHist':
                continue  # Skip the Key (LocationHist)
            formatted_summary.append(f"{key}: {value}")
        formatted_summaries.append('\n'.join(formatted_summary))

    return '\n'.join(formatted_summaries)

'''
Iterate over each dictionary in the gene_summaries list obtained earlier. 
and for each dictionary/summary iterated, creates an empty list (formated_summary) to collect the 
formatted key-value pairs of the current summary.

Then it iterates over each key-value pair in the current summary dictionary.
And ignores the key 'LocationHist', prevents its addition to the formatted summary by using the 
continue statement; is not required in the final output.

Every key:value pair is appened to the formated_summary as a string in the form key: value

Once all key-value pairs of the current summary are processed, it is merged 
into a single string with each key-value pair on a new line ('\n'.join(formatted_summary')). 
which is then appeneded to the final formatted_summaries list.

All summaries are processed and seperated by a newline before returning as the final output



'''

"""
# Example Entrez ID
entrez_id = "1921682"
  


gene_summaries = get_gene_summary(entrez_id)


if gene_summaries:
    
    formatted_summaries = format_gene_summary(gene_summaries)
    
    print(formatted_summaries)
else:
    print("No summaries found for Entrez ID:", entrez_id)
"""

def write_gene_summaries_to_file(gene_summaries, filename):
    with open(filename, "w") as file:
        for block in gene_summaries:
            file.write(block + "\n\n")  # Add a newline between blocks
        print(f"All blocks written to {filename}")

'''
Opens the specified file in write mode ("w"). 
If it does not exist, it will be created. If it does exist, its contents will be overwritten.
Iterates over each formatted summary block in the gene_summaries list and file.write(block + "\n\n") 
will write each block to the file, followed by two newline characters ("\n\n"); a space between each block
'''


file_path_A = "gene_summaries.txt"
file_path_B = "highest_expression.txt"
file_path_C = "lowest_expression.txt"


def open_files():
    try:
        os.startfile(file_path_A)
        os.startfile(file_path_B)
        os.startfile(file_path_C)
    except Exception as e: 
        print(f"An error of type {type(e).__name__} occurred: {str(e)}")
# This opens the written files from earlier

file_path_I = ["Gene Expression Heatmaps.png", "Gene Expression Heatmaps_2.png"] 



def open_images(file_path_I):
    try: 
        for file_path in file_path_I:
            img = Image.open(file_path)
            img.show() 
    except Exception as e:
        print(f"An error of type {type(e).__name__} occurred: {str(e)}")
# This opens the created images from Main python file









