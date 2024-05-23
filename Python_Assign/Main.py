'''
This file contains the instructions to be carried out on the csv data in the folders 'islets', 'liver' and 'skeletal'
Contains high thoughput RNA-seq data carried out on frozen pancreatic/liver/skeletal tissue samples from Diversity Outbred mice.

Some pre-processing carried out on the datasets including ordering certain columns, getting average expression
across genes, merging the two datasets and dropping redundant columns etc

Make gene expression heatmap plots for 20 random genes as well as their average expression across all samples
Manipulate merged dataset to obtain highest and lowest averagely expressed genes and plot heatmap plots for each
Save the plots as png files while also opening them for view once completed.

Use Entrez to obtain gene data/summary for 10 genes from each of the above plots while also processing the 
obtained text into readable .txt files that too are opened once completed.

The Entrez section takes time; introduced delay between getting results for each gene as well as between different plots
Otherwise errors are common due to traffic in the NCBI website from which the data is obtained

Hence completion of this takes time; 1-3 mins max.  

The data for the study was obtained from:

--> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266921
--> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266923
--> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266922

Contributor(s):	Churchill GA, Baker CN     

Title:  RNA-seq expression analysis of liver, adipose, skeletal muscle, heart, and pancreatic islets 
        from DO founder strain mice fed high-fat, high-sugar diet. [Islet][Liver][Skeletal]

Status:	Public on May 13, 2024

Chosen Files: mrna_annot & vst_matrix

PS. also do not carry out any data convertions on the raw csv data file lest errors show up

'''





import os # Use this to interact with the operating system/ mostly for opening files/images/text etc

import pandas as pd # For its dataframe structures and general manipulation of these tables, mainly the raw csv data 

import tkinter as tk # For the window python file; allows the creation ofa  gui for chosing which datasets to run
import matplotlib.pyplot as plt # For making plots and subplots to visualise gene expression and saving images
import seaborn as sns # For making heatmaps to visualise the gene expressions


from tkinter import ttk # For making dropdown lists in the tkinter window, etc
import subprocess # To run entire python files from another file

  
from Bio import Entrez # To obtain gene information from the entrez databases with the use of the entrez_id 
import time # To introduce delay between starting certain code snippets
from PIL import Image # To open image files




def path_file(value):
    if value == 'Islet':
        annotation_file = os.path.join('pancreatic islets', 'GSE266921_mrna_annots.csv')
        vst_file = os.path.join('pancreatic islets', 'GSE266921_vst_matrix.csv' )
    elif value == 'Liver':
        annotation_file = os.path.join('liver', 'GSE266923_mrna_annots.csv')
        vst_file = os.path.join('liver', 'GSE266923_vst_matrix.csv' )
    elif value == 'Skeletal':
        annotation_file = os.path.join('skeletal', 'GSE266922_mrna_annots.csv')
        vst_file = os.path.join('skeletal', 'GSE266922_vst_matrix.csv' )
    else:
        raise ValueError(f"Unknown value: {value}. Expected 'Islet', 'Liver', or 'Skeletal'.")
        
    data_mrna_annot = pd.read_csv(annotation_file, float_precision='round_trip')
    data_vst_matrix = pd.read_csv(vst_file)
    
    return data_mrna_annot, data_vst_matrix 

# Returns the read csv files according to the value which is selected from the dropdown menu below




        
root = tk.Tk()
root.title("RNA-seq Expression Analysis")
root.geometry("300x170+500+200")
root.config(bg='gray19')
# Initilises the tkinter window/ name, size and dimension and color

items = ["Islet", "Liver", "Skeletal"]

selected_item = tk.StringVar()
# creates a StringVar type of variable that can be used to store and retrieve string values.

combo_box = ttk.Combobox(root, textvariable=selected_item, values=items)
combo_box.pack(pady=10)
# Merge the StringVar to the dropdown combobox whose values are the item list and within the root window

##########################
# Main Analysis Function #
##########################

def run_analysis():    
    try:
        data_mrna_annot, data_vst_matrix = path_file(selected_item.get())
        # Get the vst and annot data according to the selected item from dropdown list by following the path_file function
        print(data_mrna_annot)
        print(data_vst_matrix)
    except Exception as e:
         print(f"An error of type {type(e).__name__} occurred: {str(e)}")    
    
    if selected_item.get() == "Islet":
        selected_item_name = 'Pancreatic Islet'
        # selected_item_list is to allow the name to be shown in the graphs created
        try:
            from Functions import sorting_function
            sorted_data = sorted(data_vst_matrix.columns, key=sorting_function)
            modified_data = data_vst_matrix[sorted_data]
            print(modified_data)
            # This sorts a selected number of columns (sample columns) into order 
            # Originaly it went: Sample1, Sample10, Sample100, Sample11 etc
            # This ordered it to Sample1, Sample2...Sample100 etc. 
        except Exception as e:
            print(f"An error of type {type(e).__name__} occurred: {str(e)}")
    elif selected_item.get() == "Liver":
        selected_item_name = 'Liver'
        modified_data = data_vst_matrix
        # No modification/sorting is done unlike in the Islet dataset
    elif selected_item.get() == "Skeletal":
        selected_item_name = 'Skeletal Muscle'
        modified_data = data_vst_matrix
        
    try:
        from Functions import data_conversion
        data_mrna_annot = data_conversion(data_mrna_annot, 'Entrez_ID')
        # Converts Entrez_ID values to int64
        from Functions import merging_function

        merged_dataset = merging_function(data_mrna_annot, modified_data)
    # This merged the two datasets; it automatically dropped redunant columns from what I have observed
    except Exception as e:
        print(f"An error of type {type(e).__name__} occurred: {str(e)}")
        
    if selected_item.get() == "Islet":
        valid_columns = merged_dataset.loc[:, 'Sample1':]
        # This gets the all the columns after column 'Sample1' which is the columns with the expression data
    elif selected_item.get() == "Liver" or selected_item.get() == "Skeletal":
        valid_columns = merged_dataset.loc[:, 'NZO_5':]
        # Gets all the columns after column 'NZO_5' which is the columns with the expression data (same for both skeletal and liver)
    try:
        analysis_data = pd.concat([merged_dataset[['Gene_symbol']], valid_columns], axis=1)
        print(analysis_data)
        # This creates a new dataframe containing the expression data joined to the gene name
        merged_dataset['Average_Expression'] = valid_columns.mean(axis=1)
        print(merged_dataset)
        # This adds a column containing the average gene expression across all samples for each gene
        # axis=1; mean is calculated for each row
        average_expression = pd.DataFrame(merged_dataset, columns = ['Gene_symbol', 'Average_Expression'])
        print(average_expression)
        # This creates a dataframe with only average expression and gene name, used later on


        data_analysis = analysis_data.copy()
        data_analysis_2 = merged_dataset.copy()
        # Creates copies of the above two dataframes (df) for further manipulation
        # As there is indexing and other changes done, wanted to keep the orginal data



        data_analysis.set_index('Gene_symbol', inplace=True)
        subset_genes = data_analysis.sample(20)
        print(subset_genes)
        # Sets Gene_symbol as index and modifies orginal data_analysis without creating a new df
        # Store 20 random genes from the data_analysis df in a new variable (no average expression data in this df)
        subset_genes_index = subset_genes.index
        # This stores the gene index names in another df 

        matched_rows = average_expression[average_expression['Gene_symbol'].isin(subset_genes_index)]
        matched_rows.set_index('Gene_symbol', inplace=True)
        print(matched_rows)
        # This creates a df with the average expression for the 20 random genes chosen earlier 

        highest_expression = average_expression.nlargest(20, 'Average_Expression')
        highest_expression.set_index('Gene_symbol', inplace=True)
        print(highest_expression)
        # This gets the 20 genes with the highest average expression  
        lowest_expression = average_expression.nsmallest(20, 'Average_Expression')
        lowest_expression.set_index('Gene_symbol', inplace=True)
        print(lowest_expression)
        # This gets the 20 genes with the lowest average expression 




        data_analysis_2 = pd.DataFrame(data_analysis_2, columns = ['Gene_symbol', 'Entrez_ID', 'Average_Expression'])
        print(data_analysis_2)
        # Drops all data except the columns given above
        matched_entrez = data_analysis_2[data_analysis_2['Gene_symbol'].isin(subset_genes_index)]
        print(matched_entrez)
        # Gets the entrez_ids from the randomly chosen 20 genes from before (to cross with the entrez database)

        highest_entrez = data_analysis_2.nlargest(20, 'Average_Expression')
        print(highest_entrez)
        lowest_entrez = data_analysis_2.nsmallest(20, 'Average_Expression')
        print(lowest_entrez)
        # Get the entrez id from the 20 highest and lowest expressed genes





############
# PLOTTING #
############

        fig_1, axes = plt.subplots(1, 2, figsize=(15, 8))
        # creates subplots arranged side by side of size 15'x 8'

        sns.heatmap(subset_genes, cmap='viridis', ax=axes[0])
        axes[0].set_title(f"Heatmap of Gene Expression Levels for 20 Random {selected_item_name} Genes")
        axes[0].set_xlabel('Samples')
        axes[0].set_ylabel('Gene')
        # Creates seaborn heatmap of dataset subset_genes on the axis[0]; first plot/left plot
        # Title and x/y labels given 

        sns.heatmap(matched_rows, cmap='viridis', ax=axes[1])
        axes[1].set_title(f"Heatmap of Average Gene Expression Levels for 20 Random {selected_item_name} Genes")
        axes[1].set_xlabel('')
        axes[1].set_ylabel('Gene')
        # Creates seaborn heatmap of dataset matched_rows on the axis[1]; second plot/right plot
        # x label not realy required as it takes from the dataset

        plt.tight_layout()
        plt.show()
        fig_1.savefig('Gene Expression Heatmaps.png')
        # .tight prevent overlapping between the two, shows the created plot and saves it in the same location as the 
        # python file

        fig_2, axes = plt.subplots(1, 2, figsize=(15, 8))


        sns.heatmap(highest_expression, cmap='viridis', ax=axes[0])
        axes[0].set_title(f"Heatmap of Highest Gene Expression Levels in {selected_item_name} Cells")
        axes[0].set_xlabel('Samples')
        axes[0].set_ylabel('Gene')

        heatmap = sns.heatmap(lowest_expression, cmap='viridis', ax=axes[1])
        axes[1].set_title(f"Heatmap of Lowest Gene Expression Levels in {selected_item_name} Cells")
        axes[1].set_xlabel('')
        axes[1].set_ylabel('Gene')
        heatmap.invert_yaxis()

        plt.tight_layout()
        plt.show()
        fig_2.savefig('Gene Expression Heatmaps_2.png')
        # Same as above plots, but instead compares the highest and lowest expressed genes while also saving them in 
        # same location
        

######################
# ENTREZ DATA SEARCH #
######################

        from Functions import get_gene_summary
        from Functions import format_gene_summary

        def process_gene_summaries(dataframe):
            summaries = []
            
            entrez_list = dataframe['Entrez_ID'].tolist()
            # Gets the entrez id values into a list 
            entrez_list = [value for value in entrez_list if pd.notna(value)]
            print(entrez_list)
            # Removes and NA values from the list 
            for value in entrez_list[:10]:
                # Gets first 10 entrez ids (0-9th) one at a time and carry out the following instructions
            #  if pd.notna(value):  # Check if value is not NA
                
                    
                try:
                    # Retrieve and format gene summary
                    value = str(value)
                    print(value)
                    # Transforms the int64 to a string value at every iteration
                    gene_summary = get_gene_summary(value)
                    formatted_summary = format_gene_summary(gene_summary)
                    # Gets the gene information and formats it (refer to the Functions python file)
                    # Store the formatted summary in the empty list from before
                    summaries.append(formatted_summary)
                except RuntimeError as e:
                    print(f"No results found for Entrez ID {value}: {e}")
                
                # Introduce a delay(3 second) before fetching the next summary (get a lot of runtime errors otherwise)
                time.sleep(3) 
            return summaries






        from Functions import write_gene_summaries_to_file
        # Refer to Functions.py for more information 
        try:
            summary_gene = process_gene_summaries(matched_entrez)
            filename_1 = "gene_summaries.txt"
            write_gene_summaries_to_file(summary_gene, filename_1)
            print("Summaries written to gene_summaries.txt")
        except Exception as e:
            print(e)
            
        time.sleep(15) # Introduces delay between requesting the different gene information (random/highest/lowest)

        try:
            high_gene_summaries = process_gene_summaries(highest_entrez)
            filename_2 = "highest_expression.txt"
            write_gene_summaries_to_file(high_gene_summaries, filename_2)
            print("Summaries written to highest_expression.txt")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
            
        time.sleep(15)

        try:
            low_gene_summaries = process_gene_summaries(lowest_entrez)
            filename_3 = "lowest_expression.txt"
            write_gene_summaries_to_file(low_gene_summaries, filename_3)
            print("Summaries written to lowest_expression.txt")
        except Exception as e:  
            print(f"An error of type {type(e).__name__} occurred: {str(e)}")




##############################
# RESULT OPENING AND VIEWING #
##############################


        from Functions import open_files
        from Functions import open_images
        try:
            open_files()
            file_path_I = ["Gene Expression Heatmaps.png", "Gene Expression Heatmaps_2.png"]             
            open_images(file_path_I)
        # Refer to the Functions file for more information

        # This opens both the newly written text files and created images for veiwing
        # Beware that the files (text/image) are modified everytime the different python files (Islets/Liver/Skeletal) are run
        # No new text/image files are created to prevent congestion of generated results. 
        except Exception as e:
            print(f"An error of type {type(e).__name__} occurred: {str(e)}") 
    except Exception as e:
        print(f"An error of type {type(e).__name__} occurred: {str(e)}") 
        
    
def handle_selection(event):
    #global selected_item_name
    #selected_item_name = selected_item.get()
    
    print("Selected item:", selected_item.get())
    #return selected_item_name
# Show the selected item on console  
# The event of selecting an item in the combobox results in the function being automatically called and printing 


combo_box.bind("<<ComboboxSelected>>", handle_selection)
# Binds the function from before to the combobox 

button = tk.Button(root, text="Run Script", command=run_analysis)
button.pack()
# Once the button is pressed, the run_analysis function is called and carried out giving the results

root.mainloop()
# Runs the UI application and keeps it from disappearing before the script is finished running
 








######################################################################################
# If want to visualise the individual variables at each level can run the code below #
#                       Remove """ """ from top and bottom                           #
#                           Only for Islets cells                                    #
######################################################################################





"""

data_file_mrna_annot = os.path.join('pancreatic islets', 'GSE266921_mrna_annots.csv')
data_file_vst_matrix = os.path.join('pancreatic islets', 'GSE266921_vst_matrix.csv' )
# This determines the relative file path of the data to be analysed in this program file and save it in the given variables





data_mrna_annot = pd.read_csv(data_file_mrna_annot,float_precision='round_trip')
data_vst_matrix = pd.read_csv(data_file_vst_matrix)

'''
This reads the two datasets using the pd module
The annot file contains data that needs to be kept as precise as possible (chr start/end/middle) hence use the float_precision 
to maintain precision by round-trip conversion,where it will try to convert floating-point numbers back and 
forth to strings.
'''


from Functions import data_conversion
data_file_mrna_annot = data_conversion(data_mrna_annot, 'Entrez_ID')
# The annot file also contains the entrez_id which is converted to a int64 value showcasing the full unchanged
# value as it is important to be kept as is, since it is used to identify the gene in the entrez databases
# No type of rounding can be done on this
# There was an issue with turning it into a string...I think




print(data_mrna_annot)
print(data_vst_matrix)
# Prints the two dataframes


# print("File path:", data_file_mrna_annot)

from Functions import sorting_function
sorted_data = sorted(data_vst_matrix.columns, key=sorting_function)
modified_data = data_vst_matrix[sorted_data]
print(modified_data)
# This sorts a selected number of columns (sample columns) into order 
# Originaly it went: Sample1, Sample10, Sample100, Sample11 etc
# This ordered it to Sample1, Sample2...Sample100 etc. 


from Functions import merging_function

merged_dataset = merging_function(data_mrna_annot, modified_data)
# This merged the two datasets; it automatically dropped redunant columns from what I have observed

valid_columns = merged_dataset.loc[:, 'Sample1':]
# This gets the all the columns after column 'Sample1' which is the columns with the expression data
analysis_data = pd.concat([merged_dataset[['Gene_symbol']], valid_columns], axis=1)
print(analysis_data)
# This creates a new dataframe containing the expression data joined to the gene name
merged_dataset['Average_Expression'] = valid_columns.mean(axis=1)
print(merged_dataset)
# This adds a column containing the average gene expression across all samples for each gene
# axis=1; mean is calculated for each row
average_expression = pd.DataFrame(merged_dataset, columns = ['Gene_symbol', 'Average_Expression'])
print(average_expression)
# This creates a dataframe with only average expression and gene name, used later on


data_analysis = analysis_data.copy()
data_analysis_2 = merged_dataset.copy()
# Creates copies of the above two dataframes (df) for further manipulation
# As there is indexing and other changes done, wanted to keep the orginal data



data_analysis.set_index('Gene_symbol', inplace=True)
subset_genes = data_analysis.sample(20)
print(subset_genes)
# Sets Gene_symbol as index and modifies orginal data_analysis without creating a new df
# Store 20 random genes from the data_analysis df in a new variable (no average expression data in this df)
subset_genes_index = subset_genes.index
# This stores the gene index names in another df 

matched_rows = average_expression[average_expression['Gene_symbol'].isin(subset_genes_index)]
matched_rows.set_index('Gene_symbol', inplace=True)
print(matched_rows)
# This creates a df with the average expression for the 20 random genes chosen earlier 

highest_expression = average_expression.nlargest(20, 'Average_Expression')
highest_expression.set_index('Gene_symbol', inplace=True)
print(highest_expression)
# This gets the 20 genes with the highest average expression  
lowest_expression = average_expression.nsmallest(20, 'Average_Expression')
lowest_expression.set_index('Gene_symbol', inplace=True)
print(lowest_expression)
# This gets the 20 genes with the lowest average expression 




data_analysis_2 = pd.DataFrame(data_analysis_2, columns = ['Gene_symbol', 'Entrez_ID', 'Average_Expression'])
print(data_analysis_2)
# Drops all data except the columns given above
matched_entrez = data_analysis_2[data_analysis_2['Gene_symbol'].isin(subset_genes_index)]
print(matched_entrez)
# Gets the entrez_ids from the randomly chosen 20 genes from before (to cross with the entrez database)

highest_entrez = data_analysis_2.nlargest(20, 'Average_Expression')
print(highest_entrez)
lowest_entrez = data_analysis_2.nsmallest(20, 'Average_Expression')
print(lowest_entrez)
# Get the entrez id from the 20 highest and lowest expressed genes









fig_1, axes = plt.subplots(1, 2, figsize=(15, 8))
# creates subplots arranged side by side of size 15'x 8'

sns.heatmap(subset_genes, cmap='viridis', ax=axes[0])
axes[0].set_title('Heatmap of Gene Expression Levels in Pancreatic Islets of 20 Random Genes')
axes[0].set_xlabel('Samples')
axes[0].set_ylabel('Gene')
# Creates seaborn heatmap of dataset subset_genes on the axis[0]; first plot/left plot
# Title and x/y labels given 

sns.heatmap(matched_rows, cmap='viridis', ax=axes[1])
axes[1].set_title('Heatmap of Average Gene Expression Levels in Pancreatic Islets of 20 Random Genes')
axes[1].set_xlabel('')
axes[1].set_ylabel('Gene')
# Creates seaborn heatmap of dataset matched_rows on the axis[1]; second plot/right plot
# x label not realy required as it takes from the dataset

plt.tight_layout()
plt.show()
fig_1.savefig('Gene Expression Heatmaps.png')
# .tight prevent overlapping between the two, shows the created plot and saves it in the same location as the 
# python file

fig_2, axes = plt.subplots(1, 2, figsize=(15, 8))


sns.heatmap(highest_expression, cmap='viridis', ax=axes[0])
axes[0].set_title('Heatmap of Highest Gene Expression Levels in Pancreatic Islet Cells')
axes[0].set_xlabel('Samples')
axes[0].set_ylabel('Gene')

heatmap = sns.heatmap(lowest_expression, cmap='viridis', ax=axes[1])
axes[1].set_title('Heatmap of Lowest Gene Expression Levels in Pancreatic Islet Cells')
axes[1].set_xlabel('')
axes[1].set_ylabel('Gene')
heatmap.invert_yaxis()

plt.tight_layout()
plt.show()
fig_2.savefig('Gene Expression Heatmaps_2.png')
# Same as above plots, but instead compares the highest and lowest expressed genes while also saving them in 
# same location





from Functions import get_gene_summary
from Functions import format_gene_summary

def process_gene_summaries(dataframe):
    summaries = []
    
    entrez_list = dataframe['Entrez_ID'].tolist()
    # Gets the entrez id values into a list 
    entrez_list = [value for value in entrez_list if pd.notna(value)]
    print(entrez_list)
    # Removes and NA values from the list 
    for value in entrez_list[:10]:
        # Gets first 10 entrez ids (0-9th) one at a time and carry out the following instructions
    #  if pd.notna(value):  # Check if value is not NA
        
            
        try:
            # Retrieve and format gene summary
            value = str(value)
            print(value)
            # Transforms the int64 to a string value at every iteration
            gene_summary = get_gene_summary(value)
            formatted_summary = format_gene_summary(gene_summary)
            # Gets the gene information and formats it (refer to the Functions python file)
            # Store the formatted summary in the empty list from before
            summaries.append(formatted_summary)
        except RuntimeError as e:
            print(f"No results found for Entrez ID {value}: {e}")
        
        # Introduce a delay(3 second) before fetching the next summary (get a lot of runtime errors otherwise)
        time.sleep(3) 
    return summaries






from Functions import write_gene_summaries_to_file
# Refer to Functions.py for more information 
try:
    summary_gene = process_gene_summaries(matched_entrez)
    filename_1 = "gene_summaries.txt"
    write_gene_summaries_to_file(summary_gene, filename_1)
    print("Summaries written to gene_summaries.txt")
except Exception as e:
    print(e)
    
time.sleep(15) # Introduces delay between requesting the different gene information (random/highest/lowest)

try:
    high_gene_summaries = process_gene_summaries(highest_entrez)
    filename_2 = "highest_expression.txt"
    write_gene_summaries_to_file(high_gene_summaries, filename_2)
    print("Summaries written to highest_expression.txt")
except Exception as e:
    print(f"An error occurred: {str(e)}")
    
time.sleep(15)

try:
    low_gene_summaries = process_gene_summaries(lowest_entrez)
    filename_3 = "lowest_expression.txt"
    write_gene_summaries_to_file(low_gene_summaries, filename_3)
    print("Summaries written to lowest_expression.txt")
except Exception as e:  
    print(f"An error of type {type(e).__name__} occurred: {str(e)}")









from Functions import open_files
open_files()


file_path_I = ["Gene Expression Heatmaps.png", "Gene Expression Heatmaps_2.png"] 

from Functions import open_images
open_images(file_path_I)
# Refer to the Functions file for more information

# This opens both the newly written text files and created images for veiwing
# Beware that the files (text/image) are modified everytime the different python files (Islets/Liver/Skeletal) are run
# No new text/image files are created to prevent congestion of generated results. 


"""




















