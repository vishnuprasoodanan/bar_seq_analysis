#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:40:56 2024

@author: vishnu
"""
import os
import pandas as pd

# Specify the folder path
folder_path = '/Users/vishnu/WORK/Fabiola_Barseq_Analysis/'

# List all files in the folder
file_list = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# Loop through each CSV file and print its content
for file_name in file_list:
    file_path = os.path.join(folder_path, file_name)
    fh = open(file_path, 'r', encoding = 'utf-8')
 #fh = open('/Users/vishnu/WORK/Fabiola_Barseq_Analysis/uni_KO_Analysis/uni_KO.csv', 'r', encoding = 'utf-8')
    df = pd.read_csv(fh)
    
    print(df.columns)
    
    excluded_columns = [0, 2, 3]
    subset_df = df.drop(df.columns[excluded_columns], axis=1)
    
    sum_by_column2 = subset_df.groupby('final_annotation').sum()
    
    print("Original DataFrame:")
    print(df)
    
    print("\nSubset DataFrame:")
    print(subset_df)
    
    print("\nSum by 'final_annotation':")
    print(sum_by_column2)
    
    # Create an empty DataFrame
    empty_df = pd.DataFrame()
    
    # Vector of names
    names = sum_by_column2.columns
    
    # Add a column to the empty DataFrame with the vector of names
    empty_df['Names'] = names
    
    # Display the resulting DataFrame
    print(empty_df)
    
    empty_df['Status'] = empty_df['Names'].apply(lambda x: x.split('_')[0])
    
    print(empty_df)
    
    output_file_name1 = file_name.replace('.csv', '_gene_count.txt')
    output_file_path1 = os.path.join(folder_path, output_file_name1)
    
    # Write the content to the output file
    sum_by_column2.to_csv(output_file_path1, sep='\t')
    print(f'Content written to output file {output_file_name1}\n' + '='*30 + '\n')
    
    # Construct the output file name
    output_file_name2 = file_name.replace('.csv', '_metadata.txt')
    output_file_path2 = os.path.join(folder_path, output_file_name2)
    
    # Write the content to the output file
    empty_df.to_csv(output_file_path2, sep='\t', index=False)
    print(f'Content written to output file {output_file_name2}\n' + '='*30 + '\n')
