'''
This code is a part of DraftMet
Author : Manjnatha B K
Date : 6th Feb, 2024

'''

import requests
from bs4 import BeautifulSoup
import re
import pandas as pd
from tqdm import tqdm

#take input from user
org_name = input(str("Please enter KEGG organism code: "))
organism_id = input(str("Please enter UNIPROT oragnism ID: "))

URL = f"https://rest.kegg.jp/list/pathway/{org_name}"

r = requests.get(URL)

soup = BeautifulSoup(r.content, 'html.parser')


html_pathway_string = str(soup)

# Split the string into lines
pathway_lines = html_pathway_string.strip().split('\n')

# Initialize lists to store pathway code and description
pathway_codes = []
pathway_descriptions = []

# Iterate over pathway lines to extract pathway code and description
for line in pathway_lines:
    # Split each line into code and description
    code, description = line.split('\t')
    pathway_codes.append(code.strip())
    pathway_descriptions.append(description.strip())

# Create DataFrame
df = pd.DataFrame({'Pathway_Code': pathway_codes, 'Pathway_Description': pathway_descriptions})

ec_df = pd.DataFrame(columns=['EC Numbers', 'Pathway_Code'])

process_bar = tqdm(total= len(pathway_codes), desc="processing")

for i in range(len(pathway_codes)):
    URL = "https://www.kegg.jp/entry/"+pathway_codes[i]+"/"
    r = requests.get(URL)
    soup = BeautifulSoup(r.content, 'html.parser')
    html_string = str(soup)
    regex = re.compile(r'\[EC:.*?\]')
    matches = regex.findall(html_string)
    # Iterate through each entry and apply the regex pattern
    for entry_string in matches:
        extracted_content = re.findall(r'/entry/([^"]+)', entry_string)
        # Append all extracted content to the DataFrame
        for content in extracted_content:
            new_row = pd.DataFrame({'EC Numbers': [content], 'Pathway_Code': [pathway_codes[i]]})
            ec_df = pd.concat([ec_df, new_row], ignore_index=True)
    process_bar.update(1)

process_bar.close()

main_df = pd.merge(ec_df, df, how = 'inner')
#main_df.to_csv(f"{org_name}_draft.csv",index=False)
#print(f"Output file is saved as {org_name}_draft.csv")

ECList =main_df['EC Numbers'].tolist()

df1 = []
uni_data = []
process_bar2 = tqdm(total= len(ECList), desc="Processing UNIPROT Entries")

# Loop through each EC number
for ec_num in ECList:
    URL = f"https://rest.uniprot.org/uniprotkb/search?query=(EC:{ec_num})%20AND%20(organism_id:{organism_id})&format=tsv"
    #print(URL)
    s = requests.get(URL)
    df1.append(str(s.text))
    process_bar2.update(1)

process_bar2.close()

# Iterate through data and extract UniProt entries
for i, ec_num in zip(df1, ECList):
    for line in i.rstrip().split("\n"):
        data1 = line.split("\t")
        # Append the EC number as the first column
        uni_data.append([ec_num] + data1)
        dataframe = pd.DataFrame(uni_data)

# Filter the rows that contain the substring "Entry"
substring = 'Entry'
filter = dataframe[1].str.contains(substring)
filtered_df = dataframe[~filter]

 #Rename the columns
filtered_df.columns = ['EC Numbers', 'Uniprot Entry', 'Uniprot ID', 'Reviewed', 'Enzyme name', 'Gene ID', 'Organism', 'Sequence length']  # Adjust column names as per your data

new_df = filtered_df.drop_duplicates(subset=['Uniprot Entry'])

# Create an empty dataframe to store the merged data
merged_rows = []

# Iterate through unique EC numbers
for ec_number in new_df['EC Numbers'].unique():
    # Match rows with the same EC number in both dataframes
    new_df_subset = new_df[new_df['EC Numbers'] == ec_number]
    main_df_subset = main_df[main_df['EC Numbers'] == ec_number]
    
    # Concatenate values for each column
    merged_row = {
        'EC Numbers': ec_number,
        'Uniprot Entry': ', '.join(new_df_subset['Uniprot Entry'].tolist()),
        'Uniprot ID': ', '.join(new_df_subset['Uniprot ID'].tolist()),
        'Reviewed': ', '.join(new_df_subset['Reviewed'].tolist()),
        'Enzyme name': ', '.join(new_df_subset['Enzyme name'].tolist()),
        'Gene ID': ', '.join(new_df_subset['Gene ID'].tolist()),
        'Organism': ', '.join(new_df_subset['Organism'].tolist()),
        'Sequence length': ', '.join(new_df_subset['Sequence length'].astype(str).tolist()),
    }
    
    # Append the merged row to the list
    merged_rows.append(merged_row)

# Convert the list of dictionaries into a DataFrame
merged_df = pd.DataFrame(merged_rows)

test_df = pd.merge(merged_df, main_df, how = 'right')

test_df.to_csv(f'{org_name}.csv')

print(f'DraftMet file was sucessfully saved as {org_name}.csv')


