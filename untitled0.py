#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:23:02 2025

@author: tanu
"""
import requests

sample_ids = ['SAMEA4978670', 'SAMN05223155', 'SAMEA4649982', 'SAMN02604002']

base_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"

for sample_id in sample_ids:
    url = f"{base_url}?result=read_run&accession={sample_id}&fields=accession,run_alias,study_accession,sample_accession"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.text
        if data:
            print(f"Sample {sample_id} found in ENA!")
        else:
            print(f"Sample {sample_id} not found in ENA.")
    else:
        print(f"Failed to fetch data for {sample_id}. Status Code: {response.status_code}")
        

########
# gca
https://www.ebi.ac.uk/ena/portal/api/filereport?result=assembly&accession=SAMN02604002-H 'accept: */*
import requests
import csv
from io import StringIO
# List of sample IDs to check
sample_ids = ['SAMEA4978670', 'SAMN05223155', 'SAMEA4649982', 'SAMN02604002']


# Base URL for the ENA API endpoint targeting assembly information
base_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Iterate through each BioSample ID to fetch assembly details
for sample_id in sample_ids:
    # Construct the URL with the appropriate parameters for assembly data
    url = f"{base_url}?result=assembly&accession={sample_id}&fields=accession,assembly_title,description"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.text
        if data:
            print(f"Assembly details for Sample {sample_id}:\n{data}")
        else:
            print(f"No assembly data found for Sample {sample_id}.")
    else:
        print(f"Failed to fetch data for {sample_id}. Status Code: {response.status_code}")
###############
https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=status%3Dpublic&fields=accession,assembly_set_accession,assembly_name,assembly_title,tax_id,scientific_name,assembly_level,sample_accession,study_accession,study_name,strain,genome_representation,wgs_set,sequence_accession
https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=status%3Dpublic%20AND%20tax_eq(1313)&fields=accession,sample_accession

https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=status%3Dpublic&fields=accession,assembly_set_accession,assembly_name,assembly_title,tax_id,scientific_name,assembly_level,sample_accession,study_accession,study_name,strain

############
# ncbi checkm api

# GET: /genome/accession/{accessions}/dataset_report
#https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#
##########
curl -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_000179855.1/dataset_report" \
 -H 'accept: application/json' 


# List of assembly names or accessions
#assembly_names = ["GCA_000179855.1", "GCA_000007045.1"]
assembly_names = ["GCA_000007045.1"]

accession'
'current_accession'
'paired_accession'
'status'
'assembly_name'
'bioprojects': [{'accession'}]
'genome_notes'
    "checkm_info": {
        "checkm_marker_set": 
        "checkm_species_tax_id": 
        "checkm_marker_set_rank": 
        "checkm_version": 
        "completeness":
        "completeness_percentile": ''


import requests
import pandas as pd
from pandas import json_normalize

def fetch_assembly_data(accession):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
    headers = {'accept': 'application/json'}  # Ensure we specify JSON format
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()  # Parse JSON from response directly
    else:
        return None

# Example usage:
accession = "GCA_000179855.1"
data = fetch_assembly_data(accession)["reports"][0]
print(data)  # This will print the JSON response for the specified accession


#All the following as columns from below:
# data['checkm_info']
#         "checkm_marker_set": 
#         "checkm_species_tax_id": 
#         "checkm_marker_set_rank": 
#         "checkm_version": 
#         "completeness":
#         "completeness_percentile": ''


####
foozzz={
'gca': data['accession'],
'gcf' : data['assembly_info']['paired_assembly']['accession'],
'status' : data['assembly_info']['paired_assembly']['status'],
'biosample': data['assembly_info']['biosample']['accession'],
'bioproject': data['assembly_info']['biosample']['bioprojects'][0]['accession'],
'notes' : data['assembly_info']['genome_notes'],
'checkm_info': data['checkm_info']
}
df =  json_normalize(foozzz)
####

#read data for notgca acc
n = pd.read_csv("/home/pub/Work/data_arise_proteome/spneumo_dataset/notgca_acc.txt", header = None)
accessionsL = n[0].tolist()
accessionsL = accessionsL[0:5]


dfs = []

for accession in accessionsL:
    # Fetch the data for the current accession
    data = fetch_assembly_data(accession)["reports"][0]
    
    # Normalize and create DataFrame
    normalized_data = {
        'gca': data['accession'],
        'gcf': data['assembly_info']['paired_assembly']['accession'],
        'status': data['assembly_info']['paired_assembly']['status'],
        'biosample': data['assembly_info']['biosample']['accession'],
        'bioproject': data['assembly_info']['biosample']['bioprojects'][0]['accession'],
        'notes': data['assembly_info'].get('genome_notes', ''),  # Handling missing notes
        **data['checkm_info']  # Spread checkm_info fields directly into the dictionary
    }
    
    df = json_normalize(normalized_data)
    dfs.append(df)

# Concatenate all DataFrames into one
final_df = pd.concat(dfs, ignore_index=True)

print(final_df)

##########


# Prepare to collect DataFrame for each accession
dfs = []

for accession in accessionsL:
    # Fetch the data for the current accession
    report = fetch_assembly_data(accession)["reports"][0]
    
    # Normalize and create DataFrame
    normalized_data = {
        'gca': report.get('accession', pd.NA),
        'gcf': report.get('assembly_info', {}).get('paired_assembly', {}).get('accession', pd.NA),
        'status': report.get('assembly_info', {}).get('paired_assembly', {}).get('status', pd.NA),
        'biosample': report.get('assembly_info', {}).get('biosample', {}).get('accession', pd.NA),
        'bioproject': report.get('assembly_info', {}).get('biosample', {}).get('bioprojects', [{}])[0].get('accession', pd.NA),
        'notes': report.get('assembly_info', {}).get('genome_notes', pd.NA),
        **{key: report.get('checkm_info', {}).get(key, pd.NA) for key in ['checkm_marker_set', 
                                                                          'checkm_species_tax_id', 
                                                                          'checkm_marker_set_rank', 
                                                                          'checkm_version', 
                                                                          'completeness', 
                                                                          'completeness_percentile']}
    }
    
    df = json_normalize(normalized_data)
    dfs.append(df)

# Concatenate all DataFrames into one
final_df1 = pd.concat(dfs, ignore_index=True)

print(final_df1)
##############

import pandas as pd
from pandas import json_normalize
from tqdm import tqdm  # Assuming tqdm is already installed

# List of accession numbers
accessions = ["GCA_000179855.1", "GCA_000123456.1", "GCA_000654321.1"]

# Prepare to collect DataFrame for each accession
dfs = []

for accession in tqdm(accessionsL, desc="Processing accessions"):
    # Fetch the data for the current accession
    data = fetch_assembly_data(accession)
    reports = data.get('reports', [{}])  # Use .get() to provide a default empty dict if 'reports' is missing

    if reports:  # Check if there is at least one report
        report = reports[0]  # Take the first report

        # Normalize and create DataFrame
        normalized_data = {
            'gca': report.get('accession', pd.NA),
            'gcf': report.get('assembly_info', {}).get('paired_assembly', {}).get('accession', pd.NA),
            'status': report.get('assembly_info', {}).get('paired_assembly', {}).get('status', pd.NA),
            'biosample': report.get('assembly_info', {}).get('biosample', {}).get('accession', pd.NA),
            'bioproject': report.get('assembly_info', {}).get('biosample', {}).get('bioprojects', [{}])[0].get('accession', pd.NA),
            'notes': report.get('assembly_info', {}).get('genome_notes', pd.NA),
            **{key: report.get('checkm_info', {}).get(key, pd.NA) for key in [
                'checkm_marker_set', 
                'checkm_species_tax_id', 
                'checkm_marker_set_rank', 
                'checkm_version', 
                'completeness', 
                'completeness_percentile'
            ]}
        }
        
        df = json_normalize(normalized_data)
        dfs.append(df)
    else:
        print(f"No reports available for accession: {accession}")

# Concatenate all DataFrames into one
final_df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

print(final_df)
