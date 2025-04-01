#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:23:02 2025

@author: tanu
"""
import sys, os
import requests
import pandas as pd
from pandas import json_normalize
from tqdm import tqdm
###############################################################################
# Global
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"
###############################################################################
# ncbi checkm api

# GET: /genome/accession/{accessions}/dataset_report
#https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#
# GET /genome/accession/{accessions}/dataset_report
##########
#curl -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_000179855.1/dataset_report" -H 'accept: application/json' 

# List of assembly names or accessions
#assembly_names = ["GCA_000179855.1", "GCA_000007045.1"]
#assembly_names = ["GCA_000007045.1"]

# 'accession'
# 'status'
# 'assembly_name'
# 'bioprojects': [{'accession'}]
# 'genome_notes'
# 'checkm_info': {
#         "checkm_marker_set": 
#         "checkm_species_tax_id": 
#         "checkm_marker_set_rank": 
#         "checkm_version": 
#         "completeness":
#         "completeness_percentile":}

def fetch_assembly_data(accession):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
    headers = {'accept': 'application/json'}  # Ensure we specify JSON format
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()  # Parse JSON from response directly
    else:
        return None

# Example usage: one GCA
# accession = ["GCA_000179855.1", "GCA_001134245.1"]
# for acc in accession:
#    data = fetch_assembly_data(acc)["reports"][0]
#    print(data)  # This will print the JSON response for the specified accession

# # does not work if the key does not exist for any accession:
# egD={
# 'gca': data['accession'],
# 'gcf' : data['assembly_info']['paired_assembly']['accession'],
# 'status' : data['assembly_info']['paired_assembly']['status'],
# 'biosample': data['assembly_info']['biosample']['accession'],
# #'bioproject': data['assembly_info']['biosample']['bioprojects'][0]['accession'],
# 'bioproject': data['assembly_info']['bioproject_accession'],
# 'notes' : data['assembly_info']['genome_notes'],
# 'checkm_contamination': data['checkm_info']['contamination'],
# 'checkm_info': data['checkm_info']
# }
# df_eg =  json_normalize(egD)


###############################################################################
#============
# Spneumo GCA with no checkM from file: CheckM_report_prokaryotes.txt (directly from ncbi)
#============
# List of accession numbers
#accessionsL = ["GCA_000179855.1", "GCA_000123456.1", "GCA_000654321.1"]
#accessionsL = ['GCA_000179855.1', 'GCA_001134245.1', 'GCA_000983055.1', 
#               'GCA_000983075.2',  'GCA_001083385.1', 'GCA_001083445.1']
#accessionsL = ["GCA_000179855.1"] # website contamination 0, but being returned as nan

# List of accessions: from file notgca acc.txt
infile_no_checkm = basedir + "/gca_acc_NOcheckM.txt"
print(f"\nReading Input file: {infile_no_checkm}")
no_checkm = pd.read_csv(infile_no_checkm, header = None)
accessionsL = no_checkm[0].tolist()
#accessionsL = accessionsL[0:5]
print(f"\nNo. of s.pneumo GCA ids without checkM from NCBI flat file: {len(accessionsL)}")

# output file
out_name = os.path.splitext(os.path.basename(infile_no_checkm))[0] + "_checkM_API.txt"
outfile_no_checkm = basedir + "/" + out_name


# Prepare to make a DataFrame for each accession, with the selected information from the json 
# If the key doesn't exist, it is created with value NA
# If the GCA itself doesn't exist, i.e reports is empty, create it blank
dfs = []


for accession in tqdm(accessionsL, desc="Processing accessions"):
    # Fetch the data for the current accession
    data = fetch_assembly_data(accession)
    reports = data.get('reports', [{}])  # Use .get() to provide a default empty dict if 'reports' is missing

    if reports:  # Check if there is at least one report
        report = reports[0]  # Take the first report

        # Normalize and create DataFrame
        selected_data = {
            'gca': report.get('accession', pd.NA),
            'gcf': report.get('assembly_info', {}).get('paired_assembly', {}).get('accession', pd.NA),
            'status': report.get('assembly_info', {}).get('paired_assembly', {}).get('status', pd.NA),
            'biosample': report.get('assembly_info', {}).get('biosample', {}).get('accession', pd.NA),
            #'bioproject': report.get('assembly_info', {}).get('biosample', {}).get('bioprojects', [{}])[0].get('accession', pd.NA),
            'bioproject': report.get('assembly_info', {}).get('bioproject_accession', pd.NA),
            'notes': report.get('assembly_info', {}).get('genome_notes', pd.NA),
            **{key: report.get('checkm_info', {}).get(key, pd.NA) for key in [
                'checkm_marker_set', 
                'checkm_species_tax_id', 
                'checkm_marker_set_rank', 
                'checkm_version', 
                'completeness', 
                'contamination',
                'completeness_percentile'
            ]}
        }
        
        df = json_normalize(selected_data)
        dfs.append(df)
    else:
        print(f"No reports available for accession: {accession}")

# Concatenate all DataFrames into one
final_df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

print(final_df)
print(f"\nWriting output file: {outfile_no_checkm}")
print(f"\nOutfile shape: {final_df.shape}")

final_df.to_csv(outfile_no_checkm, index=False)

###############################################################################
#============
# Spneumo GCA with checkM from file: CheckM_report_prokaryotes.txt (directly from ncbi)
#============
# List of accessions: from file notgca acc.txt
infile_checkm = basedir + "/gca_acc_checkM.txt"
print(f"\nReading Input file: {infile_checkm}")
checkm = pd.read_csv(infile_checkm, header = None)
accessionsL = checkm[0].tolist()
#accessionsL = accessionsL[0:5]
print(f"\nNo. of s.pneumo GCA ids with checkM from NCBI flat file: {len(accessionsL)}")

# output file
out_name = os.path.splitext(os.path.basename(infile_checkm))[0] + "_checkM_API.txt"
outfile_checkm = basedir + "/" + out_name

# Prepare to make a DataFrame for each accession, with the selected information from the json 
# If the key doesn't exist, it is created with value NA
# If the GCA itself doesn't exist, i.e reports is empty, create it blank
dfs = []

for accession in tqdm(accessionsL, desc="Processing accessions"):
    # Fetch the data for the current accession
    data = fetch_assembly_data(accession)
    reports = data.get('reports', [{}])  # Use .get() to provide a default empty dict if 'reports' is missing

    if reports:  # Check if there is at least one report
        report = reports[0]  # Take the first report

        # create DataFrame
        selected_data = {
            'gca': report.get('accession', pd.NA),
            'gcf': report.get('assembly_info', {}).get('paired_assembly', {}).get('accession', pd.NA),
            'status': report.get('assembly_info', {}).get('paired_assembly', {}).get('status', pd.NA),
            'biosample': report.get('assembly_info', {}).get('biosample', {}).get('accession', pd.NA),
#            'bioproject': report.get('assembly_info', {}).get('biosample', {}).get('bioprojects', [{}])[0].get('accession', pd.NA),
            'bioproject': report.get('assembly_info', {}).get('bioproject_accession', pd.NA),

            'notes': report.get('assembly_info', {}).get('genome_notes', pd.NA),
            **{key: report.get('checkm_info', {}).get(key, pd.NA) for key in [
                'checkm_marker_set', 
                'checkm_species_tax_id', 
                'checkm_marker_set_rank', 
                'checkm_version', 
                'completeness', 
                'contamination',
                'completeness_percentile'
            ]}
        }
        
        df = json_normalize(selected_data)
        dfs.append(df)
    else:
        print(f"No reports available for accession: {accession}")

# Concatenate all DataFrames into one
final_df2 = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

print(final_df2)
print(f"\nWriting output file: {outfile_checkm}")
print(f"\nOutfile shape: {final_df2.shape}")

final_df2.to_csv(outfile_checkm, index=False)

##
final_df2.groupby(['notes']).value_counts()

