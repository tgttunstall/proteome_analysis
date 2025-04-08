import os, sys
import pandas as pd

###############################################################################
###############################################################################
# Read data
homedir = os.path.expanduser("~")
#basedir =  homedir + "/Documents/arise/spneumo_dataset"
basedir =  "/home/pub/Work/data_arise_proteome/spneumo_dataset"

#=====
# ATB
#=====
input_file_atb_paths  = basedir+"/atb_proteome_paths_TT.txt"
input_file_ena_gca =  basedir+"/ena_gca_biosample.txt" # try to get  proteincount as well..

outfile_atb_checkm = basedir+"/input_list_atb_checkm.txt"


#=====
# UP
#=====
input_file_up  = basedir+"/spneumo_biosample_info.out"
up_file_path_toadd = '/nfs/research/martin/uniprot/research/spneumo_dataset/up_data/'
#up_file_path_toadd = '/home/pub/Work/data_arise_proteome/spneumo_dataset'
input_file_up_pid =  basedir+"/spneumo_proteomeids.out" 

outfile_up_checkm = basedir+"/input_list_up_checkm.txt" 

###############################################################################
#=====
# ATB Processing
#=====
#Step 1: Extract Biosample IDs and Paths from File
# Read the file into a DataFrame
df = pd.read_csv(input_file_atb_paths, header=None, names=['path'])

# Extract the biosample ID from the path
df['biosample'] = df['path'].apply(lambda x: x.split('/')[-1].split('.')[0])

print(df.head())


#Step 2: Load Genome Assembly Information
# File containing genome assembly info
# Assuming this file is tab-separated and has columns 'BiosampleID' and 'GenomeAssembly'
ena_df = pd.read_csv(input_file_ena_gca, sep='\t') #3020015,14

# subset columns
genome_df = ena_df[['assembly_set_accession', 
                    #'assembly_name',
                    'tax_id', 
                    'scientific_name', 
                    'assembly_level',
                    'sample_accession', 
                    'study_accession', 
                    'genome_representation']] #3020015,7

# rename cols
genome_df.rename(columns = {'assembly_set_accession': 'gc_set_acc',
                  #'assembly_name',
                  #tax_id': 'taxid',
                  'scientific_name': 'species_name',
                  'assembly_level':'assembly_level_assembled' ,
                  'sample_accession': 'biosample',
                  'study_accession': 'project_acc',
                  'genome_representation':'assembly_level'
                  }, inplace=True)

print(genome_df.head())
a = genome_df.head()
print(genome_df.columns)

#Step 3: Merge DataFrames
merged_df = pd.merge(df, genome_df, on='biosample', how='left')
print(merged_df.head())

# count na
c1=merged_df['gc_set_acc'].isna().sum()
merged_df['gc_set_acc'].isna().sum()

#Step 3a: Fill Na with 'GA'
#merged_df['gc_set_acc'].fillna('GA', inplace=True)
#For example, when doing 'df[col].method(value, inplace=True)', 
# try using 'df.method({col: value}, inplace=True)' or df[col] = df[col].method(value) instead, to perform the operation inplace on the original object.
merged_df.fillna({'gc_set_acc': 'GA'}, inplace=True)
merged_df.isna().sum()
merged_df['gc_set_acc'].value_counts()
(merged_df[['gc_set_acc']] == 'GA').sum() == c1

#Step 4: Add value counts as a column
b_counts = merged_df['biosample'].value_counts()
merged_df['biosample_count'] = merged_df['biosample'].map(b_counts)

#Step 5: Sort the df
merged_dfS = merged_df.sort_values(['biosample_count', 'biosample', 'gc_set_acc'], ascending=[False, True, False])

#Step 6: Drop duplicates: keeps the first (I have ordered by latest GA)
merged_dfSU  = merged_dfS.drop_duplicates(['biosample'])

#Step 7 Subset the cols we need for output
merged_dfSU_output = merged_dfSU[['biosample','gc_set_acc','path']]

#Step 8: Save the Results
merged_dfSU_output.to_csv(outfile_atb_checkm, header=False, index=False, sep="\t")
print(f"Input list for checkm ATB dataset: {outfile_atb_checkm}")

###############################################################################
#=====
# UP Processing
#=====
#Step 1: Load UP Data
up_all_df = pd.read_csv(input_file_up, sep='\t') #27622, 20
up_all_df['biosample'].value_counts()
up_all_df.nunique()
up_all_df.columns

#Step 2 : Load PID Information
# File containing genome assembly info
# Assuming this file is tab-separated and has columns 'BiosampleID' and 'GenomeAssembly'
up_pid_df = pd.read_csv(input_file_up_pid, sep='\t', header=None) #26748,4

up_pid_df.columns=['pid', 'upid', 'gca_set_acc', 'species_name']

up_pid_df.columns
up_pid_df.nunique()
up_pid_df['species_name'].value_counts()

#Step 3: Merge DataFrames
merged_df_up = pd.merge(up_all_df, up_pid_df, on='upid', how='left') #27622,23
print(merged_df_up.head())

# count na
c2=merged_df_up['gc_set_acc'].isna().sum()
merged_df_up['gc_set_acc'].isna().sum()

#Step 3a: Fill Na with 'GA'
merged_df_up.fillna({'gc_set_acc': 'GA'}, inplace=True)
merged_df_up.isna().sum()
merged_df_up['gc_set_acc'].value_counts()
(merged_df_up[['gc_set_acc']] == 'GA').sum() == c2

#Step 4: Drop duplicates UPID
merged_df_upU  = merged_df_up.drop_duplicates(['upid']) #26749,23

#Step 5: Add value counts as a column
b_counts = merged_df_upU['biosample'].value_counts()
merged_df_upU['biosample_count'] = merged_df_upU['biosample'].map(b_counts)

#Step 6: Sort the df
merged_df_upUS = merged_df_upU.sort_values(['biosample_count', 'biosample', 'gc_set_acc'], ascending=[False, True, False])

#Step 7: Drop duplicates UPID again: sanity check
merged_df_upUS_1  = merged_df_upUS.drop_duplicates(['upid']) #26749,23
merged_df_upUS_1.equals(merged_df_upUS)

#Step 8: Add  path
merged_df_upUS['path'] = up_file_path_toadd + 'proteome_' + merged_df_upUS['pid'].astype(str) + '.fa'

#Step 9 Subset the cols we need for output
merged_dfUS_output = merged_df_upUS[['pid','gc_set_acc','path']]

#Step 8: Save the Results
merged_dfUS_output.to_csv(outfile_up_checkm, header=False, index=False, sep="\t")
print(f"Input list for UP dataset: {outfile_up_checkm}")