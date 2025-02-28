#!/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from tqdm import 
pip install matplotlib
pip install matplotlib_venn
pip install seaborn
pip install seaborn[stats]
pip3 install --upgrade cmplot

import seaborn as sns


#input_file="/home/tunstall/Documents/arise/spneumo_dataset/outL.tsv"
input_file="/home/pub/Work/data_arise_proteome/spneumo_dataset/outL.tsv"
input_file ="/home/pub/Work/data_arise_proteome/spneumo_dataset/clusterizes_proteomes"
input_file = "/home/pub/Work/data_arise_proteome/spneumo_dataset/Labelled_Species_protein_cluster.tsv"
input_file ="/home/pub/Work/data_arise_proteome/spneumo_dataset/clusterizes_proteomes2"
input_file ="/home/pub/Work/data_arise_proteome/spneumo_dataset/clusterizes_proteomes_all"



chunksize=100000
#output_file="/home/tunstall/Documents/arise/spneumo_dataset/Count_outL.tsv"
#output_file="/home/tunstall/Documents/arise/spneumo_dataset/Count_Labelled_Species_protein_cluster.tsv"
output_file="/home/pub/Work/data_arise_proteome/spneumo_dataset/Count_Labelled_Species_protein_cluster.tsv"

min_proteomes=0

#########
df = pd.read_csv(input_file, sep = "\t")


a = df.groupby(['cluster_id',  'protein_id']).size()
b = df.groupby('cluster_id')['protein_id'].nunique()
c = df.groupby('protein_id')['cluster_id'].nunique()

#########


##########


def process_file(input_file, output_file, chunksize):
    """Processes the input file to count proteomes per cluster and saves the result."""
    cluster_counts = pd.Series(dtype=int)
    
    # Count total lines for progress tracking
    total_lines = sum(1 for _ in open(input_file, 'r'))
    progress_bar = tqdm(total=total_lines, desc="Processing", unit="lines")
    
    # Read in chunks and count occurrences
    for chunk in pd.read_csv(input_file, sep='\t', usecols=['cluster_id'], chunksize=chunksize):
        #cluster_counts = cluster_counts.add(chunk['cluster_id'].value_counts(), fill_value=0)
        cluster_counts = cluster_counts.add(chunk.groupby(['cluster_id', 'proteomes']).value_counts(), fill_value=0)

        progress_bar.update(len(chunk))
    
    progress_bar.close()
    
    # Save to file
    cluster_counts = cluster_counts.astype(int).reset_index()
    cluster_counts.columns = ['cluster_id', 'count_of_proteomes']
    cluster_counts.to_csv(output_file, sep='\t', index=False)
    print(f"\nCluster count saved to: {output_file}")

def plot_data(input_file, output_plot, min_proteomes=0):
    """Plots the distribution of cluster sizes, filtering by minimum proteomes if needed."""
    df = pd.read_csv(input_file, sep='\t')
    #df = pd.read_csv(output_file, sep='\t')

    # Filter based on min_proteomes threshold
    df2 = df[df['count_of_proteomes'] >= min_proteomes]

    plt.figure(figsize=(10, 6))
    plt.hist(df2['count_of_proteomes'], bins=30, edgecolor='black')

    plt.title('Distribution of Cluster Sizes')
    plt.xlabel('Cluster Size (Number of Proteomes)')
    plt.ylabel('Frequency')
    plt.grid(True)
    
    # Apply log scale to the axis
    #plt.yscale("log")
    #plt.xscale("log")
    
    plt.savefig(output_plot)
    print(f"Plot saved to: {output_plot}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process cluster files and generate plots.")
    parser.add_argument("input_file", help="Path to the input file (raw dataset or precomputed count file).")
    parser.add_argument("--process", action="store_true", help="Compute cluster proteome counts.")
    parser.add_argument("--plot", action="store_true", help="Generate cluster size distribution plot.")
    parser.add_argument("--output_file", default=None, help="Custom output file for counts (default: 'Count_<input_file>').")
    parser.add_argument("--output_plot", default=None, help="Custom output file for plot (default: 'Plot_<input_file>.png').")
    parser.add_argument("--chunksize", type=int, default=100000, help="Number of rows to process per chunk (default: 100000).")
    parser.add_argument("--min_proteomes", type=int, default=1, help="Minimum proteomes per cluster to include in plot (default: 1).")

    # Step 2: Generate plot if requested
    if args.plot:
        plot_input = args.output_file if args.process else args.input_file
        print(f"Generating plot from: {plot_input}")
        plot_data(plot_input, args.output_plot, args.min_proteomes)

###############################################################################:wq

###############################################################################
#input_file = "/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Labelled_Species_protein_cluster.tsv"
#output_file = "/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Count_Labelled_Species_protein_cluster.tsv"
#chunksize= 100000
#min_proteomes=0

#f = pd.read_csv("/home/pub/Work/data_arise_proteome/spneumo_dataset/test_ds/data/Labelled_Species_protein_cluster.tsv", sep = "\t")

#f.groupby(['cluster_id']).size()
#f.groupby(['protein_id']).size()
#f.groupby(['proteomes']).size()
#f.groupby(['cluster_id', 'proteomes'])
#f.value_counts(subset = ['cluster_id'])

#f2 = f['cluster_id'].value_counts()

######
people = ['Hannah', 'Bethany', 'Kris', 'Alex', 'Earl', 'Lori']
reputation = ['awesome', 'cool', 'brilliant', 'meh', 'awesome', 'cool']
dictionary = dict(zip(people, reputation))
df_eg = pd.DataFrame(dictionary.values(), dictionary.keys())
df_eg = df.rename(columns={0:'reputation'})
seaborn.countplot(x='reputation', data=df_eg)


##
input_file="/home/pub/Work/data_arise_proteome/spneumo_dataset/clusterizes_proteomes"
df2 = pd.read_csv(input_file_cs, header = None) #727307
df = df.rename(columns={0:'count_proteomes'})
df3 = df.groupby('count_proteomes').size()
df.groupby('count_proteomes').size().plot(kind='bar')

sns.countplot('count_proteomes', data=df)
sns.countplot(y='count_proteomes', data=df)

sns.countplot(data=df, x='count_proteomes', order=df.count_proteomes.value_counts().index)

fig,ax = plt.subplots(figsize=(10,16))
grouped=df.groupby('count_proteomes').size(). \
sort_values(ascending=False).plot(kind='barh',ax=ax)

sns.histplot(data=df, x='count_proteomes', bins = 5000)
plt.show()


###

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data from the TSV file
df = pd.read_csv(input_file, sep='\t')

# Group by 'cluster_id' and count unique 'protein_id' values
cluster_sizes = df.groupby('cluster_id')['protein_id'].nunique().reset_index()
cluster_sizes.rename(columns={'protein_id': 'num_proteins'}, inplace=True)

# Create the bar plot using Seaborn
plt.figure(figsize=(12, 6))  # Adjust figure size as needed
sns.barplot(x='cluster_id', y='num_proteins', data=cluster_sizes, palette='viridis') # You can change the pallete
plt.xlabel('Cluster ID')
plt.ylabel('Number of Unique Proteins')
plt.title('Cluster Sizes Based on Number of Unique Proteins')
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for readability
plt.tight_layout()  # Adjust layout to

## histplot

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data (consider using chunking if the file is too large to fit in memory)
df = pd.read_csv(input_file, sep='\t')

# Calculate cluster sizes
cluster_sizes = df.groupby('cluster_id')['proteomes'].nunique()

# Create a histogram of cluster sizes
plt.figure(figsize=(10, 6))
sns.histplot(cluster_sizes, bins=50, kde=True)  # Adjust the number of bins as needed
plt.xlabel('Cluster Size (Number of Unique Proteins)')
plt.ylabel('Number of Clusters')
plt.title('Distribution of Cluster Sizes')
plt.show()




#######################################################################

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm
import os

def process_file(input_file, output_file, chunksize):
    """Processes the input file to count unique proteomes per cluster and saves the result."""
    cluster_counts = {}

    total_lines = sum(1 for _ in open(input_file, 'r'))
    progress_bar = tqdm(total=total_lines, desc="Processing", unit="lines")

    for chunk in pd.read_csv(input_file, sep='\t', usecols=['cluster_id', 'proteomes'], chunksize=chunksize):
        for cluster_id, group in chunk.groupby('cluster_id'):
            unique_proteomes = group['proteomes'].nunique()
            if cluster_id in cluster_counts:
                cluster_counts[cluster_id] += unique_proteomes
            else:
                cluster_counts[cluster_id] = unique_proteomes

        progress_bar.update(len(chunk))

    progress_bar.close()

    # Convert the dictionary to a Pandas Series
    cluster_counts_series = pd.Series(cluster_counts)

    # Save the counts to a file
    cluster_counts_series.to_csv(output_file, sep='\t', header=True)  # Changed to TSV
    print(f"Cluster proteome counts saved to: {output_file}")

cluster_sizes = cluster_counts_series.copy()

#def plot_data(input_file, output_plot, min_proteomes):
def plot_data(output_file, output_plot, min_proteomes):

    """Generates a histogram of cluster sizes from a precomputed counts file."""
    try:
        #cluster_sizes = pd.read_csv(input_file, sep='\t', index_col=0, header=0)
        cluster_sizes = pd.read_csv(input_file, sep='\t', header=None)

        cluster_sizes1 = pd.read_csv(input_file, sep='\t', header=None)
        cluster_sizes2 = pd.read_csv(input_file, sep='\t', header=None)

        cluster_sizes = pd.read_csv(output_file, sep='\t', index_col=0, header=0)

    except Exception as e:
        #print(f"Error reading count file {input_file}: {e}")
        print(f"Error reading count file {output_file}: {e}")

        return

    # Apply the filtering
    cluster_sizes = cluster_sizes[cluster_sizes > min_proteomes]
    cluster_sizes=cluster_sizes.dropna()
    if len(cluster_sizes) == 0:
         print("No clusters left after filtering")
         return


    plt.figure(figsize=(10, 6))
    sns.histplot(cluster_sizes, bins=50, kde=True)  # Adjust bins as needed
    #sns.histplot(cluster_sizes1, bins=50, kde=True)  # Adjust bins as needed

    plt.xlabel('Cluster Size (Number of Unique Proteomes)')
    plt.ylabel('Number of Clusters')
    plt.title('Distribution of Cluster Sizes ()')
    plt.yscale('log')  # Consider log scale if appropriate
    plt.tight_layout()

    if output_plot:
        plt.savefig(output_plot)
        print(f"Cluster size distribution plot saved to: {output_plot}")
    else:
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process cluster files and generate plots.")
    parser.add_argument("input_file", help="Path to the input file (raw dataset or precomputed count file).")
    parser.add_argument("--process", action="store_true", help="Compute cluster proteome counts.")
    parser.add_argument("--plot", action="store_true", help="Generate cluster size distribution plot.")
    parser.add_argument("--output_file", default=None, help="Custom output file for counts (default: 'Count_<input_file>').")
    parser.add_argument("--output_plot", default=None, help="Custom output file for plot (default: 'Plot_<input_file>.png').")
    parser.add_argument("--chunksize", type=int, default=100000, help="Number of rows to process per chunk (default: 100000).")
    parser.add_argument("--min_proteomes", type=int, default=1, help="Minimum proteomes per cluster to include in plot (default: 1).")

    args = parser.parse_args()

    if args.process:
        output_file = args.output_file if args.output_file else f"Count_{os.path.basename(args.input_file)}"
        process_file(args.input_file, output_file, args.chunksize)

    if args.plot:
        plot_input = args.output_file if args.process and args.output_file else (f"Count_{os.path.basename(args.input_file)}" if args.process else args.input_file)
        output_plot = args.output_plot if args.output_plot else f"Plot_{os.path.basename(args.input_file)}.png"
        print(f"Generating plot from: {plot_input}")
        plot_data(plot_input, output_plot, args.min_proteomes)


#####
#%


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_data(input_file, output_plot, min_proteomes):
    """Generates a histogram of cluster sizes from a precomputed counts file,
    with the x-axis showing cluster size as a percentage of the maximum cluster size."""

    try:
        #cluster_sizes = pd.read_csv(input_file, sep='\t', index_col=0, header=0)
        cluster_sizes = pd.read_csv(input_file, sep='\t', header=None)

    except Exception as e:
        print(f"Error reading count file {input_file}: {e}")
        return

    # Apply the filtering (before calculating max size)
    cluster_sizes = cluster_sizes[cluster_sizes > min_proteomes]
    cluster_sizes = cluster_sizes.dropna()
    if len(cluster_sizes) == 0:
        print("No clusters left after filtering")
        return

    # Calculate the maximum cluster size
    max_cluster_size = cluster_sizes.max()

    # Calculate the percentage of maximum cluster size for each cluster
    cluster_sizes_percent = (cluster_sizes / max_cluster_size) * 100

    # Define bin edges for 5% intervals
    bins = list(range(0, 105, 5))

    # Create the histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(cluster_sizes_percent, bins=bins, kde=True)

    plt.xlabel('Cluster Size (Percentage of Maximum Cluster Size)')
    plt.ylabel('Number of Clusters')
    plt.title('Distribution of Cluster Sizes (as Percentage of Maximum Size)')
    plt.xticks(bins)  # Set x-axis ticks to bin edges
    plt.yscale('log')
    plt.tight_layout()

    if output_plot:
        plt.savefig(output_plot)
        print(f"Cluster size distribution plot saved to: {output_plot}")
    else:
        plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_data(input_file, output_plot, min_proteomes):
    """Generates a histogram of cluster sizes from a precomputed counts file,
    with the x-axis showing cluster size as a percentage of the maximum cluster size,
    and prints descriptive statistics."""

    try:
        cluster_sizes = pd.read_csv(input_file, sep='\t', index_col=0, header=0, squeeze=True)
    except Exception as e:
        print(f"Error reading count file {input_file}: {e}")
        return

    # Apply the filtering (before calculating max size)
    cluster_sizes = cluster_sizes[cluster_sizes >= min_proteomes]
    cluster_sizes = cluster_sizes.dropna()
    if len(cluster_sizes) == 0:
        print("No clusters left after filtering")
        return

    # Calculate the maximum cluster size
    max_cluster_size = cluster_sizes.max()

    # Calculate the percentage of maximum cluster size for each cluster
    cluster_sizes_percent = (cluster_sizes / max_cluster_size) * 100

    # Define bin edges for 5% intervals
    bins = list(range(0, 105, 5))

    # Create the histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(cluster_sizes_percent, bins=bins, kde=True)

    plt.xlabel('Cluster Size (Percentage of Maximum Cluster Size)')
    plt.ylabel('Number of Clusters')
    plt.title('Distribution of Cluster Sizes (as Percentage of Maximum Size)')
    plt.xticks(bins)  # Set x-axis ticks to bin edges
    plt.yscale('log')
    plt.tight_layout()

    if output_plot:
        plt.savefig(output_plot)
        print(f"Cluster size distribution plot saved to: {output_plot}")
    else:
        plt.show()

    # Calculate and print descriptive statistics
    print("Descriptive Statistics:")
    print(f"  Minimum Cluster Size: {cluster_sizes.min()}")
    print(f"  Maximum Cluster Size: {cluster_sizes.max()}")
    print(f"  Mean Cluster Size: {cluster_sizes.mean()}")
    print(f"  Median Cluster Size: {cluster_sizes.median()}")
    print(f"  Standard Deviation: {cluster_sizes.std()}")
    print(f"  Number of Clusters: {len(cluster_sizes)}")

###############################################################################
#TODO:
import pandas as pd
import subprocess
import os

def process_cluster_data_cmd(input_file, output_file):
    """
    Processes cluster data using a command-line pipeline, mimicking the following shell command:
    cut -f 1,3 input_file | sort -u | cut -f 1 | sort | uniq -c | awk '{print $1}' | sort -n > output_file

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output file to store the processed counts.
    """

    try:
        # Construct the command-line pipeline
        cmd = f"""
        cut -f 1,3 {input_file} | 
        sort -u | 
        cut -f 1 | 
        sort | 
        uniq -c | 
        awk '{{print $1}}' | 
        sort -n
        """

        # Execute the command using subprocess
        result = subprocess.run(cmd, shell=True, text=True, capture_output=True, check=True)

        # Write the output to the specified file
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        print(f"Processed data saved to {output_file}")

    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print(f"Stderr: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def process_cluster_data_pandas(input_file, output_file):
    """
    Processes cluster data using pandas to mimic the following shell command:
    cut -f 1,3 input_file | sort -u | cut -f 1 | sort | uniq -c | awk '{print $1}' | sort -n > output_file
    
    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output file to store the processed counts.
    """
    try:
        # Read the data, selecting and naming the columns
        df = pd.read_csv(input_file, sep='\t', usecols=[0, 2], names=['cluster_id', 'proteomes'], header=None)

        # Create unique combinations of cluster_id and proteomes
        df['combined'] = df['cluster_id'].astype(str) + '_' + df['proteomes'].astype(str)
        unique_combinations = df['combined'].unique()

        # Extract cluster_id from unique combinations
        cluster_ids = [combo.split('_')[0] for combo in unique_combinations]

        # Count occurrences of each cluster_id
        cluster_counts = pd.Series(cluster_ids).value_counts().sort_index()

        # Save the counts to a file
        cluster_counts.to_csv(output_file, sep='\t', header=False)
        
        print(f"Processed data saved to {output_file}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Example usage:
    input_file = "outL.tsv"  # Replace with your input file path
    output_file = "clusterizes_proteomes2"  # Replace with your desired output file path
    
    # Call the function
    process_cluster_data_pandas(input_file, output_file)  # You can choose either Pandas or CMD

    #Verify to make sure they both work properly
    process_cluster_data_cmd(input_file, "clusterizes_proteomes2_cmd")
