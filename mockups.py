#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:31:16 2024

@author: tanu
"""
import matplotlib.pyplot as plt
import pandas as pd



import matplotlib.pyplot as plt
import numpy as np

# Sample data
RPs = ['RP1', 'RP2', 'RP3']
unique_proteins = [20, 15, 25]
common_proteins = [50, 60, 55]
fragments = [30, 25, 20]

ind = np.arange(len(RPs))  # the x locations for the groups
width = 0.5  # the width of the bars

fig, ax = plt.subplots()
p1 = ax.bar(ind, unique_proteins, width, label='Unique Proteins')
p2 = ax.bar(ind, common_proteins, width, bottom=unique_proteins, label='Common Proteins')
p3 = ax.bar(ind, fragments, width, bottom=[i+j for i,j in zip(unique_proteins, common_proteins)], label='Fragments')

ax.axhline(0, color='grey', linewidth=0.8)
ax.set_ylabel('Percentage of Protein Diversity')
ax.set_title('Protein Diversity Coverage by Reference Proteomes')
ax.set_xticks(ind)
ax.set_xticklabels(RPs)
ax.legend()

plt.show()
#########

import matplotlib.pyplot as plt

# Sample data
RPs = ['RP1', 'RP2', 'RP3', 'RP4', 'RP5']
# These values should be cumulative coverage percentages
protein_coverage = [20, 45, 65, 80, 90]  # Example cumulative percentages

plt.figure(figsize=(10, 5))
plt.plot(RPs, protein_coverage, marker='o', linestyle='-', color='b')

# Adding titles and labels
plt.title('Cumulative Protein Diversity Coverage by Reference Proteomes')
plt.xlabel('Number of Reference Proteomes Considered')
plt.ylabel('Cumulative Percentage of Protein Diversity Covered')
plt.grid(True)

# Highlighting the inflection point
plt.axvline(x=3, color='r', linestyle='--')  # Assuming RP3 is the inflection point
plt.axhline(y=65, color='r', linestyle='--')  # Corresponding y value to the inflection point

plt.show()
import matplotlib.pyplot as plt

# Total identified proteins from ggCaller
total_identified_proteins = 9400  # Hypothetical total count for simplicity

# Hypothetical unique proteins covered by each RP
unique_proteins_by_rp = [4000, 3000, 2000, 400, 100]  # Each RP covers fewer new proteins

# Calculate cumulative percentage of diversity covered
cumulative_coverage = []
current_coverage = 0
for proteins in unique_proteins_by_rp:
    current_coverage += proteins
    cumulative_coverage.append(current_coverage / total_identified_proteins * 100)

# X-axis labels (number of RPs considered)
rps_considered = range(1, len(unique_proteins_by_rp) + 1)

# Plotting the line graph
plt.figure(figsize=(10, 6))
plt.plot(rps_considered, cumulative_coverage, marker='o', linestyle='-', color='b')
plt.title('Cumulative Protein Diversity Coverage by Number of Reference Proteomes')
plt.xlabel('Number of Reference Proteomes Considered')
plt.ylabel('Cumulative Percentage of Protein Diversity Covered (%)')
plt.grid(True)
plt.axhline(y=100, color='r', linestyle='--', label='100% Coverage')
plt.legend()
plt.show()
##############













import matplotlib.pyplot as plt

# Hypothetical total identified proteins from ggCaller
total_identified_proteins = 9400  # Simplified total count for demonstration

# Hypothetical unique proteins covered by each RP
unique_proteins_by_rp = [4000, 3000, 2000, 400, 100]  # Each subsequent RP covers fewer new proteins

# Calculate cumulative percentage of diversity covered
cumulative_coverage = []
current_coverage = 0
for proteins in unique_proteins_by_rp:
    current_coverage += proteins
    cumulative_coverage.append(min(100, current_coverage / total_identified_proteins * 100))  # Cap at 100%

# X-axis labels (number of RPs considered)
rps_considered = range(1, len(unique_proteins_by_rp) + 1)

# First graph: Bar graph with individual RP contributions
plt.figure(figsize=(10, 5))
plt.bar(rps_considered, [x / total_identified_proteins * 100 for x in unique_proteins_by_rp], color='skyblue')
plt.title('Protein Diversity Covered by Individual Reference Proteomes')
plt.xlabel('Reference Proteomes')
plt.ylabel('Percentage of Protein Diversity Covered (%)')
plt.ylim(0, 100)  # Ensure y-axis goes up to 100%
plt.xticks(rps_considered)
plt.grid(True)
plt.show()

# Second graph: Line graph for cumulative diversity
plt.figure(figsize=(10, 5))
plt.plot(rps_considered, cumulative_coverage, marker='o', linestyle='-', color='b')
plt.title('Cumulative Protein Diversity Coverage by Reference Proteomes')
plt.xlabel('Number of Reference Proteomes Considered')
plt.ylabel('Cumulative Percentage of Protein Diversity Covered (%)')
plt.ylim(0, 100)  # Ensure y-axis goes up to 100%
plt.xticks(rps_considered)
plt.axhline(y=100, color='red', linestyle='--', label='100% Coverage')
plt.grid(True)
plt.legend()
plt.show()

########


import matplotlib.pyplot as plt

# Sample data for RP1 to RP5
RPs = ['RP1', 'RP2', 'RP3', 'RP4', 'RP5']
unique_proteins = [4000, 3000, 2000, 400, 100]  # Counts of unique proteins for each RP
common_proteins = [2000, 1500, 1000, 300, 50]   # Counts of common proteins for each RP
fragments = [500, 300, 200, 100, 25]            # Counts of fragments for each RP

# Calculate cumulative coverage in terms of actual counts, not percentage
cumulative_coverage = []
current_coverage = 0
for i in range(len(RPs)):
    current_coverage += unique_proteins[i] + common_proteins[i] + fragments[i]
    cumulative_coverage.append(current_coverage)  # Absolute cumulative protein count

# Normalize cumulative coverage for the line graph to display as a percentage
total_proteins = cumulative_coverage[-1]  # Total proteins covered by all RPs
cumulative_percent = [x / total_proteins * 100 for x in cumulative_coverage]  # Convert to percentage

# X-axis labels (number of RPs considered)
rps_considered = range(1, len(unique_proteins) + 1)

# First graph: Bar graph with stacked components
plt.figure(figsize=(10, 5))
bar_width = 0.85
plt.bar(RPs, unique_proteins, color='blue', edgecolor='white', width=bar_width, label='Unique Proteins')
plt.bar(RPs, common_proteins, bottom=unique_proteins, color='green', edgecolor='white', width=bar_width, label='Common Proteins')
plt.bar(RPs, fragments, bottom=[u+c for u, c in zip(unique_proteins, common_proteins)], color='red', edgecolor='white', width=bar_width, label='Fragments')

plt.xlabel('Reference Proteomes')
plt.ylabel('Count of Proteins')
plt.title('Protein Diversity Coverage by Individual Reference Proteomes')
plt.legend()
plt.show()

# Second graph: Line graph for cumulative diversity (percentage)
plt.figure(figsize=(10, 5))
plt.plot(rps_considered, cumulative_percent, marker='o', linestyle='-', color='b')
plt.title('Cumulative Protein Diversity Coverage by Reference Proteomes')
plt.xlabel('Number of Reference Proteomes Considered')
plt.ylabel('Cumulative Percentage of Protein Diversity Covered (%)')
plt.axhline(y=100, color='red', linestyle='--', label='100% Coverage')
plt.xticks(rps_considered)
plt.ylim(0, 100)  # Ensure y-axis goes up to 100%
plt.legend()
plt.show()

#########

import matplotlib.pyplot as plt
import numpy as np

# Sample data: Randomly generated frequencies of protein occurrence
np.random.seed(0)
total_proteins = 1000
frequencies = np.random.randint(1, 100, total_proteins)  # Frequency percentage of each protein

# Assume 30% of these are documented in UniProt
in_uniprot = np.random.choice([True, False], total_proteins, p=[0.3, 0.7])

# Frequency bins
bins = np.linspace(0, 100, 11)  # 0%, 10%, 20%, ..., 100%

# Histogram data
hist_all, _ = np.histogram(frequencies, bins)
hist_uniprot, _ = np.histogram(frequencies[in_uniprot], bins)

# Plotting
plt.figure(figsize=(10, 6))
plt.bar(bins[:-1], hist_all, width=np.diff(bins), align='edge', alpha=0.6, color='red', label='ggCaller Only')
plt.bar(bins[:-1], hist_uniprot, width=np.diff(bins), align='edge', alpha=0.6, color='green', label='Represented in UniProt')

plt.xlabel('Frequency of Occurrence (%)')
plt.ylabel('Number of Proteins')
plt.title('Protein Frequency Distribution in ggCaller vs. UniProt Coverage')
plt.xticks(bins)
plt.legend()
plt.grid(True)
plt.show()


#########

import matplotlib.pyplot as plt
import numpy as np

# Define core and accessory based on frequency
core_threshold = 80  # percent
accessory_threshold = 20  # percent

# Create mock data
total_proteins = 1000
frequencies = np.random.randint(1, 100, total_proteins)
in_uniprot = np.random.choice([True, False], total_proteins, p=[0.3, 0.7])

# Define core and accessory
is_core = frequencies >= core_threshold
is_accessory = frequencies < accessory_threshold

# Plotting
plt.figure(figsize=(12, 8))
plt.hist(frequencies[is_core], bins=10, color='blue', alpha=0.6, label='Core Proteins')
plt.hist(frequencies[is_accessory], bins=10, color='red', alpha=0.6, label='Accessory Proteins')
plt.xlabel('Frequency of Occurrence (%)')
plt.ylabel('Number of Proteins')
plt.title('Core vs. Accessory Proteins in S. pneumoniae Pan-Proteome')
plt.legend()
plt.grid(True)
plt.show()
#########

import matplotlib.pyplot as plt
import numpy as np

# Sample data generation for demonstration
np.random.seed(0)
total_proteins = 1000
frequencies = np.random.randint(1, 100, total_proteins)  # Random frequency of each protein
in_uniprot = np.random.choice([True, False], total_proteins, p=[0.3, 0.7])  # Whether each protein is in UniProt

# Define thresholds for core and accessory proteins
core_threshold = 80  # Percent, threshold for core proteins
accessory_threshold = 20  # Percent, threshold for accessory proteins

# Determine core and accessory proteins
is_core = frequencies >= core_threshold
is_accessory = frequencies < accessory_threshold

# Data for plotting
core_in_uniprot = frequencies[is_core & in_uniprot]
core_not_in_uniprot = frequencies[is_core & ~in_uniprot]
accessory_in_uniprot = frequencies[is_accessory & in_uniprot]
accessory_not_in_uniprot = frequencies[is_accessory & ~in_uniprot]

# Bins for the histogram
bins = np.linspace(0, 100, 11)

# Plotting
plt.figure(figsize=(12, 8))

# Core proteins
plt.hist([core_in_uniprot, core_not_in_uniprot], bins=bins, stacked=True, label=['Core in UniProt', 'Core not in UniProt'], color=['blue', 'lightblue'], alpha=0.8)

# Accessory proteins
plt.hist([accessory_in_uniprot, accessory_not_in_uniprot], bins=bins, stacked=True, label=['Accessory in UniProt', 'Accessory not in UniProt'], color=['red', 'pink'], alpha=0.8)

plt.xlabel('Frequency of Occurrence (%)')
plt.ylabel('Number of Proteins')
plt.title('Core vs. Accessory Proteins in S. pneumoniae Pan-Proteome (UniProt Coverage)')
plt.legend()
plt.grid(True)
plt.show()



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates

# Data Preparation
data = {
    'RP': ['RP1', 'RP2', 'RP3', 'RP4', 'RP5'],
    'Common Proteins': [200, 180, 220, 210, 205],
    'Unique Proteins': [50, 70, 65, 60, 80],
    'Fragments': [20, 15, 25, 20, 30],
    'Shannon Index': [1.8, 1.9, 1.85, 1.88, 1.9],
    'Reciprocal Simpson': [0.9, 0.92, 0.91, 0.93, 0.95]
}

df = pd.DataFrame(data)

# Plotting
plt.figure(figsize=(10, 5))
parallel_coordinates(df, 'RP', color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
plt.title('Parallel Coordinates Plot of Proteome Diversity across Reference Proteomes')
plt.ylabel('Value')
plt.grid(True)
plt.legend(title='Reference Proteome', loc='upper right')
plt.show()


import matplotlib.pyplot as plt
import numpy as np

# Example protein frequency data
proteins = {
    "Protein A": 20,
    "Protein B": 15,
    "Protein C": 30,
    "Protein D": 10,
    "Protein E": 25
}

# Calculate frequencies as percentages of total samples
total_samples = 100  # Assume 100 total samples in the population
frequencies = np.array(list(proteins.values())) / total_samples * 100

# Sort frequencies for cumulative plot
sorted_frequencies = np.sort(frequencies)
cumulative = np.cumsum(sorted_frequencies)

# Plot cumulative distribution
plt.figure(figsize=(8, 6))
plt.plot(sorted_frequencies, cumulative, marker='o', label='Cumulative Frequency')
plt.axhline(y=95, color='r', linestyle='--', label='95% Threshold (Core Proteins)')
plt.axhline(y=5, color='g', linestyle='--', label='5% Threshold (Rare Proteins)')
plt.xlabel('Protein Frequency (%)')
plt.ylabel('Cumulative Proportion (%)')
plt.title('Cumulative Distribution of Protein Frequencies')
plt.legend()
plt.grid(True)
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Example protein presence-absence matrix
proteins = ['Protein A', 'Protein B', 'Protein C', 'Protein D', 'Protein E']
samples = ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5']

# Binary matrix (1: present, 0: absent)
presence_absence = np.array([
    [1, 1, 1, 0, 1],
    [1, 0, 1, 1, 0],
    [0, 1, 1, 1, 1],
    [1, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Create the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(presence_absence, annot=True, cmap='coolwarm', xticklabels=samples, yticklabels=proteins, cbar=False)
plt.title('Protein Presence-Absence Heatmap')
plt.xlabel('Samples')
plt.ylabel('Proteins')
plt.show()

#pip install upsetplot

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet

# Example Data: Protein presence across datasets (UniProt, ggCaller, and Panproteome)
data = {
    ('UniProt',): 10,  # Proteins unique to UniProt
    ('ggCaller',): 15,  # Proteins unique to ggCaller
    ('Panproteome',): 20,  # Proteins unique to Panproteome
    ('UniProt', 'ggCaller'): 25,  # Proteins shared between UniProt and ggCaller
    ('UniProt', 'Panproteome'): 10,  # Proteins shared between UniProt and Panproteome
    ('ggCaller', 'Panproteome'): 5,  # Proteins shared between ggCaller and Panproteome
    ('UniProt', 'ggCaller', 'Panproteome'): 30,  # Proteins shared across all datasets
}

# Convert data into a pandas Series for UpSet plot
data_series = pd.Series(data)

# Create and plot the UpSet plot
upset = UpSet(data_series, subset_size='count', show_counts=True, sort_by='degree')
plt.figure(figsize=(12, 8))
upset.plot()
plt.suptitle('UpSet Plot: Protein Overlaps Across Datasets', y=1.02, fontsize=14)
plt.show()


import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet

# Example Data: Protein presence across datasets (UniProt, ggCaller, and Panproteome)
data = {
    (True, False, False): 10,  # Proteins unique to UniProt
    (False, True, False): 15,  # Proteins unique to ggCaller
    (False, False, True): 20,  # Proteins unique to Panproteome
    (True, True, False): 25,  # Proteins shared between UniProt and ggCaller
    (True, False, True): 10,  # Proteins shared between UniProt and Panproteome
    (False, True, True): 5,   # Proteins shared between ggCaller and Panproteome
    (True, True, True): 30,   # Proteins shared across all datasets
}

# Convert data into a pandas Series for UpSet plot
data_series = pd.Series(data, index=pd.MultiIndex.from_tuples(
    data.keys(), names=['UniProt', 'ggCaller', 'Panproteome']))

# Create and plot the UpSet plot
upset = UpSet(data_series, subset_size='count', show_counts=True, sort_by='degree')
plt.figure(figsize=(12, 8))
upset.plot()
plt.suptitle('UpSet Plot: Protein Overlaps Across Datasets', y=1.02, fontsize=14)
plt.show()


######
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import entropy

# Sample data
data = {
    'Protein_A': 20,
    'Protein_B': 15,
    'Protein_C': 30,
    'Protein_D': 10,
    'Protein_E': 25
}

# Convert to numpy array for calculations
frequencies = np.array(list(data.values()))
total_proteomes = np.sum(frequencies)

# Calculate proportions
proportions = frequencies / total_proteomes

# Shannon diversity index
shannon_diversity = entropy(proportions, base=2)

# Simpson diversity index
simpson_diversity = np.sum(proportions**2)

# Inverse Simpson diversity index
inverse_simpson = 1 / simpson_diversity

# Pielou's evenness
max_shannon = np.log2(len(proportions))
pielou_evenness = shannon_diversity / max_shannon

# Coverage
coverage = {protein: freq / total_proteomes for protein, freq in data.items()}

# Print results
print(f"Shannon Diversity Index: {shannon_diversity:.4f}")
print(f"Simpson Diversity Index: {simpson_diversity:.4f}")
print(f"Inverse Simpson Diversity Index: {inverse_simpson:.4f}")
print(f"Pielou's Evenness: {pielou_evenness:.4f}")
print("\nProtein Coverage:")
for protein, cov in coverage.items():
    print(f"{protein}: {cov:.2%}")

# Plotting Diversity Indices
diversity_indices = [shannon_diversity, simpson_diversity, inverse_simpson, pielou_evenness]
index_names = ['Shannon', 'Simpson', 'Inverse Simpson', "Pielou's Evenness"]

plt.figure(figsize=(10, 5))
plt.bar(index_names, diversity_indices, color=['blue', 'orange', 'green', 'red'])
plt.title('Diversity Indices')
plt.ylabel('Index Value')
plt.ylim(0, max(diversity_indices) + 0.1)  # Adjust y-axis limit for better visualization
plt.grid(axis='y')
plt.show()

# Plotting Protein Coverage
proteins = list(coverage.keys())
coverage_values = list(coverage.values())

plt.figure(figsize=(10, 5))
plt.bar(proteins, coverage_values, color='purple')
plt.title('Protein Coverage')
plt.ylabel('Coverage (%)')
plt.ylim(0, 1)  # Coverage is a proportion (0 to 1)
plt.xticks(rotation=45)
plt.grid(axis='y')
plt.axhline(y=0.5, color='r', linestyle='--')  # Optional: Add a line at 50% coverage for reference
plt.show()


#####

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import entropy

# Sample data
data = {
    'Protein_A': 20,
    'Protein_B': 15,
    'Protein_C': 30,
    'Protein_D': 10,
    'Protein_E': 25
}

proteins = list(data.keys())
frequencies = np.array(list(data.values()))
total_proteomes = np.sum(frequencies)
proportions = frequencies / total_proteomes

# Calculate Shannon diversity contribution for each protein
shannon_contributions = -proportions * np.log2(proportions)
shannon_diversity = np.sum(shannon_contributions)

# Calculate coverage
coverage = frequencies / total_proteomes

# Create the plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Diversity contribution plot
ax1.bar(proteins, shannon_contributions)
ax1.set_title('Shannon Diversity Contribution by Protein')
ax1.set_ylabel('Diversity Contribution')
ax1.text(0.02, 0.95, f'Total Shannon Diversity: {shannon_diversity:.4f}', 
         transform=ax1.transAxes, verticalalignment='top')

# Coverage plot
ax2.bar(proteins, coverage)
ax2.set_title('Protein Coverage')
ax2.set_ylabel('Coverage')
ax2.set_ylim(0, 1)

# Add percentage labels on top of each bar
for i, v in enumerate(coverage):
    ax2.text(i, v, f'{v:.1%}', ha='center', va='bottom')

plt.tight_layout()
plt.show()

# Print additional diversity indices
simpson_diversity = np.sum(proportions**2)
inverse_simpson = 1 / simpson_diversity
max_shannon = np.log2(len(proportions))
pielou_evenness = shannon_diversity / max_shannon

print(f"Shannon Diversity Index: {shannon_diversity:.4f}")
print(f"Simpson Diversity Index: {simpson_diversity:.4f}")
print(f"Inverse Simpson Diversity Index: {inverse_simpson:.4f}")
print(f"Pielou's Evenness: {pielou_evenness:.4f}")

#####


import matplotlib.pyplot as plt
import numpy as np

# Sample data
proteins = ['A', 'B', 'C', 'D', 'E']
rps = ['RP1', 'RP2', 'RP3', 'RP4', 'RP5']

# Generate random data for demonstration
np.random.seed(42)
frequencies = np.random.randint(10, 100, size=(len(rps), len(proteins)))
coverages = frequencies / frequencies.sum(axis=1)[:, np.newaxis]
diversity_indices = -np.sum(coverages * np.log2(coverages), axis=1)

# Normalize diversity contribution for visualization
diversity_contrib = coverages * np.log2(coverages) / diversity_indices[:, np.newaxis]

# Set up the plot
fig, ax1 = plt.subplots(figsize=(15, 10))
ax2 = ax1.twinx()

# Width of each bar group
width = 0.15

# Colors for each RP
colors = plt.cm.get_cmap('Set1')(np.linspace(0, 1, len(rps)))

# Plot bars for frequencies
for i, rp in enumerate(rps):
    x = np.arange(len(proteins)) + i * width
    ax1.bar(x, frequencies[i], width, label=f'{rp} Frequency', color=colors[i], alpha=0.7)

# Plot lines for coverage and diversity
for i, rp in enumerate(rps):
    x = np.arange(len(proteins)) + i * width
    ax2.plot(x, coverages[i], 'o-', label=f'{rp} Coverage', color=colors[i])
    ax2.plot(x, -diversity_contrib[i], 's--', label=f'{rp} Diversity Contrib', color=colors[i])

# Customize the plot
ax1.set_xlabel('Proteins')
ax1.set_ylabel('Frequency')
ax2.set_ylabel('Coverage / Diversity Contribution')
ax1.set_title('Protein Analysis across Reference Proteomes')
ax1.set_xticks(np.arange(len(proteins)) + width * (len(rps) - 1) / 2)
ax1.set_xticklabels(proteins)

# Add legends
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Adjust layout and display
plt.tight_layout()
plt.show()

# Print diversity indices
for i, rp in enumerate(rps):
    print(f"{rp} Shannon Diversity Index: {diversity_indices[i]:.4f}")

#########

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import plotly.io as pio
import plotly.express as px
pio.renderers.default='browser'

# Sample data for 3 RPs and 3 proteins
rps = ['RP1', 'RP2', 'RP3']
proteins = ['A', 'B', 'C']

# Generate sample data
np.random.seed(42)
frequencies = np.random.randint(10, 100, size=(len(rps), len(proteins)))
coverages = frequencies / frequencies.sum(axis=1)[:, np.newaxis]
diversity_indices = -np.sum(coverages * np.log2(coverages), axis=1)

# Create subplots
fig = make_subplots(rows=2, cols=1, shared_xaxes=True, 
                    vertical_spacing=0.1,
                    subplot_titles=("Protein Frequencies", "Coverage and Diversity"))

# Colors for each protein
colors = ['red', 'green', 'blue']

# Add frequency bars
for i, protein in enumerate(proteins):
    fig.add_trace(go.Bar(x=rps, y=frequencies[:, i], name=f'Protein {protein} Frequency',
                         marker_color=colors[i]), row=1, col=1)

# Add coverage lines
for i, protein in enumerate(proteins):
    fig.add_trace(go.Scatter(x=rps, y=coverages[:, i], mode='lines+markers',
                             name=f'Protein {protein} Coverage', line=dict(color=colors[i])), 
                  row=2, col=1)

# Add diversity index
fig.add_trace(go.Scatter(x=rps, y=diversity_indices, mode='lines+markers',
                         name='Diversity Index', line=dict(color='purple', dash='dash')),
              row=2, col=1)

# Update layout
fig.update_layout(height=700, width=800, title_text="Protein Analysis Across RPs",
                  hovermode="x unified")
fig.update_xaxes(title_text="Reference Proteomes", row=2, col=1)
fig.update_yaxes(title_text="Frequency", row=1, col=1)
fig.update_yaxes(title_text="Coverage / Diversity", row=2, col=1)

# Show plot
fig.show()

####
import matplotlib.pyplot as plt
import numpy as np

# Sample data for proteins and their frequencies
proteins = ['A', 'B', 'C', 'D', 'E']
frequencies = [20, 15, 30, 10, 25]

# Create a bar chart
plt.figure(figsize=(10, 6))
x = np.arange(len(proteins))  # the label locations
plt.bar(x, frequencies, color='skyblue')

# Add labels and title
plt.xticks(x, proteins)
plt.xlabel('Proteins')
plt.ylabel('Frequency')
plt.title('Frequency of Protein in all proteomes in UP')

# Add frequency labels on top of each bar
for i, v in enumerate(frequencies):
    plt.text(i, v + 0.5, str(v), ha='center')

plt.tight_layout()
plt.show()

##########
import matplotlib.pyplot as plt
import numpy as np

# Sample data for diversity indices
proteome_sets = ['RP1', 'RP2', 'RP3', 'Benchmark (S. pneumoniae)']
diversity_indices = [3.5, 4.0, 2.8, 4.2]

# Create a bar chart
plt.figure(figsize=(10, 6))
x = np.arange(len(proteome_sets))  # the label locations
plt.bar(x, diversity_indices, color=['skyblue', 'lightgreen', 'salmon', 'orange'])

# Add labels and title
plt.xticks(x, proteome_sets)
plt.xlabel('Proteome Sets')
plt.ylabel('Shannon Diversity Index (OR ANY other suitable metric)')
plt.title('Comparison of Diversity Indices Across Proteome Sets')

# Add value labels on top of each bar
for i, v in enumerate(diversity_indices):
    plt.text(i, v + 0.05, f'{v:.1f}', ha='center')

plt.ylim(0, max(diversity_indices) + 1)  # Adjust y-axis limit for better visualization
plt.axhline(y=4.0, color='r', linestyle='--', label='Average RP Diversity')  # Optional: Average line for RPs

# Add legend
plt.legend()
plt.tight_layout()
plt.show()
