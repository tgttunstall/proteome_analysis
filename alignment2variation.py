#!/usr/bin/env python
# coding: utf-8

# # alignment 2 variation

# In[ ]:


# changelog
# Wed  4 Dec 2024 11:00:37 GMT 0.1 started work
# Wed  4 Dec 2024 11:52:27 GMT 0.1 vectorised implementation (align2vardf_vec)
# Wed  4 Dec 2024 12:53:18 GMT 0.2 performance tests of vectorised implementation vs for loop
# Wed  4 Dec 2024 14:00:18 GMT 0.3 find_rows_with_rare_variants


# In[ ]:


from Bio import Phylo
from Bio import AlignIO
import pandas as pd
import numpy as np
from tqdm import tqdm

import timeit  # performance evaluation


# In[ ]:


# read alignment
ALIGN_FORMAT = "fasta"
path = "/Users/insana/workEBI/proteome_alignment_idea/mmseq_labelling/"

aln_fname = "986.aln"  # small pig alignment
alignment_small = AlignIO.read(path + aln_fname, ALIGN_FORMAT)

aln_fname = "829.aln"  # large ecoli alignment
alignment_large = AlignIO.read(path + aln_fname, ALIGN_FORMAT)


# In[ ]:


def align2vardf_loop(alignment, nogaps=True):
    variation_dict = {}

    aln_length = alignment.get_alignment_length()

    if nogaps:
        for pos in range(0, aln_length):
            column = alignment[:, pos]

            prev_res = column[0]

            if prev_res == "-":
                continue

            found_var = False
            for res in column[1:]:
                if res == "-":
                    found_var = False
                    break
                if res != prev_res:
                    found_var = True
            if found_var:
                # print(f"{pos} {column}")
                variation_dict[pos] = list(column)
    else:
        for pos in range(0, aln_length):
            column = alignment[:, pos]

            prev_res = column[0]
            for res in column[1:]:
                if res != prev_res:
                    # print(f"{pos} {column}") #debug
                    variation_dict[pos] = list(column)
                    break

    # create variation df
    variation_df = pd.DataFrame.from_dict(variation_dict, orient="columns")
    variation_df.index = [alignment[i].id for i in range(len(alignment))]
    return variation_df


variation_df = align2vardf_loop(alignment_small)
variation_df


# In[ ]:


def align2df(alignment):
    return pd.DataFrame(
        np.array([list(record.seq) for record in alignment]),
        index=[record.id for record in alignment],
        columns=range(0, alignment.get_alignment_length()),
    )


# example:
alignment_df = align2df(alignment_small)


def align2vardf_vec(alignment, nogaps=True):
    """
    Create a DataFrame containing only the columns from an alignment with variation
    (i.e. where sequence identity is not 100% over all the sequences of the alignment).

    Parameters:
    alignment -- Bio.Align object
    nogaps -- bool, if True, excludes columns with gaps

    Returns:
    A pandas DataFrame where rows are sequence IDs and columns are positions with variation.
    """

    # convert to numpy
    alignment_array = np.array([list(record.seq) for record in alignment])

    # get reference row for comparison
    ref_row = alignment_array[0]

    # find columns with variation
    is_variable = np.any(alignment_array != ref_row, axis=0)

    if nogaps:
        is_gap_free = ~np.any(alignment_array == "-", axis=0)  # gap free columns
        valid_columns = is_gap_free & is_variable
    else:
        valid_columns = is_variable

    # valid_positions = np.where(valid_columns)[0]
    valid_positions = np.flatnonzero(valid_columns)  # faster?

    # create variation df
    variation_df = pd.DataFrame(
        alignment_array[:, valid_columns],
        index=[record.id for record in alignment],
        columns=valid_positions,
    )

    return variation_df


variation_df = align2vardf_vec(alignment_large, nogaps=True)
variation_df


# In[ ]:


def compute_variant_stats(variation_df):
    """
    Compute statistics for each column in an alignment.

    Parameters:
    variation_df -- pandas DataFrame with sequence IDs as rows and variable positions as columns.

    Returns:
    stats_df -- pandas DataFrame with variable positions as columns, residues as rows,
                and relative abundances as values.
    """
    # Compute frequencies for each position
    stats_dict = {}

    for col in variation_df.columns:
        counts = variation_df[col].value_counts(normalize=True)  # relative frequency
        stats_dict[col] = counts

    # Combine all counts into a DataFrame
    stats_df = pd.DataFrame(stats_dict).fillna(
        0
    )  # Fill missing values with 0 for residues not present

    return stats_df


stats_df = compute_variant_stats(variation_df)
round(stats_df, 3)


# In[ ]:


def find_rows_with_rare_variants(variation_df, stats_df, threshold=2):
    """
    Identify rows where two or more columns have rare variants.

    Parameters:
    variation_df -- pandas DataFrame with sequence IDs as rows and variable positions as columns.
    stats_df -- pandas DataFrame with variable positions as columns, residues as rows,
                and relative abundances as values.
    threshold -- int, minimum number of rare variants in a row to flag it (default: 2).

    Returns:
    rare_rows -- pandas Index of sequence IDs with two or more rare variants.
    rare_variants_dict -- dict where keys are sequence IDs and values are lists of positions
                          with rare variants.
    """
    # consensus for each position
    consensus = stats_df.idxmax(axis=0)  # most frequent residue for each position

    # create a boolean DataFrame where True indicates a rare variant
    rare_variants = variation_df != consensus

    # count the number of rare variants in each row
    rare_counts = rare_variants.sum(axis=1)

    # identify rows with >= threshold rare variants
    covariating_rows = rare_counts[rare_counts >= threshold].index

    # collect positions of rare variants for flagged rows
    rare_variants_dict = {
        row: rare_variants.columns[rare_variants.loc[row]].tolist()
        for row in covariating_rows
    }

    return covariating_rows, rare_variants_dict


variants_index, variants_dict = find_rows_with_rare_variants(variation_df, stats_df)
# variants_dict


# In[ ]:


# show first and last variant row:
display(variation_df.loc[variants_index[0:1]][variants_dict[variants_index[0]]])
display(variation_df.loc[variants_index[-1:]][variants_dict[variants_index[-1]]])

# show all:
# for variant_protein in variants_index:
#   print(variation_df.loc[variant_protein][variants_dict[variant_protein]])


# In[ ]:


# Evaluating performance for two alignments
# test_repeat_count = 10
test_repeat_count = 1
print(f"Repeating each performance test {test_repeat_count} times")
for alignment in [alignment_small, alignment_large]:
    num_seqs = len(alignment)
    aln_length = alignment.get_alignment_length()
    print(f"Testing on alignment with {num_seqs} sequences and {aln_length} length")
    for nogaps in [True, False]:
        if nogaps:
            print("\texcluding positions with gaps")
        else:
            print("\tincluding positions with gaps")
        for function_to_test in [align2vardf_loop, align2vardf_vec]:
            execution_time = timeit.timeit(
                lambda: function_to_test(alignment, nogaps=nogaps),
                number=test_repeat_count,
            )
            print(f"\t\t{function_to_test.__name__}: {execution_time:.4f} seconds")
"""
Repeating each performance test 10 times
Testing on alignment with 66 sequences and 480 length
	excluding positions with gaps
		align2vardf_np: 0.0486 seconds
		align2vardf_for: 0.2171 seconds
	including positions with gaps
		align2vardf_np: 0.0251 seconds
		align2vardf_for: 0.2071 seconds
Testing on alignment with 105522 sequences and 424 length
	excluding positions with gaps
		align2vardf_np: 38.3389 seconds
		align2vardf_for: 311.9611 seconds
	including positions with gaps
		align2vardf_np: 39.4331 seconds
		align2vardf_for: 299.2920 seconds
"""
