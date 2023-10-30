#!/usr/bin/env python3
"""
Author: Mohamed S. Sarhan
Email: mohamed.sarhan@eurac.edu; m.sabrysarhan@gmail.com
Date: September 29, 2023
Description: This Python script filters VEP annotation output file to keep the hightest imapct for each variant.
"""
import argparse
import pandas as pd
import os

def read_vep_annotation_file(file_path):
    """
    Read a VEP output annotation file into a Pandas DataFrame, excluding lines starting with '##'.
    The 'Extra' column is parsed into dictionaries for each row.

    Args:
        file_path (str): The path to the VEP output annotation file.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the data from the VEP output file with separate columns for IMPACT, DISTANCE, and STRAND, and 'Extra' parsed as dictionaries.
    """
    try:
        lines = []
        with open(file_path, 'r') as file, open(file_path + '.tmp', 'w') as output:
            for line in file:
                if not line.startswith('#'):
                    output.write(line)

        df = pd.read_csv(file_path + '.tmp', sep='\t', header=None, engine='python')
        # Add column names based on the VEP output format, you may need to adjust this
        df.columns = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type',
                      'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids',
                      'Codons', 'Existing_variation', 'Extra']

        # Parse the 'Extra' column into dictionaries
        df['Extra'] = df['Extra'].apply(lambda x: dict(item.split('=') for item in x.split(';') if '=' in item))

        os.remove(file_path + '.tmp')
        return df
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None


def select_highest_impact_rows(df):
    """
    Group rows by the 'Location' column and select rows based on the lowest distance value in the 'Extra' column,
    or append the row if 'DISTANCE' is not found.

    Args:
        df (pd.DataFrame): A Pandas DataFrame containing VEP annotation data with 'Extra' column as a dictionary.

    Returns:
        pd.DataFrame: A Pandas DataFrame with rows grouped by 'Location' and selected based on the specified conditions.
    """
    if 'Location' not in df.columns:
        raise ValueError("DataFrame must contain a 'Location' column for grouping.")

    grouped = df.groupby('Location')
    selected_rows = []

    for _, group in grouped:
        if len(group) == 1:
            selected_rows.append(group.iloc[0])
        else:
            lowest_distance_idx = None
            for index, row in group.iterrows():
                if 'DISTANCE' not in row['Extra']:
                    selected_rows.append(row)
                    break

            if lowest_distance_idx is None or row['Extra']['DISTANCE'] < group.loc[lowest_distance_idx]['Extra']['DISTANCE']:
                lowest_distance_idx = index
                selected_rows.append(group.loc[lowest_distance_idx])

    selected_df = pd.DataFrame(selected_rows)
    #selected_df['Extra'] = selected_df['Extra'].apply(frozenset)
    return selected_df.drop_duplicates(subset = ['Uploaded_variation'])


import pandas as pd

def sort_df_by_index(df):
    """
    Sort a DataFrame by its index.

    Args:
        df (pd.DataFrame): A Pandas DataFrame to be sorted.

    Returns:
        pd.DataFrame: A sorted Pandas DataFrame based on the index.
    """
    sorted_df = df.sort_index()

    return sorted_df

def rename_snps_in_df(df, prefix):
    """
    Rename SNPs in the 'uploaded_variation' column of a DataFrame using three words.

    Args:
        df (pd.DataFrame): The DataFrame containing the SNP data from VEP.
        rename_words (list): A list of three words to use for renaming.

    Returns:
        pd.DataFrame: A DataFrame with the 'uploaded_variation' column (SNPs) renamed based on the provided words.
    """

    for index, row in df.iterrows():
        var_name = row['Uploaded_variation'].split('_')
        var_new_name = prefix[0] +':' + 'contig_' + var_name[0] + ':' + var_name[1] 
        row['Uploaded_variation'] = var_new_name


    return df


def main():
    parser = argparse.ArgumentParser(description='Filter Variant Effect Predictor (VEP) annotation to keep only the highest impact effect.')
    parser.add_argument('--vep', required=True, help='Input VEP output file path', default=None)
    parser.add_argument('--rename', required=False, help='prefix', default=None)
    parser.add_argument('--out', required=True, help='Output file path', default=None)

    args = parser.parse_args()

    # Example usage:
    file_path = args.vep

    vep_df = read_vep_annotation_file(file_path)

    selected_df = select_highest_impact_rows(vep_df)


    sorted_vep_df = sort_df_by_index(selected_df)

    #print(sorted_vep_df)
    if args.rename == None:
        sorted_vep_df.to_csv(args.out, sep='\t', index=False)
    else:
        words = list(args.rename.split(","))
        #print(words)
        result_df = rename_snps_in_df(sorted_vep_df, words)
        #print(result_df)
        result_df.to_csv(args.out, sep='\t', index=False)


if __name__ == "__main__":
    main()