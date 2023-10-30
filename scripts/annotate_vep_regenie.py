#!/bin/usr/env python3
"""
Author: Mohamed S. Sarhan
Email: mohamed.sarhan@eurac.edu; m.sabrysarhan@gmail.com
Date: October 3, 2023
Description: This Python script combine VEP annotation output file and regenie output by SNP name.
"""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='VEP annotate regenie output.')
    parser.add_argument('--vep', required=True, help='Input VEP highest impact SNP', default=None)
    parser.add_argument('--regenie', required=True, help='regenie output file', default=None)
    parser.add_argument('--out', required=True, help='Output file path for annotated regenie', default=None)

    args = parser.parse_args()

    # Read the VEP annotation file into a Pandas DataFrame
    vep_df = pd.read_csv(args.vep, sep='\t', engine='python')

    # Read the regenie output file into a Pandas DataFrame
    regenie_df = pd.read_csv(args.regenie, sep='\t', engine='python')

    # Perform an inner join on the specified columns
    merged_df = vep_df.merge(regenie_df, left_on='Uploaded_variation', right_on='ID', how='inner')

    # Check if the merged DataFrame is empty
    if merged_df.empty:
        # Return an empty DataFrame
        result_df = pd.DataFrame(columns=vep_df.columns.tolist() + regenie_df.columns.tolist())
        result_df.to_csv(args.out, sep='\t', index = False)
    else:
        # Return the merged DataFrame
        result_df = merged_df
        result_df.to_csv(args.out, sep='\t', index = False)


if __name__ == "__main__":
    main()


