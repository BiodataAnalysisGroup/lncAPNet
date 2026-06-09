#!/usr/bin/env python
import argparse
import gseapy as gp
import pandas as pd
import numpy as np
import os
import openpyxl
from sklearn.model_selection import train_test_split
import gc

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Pathway Enrichment Analysis (GO only)')
    parser.add_argument('--ms_tab', required=True, help='Path to ms_tab Excel file')
    parser.add_argument('--gmt_GO', required=True, help='Path to GO GMT file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--activity', required=False, help='Path to activity.csv file for PASNet input generation')
    parser.add_argument('--metadata', required=False, help='Path to metadata.csv file for PASNet input generation')
    parser.add_argument('--p_value', type=float, default=0.05, help='P-value threshold (default: 0.05)')
    parser.add_argument('--size_threshold', type=int, default=30, help='Minimum size threshold (default: 30)')
    parser.add_argument('--cutoff', type=float, default=0.1, help='Enrichment cutoff (default: 0.1)')
    parser.add_argument('--group0', required=True, help='Control group name')
    parser.add_argument('--group1', required=True, help='Treatment group name')
    parser.add_argument('--comparison', required=True, help='Comparison name')
    return parser.parse_args()

def filter_ms_tab(ms_tab_file, group0, group1, p_value=0.05, size_threshold=30):
    """Filter ms_tab based on conditions and return sorted unique gene list"""
    ms_tab = pd.read_excel(ms_tab_file)
    
    # Dynamically construct column names based on group parameters
    adj_pval_col = f"adj.P.Val.{group1}.Vs.{group0}_DA"
    logfc_col = f"logFC.{group1}.Vs.{group0}_DA"
    
    # Filtering based on conditions
    pos = ms_tab[(ms_tab[adj_pval_col] < p_value) & 
                 (ms_tab[logfc_col] > 0) & 
                 (ms_tab["Size"] > size_threshold)]
    neg = ms_tab[(ms_tab[adj_pval_col] < p_value) & 
                 (ms_tab[logfc_col] < 0) & 
                 (ms_tab["Size"] > size_threshold)]
    
    # Row-bind (combine) the results
    ms_tab_filtered = pd.concat([pos, neg], axis=0, ignore_index=True)
    
    # Get unique values from the 'geneSymbol' column
    gene_list = list(ms_tab_filtered['geneSymbol'].unique())
    gene_list.sort()
    
    print(f"Filtered gene list contains {len(gene_list)} genes")
    return gene_list

def run_enrichment(gene_list, gmt_file, outdir, name, cutoff=0.1, p_value=0.05):
    """Run enrichment analysis with error handling for empty results"""
    print(f"Running enrichment for {name}...")
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list, 
            gene_sets=gmt_file,
            organism="human", 
            outdir=os.path.join(outdir, f"{name}_res"),
            cutoff=cutoff,
            no_plot=True
        )
        
        results = enr.results
        
        # Check if results are empty
        if results.empty:
            print(f"{name}: No enrichment terms found with cutoff={cutoff}")
            return pd.DataFrame()
        
        # Filter by p-value
        results_filt = results[results['P-value'] < p_value]
        print(f"{name}: {len(results_filt)} pathways with P-value < {p_value}")
        
        if results_filt.empty:
            print(f"{name}: No pathways passed the filtering criteria (P-value < {p_value})")
            return pd.DataFrame()
        
        return results_filt
        
    except ValueError as e:
        if "No enrich terms" in str(e):
            print(f"{name}: No enrichment terms found with cutoff={cutoff}")
            return pd.DataFrame()
        else:
            raise
    except Exception as e:
        print(f"{name}: Error during enrichment analysis: {str(e)}")
        return pd.DataFrame()

def process_results_to_matrix(results_df):
    """Convert results to gene-pathway matrix"""
    if results_df.empty:
        return pd.DataFrame()
    
    results_df = results_df.copy()
    results_df['Genes'] = results_df['Genes'].str.split(';')
    results_matrix = results_df.explode("Genes").pivot_table(
        index="Term", 
        columns="Genes", 
        aggfunc="size", 
        fill_value=0
    ).reset_index()
    results_matrix = results_matrix.set_index('Term')
    return results_matrix

def generate_pasnet_inputs(activity_file, metadata_file, pt_go, outdir, comparison, group0, group1):
    """Generate PASNet input files"""
    print("Generating PASNet input files...")
    
    # Load activity data
    activity = pd.read_csv(activity_file)
    activity['Unnamed: 0'] = activity['Unnamed: 0'].str.replace('_TF', '').str.replace('_SIG', '').str.replace('_LNC', '')
    activity = activity.drop_duplicates(subset=['Unnamed: 0']).set_index('Unnamed: 0').T
    
    # Load metadata
    metadata = pd.read_csv(metadata_file, index_col=0)
    metadata[comparison] = metadata[comparison].replace({group1: '1', group0: '0'})
    metadata = metadata[[comparison]]
    
    
    # Subset activity data to GO genes
    activity_go = activity.loc[:, activity.columns.isin(pt_go.columns)]
    
    # Merge with metadata
    activity_go = activity_go.merge(metadata, left_index=True, right_index=True)
    
    # Split into train and test sets
    pretrain_go, validation_go = train_test_split(activity_go, test_size=0.2, random_state=123, stratify=activity_go[comparison])
    train_go, validation_grinding_go = train_test_split(pretrain_go, test_size=0.2, random_state=123, stratify=pretrain_go[comparison])
    
    # Create output directory
    os.makedirs(os.path.join(outdir, "PASNet/Input/GO"), exist_ok=True)

    # Save to Excel
    train_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Training.xlsx"), index=True)
    validation_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Validation.xlsx"), index=True)
    validation_grinding_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Validation_grinding.xlsx"), index=True)
    
    # Filter pathway matrix to only include genes present in activity data
    pt_go_filt = pt_go.loc[:, pt_go.columns.isin(activity_go.columns)]
    pt_go_filt.to_excel(os.path.join(outdir, "PASNet/Input/GO/pt_fixed.xlsx"))
    
    print("PASNet input files generated successfully")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Filter ms_tab to get gene list
    gene_list = filter_ms_tab(args.ms_tab, args.group0, args.group1, args.p_value, args.size_threshold)
    
    # Run GO enrichment only
    print("\n=== Running Gene Ontology Enrichment ===")
    results_go = run_enrichment(gene_list, args.gmt_GO, args.outdir, "GO_Biological_Process", 
                                args.cutoff, args.p_value)
    
    # Process GO results only if not empty
    if not results_go.empty:
        print("\n=== Processing GO results ===")
        results_go_matrix = process_results_to_matrix(results_go)
        if not results_go_matrix.empty:
            os.makedirs(os.path.join(args.outdir, "PASNet/Input/GO"), exist_ok=True)
            results_go_matrix.to_excel(os.path.join(args.outdir, "PASNet/Input/GO/pt_fixed.xlsx"))
            print(f"GO matrix saved: {results_go_matrix.shape[0]} pathways x {results_go_matrix.shape[1]} genes")
        else:
            print("GO matrix is empty after processing")
            results_go_matrix = pd.DataFrame()
    else:
        print("\n=== Skipping GO results processing (no enrichment found) ===")
        results_go_matrix = pd.DataFrame()
    
    # Generate PASNet inputs if activity and metadata are provided and we have GO results
    if args.activity and args.metadata:
        if not results_go_matrix.empty:
            generate_pasnet_inputs(args.activity, args.metadata, 
                                   results_go_matrix, args.outdir,
                                   args.comparison, args.group0, args.group1)
        else:
            print("\nSkipping PASNet input generation (missing GO results)")
    else:
        print("\nSkipping PASNet input generation (activity and metadata not provided)")
    
    # Garbage collection
    gc.collect()
    
    print("\n=== Enrichment analysis completed successfully ===")

if __name__ == "__main__":
    main()