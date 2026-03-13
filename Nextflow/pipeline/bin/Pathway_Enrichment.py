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
    parser = argparse.ArgumentParser(description='Pathway Enrichment Analysis')
    parser.add_argument('--ms_tab', required=True, help='Path to ms_tab Excel file')
    parser.add_argument('--gmt_GO', required=True, help='Path to GO GMT file')
    parser.add_argument('--gmt_Reactome', required=True, help='Path to Reactome GMT file')
    parser.add_argument('--gmt_KEGG', required=True, help='Path to KEGG GMT file')
    parser.add_argument('--gmt_Wikipathway', required=True, help='Path to Wikipathway GMT file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--activity', required=False, help='Path to activity.csv file for PASNet input generation')
    parser.add_argument('--metadata', required=False, help='Path to metadata.csv file for PASNet input generation')
    parser.add_argument('--p_value', type=float, default=0.05, help='P-value threshold (default: 0.05)')
    parser.add_argument('--adj_p_value', type=float, default=0.05, help='Adjusted P-value threshold (default: 0.05)')
    parser.add_argument('--size_threshold', type=int, default=30, help='Minimum size threshold (default: 30)')
    parser.add_argument('--cutoff', type=float, default=0.1, help='Enrichment cutoff (default: 0.1)')
    parser.add_argument('--group0', required=True, help='Control group name')
    parser.add_argument('--group1', required=True, help='Treatment group name')
    return parser.parse_args()

def filter_ms_tab(ms_tab_file, group0, group1, p_value=0.05, size_threshold=30):
    """Filter ms_tab based on conditions"""
    prostate = pd.read_excel(ms_tab_file)
    
    # Dynamically construct column names based on group parameters
    adj_pval_col = f"adj.P.Val.{group1}.Vs.{group0}_DA"
    logfc_col = f"logFC.{group1}.Vs.{group0}_DA"
    
    # Filtering based on conditions
    pos = prostate[(prostate[adj_pval_col] < p_value) & 
                   (prostate[logfc_col] > 0) & 
                   (prostate["Size"] > size_threshold)]
    neg = prostate[(prostate[adj_pval_col] < p_value) & 
                   (prostate[logfc_col] < 0) & 
                   (prostate["Size"] > size_threshold)]
    
    return pos, neg
    
    # Row-bind (combine) the results
    ms_tab = pd.concat([pos, neg], axis=0, ignore_index=True)
    
    # Get unique values from the 'hgnc_symbol.y' column
    ms_tab = list(ms_tab['geneSymbol'].unique())
    ms_tab.sort()
    
    print(f"Filtered gene list contains {len(ms_tab)} genes")
    return ms_tab

def run_enrichment(gene_list, gmt_file, outdir, name, cutoff=0.1, p_value=0.05, adj_p_value=0.05):
    """Run enrichment analysis with error handling for empty results"""
    print(f"Running enrichment for {name}...")
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list, 
            gene_sets=gmt_file,
            organism="human", 
            outdir=os.path.join(outdir, f"{name}_res"),
            cutoff=cutoff,
            no_plot=True  # Disable plotting to avoid errors when no results
        )
        
        results = enr.results
        
        # Check if results are empty
        if results.empty:
            print(f"{name}: No enrichment terms found with cutoff={cutoff}")
            return pd.DataFrame()  # Return empty DataFrame
        
        # Filter by p-value
        results_filt_pval = results[(results['P-value'] < p_value)]
        print(f"{name}: {len(results_filt_pval)} pathways with P-value < {p_value}")
        
        # Filter by adjusted p-value
        results_filt = results[(results['Adjusted P-value'] < adj_p_value)]
        print(f"{name}: {len(results_filt)} pathways with Adjusted P-value < {adj_p_value}")
        
        if results_filt.empty:
            print(f"{name}: No pathways passed the filtering criteria (Adjusted P-value < {adj_p_value})")
            return pd.DataFrame()
        
        return results_filt
        
    except ValueError as e:
        if "No enrich terms" in str(e):
            print(f"{name}: No enrichment terms found with cutoff={cutoff}")
            return pd.DataFrame()  # Return empty DataFrame
        else:
            raise  # Re-raise if it's a different ValueError
    except Exception as e:
        print(f"{name}: Error during enrichment analysis: {str(e)}")
        return pd.DataFrame()  # Return empty DataFrame on any error

def process_results_to_matrix(results_df):
    """Convert results to gene-pathway matrix"""
    if results_df.empty:
        return pd.DataFrame()
    
    results_df['Genes'] = results_df['Genes'].str.split(';')
    results_matrix = results_df.explode("Genes").pivot_table(
        index="Term", 
        columns="Genes", 
        aggfunc="size", 
        fill_value=0
    ).reset_index()
    results_matrix = results_matrix.set_index('Term')
    return results_matrix

def generate_pasnet_inputs(activity_file, metadata_file, pt_go, pt_rest, outdir):
    """Generate PASNet input files"""
    print("Generating PASNet input files...")
    
    # Load activity data
    activity = pd.read_csv(activity_file)
    activity['Unnamed: 0'] = activity['Unnamed: 0'].str.replace('_TF', '').str.replace('_SIG', '').str.replace('_LNC', '')
    activity = activity.drop_duplicates(subset=['Unnamed: 0']).set_index('Unnamed: 0').T
    
    # Load metadata
    metadata = pd.read_csv(metadata_file, index_col=0)
    metadata['Status'] = metadata['Status'].replace({'Progressive': '1', 'Mild': '0'})
    
    # Subset activity data
    activity_go = activity.loc[:, activity.columns.isin(pt_go.columns)]
    activity_rest = activity.loc[:, activity.columns.isin(pt_rest.columns)]
    
    # Merge with metadata
    activity_go = activity_go.merge(metadata, left_index=True, right_index=True)
    activity_rest = activity_rest.merge(metadata, left_index=True, right_index=True)
    
    # Split into train and test sets
    pretrain_go, validation_go = train_test_split(activity_go, test_size=0.2, random_state=123, stratify=activity_go['Status'])
    pretrain_rest, validation_rest = train_test_split(activity_rest, test_size=0.2, random_state=123, stratify=activity_rest['Status'])
    
    train_go, validation_grinding_go = train_test_split(pretrain_go, test_size=0.2, random_state=123, stratify=pretrain_go['Status'])
    train_rest, validation_grinding_rest = train_test_split(pretrain_rest, test_size=0.2, random_state=123, stratify=pretrain_rest['Status'])
    
    # Create output directories
    os.makedirs(os.path.join(outdir, "PASNet/Input/GO"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "PASNet/Input/REST"), exist_ok=True)
    
    # Save to Excel
    train_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Training.xlsx"), index=True)
    validation_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Validation.xlsx"), index=True)
    validation_grinding_go.to_excel(os.path.join(outdir, "PASNet/Input/GO/Validation_grinding.xlsx"), index=True)
    train_rest.to_excel(os.path.join(outdir, "PASNet/Input/REST/Training.xlsx"), index=True)
    validation_rest.to_excel(os.path.join(outdir, "PASNet/Input/REST/Validation.xlsx"), index=True)
    validation_grinding_rest.to_excel(os.path.join(outdir, "PASNet/Input/REST/Validation_grinding.xlsx"), index=True)
    
    # Filter pathway matrices
    pt_go_filt = pt_go.loc[:, pt_go.columns.isin(activity_go.columns)]
    pt_rest_filt = pt_rest.loc[:, pt_rest.columns.isin(activity_rest.columns)]
    
    pt_rest_filt.to_excel(os.path.join(outdir, "PASNet/Input/REST/pt_fixed.xlsx"))
    pt_go_filt.to_excel(os.path.join(outdir, "PASNet/Input/GO/pt_fixed.xlsx"))
    
    print("PASNet input files generated successfully")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Filter ms_tab to get gene list
    gene_list = filter_ms_tab(args.ms_tab, args.group0, args.group1, args.p_value, args.size_threshold)
    
    # Run enrichment analyses for all four databases
    print("\n=== Running Gene Ontology Enrichment ===")
    results_go = run_enrichment(gene_list, args.gmt_GO, args.outdir, "GO_Biological_Process", 
                                 args.cutoff, args.p_value, args.adj_p_value)
    
    print("\n=== Running KEGG Enrichment ===")
    results_kegg = run_enrichment(gene_list, args.gmt_KEGG, args.outdir, "KEGG", 
                                   args.cutoff, args.p_value, args.adj_p_value)
    
    print("\n=== Running Reactome Enrichment ===")
    results_reactome = run_enrichment(gene_list, args.gmt_Reactome, args.outdir, "Reactome", 
                                       args.cutoff, args.p_value, args.adj_p_value)
    
    print("\n=== Running WikiPathway Enrichment ===")
    results_wiki = run_enrichment(gene_list, args.gmt_Wikipathway, args.outdir, "WikiPathway", 
                                   args.cutoff, args.p_value, args.adj_p_value)
    
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
    
    # Combine and process other pathway databases (KEGG, Reactome, WikiPathway)
    print("\n=== Processing combined pathway results (KEGG, Reactome, WikiPathway) ===")
    
    # Combine only non-empty dataframes
    df_list = []
    if not results_kegg.empty:
        df_list.append(results_kegg)
        print(f"KEGG: {len(results_kegg)} enriched pathways")
    else:
        print("KEGG: No enriched pathways")
    
    if not results_reactome.empty:
        df_list.append(results_reactome)
        print(f"Reactome: {len(results_reactome)} enriched pathways")
    else:
        print("Reactome: No enriched pathways")
    
    if not results_wiki.empty:
        df_list.append(results_wiki)
        print(f"WikiPathway: {len(results_wiki)} enriched pathways")
    else:
        print("WikiPathway: No enriched pathways")
    
    if df_list:
        df_combined = pd.concat(df_list, axis=0, ignore_index=True)
        print(f"Combined: {len(df_combined)} total enriched pathways")
        
        df_combined_matrix = process_results_to_matrix(df_combined)
        if not df_combined_matrix.empty:
            os.makedirs(os.path.join(args.outdir, "PASNet/Input/REST"), exist_ok=True)
            df_combined_matrix.to_excel(os.path.join(args.outdir, "PASNet/Input/REST/pt_fixed.xlsx"))
            print(f"Combined matrix saved: {df_combined_matrix.shape[0]} pathways x {df_combined_matrix.shape[1]} genes")
            
            # Save gene list to text file
            with open(os.path.join(args.outdir, "gene_list.txt"), "w") as file:
                for item in df_combined_matrix.columns:
                    file.write(item + "\n")
            print(f"Gene list written to {os.path.join(args.outdir, 'gene_list.txt')}")
        else:
            print("Combined matrix is empty after processing")
            df_combined_matrix = pd.DataFrame()
    else:
        print("WARNING: No enrichment results from KEGG, Reactome, or WikiPathway")
        df_combined_matrix = pd.DataFrame()
    
    # Generate PASNet inputs if activity and metadata are provided and we have results
    if args.activity and args.metadata:
        if not results_go_matrix.empty and not df_combined_matrix.empty:
            generate_pasnet_inputs(args.activity, args.metadata, 
                                   results_go_matrix, df_combined_matrix, args.outdir)
        else:
            missing = []
            if results_go_matrix.empty:
                missing.append("GO results")
            if df_combined_matrix.empty:
                missing.append("combined pathway results")
            print(f"\nSkipping PASNet input generation (missing: {', '.join(missing)})")
    else:
        print("\nSkipping PASNet input generation (activity and metadata not provided)")
    
    # Garbage collection
    gc.collect()
    
    print("\n=== Enrichment analysis completed successfully ===")

if __name__ == "__main__":
    main()