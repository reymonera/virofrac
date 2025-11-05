import sys
import pandas as pd
import argparse

# This function should be the one putting the taxa in a
# browsable list for the tree prune function to work.
def get_validation_otu_table(filepath):
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
        
        if '\t' in first_line:
            separator = '\t'
            formato = "TSV"
        elif ',' in first_line:
            separator = ','
            formato = "CSV"
        else:
            raise ValueError("Cannot detect separator (expected tab or comma)")
        
        df = pd.read_csv(filepath, sep=separator)
        otu_col = None
        
        otu_col_candidates = ['#OTU ID', 'OTU ID', 'OTU_ID', '#OTU_ID', 'OTUID']
        for candidate in otu_col_candidates:
            if candidate in df.columns:
                otu_col = candidate
                break
        
        if otu_col is None:
            raise ValueError(f"OTU table missing OTU ID column. Expected one of: {otu_col_candidates}")

        sample_cols = [col for col in df.columns if col != otu_col]
        if len(sample_cols) <= 1:
            raise ValueError("OTU table has no sample columns")
        
        df = pd.read_csv(filepath, sep=separator)
        print(f"✅ OTU table validated")
    
    except Exception as e:
        print(f"❌ Error validation OTU Table: {e}")
        sys.exit(1)

def get_validation_tax_table(filepath):
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
        
        if '\t' in first_line:
            separator = '\t'
            formato = "TSV"
        elif ',' in first_line:
            separator = ','
            formato = "CSV"
        else:
            raise ValueError("Cannot detect separator (expected tab or comma)")
        
        df = pd.read_csv(filepath, sep=separator)
        otu_col = None
        
        otu_col_candidates = ['#OTU ID', 'OTU ID', 'OTU_ID', '#OTU_ID', 'OTUID']
        for candidate in otu_col_candidates:
            if candidate in df.columns:
                otu_col = candidate
                break
        
        if otu_col is None:
            raise ValueError(f"Tax table missing OTU ID column. Expected one of: {otu_col_candidates}")

        sample_cols = [col for col in df.columns if col != otu_col]
        if len(sample_cols) <= 1:
            raise ValueError("Tax table has no taxonomic columns")
        
        df = pd.read_csv(filepath, sep=separator)
        print(f"✅ Tax table validated")
    
    except Exception as e:
        print(f"❌ Error validation Tax Table: {e}")
        sys.exit(1)

def get_validation_reads_file(filepath):
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
        
        if first_line.startswith("@") == False:
            raise ValueError(".fastq file is not valid")
        
        print(f"✅ Reads file validated")
        
    except Exception as e:
        print(f"❌ Error validation reads file: {e}")
        sys.exit(1)

def get_validation_newick_file(filepath):
    try:
        with open(filepath, 'r') as f:
            tree_string = f.read().strip()
        
        if not tree_string.endswith(';'):
            raise ValueError("Newick string must end with semicolon (;)")
        
        open_parens = tree_string.count('(')
        close_parens = tree_string.count(')')
        
        if open_parens != close_parens:
            raise ValueError(
                f"Unbalanced parentheses: {open_parens} '(' vs {close_parens} ')'"
            )
        
        print(f"✅ Newick file validated")
        
    except Exception as e:
        print(f"❌ Error validation newick file: {e}")
        sys.exit(1)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--otu-table',
        type=str
    )
    
    parser.add_argument(
        '--reads',
        type=str
    )
    
    parser.add_argument(
        '--phylo-tree',
        type=str
    )

    parser.add_argument(
        '--tax-table',
        type=str
    )
    
    arguments = parser.parse_args()
    
    if arguments.otu_table:
        get_validation_otu_table(arguments.otu_table)
    
    if arguments.tax_table:
        get_validation_tax_table(arguments.tax_table)
    
    if arguments.reads:
        get_validation_reads_file(arguments.reads)
    
    if arguments.phylo_tree:
        get_validation_newick_file(arguments.phylo_tree)
        
    if not any([arguments.otu_table, arguments.tax_table, arguments.reads, arguments.phylo_tree]):
        sys.exit(1)
        
if __name__ == "__main__":
    main()