import sys
import pandas as pd
import argparse
from .utils import GlobalTimer

# This function should be the one putting the taxa in a
# browsable list for the tree prune function to work. This
# also contains a normalization process, in which the sum of
# sample columns is evaluated. Every sample column should sum 1.
def get_validation_otu_table(filepath):
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if '\t' in first_line:
                separator = '\t'
            elif ',' in first_line:
                separator = ','
            else:
                raise ValueError("Cannot detect separator (expected tab or comma)")

        df = pd.read_csv(filepath, sep=separator, index_col=0)

        sample_cols = df.columns.tolist()

        if len(sample_cols) <= 1:
            raise ValueError("OTU table has no sample columns")

        # Normalization process
        GlobalTimer.log("Checking if OTU table is normalized...")
        col_sums = df[sample_cols].sum(axis=0)
        #print(col_sums)
        if not (abs(col_sums - 1.0) < 1e-6).all():
            GlobalTimer.log("WARNING: OTU table columns do not sum to 1, normalizing...")
            df[sample_cols] = df[sample_cols].div(col_sums, axis=1)
            df.to_csv(filepath, sep=separator, index=True)
            GlobalTimer.log("✓ OTU table normalized and saved")

        GlobalTimer.log("✓ OTU table validated")

    except Exception as e:
        GlobalTimer.log(f"✕ Error validation OTU Table: {e}")
        sys.exit(1)

# This function evaluates the OTU table. First, it identifies the
# first column as an id column, and then proceeds to evaluate if
# the values are empty or numeric. Then, i determines if the
# column values are numeric or non-numeric.
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
        if df.shape[1] < 2:
            raise ValueError("Tax table needs an ID column and at least one numeric column")

        otu_col = df.columns[0]
        ids = df[otu_col]

        if ids.isnull().any():
            raise ValueError(f"ID column '{otu_col}' contains empty values")
        if pd.to_numeric(ids, errors='coerce').notnull().all():
            raise ValueError(f"First column '{otu_col}' looks numeric; expected text IDs")

        for col in df.columns[1:]:
            coerced = pd.to_numeric(df[col], errors='coerce')
            bad = coerced.isnull() & df[col].notnull()
            if bad.any():
                rows = df.index[bad].tolist()[:5]
                raise ValueError(f"Column '{col}' has non-numeric values at rows {rows}")

        GlobalTimer.log("✓ Tax table validated")
    
    except Exception as e:
        GlobalTimer.log(f"✕ Error validation Tax Table: {e}")
        sys.exit(1)

def get_validation_fasta_file(filepath):
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
        
        if first_line.startswith(">") == False:
            raise ValueError("fasta file is not valid, check if all assemblies are in the same file")
        
        GlobalTimer.log("✓ Fasta file validated")
        #print("✓ Reads file validated")
        
    except Exception as e:
        GlobalTimer.log("✕ Error validation fasta file: {e}")
        #print(f"✕ Error validation reads file: {e}")
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
        
        GlobalTimer.log("✓ Newick file validated")
        #print("✓ Newick file validated")
        
    except Exception as e:
        GlobalTimer.log("✕ Error validation newick file: {e}")
        #print(f"✕ Error validation newick file: {e}")
        sys.exit(1)

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--otu-table',
        type=str
    )
    
    parser.add_argument(
        '--fasta',
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
    
    if arguments.fasta:
        get_validation_fasta_file(arguments.fasta)
    
    if arguments.phylo_tree:
        get_validation_newick_file(arguments.phylo_tree)
        
    if not any([arguments.otu_table, arguments.tax_table, arguments.fasta, arguments.phylo_tree]):
        sys.exit(1)
        
if __name__ == "__main__":
    main()