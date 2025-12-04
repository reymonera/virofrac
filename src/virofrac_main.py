import argparse
from pathlib import Path
import src.tree_control as treetr
import pandas as pd
import src.output_control as oc
#Provisional
from src.utils import GlobalTimer

def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input files per argument
    parser.add_argument(
        '--otu-table',
        type=Path
    )

    parser.add_argument(
        '--fasta',
        type=Path
    )

    parser.add_argument(
        '--tax-table',
        type=Path
    )

    parser.add_argument(
        '--tree-file',
        type=Path
    )

    parser.add_argument(
        '--tree-type',
        type=str
    )

    parser.add_argument(
        '--unifrac-type',
        type=str
    )

    # Network options
    parser.add_argument(
        '--network-method',
        type=str,
        choices=['ani', 'gene-sharing'],
        required=False
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.70,
        required=False
    )
    parser.add_argument(
        '--vclust-output-dir',
        type=Path,
        required=False
    )

    # Manage metadata file and its options
    parser.add_argument(
        '--metadata', 
        type=str, 
        required=False
        )

    parser.add_argument('--legend-column',
                   action='append',
                   type=str,
                   required=False)

    parser.add_argument('--color-column',
                    action='append',
                    type=str,
                    required=False)
    
    parser.add_argument('--color-gradient',
                   type=str)

    args = parser.parse_args()

    return args

def get_input_otu_table():
    args = parse_arguments()
    input_otu = pd.read_csv(args.otu_tablesharing, sep=None, engine='python')

    #print("✅ OTU Table correctly loaded")
    GlobalTimer.log("✅ OTU Table correctly loaded")

    return input_otu

def get_input_tax_table():
    args = parse_arguments()
    input_tax = pd.read_csv(args.tax_table, sep=None, engine='python',
                            na_values=['', ' ', 'NA', 'N/A', 'nan', 'NaN', 'null'], 
                            keep_default_na=True)
    
    input_tax = input_tax.replace(r'^\s*$', pd.NA, regex=True)

    #print("✅ Taxonomic Table correctly loaded")
    GlobalTimer.log("✅ Taxonomic Table correctly loaded")

    return input_tax

def get_input_otu_tree():
    args = parse_arguments()
    input_tree = None

    if args.tree_type.strip() == 'taxonomic':
        otu_tax = get_input_tax_table()
        input_tree = treetr.get_otu_tree(otu_tax)
    
    #print("✅ Tree correctly loaded")
    GlobalTimer.log("✅ Tree correctly loaded")

    return input_tree

def get_input_fasta_file():
    args = parse_arguments()
    fasta_file = args.fasta
    
    #print("✅ Tree correctly loaded")
    GlobalTimer.log("✅ Fasta file correctly loaded")

    return fasta_file

def main():
    #print("Managing inputs...")
    GlobalTimer.log("Managing inputs...")

    args = parse_arguments()
    otu_table =  get_input_otu_table()
    
    if args.tree_type.strip() == 'phylogenetic' or args.tree_type.strip() == 'taxonomic':
        GlobalTimer.log("Implementing the tree option...")
        otu_tree = get_input_otu_tree()
        matrix = oc.get_frac_matrix_output(otu_table, otu_tree, args.unifrac_type)

    elif args.tree_type.strip() == 'network':
        if not args.network_method:
            raise ValueError("Network tree type requires --network-method (ani or gene-sharing)")
        GlobalTimer.log("Implementing the network option...")

        fasta_file = get_input_fasta_file()

        if args.network_method == 'ani':
            matrix = oc.get_net_frac_matrix_output(
                otu_table, 
                fasta_file, 
                args.unifrac_type,
                args.network_method,
                args.threshold,
                args.vclust_output_dir
            )
        elif args.network_method == 'gene-sharing':
            matrix = oc.get_net_frac_matrix_output(
                otu_table, 
                fasta_file, 
                args.unifrac_type,
                args.network_method
            )

    oc.get_heatmap_output(
        matrix, 
        otu_table, 
        args.unifrac_type,
        args.color_gradient,
        args.metadata,
        args.legend_column,
        args.color_column)

    return 0

if __name__ == "__main__":
    main()
