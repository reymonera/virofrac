#!/bin/bash

# --------------------------------------
#             PARAMETERS
# --------------------------------------
if [[ $# -eq 0 ]]; then
    echo ' _  _  __  ___   __  ___  ___    __   __ '
    echo '( )( )(  )(  ,) /  \(  _)(  ,)  (  ) / _)'
    echo ' \\//  )(  )  \( () )) _) )  \  /__\( (_ '
    echo ' (__) (__)(_)\_)\__/(_)  (_)\_)(_)(_)\__)'
    echo ""
    echo "Virofrac v. 1.0.0"
    echo "Usage: $(basename "$0") [options]"
    echo ""
    echo "Made with elegance 🍷"
    echo ""
    echo "General options:"
    echo "  -h, --help                              Shows this helpful text :)"
    echo ""
    echo "Input options:"
    echo "  -f, --fasta [file]                      Selects .fasta file"
    echo "  -o, --otu-table [file]                  Selects OTU table [.tsv|.csv|.tab|.tabular]"
    echo "  -t, --tax-table [file]                  Selects taxonomic table [.tsv|.csv|.tab|.tabular]"
    echo ""
    echo "UniFrac distance options:"
    echo "  -uu, --unweighted-unifrac               Selects the unweighted unifrac option"
    echo "  -uw, --unnormalized-weighted-unifrac    Selects the unnormalized weighted unifrac option"
    echo "  -nw, --normalized-weighted-unifrac      Selects the normalized weighted unifrac option"
    echo ""
    echo "Tree options:"
    echo "  -tt, --tax-tree                         Selects the taxonomic tree mode"
    echo "  -pt, --phy-tree [file]                  Selects the phylogenetic tree option [.newick]"
    echo ""
    echo "Network options:"
    echo "  -nt, --network                          Selects network option"
    echo "  --ani                                   Selects ANI based network option (works wth vclust)"
    echo "  --gene-sharing                          Selects gene-sharing based network option (works wth vcontact3)"
    echo "  --threshold [value]                     Selects a threshold (XXX!!!!!) for the network-based clustering. Default: 0.70"
    echo ""
    echo "Metadata/Plot options:"
    echo "  -m, --metadata [file]                   A metadata file used for the final heatmap plot. If not used, it will output a default plot. [.tsv|.csv|.tab|.tabular]"
    echo "  -l, --legend-column [column_name]       This will apply color strips and a legend to the final heatmap plot. Requires a specified color column."
    echo "  -c, --color-column [column_name]        If used, the script will apply the colors in the selected column. Requires a specified legend column."
    echo "                                          NOTE: Specify the legend columns and then the color columns. The script will use the corresponding order."
    echo "  -g, --color-gradient [string]           Deafult: 'coolwarm'. When used, user can specify another gradient or a custom gradient using hexacode"
    echo ""
    exit 0
fi

# Tree multiple flags detection
TREE_TAX_USED=false
TREE_PHY_USED=false

# UniFrac flags detection
UNIFRAC_UU_USED=false
UNIFRAC_UW_USED=false
UNIFRAC_NW_USED=false

# Arrays for metadata
LEGEND_COLUMNS=()
COLOR_COLUMNS=()

# Network flags detection
TREE_NET_USED=false
ANI_USED=false
GENE_SHARING_USED=false
NETWORK_THRESHOLD="0.70"
VCLUST_OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        # Input flags
        -o|--otu-table)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -o/--otu-table requires an OTU table"
                exit 1
            fi
            OTU_TABLE_FILE="$2"
            shift 2
            ;;
        -t|--tax-table)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -t/--tax-table requires a taxonomic table"
                exit 1
            fi
            TAX_TABLE_FILE="$2"
            shift 2
            ;;
        -f|--fasta)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -f/--fasta requires fasta file"
                exit 1
            fi
            FASTA_FILE="$2"
            shift 2
            ;;
        -tt|--tax-tree)
            TREE_TAX_USED=true
            TREE_TYPE="taxonomic"
            shift
            ;;
        -pt|--phy-tree)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -pt/--phy_tree requires a phylogenetic tree file"
                exit 1
            fi
            PHY_TREE_FILE="$2"
            TREE_PHY_USED=true
            TREE_TYPE="phylogenetic"
            shift 2
            ;;
        # UniFrac/Distance metrics flags
        -uu|--unweighted-unifrac)
            UNIFRAC_UU_USED=true
            UNIFRAC_TYPE="unweighted unifrac"
            shift
            ;;
        -uw|--unnormalized-weighted-unifrac)
            UNIFRAC_UW_USED=true
            UNIFRAC_TYPE="unnormalized weighted unifrac"
            shift
            ;;
        -nw|--normalized-weighted-unifrac)
            UNIFRAC_NW_USED=true
            UNIFRAC_TYPE="normalized weighted unifrac"
            shift
            ;;
        # Network flags
        -nt|--network)
            TREE_NET_USED=true
            TREE_TYPE="network"
            shift
            ;;
        --ani)
            ANI_USED=true
            NETWORK_METHOD="ani"
            shift
            ;;
        --gene-sharing)
            GENE_SHARING_USED=true
            NETWORK_METHOD="gene-sharing"
            shift
            ;;
        --threshold)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --threshold requires a numeric value"
                exit 1
            fi
            NETWORK_THRESHOLD="$2"
            shift 2
            ;;
        # Metadata - Plot Output flags
        -m|--metadata)
            METADATA_FILE="$2"
            shift 2
            ;;
        -l|--legend-column)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -l/--legend-column requires a column name"
                exit 1
            fi
            LEGEND_COLUMNS+=("$2")
            shift 2
            ;;

        -c|--color-column)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -c/--color-column requires a column name"
                exit 1
            fi
            COLOR_COLUMNS+=("$2")
            shift 2
            ;;
        -g|--color-gradient)
            COLOR_GRADIENT="$2"
            shift 2
            ;;
        -h|--help)
                echo ' _  _  __  ___   __  ___  ___    __   __ '
                echo '( )( )(  )(  ,) /  \(  _)(  ,)  (  ) / _)'
                echo ' \\//  )(  )  \( () )) _) )  \  /__\( (_ '
                echo ' (__) (__)(_)\_)\__/(_)  (_)\_)(_)(_)\__)'
                echo ""
                echo "Virofrac v. 1.0.0"
                echo "Usage: $(basename "$0") [options]"
                echo ""
                echo "Made with elegance 🍷"
                echo ""
                echo "General options:"
                echo "  -h, --help                              Shows this helpful text :)"
                echo ""
                echo "Input options:"
                echo "  -f, --fasta [file]                      Selects .fasta file"
                echo "  -o, --otu-table [file]                  Selects OTU table [.tsv|.csv|.tab|.tabular]"
                echo "  -t, --tax-table [file]                  Selects taxonomic table [.tsv|.csv|.tab|.tabular]"
                echo ""
                echo "UniFrac distance options:"
                echo "  -uu, --unweighted-unifrac               Selects the unweighted unifrac option"
                echo "  -uw, --unnormalized-weighted-unifrac    Selects the unnormalized weighted unifrac option"
                echo "  -nw, --normalized-weighted-unifrac      Selects the normalized weighted unifrac option"
                echo ""
                echo "Tree options:"
                echo "  -tt, --tax-tree                         Selects the taxonomic tree mode"
                echo "  -pt, --phy-tree [file]                  Selects the phylogenetic tree option [.newick]"
                echo ""
                echo "Network options:"
                echo "  -nt, --network                          Selects network option"
                echo "  --ani                                   Selects ANI based network option (works wth vclust)"
                echo "  --gene-sharing                          Selects gene-sharing based network option (works wth vcontact3)"
                echo "  --threshold [value]                     Selects a threshold (XXX!!!!!) for the network-based clustering. Default: 0.70"
                echo ""
                echo "Metadata/Plot options:"
                echo "  -m, --metadata [file]                   A metadata file used for the final heatmap plot. If not used, it will output a default plot. [.tsv|.csv|.tab|.tabular]"
                echo "  -l, --legend-column [column_name]       This will apply color strips and a legend to the final heatmap plot. Requires a specified color column."
                echo "  -c, --color-column [column_name]        If used, the script will apply the colors in the selected column. Requires a specified legend column."
                echo "                                          NOTE: Specify the legend columns and then the color columns. The script will use the corresponding order."
                echo "  -g, --color-gradient [string]           Deafult: 'coolwarm'. When used, user can specify another gradient or a custom gradient using hexacode"
                echo ""
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'"
            echo "Use -h/--help to get the help you need, traveler!"
            exit 1
            ;;
    esac
done

# --------------------------------------
#             CITATION
# --------------------------------------
echo "============================================================"
echo "Thank you for using Virofrac v.1.0.0!"
echo ""
echo "If you use this tool in your work, please cite as following:"
echo "(Also never cite that stupid market chicken paper, please)"
echo ""
echo "Castillo-Vilcahuaman et al. Virofrac loves Timoteo"
echo "The UDC Journal, 25(16):2078–2079 (2026)." 
echo "doi:10.1093/judc/btp352"
echo "============================================================"

# Starts validation with Bash.
echo "Checking if files are valid..."

# --------------------------------------
#    TAX-TABLE REQUIRED WITH OTU-TABLE
# --------------------------------------
if [[ -z "$OTU_TABLE_FILE" ]] && [[ -z "$TAX_TABLE_FILE" ]]; then
    echo "Error: --tax-table is required when using --otu-table"
    echo "Usage: virofrac.sh --otu-table [file] --tax-table [file] --[tree option]"
    exit 1
fi

if [[ -z "$OTU_TABLE_FILE" ]] && [[ -z "$TAX_TABLE_FILE" ]]; then
    echo "Warning: --tax-table provided but --otu-table not selected"
    echo "         --tax-table is only used with --otu-table"
fi

# # --------------------------------------
# #    CHECK IF OTU_TABLE OR READS
# # --------------------------------------

# INPUT_COUNT=0
# [[ -n "$OTU_TABLE_FILE" ]] && ((INPUT_COUNT++))
# [[ -n "$FASTA_FILE" ]] && ((INPUT_COUNT++))

# if [[ $INPUT_COUNT -eq 0 ]]; then
#     echo "Error: Must provide either -o/--otu-table OR -f/--fasta"
#     exit 1#17becf
# fi

# if [[ $INPUT_COUNT -gt 1 ]]; then
#     echo "Error: Cannot use both -o/--otu-table and -f/--fasta together"
#     echo "       Please choose only one input type"
#     exit 1
# fi

# --------------------------------------
#      MUTUALLY EXCLUSIVE TREES
# --------------------------------------
# User will only be able to select one type of tree.
TREE_OPTIONS=0
[[ "$TREE_TAX_USED" == true ]] && ((TREE_OPTIONS++))
[[ "$TREE_PHY_USED" == true ]] && ((TREE_OPTIONS++))
[[ "$TREE_NET_USED" == true ]] && ((TREE_OPTIONS++))

if [[ $TREE_OPTIONS -gt 1 ]]; then
    echo "Error: Only one tree type can be selected"
    echo ""
    echo "You used multiple tree options:"
    [[ "$TREE_TAX_USED" == true ]] && echo "  • -tt/--tax-tree (taxonomic)"
    [[ "$TREE_PHY_USED" == true ]] && echo "  • -uw/--phy-tree (phylogenetic)"
    [[ "$TREE_NET_USED" == true ]] && echo "  • -nw/--network (network)"
    echo ""
    echo "Please select only ONE tree type."
    exit 1
fi

# --------------------------------------
#      MUTUALLY UNIFRAC TYPES
# --------------------------------------
UNIFRAC_OPTIONS=0
[[ "$UNIFRAC_UU_USED" == true ]] && ((UNIFRAC_OPTIONS++))
[[ "$UNIFRAC_UW_USED" == true ]] && ((UNIFRAC_OPTIONS++))
[[ "$UNIFRAC_NW_USED" == true ]] && ((UNIFRAC_OPTIONS++))

if [[ $UNIFRAC_OPTIONS -gt 1 ]]; then
    echo "Error: Only one UniFrac type can be selected"
    echo ""
    echo "You used multiple UniFrac options:"
    [[ "$UNIFRAC_UU_USED" == true ]] && echo "  • -uu/--unweighted-unifrac"
    [[ "$UNIFRAC_UW_USED" == true ]] && echo "  • -uw/--unnormalized-weighted-unifrac"
    [[ "$UNIFRAC_NW_USED" == true ]] && echo "  • -nw/--normalized-weighted-unifrac"
    echo ""
    echo "Please select only ONE UniFrac type."
    exit 1
fi

if [[ $UNIFRAC_OPTIONS -eq 0 ]]; then
    echo "Error: Must provide a UniFrac distance type"
    echo ""
    echo "Choose one:"
    echo "  • -uu/--unweighted-unifrac"
    echo "  • -uw/--unnormalized-weighted-unifrac"
    echo "  • -nw/--normalized-weighted-unifrac"
    exit 1
fi

# --------------------------------------
#           HEATMAP CONTROLS
# --------------------------------------
# If the metadata option is used, then the legend column needs to be specified.
if [[ -n "$METADATA_FILE" ]] && [[ ${#LEGEND_COLUMNS[@]} -eq 0 ]]; then
    echo "Error: --metadata requires at least one --legend-column"
    echo "Usage: --metadata [file] --legend-column [column_name]"
    exit 1
fi

# Color column option requires the metadata option.
if [[ ${#COLOR_COLUMNS[@]} -gt 0 ]] && [[ -z "$METADATA_FILE" ]]; then
    echo "Error: --color-column requires --metadata"
    echo "Usage: --metadata [file] --legend-column [column] --color-column [column]"
    exit 1
fi

# Same number of color columns and legend columns.
if [[ ${#COLOR_COLUMNS[@]} -gt 0 ]] && [[ ${#COLOR_COLUMNS[@]} -ne ${#LEGEND_COLUMNS[@]} ]]; then
    echo "Error: Number of --color-column must match number of --legend-column"
    echo "You provided:"
    echo "  ${#LEGEND_COLUMNS[@]} legend columns: ${LEGEND_COLUMNS[*]}"
    echo "  ${#COLOR_COLUMNS[@]} color columns: ${COLOR_COLUMNS[*]}"
    exit 1
fi

# --------------------------------------
#           NETWORK OPTIONS
# --------------------------------------
# If network mode is selected, a network method must be chosen
if [[ "$TREE_NET_USED" == true ]]; then
    NETWORK_METHODS=0
    [[ "$ANI_USED" == true ]] && ((NETWORK_METHODS++))
    [[ "$GENE_SHARING_USED" == true ]] && ((NETWORK_METHODS++))
    
    if [[ $NETWORK_METHODS -eq 0 ]]; then
        echo "Error: --network requires a network method"
        echo ""
        echo "Choose one:"
        echo "  • --ani (ANI-based clustering with vclust)"
        echo "  • --gene-sharing (gene-sharing based clustering with vcontact3)"
        exit 1
    fi
    
    if [[ $NETWORK_METHODS -gt 1 ]]; then
        echo "Error: Only one network method can be selected"
        echo ""
        echo "You used multiple network methods:"
        [[ "$ANI_USED" == true ]] && echo "  • --ani"
        [[ "$GENE_SHARING_USED" == true ]] && echo "  • --gene-sharing"
        echo ""
        echo "Please select only ONE network method."
        exit 1
    fi
    
    echo "Using network threshold: $NETWORK_THRESHOLD"
fi

# Threshold only makes sense with network mode
if [[ "$NETWORK_THRESHOLD" != "0.70" ]] && [[ "$TREE_NET_USED" != true ]]; then
    echo "Error: --threshold requires -nt/--network"
    echo "Usage: -nt --ani --threshold [value]"
    exit 1
fi

# Validate threshold is a valid number between 0 and 1
if ! [[ "$NETWORK_THRESHOLD" =~ ^[0-1]?\.?[0-9]+$ ]]; then
    echo "Error: --threshold must be a number between 0 and 1"
    echo "Got: $NETWORK_THRESHOLD"
    exit 1
fi

# Additional range check using bc
if (( $(echo "$NETWORK_THRESHOLD < 0 || $NETWORK_THRESHOLD > 1" | bc -l) )); then
    echo "Error: --threshold must be between 0 and 1"
    echo "Got: $NETWORK_THRESHOLD"
    exit 1
fi

# ANI and gene-sharing shouldn't be used without network flag
if [[ "$ANI_USED" == true || "$GENE_SHARING_USED" == true ]] && [[ "$TREE_NET_USED" != true ]]; then
    echo "Error: --ani and --gene-sharing require -nt/--network"
    exit 1
fi

# Threshold only applies to ANI method
if [[ "$NETWORK_THRESHOLD" != "0.70" ]] && [[ "$GENE_SHARING_USED" == true ]]; then
    echo "Warning: --threshold is ignored with --gene-sharing method"
    echo "         Threshold only applies to --ani method"
fi

# Create vclust output directory if using ANI network method
if [[ "$ANI_USED" == true ]]; then
    VCLUST_OUTPUT_DIR="$(pwd)/vclust_output"
    
    if [[ -d "$VCLUST_OUTPUT_DIR" ]]; then
        echo "Warning: vclust output directory already exists. Cleaning up..."
        rm -rf "$VCLUST_OUTPUT_DIR"
    fi
    
    mkdir -p "$VCLUST_OUTPUT_DIR"
    
    if [[ ! -d "$VCLUST_OUTPUT_DIR" ]]; then
        echo "Error: Failed to create vclust output directory"
        exit 1
    fi
    
    echo "Created vclust output directory: $VCLUST_OUTPUT_DIR"
fi

# --------------------------------------
#            MISSING FILES
# --------------------------------------
# Handling a missing OTU_TABLE_FILE
if [[ -n "$OTU_TABLE_FILE" ]]; then
    if [[ ! -f "$OTU_TABLE_FILE" ]]; then
        echo "Error: OTU table file not found: $OTU_TABLE_FILE"
        exit 1
    fi
fi

# Handling a missing FASTA_FILE
if [[ -n "$FASTA_FILE" ]]; then
    if [[ ! -f "$FASTA_FILE" ]]; then
        echo "Error: Reads file not found: $FASTA_FILE"
        exit 1
    fi
fi

# Handling a missing PHY_TREE_FILE
if [[ -n "$PHY_TREE_FILE" ]]; then
    if [[ ! -f "$PHY_TREE_FILE" ]]; then
        echo "Error: Phylogenetic tree file not found: $PHY_TREE_FILE"
        exit 1
    fi
fi

# Metadata file should exist.
if [[ -n "$METADATA_FILE" ]] && [[ ! -f "$METADATA_FILE" ]]; then
    echo "Error: Metadata file not found: $METADATA_FILE"
    exit 1
fi

# --------------------------------------
#            EMPTY FILES
# --------------------------------------
# Handling an empty OTU_TABLE_FILE
if [[ -n "$OTU_TABLE_FILE" ]]; then
    if [[ ! -s "$OTU_TABLE_FILE" ]]; then
        echo "Error: OTU table file is empty"
        exit 1
    fi
fi

# Handling an empty FASTA_FILE
if [[ -n "$FASTA_FILE" ]]; then
    if [[ ! -s "$FASTA_FILE" ]]; then
        echo "Error: Reads file is empty"
        exit 1
    fi
fi

# Handling an empty PHY_TREE_FILE
if [[ -n "$PHY_TREE_FILE" ]]; then
    if [[ ! -s "$PHY_TREE_FILE" ]]; then
        echo "Error: Phylogenetic tree file is empty"
        exit 1
    fi
fi

# Handling an empty METADATA_FILE
if [[ -n "$METADATA_FILE" ]] && [[ ! -s "$METADATA_FILE" ]]; then
    echo "Error: Metadata file is empty"
    exit 1
fi
# --------------------------------------
#            FILE EXTENSIONS
# --------------------------------------
# Handling file extensions in FASTA_FILE
if [[ -n "$FASTA_FILE" ]]; then
    case "$FASTA_FILE" in
        *.fastq|*.fq|*.fastq.gz|*.fq.gz)
            ;;
        *)
            echo "Error: Reads file doesn't have expected extension"
            echo "Expected: .fastq, .fq, .fasta, .fa, or .gz compressed versions" # Pretty sure it is just .fastq or .fq, but we haven't reached this stage yet.
            echo "Got: $FASTA_FILE"
            exit 1
            ;;
    esac
fi

# Handling file extensions in OTU_TABLE_FILE
if [[ -n "$OTU_TABLE_FILE" ]]; then
    case "$OTU_TABLE_FILE" in
        *.csv|*.tsv|*.tab|*.tabular|*.txt)
            ;;
        *)
            echo "Error: OTU table doesn't have expected extension"
            echo "Expected: .csv, .tsv, .tab, .tabular, or .txt"
            echo "Got: $OTU_TABLE_FILE"
            exit 1
            ;;
    esac
fi

# Handling file extensions in PHY_TREE_FILE
if [[ -n "$PHY_TREE_FILE" ]]; then
    case "$PHY_TREE_FILE" in
        *.newick)
            ;;
        *)
            echo "Error: Phylogenetic tree file doesn't have expected extension (.newick)"
            exit 1
            ;;
    esac
fi

# Starts validation with Python.
echo "Checking if files are correctly formatted..."

# --------------------------------------
#     PYTHON -  VALIDATION
# --------------------------------------
python3 -m src.file_validation --otu-table "$OTU_TABLE_FILE" || exit 1
python3 -m src.file_validation --tax-table "$TAX_TABLE_FILE" || exit 1

if [[ -f $FASTA_FILE ]]; then
    python3 -m src.file_validation --reads "$FASTA_FILE" || exit 1
fi

if [[ $TREE_TYPE == 'phylogenetic' ]]; then
    python3 -m src.file_validation --phylo-tree "$PHY_TREE_FILE" || exit 1
fi

# Starts validation with workflow.
echo "Starting with workflow..."

# --------------------------------------
#     PYTHON - WORKFLOW STARTS
# --------------------------------------

PYTHON_CMD="python3 -m src.virofrac_main"

# Deals with the OTU Table inpt
if [[ -n "$OTU_TABLE_FILE" ]] && [[ -n "$TAX_TABLE_FILE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --otu-table \"$OTU_TABLE_FILE\" --tax-table \"$TAX_TABLE_FILE\""
elif [[ -n "$FASTA_FILE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --reads \"$FASTA_FILE\""
fi

# Deals with the used tree
if [[ -n "$TREE_TYPE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --tree-type \"$TREE_TYPE\""
    
    if [[ "$TREE_TYPE" == "phylogenetic" ]] && [[ -n "$PHY_TREE_FILE" ]]; then
        PYTHON_CMD="$PYTHON_CMD --tree-file \"$PHY_TREE_FILE\""
    fi
fi

# Deals with the UniFrac type
if [[ -n "$UNIFRAC_TYPE" ]] && [[ -n "$OTU_TABLE_FILE" ]] && [[ -n "$TAX_TABLE_FILE" ]] && [[ -n "$TREE_TYPE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --unifrac-type \"$UNIFRAC_TYPE\""
fi

# Deals with metadata and heatmap options
if [[ -n "$METADATA_FILE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --metadata \"$METADATA_FILE\""
    
    if [[ ${#LEGEND_COLUMNS[@]} -gt 0 ]]; then
        for col in "${LEGEND_COLUMNS[@]}"; do
            PYTHON_CMD="$PYTHON_CMD --legend-column \"$col\""
        done
    fi
    
    if [[ ${#COLOR_COLUMNS[@]} -gt 0 ]]; then
        for col in "${COLOR_COLUMNS[@]}"; do
            PYTHON_CMD="$PYTHON_CMD --color-column \"$col\""
        done
    fi
fi

# Deals with network options
if [[ "$TREE_NET_USED" == true ]]; then
    if [[ -n "$NETWORK_METHOD" ]]; then
        PYTHON_CMD="$PYTHON_CMD --network-method \"$NETWORK_METHOD\""
    fi
    
    # Only pass threshold and vclust output dir for ANI method
    if [[ "$ANI_USED" == true ]]; then
        PYTHON_CMD="$PYTHON_CMD --threshold \"$NETWORK_THRESHOLD\""
        PYTHON_CMD="$PYTHON_CMD --vclust-output-dir \"$VCLUST_OUTPUT_DIR\""
    fi
fi

eval $PYTHON_CMD