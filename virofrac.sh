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
    echo "  -h, --help                             Shows this helpful text :)"
    echo ""
    echo "Input options:"
    echo "  -f, --fasta [file]                     Selects FASTA file for network options [.fasta|.fna|.fa]"
    echo "  -o, --otu-table [file]                 Selects OTU table [.tsv|.csv|.tab|.tabular]"
    echo "  -t, --tax-table [file]                 Selects taxonomic table [.tsv|.csv|.tab|.tabular]"
    echo ""
    echo "Distance options:"
    echo "  -u, --unweighted-distance              Selects the unweighted distance option."
    echo "                                         • When paired with a tree option it will perform an unweighted UniFrac distance." 
    echo "                                         • When paired with the network option it will perform an unweighted NetUniFrac distance."
    echo "  -z, --unnormalized-weighted-distance   Selects the unnormalized weighted distance option."
    echo "  -w, --normalized-weighted-distance     Selects the normalized weighted unifrac option. When paired with the network option it will perform a Weighted NetUniFrac distance."
    echo "                                         • When paired with a tree option it will perform a normalized weighted UniFrac distance." 
    echo "                                         • When paired with the network option it will perform an weighted NetUniFrac distance."
    echo "  -s, --spectre                          Selects the Spectre distance (Only available with Network option)"
    echo "  -y, --fuzzy-spectre                    Selects the Fuzzy Spectre distance for abundance influence on community distances (Only available with Network option)"
    echo ""
    echo "Tree options:"
    echo "  -x, --tax-tree                         Selects the taxonomic tree mode"
    echo "  -i, --ani-tree [file]                  Selects the ANI tree option"
    echo ""
    echo "Network options:"
    echo "  -n, --network                          Selects network option"
    echo "  -a, --ani                              Selects ANI based network option (works wth vClust)"
    echo "  -g, --gene-sharing                     Selects gene-sharing based network option (works wth DIAMOND + PyRodigal)"
    echo "  -b, --threshold [value]                Selects a threshold for the network-based clustering. Default: 0"
    echo ""
    echo "Metadata/Plot options:"
    echo "  -m, --metadata [file]                  A metadata file used for the final heatmap plot. If not used, it will output a default plot. [.tsv|.csv|.tab|.tabular]"
    echo "  -l, --legend-column [column_name]      This will apply color strips and a legend to the final heatmap plot. Requires a specified color column."
    echo "  -c, --color-column [column_name]       If used, the script will apply the colors in the selected column. Requires a specified legend column."
    echo "                                         NOTE: Specify the legend columns and then the color columns. The script will use the corresponding order."
    echo "  -r, --color-gradient [string]          Default: 'coolwarm'. When used, user can specify another gradient or a custom gradient using hexacode"
    echo "  -e, --env-column [column_name]         Environmental variable(s) for envfit projection in PCoA."
    echo ""
    echo "Output options:"
    echo "  -d, --output-directory [directory]     Defines an output directory."
    echo ""
    exit 0
fi

# Tree multiple flags detection
TREE_TAX_USED=false
TREE_ANI_USED=false

# UniFrac flags detection
UNIFRAC_UU_USED=false
UNIFRAC_UW_USED=false
UNIFRAC_NW_USED=false
SPECTRAL_CLUSTERING_USED=false
FUZZY_SPECTRAL_CLUSTERING_USED=false

# Arrays for metadata
LEGEND_COLUMNS=()
COLOR_COLUMNS=()
ENV_COLUMNS=()

# Network flags detection
TREE_NET_USED=false
ANI_USED=false
GENE_SHARING_USED=false
NETWORK_THRESHOLD="0"
THRESHOLD_USED=false
OUTPUT_DIR=""

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
        -x|--tax-tree)
            TREE_TAX_USED=true
            TREE_TYPE="taxonomic"
            shift
            ;;
        -i|--ani-tree)
            TREE_ANI_USED=true
            TREE_TYPE="ani"
            shift
            ;;
        # UniFrac/Distance metrics flags
        -u|--unweighted-unifrac)
            UNIFRAC_UU_USED=true
            DISTANCE_TYPE="unweighted unifrac"
            shift
            ;;
        -z|--unnormalized-weighted-unifrac)
            UNIFRAC_UW_USED=true
            DISTANCE_TYPE="unnormalized weighted unifrac"
            shift
            ;;
        -w|--normalized-weighted-unifrac)
            UNIFRAC_NW_USED=true
            DISTANCE_TYPE="normalized weighted unifrac"
            shift
            ;;
        -s|--spectre)
            SPECTRAL_CLUSTERING_USED=true
            DISTANCE_TYPE="spectre"
            shift
            ;;
        -y|--fuzzy-spectre)
            FUZZY_SPECTRAL_CLUSTERING_USED=true
            DISTANCE_TYPE="fuzzy spectre"
            shift
            ;;
        # Network flags
        -n|--network)
            TREE_NET_USED=true
            TREE_TYPE="network"
            shift
            ;;
        -a|--ani)
            ANI_USED=true
            NETWORK_METHOD="ani"
            shift
            ;;
        -g|--gene-sharing)
            GENE_SHARING_USED=true
            NETWORK_METHOD="gene-sharing"
            shift
            ;;
        -b|--threshold)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --threshold requires a numeric value"
                exit 1
            fi
            NETWORK_THRESHOLD="$2"
            THRESHOLD_USED=true
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
        -r|--color-gradient)
            COLOR_GRADIENT="$2"
            shift 2
            ;;
        -e|--env-column)
            if [[ -z "${2:-}" ]]; then
                echo "Error: -e/--env-column requires a column name"
                exit 1
            fi
            ENV_COLUMNS+=("$2")
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
            echo "  -h, --help                             Shows this helpful text :)"
            echo ""
            echo "Input options:"
            echo "  -f, --fasta [file]                     Selects FASTA file for network options [.fasta|.fna|.fa]"
            echo "  -o, --otu-table [file]                 Selects OTU table [.tsv|.csv|.tab|.tabular]"
            echo "  -t, --tax-table [file]                 Selects taxonomic table [.tsv|.csv|.tab|.tabular]"
            echo ""
            echo "Distance options:"
            echo "  -u, --unweighted-distance              Selects the unweighted distance option."
            echo "                                         • When paired with a tree option it will perform an unweighted UniFrac distance." 
            echo "                                         • When paired with the network option it will perform an unweighted NetUniFrac distance."
            echo "  -z, --unnormalized-weighted-distance   Selects the unnormalized weighted distance option."
            echo "  -w, --normalized-weighted-distance     Selects the normalized weighted unifrac option. When paired with the network option it will perform a Weighted NetUniFrac distance."
            echo "                                         • When paired with a tree option it will perform a normalized weighted UniFrac distance." 
            echo "                                         • When paired with the network option it will perform an weighted NetUniFrac distance."
            echo "  -s, --spectre                          Selects the Spectre distance (Only available with Network option)"
            echo "  -y, --fuzzy-spectre                    Selects the Fuzzy Spectre distance for abundance influence on community distances (Only available with Network option)"
            echo ""
            echo "Tree options:"
            echo "  -x, --tax-tree                         Selects the taxonomic tree mode"
            echo "  -i, --ani-tree [file]                  Selects the ANI tree option"
            echo ""
            echo "Network options:"
            echo "  -n, --network                          Selects network option"
            echo "  -a, --ani                              Selects ANI based network option (works wth vClust)"
            echo "  -g, --gene-sharing                     Selects gene-sharing based network option (works wth DIAMOND + PyRodigal)"
            echo "  -b, --threshold [value]                Selects a threshold for the network-based clustering. Default: 0"
            echo ""
            echo "Metadata/Plot options:"
            echo "  -m, --metadata [file]                  A metadata file used for the final heatmap plot. If not used, it will output a default plot. [.tsv|.csv|.tab|.tabular]"
            echo "  -l, --legend-column [column_name]      This will apply color strips and a legend to the final heatmap plot. Requires a specified color column."
            echo "  -c, --color-column [column_name]       If used, the script will apply the colors in the selected column. Requires a specified legend column."
            echo "                                         NOTE: Specify the legend columns and then the color columns. The script will use the corresponding order."
            echo "  -r, --color-gradient [string]          Default: 'coolwarm'. When used, user can specify another gradient or a custom gradient using hexacode"
            echo "  -e, --env-column [column_name]         Environmental variable(s) for envfit projection in PCoA."
            echo ""
            echo "Output options:"
            echo "  -d, --output-directory [directory]     Defines an output directory."
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
#             TREE OPTIONS
# --------------------------------------
# A taxonomic table is required when using the tree option
if [[ "$TREE_TAX_USED" == true ]] && { [[ -z "$OTU_TABLE_FILE" ]] || [[ -z "$TAX_TABLE_FILE" ]]; }; then
    echo "Error: --tax-tree requires both --otu-table and --tax-table"
    exit 1
fi

# A otu table and FASTA file is required when using the ANI tree option
if [[ "$TREE_ANI_USED" == true ]] && { [[ -z "$OTU_TABLE_FILE" ]] }; then
    echo "Error: --ani-tree requires --otu-table"
    exit 1
fi

# --------------------------------------
#      MUTUALLY EXCLUSIVE TREE AND NETWORK OPTION
# --------------------------------------
# User will only be able to select one type of tree or the network option
TREE_OPTIONS=0
[[ "$TREE_TAX_USED" == true ]] && ((TREE_OPTIONS++))
[[ "$TREE_ANI_USED" == true ]] && ((TREE_OPTIONS++))
[[ "$TREE_NET_USED" == true ]] && ((TREE_OPTIONS++))

if [[ $TREE_OPTIONS -gt 1 ]]; then
    echo "Error: Only one tree type can be selected"
    echo ""
    echo "You used multiple tree options:"
    [[ "$TREE_TAX_USED" == true ]] && echo "  • --tax-tree (taxonomic)"
    [[ "$TREE_ANI_USED" == true ]] && echo "  • --ani-tree (ANI-based tree)"
    [[ "$TREE_NET_USED" == true ]] && echo "  • --network (network)"
    echo ""
    echo "Please select only ONE tree type."
    exit 1
fi

# --------------------------------------
#      MUTUALLY DISTANCE TYPES
# --------------------------------------
DISTANCE_OPTIONS=0
[[ "$UNIFRAC_UU_USED" == true ]] && ((DISTANCE_OPTIONS++))
[[ "$UNIFRAC_UW_USED" == true ]] && ((DISTANCE_OPTIONS++))
[[ "$UNIFRAC_NW_USED" == true ]] && ((DISTANCE_OPTIONS++))
[[ "$SPECTRAL_CLUSTERING_USED" == true ]] && ((DISTANCE_OPTIONS++))
[[ "$FUZZY_SPECTRAL_CLUSTERING_USED" == true ]] && ((DISTANCE_OPTIONS++))

if [[ $DISTANCE_OPTIONS -gt 1 ]]; then
    echo "Error: Only one distance type can be selected"
    echo ""
    echo "You used multiple distance options:"
    [[ "$UNIFRAC_UU_USED" == true ]] && echo "  • -u/--unweighted-unifrac"
    [[ "$UNIFRAC_UW_USED" == true ]] && echo "  • -z/--unnormalized-weighted-unifrac"
    [[ "$UNIFRAC_NW_USED" == true ]] && echo "  • -w/--normalized-weighted-unifrac"
    [[ "$SPECTRAL_CLUSTERING_USED" == true ]] && echo "  • -s/--spectre"
    [[ "$FUZZY_SPECTRAL_CLUSTERING_USED" == true ]] && echo "  • -y/--fuzzy-spectre"
    echo ""
    echo "Please select only ONE distance type."
    exit 1
fi

if [[ $DISTANCE_OPTIONS -eq 0 ]]; then
    echo "Error: Must provide a distance type"
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
if [[ "$THRESHOLD_USED" == true ]] && [[ "$TREE_NET_USED" != true ]]; then
    echo "Error: --threshold requires -nt/--network"
    echo "Usage: -nt --ani --threshold [value]"
    exit 1
fi

# Spectral clustering only makes sense with network mode
if [[ "$SPECTRAL_CLUSTERING_USED" == true ]] && [[ "$TREE_NET_USED" != true ]]; then
    echo "Error: --spectre requires -nt/--network"
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
if [[ -n "$NETWORK_THRESHOLD" ]] && [[ "$GENE_SHARING_USED" == true ]]; then
    echo "Warning: --threshold is ignored with --gene-sharing method"
    echo "         Threshold only applies to --ani method"
fi

# Create output directory if using ANI or gene sharing network method
if [[ "$ANI_USED" == true ]] || [[ "$GENE_SHARING_USED" == true ]]; then
    OUTPUT_DIR="$(pwd)/network_output"
    
    if [[ -d "$OUTPUT_DIR" ]]; then
        echo "Warning: Output directory already exists. Cleaning up..."
        rm -rf "$OUTPUT_DIR"
    fi
    
    mkdir -p "$OUTPUT_DIR"
    
    if [[ ! -d "$OUTPUT_DIR" ]]; then
        echo "Error: Failed to create output directory"
        exit 1
    fi
    
    echo "Created output directory: $OUTPUT_DIR"
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
        *.fasta|*.fna|*.fa)
            ;;
        *)
            echo "Error: Fasta file doesn't have expected extension"
            echo "Expected: .fasta, .fa, or .fna"
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

if [[ -n "$TAX_TABLE_FILE" ]]; then
    python3 -m src.file_validation --tax-table "$TAX_TABLE_FILE" || exit 1
fi

if [[ -n "$FASTA_FILE" ]]; then
    python3 -m src.file_validation --fasta "$FASTA_FILE" || exit 1
fi

if [[ "$TREE_TYPE" == 'phylogenetic' ]] && [[ -n "$PHY_TREE_FILE" ]]; then
    python3 -m src.file_validation --phylo-tree "$PHY_TREE_FILE" || exit 1
fi

# --------------------------------------
#     PYTHON - WORKFLOW STARTS
# --------------------------------------

PYTHON_CMD="python3 -m src.virofrac_main"

# Handle inputs based on tree type
if [[ "$TREE_TYPE" == "network" ]]; then
    # Network mode: needs OTU table and fasta
    if [[ -n "$OTU_TABLE_FILE" ]]; then
        PYTHON_CMD="$PYTHON_CMD --otu-table \"$OTU_TABLE_FILE\""
    fi
    if [[ -n "$FASTA_FILE" ]]; then
        PYTHON_CMD="$PYTHON_CMD --fasta \"$FASTA_FILE\""
    fi
else
    # Taxonomic tree: needs OTU table and tax table
    if [[ -n "$OTU_TABLE_FILE" ]] && [[ -n "$TAX_TABLE_FILE" ]]; then
        PYTHON_CMD="$PYTHON_CMD --otu-table \"$OTU_TABLE_FILE\" --tax-table \"$TAX_TABLE_FILE\""
    fi
    # ANI tree: needs OTU table and FASTA file
    if [[ -n "$OTU_TABLE_FILE" ]] && [[ -n "$FASTA_FILE" ]]; then
        PYTHON_CMD="$PYTHON_CMD --otu-table \"$OTU_TABLE_FILE\" --fasta \"$FASTA_FILE\""
    fi
fi

# Tree type
if [[ -n "$TREE_TYPE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --tree-type \"$TREE_TYPE\""
    
    if [[ "$TREE_TYPE" == "ani" ]]; then
        OUTPUT_DIR="$(pwd)/ani_tree_output"
        mkdir -p "$OUTPUT_DIR"
        PYTHON_CMD="$PYTHON_CMD --output-dir \"$OUTPUT_DIR\""
    fi
fi

# Distance type
if [[ -n "$DISTANCE_TYPE" ]]; then
    PYTHON_CMD="$PYTHON_CMD --distance-type \"$DISTANCE_TYPE\""
fi

# Network options
if [[ "$TREE_NET_USED" == true ]]; then
    if [[ -n "$NETWORK_METHOD" ]]; then
        PYTHON_CMD="$PYTHON_CMD --network-method \"$NETWORK_METHOD\""
    fi
    
    if [[ "$ANI_USED" == true ]]; then
        PYTHON_CMD="$PYTHON_CMD --threshold \"$NETWORK_THRESHOLD\""
        PYTHON_CMD="$PYTHON_CMD --output-dir \"$OUTPUT_DIR\""
    fi

    if [[ "$GENE_SHARING_USED" == true ]]; then
        #PYTHON_CMD="$PYTHON_CMD --threshold \"$NETWORK_THRESHOLD\""
        PYTHON_CMD="$PYTHON_CMD --output-dir \"$OUTPUT_DIR\""
    fi
fi

# Metadata and heatmap options
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

if [[ -n "$COLOR_GRADIENT" ]]; then
    PYTHON_CMD="$PYTHON_CMD --color-gradient \"$COLOR_GRADIENT\""
fi

if [[ ${#ENV_COLUMNS[@]} -gt 0 ]]; then
    for col in "${ENV_COLUMNS[@]}"; do
        PYTHON_CMD="$PYTHON_CMD --env-column \"$col\""
    done
fi

eval $PYTHON_CMD