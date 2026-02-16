# <img src="./img/virofrac_min_logo.png" alt="ViroFrac logo" width="40" /> ViroFrac

The ViroFrac pipeline is a software package for estimating distances between viral communities using UniFrac-based approaches.

## Features

### Relatedness methods

ViroFrac offers 3 approaches for establihing relationships between viral OTUs:

#### Taxonomic tree

Builds a hierarchical tree based on the ICTV taxonomy. For this, the tool uses a background reference tree containing all the ICTV-recognized taxa (v. ).

#### Network

Builds a similarity network between contigs. Networks can be built using vClust (v. ) and vConTACT (v. ).

#### Custom tree

### Output

ViroFrac generates outputs for comprehensive analysis and visualization:

#### Matrix

Distance values are exported in a matrix output as a tabular-separated format (`matrix_as_dataframe.tsv`). This matrix can be imported into other analyticas or visualization packaged for downstream analysis.

**Format:**
```
        Sample1    Sample2    Sample3
Sample1 0.000      0.245      0.678
Sample2 0.245      0.000      0.543
Sample3 0.678      0.543      0.000
```

#### Tree/Network

When the taxonomic tree option is selected, ViroFrac outputs the pruned phylogenetic tree used for calculations (`otu_tree.newick`). This tree is derived from the ICTV hierarchical classification and only includes the viral taxa present in your dataset.

#### Heatmap

A heatmap is the graphical output provided by ViroFrac. The heatmap function uses the seaborn package for visualization.

## Installation

## Usage
```
bash virofrac.sh --fasta uhgv/subset_votus_renamed.fna --otu-table uhgv/subset_count_votus.tsv --network --ani --threshold 0 -uu -m uhgv/country_metadata.tsv -l continent -c color_continent -g viridis
```

## Documentation

## Citation
