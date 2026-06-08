# <img src="./img/virofrac_min_logo.png" alt="ViroFrac logo" width="40" /> ViroFrac

ViroFrac is a software package for estimating distances between viral communities using genome similarity information.

## Features

### Relatedness methods

ViroFrac offers 2 approaches for establishing relationships between viral OTUs:

#### Taxonomic tree

Builds a hierarchical tree based on the ICTV taxonomy. For this, the tool uses a background reference tree containing all the ICTV-recognized taxa (VMR MSL40).

#### Network

Builds a similarity network between contigs. Networks can be built using ANI-based distance with vClust (v1.3.1) and gene-sharing based distance with Pyrodigal (v3.7.0) and DIAMOND (v2.1.15.169).

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

#### Heatmap

A heatmap is the graphical output provided by ViroFrac. The heatmap function uses the seaborn package for visualization.

#### PCoA

A PCoA plot is also provided as a standard ViroFrac output. If included, an envfit vector can be drawn on the PCoA plots.

## Installation

## Usage
```
bash virofrac.sh --fasta uhgv/subset_votus_renamed.fna --otu-table uhgv/subset_count_votus.tsv --network --ani --threshold 0 -u -m uhgv/country_metadata.tsv -l continent -c color_continent -g viridis
```

## Documentation

## Citation
