import subprocess
import networkx as nx
import numpy as np
import pandas as pd
#import glob
from src.utils import GlobalTimer
#import shutil
from pathlib import Path
import pyrodigal

# Esta función debería de acepar el network basado en ani y en gene sharing
def get_count_table(count_table_path):
    ## Creo que no necestamos esto por el tipo de input
    table = pd.read_csv(count_table_path)
    indexed_count_table = table.set_index([0])

    return indexed_count_table

##CHECK
# Convert vclust edge list output to a square adjacency matrix.
def set_edge_list_to_matrix(edges_df):
    # Detect query column
    if 'query' in edges_df.columns:
        query_col = 'query'
    else:
        raise ValueError(f"Could not find query column. Available: {edges_df.columns.tolist()}")
    
    # Detect target/reference column
    if 'reference' in edges_df.columns:
        target_col = 'reference'
    elif 'target' in edges_df.columns:
        target_col = 'target'
    else:
        raise ValueError(f"Could not find target/reference column. Available: {edges_df.columns.tolist()}")
    
    # Detect ANI column
    if 'tani' in edges_df.columns:
        ani_col = 'tani'
    elif 'gani' in edges_df.columns:
        ani_col = 'gani'
    elif 'ani' in edges_df.columns:
        ani_col = 'ani'
    else:
        raise ValueError(f"Could not find ANI column. Available: {edges_df.columns.tolist()}")

    #print(f"DEBUG: Using columns - query: {query_col}, target: {target_col}, ani: {ani_col}")

    # Get all unique sequence IDs
    all_seqs = pd.unique(edges_df[[query_col, target_col]].values.ravel())
    n = len(all_seqs)
    
    seq_to_idx = {seq: idx for idx, seq in enumerate(all_seqs)}
    
    matrix = np.zeros((n, n))
    
    for _, row in edges_df.iterrows():
        i = seq_to_idx[row[query_col]]
        j = seq_to_idx[row[target_col]]
        ani_value = float(row[ani_col])
        matrix[i, j] = ani_value
        matrix[j, i] = ani_value

    np.fill_diagonal(matrix, 1.0)
        
    return matrix

def get_ani_matrix_vclust(input_fasta, output_dir):
    # Prefilter
    GlobalTimer.log("Prefiltering sponsored by vclust...")
    cmd_prefilter = [
        'vclust', 'prefilter',
        '-i', input_fasta,
        '-o', f'{output_dir}/filter.txt'
    ]
    subprocess.run(cmd_prefilter, check=True)

    # Alignment
    GlobalTimer.log("Alignment sponsored by vclust...")
    cmd_align = [
        'vclust', 'align',
        '-i', input_fasta,
        '-o', f'{output_dir}/ani.tsv',
        '--filter', f'{output_dir}/filter.txt'
    ]
    subprocess.run(cmd_align, check=True)

    edges_df = pd.read_csv(f'{output_dir}/ani.tsv', sep='\t')
    matrix = set_edge_list_to_matrix(edges_df)

    return matrix

# Contrary to vClust, vConTACT3 will work better
# producing the network instead of doing a matrix and
# then producing the network.
# def get_gene_sharing_network_vcontact3(input_fasta, output_dir):
#     GlobalTimer.log("Building network sponsored by vConTACT3...")
#     GlobalTimer.log("vConTACT3 is now searching for its database...")

#     vcontact3_path = Path(shutil.which('vcontact3')) #Acá es para buscar vcontact3
#     db_path = vcontact3_path.parent.parent / 'db' / 'vcontact3' #Sube 2 niveles, porque aquí están las DB normalmente
#     db_path.mkdir(parents=True, exist_ok=True)

#     existing_versions = list(db_path.glob('*.json'))
#     if existing_versions:
#         GlobalTimer.log(f"vConTACT3 database found, skipping download.")
#     else:
#         GlobalTimer.log("vConTACT3 database not found, downloading latest version...")
#         subprocess.run([
#             'vcontact3', 'prepare_databases',
#             '--get-version', 'latest',
#             '--set-location', str(db_path)
#         ], check=True)

#     subprocess.run([
#         'vcontact3', 'run',
#         '--nucleotide', input_fasta,
#         '--output', output_dir,
#         '--db-path', str(db_path),
#         '--exports', 'graphml'
#     ], check=True)

#     network_files = list(Path(output_dir).rglob('exports/networks/part*.graphml'))
#     graphs = [nx.read_graphml(str(f)) for f in network_files]
#     network = nx.compose_all(graphs)

#     # Normalization before using as an input
#     gene_sharing_weights = nx.get_edge_attributes(network, 'weight')
#     total_weight = sum(gene_sharing_weights.values())
#     normalized_weights = {edge: w / total_weight for edge, w in gene_sharing_weights.items()}
#     nx.set_edge_attributes(network, normalized_weights, 'weight')

#     return network

def get_protein_prediction(header, sequence, faa, gene_finder, genome_proteins):
    
    if not header or not sequence:
        return 0

    genome_id = header.split()[0]
    full_seq = ''.join(sequence)
    genes = gene_finder.find_genes(full_seq.encode())
    genome_proteins[genome_id] = set()
    protein_count = 0

    for i, gene in enumerate(genes):
        protein_id = f"{genome_id}_{i+1}"
        genome_proteins[genome_id].add(protein_id)
        faa.write(f">{protein_id}\n{gene.translate()}\n")
        protein_count += 1
    
    return protein_count

def get_gene_sharing_network(input_fasta, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    protein_fasta = output_dir / 'proteins.faa'
    cluster_out = output_dir / 'diamond_clusters.tsv'

    # Step 1: Predict proteins with Pyrodigal (meta mode for metagenomes)
    GlobalTimer.log("Predicting proteins with Pyrodigal...")
    gene_finder = pyrodigal.GeneFinder(meta=True)
    genome_proteins = {}  # genome_id -> set of protein_ids
    total_proteins = 0

    with open(protein_fasta, 'w') as faa:
        current_header = None
        current_seq = []
        
        with open(input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    total_proteins += get_protein_prediction(current_header, current_seq, faa, gene_finder, genome_proteins)
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)

            total_proteins += get_protein_prediction(current_header, current_seq, faa, gene_finder, genome_proteins)

    GlobalTimer.log(f"Predicted {total_proteins} proteins from {len(genome_proteins)} genomes.")

    # Step 2: Cluster proteins with Diamond linclust (linear scaling, fast, low RAM)
    GlobalTimer.log("Clustering proteins with Diamond linear scaling...")
    subprocess.run([
        'diamond', 'linclust',
        '-d', str(protein_fasta),
        '-o', str(cluster_out),
        '--approx-id', '30', #30% AA identity
        '--member-cover', '80', #80% of coverage
        '-M', '8G',
        '--header'
    ], check=True)

    # Step 3: Parse clusters into a DataFrame
    GlobalTimer.log("Building gene sharing network from protein clusters...")

    clusters = pd.read_csv(cluster_out, sep='\t', comment='#', header=None, names=['representative', 'member'])
    clusters['genome'] = clusters['member'].str.rsplit('_', n=1).str[0]

    # Build genome-pair edge weights using merge (self-join on cluster)
    genome_clusters = clusters[['representative', 'genome']].drop_duplicates()
    pairs = genome_clusters.merge(genome_clusters, on='representative', suffixes=('_a', '_b'))
    pairs = pairs[pairs['genome_a'] < pairs['genome_b']]
    edge_weights = pairs.groupby(['genome_a', 'genome_b']).size().reset_index(name='weight')

    # Build the graph
    network = nx.Graph()
    network.add_nodes_from(clusters['genome'].unique())
    for _, row in edge_weights.iterrows():
        network.add_edge(row['genome_a'], row['genome_b'], weight=float(row['weight']))

    del clusters
    del genome_clusters
    del pairs
    del edge_weights

    # Report network structure
    n_components = nx.number_connected_components(network)
    isolated = sum(1 for n in network.nodes() if network.degree(n) == 0)
    GlobalTimer.log(
        f"Final network: {network.number_of_nodes()} nodes, "
        f"{network.number_of_edges()} edges, "
        f"{n_components} connected components, {isolated} isolated nodes."
    )

    # Normalize weights
    weights = nx.get_edge_attributes(network, 'weight')
    total_weight = sum(weights.values())
    if total_weight > 0:
        normalized = {edge: w / total_weight for edge, w in weights.items()}
        nx.set_edge_attributes(network, normalized, 'weight')

    return network

def get_network(matrix):
    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.values
    network = nx.from_numpy_array(matrix)
    # Meter grafo en versión plot

    return network

def get_network_with_threshold(threshold, matrix, idx_to_seq):
    network = get_network(matrix, idx_to_seq)
    # remove_edges_from modifies in place, returns None
    network.remove_edges_from([
        (n1, n2) for n1, n2, weight in network.edges(data="weight") 
        if weight < threshold
    ])
    return network

## CHECK
# Build network directly from vclust edge list.
# Nodes will be sequence IDs automatically.
def get_network_from_edges(edges_df, threshold):
    query_col = 'query'
    target_col = 'reference'
    ani_col = 'tani'
    
    network = nx.Graph()
    
    # Add edges with ANI >= threshold
    for _, row in edges_df.iterrows():
        ani_value = float(row[ani_col])
        if ani_value >= threshold:
            network.add_edge(row[query_col], row[target_col], weight=ani_value)
    
    # Add self-loops with weight 1.0 (optional)
    for node in network.nodes():
        network.nodes[node]['self_ani'] = 1.0
    
    return network

def set_community_atribute_on_nodes(network, count_table):
    # Ensure OTU IDs are the index
    if count_table.index.dtype == 'int64':
        count_table = count_table.set_index(count_table.columns[0])
    
    # Iterate over NETWORK nodes, not count_table
    for node in network.nodes():
        if node in count_table.index:
            row = count_table.loc[node]
            communities = row[row > 0].to_dict()
            network.nodes[node]['communities'] = communities
        else:
            # Node exists in network but not in count_table
            network.nodes[node]['communities'] = {}
    
    return network

def get_input_ani_network(input_fasta, output_dir, threshold, count_table):
    get_ani_matrix_vclust(input_fasta, output_dir)

    edges_df = pd.read_csv(f'{output_dir}/ani.tsv', sep='\t')
    
    network = get_network_from_edges(edges_df, threshold)
    
    # Add community attributes
    input_network = set_community_atribute_on_nodes(network, count_table)

    return input_network

def get_input_gene_sharing_network(input_fasta, output_dir, count_table):

    network = get_gene_sharing_network(input_fasta, output_dir)
    print("==AQUI ENTRA EL NETWORK==")
    print(network)
    input_network = set_community_atribute_on_nodes(network, count_table)
    
    return input_network
