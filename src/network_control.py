import subprocess
import networkx as nx
import numpy as np
import pandas as pd
import glob
from src.utils import GlobalTimer
import shutil
from pathlib import Path

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
def get_gene_sharing_network_vcontact3(input_fasta, output_dir):
    GlobalTimer.log("Building network sponsored by vConTACT3...")
    GlobalTimer.log("vConTACT3 is now searching for its database...")

    vcontact3_path = Path(shutil.which('vcontact3')) #Acá es para buscar vcontact3
    db_path = vcontact3_path.parent.parent / 'db' / 'vcontact3' #Sube 2 niveles, porque aquí están las DB normalmente
    db_path.mkdir(parents=True, exist_ok=True)

    existing_versions = list(db_path.glob('*.json'))
    if existing_versions:
        GlobalTimer.log(f"vConTACT3 database found, skipping download.")
    else:
        GlobalTimer.log("vConTACT3 database not found, downloading latest version...")
        subprocess.run([
            'vcontact3', 'prepare_databases',
            '--get-version', 'latest',
            '--set-location', str(db_path)
        ], check=True)

    subprocess.run([
        'vcontact3', 'run',
        '--nucleotide', input_fasta,
        '--output', output_dir,
        '--db-path', str(db_path),
        '--exports', 'graphml'
    ], check=True)
    
    #network_file = glob.glob(f'{output_dir}/*.graphml')[0]
    #network = nx.read_graphml(network_file)

    network_files = list(Path(output_dir).rglob('exports/networks/part*.graphml'))
    graphs = [nx.read_graphml(str(f)) for f in network_files]
    network = nx.compose_all(graphs)

    # Normalization before using as an input
    gene_sharing_weights = nx.get_edge_attributes(network, 'weight')
    total_weight = sum(gene_sharing_weights.values())
    normalized_weights = {edge: w / total_weight for edge, w in gene_sharing_weights.items()}
    nx.set_edge_attributes(network, normalized_weights, 'weight')

    return network

def get_gene_sharing_network_vcontact2(input_fasta, output_dir):
    protein_fasta = f'{output_dir}/proteins.faa'
    diamond_out = f'{output_dir}/diamond.tsv'

    # Predecir proteínas
    GlobalTimer.log("Predicting proteins with Prodigal...")
    subprocess.run([
        'conda', 'run', '-n', 'vcontact2_env',
        'prodigal', '-p', 'meta',
        '-i', input_fasta,
        '-a', protein_fasta,
        '-f', 'gff'
    ], check=True)

    # All-vs-all con Diamond
    GlobalTimer.log("Running Diamond all-vs-all...")
    subprocess.run([
        'conda', 'run', '-n', 'vcontact2_env',
        'diamond', 'blastp',
        '-q', protein_fasta,
        '-d', protein_fasta,
        '-o', diamond_out,
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident',
        '--sensitive'
    ], check=True)

    # Construir network contando genes compartidos por par de genomas
    GlobalTimer.log("Building gene sharing network...")
    hits = pd.read_csv(diamond_out, sep='\t', names=['query', 'subject', 'pident'])
    
    # Extraer nombre del genoma de cada proteína (asume formato prodigal: genoma_1, genoma_2...)
    hits['genome1'] = hits['query'].str.rsplit('_', n=1).str[0]
    hits['genome2'] = hits['subject'].str.rsplit('_', n=1).str[0]
    
    # Filtrar self-hits y contar genes compartidos por par
    hits = hits[hits['genome1'] != hits['genome2']]
    edges = hits.groupby(['genome1', 'genome2']).size().reset_index(name='weight')
    
    G = nx.from_pandas_edgelist(edges, source='genome1', target='genome2', edge_attr='weight')
    
    # Normalizar pesos
    weights = nx.get_edge_attributes(G, 'weight')
    total_weight = sum(weights.values())
    normalized_weights = {edge: w / total_weight for edge, w in weights.items()}
    nx.set_edge_attributes(G, normalized_weights, 'weight')

    return G

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

def get_input_vcontact3_network(input_fasta, output_dir, count_table):

    network = get_gene_sharing_network_vcontact2(input_fasta, output_dir)
    print("==AQUI ENTRA EL NETWORK==")
    print(network)
    input_network = set_community_atribute_on_nodes(network, count_table)
    
    return input_network
