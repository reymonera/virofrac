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
    import xml.etree.ElementTree as ET

    GlobalTimer.log("Building network sponsored by vConTACT3...")
    GlobalTimer.log("vConTACT3 is now searching for its database...")

    vcontact3_path = Path(shutil.which('vcontact3'))
    db_path = vcontact3_path.parent.parent / 'db' / 'vcontact3'
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

    # Identify query genomes from the input fasta
    query_genomes = set()
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                query_genomes.add(line.strip().lstrip('>').split()[0])

    network_files = list(Path(output_dir).glob('exports/networks/part*.graphml'))
    GlobalTimer.log(f"Found {len(network_files)} partition files to stream.")

    ns = {'g': 'http://graphml.graphdrawing.org/xmlns'}

    def stream_edges(graphml_path):
        """Yield (source, target, weight) tuples from a GraphML file without loading it fully."""
        weight_key_id = None
        for event, elem in ET.iterparse(str(graphml_path), events=('end',)):
            tag = elem.tag.split('}')[-1]
            if tag == 'key' and elem.get('attr.name') == 'weight':
                weight_key_id = elem.get('id')
                elem.clear()
                break
            if tag == 'key':
                elem.clear()
            if tag == 'graph':
                break

        context = ET.iterparse(str(graphml_path), events=('end',))
        for event, elem in context:
            tag = elem.tag.split('}')[-1]
            if tag == 'edge':
                src = elem.get('source')
                tgt = elem.get('target')
                w = 1.0
                for data in elem.findall('g:data', ns):
                    if weight_key_id is None or data.get('key') == weight_key_id:
                        try:
                            w = float(data.text)
                        except (TypeError, ValueError):
                            pass
                        if weight_key_id is not None:
                            break
                yield src, tgt, w
                elem.clear()
            elif tag in ('node', 'key', 'data'):
                elem.clear()

    # Accumulate only direct query-to-query edges
    query_to_query_edges = {}

    for nf in network_files:
        GlobalTimer.log(f"Streaming {nf.name}...")
        edge_count = 0
        kept = 0
        for src, tgt, w in stream_edges(nf):
            edge_count += 1
            if src in query_genomes and tgt in query_genomes:
                key = (src, tgt) if src < tgt else (tgt, src)
                if key not in query_to_query_edges or w > query_to_query_edges[key]:
                    query_to_query_edges[key] = w
                kept += 1
        GlobalTimer.log(f"  {edge_count} edges scanned, {kept} retained.")

    # Build the final graph
    query_network = nx.Graph()
    query_network.add_nodes_from(query_genomes)
    for (a, b), w in query_to_query_edges.items():
        query_network.add_edge(a, b, weight=w)

    del query_to_query_edges

    # Report network structure
    n_components = nx.number_connected_components(query_network)
    isolated = sum(1 for n in query_network.nodes() if query_network.degree(n) == 0)
    GlobalTimer.log(
        f"Final network: {query_network.number_of_nodes()} nodes, "
        f"{query_network.number_of_edges()} edges, "
        f"{n_components} connected components, {isolated} isolated nodes."
    )

    # Normalize weights
    gene_sharing_weights = nx.get_edge_attributes(query_network, 'weight')
    total_weight = sum(gene_sharing_weights.values())
    if total_weight > 0:
        normalized_weights = {edge: w / total_weight for edge, w in gene_sharing_weights.items()}
        nx.set_edge_attributes(query_network, normalized_weights, 'weight')

    return query_network

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
