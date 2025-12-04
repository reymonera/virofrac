import subprocess
import networkx as nx
import numpy as np
import pandas as pd

# Esta función debería de acepar el network basado en ani y en gene sharing
def get_count_table(count_table_path):
    ## Creo que no necestamos esto por el tipo de input
    table = pd.read_csv(count_table_path)
    indexed_count_table = table.set_index([0])

    return indexed_count_table

##CHECK
def set_edge_list_to_matrix(edges_df):
    """
    Convert vclust edge list output to a square adjacency matrix.
    """
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

    print(f"DEBUG: Using columns - query: {query_col}, target: {target_col}, ani: {ani_col}")

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
    
    print(f"DEBUG: Adjacency matrix shape = {matrix.shape}")
    
    return matrix

def get_ani_matrix_vclust(input_fasta, output_dir):
    # Prefilter
    print("Prefiltering sponsored by vclust...")
    cmd_prefilter = [
        'vclust', 'prefilter',
        '-i', input_fasta,
        '-o', f'{output_dir}/filter.txt'
    ]
    subprocess.run(cmd_prefilter, check=True)

    # Alignment
    print("Alignment sponsored by vclust...")
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

def get_gene_sharing_matrix_vcontact():
    matrix = 0
    return matrix

def get_network(matrix):
    #matrix = np.loadtxt(matrix_file, delimiter=",", skiprows=1)
    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.values
    network = nx.from_numpy_array(matrix)
    # Meter grafo en versión plot

    return network

def get_network_with_threshold(threshold, matrix):
    network = get_network(matrix)
    threshold_network = network.remove_edges_from([(n1, n2) for n1, n2, 
                               weight in network.edges(data="weight") 
                               if weight < threshold])
    
    return threshold_network

def set_community_atribute_on_nodes(network, count_table):
    #count_table = get_count_table(count_table_path)
    
    for node in count_table.index:
        row = count_table.loc[node]
        communities = row[row > 0].to_dict() # Diccionario de comunidades y cantidades
        network.nodes[node]['communities'] = communities 
        # Ahora G.nodes[0]['comunidades'] = {'comunidad_A': 5, 'comunidad_B': 2}
    
    return network

def get_input_ani_network(input_fasta, output_dir, threshold, count_table):
    matrix = get_ani_matrix_vclust(input_fasta, output_dir)
    pruned_network = get_network_with_threshold(threshold, matrix)
    input_network = set_community_atribute_on_nodes(pruned_network, count_table)

    return input_network


