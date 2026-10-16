import os
import itertools
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pyrodigal
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from src.utils import GlobalTimer

# Network used by the frac calculations: one entry per edge in integer arrays, plus node names.
# edges_by_source (CSR nodes x edges) returns the edges of a set of nodes without scanning them all.
# abundance (CSC nodes x samples) replaces the 'communities' attribute of the nx.Graph.
@dataclass
class SparseNetwork:
    node_names: np.ndarray
    edge_source: np.ndarray
    edge_target: np.ndarray
    edge_weight: np.ndarray
    edges_by_source: sparse.csr_matrix
    abundance: sparse.csc_matrix = None
    sample_names: list = None

# Reads the count table and indexes it by its first column (the sequence IDs).
def get_count_table(count_table_path):
    table = pd.read_csv(count_table_path)
    indexed_count_table = table.set_index(table.columns[0])

    return indexed_count_table

# Runs vclust prefilter + align and returns the path to the ANI edge list.
# Without filtering: ani.tsv is shared by the network path, the tree path and external scripts.
def get_vclust_run(input_fasta, output_dir):
    GlobalTimer.log("Prefiltering sponsored by vclust...")
    cmd_prefilter = [
        'vclust', 'prefilter',
        '-i', input_fasta,
        '-o', f'{output_dir}/filter.txt'
    ]
    subprocess.run(cmd_prefilter, check=True)

    GlobalTimer.log("Alignment sponsored by vclust...")
    cmd_align = [
        'vclust', 'align',
        '-i', input_fasta,
        '-o', f'{output_dir}/ani.tsv',
        '--filter', f'{output_dir}/filter.txt'
    ]
    subprocess.run(cmd_align, check=True)

    return f'{output_dir}/ani.tsv'

# Runs vclust and converts the edge list into the dense matrix (tree path).
def get_ani_matrix_vclust(input_fasta, output_dir):
    ani_path = get_vclust_run(input_fasta, output_dir)
    edges_df = pd.read_csv(ani_path, sep='\t')
    matrix = set_edge_list_to_matrix(edges_df)

    return matrix

# Returns the row positions to keep so each pair (A-B or B-A) appears 
# once, keeping its last occurrence. A repeated edge collapses into one 
# and keeps the last weight.
def get_last_unique_edges(edge_source, edge_target, n_nodes):
    pair_low = np.minimum(edge_source, edge_target).astype(np.int64)
    pair_high = np.maximum(edge_source, edge_target)
    pair_keys = pair_low * n_nodes + pair_high

    _, first_idx_in_reversed = np.unique(pair_keys[::-1], return_index=True)

    return len(pair_keys) - 1 - first_idx_in_reversed

# Converts vclust edge list output to a square adjacency matrix.
# For this, the function is using tANI as the ANI value. This is
# controlled by the ani_col value. Vectorized: names -> indices
# in one path. If a pair appears more than once (A-B or B-A),
# the last row wins.
def set_edge_list_to_matrix(edges_df):
    if 'query' in edges_df.columns:
        query_col = 'query'
    else:
        raise ValueError(f"Could not find query column. Available: {edges_df.columns.tolist()}")

    if 'reference' in edges_df.columns:
        target_col = 'reference'
    elif 'target' in edges_df.columns:
        target_col = 'target'
    else:
        raise ValueError(f"Could not find target/reference column. Available: {edges_df.columns.tolist()}")

    if 'tani' in edges_df.columns:
        ani_col = 'tani'
    else:
        raise ValueError(f"Could not find ANI column. Available: {edges_df.columns.tolist()}")

    all_seqs = pd.unique(edges_df[[query_col, target_col]].values.ravel())
    n = len(all_seqs)
    seq_index = pd.Index(all_seqs)

    query_idx = seq_index.get_indexer(edges_df[query_col])
    target_idx = seq_index.get_indexer(edges_df[target_col])
    ani_values = edges_df[ani_col].to_numpy(dtype=float)

    kept_rows = get_last_unique_edges(query_idx, target_idx, n)
    unique_low = np.minimum(query_idx, target_idx)[kept_rows]
    unique_high = np.maximum(query_idx, target_idx)[kept_rows]
    unique_ani = ani_values[kept_rows]

    matrix = np.zeros((n, n))
    matrix[unique_low, unique_high] = unique_ani
    matrix[unique_high, unique_low] = unique_ani

    np.fill_diagonal(matrix, 1.0)

    return matrix

# Builds the SparseNetwork from the edge arrays 
# (node indices and weights). edges_by_source stores each edge in the
# row of its source node, so each edge appears only once.
def get_sparse_network(node_names, edge_source, edge_target, edge_weight):
    n_edges = len(edge_source)
    edges_by_source = sparse.csr_matrix(
        (np.ones(n_edges, dtype=np.int8), (edge_source, np.arange(n_edges))),
        shape=(len(node_names), n_edges)
    )

    return SparseNetwork(
        node_names=np.asarray(node_names),
        edge_source=np.asarray(edge_source, dtype=np.int32),
        edge_target=np.asarray(edge_target, dtype=np.int32),
        edge_weight=np.asarray(edge_weight, dtype=float),
        edges_by_source=edges_by_source
    )

# def get_ani_matrix_vclust(input_fasta, output_dir):
#     # Prefilter
#     GlobalTimer.log("Prefiltering sponsored by vclust...")
#     cmd_prefilter = [
#         'vclust', 'prefilter',
#         '-i', input_fasta,
#         '-o', f'{output_dir}/filter.txt'
#     ]
#     subprocess.run(cmd_prefilter, check=True)

#     # Alignment
#     GlobalTimer.log("Alignment sponsored by vclust...")
#     cmd_align = [
#         'vclust', 'align',
#         '-i', input_fasta,
#         '-o', f'{output_dir}/ani.tsv',
#         '--filter', f'{output_dir}/filter.txt'
#     ]
#     subprocess.run(cmd_align, check=True)

#     edges_df = pd.read_csv(f'{output_dir}/ani.tsv', sep='\t')
#     matrix = set_edge_list_to_matrix(edges_df)

#     return matrix

# get_protein_prediction takes a FASTA file and predicts the genes,
# transcribing how many proteins it found and then returning the
# number of proeins it found. 
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

# Reads a FASTA file and yields (genome_id, sequence) one by one 
# without loading the whole file. Skips records with an empty 
# header or sequence, like the original.
def get_read_fasta_records(input_fasta):
    current_header = None
    current_seq = []

    with open(input_fasta, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    yield current_header.split()[0], ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

    if current_header and current_seq:
        yield current_header.split()[0], ''.join(current_seq)

# Groups an iterator of records into lists of batch_size, to bound memory.
def get_record_batches(records, batch_size):
    while True:
        batch = list(itertools.islice(records, batch_size))
        if not batch:
            return
        yield batch

# GeneFinder for each worker process; each process creates its own so they don't compete with each other.
process_gene_finder = None

# Runs once when each process starts: creates that process's GeneFinder.
def init_gene_finder():
    global process_gene_finder
    process_gene_finder = pyrodigal.GeneFinder(meta=True)

# Predicts the genes of one genome inside a worker process and returns its proteins
# already formatted as FASTA, plus how many there are.
def predict_genome_proteins(genome_id, sequence):
    genes = process_gene_finder.find_genes(sequence.encode())
    fasta_text = ''.join(f">{genome_id}_{i + 1}\n{gene.translate()}\n" for i, gene in enumerate(genes))
    return fasta_text, len(genes)

# Predicts proteins (Pyrodigal in parallel), clusters them with Diamond, and builds the gene-sharing
# network: weight = number of shared clusters / maximum in the network.
# Shared clusters are computed as B @ B.T (sparse genome x cluster matrix) instead of the merge.
# Returns a SparseNetwork.
def get_gene_sharing_network(input_fasta, output_dir, threads=None, batch_size=10000):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    protein_fasta = output_dir / 'proteins.faa'
    cluster_out = output_dir / 'diamond_clusters.tsv'

    # Step 1: Predict proteins with Pyrodigal in parallel processes (meta mode for metagenomes)
    GlobalTimer.log("Predicting proteins with Pyrodigal...")
    seen_genomes = set()
    total_proteins = 0
    n_workers = threads or os.cpu_count()

    fasta_records = get_read_fasta_records(input_fasta)
    with open(protein_fasta, 'w') as faa, ProcessPoolExecutor(max_workers=n_workers, initializer=init_gene_finder) as executor:
        for record_batch in get_record_batches(fasta_records, batch_size):
            genome_ids = [genome_id for genome_id, _ in record_batch]
            sequences = [sequence for _, sequence in record_batch]
            predictions = executor.map(predict_genome_proteins, genome_ids, sequences, chunksize=20)

            fasta_texts, protein_counts = zip(*predictions)
            faa.write(''.join(fasta_texts))
            total_proteins += sum(protein_counts)
            seen_genomes.update(genome_ids)

    GlobalTimer.log(f"Predicted {total_proteins} proteins from {len(seen_genomes)} genomes.")
    del seen_genomes

    # Step 2: Cluster proteins with Diamond linclust (linear scaling, fast, low RAM)
    GlobalTimer.log("Clustering proteins with Diamond linear scaling...")
    subprocess.run([
        'diamond', 'linclust',
        '-d', str(protein_fasta),
        '-o', str(cluster_out),
        '--approx-id', '30',  # 30% AA identity
        '--member-cover', '80',  # 80% of coverage
        '-M', '8G',
        '--header'
    ], check=True)

    # Step 3: Build the gene-sharing network from the clusters
    GlobalTimer.log("Building gene sharing network from protein clusters...")

    clusters = pd.read_csv(cluster_out, sep='\t', comment='#', header=None, names=['representative', 'member'])
    clusters['genome'] = clusters['member'].str.rsplit('_', n=1).str[0]

    genome_codes, genome_names = pd.factorize(clusters['genome'], sort=True)
    cluster_codes, cluster_names = pd.factorize(clusters['representative'])
    del clusters

    n_genomes = len(genome_names)
    genome_by_cluster = sparse.csr_matrix(
        (np.ones(len(genome_codes), dtype=np.int32), (genome_codes, cluster_codes)),
        shape=(n_genomes, len(cluster_names))
    )
    genome_by_cluster.data[:] = 1  # paralogs get summed; collapse to presence/absence
    del genome_codes, cluster_codes, cluster_names

    # shared_clusters[a, b] = number of clusters genomes a and b share
    shared_clusters = genome_by_cluster @ genome_by_cluster.T
    genome_pairs = sparse.triu(shared_clusters, k=1).tocsr()
    genome_pairs.sort_indices()
    del genome_by_cluster, shared_clusters

    pair_rows = np.repeat(np.arange(n_genomes), np.diff(genome_pairs.indptr))
    pair_cols = genome_pairs.indices
    shared_counts = genome_pairs.data.astype(float)

    # Report network structure
    n_components, _ = connected_components(genome_pairs, directed=False)
    node_degrees = np.bincount(pair_rows, minlength=n_genomes) + np.bincount(pair_cols, minlength=n_genomes)
    isolated = int(np.sum(node_degrees == 0))
    GlobalTimer.log(
        f"Final network: {n_genomes} nodes, "
        f"{len(shared_counts)} edges, "
        f"{n_components} connected components, {isolated} isolated nodes."
    )

    # Normalize edge weights: Edge / Max edge in network
    # Normalize edge weights

    # Total weight
    # weights = nx.get_edge_attributes(network, 'weight')
    # total_weight = sum(weights.values())
    # if total_weight > 0:
    #     normalized = {edge: w / total_weight for edge, w in weights.items()}
    #     nx.set_edge_attributes(network, normalized, 'weight')
    edge_weights = shared_counts
    if len(shared_counts) > 0 and shared_counts.max() > 0:
        edge_weights = shared_counts / shared_counts.max()

    return get_sparse_network(genome_names, pair_rows, pair_cols, edge_weights)

# def get_gene_sharing_network(input_fasta, output_dir):
#     output_dir = Path(output_dir)
#     output_dir.mkdir(parents=True, exist_ok=True)

#     protein_fasta = output_dir / 'proteins.faa'
#     cluster_out = output_dir / 'diamond_clusters.tsv'

#     # Step 1: Predict proteins with Pyrodigal (meta mode for metagenomes)
#     GlobalTimer.log("Predicting proteins with Pyrodigal...")
#     gene_finder = pyrodigal.GeneFinder(meta=True)
#     genome_proteins = {}  # genome_id -> set of protein_ids
#     total_proteins = 0

#     with open(protein_fasta, 'w') as faa:
#         current_header = None
#         current_seq = []
        
#         with open(input_fasta, 'r') as f:
#             for line in f:
#                 line = line.strip()
#                 if line.startswith('>'):
#                     total_proteins += get_protein_prediction(current_header, current_seq, faa, gene_finder, genome_proteins)
#                     current_header = line[1:]
#                     current_seq = []
#                 else:
#                     current_seq.append(line)

#             total_proteins += get_protein_prediction(current_header, current_seq, faa, gene_finder, genome_proteins)

#     GlobalTimer.log(f"Predicted {total_proteins} proteins from {len(genome_proteins)} genomes.")

#     # Step 2: Cluster proteins with Diamond linclust (linear scaling, fast, low RAM)
#     GlobalTimer.log("Clustering proteins with Diamond linear scaling...")
#     subprocess.run([
#         'diamond', 'linclust',
#         '-d', str(protein_fasta),
#         '-o', str(cluster_out),
#         '--approx-id', '30', #30% AA identity
#         '--member-cover', '80', #80% of coverage
#         '-M', '8G',
#         '--header'
#     ], check=True)

#     # Step 3: Parse clusters into a DataFrame
#     GlobalTimer.log("Building gene sharing network from protein clusters...")

#     clusters = pd.read_csv(cluster_out, sep='\t', comment='#', header=None, names=['representative', 'member'])
#     clusters['genome'] = clusters['member'].str.rsplit('_', n=1).str[0]

#     # Build genome-pair edge weights using merge (self-join on cluster)
#     genome_clusters = clusters[['representative', 'genome']].drop_duplicates()
#     pairs = genome_clusters.merge(genome_clusters, on='representative', suffixes=('_a', '_b'))
#     pairs = pairs[pairs['genome_a'] < pairs['genome_b']]
#     edge_weights = pairs.groupby(['genome_a', 'genome_b']).size().reset_index(name='weight')

#     # Build the graph
#     network = nx.Graph()
#     network.add_nodes_from(clusters['genome'].unique())
#     for _, row in edge_weights.iterrows():
#         network.add_edge(row['genome_a'], row['genome_b'], weight=float(row['weight']))

#     del clusters
#     del genome_clusters
#     del pairs
#     del edge_weights

#     # Report network structure
#     n_components = nx.number_connected_components(network)
#     isolated = sum(1 for n in network.nodes() if network.degree(n) == 0)
#     GlobalTimer.log(
#         f"Final network: {network.number_of_nodes()} nodes, "
#         f"{network.number_of_edges()} edges, "
#         f"{n_components} connected components, {isolated} isolated nodes."
#     )

#     # Normalize edge weights

#     # Total weight
#     # weights = nx.get_edge_attributes(network, 'weight')
#     # total_weight = sum(weights.values())
#     # if total_weight > 0:
#     #     normalized = {edge: w / total_weight for edge, w in weights.items()}
#     #     nx.set_edge_attributes(network, normalized, 'weight')

#     # Edge / Max edge in network
#     weights = nx.get_edge_attributes(network, 'weight')
#     max_weight = max(weights.values())
#     if max_weight > 0:
#         normalized = {edge: w / max_weight for edge, w in weights.items()}
#         nx.set_edge_attributes(network, normalized, 'weight')

#     return network

# Converts a dense matrix into an nx.Graph (unused code, kept for compatibility).
# def get_network(matrix):
#     if isinstance(matrix, pd.DataFrame):
#         matrix = matrix.values
#     network = nx.from_numpy_array(matrix)
#     # Meter grafo en versión plot

#     return network

# Builds the nx.Graph from the matrix keeping only weights >= threshold (unused code).
# idx_to_seq is kept in the signature but not used, same as before.
# def get_network_with_threshold(threshold, matrix, idx_to_seq):
#     if isinstance(matrix, pd.DataFrame):
#         matrix = matrix.values

#     thresholded_matrix = np.where(matrix >= threshold, matrix, 0.0)
#     network = get_network(thresholded_matrix)

#     return network

# Builds the ANI network from the vclust edge list, keeping edges with tANI >= threshold.
# Repeated edges (A-B or B-A) collapse into one and keep the last weight.
def get_network_from_edges(edges_df, threshold):
    query_col = 'query'
    target_col = 'reference'
    ani_col = 'tani'

    ani_values = edges_df[ani_col].to_numpy(dtype=float)
    passes_threshold = ani_values >= threshold
    passing_edges = edges_df.loc[passes_threshold]
    n_passing = len(passing_edges)

    all_names = pd.concat([passing_edges[query_col], passing_edges[target_col]], ignore_index=True)
    node_codes, node_names = pd.factorize(all_names)
    edge_source = node_codes[:n_passing]
    edge_target = node_codes[n_passing:]

    kept_edges = get_last_unique_edges(edge_source, edge_target, len(node_names))

    return get_sparse_network(
        node_names,
        edge_source[kept_edges],
        edge_target[kept_edges],
        ani_values[passes_threshold][kept_edges]
    )

# Aligns the count table to the network's node order and stores it as a sparse nodes x samples matrix.
# Only counts > 0 are stored (equivalent to the old 'communities' dicts).
# Nodes missing from the table are left without abundance. Raises an error if there are duplicated IDs.
def set_community_atribute_on_nodes(network, count_table):
    if count_table.index.dtype == 'int64':
        count_table = count_table.set_index(count_table.columns[0])

    if not count_table.index.is_unique:
        duplicated_ids = count_table.index[count_table.index.duplicated()].unique().tolist()
        raise ValueError(f"Count table has {len(duplicated_ids)} duplicated IDs, e.g.: {duplicated_ids[:10]}")

    table_rows = count_table.index.get_indexer(network.node_names)
    nodes_in_table = np.flatnonzero(table_rows >= 0)
    counts = count_table.iloc[table_rows[nodes_in_table]].to_numpy(dtype=float)

    count_rows, sample_columns = np.nonzero(counts > 0)
    network.abundance = sparse.csc_matrix(
        (counts[count_rows, sample_columns], (nodes_in_table[count_rows], sample_columns)),
        shape=(len(network.node_names), count_table.shape[1])
    )
    network.sample_names = count_table.columns.tolist()

    return network

# Runs vclust once and builds the ANI network with the abundances aligned.
# No longer builds the dense matrix or reads ani.tsv twice.
def get_input_ani_network(input_fasta, output_dir, threshold, count_table):
    ani_path = get_vclust_run(input_fasta, output_dir)

    edges_df = pd.read_csv(ani_path, sep='\t', usecols=['query', 'reference', 'tani'])
    network = get_network_from_edges(edges_df, threshold)
    del edges_df

    input_network = set_community_atribute_on_nodes(network, count_table)
    GlobalTimer.log("Input: ANI network")

    return input_network

# Builds the gene-sharing network with the abundances aligned.
def get_input_gene_sharing_network(input_fasta, output_dir, count_table):
    network = get_gene_sharing_network(input_fasta, output_dir)
    GlobalTimer.log("Input: Gene-sharing network")
    input_network = set_community_atribute_on_nodes(network, count_table)

    return input_network
