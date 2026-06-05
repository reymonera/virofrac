import numpy as np
import pandas as pd
from ete3 import Tree
import src.network_control as net_control
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
#Provisional
import time
from src.utils import GlobalTimer

# This function gets the tree stored in the database, and
# outputs the original .newick
def get_taxa_text_tree():
    file = 'tree_db/tree_ictv.newick'
    newick_tree = open(file, 'r').read().replace('()', '')

    return newick_tree

# This function should be the one putting the taxa in a
# browsable list for the tree prune function to work.
def get_otus_taxa_list(otu_tax_table):
    all_taxa = set()
    
    for col in otu_tax_table.columns[1:]:
        taxa = otu_tax_table[col].dropna().unique()
        all_taxa.update(taxa)
    
    tax_otus = list(all_taxa)

    return tax_otus

# This function returns a dictionary based on the OTUs
# taxonomic ranks. This is done through the 
def get_otu_dictionary(otu_tax_table):
    id_column = otu_tax_table.columns[0]
    otu_dict = otu_tax_table.set_index(id_column).to_dict('index')

    return otu_dict

# This function outputs the pruned tree. For this, it
# will first extract the names in the taxonomy table,
# then it will search for matches in the tree (node
# matches), it will save the ancestors and output a
# pruned tree with its ancestors.
def get_pruned_tree_1(otu_tax_table):
    nodes_to_keep = get_otus_taxa_list(otu_tax_table)
    tree = Tree(get_taxa_text_tree(), format=1)

    node_index = {node.name: node for node in tree.traverse() if node.name}

    ancestors = []
    ancestors_and_nodes = list(nodes_to_keep)

    for i, node_name in enumerate(nodes_to_keep):
        node = node_index.get(node_name)
        
        if node:
            ancestors.extend(node.get_ancestors())
    
    for ancestor in ancestors:
        if ancestor.name:
            ancestors_and_nodes.append(ancestor.name)

    unique_nodes = list(set(ancestors_and_nodes))    
    tree.prune(unique_nodes, preserve_branch_length=True)
    tree.write(outfile="pruned_tree.newick", format=1)
    
    return tree

# This function gets the count of nan and list of 
# non-nan values in a dictionary. This is mainly 
# for the get_otu_tree function.
def get_nan_non_nan_values(taxonomy_in_dictionary):
    all_values = list(taxonomy_in_dictionary.values())
    non_nan_values = []
    
    last_valid_idx = -1
    for i, v in enumerate(all_values):
        if pd.notna(v) and str(v).strip():
            non_nan_values.append(v)
            last_valid_idx = i
    
    if last_valid_idx == -1:
        return [], 0
    
    nan_count = sum(1 for v in all_values[last_valid_idx+1:] 
                    if pd.isna(v) or not str(v).strip())
    
    return non_nan_values, nan_count

# This function has as an output the final OTU tree, which
# is then stored in a newick. For OTUs, this function will
# calculate the maximum depth of. a tree, so that any NA
# OTUs will be placed correctly with the optimum distance
# on the tree. Every OTU should be on the species level.
def get_otu_tree_1(otu_tax_table):
    tree = get_pruned_tree(otu_tax_table)
    otu_dictionary = get_otu_dictionary(otu_tax_table)

    node_index = {node.name: node for node in tree.traverse() if node.name}
    max_depth = max(tree.get_distance(leaf) for leaf in tree.iter_leaves())

    for (otu_id, taxonomy) in otu_dictionary.items():
        
        non_nan_values, nan_count = get_nan_non_nan_values(taxonomy)
        
        if not non_nan_values:
            continue
        
        most_specific = non_nan_values[-1]
        node = node_index.get(most_specific)
        
        if node:
            current_depth = tree.get_distance(node)
            missing_depth = max_depth - current_depth
            
            node.add_child(name=otu_id, dist=missing_depth)
    
    tree.write(outfile="otu_tree.newick", format=1)

    return tree

def get_pruned_tree(otu_tax_table):
    nodes_to_keep = get_otus_taxa_list(otu_tax_table)
    tree = Tree(get_taxa_text_tree(), format=1)

    node_index = {node.name: node for node in tree.traverse() if node.name}
    tree_nodes_set = set(node_index.keys())

    # Filtrar solo nodos que existen en el árbol
    valid_nodes = [n for n in nodes_to_keep if n in tree_nodes_set]
    missing_nodes = [n for n in nodes_to_keep if n not in tree_nodes_set]
    
    if missing_nodes:
        print(f"⚠️  {len(missing_nodes)} nodos no encontrados en árbol ICTV:")
        print(f"   {missing_nodes[:10]}...")

    ancestors = []
    ancestors_and_nodes = list(valid_nodes)

    for node_name in valid_nodes:
        node = node_index.get(node_name)
        if node:
            ancestors.extend(node.get_ancestors())
    
    for ancestor in ancestors:
        if ancestor.name:
            ancestors_and_nodes.append(ancestor.name)

    unique_nodes = list(set(ancestors_and_nodes))
    
    if not unique_nodes:
        raise ValueError("No se encontró ningún nodo válido en el árbol ICTV")
    
    tree.prune(unique_nodes, preserve_branch_length=True)
    tree.write(outfile="pruned_tree.newick", format=1)
    
    return tree


def get_otu_tree(otu_tax_table):
    # Obtener max_depth del árbol ICTV original
    full_tree = Tree(get_taxa_text_tree(), format=1)
    max_depth = max(full_tree.get_distance(leaf) for leaf in full_tree.iter_leaves())
    
    # Ahora sí podar
    tree = get_pruned_tree(otu_tax_table)
    otu_dictionary = get_otu_dictionary(otu_tax_table)

    node_index = {node.name: node for node in tree.traverse() if node.name}

    otus_not_placed = []

    for (otu_id, taxonomy) in otu_dictionary.items():
        
        non_nan_values, nan_count = get_nan_non_nan_values(taxonomy)
        
        if not non_nan_values:
            continue
        
        # Buscar desde el más específico hacia arriba hasta encontrar nodo válido
        node = None
        for tax_name in reversed(non_nan_values):
            node = node_index.get(tax_name)
            if node:
                break
        
        if node:
            current_depth = tree.get_distance(node)
            missing_depth = max_depth - current_depth
            node.add_child(name=otu_id, dist=missing_depth)
        else:
            otus_not_placed.append(otu_id)
    
    # if otus_not_placed:
    #     print(f"WARNING  {len(otus_not_placed)} OTUs no pudieron colocarse en el árbol")
    
    tree.write(outfile="otu_tree.newick", format=1)

    return tree

def get_ani_tree(fasta_file, output_dir):
    net_control.get_ani_matrix_vclust(fasta_file, output_dir)

    ani_results = pd.read_csv(f'{output_dir}/ani.tsv', sep='\t')
    labels = list(ani_results['query'].unique())
    ani_matrix = net_control.set_edge_list_to_matrix(ani_results)

    #labels = list(ani_matrix.index)
    distance = 1.0 - ani_matrix
    linkage_matrix = linkage(squareform(distance), method="average")
    num_assemblies = len(labels)
    
    nodes = {i: Tree(name=labels[i]) for i in range(num_assemblies)}

    for i, (left, right, dist_val, _) in enumerate(linkage_matrix):
        tree = Tree()
        tree.dist = dist_val
        tree.add_child(nodes[int(left)])
        tree.add_child(nodes[int(right)])
        nodes[num_assemblies + i] = tree
    
    return tree