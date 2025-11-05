import pandas as pd
from ete3 import Tree
#Provisional
import time

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
def get_pruned_tree(otu_tax_table):
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

# def get_pruned_tree(otu_tax_table):
#     nodes_to_keep = get_otus_taxa_list(otu_tax_table)
#     tree = Tree(get_taxa_text_tree(), format=1)
    
#     print(f"Creando índice del árbol...")
#     # CREAR ÍNDICE UNA SOLA VEZ - O(n)
#     node_index = {node.name: node for node in tree.traverse() if node.name}
#     print(f"✓ Índice creado: {len(node_index)} nodos")
    
#     ancestors = []
#     ancestors_and_nodes = list(nodes_to_keep)
    
#     print(f"Buscando ancestros para {len(nodes_to_keep)} nodos...")
#     for i, node_name in enumerate(nodes_to_keep):
#         if i % 1000 == 0:
#             print(f"  Procesando {i}/{len(nodes_to_keep)}...")
        
#         # BÚSQUEDA O(1) en lugar de O(n)
#         node = node_index.get(node_name)
        
#         if node:
#             ancestors.extend(node.get_ancestors())
    
#     for ancestor in ancestors:
#         if ancestor.name:
#             ancestors_and_nodes.append(ancestor.name)
    
#     unique_nodes = list(set(ancestors_and_nodes))
#     print(f"✓ Nodos únicos a mantener: {len(unique_nodes)}")
    
#     tree.prune(unique_nodes, preserve_branch_length=True)
#     tree.write(outfile="pruned_tree.newick", format=1)
    
#     return tree

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
def get_otu_tree(otu_tax_table):
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

# def get_otu_tree(otu_tax_table):
#     start = time.time()
#     tree = get_pruned_tree(otu_tax_table)
    
#     print(f"Creando índice del árbol podado...")
#     node_index = {node.name: node for node in tree.traverse() if node.name}
#     print(f"✓ Índice creado: {len(node_index)} nodos")
    
#     # CALCULAR PROFUNDIDAD MÁXIMA DEL ÁRBOL (UNA VEZ)
#     print(f"Calculando profundidad máxima del árbol...")
#     max_depth = max(tree.get_distance(leaf) for leaf in tree.iter_leaves())
#     print(f"✓ Profundidad máxima: {max_depth}")
    
#     otu_dictionary = get_otu_dictionary(otu_tax_table)
#     print(f"Procesando {len(otu_dictionary)} OTUs...")
    
#     for i, (otu_id, taxonomy) in enumerate(otu_dictionary.items()):
#         if i % 5000 == 0:
#             print(f"  OTU {i}/{len(otu_dictionary)}...")
        
#         non_nan_values, nan_count = get_nan_non_nan_values(taxonomy)
        
#         if not non_nan_values:
#             continue
        
#         most_specific = non_nan_values[-1]
#         node = node_index.get(most_specific)
        
#         if node:
#             # Calcular distancia faltante hasta el nivel máximo
#             current_depth = tree.get_distance(node)
#             missing_depth = max_depth - current_depth
            
#             # SIEMPRE agregar como hijo (sin if/else)
#             node.add_child(name=otu_id, dist=missing_depth)
    
#     print(f"✓ OTUs agregados: {time.time() - start:.2f}s")
#     tree.write(outfile="otu_tree.newick", format=1)
#     return tree
