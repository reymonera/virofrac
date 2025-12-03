import networkx as nx
import numpy as np
import pandas as pd

# Esta función debería de acepar el network basado en ani y en gene sharing
def get_count_table(count_table_path):
    table = pd.read_csv(count_table_path)
    indexed_count_table = table.set_index([0])

    return indexed_count_table

def get_network(matrix_file):
    matrix = np.loadtxt(matrix_file, delimiter=",", skiprows=1)
    network = nx.from_numpy_array(matrix)

    return network

def get_network_with_threshold(threshold, network):
    threshold_network = network.remove_edges_from([(n1, n2) for n1, n2, 
                               weight in network.edges(data="weight") 
                               if weight < threshold])
    
    return threshold_network

def set_community_atribute_on_nodes(network, count_table_path):
    count_table = get_count_table(count_table_path)
    
    for node in count_table.index:
        row = count_table.loc[node]
        communities = row[row > 0].to_dict() # Diccionario de comunidades y cantidades
        network.nodes[node]['communities'] = communities 
        # Ahora G.nodes[0]['comunidades'] = {'comunidad_A': 5, 'comunidad_B': 2}
    
    return network

def get_all_weights_from_edges(network, community_A, community_B):            
    all_weights = [weight for node1, node2, weight 
                   in network.edges(data='weight')
                   if (community_A in network.nodes[node1]['communities'] or
                   community_B in network.nodes[node1]['communities']) and
                   (community_B in network.nodes[node2]['communities'] or
                   community_B in network.nodes[node2]['communities'])]

    return all_weights
# Debería de retornar una lista [0.5, 1.2, 0.8]
# Solo ejes de esa comparación de pares

def get_monochromatic_edges(network, community_A, community_B):
    monochromatic_edges = []
    for node1, node2 in network.edges:
        in_A = (community_A in network.nodes[node1]['communities'] and 
                community_A in network.nodes[node2]['communities'])
        in_B = (community_B in network.nodes[node1]['communities'] and 
                community_B in network.nodes[node2]['communities'])
        
        if in_A != in_B:
            monochromatic_edges.append(network[node1][node2]['weight'])
    
    return monochromatic_edges
