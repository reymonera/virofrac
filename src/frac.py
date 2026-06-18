from ete3 import Tree
import networkx as nx
from decimal import Decimal, getcontext

# Gives a precise number of decimals for the final
# result. 20 digits of precision are determined.
getcontext().prec = 20

####################################################
#           PHYLOGENETIC - BASED                   #
####################################################

# This function processes the lines from the OTU
# table that is included in the input files. Mainly
# it works by determining community assignment per
# taxa. It ignores the header.
def get_processed_otu_table(lines):
    taxon_communities = {}
    header = None
    
    for i, line in enumerate(lines):
        columns = line.split('\t')
        
        if i == 0:
            header = columns
            continue
        
        taxon_name = columns[0]
        count_community_1 = float(columns[1])
        count_community_2 = float(columns[2])
        
        if count_community_1 > 0 and count_community_2 > 0:
            taxon_communities[taxon_name] = (header[1], header[2])
        elif count_community_1 > 0 and count_community_2 == 0:
            taxon_communities[taxon_name] = (header[1],)
        elif count_community_2 > 0:
            taxon_communities[taxon_name] = (header[2],)
    
    return taxon_communities

# This function puts a new feature on the tree. In this
# case, it refers to the community feature. The output
# is a tree with leaves and the featured community
# they belong to.

# def put_community_feature_on_tree(tree, lines):
#     taxon_communities = get_processed_otu_table(lines)
    
#     for leaf in tree.iter_leaves():
#         if leaf.name in taxon_communities:
#             leaf.community = taxon_communities[leaf.name]
    
#     return tree

def put_community_feature_on_tree(tree, lines):
    taxon_communities = get_processed_otu_table(lines)
    for leaf in tree.iter_leaves():
        leaf.community = taxon_communities.get(leaf.name, ())
    return tree

# This function gets all the branches that are unique.
# It generates a Python list that will be filled with
# the lengths of said unique branches. This function
# is not considering the root of the tree.
def get_unique_branches(tree):
    unique_branches = []
    community_list = []

    # for node in tree.traverse():
    #     if node.is_root(): 
    #         continue

    for node in tree.traverse():
        if node.is_root(): 
            continue
        descendants = node.get_leaves()
        for leaf in descendants:
            community_list.append(leaf.community)
        if len(set(community_list)) == 1 and len(community_list[0]) == 1:
            unique_branches.append(node)
        community_list = []
    return(unique_branches)

# getSumUniqueBranchLengths just sums all the 
# unique lengths that were in the previous
# created list.
def get_sum_unique_branch_lengths(tree):
    unique_branches = get_unique_branches(tree)
    unique_lengths = []
    for node in unique_branches:
        unique_lengths.append(node.dist)
    return(sum(unique_lengths))

# getSumAllBranchLengths just sums all the lengths.
# This function does not consider the root.
def get_sum_all_branch_lengths(tree):
    all_lengths = []

    for node in tree.traverse():
        if node.is_root():
            continue
        all_lengths.append(node.dist)

    return(sum(all_lengths))

# getCommunityDictionaries is a function that will
# create dictionaries for the communities that are
# being compared with each other.
def get_community_dictionaries(lines):
    first_community = {}
    second_community = {}
    community_dictionaries = []
    for i, line in enumerate(lines):

        if i == 0:
            continue
        
        columns = line.split('\t')
        
        first_community[columns[0]] = float(columns[1])
        second_community[columns[0]] = float(columns[2])
    
    community_dictionaries.append(first_community)
    community_dictionaries.append(second_community)

    return community_dictionaries

# getCommunityProportionsDictionaries is a function
# that is assigning relative proportions per taxa 
# in each community. 
def get_community_proportions_dictionaries(lines):
    communities_proportions = []
    community_dictionaries = get_community_dictionaries(lines)

    first_community = community_dictionaries[0]
    second_community = community_dictionaries[1]
    
    all_first_community = sum(first_community.values())
    all_second_community = sum(second_community.values())
    
    first_community.update((x, y/all_first_community) for x, y in first_community.items())
    second_community.update((x, y/all_second_community) for x, y in second_community.items())
    
    communities_proportions.append(first_community)
    communities_proportions.append(second_community)
    
    return communities_proportions

# This function outputs the Unweighted Unifrac value
# between 2 communities.
def get_unweighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    sum_unique_branches = Decimal(get_sum_unique_branch_lengths(tree))
    sum_all_branches = Decimal(get_sum_all_branch_lengths(tree))

    return sum_unique_branches/sum_all_branches

# This function outputs the Normalized Weighted Unifrac 
# value between 2 communities.
def get_unnormalized_weighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    first_community = get_community_proportions_dictionaries(lines)[0]
    second_community = get_community_proportions_dictionaries(lines)[1]

    sum_weighted_branches = []
    
    for node in tree.traverse():
        if node.is_root():
            continue
            
        descendants = node.get_leaves()
        sum_first_proportions = []
        sum_second_proportions = []
        
        for leaf in descendants:
            sum_first_proportions.append(first_community[leaf.name])
            sum_second_proportions.append(second_community[leaf.name])
        
        p_i = sum(sum_first_proportions)
        q_i = sum(sum_second_proportions)
        
        sum_weighted_branches.append(node.dist * abs(p_i - q_i))

    return(Decimal(sum(sum_weighted_branches)))

# This function outputs the Normalized Weighted Unifrac 
# value between 2 communities.
def get_normalized_weighted_unifrac(community_tree, lines):
    tree = put_community_feature_on_tree(community_tree, lines)

    first_community = get_community_proportions_dictionaries(lines)[0]
    second_community = get_community_proportions_dictionaries(lines)[1]

    sum_weighted_branches = []
    sum_all_weighted_branches = []
    
    for node in tree.traverse():
        if node.is_root():
            continue
            
        descendants = node.get_leaves()
        sum_first_proportions = []
        sum_second_proportions = []
        
        for leaf in descendants:
            sum_first_proportions.append(first_community[leaf.name])
            sum_second_proportions.append(second_community[leaf.name])
        
        p_i = sum(sum_first_proportions)
        q_i = sum(sum_second_proportions)
        
        sum_weighted_branches.append(node.dist * abs(p_i - q_i))
        sum_all_weighted_branches.append(node.dist * (p_i + q_i))
    
    return(Decimal(sum(sum_weighted_branches)/sum(sum_all_weighted_branches)))

####################################################
#              NETWORK - BASED
####################################################
def get_distances_from_edges(edge_list):
    distance_list = [1 - x for x in edge_list]
    return distance_list

def get_similarities_from_edges(edge_list):
    similarity_list = [x for x in edge_list]
    return similarity_list

def get_all_weights_from_edges(network, community_A, community_B):            
    all_weights = []
    for node1, node2, weight in network.edges(data='weight'):
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        # Edge is relevant if BOTH nodes belong to at least one of the communities
        if (n1_in_A or n1_in_B) and (n2_in_A or n2_in_B):
            all_weights.append(weight)
    
    return all_weights

# Get volume for any imput community, which is always the first one
def get_volume_of_first_sample(network, community_A, community_B):            
    edge_weights_of_sample_A = []
    for node1, node2, weight in network.edges(data='weight'):
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        # Edge must be relevant (both nodes belong to at least one community)
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        # Volume of A = edges where at least one node belongs to A
        if n1_in_A or n2_in_A:
            edge_weights_of_sample_A.append(weight)
    
    return sum(edge_weights_of_sample_A)

def get_monochromatic_edges(network, community_A, community_B):
    monochromatic_weights = []
    for node1, node2, weight in network.edges(data='weight'):
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        # First: Edge must be relevant (same filter as all_edges)
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        # Second: Check if monochromatic
        # Both nodes in A only (neither in B)
        both_A_only = (n1_in_A and n2_in_A and not n1_in_B and not n2_in_B)
        # Both nodes in B only (neither in A)
        both_B_only = (n1_in_B and n2_in_B and not n1_in_A and not n2_in_A)
        
        #if (n1_in_A and n2_in_A) or (n1_in_B and n2_in_B):
        if both_A_only or both_B_only:
            monochromatic_weights.append(weight)
    
    return monochromatic_weights

def get_bichromatic_edges(network, community_A, community_B):
    bichromatic_weights = []
    for node1, node2, weight in network.edges(data='weight'):
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        # First: Edge must be relevant (same filter as all_edges)
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue

        both_in_both = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
        original_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B)
        
        if original_bichromatic or both_in_both:
            bichromatic_weights.append(weight)
    
    return bichromatic_weights

def get_total_for_first_community(network, community_A, community_B):
    node_abundances = {}
    
    for node1, node2 in network.edges():
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        node_abundances[node1] = network.nodes[node1]['communities'].get(community_A, 0)
        node_abundances[node2] = network.nodes[node2]['communities'].get(community_A, 0)
    
    return sum(node_abundances.values())

def get_unweighted_netunifrac(network, community_A, community_B):
    monochromatic_edges = get_monochromatic_edges(network, community_A, community_B)
    all_edges = get_all_weights_from_edges(network, community_A, community_B)

    monochromatic_distances = get_distances_from_edges(monochromatic_edges)
    all_distances = get_distances_from_edges(all_edges)

    sum_monochromatic_distances = sum(monochromatic_distances)
    sum_all_distances = sum(all_distances)

    if sum_all_distances == 0:
        # No shared structure = maximum distance (completely different communities)
        return 1.0

    # NetUniFrac
    unweighted_netunifrac = sum_monochromatic_distances/sum_all_distances
    unweighted_netunifrac = len(monochromatic_distances)/len(all_distances)

    return unweighted_netunifrac

# def get_weighted_netunifrac(network, community_A, community_B):
#     list_num = []
#     list_dem = []

#     total_community_A = get_total_for_first_community(network, community_A, community_B)
#     total_community_B = get_total_for_first_community(network, community_B, community_A)

#     for node1, node2, weight in network.edges(data='weight'):
#         n1_in_A = community_A in network.nodes[node1]['communities']
#         n1_in_B = community_B in network.nodes[node1]['communities']
#         n2_in_A = community_A in network.nodes[node2]['communities']
#         n2_in_B = community_B in network.nodes[node2]['communities']

#         node1_A = network.nodes[node1]['communities'].get(community_A, 0)
#         node1_B = network.nodes[node1]['communities'].get(community_B, 0)
#         node2_A = network.nodes[node2]['communities'].get(community_A, 0)
#         node2_B = network.nodes[node2]['communities'].get(community_B, 0)

#         if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
#             continue

#     # Go through all distances
#     # Multiply them by the proportions of node wieghts
#         if total_community_A == 0 or total_community_B == 0:
#             continue
#         else:
#             # The (1 - weight) here is doing the distance transformation, instead of using the function assigned
#             # for this transformation as in the get_weighted_netunifrac() function
#             num = (1 - weight)*abs((node1_A + node2_A)/total_community_A - (node1_B + node2_B)/total_community_B)
#             dem = (1 - weight)*abs((node1_A + node2_A)/total_community_A + (node1_B + node2_B)/total_community_B)
#         # Save this on a list
#             list_num.append(num)
#             list_dem.append(dem)
#     # Return the division
#     # Handle empty lists
#     if sum(list_dem) == 0:
#         weighted_netunifrac = 1
#     else:
#         weighted_netunifrac = sum(list_num) / sum(list_dem)
    
#     return weighted_netunifrac

# Weighted but with monochromatic edges
def get_weighted_netunifrac(network, community_A, community_B):
    list_num = []
    list_dem = []
    
    total_community_A = get_total_for_first_community(network, community_A, community_B)
    total_community_B = get_total_for_first_community(network, community_B, community_A)
    
    if total_community_A == 0 or total_community_B == 0:
        return 1.0
    
    for node1, node2, weight in network.edges(data='weight'):
        n1_in_A = community_A in network.nodes[node1]['communities']
        n1_in_B = community_B in network.nodes[node1]['communities']
        n2_in_A = community_A in network.nodes[node2]['communities']
        n2_in_B = community_B in network.nodes[node2]['communities']
        
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        node1_A = network.nodes[node1]['communities'].get(community_A, 0)
        node1_B = network.nodes[node1]['communities'].get(community_B, 0)
        node2_A = network.nodes[node2]['communities'].get(community_A, 0)
        node2_B = network.nodes[node2]['communities'].get(community_B, 0)
        
        prop_diff = abs((node1_A + node2_A)/total_community_A - 
                        (node1_B + node2_B)/total_community_B)
        prop_sum  = abs((node1_A + node2_A)/total_community_A + 
                        (node1_B + node2_B)/total_community_B)
        
        # Monochromatic check inline
        both_A_only = (n1_in_A and n2_in_A and not n1_in_B and not n2_in_B)
        both_B_only = (n1_in_B and n2_in_B and not n1_in_A and not n2_in_A)
        is_mono     = both_A_only or both_B_only
        
        if is_mono:
            list_num.append((1 - weight) * prop_diff)
        
        list_dem.append((1 - weight) * prop_sum)
    
    if sum(list_dem) == 0:
        return 1.0
    
    return sum(list_num) / sum(list_dem)
    
def put_edge_community_data(network):
    edge_data = []
    for node1, node2, weight in network.edges(data='weight'):
        communities_1 = network.nodes[node1].get('communities', {})
        communities_2 = network.nodes[node2].get('communities', {})
        edge_data.append((communities_1, communities_2, weight))

    return edge_data

# This function performs a distance based on spectral clustering. It
# requires similarity instead of distance as an input. If the denominator
# is equal to 0, then it will automatically return 1.
def get_based_spectral_clustering_1(edge_data, community_A, community_B):
    sum_bichromatic = 0.0
    volume_A = 0.0
    volume_B = 0.0

    for communities_1, communities_2, weight in edge_data:
        n1_in_A = community_A in communities_1
        n1_in_B = community_B in communities_1
        n2_in_A = community_A in communities_2
        n2_in_B = community_B in communities_2

        # Edge must be relevant
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue

        # Volume A: at least one node in A
        if n1_in_A or n2_in_A:
            volume_A += weight

        # Volume B: at least one node in B
        if n1_in_B or n2_in_B:
            volume_B += weight

        # Bichromatic: nodes differ in membership, or both in both
        both_in_both = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
        original_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B)
        if original_bichromatic or both_in_both:
            sum_bichromatic += weight

    if volume_A == 0 or volume_B == 0:
        return 1.0

    return 1 - ((sum_bichromatic / volume_A + sum_bichromatic / volume_B) * 0.5)

def get_based_spectral_clustering(edge_data, community_A, community_B):
    
    # Calculate totals directly from edge_data
    total_A = 0.0
    total_B = 0.0
    for communities_1, communities_2, weight in edge_data:
        for communities in [communities_1, communities_2]:
            if community_A in communities:
                total_A += communities[community_A]
            if community_B in communities:
                total_B += communities[community_B]
    
    if total_A == 0 or total_B == 0:
        return 1.0
    
    sum_bichromatic = 0.0
    volume_A = 0.0
    volume_B = 0.0
    
    for communities_1, communities_2, weight in edge_data:
        n1_in_A = community_A in communities_1
        n1_in_B = community_B in communities_1
        n2_in_A = community_A in communities_2
        n2_in_B = community_B in communities_2
        
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        node1_A = communities_1.get(community_A, 0)
        node1_B = communities_1.get(community_B, 0)
        node2_A = communities_2.get(community_A, 0)
        node2_B = communities_2.get(community_B, 0)
        
        delta1 = (node1_A - node1_B)
        delta2 = (node2_A - node2_B)

        factor = weight * ((1 - delta1 * delta2) / 2)
        #print("delta1 is: ", delta1)
        #print("delta2 is: ", delta2)
        #print("factor is: ", factor)
        
        if n1_in_A or n2_in_A:
            volume_A += weight
        if n1_in_B or n2_in_B:
            volume_B += weight
        
        both_in_both   = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
        is_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B) or both_in_both
        
        if is_bichromatic:
            sum_bichromatic += factor
    
    vol_A_adj = volume_A #* total_A
    vol_B_adj = volume_B #* total_B
    
    if vol_A_adj == 0 or vol_B_adj == 0:
        return 1.0
    
    return 1 - 0.5 * (sum_bichromatic / vol_A_adj + sum_bichromatic / vol_B_adj)

def get_based_spectral_clustering_2(edge_data, community_A, community_B):
    
    # Calculate totals directly from edge_data
    total_A = 0.0
    total_B = 0.0
    for communities_1, communities_2, weight in edge_data:
        for communities in [communities_1, communities_2]:
            if community_A in communities:
                total_A += communities[community_A]
            if community_B in communities:
                total_B += communities[community_B]
    
    if total_A == 0 or total_B == 0:
        return 1.0
    
    sum_bichromatic = 0.0
    volume_A = 0.0
    volume_B = 0.0
    
    for communities_1, communities_2, weight in edge_data:
        n1_in_A = community_A in communities_1
        n1_in_B = community_B in communities_1
        n2_in_A = community_A in communities_2
        n2_in_B = community_B in communities_2
        
        if not ((n1_in_A or n1_in_B) and (n2_in_A or n2_in_B)):
            continue
        
        node1_A = communities_1.get(community_A, 0)
        node1_B = communities_1.get(community_B, 0)
        node2_A = communities_2.get(community_A, 0)
        node2_B = communities_2.get(community_B, 0)
        
        if n1_in_A or n2_in_A:
            volume_A += weight
        if n1_in_B or n2_in_B:
            volume_B += weight
        
        both_in_both   = (n1_in_A and n1_in_B and n2_in_A and n2_in_B)
        is_bichromatic = (n1_in_A != n2_in_A or n1_in_B != n2_in_B) or both_in_both
        
        if is_bichromatic:
            # Numerator: edge weight × sum of node abundances in both communities
            node_abund_sum = node1_A + node1_B + node2_A + node2_B
            sum_bichromatic += weight * node_abund_sum
    
    # Denominator: volume × |total_A - total_B|
    diff = abs(total_A - total_B)
    vol_A_adj = volume_A * diff
    vol_B_adj = volume_B * diff
    
    if vol_A_adj == 0 or vol_B_adj == 0:
        return 1.0
    
    return 1 - 0.5 * (sum_bichromatic / vol_A_adj + sum_bichromatic / vol_B_adj)