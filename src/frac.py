from ete3 import Tree
from decimal import Decimal, getcontext

# Gives a precise number of decimals for the final
# result. 20 digits of precision are determined.
getcontext().prec = 20

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
        count_community_1 = int(columns[1])
        count_community_2 = int(columns[2])
        
        if count_community_1 > 0 and count_community_2 > 0:
            taxon_communities[taxon_name] = (header[1], header[2])
        elif count_community_1 > 0 and count_community_2 == 0:
            taxon_communities[taxon_name] = (header[1],)
        else:
            taxon_communities[taxon_name] = (header[2],)
    
    return taxon_communities

# This function puts a new feature on the tree. In this
# case, it refers to the community feature. The output
# is a tree with leaves and the featured community
# they belong to.
def put_community_feature_on_tree(tree, lines):
    taxon_communities = get_processed_otu_table(lines)
    
    for leaf in tree.iter_leaves():
        if leaf.name in taxon_communities:
            leaf.community = taxon_communities[leaf.name]
    
    return tree

# This function gets all the branches that are unique.
# It generates a Python list that will be filled with
# the lengths of said unique branches. This function
# is not considering the root of the tree.
def get_unique_branches(tree):
    unique_branches = []
    community_list = []

    for node in tree.traverse():
        if node.is_root(): 
            continue

    for node in tree.traverse():
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
        
        first_community[columns[0]] = int(columns[1])
        second_community[columns[0]] = int(columns[2])
    
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
