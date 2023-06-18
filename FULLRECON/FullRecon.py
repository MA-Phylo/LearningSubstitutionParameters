# Import modules and functions
import networkx as nx
import numpy as np
import pylab
import random
import sys
from Bio import SeqIO, Phylo
from collections import Counter
from LeafRecon import LeafRecon
from leaves_data import leaf_distribution
from io import StringIO
from Normalization import normalize_matrix
from Normalization import normalize_distribution
from SepRecon import SepRecon
from sklearn.preprocessing import normalize

def obtain_matrix(Results, source, target):
    for matrix in Results:
        if matrix.source == source and matrix.target==target:
            return matrix


def FullRecon(sequences, tree):
    ###--------------------------------------------------------------------------###
    ### Preparing the data

    # Read sequences (FASTA file)
    path_s = sequences
    sequences = [i for i in SeqIO.parse(path_s, 'fasta')]
    n_sequences = len(sequences) # Number of sequences

    # Read phylogenetic tree (Newick tree format)
    path_t = tree
    tree_file = open(path_t, "r")
    tree = tree_file.read()
    tree = Phylo.read(StringIO(tree), "newick")

    # Change nodes names
    for idx, clade in enumerate(tree.get_nonterminals()):
        clade.name = "Node_" + str(idx) if idx > 0 else "Root"
    # Change leaves names
    for idx, clade in enumerate(tree.get_terminals()):
        clade.name = "Leaf_" + clade.name

    # Map the tree's leaves to their corresponding DNA sequences
    leaves_sequences = {}
    for i in range (n_sequences):
        leaves_sequences[list(tree.get_terminals())[i]] = sequences[i].seq

    net = Phylo.to_networkx(tree) # Tree to network

    ###--------------------------------------------------------------------------###
    ### Phase 0: Initialization

    L = set(tree.get_terminals()) # Set of tree's leaves
    U = set() # Set of uncovered leaves
    S = set() # Set of separator edges
    w = random.choice(tuple(L)) #Pick randomly one leaf
    t = 0
    a = [0] * len(tree.get_terminals())#Vector of leaves
    a[t] = w
    Results = [] # Storing the mutation matrices
    Internal_Nodes_Distr = {}

    ###--------------------------------------------------------------------------###
    ### Phase 1: Subtree reconstruction

    SR_input = [] # Phase 2 input
    while U != L:
        U, Results, SR_input, Internal_Nodes_Distr = LeafRecon(leaves_sequences, tree, a[t], U, S, net, Results, SR_input, Internal_Nodes_Distr)
        t = t + 1
        for leaf in SR_input:
            if leaf.t == t:
                a[t] = leaf.a
    ###--------------------------------------------------------------------------###
    ### Phase 2: Patching up
    for nodes in SR_input:
        Results = SepRecon(leaves_sequences, nodes, Results, Internal_Nodes_Distr, net)

    ###--------------------------------------------------------------------------###

    ###--------------------------------------------------------------------------###
    ### Phase ∞: Normalization
    for MutationMatrix in Results:
        MutationMatrix.matrix = normalize_matrix(MutationMatrix.matrix)

    ###--------------------------------------------------------------------------###

    final_nodes_distributions = {}
    for node in net.nodes():
        if (net.degree(node) == 1):
            final_nodes_distributions[node] = leaf_distribution(leaves_sequences, node)[0]
        else:
            final_nodes_distributions[node] = Internal_Nodes_Distr[node]

    for k, v in final_nodes_distributions.items():
        final_nodes_distributions[k] = normalize_distribution(v)

    # Filtering
    Final_Results = []
    for edge in net.edges():
        res = obtain_matrix(Results, edge[0], edge[1])
        Final_Results.append(res)

    return Final_Results, final_nodes_distributions, net
