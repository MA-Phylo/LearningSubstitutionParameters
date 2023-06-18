import copy
import numpy as np
import networkx as nx
import operator
import random
from eigenvalue_decomposition import eigenvalue_decomposition
from leaves_data import leaf_distribution
from leaves_data import matrix_from_leaves
from my_classes import MM
from my_classes import SepReconInput
from typing import NamedTuple

# Recover a mutation matrix previously computed
def obtain_matrix(Results, source, target):
    for matrix in Results:
        if matrix.source == source and matrix.target==target:
            return matrix.matrix

# Subroutine: LeafRecon
def LeafRecon(leaves_sequences, tree, a, U, S, net, Results, SR_input, Internal_Nodes_Distr):
    W = {a} #Current leaf
    U = U.union({a}) # Set of uncovered leaves

    #Tree's depth
    Delta = 0
    for clade in tree.get_terminals():
        p=tree.get_path(clade)
        Delta=max(Delta, len(p))

    # Subset of nodes at distance at most Delta from leaf a
    B = set()
    for node in net.nodes():
        sp=nx.shortest_path_length(net, source=a, target=node)
        if (sp > 0 and sp <= Delta):
            B.add(node)

    # Outside neighbors of W without using edges in S
    N = set()
    aux_net = nx.Graph(net)
    for edge in S:
        aux_net.remove_edge(edge[0], edge[1])
    for node in nx.node_boundary(aux_net, W):
        N.add(node)

    while B.intersection(N) != set():
        intersect = B.intersection(N)

        # Pick any vertex in the intersection computed
        r = random.choice(tuple(intersect))

        # Let (ro, r) be the edge with endpoint r in the path from a to r
        path = nx.shortest_path(net,a,r)
        r0 = path[-2]

        # If r0 = a, we set P_a_ro to be the identity matrix
        if (r0 == a):
            Results.append(MM(a, r0, np.identity(4)))

        # If r is not a leaf
        if (net.degree(r) > 1):

            # let (r,r1) and (r, r2) be the other two edges out of r
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r,r0)
            r1 = list(aux_net.edges(r))[0][1]
            r2 = list(aux_net.edges(r))[1][1]

            # Find the closest leaf b connected to r1 by a path not using (r0, r)
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r0,r)
            sp = nx.shortest_path_length(aux_net, r1)
            sp = {key:val for key, val in sp.items() if aux_net.degree(key) == 1}
            sp = sorted(sp.items(), key=operator.itemgetter(1))
            b = sp[0][0]

            # Find the closest leaf c connected to r2 by a path not using (r0, r)
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r0,r)
            sp = nx.shortest_path_length(aux_net, r2)
            sp = {key:val for key, val in sp.items() if aux_net.degree(key) == 1}
            sp = {key:val for key, val in sp.items() if key != b}
            sp = sorted(sp.items(), key=operator.itemgetter(1))
            c = sp[0][0]

            # Eigenvalue decomposition to obtain to P_a_r
            P_a_r = eigenvalue_decomposition(leaves_sequences, a, b, c)
            Results.append(MM(a, r, P_a_r))
            P_a_r = obtain_matrix(Results, a, r)

        # Otherwise r is a leaf
        else:
            P_a_r=matrix_from_leaves(leaves_sequences, a, r)
            U = U.union({r})

        # Compute P_ro_r
        P_a_r0 = obtain_matrix(Results, a, r0)
        P_r0_r = np.matmul(np.linalg.inv(P_a_r0), P_a_r)
        Results.append(MM(r0, r, P_r0_r))

        # Compute  π_r
        π_a = leaf_distribution(leaves_sequences, a)[0]
        π_r = np.matmul(π_a, P_a_r)
        if (net.degree(r) > 1):
            Internal_Nodes_Distr[r] = π_r

        # Obtain P_r_ro
        if (net.degree(r0) == 1):
            π_r0 = leaf_distribution(leaves_sequences, r0)[0]
        else:
            π_r0 = Internal_Nodes_Distr[r0]
        P_r_r0 = np.matmul(np.matmul(np.linalg.inv(np.diag(π_r)),np.matrix.transpose(P_r0_r)), np.diag(π_r0))
        Results.append(MM(r, r0, P_r_r0))

        # Obtain P_r_a (if r0 != a)
        if (r0 != a):
            P_r_a = np.matmul(np.matmul(np.linalg.inv(np.diag(π_r)),np.matrix.transpose(P_a_r)), np.diag(π_a))
            Results.append(MM(r, a, P_r_a))

        # Update W
        W = W.union({r})

        # If r is not a leaf and the distance between a and r is Delta
        if (net.degree(r) > 1 and nx.shortest_path_length(net, source=a, target=r) == Delta):
            if ((r, r1) not in S):
                t_prime = 1 + len(S)
                SR_input.append(SepReconInput(t_prime, b, a, r1, r))
                S.add((r,r1))
            if ((r, r2) not in S):
                t_prime = 1 + len(S)
                SR_input.append(SepReconInput(t_prime, c, a, r2, r))
                S.add((r,r2))

        # Recompute B (subset of nodes at distance at most Delta from leaf a)
        B = set()
        for node in net.nodes():
            sp=nx.shortest_path_length(net, source=a, target=node)
            if (sp > 0 and sp <= Delta):
                B.add(node)

        # Recompute N (Outside neighbors of W without using edges in S)
        N = set()
        aux_net = nx.Graph(net)
        for edge in S:
            aux_net.remove_edge(edge[0], edge[1])
        for node in nx.node_boundary(aux_net, W):
            N.add(node)
        
    return U, Results, SR_input, Internal_Nodes_Distr
