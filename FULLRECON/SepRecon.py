import numpy as np
from leaves_data import matrix_from_leaves
from my_classes import MM
from leaves_data import leaf_distribution

# Recover a transition matrix previously computed
def obtain_matrix(Results, source, target):
    for matrix in Results:
        if matrix.source == source and matrix.target==target:
            return matrix.matrix

# Subroutine: SepRecon
def SepRecon(leaves_sequences, nodes, Results, Internal_Nodes_Distr, net):

    # Estimate Paa'
    P_a_a_prime = matrix_from_leaves(leaves_sequences, nodes.a, nodes.a_prime)

    # Recover Paw
    if (nodes.a == nodes.w):
        P_a_w = matrix_from_leaves(leaves_sequences, nodes.a, nodes.w)
    else:
        P_a_w = obtain_matrix(Results, nodes.a, nodes.w)

    # Recover Pw'a'
    P_w_prime_a_prime = obtain_matrix(Results, nodes.w_prime, nodes.a_prime)

    # Compute Pww'
    P_w_w_prime = np.matmul(np.matmul(np.linalg.inv(P_a_w),P_a_a_prime), np.linalg.inv(P_w_prime_a_prime))
    Results.append(MM(nodes.w, nodes.w_prime, P_w_w_prime))

    # Compute Pw'w
    if (net.degree(nodes.w) == 1):
        π_w = leaf_distribution(leaves_sequences, nodes.w)[0]
    else:
        π_w = Internal_Nodes_Distr[nodes.w]
    π_w_prime = Internal_Nodes_Distr[nodes.w_prime]
    P_w_prime_w = np.matmul(np.matmul(np.linalg.inv(np.diag(π_w_prime)),np.matrix.transpose(P_w_w_prime)), np.diag(π_w))
    Results.append(MM(nodes.w_prime, nodes.w, P_w_prime_w))

    return Results
