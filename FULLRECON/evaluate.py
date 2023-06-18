import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import sys
from scipy.linalg import logm, expm
from FullRecon import FullRecon
from simulate import matrix_generation
from Bio import Phylo, SeqIO
from io import StringIO
from my_classes import MM
from sympy import symbols, Eq, solve
from scipy.linalg import expm
from collections import Counter, OrderedDict
import itertools
import math

# Hyperparameter r
repetitions = 30

tree = sys.argv[1]
length = sys.argv[2]

# Simulation
real_matrices, sequences_file_name, tree_d = matrix_generation(tree, length)

# Matrices estimation with FullRecon
for r in range(repetitions):

      if r == 0:
          Results, final_nodes_distributions, net = FullRecon(sequences_file_name, tree)

      else:
          Res, distr, _ = FullRecon(sequences_file_name, tree)
          for MutationMatrixA in Res:
              for MutationMatrixB in Results:
                  if (str(MutationMatrixA.source) == str(MutationMatrixB.source) and str(MutationMatrixA.target)==str(MutationMatrixB.target)):
                      MutationMatrixB.matrix = MutationMatrixB.matrix + MutationMatrixA.matrix
          for k in final_nodes_distributions.keys():
              for k_ in distr.keys():
                  if (str(k) == str(k_)):
                    final_nodes_distributions[k] += distr[k_]

# Averaging the node distribution estimations
for node in final_nodes_distributions.keys():
    final_nodes_distributions[node] = 1/repetitions * final_nodes_distributions[node]

length_results = len(Results)
average_frob_error = 0
average_branch_absolute_error = 0

# Results for each transition matrix
for i in range(length_results):
    # Averaging the estimations
    Results[i].matrix = 1/repetitions * Results[i].matrix
    assert(str(real_matrices[i].source) == str(Results[i].source))
    assert(str(real_matrices[i].target) == str(Results[i].target))
    # Simulated matrix
    print(real_matrices[i].source, real_matrices[i].target, "Simulated Matrix")
    print(" ")
    print(real_matrices[i].matrix)
    print(" ")
    # Estimated matrix
    print(Results[i].source, Results[i].target, "Estimated Matrix")
    print(" ")
    print(Results[i].matrix)
    print(" ")
    # Error matrix
    error_matrix = real_matrices[i].matrix - Results[i].matrix
    print("Error matrix")
    print(" ")
    print(error_matrix)
    print(" ")
    # Frobenius norm
    error_frob_matrix = np.linalg.norm(error_matrix, 'fro')
    print("Error matrix Frobenius norm", error_frob_matrix)
    π = final_nodes_distributions[Results[i].source]
    π_ = final_nodes_distributions[Results[i].target]
    Dπ = np.diag(π)
    Dπ_ = np.diag(π_)
    M = Results[i].matrix
    print(" ")
    print("******************************************************")
    print(" ")
    # Branch length estimation
    estimated_l = 1/4*(-np.log(((np.sqrt(np.linalg.det(Dπ)))*(np.linalg.det(M)))/(np.sqrt(np.linalg.det(Dπ_)))))
    current_edge = list(net.edges())[i]
    branch_error = np.abs(estimated_l - current_edge[1].branch_length)
    print("True branch length: ", current_edge[1].branch_length)
    print("Estimated branch length: ", estimated_l)
    average_frob_error += error_frob_matrix
    average_branch_absolute_error += branch_error
    print(" ")
    print("------------------------------------------------------------------------------------------------------------------")
    print(" ")


print("Average branch length absolute error", (average_branch_absolute_error / length_results))
print("Average error matrix Frobenius norm", average_frob_error / length_results)

#--------------------------------------------------------------------------------------------------------------

# Computing the log-likelihood of the phylogenetic tree given the estimated parameters
# (code provided by Martí Cortada)

path_s = sequences_file_name
sequences = [i for i in SeqIO.parse(path_s, 'fasta')]
EM_leave_sequences = []
for i in sequences:
    seq = ""
    for j in i.seq:
        if j == "A":
            seq += "0"
        elif j == "G":
            seq += "1"
        elif j == "C":
            seq += "2"
        elif j == "T":
            seq += "3"
        else:
            raise ValueError("Invalid nucleotide")
    EM_leave_sequences.append(seq)

sequence_length = len(EM_leave_sequences[0])
number_of_sequences = len(EM_leave_sequences)
occurrences = []
for i in range(sequence_length):
    _ = ''
    for j in range(number_of_sequences):
        _ += EM_leave_sequences[j][i]
    occurrences.append(_)
c = Counter(occurrences)
u_i = OrderedDict(sorted(c.items()))


n_int = 2
hidden_combinations = itertools.product([0,1,2,3], repeat=n_int)
sorted_hidden_combinations = sorted(hidden_combinations)
for i in range(len(sorted_hidden_combinations)):
    sorted_hidden_combinations[i] = ''.join([str(s) for s in sorted_hidden_combinations[i]])

states = list(itertools.product(list(sorted_hidden_combinations), list(u_i.keys())))
states = [i[0]+i[1] for i in states]

r_distr = final_nodes_distributions[Results[0].source]
print(r_distr)

def map_node(n):
    if n.name == "Root":
        return 0
    if n.name == "Node_1":
        return 1
    if n.name == "Leaf_1":
        return 2
    if n.name == "Leaf_2":
        return 3
    if n.name == "Leaf_3":
        return 4
    if n.name == "Leaf_4":
        return 5


def log_likelihood(states, u_i, params, root_distribution, n_int):
    """
    #Compute the log-likehood of the tree
    """
    logL = 0
    for obs in states:
        p_i = 0
        observed_data = obs[n_int:]
        if observed_data in u_i:
            for state in states:
                if state[n_int:] == observed_data:
                    pi = root_distribution[int(state[0])]
                    for p in params:
                        u, v = p.source, p.target
                        u = map_node(u)
                        v = map_node(v)
                        pi *= p.matrix[int(state[u]), int(state[v])]
                    p_i += pi
        logL += (u_i[observed_data] * math.log(p_i))
    return logL

print("Log-likelihood: ", log_likelihood(states, u_i, Results, r_distr, 2))

Phylo.draw(tree_d)
