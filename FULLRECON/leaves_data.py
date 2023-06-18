import numpy as np
from collections import Counter
np.set_printoptions(suppress=True)

# Given a leaf and its sequence, it returns its nucleotide distribution.
def leaf_distribution(leaves_sequences, a):
    seq_a = leaves_sequences[a]
    #Counter for sequence A
    len_sequence = len(seq_a) # Length of sequences
    v = []
    for i in range (len_sequence):
        str = '' + seq_a[i]
        v.append(str)
    counter_a = Counter(v) # Count ocurrences

    distribution = np.ones((1,4))
    distribution[0, 0] = counter_a.get("A", 0) / len_sequence
    distribution[0, 1] = counter_a.get("G", 0) / len_sequence
    distribution[0, 2] = counter_a.get("C", 0) / len_sequence
    distribution[0, 3] = counter_a.get("T", 0) / len_sequence
    return distribution

# Given two leaves a,b and its sequences, it returns the matrix Pab
def matrix_from_leaves(leaves_sequences, leaf_a, leaf_b):
    seq_a = leaves_sequences[leaf_a]
    seq_b = leaves_sequences[leaf_b]

    #Counter for sequence A
    len_sequence = len(seq_a) # Length of sequences
    v = []
    for i in range (len_sequence):
        str = '' + seq_a[i]
        v.append(str)
    counter_a = Counter(v) # Count ocurrences

    # Counter for sequences A and B
    len_sequence = len(seq_a) # Length of sequences
    v = []
    for i in range (len_sequence):
        str = '' + seq_a[i] + seq_b[i]
        v.append(str)
    counter_ab = Counter(v) # Count ocurrences

    # P_ab matrix estimation
    P_a_b = np.zeros((4,4))
    P_a_b[0,0] = counter_ab.get("AA", 0) / counter_a.get("A", 0)
    P_a_b[0,1] = counter_ab.get("AG", 0) / counter_a.get("A", 0)
    P_a_b[0,2] = counter_ab.get("AC", 0) / counter_a.get("A", 0)
    P_a_b[0,3] = counter_ab.get("AT", 0) / counter_a.get("A", 0)
    P_a_b[1,0] = counter_ab.get("GA", 0) / counter_a.get("G", 0)
    P_a_b[1,1] = counter_ab.get("GG", 0) / counter_a.get("G", 0)
    P_a_b[1,2] = counter_ab.get("GC", 0) / counter_a.get("G", 0)
    P_a_b[1,3] = counter_ab.get("GT", 0) / counter_a.get("G", 0)
    P_a_b[2,0] = counter_ab.get("CA", 0) / counter_a.get("C", 0)
    P_a_b[2,1] = counter_ab.get("CG", 0) / counter_a.get("C", 0)
    P_a_b[2,2] = counter_ab.get("CC", 0) / counter_a.get("C", 0)
    P_a_b[2,3] = counter_ab.get("CT", 0) / counter_a.get("C", 0)
    P_a_b[3,0] = counter_ab.get("TA", 0) / counter_a.get("T", 0)
    P_a_b[3,1] = counter_ab.get("TG", 0) / counter_a.get("T", 0)
    P_a_b[3,2] = counter_ab.get("TC", 0) / counter_a.get("T", 0)
    P_a_b[3,3] = counter_ab.get("TT", 0) / counter_a.get("T", 0)

    return P_a_b
    
