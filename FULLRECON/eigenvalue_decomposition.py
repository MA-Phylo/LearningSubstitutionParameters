import numpy as np
from collections import Counter
from relabeling import permutation
np.set_printoptions(suppress=True)

def eigenvalue_decomposition(leaves_sequences, a, b, c):
    seq_a = leaves_sequences[a]
    seq_b = leaves_sequences[b]
    seq_c = leaves_sequences[c]

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

    # Counter for sequences A, B and C
    len_sequence = len(seq_a) # Length of sequences
    v = []
    for i in range (len_sequence):
        str = '' + seq_a[i] + seq_b[i] + seq_c[i]
        v.append(str)
    counter_abc = Counter(v) # Count ocurrences

    # P_a_b matrix estimation
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

    #P_a_b_U matrix estimation
    P_a_b_U = np.zeros((4,4))
    U = np.random.normal(0, 1, size=(4, 1))
    P_a_b_U[0, 0] = U[0]*(counter_abc.get("AAA", 0)/counter_a.get("A", 0)) + U[1]*(counter_abc.get("AAG", 0)/counter_a.get("A", 0)) + U[2]*(counter_abc.get("AAC", 0)/counter_a.get("A", 0)) + U[3]*(counter_abc.get("AAT", 0)/counter_a.get("A", 0))
    P_a_b_U[0, 1] = U[0]*(counter_abc.get("AGA", 0)/counter_a.get("A", 0)) + U[1]*(counter_abc.get("AGG", 0)/counter_a.get("A", 0)) + U[2]*(counter_abc.get("AGC", 0)/counter_a.get("A", 0)) + U[3]*(counter_abc.get("AGT", 0)/counter_a.get("A", 0))
    P_a_b_U[0, 2] = U[0]*(counter_abc.get("ACA", 0)/counter_a.get("A", 0)) + U[1]*(counter_abc.get("ACG", 0)/counter_a.get("A", 0)) + U[2]*(counter_abc.get("ACC", 0)/counter_a.get("A", 0)) + U[3]*(counter_abc.get("ACT", 0)/counter_a.get("A", 0))
    P_a_b_U[0, 3] = U[0]*(counter_abc.get("ATA", 0)/counter_a.get("A", 0)) + U[1]*(counter_abc.get("ATG", 0)/counter_a.get("A", 0)) + U[2]*(counter_abc.get("ATC", 0)/counter_a.get("A", 0)) + U[3]*(counter_abc.get("ATT", 0)/counter_a.get("A", 0))
    P_a_b_U[1, 0] = U[0]*(counter_abc.get("GAA", 0)/counter_a.get("G", 0)) + U[1]*(counter_abc.get("GAG", 0)/counter_a.get("G", 0)) + U[2]*(counter_abc.get("GAC", 0)/counter_a.get("G", 0)) + U[3]*(counter_abc.get("GAT", 0)/counter_a.get("G", 0))
    P_a_b_U[1, 1] = U[0]*(counter_abc.get("GGA", 0)/counter_a.get("G", 0)) + U[1]*(counter_abc.get("GGG", 0)/counter_a.get("G", 0)) + U[2]*(counter_abc.get("GGC", 0)/counter_a.get("G", 0)) + U[3]*(counter_abc.get("GGT", 0)/counter_a.get("G", 0))
    P_a_b_U[1, 2] = U[0]*(counter_abc.get("GCA", 0)/counter_a.get("G", 0)) + U[1]*(counter_abc.get("GCG", 0)/counter_a.get("G", 0)) + U[2]*(counter_abc.get("GCC", 0)/counter_a.get("G", 0)) + U[3]*(counter_abc.get("GCT", 0)/counter_a.get("G", 0))
    P_a_b_U[1, 3] = U[0]*(counter_abc.get("GTA", 0)/counter_a.get("G", 0)) + U[1]*(counter_abc.get("GTG", 0)/counter_a.get("G", 0)) + U[2]*(counter_abc.get("GTC", 0)/counter_a.get("G", 0)) + U[3]*(counter_abc.get("GTT", 0)/counter_a.get("G", 0))
    P_a_b_U[2, 0] = U[0]*(counter_abc.get("CAA", 0)/counter_a.get("C", 0)) + U[1]*(counter_abc.get("CAG", 0)/counter_a.get("C", 0)) + U[2]*(counter_abc.get("CAC", 0)/counter_a.get("C", 0)) + U[3]*(counter_abc.get("CAT", 0)/counter_a.get("C", 0))
    P_a_b_U[2, 1] = U[0]*(counter_abc.get("CGA", 0)/counter_a.get("C", 0)) + U[1]*(counter_abc.get("CGG", 0)/counter_a.get("C", 0)) + U[2]*(counter_abc.get("CGC", 0)/counter_a.get("C", 0)) + U[3]*(counter_abc.get("CGT", 0)/counter_a.get("C", 0))
    P_a_b_U[2, 2] = U[0]*(counter_abc.get("CCA", 0)/counter_a.get("C", 0)) + U[1]*(counter_abc.get("CCG", 0)/counter_a.get("C", 0)) + U[2]*(counter_abc.get("CCC", 0)/counter_a.get("C", 0)) + U[3]*(counter_abc.get("CCT", 0)/counter_a.get("C", 0))
    P_a_b_U[2, 3] = U[0]*(counter_abc.get("CTA", 0)/counter_a.get("C", 0)) + U[1]*(counter_abc.get("CTG", 0)/counter_a.get("C", 0)) + U[2]*(counter_abc.get("CTC", 0)/counter_a.get("C", 0)) + U[3]*(counter_abc.get("CTT", 0)/counter_a.get("C", 0))
    P_a_b_U[3, 0] = U[0]*(counter_abc.get("TAA", 0)/counter_a.get("T", 0)) + U[1]*(counter_abc.get("TAG", 0)/counter_a.get("T", 0)) + U[2]*(counter_abc.get("TAC", 0)/counter_a.get("T", 0)) + U[3]*(counter_abc.get("TAT", 0)/counter_a.get("T", 0))
    P_a_b_U[3, 1] = U[0]*(counter_abc.get("TGA", 0)/counter_a.get("T", 0)) + U[1]*(counter_abc.get("TGG", 0)/counter_a.get("T", 0)) + U[2]*(counter_abc.get("TGC", 0)/counter_a.get("T", 0)) + U[3]*(counter_abc.get("TGT", 0)/counter_a.get("T", 0))
    P_a_b_U[3, 2] = U[0]*(counter_abc.get("TCA", 0)/counter_a.get("T", 0)) + U[1]*(counter_abc.get("TCG", 0)/counter_a.get("T", 0)) + U[2]*(counter_abc.get("TCC", 0)/counter_a.get("T", 0)) + U[3]*(counter_abc.get("TCT", 0)/counter_a.get("T", 0))
    P_a_b_U[3, 3] = U[0]*(counter_abc.get("TTA", 0)/counter_a.get("T", 0)) + U[1]*(counter_abc.get("TTG", 0)/counter_a.get("T", 0)) + U[2]*(counter_abc.get("TTC", 0)/counter_a.get("T", 0)) + U[3]*(counter_abc.get("TTT", 0)/counter_a.get("T", 0))

    # Product
    F=np.matmul(P_a_b_U,np.linalg.inv(P_a_b))

    # Get eigenvectors
    [w,v]=np.linalg.eig(F)

    # Eigenvectors with norm l1 = 1
    v[:, 0] = v[:, 0] / np.linalg.norm(v[:, 0], ord=1)
    v[:, 1] = v[:, 1] / np.linalg.norm(v[:, 1], ord=1)
    v[:, 2] = v[:, 2] / np.linalg.norm(v[:, 2], ord=1)
    v[:, 3] = v[:, 3] / np.linalg.norm(v[:, 3], ord=1)

    v = permutation(v)

    #Fixing 'Error on estimated eigenvectors'
    et=np.linalg.inv(v)
    ones=np.ones((4,1))
    et=np.matmul(et,ones)
    v[:, 0]=et[0]*v[:, 0]
    v[:, 1]=et[1]*v[:, 1]
    v[:, 2]=et[2]*v[:, 2]
    v[:, 3]=et[3]*v[:, 3]

    P_a_r=v

    return P_a_r
