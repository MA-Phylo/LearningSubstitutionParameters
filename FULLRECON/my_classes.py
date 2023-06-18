# Transition Matrix class. Defined by an ancestor node, a descendant node
# and the corresponding matrix.
class MM:
    def __init__(self, source, target, matrix):
        self.source = source
        self.target = target
        self.matrix = matrix

# SepRecon input class. Defined by a t value and the corresponding leaves a, a'
# and separator edge (w, w')
class SepReconInput:
    def __init__(self, t, a, a_prime, w, w_prime):
        self.t = t
        self.a = a
        self.a_prime = a_prime
        self.w = w
        self.w_prime = w_prime
