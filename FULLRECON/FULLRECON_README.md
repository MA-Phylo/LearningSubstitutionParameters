## FULLRECON - Substitution parameters estimation on phylogenetic trees

It includes all the scripts needed to the reproduce the pipeline described in the report. Given a phylogenetic tree with a pre-fixed number of substitutions at the edges (branch lengths), 
generates its transition matrices and an alignment of the leaf DNA sequences. Then, FULLRECON is ran using as input these sequences to estimate the generated matrices from which they come.

Command run example: if we want to study the phylogenetic tree in tree_4L.txt (in a Newick format) and use an alignment of length 1000, we should run  `python3 evaluate.py tree_4L.txt 1000`

The output displayed includes, for each branch, the simulated and estimated transition matrix and the branch length estimation. The likelihood of the phylogenetic tree given the estimated parameters it is also provided.
