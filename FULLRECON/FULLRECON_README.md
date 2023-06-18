## FULLRECON - Substitution parameters estimation on phylogenetic trees

It includes all the scripts needed to the reproduce the pipeline described in the report. Given a phylogenetic tree with a pre-fixed number of substitutions at the edges (branch lengths), 
generates its transition matrices and the corresponding alignment of the leaf DNA sequences. Then, FULLRECON is ran using as input these sequences to compare the generated matrices from which they come with the estimated ones.

Command run example:  `python3 evaluate.py tree_4L.txt 1000`. tree_4L.txt is the studied phylogenetic tree (in a Newick format) and 1000 the alignment length.
