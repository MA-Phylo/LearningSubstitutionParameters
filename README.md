---------------------------------------------HOW TO GENERATE ALIGNMENTS?---------------------------------------------

Option 1: Generate t FASTA files of alignments of length L given a tree in a newick format:

For example, t = 5 and L = 1000;  
	`python3 GenGM.py tree_4L.txt 5 1000` 

Option 2: Generate FASTA files of alignments of given lengths L1...Ld. **Note that in this case, lenghts sequences are preceeded by an L**

For example, L1 = 500, L2 = 1000 and L3 = 10000;  
	`python3 GenGM.py tree_4L.txt L500 L1000 L10000` 
