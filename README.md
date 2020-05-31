# Sample-Code
This repository includes a number of scripts that I have written in python during my MSc course Algorithms in Bioinformatics. These scripts each implement a common bioinformatics tool. They are not designed to be used as standalone tools but to showcase their capabilities through function calls and printed results.

When downloading and running, keep in mind that SSAHA_sequence_search.py and Kmeans_clustering.py need files in misc to run the functions in main().

# DeBruijn_assembler.py
This script implements a DeBruijn assembler, an algorithm that splits reads into k-mers and creates a graph of consecutive k-mers. By finding a Eulerian path through the graph, an assembly can be constructed.

# Kmeans_clustering.py
This script implements a K-means clustering algorithm and evaluates between- and within-group variance.

# Needleman_wunsch_alignment.py
This script implements a dynamic programming algorithm that finds an optimal global sequence alignment.

# SSAHA_sequence_search.py
This script implements SSAHA (Sequence Search and Alignment by Hashing Algorithm) published at https://pubmed.ncbi.nlm.nih.gov/11591649/ to find closely related sequences within a hashed sequence database.
