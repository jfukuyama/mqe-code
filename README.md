# Code for paper on phylogenetic distances

This repository contains the code to reproduce the analyses in the paper "Uncovering the features constructed by a class of phylogenetically-informed distances."
- The file `eigen-fns.R` contains some helper functions related to the theoretical eigenvectors of the balanced binary tree.
- The file `evec-plots.R` contains the code that creates the plots of the theoretical eigenvectors/eigenvalues of the systematic trees discussed in the paper.
- The file `tsimane.R` contains the code that performs the analysis on the Tsimane microbiota data. Before running this script, you need to download the data from https://purl.stanford.edu/tv993xn7633 and place the file `Sprockett_Tsimane_ps.rds` in the same directory as `tsimane.R`.
