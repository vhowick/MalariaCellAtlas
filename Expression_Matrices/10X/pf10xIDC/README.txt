pf10xIDC: 6737 P. falciparum cells that passed QC (>100 genes per cell).

important pheno columns


topcell_sp: the cell prediction cluster based on the top nearest neighbour top cell, same as topcell_ac but if cosine similarity was less than 0.3 it was 'unassigned'

topsim: cosine similarity of the top nearest neighbour

all other scmap columns are from the top 3NN method:

pb_cell: the P. berghei reference cell assigned to this P. falciparum cell based on pc coordinate space method.

pb_bulk: bulk assignment of the pb cell

pbpt: psuedotime value of assigned pb cell

PC1mean, PC2mean: mean of the PCs of the top 3NN
X_snap, Y_snap: PC coordinates of the nearest cell to PC1mean, PC2mean

absclust: CCA cluster of assigned cell

topsim1-3: cosine similarity of NN 1-3

stage_pred: the CCA cluster of assigned cell with cells with cosine similarity of less than 0.3 being unassigned. 

clock_pseudotime: pseudotime calculated from pf PCA using clock method
PC1 and PC2 PCA coordinates of Pf used to calculate clock pseudotime and shown in fig 3

bulk: highest cor with Lopez et al

