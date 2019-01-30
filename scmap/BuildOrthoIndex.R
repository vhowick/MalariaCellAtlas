###Making an ortholog reference index

#library(Seurat)
library(scater)
library(scmap)


###Read in SCE object and ortholog table
pb_filtered_sce <- readRDS("/Users/vh3/Documents/MCA/FieldDATA/pb_filtered_sce_20181031.rds")
newot <- read.csv("/Users/vh3/Documents/MCA/SCmap2/ortho2.csv")

#merge gene info (rowData) with ortholog table, reorder genes back to original rowData, and replace rowData with the merged df
rd <- rowData(pb_filtered_sce)
newot$Gene <- as.character(newot$Gene)
mrg <- merge(rd, newot, by.x="gene", by.y="Gene", all.x=TRUE, all.y=FALSE)

genenames <- rd$gene
mrg <- mrg[match(genenames, mrg$gene), ]
head(mrg$gene, 20)
rowData(pb_filtered_sce) <- mrg


#Remove genes that have no orth
pb_filtered_sce_orth <- pb_filtered_sce[!is.na(rowData(pb_filtered_sce)$orth_name), ]
#Replace feature_symbol and rownames with orth_name
rowData(pb_filtered_sce_orth)$feature_symbol <- rowData(pb_filtered_sce_orth)$orth_name
rownames(pb_filtered_sce_orth) <- rowData(pb_filtered_sce_orth)$orth_name

#prep the SCE, if was originally a Suerat object need the dfs to be regular matrices
pb_filtered_sce_orth <- pb_filtered_sce_orth[, colData(pb_filtered_sce_orth)$absclust3 != "8"]
sce <- pb_filtered_sce_orth
pca <- plotPCA(sce)
pcs <- pca$data
table(rownames(pcs)==colnames(sce))
colData(sce) <- cbind(colData(sce), pcs)
rowData(sce)$feature_symbol <- rowData(sce)$gene
counts(sce) <- as.matrix(counts(sce))
logcounts(sce) <- as.matrix(logcounts(sce))

#build scmap-cell reference index, save this rds
sce <- selectFeatures(sce, suppress_plot = FALSE, n_features = 500)
table(rowData(sce)$scmap_features)

set.seed(1)
sce <- indexCell(sce)
names(metadata(sce)$scmap_cell_index)
length(metadata(sce)$scmap_cell_index$subcentroids)
dim(metadata(sce)$scmap_cell_index$subcentroids[[1]])
metadata(sce)$scmap_cell_index$subcentroids[[1]][,1:5]

saveRDS(pb_filtered_sce_orth, file="pb_filtered_sce_orthindex_20181109.rds")

