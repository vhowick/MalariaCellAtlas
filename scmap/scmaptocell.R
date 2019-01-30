setwd("/Users/vh3/Documents/MCA/SCmap")
library(scmap)
library(SingleCellExperiment)
library(plyr)
library(ggplot2)
library(viridis)
library(scater)

#Read in ortho index and query data sets. sce can be built using BuildOrthoIndex.R, follow the same first half of script to get an ortho SCE for your query data set. Also read in PCA coordinates if not already in your colData
sce <- readRDS("pb_filtered_sce_orthindex_20181109.rds")
pfobj <- readRDS("/Users/vh3/Documents/MCA/10x/24612/qc3d7100_ORTHO_SCEobj_20190103.rds")
pfpc <- read.csv("/Users/vh3/Documents/MCA/10x/24612/tmmfilt_PCS_qc3d7100_20190103.csv", row.names = 1)
table(rownames(pfpc)==colnames(pfobj))
colData(pfobj) <- cbind(colData(pfobj), pfpc)


#Project query data set onto cell index
scmapCell_results <- scmapCell(
  pfobj, 
  list(
    yan = metadata(sce)$scmap_cell_index
  )
)

##Look into the results
# For each dataset there are two matricies. cells matrix contains the top 10 (scmap default) cell IDs of the cells of the reference dataset that a given cell of the projection dataset is closest to:
#   
#   Give assignments in two ways:
#   1. Take the top cell assignment abs clust, if cosine similarity is less than 0.4 (or adjust if needed) mark as unassigned
# 2. For the top 3 nearest neighbors, get a mean of the PCA coordinates and snap to the nearest cell of those coordinates. If any of the top three cells are sim below 0.4 then mark as unassigned.


##Top cell assignment method
scmapCell_results$yan$cells[, 1:3]
getcells <- scmapCell_results$yan$cells[1, ]
cdsce <- colData(sce)[getcells, ]
topsim <- scmapCell_results$yan$similarities[1, ]

pfobj$topcell <- cdsce$sample_id
pfobj$topcell_ac <- cdsce$absclust
pfobj$indexPC1 <- cdsce$PC1
pfobj$indexPC2 <- cdsce$PC2
pfobj$pbpt <- cdsce$pseudotime
pfobj$pbbulk <- cdsce$bulk
pfobj$topcell_sp <- pfobj$topcell_ac
pfobj$topsim <- topsim
pfobj$topcell_sp[pfobj$topsim < 0.4] <- "unassigned"
table(pfobj$topcell_sp)


#### TOP 3NN method

#This function makes a list of the PC means for each cell and then do.call below rbinds them into a dataframe called big_data


datalist = list()

for (i in colnames(scmapCell_results$yan$cells)) {
  
  getcellstest <- scmapCell_results$yan$cells[1:3, i]
  cdscetest <- colData(sce)[getcellstest, ]
  PC1mean <- mean(cdscetest$PC1)
  PC2mean <- mean(cdscetest$PC2)
  # ... make some data
  dat <- data.frame(i, PC1mean, PC2mean)
  dat$i <- i  # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
# or big_data <- dplyr::bind_rows(datalist)
# or big_data <- data.table::rbindlist(datalist)

test <- big_data[1, ]

df <- data.frame(X=colData(sce)$PC1, Y=colData(sce)$PC2, row.names = rownames(colData(sce)))

#the snap function snaps to the nearest cell in PC coordiantes
snap <- function(df, test){
  require(Biobase)
  d <- matchpt(as.matrix(df),
               as.matrix(data.frame(X=test$PC1mean,Y=test$PC2mean)))
  
  min_row <- rownames(d[d$distance==min(d$distance),])
  
  test$X_snap <- unique(df[min_row,"X"])
  test$Y_snap <- unique(df[min_row,"Y"])
  test$pb_cell <- min_row
  
  test
}

#this loops through each cell and in big_data and runs the snap function
datalist2 = list()
colnames(big_data) <- c("sample_id", "PC1mean", "PC2mean")
for (i in rownames(big_data)) {
  test <- big_data[i, ]
  coord <- snap(df, test)
  coord$i <- i
  datalist2[[i]] <- coord
}
big_data2 = do.call(rbind, datalist2)


table(rownames(big_data2)==rownames(colData(pfobj)))

allpbcd <- colData(sce)
allpbcd <- as.data.frame(allpbcd)
pbabsclust <- allpbcd[, c("absclust", "pseudotime"), drop=FALSE]
pbabsclust$pb_sample_id <- rownames(pbabsclust)

# Now merge the pc cell asignments with their abs clust and get in the right order
big_data3 <- merge(big_data2, pbabsclust, by.x = "pb_cell", by.y = "pb_sample_id", all.x=TRUE, all.y=FALSE)
big_data4 <- big_data3[match(rownames(big_data2), big_data3$sample_id), ]


colors <- c("6"="#78C679",
            "2"="#D1EC9F",
            "0"="#FEB24C",
            "1"="#F4CF63",
            "3"="#FEEEAA",
            "4"="#85B1D3",
            "7"="#9ecae1",
            "5"="#C9E8F1",
            "M"= "#B7B7D8",
            "F"="#9C96C6",
            "unassigned"="black")

ggplot(big_data4, aes(PC1mean, PC2mean)) + geom_point(aes(colour=factor(big_data4$absclust))) + scale_color_manual(values = colors) 

ggplot(big_data4, aes(X_snap, Y_snap)) + geom_point(aes(colour=factor(big_data4$absclust))) + scale_color_manual(values = colors) 

##add info to SCE and save colData, to be an assigned cell all 3NN must have a cos sim >0.4
scmapCell_results$yan$similarities[, 1:3]
topsim1 <- scmapCell_results$yan$similarities[1, ]
topsim2 <- scmapCell_results$yan$similarities[2, ]
topsim3 <- scmapCell_results$yan$similarities[3, ]

table(big_data4$sample_id==rownames(colData(pfobj)))
#pfobj$pb_cell <- big_data4$pb_cell
#pfobj$PC1mean <- big_data4$PC1mean
bd4 <- big_data4[, c("pb_cell", "sample_id", "PC1mean", "PC2mean", "X_snap", "Y_snap", "absclust", "pseudotime" )]

colData(pfobj) <- cbind(colData(pfobj), bd4)

pfobj$topsim1 <- topsim1
pfobj$topsim2 <- topsim2
pfobj$topsim3 <- topsim3
pfobj$stage_pred <- pfobj$absclust
pfobj$stage_pred[pfobj$topsim1 < 0.4 | pfobj$topsim2 < 0.4 | pfobj$topsim3 < 0.4] <- "unassigned"
table(pfobj$stage_pred)

#write.csv(bd4, "pfcellassignmentswithmean3nn_20181029.csv")
write.csv(colData(pfobj), "pf3d7100scmapclusts2methodindexn100_20190107.csv")

