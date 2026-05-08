##### Packages #####
library(Seurat)
library(tidyverse)

set.seed(666)
#
##### Loading In #####
sc2week<- Read10X_h5("../../results/02_cellrangercount/2wkold/outs/filtered_feature_bc_matrix.h5")
sc6week<- Read10X_h5("../../results/02_cellrangercount/6wkold/outs/filtered_feature_bc_matrix.h5")

sdata.sc2week<- CreateSeuratObject(sc2week)
sdata.sc6week<- CreateSeuratObject(sc6week)

alldata<- merge(sdata.sc2week, c(sdata.sc6week), add.cell.ids = c("2Week", "6Week"))

rm(sc2week, sc6week, sdata.sc2week, sdata.sc6week)

alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "Mitochondria")
alldata$mitoPercent<- PercentageFeatureSet(alldata, pattern = '^mt-')
# neither of the above show any mitochondrial %, cannot filter out what is written in the paper
# no mitochondrial percent in this data 
data.filt<- subset(alldata, subset = nCount_RNA > 250 & 
                     nFeature_RNA < 2500 & 
                     mitoPercent < 5)



data.filt
# An object of class Seurat 
# 32285 features across 9775 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 2 layers present: counts.1, counts.2

rm(alldata)

data.filt<- JoinLayers(data.filt)
#
##### QC Metrics #####
VlnPlot(data.filt, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), pt.size = 0.25, ncol = 3)
# FeatureScatter(data.filt, feature1 = "nCount_RNA", feature2 = "mitoPercent")
FeatureScatter(data.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dim(data.filt)
# 32285  9775
ncol(data.filt)
# 9775

data.filt <- NormalizeData(data.filt, normalization.method = "LogNormalize", scale.factor = 10000)

data.filt$sample<- rownames(data.filt@meta.data)
data.filt@meta.data<- separate(data.filt@meta.data, col = 'sample', into = c('condition', 'barcode'), sep = '_')
data.filt = FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = 2000, assay = "RNA")

# show top 10 variable genes (features)
top10vf<- head(VariableFeatures(data.filt), 10)
vfplot<- VariableFeaturePlot(data.filt)
LabelPoints(plot = vfplot, points = top10vf, repel = TRUE)

data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "mitoPercent"),verbose = TRUE)

data.filt <- RunPCA(data.filt, features = VariableFeatures(object = data.filt))
print(data.filt[["pca"]], dims = 1:5, nfeatures = 5)

# check the elbow plot
ElbowPlot(data.filt)

data.filt = RunUMAP(data.filt, dims = 1:9, verbose = T)

data.filt<- JoinLayers(data.filt)


rm(vfplot, top10vf)
#
##### DOUBLET FINDER #####
library(DoubletFinder)

# if error in getGlobalsAndPackages occurs:
options(future.globals.maxSize = 8000 * 1024^2)

sweep.res<- paramSweep(data.filt, PCs = 1:9, sct = TRUE)
sweep.stats<- summarizeSweep(sweep.res, GT = FALSE)
bcmvn<- find.pK(sweep.stats)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()
pK<- bcmvn%>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK<- as.numeric(as.character(pK[[1]]))
# assume 5% of cells are heterotypic doublets
nExp<- round(ncol(data.filt) * 0.05) 
data.filt<- doubletFinder(data.filt, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:9)
table(data.filt@meta.data$DF.classifications_0.25_0.07_489)
# Doublet Singlet 
#   489    9286
DF.name = colnames(data.filt@meta.data)[grepl("DF.classifications_0.25_0.07_489", colnames(data.filt@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "condition") + NoAxes(), DimPlot(data.filt, group.by = DF.name) + NoAxes())
VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
dim(data.filt)
# 32285  9286

rm(sweep.res, sweep.stats, nExp, pK, bcmvn, DF.name)
data.filt@meta.data$DF.classifications_0.25_0.07_489 <- NULL
data.filt@meta.data$pANN_0.25_0.07_489 <- NULL
#
##### Data Integration #####

data.filt.split<- SplitObject(data.filt, split.by = "condition")

for (i in 1:length(data.filt.split)) {
  data.filt.split[[i]]<- NormalizeData(data.filt.split[[i]])
  data.filt.split[[i]]<- FindVariableFeatures(data.filt.split[[i]], selection.method = "vst", nfeatures = 2000)
}

features<- SelectIntegrationFeatures(data.filt.split)
for (i in seq_along(along.with = data.filt.split)) {
  data.filt.split[[i]]<- ScaleData(data.filt.split[[i]], features = features)
  data.filt.split[[i]]<- RunPCA(data.filt.split[[i]], features = features)
}

anchors<- FindIntegrationAnchors(object.list = data.filt.split)
data.filt.integrated<- IntegrateData(anchorset = anchors)
data.filt.integrated<- ScaleData(data.filt.integrated)
data.filt.integrated<- RunPCA(data.filt.integrated)
data.filt.integrated<- RunUMAP(data.filt.integrated, dims = 1:11, reduction.name = "UMAP")
data.filt.integrated<- FindNeighbors(data.filt.integrated, dims = 1:11)
data.filt.integrated<- FindClusters(data.filt.integrated)
# Number of nodes: 9286
# Number of edges: 305622
# Maximum modularity in 10 random starts: 0.8677
# Number of communities: 14

DimPlot(data.filt.integrated, label = TRUE)
DimPlot(data.filt.integrated, split.by = "condition", label = TRUE)

rm(data.filt, data.filt.split, anchors, features, i)

DefaultAssay(data.filt.integrated)<- "RNA"
#
##### Typing Dot Plot #####
# from figure 1 to cell type the clusters
markers.to.plot.typing<- c("Col1a1", "Fmod", "Col3a1", "Gsn", "Sox9", "AY036118", "Gm42418", "Ly6a", "Tppp3", "Pdgfra", "Cd34", "Dpt", "Lyz2", "Pecam1", "Hba-a1", "Myl9", "Mki67", "Stmn1", "C1qa", "Mkx", "Scx", "Tnmd", "Postn", "Mbp", "Mmrn1", "Ptprc")

DotPlot(data.filt.integrated, features = c(markers.to.plot.typing), 
        dot.scale = 7,
        scale = FALSE,
        group.by = "seurat_clusters") + RotatedAxis()

rm(markers.to.plot.typing)

# two major genes of interest
FeaturePlot(data.filt.integrated, c("Cd55", "Cd248"))

##### Rename Cell Idents #####
# dot plot looks completely different than what is shown in the paper
# whatever I can't match will be called 'tendon'
# SP1/2: tendon stem/progenitor cell
# RBC: red blood cell
# LC: lymphocyte
data.filt.integrated<- RenameIdents(data.filt.integrated, c(`0` = "0_Tendon", `1` = "1_Tendon", `2` = "2_SP1", `3` = "3_RBC", `4` = "4_SP2", `5` = "5_Tendon", `6` = "6_Tendon", `7` = "7_Tendon", `8` = "8_LC", `9` = "9_Tendon", `10` = "10_Tendon", `11` = "11_Tendon", `12` = "12_Tendon", `13` = "13_Tendon", `14` = "14_Tendon"))

data.filt.integrated$CellType <- Idents(data.filt.integrated)

# UMAP with annotated clusters
DimPlot(data.filt.integrated, label = TRUE)

DimPlot(data.filt.integrated, split.by = "condition", label = TRUE)

# remake the dot plot with annotated clusters 
markers.to.plot.typing<- c("Col1a1", "Fmod", "Col3a1", "Gsn", "Sox9", "AY036118", "Gm42418", "Ly6a", "Tppp3", "Pdgfra", "Cd34", "Dpt", "Lyz2", "Pecam1", "Hba-a1", "Myl9", "Mki67", "Stmn1", "C1qa", "Mkx", "Scx", "Tnmd", "Postn", "Mbp", "Mmrn1", "Ptprc")

DotPlot(data.filt.integrated, features = c(markers.to.plot.typing), 
        dot.scale = 7,
        scale = FALSE,
        group.by = "CellType") + RotatedAxis()

rm(markers.to.plot.typing)


##### Barchart #####
library(plyr)

# this chart visualized the differences in clusters better than the original and it shows how different the cell recovery is between the groups 

barchart_10<-table(Idents(data.filt.integrated), data.filt.integrated$condition)
barchart_10<-as.data.frame(barchart_10)
barchart_10$Cluster<-as.character(barchart_10$Var1)
print(barchart_10)
barchart_10_sorted<-arrange(barchart_10, Var2, Var1)
head(barchart_10_sorted)
barchart_10_cumsum<-ddply(barchart_10_sorted, "Var2", transform, label_ypos=cumsum(Freq)-0.5*Freq)
head(barchart_10_cumsum)
mt10_barchart<-ggplot(barchart_10_cumsum[which(barchart_10$Freq>0),], aes(x=Var2, y=Freq, fill=Var1))+
  geom_bar(position= 'stack', stat = 'identity')+
  geom_text(aes(label=Freq), size=1.5, position = position_stack(vjust=0.5, reverse = FALSE))+
  labs(x=NULL, y=NULL, title = 'Number of Cells for Each Cluster of Each Condition')+
  theme_bw()
print(mt10_barchart)

rm(barchart_10, barchart_10_cumsum, barchart_10_sorted, mt10_barchart)

##### Violin Plots of Genes #####
VlnPlot(data.filt.integrated, features = c("Pdgfra", "Cd34","Cd55", "Cd248"), idents = c("2_SP1", "4_SP2"), split.by = "condition",ncol = 2)


##### Trajectory and Pseudotime #####
# to see changes between groups, the object must be split
library(monocle3)
library(SeuratWrappers)

splitdata<- SplitObject(data.filt.integrated, split.by = "condition")

sc2week<- splitdata$`2Week`
sc6week<- splitdata$`6Week`
sc2week<- JoinLayers(sc2week, assay = 'RNA')
sc6week<- JoinLayers(sc6week, assay = 'RNA')
sc2week$CellType <- Idents(sc2week)
sc6week$CellType <- Idents(sc6week)
DefaultAssay(sc2week)<- "RNA"
DefaultAssay(sc6week)<- "RNA"

cds2wk<- as.cell_data_set(sc2week)
cds2wk<- cluster_cells(cds2wk, cluster_method = 'louvain')
cds2wk<- learn_graph(cds2wk)
plot_cells(cds2wk, 
           color_cells_by = "ident", 
           label_groups_by_cluster = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE)
# picking root node cluster 2
cds2wk<- order_cells(cds2wk)
plot_cells(cds2wk, 
           color_cells_by = "pseudotime", 
           group_cells_by = "cluster", 
           label_cell_groups = FALSE, 
           label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           trajectory_graph_color = "black")

# make a box plot for visualization
cds2wk$monocle3_pseudotime<- pseudotime(cds2wk)
data.pseudo2wk<- as.data.frame(colData(cds2wk))
ggplot(data.pseudo2wk, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime), fill = ident)) + geom_boxplot()


cds6wk<- as.cell_data_set(sc6week)
cds6wk<- cluster_cells(cds6wk, cluster_method = 'louvain')
cds6wk<- learn_graph(cds6wk)
plot_cells(cds6wk, 
           color_cells_by = "ident", 
           label_groups_by_cluster = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE)
# picking root node cluster 2
cds6wk<- order_cells(cds6wk)
plot_cells(cds6wk, 
           color_cells_by = "pseudotime", 
           group_cells_by = "cluster", 
           label_cell_groups = FALSE, 
           label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, 
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           trajectory_graph_color = "black")

# make a box plot for visualization
cds6wk$monocle3_pseudotime<- pseudotime(cds6wk)
data.pseudo6wk<- as.data.frame(colData(cds6wk))
ggplot(data.pseudo6wk, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime), fill = ident)) + geom_boxplot()

rm(cds2wk, data.pseudo2wk, cds6wk, data.pseudo6wk, splitdata, sc2week, sc6week)
#
