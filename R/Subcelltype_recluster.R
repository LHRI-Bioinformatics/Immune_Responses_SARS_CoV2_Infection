
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggbreak)

library(SingleR)
library(scran)
library(scRNAseq)
library(scuttle)
library(celldex)
library(singleCellNet)
source("./Function.R")
source("./doublet_function")

Subecells="sub_cellselial"

sub_cells_cells=subset(all_Samples_f.list_NFSR.inted_annoted, subset =major_cell_type== c(Subecells))

DefaultAssay(sub_cells_cells) <- "RNA"

###As some Samples have few sub_cellselial. So we would filter our samples with cell number at least 20
Cellnamubers=as.data.frame(table(sub_cells_cells$orig.ident,sub_cells_cells$state))
Cellnamubers1=Cellnamubers[Cellnamubers$Freq>20,]
sub_cells_cells1=subset(sub_cells_cells, subset =orig.ident %in% Cellnamubers1$Var1)
sub_cells_cells2=doublet(sub_cells_cells1)

##as the doublet cells will affected the sub_cellselial annotation in downsream study, we will remove doulbet if there is the clusters are all of doublet.
sub_cells_cells3=subset(sub_cells_cells2,subset = DoubletScore_n == "singlet")
OVER20CELLS=names(table(sub_cells_cells3$orig.ident)[table(sub_cells_cells3$orig.ident)>20])
sub_cells_cells4=subset(sub_cells_cells3,subset = orig.ident %in% OVER20CELLS)


Subsesub_cells.list <- SplitObject(sub_cells_cells4, split.by = "orig.ident")
for (i in 1:length(Subsesub_cells.list)) {
  Subsesub_cells.list[[i]] <- NormalizeData(Subsesub_cells.list[[i]], verbose = FALSE)
  Subsesub_cells.list[[i]] <- FindVariableFeatures(Subsesub_cells.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}


features <- SelectIntegrationFeatures(object.list = Subsesub_cells.list)
Subsesub_cells.list_NFSR <- lapply(X = Subsesub_cells.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs=15,approx = FALSE)
})



###so chose 6,8 as reference integation achors
###adjust parameter based on the subg cell types, Epithelial and Macrophage are different
anchors <- FindIntegrationAnchors(object.list = Subsesub_cells.list_NFSR, reference = c(15, 16), reduction = "rpca",dims = 1:10,k.anchor = 5,
                                  k.filter = k.filter,
                                  k.score = k.score,
                                  max.features = max.features,
                                  nn.method = "annoy",
                                  n.trees = n.trees)

##small cell number can not do integratedDana, reduce K.weight to sovle this problem.Mar2nd2023.
all_Samples_f.list_NFSR.inted <- IntegrateData(anchorset = anchors, dims = 1:10,k.weight=18)


DefaultAssay(all_Samples_f.list_NFSR.inted)


###first generate data and scaledata in RNA assay
DefaultAssay(all_Samples_f.list_NFSR.inted) <- "RNA"
all_Samples_f.list_NFSR.inted[['percent.mito']] <- PercentageFeatureSet(all_Samples_f.list_NFSR.inted, pattern = "^MT-")
sub_cells.Integrated <- NormalizeData(object = all_Samples_f.list_NFSR.inted, normalization.method = "LogNormalize", scale.factor = 1e4)
sub_cells.Integrated <- FindVariableFeatures(object = sub_cells.Integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
sub_cells.Integrated <- ScaleData(sub_cells.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))


##change to integrated assay
DefaultAssay(sub_cells.Integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
sub_cells.Integrated <- ScaleData(sub_cells.Integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
sub_cells.Integrated <- RunPCA(sub_cells.Integrated, verbose = FALSE)


ElbowPlot(object = sub_cells.Integrated,ndims = 10)

# 
 sub_cells.Integrated <- FindNeighbors(object = sub_cells.Integrated, dims = 1:10)

sub_cells.Integrated <- FindClusters(object = sub_cells.Integrated, resolution = 0.4)

##tsne and umap
sub_cells.Integrated <- RunTSNE(object = sub_cells.Integrated, dims = 1:10)
sub_cells.Integrated <- RunUMAP(sub_cells.Integrated, reduction = "pca", dims = 1:10)





###annotation 
annotation_sub_cells=read.csv("sub_cells_cells_annotation_May1st2023.csv")
uniq_annoted= unique.data.frame(annotation_sub_cells[,c("cluster","annoted")])
uniq_annoted_cl=uniq_annoted$annoted
names(uniq_annoted_cl) <- levels(sub_cells.Integrated)
Idents()
sub_cells.Integrated_annoted <- RenameIdents(sub_cells.Integrated, uniq_annoted_cl)
sub_cells.Integrated_annoted$annotation_summary=Idents(sub_cells.Integrated_annoted)




