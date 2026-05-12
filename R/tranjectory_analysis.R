


#####MONCLE 
source("../../Function.R")
library(Seurat)
library(monocle3)

library(ggplot2)
library(stringr)
library(DDRTree)
library(reticulate)
library(devtools)
library(VGAM)
library(Matrix.utils)
library(DelayedMatrixStats)

#### keep data same to using  Subcells_Trajectory_analysis_add_clonotypes_05_18_2023.rds data

Subcells_Inted_Annoted=Subcells.Integrated_annoted


DefaultAssay(Subcells_Inted_Annoted)="RNA"

data <- as(as.matrix(Subcells_Inted_Annoted@assays$RNA@data), 'sparseMatrix')
Subcellss_meta=Subcells_Inted_Annoted@meta.data

Subcellss_gene_meta=data.frame(id=rownames(data),gene_short_name=rownames(data))
rownames(Subcellss_gene_meta)=rownames(data)


Subcellss_inte_cds <- new_cell_data_set(data,
                                     cell_metadata = Subcellss_meta,
                                     gene_metadata = Subcellss_gene_meta)
##
Subcellss_inte_cds <- preprocess_cds(Subcellss_inte_cds, num_dim = 50,norm_method = "none")
plot_pc_variance_explained(Subcellss_inte_cds)
Subcellss_inte_cds_a <- align_cds(Subcellss_inte_cds,alignment_k = 50,alignment_group = "orig.ident")

Subcellss_inte_cds_ad= reduce_dimension(Subcellss_inte_cds_a)


Subcellss_inte_cds_ac <- cluster_cells(Subcellss_inte_cds_ad,k=100)

Subcellss_inte_cds_ac@clusters$UMAP$partitions[Subcellss_inte_cds_ac@clusters$UMAP$partitions != "1"] <- "1"

Subcellss_inte_cds_acl <- learn_graph(Subcellss_inte_cds_ac)



Subcellss_inte_cds_acl <- order_cells(Subcellss_inte_cds_acl, root_pr_nodes=get_earliest_principal_node_Subcellss(Subcellss_inte_cds_acl))



Subcells_test1=Subcellss_inte_cds_acl
#### pseudotime summary
### using ggboxplot Jan9th2024 



