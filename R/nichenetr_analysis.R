library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)


organism = "human"


#lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
lr_network=readRDS("lr_network_human_21122021.rds")
ligand_target_matrix=readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

#seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
###read seuratOBJ
All_inte_Anned=all_Samples_f.list_NFSR.inted_annoted_sub
All_inte_Anned$state=str_replace_all(All_inte_Anned$state,"unvac","Unvaxed")
All_inte_Anned$state=str_replace_all(All_inte_Anned$state,"vacced","Vaxed")
All_inte_Anned$celltype=All_inte_Anned$major_cell_type


DefaultAssay(All_inte_Anned)="RNA"
nichenet_output1 = nichenet_seuratobj_aggregate(
  seurat_obj = All_inte_Anned, 
  receiver = c("mDC","pDC","Macrophage" ) ,
  condition_colname = "state", condition_oi ="Vaxed_R"  , condition_reference = "Unvaxed_R"  , 
  sender = c( "CD4","CD8","B_cell","NK"), 
  ligand_target_matrix = ligand_target_matrix[,c("CLU","ANXA1","CRTAM","COL19A1","TNFSF4")],
  lr_network = lr_network[lr_network$from %in% c("CLU","ANXA1","CRTAM","COL19A1","TNFSF4"),],
  weighted_networks = weighted_networks,filter_top_ligands=T,
  top_n_ligands=500,top_n_targets=1000,assay_oi="RNA",cutoff_visualization=0.1,lfc_cutoff=0.1,expression_pct=0.03
  
)
nichenet_output1$top_ligands [nichenet_output1$top_ligands %in% c("ANXA1","CRTAM") ]
nichenet_output1$top_receptors [nichenet_output1$top_receptors %in% c("FPR1","CADM1") ]

saveRDS(nichenet_output1,"TNK_to_DC_Macrophage__Vaxed_R_vs_UnvaxedR_targets_Jan23th2024.rds")

