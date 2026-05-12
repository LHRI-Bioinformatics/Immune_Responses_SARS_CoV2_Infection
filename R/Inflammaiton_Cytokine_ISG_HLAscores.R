
library(Seurat)
library(stringr)
library(ggplot2)
#####Add module score 
##interferon

All_inte_Anned=readRDS("all_samples_annotated.rds")


Cytokines=read_excel("Cytokine_Inflammation_Score/Cytokines_genes.xlsx", sheet = 2)
Cytokines=as.matrix(as.data.frame(Cytokines))
Cytokines_name=sapply(Cytokines, as.character) 
Cytokines_name=unique(Cytokines_name[!is.na(Cytokines_name)])


Cytokines_name1=list(c(Cytokines_name,'CXCL8','CCL22','CCL17','CCL19','CXCL12','IL33','CXCL12'))


All_inte_Anned1_mean=as.data.frame(AverageExpression(All_inte_Anned,assays = "RNA", features =unlist(Cytokines_name1) ))

##select top30 high expresse cytokines genes
Inflam_cyto_h=rownames(All_inte_Anned1_mean)[order(rowMeans(All_inte_Anned1_mean),decreasing = T)][1:30]

DefaultAssay(All_inte_Anned)="RNA"
All_inte_Anned2 <- AddModuleScore(
  object = All_inte_Anned,
  features = list(Inflam_cyto_h),
  ctrl = 5,
  name = 'Cytokine_score'
)


#### Interferon genes 
inferon_gene=read.table("interferon_genes_2.txt",sep="\t",row.names=NULL)

inferon_genes=unique(c(inferon_gene[,2],colnames(inferon_gene)[2]))



##check the average expressin of the inferon genes and select the top30 genes

All_inte_Anned1_mean=as.data.frame(AverageExpression(All_inte_Anned,assays = "RNA", features =(inferon_genes) ))
inferon_genes_f=rownames(All_inte_Anned1_mean)[order(rowMeans(All_inte_Anned1_mean),decreasing = T)][1:30]

All_inte_Anned_cyto <- AddModuleScore(
  object = All_inte_Anned2,
  features = list(inferon_genes_f),
  ctrl = 5,
  name = 'Inferon_score'
)

#### add inflammtion
Inflam=read_excel("Cytokine_Inflammation_Score/Cytokines_genes.xlsx", sheet = 1)
Inflam=as.matrix(as.data.frame(Inflam))
Inflam_na=sapply(Inflam, as.character) 
Inflam_na1=unique(Inflam_na[!is.na(Inflam_na)])

Inflam_na1=list(c(Inflam_na1))


All_inte_Anned1_mean=AverageExpression(All_inte_Anned_cyto,assays = "RNA", features =unlist(Inflam_na1) )
All_inte_Anned1_mean=(as.data.frame(All_inte_Anned1_mean))

Inflam_na1_h=rownames(All_inte_Anned1_mean)[order(rowMeans(All_inte_Anned1_mean),decreasing = T)][1:30]





###Use high expression gene with mean high than 0.5
All_inte_Anned3 <- AddModuleScore(
  object = All_inte_Anned_cyto,
  features = list(Inflam_na1_h),
  ctrl = 5,
  name = 'Inflammation_score_high'
)




####HLA_I gene score


Salina_HLAI1=rownames(All_inte_Anned)[grep("^HLA",rownames(All_inte_Anned))]
Salina_HLAI2=c("HLA-F","HLA-G","HLA-A","HLA-E","HLA-C","HLA-B")

All_inte_Anned_HLAI <- AddModuleScore(
  object = All_inte_Anned3,
  features = list(Salina_HLAI2),
  ctrl = 5,
  name = 'HLAI_score_all'
)

Salina_HLAII1=rownames(All_inte_Anned)[grep("^HLA",rownames(All_inte_Anned))]

Salina_HLAII2=Salina_HLAII1[!Salina_HLAII1 %in% c("HLA-F","HLA-G","HLA-A","HLA-E","HLA-C","HLA-B")]

All_inte_Anned_HLAII <- AddModuleScore(
  object = All_inte_Anned_HLAI,
  features = list(Salina_HLAII2),
  ctrl = 5,
  name = 'HLAII_score_all'
)

saveRDS(All_inte_Anned_HLAII,"all_Samples_Cytokine_Inflammat_ISG_HLAscores.rds")


