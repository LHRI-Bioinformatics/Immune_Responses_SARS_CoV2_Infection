
library(Seurat)

library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)

all_Samples_f.list_NFSR.inted_annoted=readRDS(Seuratdataset)

all_Samples_f.list_NFSR.inted_annoted$celltype.stim<-paste((all_Samples_f.list_NFSR.inted_annoted$major_cell_type),
                                                           all_Samples_f.list_NFSR.inted_annoted$state,sep = '_')

Idents(all_Samples_f.list_NFSR.inted_annoted)<-'celltype.stim'


All_idents<-as.character(unique(Idents(all_Samples_f.list_NFSR.inted_annoted)))

library(DESeq2)
library(stringr)
library(Seurat)
library(limma)
obj="All_sample"


all_Samples_f.list_NFSR.inted_annoted$state=str_replace_all(all_Samples_f.list_NFSR.inted_annoted$state,"unvac","Unvaxed")
all_Samples_f.list_NFSR.inted_annoted$state=str_replace_all(all_Samples_f.list_NFSR.inted_annoted$state,"vacced","Vaxed")

selected_clustered=selected_clustered
###change back to clusters
all_Samples_f.list_NFSR.inted_annoted=subset(all_Samples_f.list_NFSR.inted_annoted,subset= (seurat_clusters %in% selected_clustered) & (state %in% c( "Unvaxed_A","Unvaxed_R","Vaxed_A","Vaxed_R")))

all_Samples_f.list_NFSR.inted_annoted_ps= AggregateExpression(all_Samples_f.list_NFSR.inted_annoted1, assays = "RNA", return.seurat = T, group.by = c("state", "orig.ident", "seurat_clusters"))

all_Samples_f.list_NFSR.inted_annoted_ps$celltype.stim<-paste((all_Samples_f.list_NFSR.inted_annoted_ps$seurat_clusters),
                                                              all_Samples_f.list_NFSR.inted_annoted_ps$state,sep = '_')
                                                            


Date_counts=(all_Samples_f.list_NFSR.inted_annoted_ps@assays$RNA$counts)



library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = Date_counts,
                              colData = all_Samples_f.list_NFSR.inted_annoted_ps@meta.data,
                              design= ~ celltype.stim)

dds <- DESeq(dds)
All_idents<-as.character(unique((dds$celltype.stim)))

saveRDS(dds,"cluster_state_deseqdataset_Mar9th2026.rds")

library(DESeq2)
library(stringr)


infname=c("Unvaxed-A","Unvaxed-R",    "Unvaxed-A",  "Vaxed-A")
samplename=c("Vaxed-A",  "Vaxed-R",    "Unvaxed-R",  "Vaxed-R")
for (i in 1:4){
  ass2<-data.frame()
  for (xx in c( 0:23,25:29)){
    
    
    ident.2=All_idents[grep(paste("^",xx,"_",samplename[i],sep = ""),All_idents)]
    ident.1=All_idents[grep(paste("^",xx,"_",infname[i],sep = ""),All_idents)]
    print(ident.1)
    
      ass1 <- results(dds, 
                     contrast=c("celltype.stim",ident.1 ,ident.2),name = paste(ident.1,ident.2,sep = "_vs_")) 
      ass1$clusters<-xx
      
      ass1$symbol<-rownames(ass1)
      ass2<-rbind(ass2,as.data.frame(ass1))
      print("3")
  }
  ass2=ass2[!is.na(ass2$padj),]
  
  
  OriUAVA=read.csv(paste("cell/based/DEGs/",infname_o[i],"_vs_",samplename_o[i],"_DEGs.csv",sep = ""))
  pois_glmm_df1=na.omit(ass2)
  pois_glmm_df1=merge(pois_glmm_df1,OriUAVA,by=c( "clusters",'symbol'),all.x=T)
  
  write.csv(pois_glmm_df1,file = paste(infname[i],'_vs_',samplename[i],'_altered_markers_Gene_Pseudobulk_DESEQ2',obj,'.csv',sep = ''),row.names = F)
  
  ass3<-pois_glmm_df1[(pois_glmm_df1$padj<0.05)&((pois_glmm_df1$log2FoldChange>0.59)|(pois_glmm_df1$log2FoldChange<(-0.59))),]
  write.csv(ass3,file = paste(infname[i],'_vs_',samplename[i],'_altered_markersFC1.5padj0.05PCT25_Gene_Pseudobulk_mDESEQ2_',obj,'.csv',sep = ''),row.names = F)
  
  }

#####Voom by group and Voom with quality weight

####can we have differential gene and Gene expression on all ciliated cells and inflammatory  cells  July13th 2023
library(Seurat)
library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)
library(DESeq2)
library(stringr)
library(Seurat)
library(limma)

#####all_Samples_f.list_NFSR.inted_annoted can be all samples seurat dataset or subcell type sell dataset


### all Samples and or subcell types
obj="Sample datase All_sample/TNK/Macrophage/Epithelial"


all_Samples_f.list_NFSR.inted_annoted$state=str_replace_all(all_Samples_f.list_NFSR.inted_annoted$state,"unvac","Unvaxed")
all_Samples_f.list_NFSR.inted_annoted$state=str_replace_all(all_Samples_f.list_NFSR.inted_annoted$state,"vacced","Vaxed")

all_Samples_f.list_NFSR.inted_annoted_ps= AggregateExpression(all_Samples_f.list_NFSR.inted_annoted, assays = "RNA", return.seurat = T, group.by = c("state", "orig.ident", "major_cell_type"))

all_Samples_f.list_NFSR.inted_annoted_ps$celltype.stim<-paste((all_Samples_f.list_NFSR.inted_annoted_ps$major_cell_type),
                                                              all_Samples_f.list_NFSR.inted_annoted_ps$state,sep = '_')

Idents(all_Samples_f.list_NFSR.inted_annoted)<-'celltype.stim'


All_idents<-as.character(unique(Idents(all_Samples_f.list_NFSR.inted_annoted)))
Date_counts=(all_Samples_f.list_NFSR.inted_annoted_ps@assays$RNA$counts)

cluster_celltye=unique.data.frame(all_Samples_f.list_NFSR.inted_annoted1@meta.data[,c("seurat_clusters","annotated name")])

meata_data=all_Samples_f.list_NFSR.inted_annoted_ps@meta.data
meata_data$celltype.stim=str_replace_all(meata_data$celltype.stim,"-","_")

library(stringr)
keep_gene <- (rowSums(Date_counts > 0) > 50 & rowMeans(Date_counts) > 0.4)


All_edgeR_results=data.frame()

cell_dataset_s=Date_counts[keep_gene,colnames(Date_counts)[grepl(paste("-A_.+","_",sep = ""),colnames(Date_counts))]]
meata_data_cells=meata_data[grepl(paste("_.+_A$",sep = ""),meata_data$celltype.stim),]

pb_dge <- DGEList(
  counts = cell_dataset_s,
  samples = meata_data_cells,
  group = meata_data_cells$celltype.stim)

pb_dge <- calcNormFactors(pb_dge)





group_edgeR <- as.character(pb_dge$samples$group)
model.matrix(~0+group_edgeR) -> design

###manuuall add the comparasion groups
contr.matrix <- makeContrasts(
  UnvsV0 = group_edgeR0_Unvaxed_A - group_edgeR0_Vaxed_A,
  UnvsV1 = group_edgeR1_Unvaxed_A - group_edgeR1_Vaxed_A,
  UnvsV2 = group_edgeR2_Unvaxed_A - group_edgeR2_Vaxed_A,
  UnvsV3 = group_edgeR3_Unvaxed_A - group_edgeR3_Vaxed_A,
  UnvsV4 = group_edgeR4_Unvaxed_A - group_edgeR4_Vaxed_A,
  UnvsV5 = group_edgeR5_Unvaxed_A - group_edgeR5_Vaxed_A,
  UnvsV6 = group_edgeR6_Unvaxed_A - group_edgeR6_Vaxed_A,
  UnvsV7 = group_edgeR7_Unvaxed_A - group_edgeR7_Vaxed_A,
  UnvsV8 = group_edgeR8_Unvaxed_A - group_edgeR8_Vaxed_A,
  UnvsV9 = group_edgeR9_Unvaxed_A - group_edgeR9_Vaxed_A,
  UnvsV10 = group_edgeR10_Unvaxed_A - group_edgeR10_Vaxed_A,
  UnvsV11 = group_edgeR11_Unvaxed_A - group_edgeR11_Vaxed_A,
  
  UnvsV13 = group_edgeR13_Unvaxed_A - group_edgeR13_Vaxed_A,
  UnvsV14 = group_edgeR14_Unvaxed_A - group_edgeR14_Vaxed_A,
  UnvsV15 = group_edgeR15_Unvaxed_A - group_edgeR15_Vaxed_A,
  UnvsV16 = group_edgeR16_Unvaxed_A - group_edgeR16_Vaxed_A,
  UnvsV17 = group_edgeR17_Unvaxed_A - group_edgeR17_Vaxed_A,
  UnvsV18 = group_edgeR18_Unvaxed_A - group_edgeR18_Vaxed_A,
  UnvsV19 = group_edgeR19_Unvaxed_A - group_edgeR19_Vaxed_A,
  UnvsV20 = group_edgeR20_Unvaxed_A - group_edgeR20_Vaxed_A,
  UnvsV21 = group_edgeR21_Unvaxed_A - group_edgeR21_Vaxed_A,
  UnvsV22 = group_edgeR22_Unvaxed_A - group_edgeR22_Vaxed_A,
  
  UnvsV23 = group_edgeR23_Unvaxed_A - group_edgeR23_Vaxed_A,
  UnvsV24 = group_edgeR24_Unvaxed_A - group_edgeR24_Vaxed_A,
  UnvsV25 = group_edgeR25_Unvaxed_A - group_edgeR25_Vaxed_A,
  
  levels = colnames(design))


cols=as.character(pb_dge$samples$group)
table(cols)


### download https://github.com/YOU-k/voomByGroup/blob/main/voomByGroup.R 
source("voomByGroup.R")

###Both of VoomByGroup and voomWithQualityWeights were used  for analysis
#voomByGroup(pb_dge,design = design, group = group_edgeR, plot = "combine") -> y_v

voomWithQualityWeights(pb_dge, design = design, plot = TRUE) -> y_v



fit <- lmFit(y_v, design)
vfit <- contrasts.fit(fit, contrasts=contr.matrix)
tfit <- treat(vfit, lfc=0)



coef_names <- colnames(tfit$coefficients)
coef_names <- coef_names[coef_names != "Intercept"]

results_list <- lapply(coef_names, function(coef_name) {
  topTreat(tfit, coef = coef_name, number = Inf)
})


names(results_list) <- coef_names

all_results=data.frame()

for (name in names(results_list)) {
  
  c_results_list=results_list[[name]]
  c_results_list$clusters=str_replace_all(name,".+vsV","")
  c_results_list$symbol=rownames(c_results_list)
  all_results=rbind(all_results, c_results_list)
  
}

degso=read.csv("Cell/based_DEGs")
all_results1=merge(all_results,degso,by=c("clusters","symbol"))
cluster_c=unique.data.frame(all_Samples_f.list_NFSR.inted_annoted1@meta.data[,c("seurat_clusters","annoated_celltypes")])
all_results1=merge(all_results1,cluster_c,by.x="clusters",by.y="seurat_clusters")
write.csv(all_results1,"voom_DEGs.csv")


all_results2=(all_results1[all_results1$adj.P.Val<0.05 & (abs(all_results1$logFC)>0.59),])
write.csv(all_results2,"VoomDEGs_FC1.5FDR005.csv")
