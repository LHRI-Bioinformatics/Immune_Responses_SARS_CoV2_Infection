library(Platypus)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

library(dplyr)
library(ggplot2)

source("./Function.R")

list.files()

library(readxl)
library(Seurat)
library(stringr)
source("Functions_Nov.R")

Genomex_path="COUNTS READS after 10x GENOMICS pipeline"
setwd(Genomex_path)




#### get meta data
sample_meta1=read.csv("Samples_original_meta.csv")

data_list <- NULL
for (sample_cov in list.files("../couts/",pattern = "[0-9]_")){
  
  dir=paste('SARS_COvid/couts/',
            sample_cov,'/outs/multi/count/raw_feature_bc_matrix/',sep = "")
  ##featured antibody will list as another assay
  Trraw_data<-Read10X(data.dir=dir)
  TrSample <- CreateSeuratObject(counts = Trraw_data[[1]], project = sample_cov, min.cells = 3, min.features = 200)
  adt_assay <- CreateAssayObject(counts = Trraw_data[[2]][,names(TrSample$orig.ident)])
  TrSample[["ADT"]] <- adt_assay
  
  tcr_folder=paste('SARS_COvid/couts/',
                   sample_cov,'/outs/per_sample_outs/',sample_cov,sep = "")
  TrSample1=add_clonotype(tcr_folder,TrSample)
  TrSample1$group=sample_meta1$group[sample_meta1$sample_name==sample_cov]
  TrSample1$type=  sample_meta1$type[sample_meta1$sample_name==sample_cov]
  
  
  data_list=append(data_list,TrSample1)
}
COvid_data_list=data_list

###add public HC control data, download from GSE128033 and GSM3660650 and add them in meta data

sample_meta1[41,]=c("41","HC","C51","unvac","20")
sample_meta1[42,]=c("42","HC","C52","unvac","21")
sample_meta1[43,]=c("43","HC","C100","unvac","22")
sample_meta1[44,]=c("44","HC","GSM3660650","unvac","23")


sample_meta1$SRA="unknown"

sample_meta1$SRA[41:43]=c("SRR11537946" ,"SRR11537947", "SRR11537948")


##public data GSM3660650
Trraw_data<-Read10X(data.dir="Pulic_d_hc/GSM3660650_SC249NORbal_fresh")
TrSample <- CreateSeuratObject(counts = Trraw_data, project = "GSM3660650", min.cells = 3, min.features = 200)


data_list <- NULL
for (sample_cov in list.files(pattern = "SRR11")){
  
  dir=paste('refer_data/',
            sample_cov,'/outs/raw_feature_bc_matrix/',sep = "")
  ##featured antibody will list as another assay
  Trraw_data1<-Read10X(data.dir=dir)
  TrSample1 <- CreateSeuratObject(counts = Trraw_data1, project = sample_meta1$sample_name[sample_meta1$SRA==sample_cov], min.cells = 3, min.features = 200)
  
  
  TrSample1$group=sample_meta1$group[sample_meta1$SRA==sample_cov]
  TrSample1$type=  sample_meta1$type[sample_meta1$SRA==sample_cov]
  
 
  data_list=append(data_list,TrSample1)
}

HC_data=append(data_list,TrSample)


####combind public health data to our data
COvid_data_list1=append(COvid_data_list,HC_data)

data_lists=c()
sample_des=c(as.character(COvid_data_list1[[1]]$orig.ident[1]))

for (i in 2:44){
  data_l=COvid_data_list1[i]
  sample_des1=as.character(COvid_data_list1[[i]]$orig.ident[1])
  data_lists=append(data_lists,data_l)
  sample_des=append(sample_des,sample_des1)
  
  
}
length(data_lists)

all_Samples <- merge(COvid_data_list1[[1]], y = data_lists, add.cell.ids = sample_des, project = "SARS_samples")





all_Samples[["percent.mt"]] <- PercentageFeatureSet(all_Samples, pattern = "^MT-")
all_Samples <- PercentageFeatureSet(all_Samples, "^RP[SL]", col.name = "percent_ribo")

all_Samples <- PercentageFeatureSet(all_Samples, "^HB[^(P)]", col.name = "percent_hb")

all_Samples <- PercentageFeatureSet(all_Samples, "PECAM1|PF4", col.name = "percent_plat")

all_Samples@meta.data$group[all_Samples@meta.data$orig.ident=="GSM3660650"]="unvac"
all_Samples@meta.data$type[all_Samples@meta.data$orig.ident=="GSM3660650"]="HC"


### add patient info
all_Samples1=all_Samples
all_Samples1$patient <- 20
for (sample_name in sample_meta1$sample_name){
  
  patient_id=sample_meta1$patient[sample_meta1$sample_name == sample_name]
  all_Samples1$patient[all_Samples1$orig.ident == sample_name] <- patient_id
  
}



all_Samples1$state=paste(all_Samples1$group,all_Samples1$type,sep = "_")





all_Samples1$log10GenesPerUMI <- log10(all_Samples1$nFeature_RNA) / log10(all_Samples1$nCount_RNA)

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo", "percent_hb")

###filter QC,as we have the reference HC control, we used their filter method, ncouts>1000, gene >200, remove mt gene over 15%

all_Samples_f <- subset(all_Samples1, subset = nFeature_RNA > 200 &  nCount_RNA > 1000 & percent.mt < 15)


all_Samples_f.list <- SplitObject(all_Samples_f, split.by = "orig.ident")
all_Samples_f.list_NF <- lapply(X = all_Samples_f.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})


features <- SelectIntegrationFeatures(object.list = all_Samples_f.list_NF)
all_Samples_f.list_NFSR <- lapply(X = all_Samples_f.list_NF, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})




###find integration reference male and female unvacinced R,sample 2 and 6, their position at 12,37


anchors <- FindIntegrationAnchors(object.list = all_Samples_f.list_NFSR, reference = c(12, 37), reduction = "rpca",dims = 1:50)

all_Samples_f.list_NFSR.inted <- IntegrateData(anchorset = anchors, dims = 1:50)

all_Samples_f.list_NFSR.inted <- ScaleData(all_Samples_f.list_NFSR.inted, verbose = FALSE)

all_Samples_f.list_NFSR.inted <- RunPCA(all_Samples_f.list_NFSR.inted, verbose = FALSE)

all_Samples_f.list_NFSR.inted=readRDS("all_Samples_f.list_NFSR.inted_p.rds")
all_Samples_f.list_NFSR.inted <- RunUMAP(all_Samples_f.list_NFSR.inted, dims = 1:50)

all_Samples_f.list_NFSR.inted=readRDS("all_Samples_f.list_NFSR.inted_u.rds")

all_Samples_f.list_NFSR.inted <- RunTSNE(all_Samples_f.list_NFSR.inted,dims = 1:50)


all_Samples_f.list_NFSR.inted <- FindNeighbors(all_Samples_f.list_NFSR.inted, reduction = "pca", dims = 1:50)

all_Samples_f.list_NFSR.inted <- FindClusters(all_Samples_f.list_NFSR.inted, resolution = 0.5)













uniq_annoted=read.csv("Annotation_table.csv")
uniq_annoted=uniq_annoted[order(uniq_annoted$seurat_clusters,decreasing = F),]
uniq_annoted_cl=uniq_annoted$major_cell_type
names(uniq_annoted_cl) <- levels(all_Samples_f.list_NFSR.inted)
all_Samples_f.list_NFSR.inted_annoted <- RenameIdents(all_Samples_f.list_NFSR.inted, uniq_annoted_cl)


all_Samples_f.list_NFSR.inted_annoted$major_cell_type=Idents(all_Samples_f.list_NFSR.inted_annoted)

saveRDS(all_Samples_f.list_NFSR.inted_annoted,"all_samples_annotated.rds")








library(dplyr)
library(ggplot2)








### calculate the cell population of each cell type 

CELL_POPULATION=all_Samples_f.list_NFSR.inted_annoted@meta.data[,c("orig.ident","seurat_clusters","state","major_cell_type" ,"patient","group","type")]

CELL_POPULATION$cluster_state=paste(CELL_POPULATION$seurat_clusters,CELL_POPULATION$state,sep = "_")
CELL_POPULATION$cluster_orig_name=paste(CELL_POPULATION$seurat_clusters,CELL_POPULATION$orig.ident,sep = "_")
cluster_state_cells=as.data.frame(table(CELL_POPULATION$cluster_orig_name))
CELL_POPULATION_cluster_Origname=unique.data.frame(CELL_POPULATION[,c(1,2,3,4,8,9)])
colnames(cluster_state_cells)=c("cluster_orig_name" ,"Freq")

cluster_state_cells=merge(cluster_state_cells,CELL_POPULATION_cluster_Origname,by="cluster_orig_name")

Idents(all_Samples_f.list_NFSR.inted_annoted)[1:3]
for (x in 0:35){
  cluster_state_cells$total_cells[cluster_state_cells$seurat_clusters==x]=sum(cluster_state_cells$Freq[cluster_state_cells$seurat_clusters==x])
}
cluster_state_cells$Freq_Percent=(cluster_state_cells$Freq/cluster_state_cells$total_cells)*100

cluster_state_cells$seurat_clusters=factor(cluster_state_cells$seurat_clusters,levels = unique(cluster_state_cells$seurat_clusters[order(cluster_state_cells$major_cell_type)]))
cluster_state_cells$state=factor(cluster_state_cells$state,levels =c("unvac_A","unvac_R","unvac_C", "vacced_A","vacced_R","vacced_C","unvac_HC"))
##as cluster32-35 only have 2 cells, so remove these clusters
cluster_state_cells1= cluster_state_cells[!cluster_state_cells$seurat_clusters %in% 32:35,]



for (x in unique(cluster_state_cells$orig.ident)){
  cluster_state_cells$total_cells_origname[cluster_state_cells$orig.ident==x]=sum(cluster_state_cells$Freq[cluster_state_cells$orig.ident==x])
}


cluster_state_cells$cluster_statetotalcell_Freq=cluster_state_cells$Freq/cluster_state_cells$total_cells_origname

####save cluster_state_cells for population plotting


