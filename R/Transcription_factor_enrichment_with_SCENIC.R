##download human databases
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs

###### add 500 bp enrichments, Oct16th2023
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles){
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

library(Seurat)
library(SCENIC)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(RColorBrewer)

library(AUCell)
library(GENIE3)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)

library(ComplexHeatmap)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
#### data set


All_inte_Anned=all_Samples_f.list_NFSR.inted_annoted


All_inte_Anned_sub=subset(All_inte_Anned,downsample = 10000)

cellInfo <- All_inte_Anned@meta.data
exprMat <- as.matrix(All_inte_Anned@assays$RNA@data)
cellInfo=readRDS("int/cellInfo.Rds")


exprMat_1=exprMat[,colnames(exprMat) %in% rownames(cellInfo)]


list.files("int")
library(feather)
### Initialize settings


scenicOptions <- initializeScenic(org="hgnc" , dbDir="cisTarget_databases", nCores=10)

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo_1.Rds"



### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)



cellInfo_1=cellInfo[,c("state","major_cell_type","seurat_clusters")]

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo_1.Rds"
# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings


export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

