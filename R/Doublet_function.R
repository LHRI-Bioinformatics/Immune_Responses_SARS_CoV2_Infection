##check the doublelet

#BiocManager::install("scDblFinder")
library(scDblFinder)
library(scater)
library(scran)



doublet=function(all_Samples_f.list_NFSR.inted_annoted_sub){
  All.Integrated_sce=as.SingleCellExperiment(all_Samples_f.list_NFSR.inted_annoted_sub)
  set.seed(00010101)
  dec.mam <- modelGeneVarByPoisson(All.Integrated_sce)
  top.mam <- getTopHVGs(All.Integrated_sce, prop=0.1)
  library(BiocSingular)
  library(scRNAseq)
  set.seed(101010011)
  All.Integrated_sce1 <- denoisePCA(All.Integrated_sce, technical=dec.mam, subset.row=top.mam)
  All.Integrated_sce2 <- runTSNE(All.Integrated_sce1, dimred="PCA")
  snn.gr <- buildSNNGraph(All.Integrated_sce2, use.dimred="PCA", k=25)
  colLabels(All.Integrated_sce2) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
  saveRDS(All.Integrated_sce2,"All.Integrated_sce2.rds")
  All.Integrated_sce2=readRDS("All.Integrated_sce2.rds")
  dbl.out <- findDoubletClusters(All.Integrated_sce2)
  
  dbl.out=readRDS("All_samples03_14_2023dbl.out.rds")
  library(scater)
  chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
  chosen.doublet
  
  library(BiocSingular)
  set.seed(100)
  
  # Setting up the parameters for consistency with denoisePCA();
  # this can be changed depending on your feature selection scheme.
  dbl.dens <- computeDoubletDensity(All.Integrated_sce2, subset.row=top.mam, 
                                    d=ncol(reducedDim(All.Integrated_sce2)))
  
  All.Integrated_sce2$DoubletScore <- dbl.dens
  
  
  All.Integrated_sce3 <- runUMAP(All.Integrated_sce2, dimred="PCA")
  
  
  dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                   method="griffiths", returnType="call")
  
  
  
  
  All.Integrated_sce3$Doubletdefine <- dbl.calls
  All.Integrated_sce3$DoubletScore <- dbl.dens
  
  all_Samples_f.list_NFSR.inted_annoted_t=all_Samples_f.list_NFSR.inted_annoted
  All.Integrated_sce3$cellname=rownames(All.Integrated_sce3@colData)
  
  
  all_Samples_f.list_NFSR.inted_annoted_t$cellname=rownames(all_Samples_f.list_NFSR.inted_annoted_t@meta.data)
  
  
  ssw=merge(all_Samples_f.list_NFSR.inted_annoted_t@meta.data,All.Integrated_sce3@colData[,c("cellname","DoubletScore", "Doubletdefine")],
            by="cellname",all.x=T)
  
  ssw=as.data.frame(ssw)
  ssw$Doubletdefine=as.character(ssw$Doubletdefine)
  ssw$Doubletdefine[is.na(ssw$Doubletdefine)]="no_detected"
  rownames(ssw)=ssw$cellname
  all_Samples_f.list_NFSR.inted_annoted_t@meta.data <-ssw
  
  all_Samples_f.list_NFSR.inted_annoted_t
}




