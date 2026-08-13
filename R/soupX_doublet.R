#install.packages("SoupX")
library(SoupX)
library(Seurat)
library(scDblFinder)

contamination_data=data.frame()
seurat_list=list()
####All_the_samples is the name of all our samples
for (sample_id in All_the_samples){
  row_dir=paste('/counts_path/couts/',sample_id,'/outs/multi/count/raw_feature_bc_matrix/',sep = "")
  
  
  fil_dir=paste('/counts_path/couts/',sample_id,'/outs/per_sample_outs/',sample_id,'/count/sample_filtered_feature_bc_matrix/',sep = "")
  
  raw_matrix  <- Read10X(row_dir)
  
  filt_matrix <- Read10X(fil_dir)
  
  
  
  names(raw_matrix)
  names(filt_matrix)
  
  raw_matrix  <- raw_matrix$`Gene Expression`
  filt_matrix <- filt_matrix$`Gene Expression`
  
  #sc <- SoupChannel(raw_matrix, filt_matrix)
  
  srt <- CreateSeuratObject(filt_matrix)
  srt <- RenameCells(srt, add.cell.id = sample_id)  # add prefix consistently
  srt <- NormalizeData(srt)
  srt <- FindVariableFeatures(srt)
  srt <- ScaleData(srt)
  srt <- RunPCA(srt)
  srt <- FindNeighbors(srt, dims = 1:20)
  srt <- FindClusters(srt)
  srt <- RunUMAP(srt, dims = 1:20)
  
  sc <- SoupChannel(raw_matrix, filt_matrix)
  
  
  
  clusters <- srt$seurat_clusters
  sc <- setClusters(sc, setNames(srt$seurat_clusters, colnames(srt)))
  
  sc <- setDR(sc, srt@reductions$umap@cell.embeddings)
  
  # Estimate contamination and correct
  sc <- autoEstCont(sc,contaminationRange = c(0.00001, 0.8),forceAccept = TRUE)
  adj_matrix <- adjustCounts(sc)
  
  # Create Seurat object with corrected counts, tag with sample name
  srt_clean <- CreateSeuratObject(counts = adj_matrix, project = sample_id)
  srt_clean$sample <- sample_id
  
  seurat_list[[sample_id]] <- srt_clean
  #colnames(srt_filtered)[!colnames(srt_filtered) %in% rownames(sc$metaData)]
  
  
  
  
  
  print(paste(sample_id,'the unique the RHO is',unique(sc$metaData$rho)))
 rownames( sc$metaData)=paste(sample_id,rownames(sc$metaData),sep = "_")
 contamination_data=rbind(contamination_data,sc$metaData)
  
}
 


library(stringr)
pdf("../../soupX/Rho_values_wholes_samples_Aug4th2026.pdf",width = 10,height = 8)

contamination_data$sample_id_old=contamination_data$sample_id
contamination_data$sample_id=str_replace_all(contamination_data$sample_id,"_.+","")
p=ggplot(unique.data.frame(contamination_data[,c("rho","sample_id")]), aes(x = reorder(sample_id, rho), y = rho)) +
  geom_col(fill = "steelblue") +
  #geom_hline(yintercept = 0.03, linetype = "dashed", color = "red") +
  labs(x = "Sample", y = "SoupX estimated rho (contamination fraction)",
       title = "Estimated ambient RNA contamination across samples") +
  #theme_minimal() +
  theme(axis.text = element_text(angle = 90, hjust = 1, size = 15),panel.grid = element_blank(), panel.background =  element_blank(), 
        axis.title  = element_text( size = 20),axis.line = element_line(size = 1),title = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 
print(p)
dev.off()
#ggsave("Rho_values_wholes_samples.pdf", width = 10, height = 8)

### doulet checking

seurat.list <- lapply(seurat_list, function(obj) {
  DefaultAssay(obj) <- "RNA"
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, BPPARAM = SerialParam())  # no samples= needed, each list element IS one sample
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  obj
})
saveRDS(seurat.list,"//Alls_samples_SoupX__dobuleted_wholedata.rds")
