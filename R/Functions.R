library("plotrix")
add_clonotype <- function(tcr_folder, seurat_obj){
  vdj_data=data.frame()
  for (vdj_type in c("/vdj_b/","/vdj_t/")){
    if ((file.size(paste(tcr_folder,vdj_type,"filtered_contig_annotations.csv", sep=""))>0) & (file.size(
      paste(tcr_folder,vdj_type,"clonotypes.csv", sep=""))>0)){
      tcr <- read.csv(paste(tcr_folder,vdj_type,"filtered_contig_annotations.csv", sep=""))
      # Remove the -1 at the end of each barcode.
      # Subsets so only the first line of each barcode is kept,
      # as each entry for given barcode will have same clonotype.
      #no need to rm -1 in names as it is keep in seurat object
      #tcr$barcode <- gsub("-1", "", tcr$barcode)
      tcr <- tcr[!duplicated(tcr$barcode), ]
      # Only keep the barcode and clonotype columns. 
      # We'll get additional clonotype info from the clonotype table.
      tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene")]
      names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
      # Clonotype-centric info.
      clono <- read.csv(paste(tcr_folder,vdj_type,"clonotypes.csv", sep=""))
      # Slap the AA sequences onto our original table by clonotype_id.
      tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
      # Reorder so barcodes are first column and set them as rownames.
      tcr <- tcr[, c(2,1,3:7)]
      rownames(tcr) <- tcr[,1]
      tcr[,1] <- NULL
      vdj_data=rbind(vdj_data,tcr)
      
    }
  }
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=vdj_data)
  return(clono_seurat)
}











## unvac_C only have one sample

data_pvalue=function(aqw,column,value){
  
  
  unvac=c("unvac_A", "unvac_R", "unvac_A", "vacced_A","unvac_C")
  vacced=c("vacced_A", "vacced_R", "unvac_R","vacced_R","vacced_C")
  stat_data=data.frame()
  
  for (i in 1:length(unvac)){
    for (x_cell in unique(aqw[,column])){
      x_v=aqw[,value][(aqw[,column]==x_cell)&(aqw$state==unvac[i])]
      y_v=aqw[,value][(aqw[,column]==x_cell)&(aqw$state==vacced[i])]
      if ((length(x_v)==0)|(length(y_v)==0)) next
      else if (length(x_v)==1 & length(y_v)==1) next
      else if (length(x_v)==1){
        se=t.test(mu=x_v,x=y_v)
        se_df=data.frame(state=paste(unvac[i],vacced[i],sep = " vs "),
                         celltype=x_cell,
                         mean_x=x_v,
                         se_x=0,
                         se_y=std.error(y_v),
                         mean_y=se$estimate[[1]],
                         p_value=se$p.value,
                         stderr=se$stderr)
      }else if (length(y_v)==1){
        se=t.test(mu=y_v,x=x_v)
        se_df=data.frame(state=paste(unvac[i],vacced[i],sep = " vs "),
                         celltype=x_cell,
                         mean_x=se$estimate[[1]],
                         mean_y=y_v,
                         se_x=std.error(x_v),
                         se_y=0,
                         p_value=se$p.value,
                         stderr=se$stderr)
      } else {
        se=t.test(x_v,y_v)
        se_df=data.frame(state=paste(unvac[i],vacced[i],sep = " vs "),
                         celltype=x_cell,
                         mean_x=se$estimate[[1]],
                         mean_y=se$estimate[[2]],
                         se_x=std.error(x_v),
                         se_y=std.error(y_v),
                         p_value=se$p.value,
                         stderr=se$stderr)}
      
      
      stat_data=rbind(stat_data,se_df)
      
    }
    
  }
  
  return(stat_data)
}

library(data.table)
library(Matrix)
library(readr)



Covid_gene_matrix=function(dir){
  setwd(dir)
  Covid.d=data.frame(c('CD8A',"orf1ab","S","ORF3a", "E","M","ORF6","ORF7a", "ORF8","N","ORF10" ))
  
    
  for (x in list.files()){
    mat <- readMM(file = paste(x,'/outs/multi/count/raw_feature_bc_matrix/matrix.mtx.gz',sep = ''))
    feature.names<- read.delim(paste(x,"/outs/multi/count/raw_feature_bc_matrix/features.tsv.gz",sep = ""),
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(paste(x,"/outs/multi/count/raw_feature_bc_matrix/barcodes.tsv.gz",sep = ""),
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V2
    
    
    HIVmat<-as.data.frame(rowSums(mat[rownames(mat) %in% 
                                        (c('CD8A',"orf1ab","S","ORF3a", "E","M","ORF6","ORF7a", "ORF8","N","ORF10" )),])
    )
    colnames(HIVmat)=x
    Covid.d=cbind(Covid.d,as.data.frame(HIVmat))
    
  }
  
  
  return(t(Covid.d))
  
}




###plot the population of each cell types Jan16th 2023
PUPULATION_POLOT=function(data,celltype_cluster,cluster_statetotalcell_Freq,geom_text,ncol,y_title){
  cluster_state_cells_stats=data_pvalue(data,celltype_cluster,cluster_statetotalcell_Freq)
  
  cluster_state_cells_stats_d=cluster_state_cells_stats[cluster_state_cells_stats$state %in% c("unvac_A vs vacced_A","unvac_R vs vacced_R","unvac_C vs vacced_C"),]
  
  cluster_state_cells_stats_d$p_value_d=round(cluster_state_cells_stats_d$p_value,3)
  cluster_state_cells_stats_d$group=factor(str_replace_all(cluster_state_cells_stats_d$state,".+_",""),levels = c("A","R","C"))
  colors <- c("Unvaxed" = "blue", "Vaxed" = "red")
  
  
  
  p=ggplot(data = cluster_state_cells_stats_d,aes(x = group)) + 
    geom_line( aes( y = mean_x,group = 1, color = "Unvaxed"))+
    geom_point( aes( y = mean_x,group = 1, color = "Unvaxed")) +
    geom_errorbar(aes(ymin=mean_x-se_x, ymax=mean_x+se_x, color = "Unvaxed"), width=.2,position=position_dodge(0.05))+
    
    geom_line( aes(y = mean_y,group = 1, color = "Vaxed"))+
    geom_point( aes( y = mean_y,group = 1, color = "Vaxed")) +
    geom_errorbar(aes(ymin=mean_y-se_y, ymax=mean_y+se_y,color = "Vaxed"), width=.2,position=position_dodge(0.05))+
    
    ylab(y_title) + facet_wrap(~ celltype,scales = "free",ncol = ncol)
  
  
  if (geom_text){
    p=p+geom_text(aes(x = group, y = Inf,label=paste("p:",p_value_d),hjust = 0.5, vjust = 1,group = 1))
  }
  
  print(p)
    
  
  
 
  
}








# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="CD4+_N"){
  cell_ids <- which(colData(cds)[, "annotation_summary"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
# a helper function to identify the root principal points:
get_earliest_principal_node_CD8 <- function(cds, time_bin="CD8_N"){
  cell_ids <- which(colData(cds)[, "annotation_summary"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

##MACROPHAGE CELLS, FOLLOW THE REFERENCE PAPER, FCN1 HI AS THE ORIGINAL CLUSTER, WE USE THE Mac_CSF1R WHICH HAS HIGH FCN1
get_earliest_principal_node_Macros <- function(cds, time_bin="Mac_CSF1R"){
  cell_ids <- which(colData(cds)[, "Annotation_summary"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}




#### using FindMarker function to call DEG in Covid samples





CallDEG=function(all_Samples_f.list_NFSR.inted_annoted,obj){
  
  ##FIDN DIFFERENT GENES BETWEEN TRO GROUP
  
  
  
  all_Samples_f.list_NFSR.inted_annoted$celltype.stim<-paste((all_Samples_f.list_NFSR.inted_annoted$seurat_clusters),
                                                             all_Samples_f.list_NFSR.inted_annoted$state,sep = '_')
  
  Idents(all_Samples_f.list_NFSR.inted_annoted)<-'celltype.stim'
  
  
  All_idents<-as.character(unique(Idents(all_Samples_f.list_NFSR.inted_annoted)))
  
  library(DESeq2)
  library(stringr)
  
  library(limma)
  
  infname=c("unvac_A","unvac_R",  "unvac_C",  "unvac_A",  "vacced_A", "unvac_R","vacced_R")
  samplename=c("vacced_A",  "vacced_R",  "vacced_C",  "unvac_R",  "vacced_R", "unvac_HC","unvac_HC")
  
  for (i in 1:7){
    ass2<-data.frame()
    for (xx in unique(all_Samples_f.list_NFSR.inted_annoted$seurat_clusters)){
      
      
      ident.2=All_idents[grep(paste("^",xx,"_",samplename[i],sep = ""),All_idents)]
      ident.1=All_idents[grep(paste("^",xx,"_",infname[i],sep = ""),All_idents)]
      print(ident.1)
      if (length(all_Samples_f.list_NFSR.inted_annoted$celltype.stim[all_Samples_f.list_NFSR.inted_annoted$celltype.stim==ident.1])>3 &
          length(all_Samples_f.list_NFSR.inted_annoted$celltype.stim[all_Samples_f.list_NFSR.inted_annoted$celltype.stim==ident.2])>3){
        ass1<-FindMarkers(all_Samples_f.list_NFSR.inted_annoted,ident.2 = ident.2,ident.1 = ident.1, 
                          verbose = FALSE,pseudocount.use = 0.0001)    
        ass1$clusters<-xx
        
        ass1$symbol<-rownames(ass1)
        
        
        
        ass2<-rbind(ass2,ass1)
        print("3")
      }
      
      
    }
    ass2=ass2[order(ass2$clusters),]
    write.csv(ass2,file = paste(infname[i],'_vs_',samplename[i],'_altered_markers_',obj,'.csv',sep = ''),row.names = F)
    ass3<-ass2[(ass2$p_val<0.05)&((ass2$avg_logFC>0.59)|(ass2$avg_logFC<(-0.59)))&((ass2$pct.1>0.25)|(ass2$pct.2>0.25)),]
    write.csv(ass3,file = paste(infname[i],'_vs_',samplename[i],'_altered_markersFC1.5p0.05PCT25_',obj,'.csv',sep = ''),row.names = F)
    
  }
  
}






CallDEP=function(all_Samples_f.list_NFSR.inted_annoted,obj){
  ### as protein have 136, CLR was used in the normalization method
  ##FIDN DIFFERENT GENES BETWEEN TwO GROUP
  
  DefaultAssay(all_Samples_f.list_NFSR.inted_annoted)="ADT"
  
  all_Samples_f.list_NFSR.inted_annoted <- NormalizeData(all_Samples_f.list_NFSR.inted_annoted, assay = "ADT", normalization.method = "CLR",margin = 2)
  all_Samples_f.list_NFSR.inted_annoted <- ScaleData(all_Samples_f.list_NFSR.inted_annoted, assay = "ADT")
  
  all_Samples_f.list_NFSR.inted_annoted$celltype.stim<-paste((all_Samples_f.list_NFSR.inted_annoted$seurat_clusters),
                                                             all_Samples_f.list_NFSR.inted_annoted$state,sep = '_')
  
  Idents(all_Samples_f.list_NFSR.inted_annoted)<-'celltype.stim'
  
  
  All_idents<-as.character(unique(Idents(all_Samples_f.list_NFSR.inted_annoted)))
  
  library(DESeq2)
  library(stringr)
  
  library(limma)
  
  infname=c("unvac_A","unvac_R",  "unvac_C",  "unvac_A",  "vacced_A", "unvac_R","vacced_R")
  samplename=c("vacced_A",  "vacced_R",  "vacced_C",  "unvac_R",  "vacced_R", "unvac_HC","unvac_HC")
  
  for (i in 1:7){
    ass2<-data.frame()
    for (xx in unique(all_Samples_f.list_NFSR.inted_annoted$seurat_clusters)){
      
      
      ident.2=All_idents[grep(paste("^",xx,"_",samplename[i],sep = ""),All_idents)]
      ident.1=All_idents[grep(paste("^",xx,"_",infname[i],sep = ""),All_idents)]
      print(ident.1)
      if (length(all_Samples_f.list_NFSR.inted_annoted$celltype.stim[all_Samples_f.list_NFSR.inted_annoted$celltype.stim==ident.1])>3 &
          length(all_Samples_f.list_NFSR.inted_annoted$celltype.stim[all_Samples_f.list_NFSR.inted_annoted$celltype.stim==ident.2])>3){
        ass1<-FindMarkers(all_Samples_f.list_NFSR.inted_annoted,ident.2 = ident.2,ident.1 = ident.1, 
                          verbose = FALSE,pseudocount.use = 0.0001)    
        ass1$clusters<-xx
        
        ass1$symbol<-rownames(ass1)
        
        
        
        ass2<-rbind(ass2,ass1)
        
      }
      
      
    }
    ass2=ass2[order(ass2$clusters),]
    write.csv(ass2,file = paste(infname[i],'_vs_',samplename[i],'_altered_markers_Protein',obj,'.csv',sep = ''),row.names = F)
    ass3<-ass2[(ass2$p_val<0.05)&((ass2$avg_logFC>0.59)|(ass2$avg_logFC<(-0.59)))&((ass2$pct.1>0.25)|(ass2$pct.2>0.25)),]
    write.csv(ass3,file = paste(infname[i],'_vs_',samplename[i],'_altered_markersFC1.5p0.05PCT25_Protein',obj,'.csv',sep = ''),row.names = F)
    
  }
  
}