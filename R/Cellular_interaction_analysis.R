library(sceasy)
library(reticulate)
use_condaenv('your_path_software/Anaconda3/2022.05/envs/R-4.1.3')
loompy <- reticulate::import('loompy')
library(Seurat)
library(sceasy)
library(reticulate)

#loompy <- reticulate::import('loompy')
library(Seurat)
library(stringr)
library(ktplots)
library(SingleCellExperiment)

All_inte_Anned=readRDS("all_Samples.rds")


sed=sceasy::convertFormat(All_inte_Anned, from="seurat", to="anndata",
                     outFile='all_Samples_seura.h5ad')

All_inte_Anned_UnA=subset(All_inte_Anned,subset=state == "Unvax_A")

All_sce <- as.SingleCellExperiment(seurat_obj)


##make Cellphone data June2nd2023



###same method to generate VaxedA, Unvaxed_R and Vaxed_R input data





##### Acute downstream analysis


naive_means <- read.delim("UnvaxedA/method2/statistical_analysis_means_09_20_2023_12:37:17.txt", check.names = FALSE)
naive_pvals <- read.delim("UnvaxedA/method2/statistical_analysis_pvalues_09_20_2023_12:37:17.txt", check.names = FALSE)
naive_decon <- read.delim("UnvaxedA/method2/statistical_analysis_deconvoluted_09_20_2023_12:37:17.txt", check.names = FALSE)
naive_pvals_o <- read.delim("../All_sample_results/UnvaxedA/degs_analysis_relevant_interactions_result_06_07_2023_14:42:46.txt", check.names = FALSE)


treated_means <- read.delim("VaxedA/method2/statistical_analysis_means_09_20_2023_13:00:10.txt", check.names = FALSE)
treated_pvals <- read.delim("VaxedA/method2/statistical_analysis_pvalues_09_20_2023_13:00:10.txt", check.names = FALSE)
treated_decon <- read.delim("VaxedA/method2/statistical_analysis_deconvoluted_09_20_2023_13:00:10.txt", check.names = FALSE)
treated_pvals_In <- read.delim("../All_sample_results/VaxedA//degs_analysis_relevant_interactions_result_06_07_2023_14:49:53.txt", check.names = FALSE)

means <- combine_cpdb(naive_means, treated_means)
pvals <- combine_cpdb(naive_pvals, treated_pvals)
decon <- combine_cpdb(naive_decon, treated_decon)




All_sce$state=str_replace_all(All_sce$state,"unvac","Unvax")
All_sce$state=str_replace_all(All_sce$state,"vacced","Vaxed")




breaksList = seq(0, 150, by =0.15)
symmetrical=TRUE
naive_pvals1=naive_pvals[,colnames(naive_pvals)[!grepl("unknow",colnames(naive_pvals))]]
colnames(naive_pvals1)=str_replace_all(colnames(naive_pvals1),"Unvax_A_","")


pdf(file = (paste("UnvaxedA_interactions_",  format(Sys.Date(), "%m_%d_%Y"),".PDF",
                  sep = "")),height =8,width = 10)
p1=plot_cpdb_heatmap(pvals =  naive_pvals1,  symmetrical = symmetrical,angle_col = 90,
                     cluster_rows = FALSE,cluster_cols = FALSE,breaks = breaksList,fontsize_row=35,fontsize_col =35, 
                     fontsize = 25,low_col = "blue",high_col = "red")
print(p1)
dev.off()


treated_pvals1=treated_pvals[,colnames(treated_pvals)[!grepl("unknow",colnames(treated_pvals))]]
colnames(treated_pvals1)=str_replace_all(colnames(treated_pvals1),"Vaxed_A_","")



pdf(file = (paste("VaxedA_interactions_",  format(Sys.Date(), "%m_%d_%Y"),".pdf",
                  sep = "")),height =8,width = 10)
p2=plot_cpdb_heatmap(pvals = treated_pvals1,  symmetrical = symmetrical,angle_col = 90,
                     cluster_rows = FALSE,cluster_cols = FALSE,breaks = breaksList,fontsize_row=35,fontsize_col =35, 
                     fontsize = 25,low_col = "blue",high_col = "red")
print(p2)
dev.off()

#### Epithelial receptors Jan14th2025
keep_significant_only = "Significant"
dir.create("Acute")
setwd("Acute/")### Modified on Oct30th2024

pdf(file = (paste("Acute_Epithelial_Recptors_Levles_",  format(Sys.Date(), "%m_%d_%Y"),".pdf",
                  sep = "")),height =7,width = 15)
G=plot_cpdb(cell_type1 = 'Epithelial', cell_type2 = 'Macrophage|CD4|CD8|B_cell|Mast|mDC|pDC', scdata = All_sce,
            celltype_key = 'major_cell_type', # column name where the cell ids are located in the metadata
            splitby_key = 'state', # column name where the grouping column is. Optional.
            means = means, pvals = pvals,degs_analysis = FALSE,keep_significant_only=TRUE,
            genes = c("EREG","ACE2","TMPRSS2","AXL","NRP1","CTSL","CTSB","SCARB1","BSG","ADAM17","SLC6A19","RPS2","CLEC4M")) +
  small_axis(fontsize = 20) +  small_legend(fontsize = 20) +
  theme(legend.key.size = unit(0.7, 'cm'))+theme(axis.text.x = element_text( size = 20,color="black",face = "bold"),axis.text.y = element_text(
    size = 20,color="black",face = "bold"),axis.title =element_text(size = 20,color="black",face = "bold"),
    legend.text =element_text(size = 20) ,#strip.text=element_text(size=20),
    
   
  )+
  scale_size(range = c(1.2, 8))+
  theme(legend.direction = "vertical", legend.box = "vertical")
G$data$Var2=str_replace_all(G$data$Var2,"-.+_A_","-")
G$data$Var2=str_replace_all(G$data$Var2,"(axed_A)|(vax_A)","")
G$data$Var2=str_replace_all(G$data$Var2,"(thelial)|(rophage)","")
S1=as.character(unique(G$data$Var2))

G$data$Var2 <- factor(G$data$Var2, levels = S1[order(S1)])
print(G)

dev.off()




#### Means files and pvals file should keep in same order!!!


naive_means1=naive_means[,!grepl("(Mast)|(Plasma)|(unknown)|(B_cell)",colnames(naive_means))]
naive_pvals1=naive_pvals[,!grepl("(Mast)|(Plasma)|(unknown)|(B_cell)",colnames(naive_pvals))]
colnames(naive_means1)=str_replace_all(colnames(naive_means1),"Unvax_A_","")
colnames(naive_pvals1)=str_replace_all(colnames(naive_pvals1),"Unvax_A_","")
naive_decon1=naive_decon
colnames(naive_decon1)=str_replace_all(colnames(naive_decon1),"Unvax_A_","")
All_sce_UNA=All_sce[,(colData(All_sce)$state=="Unvax_A")]
library(ComplexHeatmap)

cell_type1 = 'Unvax_A_Epithelial'
cell_type2 = 'Unvax_A_CD4|Unvax_A_CD8|Unvax_A_Macrophage|Unvax_A_NK|Unvax_A_mDC|Unvax_A_pDC'
cell_type1=unlist(str_split(cell_type1,"\\|"))

cell_type2=unlist(str_split(cell_type2,"\\|"))
paste(cell_type1,cell_type2,sep = "|")

in_f=outer(cell_type1, cell_type2, FUN = "paste",sep="|")
dim(in_f) =NULL
in_r=outer(cell_type2, cell_type1, FUN = "paste",sep="|")
dim(in_r) =NULL

all_Uvc=c(in_r,in_f)
all_Vc=str_replace_all(all_Uvc,"Unvax","Vaxed")

uniqe_n = naive_pvals[(naive_pvals$interacting_pair %in% naive_pvals_o$interacting_pair) &
                        (rowMins(as.matrix(naive_pvals[,all_Uvc]))<=0.05),]

uniqe_t=treated_pvals[(treated_pvals$interacting_pair %in% treated_pvals_In$interacting_pair)&
                        (rowMins(as.matrix(treated_pvals[,all_Vc]))<=0.05),]
test_list2 <- c("#e8e8e8",   "#342020",    "#b7b7b7",     "#606060",       "black",     "#404040",   "#525252")
names(test_list2) <- unique(All_sce_UNA$major_cell_type)[!grepl("(Mast)|(Plasma)|(unknown)",unique(All_sce_UNA$major_cell_type))][c(1:6,8)]

#### UnvaxedA vs VaxA uniques
interactions=str_replace_all(naive_means$interacting_pair[naive_means$interacting_pair %in% uniqe_n$interacting_pair[!uniqe_n$interacting_pair %in% uniqe_t$interacting_pair]],"_","-")


interactions=str_replace_all(interactions,"_","-")

pdf(file = (paste( "Epithelial_Recptors_Levles_UnvaxedAvsVaxed_Circleplot_Uinques",format(Sys.Date(), "%m_%d_%Y"),".pdf",sep="")),height =3,width = 9.5)

p <- plot_cpdb4(interaction =interactions,cell_type1='Epithelial', cell_type2 = 'CD4|CD8|Macrophage|NK|mDC|pDC',
                scdata = All_sce_UNA,
                celltype_key  = 'major_cell_type',
              
                means = naive_means1[naive_means1$interacting_pair %in% naive_pvals_o$interacting_pair,],
                pvals = naive_pvals1[naive_pvals1$interacting_pair %in% naive_pvals_o$interacting_pair,],
                desiredInteractions = list(c('Epithelial','CD8'),
                                           c('Epithelial','CD4'),
                                           c('Epithelial','Macrophage'),
                                           c('Epithelial','NK'),
                                           c('Epithelial','mDC'),
                                           c('Epithelial','pDC'),
                                           c('CD8','Epithelial'),
                                           c('CD4','Epithelial'),
                                           c('Macrophage','Epithelial'),
                                           c('NK','Epithelial'),
                                           c('mDC','Epithelial'),
                                           c('pDC','Epithelial')
                ),
                deconvoluted =  naive_decon1,keep_significant_only = TRUE,
                standard_scale = TRUE,
                grid_colors	=test_list2,
                
                alpha = 1,frac = 0.5,
                
                remove_self = TRUE) +guides(fill=guide_legend(ncol=2))
print(p)
dev.off()



interactions=str_replace_all(treated_means$interacting_pair[treated_means$interacting_pair %in% uniqe_t$interacting_pair[!uniqe_t$interacting_pair %in% uniqe_n$interacting_pair]],"_","-")
treated_means1=treated_means[,!grepl("(Mast)|(Plasma)|(unknown)|(B_cell)",colnames(treated_means))]
treated_pvals1=treated_pvals[,!grepl("(Mast)|(Plasma)|(unknown)|(B_cell)",colnames(treated_pvals))]
colnames(treated_means1)=str_replace_all(colnames(treated_means1),"Vaxed_A_","")
colnames(treated_pvals1)=str_replace_all(colnames(treated_pvals1),"Vaxed_A_","")
treated_decon1=treated_decon
colnames(treated_decon1)=str_replace_all(colnames(treated_decon1),"Vaxed_A_","")
All_sce_VA=All_sce[,(colData(All_sce)$state=="Vaxed_A")]
test_list2 <- c("#e8e8e8",   "#342020",    "#b7b7b7",     "#606060",       "black",     "#404040",   "#525252")
names(test_list2) <- unique(All_sce_UNA$major_cell_type)[!grepl("(Mast)|(Plasma)|(unknown)",unique(All_sce_UNA$major_cell_type))][c(1:6,8)]



pdf(file = (paste("Epithelial_Recptors_Levles_vaxedAvSUnvaxed_Circleplot_Uniques_",  format(Sys.Date(), "%m_%d_%Y"),".pdf",
                  sep = "")),height =5.5,width = 13)

p <- plot_cpdb4(interaction =interactions,cell_type1 = 'Epithelial', cell_type2 = 'CD4|CD8|Macrophage|NK|mDC|pDC',
                scdata = All_sce_VA,
                celltype_key  = 'major_cell_type',
                #split.by = 'state',# column name where the cell ids are located in the metadata
                means = treated_means1,
                pvals = treated_pvals1,desiredInteractions = list(c('Epithelial','CD8'),
                                                                  c('Epithelial','CD4'),
                                                                  c('Epithelial','Macrophage'),
                                                                  c('Epithelial','NK'),
                                                                  c('Epithelial','mDC'),
                                                                  c('Epithelial','pDC'),
                                                                  c('CD8','Epithelial'),
                                                                  c('CD4','Epithelial'),
                                                                  c('Macrophage','Epithelial'),
                                                                  c('NK','Epithelial'),
                                                                  c('mDC','Epithelial'),
                                                                  c('pDC','Epithelial')
                ),
                deconvoluted = treated_decon1,keep_significant_only = TRUE,
                standard_scale = TRUE,
                remove_self = TRUE,grid_colors	=test_list2,
                
                alpha = 1
)
p
dev.off()



###### Recovery sample using same method to analysis 






