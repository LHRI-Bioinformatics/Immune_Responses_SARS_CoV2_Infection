
####conver Seurat to loupe browser
#install.packages("hdf5r")
#install.packages("remotes")
#remotes::install_github("10XGenomics/loupeR")
#loupeR::setup()
library(loupeR)
library(hdf5r)
library(Seurat)
testrds=readRDS("T_NK.Integrated_annoted_Seurat.rds")

DefaultAssay(testrds)="RNA"

testrds2=subset(testrds,subset = cdr3s_aa %in% testrds$cdr3s_aa[!is.na(testrds$cdr3s_aa)])

testrds2@meta.data$d_gene[ nchar(testrds2@meta.data$d_gene)==0]="NA"
testrds2@meta.data$v_gene[ nchar(testrds2@meta.data$v_gene)==0]="NA"
testrds2@meta.data$j_gene[ nchar(testrds2@meta.data$j_gene)==0]="NA"
testrds2@meta.data$cdr3s_aa[ nchar(testrds2@meta.data$cdr3s_aa)==0]="NA"


testrds4=testrds2

VDJ_order=read.csv('integration_vloupe_files_Sample_sheet' )



library(stringr)

VDJ_order1=as.data.frame(VDJ_order$sample_id)
VDJ_order1$Samplesorder=rownames(VDJ_order1)
colnames(VDJ_order1)=c("orig.ident","Samplesorder")
for (x in VDJ_order1$orig.ident){
  testrds4@meta.data$Sample_order[testrds4@meta.data$orig.ident==x]=VDJ_order1$Samplesorder[VDJ_order1$orig.ident==x]
}

testrds4@meta.data$rownames=rownames(testrds4@meta.data)

testrds4@meta.data$rownames=paste(str_replace_all(testrds4@meta.data$rownames,"(.+_)|(-.+)",""),testrds4@meta.data$Sample_order,sep = "-")


testrds4 <- RenameCells(testrds4, new.names = testrds4@meta.data$rownames)


all_barcode=read.csv("VDJ_View_AllBarcode.csv")

testrds5=subset(testrds4,subset =  rownames %in%  all_barcode$Barcode)



##remove blank value
#testrds2@meta.data$d_gene[ is.na(testrds2@meta.data$d_gene)]="NA"
dir_files="VDJ_vloupe/"


create_loupe_from_seurat(
  obj = testrds5,
  output_dir = dir_files,  # Optional: specify output directory
  output_name = "All_samples_RNA_VDJ_annoatedTNK_changeID_simplecolumn",          # Optional: specify file name (without .cloupe extension)
  metadata_cols = c("seurat_clusters","annotation_summary","state",'orig.ident'), ### must have this parametere, otherwise have error
  #feature_ids = feature_ids,               # Optional: feature IDs from step 3
  dedup_clusters = FALSE,                  # Set to TRUE to deduplicate cluster names
  force = FALSE                            # Set to TRUE to overwrite existing files
)

