##### Single cell analysis ####

##### Import libraries #####
library(rliger)
library(liger)
library(Seurat)
library(Matrix)
library(SeuratData)
library(SeuratWrappers)
library(dplyr)
library(plyr)
library(ggplot2)
library(quanteda)
library(cowplot)
library(edgeR)

###snRNAseq analysis from multi-omics dataset
###We extracted count data with july option 
###scRNAseq data
iNSC_P<-Read10X("./single-cell-intron/iNSC-P-intron/outs/filtered_feature_bc_matrix/") #DR5a
iNSC_S<-Read10X("./single-cell-intron/iNSC-S-intron/outs/filtered_feature_bc_matrix/") #103
iNSC_P0<-Read10X("./single-cell-intron/IB4-111-intron/outs/filtered_feature_bc_matrix/")
iNSC_S0<-Read10X("./single-cell-intron/IB4-113-intron/outs/filtered_feature_bc_matrix/")
###snRNAseq data batch1
iNSC_P1<-Read10X("./iNSC-scRNA-scATAC/example2/outs/filtered_feature_bc_matrix/") #DR5a-Healthy
iNSC_S1<-Read10X("./iNSC-scRNA-scATAC/example5/outs/filtered_feature_bc_matrix/") # 103-PMS
iNSC_S2<-Read10X("./iNSC-scRNA-scATAC/example4/outs/filtered_feature_bc_matrix/") #113-PMS
###snRNAseq data batch2
iNSC_P2<-Read10X("./iNSC-scRNA-scATAC/IB6-111-intron/outs/filtered_feature_bc_matrix/") #111-Healthy
iNSC_S3<-Read10X("./iNSC-scRNA-scATAC/IB6-102-intron/outs/filtered_feature_bc_matrix/") #102-PMS2

# Basic filterring process
iNSC_P <- CreateSeuratObject(counts = iNSC_P, project = "iNSC",min.features=1000,min.cells=3)   
iNSC_S <- CreateSeuratObject(counts = iNSC_S, project = "iNSC",min.features=1000,min.cells=3)   
iNSC_P0 <- CreateSeuratObject(counts = iNSC_P0, project = "iNSC",min.features=1000,min.cells=3)   
iNSC_S0 <- CreateSeuratObject(counts = iNSC_S0, project = "iNSC",min.features=1000,min.cells=3)   
iNSC_P1 <- CreateSeuratObject(counts = iNSC_P1[["Gene Expression"]], project = "iNSC",min.features=1000,min.cells=3)   
iNSC_S1 <- CreateSeuratObject(counts = iNSC_S1[["Gene Expression"]], project = "iNSC",min.features=1000,min.cells=3)   
iNSC_S2 <- CreateSeuratObject(counts = iNSC_S2[["Gene Expression"]], project = "iNSC",min.features=1000,min.cells=3)
iNSC_P2 <- CreateSeuratObject(counts = iNSC_P2, project = "iNSC",min.features=1000,min.cells=3)   
iNSC_S3 <- CreateSeuratObject(counts = iNSC_S3, project = "iNSC",min.features=1000,min.cells=3)   

# Meta data
iNSC_P@meta.data[,"Treatment"]<-"Healthy"
iNSC_S@meta.data[,"Treatment"]<-"MS"
iNSC_P0@meta.data[,"Treatment"]<-"Healthy"
iNSC_S0@meta.data[,"Treatment"]<-"MS"
iNSC_P1@meta.data[,"Treatment"]<-"Healthy"
iNSC_S1@meta.data[,"Treatment"]<-"MS"
iNSC_P2@meta.data[,"Treatment"]<-"Healthy"
iNSC_S2@meta.data[,"Treatment"]<-"MS"
iNSC_S3@meta.data[,"Treatment"]<-"MS"

iNSC_P@meta.data[,"Replicate"]<-"1"
iNSC_S@meta.data[,"Replicate"]<-"1"
iNSC_P0@meta.data[,"Replicate"]<-"1"
iNSC_S0@meta.data[,"Replicate"]<-"1"
iNSC_P1@meta.data[,"Replicate"]<-"1"
iNSC_S1@meta.data[,"Replicate"]<-"1"
iNSC_P2@meta.data[,"Replicate"]<-"1"
iNSC_S2@meta.data[,"Replicate"]<-"1"
iNSC_S3@meta.data[,"Replicate"]<-"1"

iNSC_P@meta.data[,"ID"]<-"scDR5a"
iNSC_S@meta.data[,"ID"]<-"sc103_8"
iNSC_P0@meta.data[,"ID"]<-"sc111_v1"
iNSC_S0@meta.data[,"ID"]<-"sc113_7g"
iNSC_P1@meta.data[,"ID"]<-"snDR5a"
iNSC_S1@meta.data[,"ID"]<-"sn103_8"
iNSC_P2@meta.data[,"ID"]<-"sn111_v2"
iNSC_S2@meta.data[,"ID"]<-"sn113_7g"
iNSC_S3@meta.data[,"ID"]<-"sn102_3"

iNSC_P@meta.data[,"Method"]<-"scRNAseq"
iNSC_S@meta.data[,"Method"]<-"scRNAseq"
iNSC_P0@meta.data[,"Method"]<-"scRNAseq"
iNSC_S0@meta.data[,"Method"]<-"scRNAseq"
iNSC_P1@meta.data[,"Method"]<-"snRNAseq"
iNSC_S1@meta.data[,"Method"]<-"snRNAseq"
iNSC_P2@meta.data[,"Method"]<-"snRNAseq"
iNSC_S2@meta.data[,"Method"]<-"snRNAseq"
iNSC_S3@meta.data[,"Method"]<-"snRNAseq"

iNSC_P@meta.data[,"ID2"]<-"DR5a"
iNSC_S@meta.data[,"ID2"]<-"103_8"
iNSC_P0@meta.data[,"ID2"]<-"111"
iNSC_S0@meta.data[,"ID2"]<-"113_7g"
iNSC_P1@meta.data[,"ID2"]<-"DR5a"
iNSC_S1@meta.data[,"ID2"]<-"103_8"
iNSC_P2@meta.data[,"ID2"]<-"111"
iNSC_S2@meta.data[,"ID2"]<-"113_7g"
iNSC_S3@meta.data[,"ID2"]<-"102_3"

WD<-merge(x=iNSC_P,y=list(iNSC_S,iNSC_P0,iNSC_S0,iNSC_P1,iNSC_S1,iNSC_P2,iNSC_S2,iNSC_S3), add.cell.ids=c("iNSC_P","iNSC_S","iNSC_P0","iNSC_S0","iNSC_P1","iNSC_S1","iNSC_P2","iNSC_S2","iNSC_S3"),project="BRAIN")
###Normalize and find variable features

#QC clusters
FeatureScatter(WD, "nCount_RNA", "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0.5)
FeatureScatter(WD, "nCount_RNA", "nFeature_RNA", group.by = "ID", pt.size = 0.5)

###Switch to the integrated data for downstream analyses
DefaultAssay(object = WD) <- "integrated"
DefaultAssay(object = WD) <- "RNA"

WD <- NormalizeData(WD, normalization.method = "LogNormalize", scale.factor = 10000)
WD <- FindVariableFeatures(WD, selection.method = "vst", nfeatures = 2000)
WD <- ScaleData(object = WD, verbose = FALSE)

VlnPlot(WD, "nFeature_RNA", group.by = "ID", pt.size=0)
VlnPlot(WD, "mitoRatio", group.by = "ID", pt.size=0.1)

WD <- RunPCA(object = WD, npcs = 30, verbose = FALSE)
WD <- FindNeighbors(object = WD, dims = 1:20)
WD <- FindClusters(object = WD)
WD <- RunUMAP(object = WD, reduction = "pca", 
              dims = 1:30)
p1 <- DimPlot(object = WD, reduction = "umap", group.by = "ID")
p2 <- DimPlot(object = WD, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2)

summary(factor(WD$Treatment))

View(WD@meta.data)
WD$log10GenesPerUMI <- log10(WD$nFeature_RNA) / log10(WD$nCount_RNA)

# Compute percent mito ratio
WD$mitoRatio <- PercentageFeatureSet(object = WD, pattern = "^MT-")
WD$mitoRatio <- WD@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- WD@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Create sample column
View(metadata)

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=Treatment, fill=Treatment)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

library(tidyverse)

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=Treatment, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Add metadata back to Seurat object
WD@meta.data <- metadata

write.csv(metadata,file="metadata-integrated-2022-scRNAseq.csv")
saveRDS(WD,file="metadata-integrated-2022-scRNAseq-raw.rds")
WD <- readRDS(file="metadata-integrated-2022-scRNAseq-raw.rds")

# Filter out low quality reads using selected thresholds - these will change with experiment
# Python package SCVI (generating the doublet ranking top to bottom, exclude only 5% top doublet candidate cells
metadata_before_filter <- read.table(file="metadata-integrated-2022-scRNAseq.csv-doublet.csv",sep=",",header=TRUE)
rownames(metadata_before_filter) <- metadata_before_filter$name
WD@meta.data <- metadata_before_filter

# stringent filter test #1
filtered_seurat <- subset(x = WD, 
                          subset= (nUMI >= 10000) & 
                            (nGene >= 1000) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# relaxed filter 10/31/2022, current version (to increase # of downstream cells)
filtered_seurat <- subset(x = WD, 
                          subset= (nUMI >= 8000) & 
                            (nGene >= 1000) & 
                            (log10GenesPerUMI > 0.75) & 
                            (mitoRatio < 0.25))

#relaxed result
#An object of class Seurat 
#30525 features across 32446 samples within 1 assay 
#Active assay: RNA (30525 features, 2000 variable features)
#2 dimensional reductions calculated: pca, umap
# Python package SCVI (generating the doublet ranking top to bottom, exclude only 5% top doublet candidate cells
metadata <- read.table(file="metadata-integrated-harmony2-doublet2.csv",sep=",",header=TRUE)
rownames(metadata) <- metadata$name

filtered_seurat@meta.data <- metadata

Idents(object=filtered_seurat)<-"doublet"
#SetIdent(pbmc, value="doublet")
filtered_seurat <- subset(x = filtered_seurat, idents = "singlet")
#An object of class Seurat 
#30525 features across 30227 samples within 1 assay 
#Active assay: RNA (30525 features, 2000 variable features)
#2 dimensional reductions calculated: pca, umap

Idents(object=filtered_seurat)<-"seurat_clusters"
DimPlot(object = filtered_seurat, reduction = "umap", group.by = "seurat_clusters", 
              label = TRUE, repel = TRUE)

WD <- filtered_seurat
# 10/31/2022 updated
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. 
# RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, 
# so we can make sure that the Harmony objective function gets better with each round.
library(cowplot)
library(harmony)
options(repr.plot.height = 2.5, repr.plot.width = 6)
integrated <- WD %>% 
  RunHarmony("ID", plot_convergence = TRUE)

# To directly access the new Harmony embeddings, use the Embeddings command.
harmony_embeddings <- Embeddings(integrated, 'harmony')
harmony_embeddings[1:5, 1:5]

# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = integrated, reduction = "harmony", pt.size = .1, group.by = "ID")
p2 <- VlnPlot(object = integrated, features = "harmony_1", group.by = "ID", pt.size = .0)
plot_grid(p1,p2)

integrated <- integrated %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(integrated, reduction = "umap", group.by = "ID2", pt.size = .5, split.by = 'ID2')

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(integrated, reduction = "umap", group.by = "Treatment", pt.size = .5, split.by = 'Treatment')

options(repr.plot.height = 4, repr.plot.width = 6)
p1 <- DimPlot(object = pintegrated, reduction = "umap", pt.size = .5, group.by = "ID")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = .5)
plot_grid(p1,p2)

# adding UMAP information
metadata <- integrated@meta.data
metadata2 <- data.frame(integrated[["umap"]]@cell.embeddings)
integrated@meta.data <- cbind(metadata, metadata2)

# 10/31/2022 updated/ removal of doublet cells (<10%)
# Ver2 second relaxed version
saveRDS(integrated_singlet,file="metadata-integrated-2022-scRNAseq-harmony-singlet.rds")
write.csv(iintegrated_singlet@meta.data, file="metadata-integrated-2022-scRNAseq-harmony-singlet.csv")
integrated_singlet <- readRDS(file="metadata-integrated-2022-scRNAseq-harmony-singlet.rds")

###################################################################3
#Downstream visualisation doublt/singlet test old version
###################################################################3
# this error fixed after importing tidyverse library
# Error: Must request at least one colour from a hue palette.
library(tidyverse)
library(ggplot2)

DimPlot(integrated, reduction = "umap", group.by = "ID2", pt.size = .5, split.by = 'ID2')
DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", pt.size = .5, split.by ="seurat_clusters")
DimPlot(integrated, reduction = "umap", group.by = "doublet", pt.size = .1, split.by = 'doublet')

Idents(object=integrated)
levels(x = integrated)

###Visualize markers by cluster
Idents(object=integrated)<-"doublet"
Idents(object=integrated)
levels(x = integrated)

integrated_singlet <- subset(x = integrated, idents = "singlet")
integrated_singlet@meta.data <- read.table(file="metadata-integrated-harmony2-singlet2.csv",sep=",",header=TRUE)

Idents(object=integrated_singlet)<-"seurat_clusters"

DimPlot(integrated, reduction = "umap", group.by = "ID2", pt.size = .5, split.by = 'ID2')
DimPlot(integrated_singlet, reduction = "umap", group.by = "doublet", pt.size = .1, split.by = 'doublet')
DimPlot(integrated, reduction = "umap", label = TRUE)

options(repr.plot.height = 4, repr.plot.width = 6)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "ID2", label = TRUE, pt.size = .5)
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "ID2", label = TRUE, pt.size = .5)
plot_grid(p1,p2)

###################################################################3
#Downstream visualisation doublt/singlet test old version
###################################################################3
# Very nice feature RidgePlot generation key genes
features <- c( "NES", "SOX2", "HES6", "YBX1","PAX6", "DACH1","IFIT1","IFIT3")
RidgePlot(integrated_singlet, features = features, ncol = 3)

library(shiny)
library(Seurat)
library(ShinyCell)

integrated_singlet[["umap"]]@cell.embeddings

#DefaultAssay(pbmc) <- 'RNA'
Idents(object=integrated_singlet)<-"seurat_clusters"
integrated_singlet <- subset(x = integrated_singlet, idents = c("0","1","2","3","4","5","6","7","8","9"))

# excluding cluster 10,11 small cluster, specific to the one sample 
scConf = createConfig(integrated_singlet)
makeShinyApp(integrated_singlet, scConf, gene.mapping = TRUE,
             shiny.title = "iNSC (MS vs Healthy)")

WD.brain <- readRDS(file="metadata-integrated-2022-scRNAseq-harmony-singlet.rds")

###Normalize and find variable features
table(WD.brain@meta.data$ID,WD.brain@meta.data$seurat_clusters) 
t(round(prop.table(table(WD.brain@meta.data$Treatment,WD.brain@meta.data$seurat_clusters),margin=1),4))
cbind(t(round(prop.table(table(WD.brain@meta.data$Treatment,WD.brain@meta.data$seurat_clusters),margin=1),4)))
test<-chisq.test(table(WD.brain@meta.data$seurat_clusters==0, WD.brain@meta.data$Treatment))
str(test)
chi.res<-data.frame(mat.or.vec(17,3))
names(chi.res)<-c("cluster","chi_seq_stat","p-value")
for (i in 0:17){
  test<-chisq.test(table(WD.brain@meta.data$seurat_clusters==i, WD.brain@meta.data$Treatment))
  chi.res[i+1,1]<-i
  chi.res[i+1,2]<-round(test$statistic,3)
  chi.res[i+1,3]<-test$p.value
}
cbind(t(round(prop.table(table(WD.brain@meta.data$Treatment,WD.brain@meta.data$seurat_clusters),margin=1),4)),chi.res$'p-value')
a <- cbind(t(round(prop.table(table(WD.brain@meta.data$Treatment,WD.brain@meta.data$seurat_clusters),margin=1),4)),chi.res$'p-value')
write.csv(a, file="cell-proportion.txt")
a <- structure(list(`X$Days` = c("10", "38", "66", "101", "129", "185", "283", 
                                 "374")), .Names = "X$Days")

#this doesn't work need to fix it - bargarph
library(gcookbook) # Load gcookbook for the cabbage_exp data set
ggplot(b, aes(x = colnames(b), y = colnames(b), fill = Cultivar)) +
  geom_col(position = "fill")

###Switch to the integrated data for downstream analyses
DefaultAssay(object = WD) <- "integrated"
DefaultAssay(object = WD.brain) <- "RNA"

par(mfrow=c(2,2))
VlnPlot(WD.brain,features=c("SOX2","NES"),pt.size=0)
VlnPlot(WD.brain,features=c("HMGB1","HMGB2"),pt.size=0)
VlnPlot(WD.brain,features=c("CDKN1A","CDKN1C"),pt.size=0)
VlnPlot(WD.brain,features=c("EBF1","DACH1"),pt.size=0)
VlnPlot(WD.brain,features=c("IFIT1","IFIT3"),pt.size=0)

###Visualize markers by cluster
WD.brain <- pbmc # when readRDS from file
Idents(object=WD.brain)<-"seurat_clusters"
DefaultAssay(object = WD.brain) <- "RNA"
WD.brain <- NormalizeData(WD.brain)

DotPlot(WD.brain,col.min=0,features=c("MSI1", "VIM", "ASCL1","SLC1A3", "OLIG2"))+labs(title="")

DotPlot(WD.brain,col.min=0,features=c("SOX2", "NES", "DACH1","PLZF", "ZO1", "PAX6", "AP2A","HNK1","SOX10","HES6"))+labs(title="")

DotPlot(WD.brain,col.min=0,features=c("ANGPT2",
                                      "APOE",
                                      "B2M",
                                      "C1orf54",
                                      "C1QBP",
                                      "CHRDL1",
                                      "CLU",
                                      "COL11A1",
                                      "COL15A1",
                                      "CRTAP",
                                      "CST3",
                                      "EDA",
                                      "FGF2",
                                      "FN1",
                                      "FST",
                                      "FSTL1",
                                      "FUCA2",
                                      "IGFBP2",
                                      "IGFBP5",
                                      "IGFBP6",
                                      "IGSF10",
                                      "ISG15",
                                      "LGALS1",
                                      "LUZP2",
                                      "MANF",
                                      "MIF",
                                      "NPC2",
                                      "NRG2",
                                      "OLFML3",
                                      "PPIB",
                                      "PRSS23",
                                      "RSPO2",
                                      "SEMA3D",
                                      "SEMA3E",
                                      "ST6GAL1",
                                      "TCN2",
                                      "TIMP1",
                                      "TIMP3",
                                      "VASH2"))

# 3/14/2023 Alex's request
FeaturePlot(pbmc,  features = c("YBX1", "DTYMK", "SOX11"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

# 11/21/2022 Isabel's request
FeaturePlot(WD.brain,  features = c("MSI1", "VIM", "ASCL1","SLC1A3", "OLIG2"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("SOX2", "NES", "DACH1","PLZF", "ZO1", "PAX6", "AP2A","HNK1","PROM1","SOX11","HES6","YBX1","ETNPPL"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("HMGB1","HMGB2","SOX2", "NES", "DACH1","PLZF", "ZO1", "PAX6", "AP2A","HNK1","PROM1","SOX10","CDKN1A","HES6","YBX1","ETNPPL"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("HMGB1","HMGB2","SOX2", "NES", "DACH1","PLZF", "ZO1", "PAX6","IFIT1","IFIT3","SH2D4B"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("DACH1","ZBTB16", "ZBT16", "TJP1", "PLZF", "ZO1", "TFAP2A","B3GAT1","SOX10","CDKN1C"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("ZBTB16", "ZBT16", "TJP1", "AP2A","HNK1","SOX10"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("B2M",
                                    "TIMP1",
                                    "EDA",
                                    "FN1",
                                    "FGF2",
                                    "MIF",
                                    "IGS15"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("FGF2",
                                    "IGFBP5",
                                    "MIF",
                                    "TUBGCP2"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("FN1","IGFBP5",
                                    "FUCA2",
                                    "FSTL1","CHRDL1","CST3"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c(
                                    "IGFBP5",
                                    "NADK",
                                    "POU5F1",
                                    "TXNIP"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")


FeaturePlot(WD.brain,  features = c("IFIT1","IFIT3","EBF1","NES"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            reduction = "umap")

FeaturePlot(WD.brain,  features = c("SOX2","NES","DACH1","PAX6","PROM1","YBX1","SOX1","NUMB","BMI1","PSEN1"), 
            cols =c("lightgrey", "red"),
            min.cutoff = "q9", 
            keep.scale = "all",
            reduction = "umap")

WD.brain.markers <- FindAllMarkers(WD.brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WD.brain.markers.all <- FindAllMarkers(WD.brain, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(WD.brain.markers,file="NIH-NIA-scRNAseq-and-snRNAseq.markers.rds")
WD.brain.markers = readRDS(file="NIH-NIA-scRNAseq-and-snRNAseq.markers.rds")

# R version4 didn't work
WD.brain.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- WD.brain.markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
DoHeatmap(WD.brain, features = top10$gene) 

gene_of_interest = c("SOX2", "NES", "DACH1","PLZF", "ZO1", "PAX6", "AP2A","HNK1","SOX10","HES6")
DoHeatmap(WD.brain, features = gene_of_interest)

VariableFeatures(WD) <- var.genes
WD.brain <- ScaleData(WD.brain, features=VariableFeatures(WD.brain))

write.csv(WD.brain.markers,file="NIH-NIA-scRNAseq-and-snRNAseq.markers.csv")
write.csv(WD.brain.markers,file="NIH-NIA-scRNAseq-and-snRNAseq.markers-singlet.csv")
write.csv(WD.brain.markers.all,file="NIH-NIA-scRNAseq-and-snRNAseq.markers-singlet-all.csv")

##############3 Velocity Analysis Preparation ##################

#Filtered
#> summary(factor(filtered_seurat$Treatment))
#Old Young 
#19548 12987 

WD<-filtered_seurat
WD$barcode <- colnames(WD)
WD$UMAP_1 <- WD@reductions$umap@cell.embeddings[,1]
WD$UMAP_2 <- WD@reductions$umap@cell.embeddings[,2]
WD$seurat_clusters_name <- paste("Cluster",WD$seurat_clusters,sep="_")
write.csv(WD@meta.data, file='metadata_brain.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(WD, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts_brain.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(WD@reductions$pca@cell.embeddings, file='pca_brain.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names_brain.csv',
  quote=F,row.names=F,col.names=F
)

# Finding sub-groups here MS vs Healthy
DefaultAssay(pbmc) <- 'RNA'
Idents(object=pbmc)<-"Treatment"
healthy <- subset(x = pbmc, idents = c("Healthy"))
msdrived <- subset(x = pbmc, idents = c("MS"))
Idents(object=pbmc)<-"ID2"
msdrived2 <- subset(x = pbmc, idents = c("103_8","113_7g"))

#healthy <- obj.list[1]
#msdrived <- obj.list[2]

#healthy = as.Seurat(healthy)
#msdrived = as.Seurat(msdrived)

healthy$barcode <- colnames(healthy)
healthy$UMAP_1 <- healthy@reductions$umap@cell.embeddings[,1]
healthy$UMAP_2 <- healthy@reductions$umap@cell.embeddings[,2]
healthy$seurat_clusters_name <- paste("Cluster",healthy$seurat_clusters,sep="_")
write.csv(healthy@meta.data, file='metadata_healthy.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(healthy, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts_healthy.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(healthy@reductions$pca@cell.embeddings, file='pca_healthy.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names_healthy.csv',
  quote=F,row.names=F,col.names=F
)

#MS Derived iNSC -processing
msdrived$barcode <- colnames(msdrived)
msdrived$UMAP_1 <- msdrived@reductions$umap@cell.embeddings[,1]
msdrived$UMAP_2 <- msdrived@reductions$umap@cell.embeddings[,2]
msdrived$seurat_clusters_name <- paste("Cluster",msdrived$seurat_clusters,sep="_")
write.csv(msdrived@meta.data, file='metadata_msdrived.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(msdrived, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts_msdrived.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(msdrived@reductions$pca@cell.embeddings, file='pca_msdrived.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names_msdrived.csv',
  quote=F,row.names=F,col.names=F
)

#MS Derived iNSC -processing ver2
msdrived$barcode <- colnames(msdrived)
msdrived$UMAP_1 <- msdrived@reductions$umap@cell.embeddings[,1]
msdrived$UMAP_2 <- msdrived@reductions$umap@cell.embeddings[,2]
msdrived$seurat_clusters_name <- paste("Cluster",msdrived$seurat_clusters,sep="_")
write.csv(msdrived@meta.data, file='metadata_msdrived2.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(msdrived, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('./', 'counts_msdrived2.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(msdrived@reductions$pca@cell.embeddings, file='pca_msdrived2.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names_msdrived2.csv',
  quote=F,row.names=F,col.names=F
)

# Testing SCP library 11/1/2022
# Re-install SCP library 1/24/2023
library(SCP)

# Finding sub-groups here MS vs Healthy
DefaultAssay(pbmc) <- 'RNA'
Idents(object=pbmc)<-"seurat_clusters"
pbmc_cluster8 <- subset(x = pbmc, idents = c("0","1","2","3","4","8"))
pbmc_cluster8 <- subset(x = pbmc, idents = c("3"))
pbmc_cluster8 <- Standard_SCP(srt = pbmc_cluster8, nonliner_reduction="umap")

# Very nice feature RidgePlot generation key genes
stem_features = c(
  "SOX2","NES", "DACH1","PAX6", "PROM1",# Pre-endocrine
  "YBX1", 
  "SOX1","NUMB", # Ductal
  "BMI1","PSEN1" # EPs
  # Beta, Alpha, Delta, Epsilon
)

glial_features = c(
  "ID4","SOX9", "SOX10","FGF2", "HES1",
  "CD44","FABP7","SLC1A3","CLU", 
  "LGALS3","PPM1K", # Ductal
  "DTYMK","TPX2","KLF15"
)

neural_features = c("SOX11","HES6","ASCL1","ID2","DCX")

features <- c( "NES", "SOX2", "HES6", "YBX1","PAX6", "DACH1","IFIT1","IFIT3")
RidgePlot(pbmc_cluster8, features = stem_features, ncol = 3)
DotPlot(pbmc_cluster8,col.min=0,features=stem_features)+labs(title="")


ClassDimPlot(
  srt = pancreas1k, group.by = c("CellType", "SubCellType"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
)

ClassDimPlot(
  srt = pbmc_cluster8, group.by = c("ID2", "Treatment"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
)

ClassDimPlot(
  srt = pbmc, group.by = c("Treatment","ID2"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank", palcolor=c("#27AE60","orange","#5DADE2","red")
)

ClassDimPlot(
  srt = pbmc_cluster8, group.by = c("Treatment","ID2"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank", palcolor=c("#27AE60","orange","#5DADE2","red")
)

ClassDimPlot(
  srt = pbmc_cluster8, group.by = c("Treatment"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
)

# it doesn't work 1/24/2023
GroupHeatmap(
  srt = pbmc_cluster8,
  features = c(
    "RYR2", "IFIT2", # Ductal
    "OASL","IFIT1", "IFIT3", # EPs
    "HERC5", "CDKN2A", # Pre-endocrine
    "KCNIP4", "NEAT1", # Endocrine
    "NES", "SOX2", "DACH1", "EBF1" # Beta, Alpha, Delta, Epsilon
  ),
  group.by = c("Treatment", "seurat_clusters"),
  heatmap_palette = "YlOrRd",
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)

# TF analysis
GroupHeatmap(
  srt = pbmc_cluster8,
  features = c(
    "SOX2","IRF3", "STAT2","BHLHE22", "SOX10",# Pre-endocrine
    "KLF6", 
    "SOX4","RFX3", # Ductal
    "SOX8","MEIS1", "PDX1" # EPs
     # Beta, Alpha, Delta, Epsilon
  ),
  group.by = c("Treatment", "seurat_clusters"),
  heatmap_palette = "YlOrRd",
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)

par(mfrow=c(2,2))
VlnPlot(pbmc_cluster8,features=c("SOX2","NES","DACH1","PAX6","YBX1","PROM1"),pt.size=0)

# Stemness analysis
ht <- GroupHeatmap(
  srt = pbmc_cluster8,
  features = c("SOX11","HES6","ASCL1","DCX"),
  group.by = c("seurat_clusters"),
  heatmap_palette = "YlOrRd",
  exp_method = c("log2fc"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = FALSE, add_reticle = FALSE
)
print(ht$plot)

# Stemness analysis
ht <- GroupHeatmap(
  srt = pbmc_cluster8,
  features = c(stem_features,glial_features),
  group.by = c("seurat_clusters"),
  heatmap_palette = "YlOrRd",
  exp_method = c("log2fc"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = FALSE, add_reticle = FALSE
)
print(ht$plot)

ht <- FeatureHeatmap(
  srt = pbmc_cluster8, group.by = "seurat_clusters", features = stem_features, 
  height = 5, width = 4
)
print(ht$plot)

# DEG+GO terms
pbmc_cluster8 <- RunDEtest(srt = pbmc_cluster8, group_by = "Treatment", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = pbmc_cluster8, group_by = "Treatment")
DEGs <- pbmc_cluster8@tools$DEtest_Treatment$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "SP"))
ht <- ExpHeatmap(
  srt = pancreas_sub, group.by = "CellType", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Mus_musculus", db = c("GO_BP", "KEGG", "WikiPathway"), anno_terms = TRUE,
  feature_annotation = c("TF", "SP"), feature_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)

###Looks like there's significant differences in cluster membership by treatment
# SCP tutorial scVELO test 11/1/2022
pbmc <- RunSCVELO(
  srt = pbmc_cluster3, group_by = "seurat_clusters",
  liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE
)
VelocityPlot(srt = pbmc, reduction = "UMAP", group_by = "seurat_clusters")

VelocityPlot(srt = pbmc, reduction = "UMAP", plot_type = "stream")


library(data.table)
data_to_write_out <- as.data.frame(as.matrix(Breast_1@scale.data))
fwrite(x = data_to_write_out, file = "breast_1.outfile.csv")

# Extract subset (for example cluster 7) from the analyssi
cells.use <- WhichCells(object = WD, ident = 0)
expr <- GetAssayData(object = WD, assay.type = "RNA", slot = "counts")[, cells.use] # slot = "data" normalized one
expr <- as(Class = 'matrix', object = expr)
write.csv(x = expr, file = "./sub_clusters/sub_cluster_0_v3/expression_merged_cluster0.csv", quote = FALSE)

for (i in 1:10) {
  cells.use <- WhichCells(object = WD, ident = i)
  expr <- GetAssayData(object = WD, assay.type = "RNA", slot = "counts")[, cells.use] # slot = "data" normalized one
  expr <- as(Class = 'matrix', object = expr)
  write.csv(x = expr, file = paste(paste("./sub_clusters/expression_merged_cluster_v3_",i,sep=""),"csv",sep="."), quote = FALSE)
}

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)

WD <-pbmc

DimPlot(WD, label = TRUE) + NoLegend()
erythroid <- WD[, WD$seurat_clusters %in% c(0,1,2,3,4)]
lymphoid <- WD[, WD$seurat_clusters %in% c(6,7)]
all <- WD[, WD$seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9)]

erythroid.cds <- as.cell_data_set(erythroid)
erythroid.cds <- cluster_cells(cds = erythroid.cds, reduction_method = "UMAP")
erythroid.cds <- learn_graph(erythroid.cds, use_partition = TRUE)

lymphoid.cds <- as.cell_data_set(lymphoid)
lymphoid.cds <- cluster_cells(cds = lymphoid.cds, reduction_method = "UMAP")
lymphoid.cds <- learn_graph(lymphoid.cds, use_partition = TRUE)

cluster_cells(all.cds, reduction_method="UMAP", cluster_method="leiden")
all.cds <- as.cell_data_set(all)
all.cds <- cluster_cells(cds = all.cds, reduction_method = "UMAP", k=28)
all.cds <- learn_graph(all.cds, use_partition = TRUE)

cells.use <- WhichCells(object = WD, ident = c(0,1,2,3))
cells.use <- WhichCells(object = WD, ident = c(3))

order_cells(erythroid.cds, reduction_method="UMAP", root_cells=cells.use)
order_cells(lymphoid.cds, reduction_method="UMAP", root_cells=cells.use)
order_cells(all.cds, reduction_method="UMAP", root_cells=cells.use)
# plot trajectories colored by pseudotime
plot_cells(
  cds = all.cds,
  color_cells_by = "cluster",
  reduction_method = c("UMAP"),
  show_trajectory_graph = TRUE,
  cell_size = 0.5,
  graph_label_size=5,
)

ciliated_genes <- c("AL627309.1")

plot_cells(all.cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

plot_cells_3d(
  cds = all.cds,
  color_cells_by = "cluster",
  reduction_method = c("UMAP"),
  show_trajectory_graph = TRUE
)

# plot trajectories colored by pseudotime
plot_cells(
  cds = all.cds,
  color_cells_by = "seurat_clusters",
  reduction_method = c("UMAP"),
  show_trajectory_graph = TRUE,
  graph_label_size=3.5,
  cell_size = 0.5
)

all.cds <- order_cells(all.cds, reduction_method="UMAP", root_cells=cells.use)

plot_cells(
  cds = all.cds,
  color_cells_by = "pseudotime",
  label_cell_groups=FALSE,
  label_leaves=FALSE,
  label_branch_points=FALSE,
  graph_label_size=0.0,
  cell_size=0.5
) + scale_color_viridis_c(option = "plasma",begin = 0, end = 0.6) 


