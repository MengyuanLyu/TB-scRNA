library(Seurat)
library(celldex)
library(SingleR)
library(ggplot2)
library(paletteer) 
library(metap)
library(multtest)
library(dplyr)

#read rawdata
TB <- readRDS("TB_W_Ab_ST_Samplesonly.rds")
HC <- readRDS("HC_W_Ab_ST_Samplesonly.rds")

#quality control
HC[["percent.mt"]] <- PercentageFeatureSet(HC, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(HC@assays$RNA)) 
HB.genes <- rownames(HC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
HC[["percent.HB"]]<-PercentageFeatureSet(HC, features=HB.genes) 
minGene=500
maxGene=4000
pctMT=15
HC <- subset(HC, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
TB <- subset(TB, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

#combine
scRNA<-list(HC,TB)
names(scRNA)<-c("HC","TB")
scRNA <- lapply(X = scRNA, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = scRNA)
immune.anchors <- FindIntegrationAnchors(object.list = scRNA, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

#clustering
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE) 
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)

#cell annonation
allmarkers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct = 0.25, 
                             thresh.use = 0.25)
DimPlot(immune.combined, group.by="celltype", label=F, label.size=4, reduction='umap',cols=col)+labs(x = "UMAP_1", y = "UMAP_2")+ theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 

#Marker genes
FeaturePlot(immune.combined, features = c("STMN1"), min.cutoff = "q9")

#differentially expressed genes
NK <- FindMarkers(immune.combined, ident.1 = "NK cells_HC", ident.2 = "NK cells_TB", verbose = FALSE)
NK <- subset(NK, p_val_adj < 0.05)

#re-reduction for T cells
Tcells <- subset(immune.combined@meta.data,celltype=='T cells')
Tcells  <- subset(immune.combined, cells=row.names(Tcells))
DefaultAssay(Tcells) <- "integrated"
Tcells <- ScaleData(Tcells,  verbose = FALSE)
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)
Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:30)
Tcells <- FindNeighbors(Tcells, reduction = "pca", dims = 1:30) 
Tcells <- FindClusters(Tcells, resolution = 1)
DimPlot(Tcells, reduction = "umap") 

#primery functional enrichment
T_cells <- irGSEA.score(object = T_cells, assay = "RNA",
                        slot = "data", seeds = 123, ncores = 1,min.cells = 3, min.feature = 0,
                        custom = F, geneset = NULL, msigdb = T,
                        species = "Homo sapiens", category = "H",  
                        subcategory = NULL, geneid = "symbol",
                        method = c("AUCell", "UCell", "singscore","ssgsea"),
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                        kcdf = 'Gaussian')
Seurat::Assays(T_cells)
result.dge <- irGSEA.integrate(object = T_cells ,
                               group.by = "celltype.group",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore","ssgsea"))
irGSEA.heatmap.plot <- irGSEA.bubble(object = result.dge,
                                     method = "ssgsea",
                                     top = 50,
                                     show.geneset = NULL,
                                     cluster.color =col,
                                     significance.color = c('#D9D9D9','#FBB4AE'),
                                     direction.color = c('#80CDC1','#FFD92F'))

#functional assessment
T_cells <- ScaleData(T_cells,  verbose = FALSE)#给RNA计算scale.data

naive  <- subset(T_cells, cells=row.names(naive))
gene_data <- GetAssayData(naive, slot = 'scale.data')
gene_specific <- gene_data['CCR7',]
result <- NULL 
for(i in 1:nrow(gene_data)){
  
  result <-c(result,cor.test(gene_specific, gene_data[i,],method="spearman")$estimate)
  cat(i,'out of',nrow(gene_data),'\n')
  
}
result<-data.frame (genes = rownames(gene_data),cor=result)

a <- c('CCR7','RPS14','EEF1A1','RPS28','HLA-E','RPL29','RPL13','RPS18','RPS19','RPS7')
gene_data1 <- GetAssayData(T_cells, slot = 'scale.data')
gene_specific <- gene_data1[a,]
y<-as.character(Idents(T_cells))
y[as.character(Idents(T_cells))=="CD4+ naive T cells"|as.character(Idents(T_cells))=="CD8+ naive T cells"]<-1
y[as.character(Idents(T_cells))=="CD4+ T helper cells"|as.character(Idents(T_cells))=="CD8+ effector T cells"
  |as.character(Idents(T_cells))=="CD8+ effector memory T cells"|as.character(Idents(T_cells))=="Natural killer T cells"
  |as.character(Idents(T_cells))=="MAIT cells"|as.character(Idents(T_cells))=="CD4+ effector T cells"
  |as.character(Idents(T_cells))=="CD4+ regulatory T cells"]<-0
y <- as.numeric(y)
gene_specificx<-data.frame(cbind(y,gene_specific))
naive <-glm(y~., data=gene_specificx,family = "binomial")

#The same code was used to process the data of myeloid cells, B cells and NK cells. 

