### libraries

library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scCustomize)
library(ggalluvial)
library(ggrepel)
library(RColorBrewer)
library(reshape2)

### configuration
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

### files
idir12h<-"12h/outs/filtered_feature_bc_matrix"
dir10h<-"10h/outs/filtered_feature_bc_matrix"
dir8h<-"8h/outs/filtered_feature_bc_matrix"
annoFile<-"nv2.func-04.04.23.tsv"
g2mFile<-"candidate_g2m.txt"
sFile<-"candidate_s.txt"
ribFile<-"ribosomal.putative.txt"

### functions
renameFeaturePlot <- function(x, y){
  for (i in 1:length(x)){
    nv2<-x[[i]]$labels$title
    newName<-y %>% filter(geneID==nv2) %>% select(gene_short_name) %>% unlist() %>% unname()
    x[[i]]$labels$title <- newName
  }
  return(x)
}

### load data
raw12h<-Read10X(data.dir=dir12h)
raw10h<-Read10X(data.dir=dir10h)
raw8h<-Read10X(data.dir=dir8h)
annot<-read_tsv(annoFile) %>% select(geneID, nve, gene_short_name, cdsName, TFFam, PFAMs, lncrna, GO)
sGenes<-read_tsv(sFile, col_names=F) %>% select(c(1)) %>% unlist() %>% unname()
g2mGenes<-read_tsv(g2mFile, col_names=F) %>% select(c(1)) %>% unlist() %>% unname()
ribGenes<-read_tsv(ribFile, col_names=F) %>% select(c(1)) %>% unlist() %>% unname()
manuMarkers <- c("NV2.6608", "NV2.15833", "NV2.2150", "NV2.9419", "NV2.8483", "NV2.22508", "NV2.234", "NV2.10891", "NV2.15303", "NV2.10624", "NV2.11441", "NV2.472")

### Create Seurat objects
dataObj<-list(
  data8h = CreateSeuratObject(counts=raw8h, min.cells=3, min.features=200, project="8h"),
  data10h = CreateSeuratObject(counts=raw10h, min.cells=3, min.features=200, project="10h"),
  data12h = CreateSeuratObject(counts=raw12h, min.cells=3, min.features=200, project="12h")
)

### plot QC
# distribution of genes and reads per cell for filtering
rbind(dataObj[[1]]@meta.data, dataObj[[2]]@meta.data,dataObj[[3]]@meta.data) %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) +
    theme_bw() +
    geom_density_2d_filled(aes(group=orig.ident)) +
    facet_wrap(~orig.ident, scales="free") +
    labs(x="UMI per cell", y="Genes per cell")

### Filtering
dataFilt<-list(
  data8h = subset(dataObj[["data8h"]], subset = nFeature_RNA>=1500 & nFeature_RNA <=7000 & nCount_RNA>=4000 & nCount_RNA<=30000),
  data10h = subset(dataObj[["data10h"]], subset = nFeature_RNA>=1500 & nFeature_RNA <=7000 & nCount_RNA>=5000 & nCount_RNA<=30000),
  data12h = subset(dataObj[["data12h"]], subset = nFeature_RNA>=1000 & nFeature_RNA <=4000 & nCount_RNA>=4000 & nCount_RNA<=30000)
)

### function for several steps:
normAndCC <- function(x){
	x <- SCTransform(x)
	x <- RunPCA(x)
	x <- RunUMAP(x, dims=1:30)
	s <- sGenes[sGenes %in% rownames(x)]
	g2m <- g2mGenes[g2mGenes %in% rownames(x)]
	x<- CellCycleScoring(x, s.features=s, g2m.features=g2m)
	x$CCDiff <- x$S.Score - x$G2M.Score
	x$percent.cox1 <- PercentageFeatureSet(x, features="NV2.25931")
	rib<-ribGenes[ribGenes %in% rownames(x)]
	x$rib <- PercentageFeatureSet(x, features=rib)
	x@meta.data$orig.ident <- factor(x@meta.data$orig.ident, levels=c("8h", "10h", "12h"))
	x@meta.data$Phase <- factor(x@meta.data$Phase, levels=c("G1", "S", "G2M"))
	return(x)
}

### normalize
dataNorm <- lapply(dataFilt, function(x) {normAndCC(x)})

### regress cell cycle
dataCCReg <- lapply(dataNorm, function(x){SCTransform(x, vars.to.regress="CCDiff", assay = "RNA")})
dataCCReg <- lapply(dataCCReg, function(x){RunPCA(x, reduction.name="pca.ccregress")})
dataCCReg <- lapply(dataCCReg, function(x){RunUMAP(x, reduction="pca.ccregress", reduction.name="umap.ccregress", dims=1:30)})

### merge and integrate
nv2IntSC<-merge(dataCCReg[[1]], y=c(dataCCReg[[2]], dataCCReg[[3]]), add.cell.ids = c("8h", "10h", "12h"), project="Nv2_dev", merge.data = TRUE, merge.dr=TRUE)
nv2IntSC@meta.data$orig.ident <- factor(nv2IntSC@meta.data$orig.ident, levels=c("8h", "10h", "12h"))
nv2IntSC@meta.data$Phase <- factor(nv2IntSC@meta.data$Phase, levels=c("G1", "S", "G2M"))

### plotting post norm QC
VlnPlot(nv2IntSC, split.by="orig.ident", features=c("nCount_RNA", "nFeature_RNA", "percent.cox1", "rib"))

### plotting cell phase by time point
DimPlot(nv2IntSC, split.by="orig.ident", group.by="Phase", reduction="umap.ccregress", pt.size=2, alpha=0.5)

### plotting molecular markers
# first check the gene names
annot %>% filter(geneID %in% manuMarkers)

# plot markers in paper
FeaturePlot_scCustom(nv2IntSC, split.by="orig.ident", colors_use=rev(brewer.pal(11, "Spectral")), features=manuMarkers[1:3], reduction="umap.ccregress", alpha_exp=0.5, order=TRUE, pt.size=1.5)
FeaturePlot_scCustom(nv2IntSC, split.by="orig.ident", colors_use=rev(brewer.pal(11, "Spectral")), features=manuMarkers[4:6], reduction="umap.ccregress", alpha_exp=0.5, order=TRUE, pt.size=1.5)
FeaturePlot_scCustom(nv2IntSC, split.by="orig.ident", colors_use=rev(brewer.pal(11, "Spectral")), features=manuMarkers[7:9], reduction="umap.ccregress", alpha_exp=0.5, order=TRUE, pt.size=1.5)
FeaturePlot_scCustom(nv2IntSC, split.by="orig.ident", colors_use=rev(brewer.pal(11, "Spectral")), features=manuMarkers[10:12], reduction="umap.ccregress", alpha_exp=0.5, order=TRUE, pt.size=1.5)
