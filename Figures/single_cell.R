## check the eQTM genes in published dataset PMID30135581

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tidyr)

setwd("/groups/umcg-griac/tmp01/projects/umcg-cqi/PIAMA/allergy_prediction/data/scRNA_polyp")
meta_data<-read.table("20180822_PolypAll_cleaned_metadata.txt",header=T)
raw_counts<-read.table("20180822_PolypAll_cleaned_rawdata.txt",header=T)
rownames(raw_counts)<-raw_counts$X
raw_counts<-raw_counts[,-1]

polyp <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.genes = 200, project = "mydata_scRNAseq")
AddMetaData(object=polyp,metadata=meta_data)

polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp1TOT'] <- 'YES'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp2TOT'] <- 'YES'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp3TOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp4TOT'] <- 'YES'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp5TOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp6ATOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp6BTOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp7TOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp8TOT'] <- 'NO'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp9TOT'] <- 'YES'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp11TOT'] <- 'YES'
polyp@meta.data$polyp[polyp@meta.data$orig.ident == 'Polyp12TOT'] <- 'YES'

#Filter cells with UMI counts greater than 12000
polyp <- subset(polyp, subset = nCount_RNA<=12000)

# nomalize
polyp <- NormalizeData(polyp, normalization.method = "LogNormalize", scale.factor = 10000)

#Determine variable genes for input into PCA
polyp <- FindVariableFeatures(polyp, selection.method = "vst", nfeatures = 2000)

# scaling
all.genes <- rownames(polyp)
polyp <- ScaleData(polyp, features = all.genes)

#Run PCA
polyp <- RunPCA(polyp, features = VariableFeatures(object = polyp))

#Visualize PCA gene loadings and PCA space
#VizPCA(polyp,1, num.genes= 60, do.balanced = TRUE)
#VizPCA(polyp,2, num.genes= 60, do.balanced = TRUE)
#VizPCA(polyp,3, num.genes= 60, do.balanced = TRUE)
#VizPCA(polyp,4, num.genes= 60, do.balanced = TRUE)

#PCAPlot(polyp, 1, 2)

#Threshold which PCs to use for further dimensionality reduction
pdf("polyp_elbow_PCA.pdf")
ElbowPlot(polyp)
dev.off()

#Find clusters and run TSNE over the first 12 PCs
polyp <- FindNeighbors(polyp, dims = 1:20)
polyp <- FindClusters(polyp, resolution = 0.8)

polyp <- RunTSNE(polyp, dims.use = 1:20, do.fast = T)
polyp <- RunUMAP(polyp, dims= 1:20)

#Visualize results by clusters identified, original identity of sample, or polyp status
pdf("polyp_tsne.pdf")
DimPlot(polyp, reduction = "tsne")
dev.off()

pdf("polyp_umap.pdf")
DimPlot(polyp, reduction = "umap")
dev.off()

#Display some of the top marker genes for clusters used in Figure 1
pdf("polyp_marker1.pdf")
FeaturePlot(polyp, features = c('KRT5','KRT8','LTF','FOXJ1','TRBC2','CD79A','COL1A2','SPARCL1','DARC'))
dev.off()

pdf("polyp_marker2.pdf")
FeaturePlot(polyp, features = c('HLA-DRA','TYROBP','TPSAB1','CLC'))      
dev.off()

#Identify marker genes for clusters using a ROC test and print table
polyp.markers.roc= FindAllMarkers(polyp,return.thresh = 0.65, only.pos = TRUE, test.use = "roc", verbose = T)
write.table(polyp.markers.roc, 'PolypALLTOT_markers_roc.txt', sep = '\t', col.names= NA)

top.markers<-polyp.markers.roc %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = myAUC)


#Name cells for Figure 1
polyp <- AddMetaData(polyp, NA, col.name = 'subset')
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 0] <- 'Basal'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 1] <- 'Apical'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 2] <- 'Glandular'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 3] <- 'Apical'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 4] <- 'Fibroblast'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 5] <- 'Endothelial'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 6] <- 'PlasmaCell'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 7] <- 'Basal'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 8] <- 'PlasmaCell'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 9] <- 'Apical'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 10] <- 'Myeloid'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 11] <- 'TCell'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 12] <- 'PlasmaCell'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 13] <- 'Ciliated'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 14] <- 'Glandular'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 15] <- 'Basal'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 16] <- 'MastCell'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 17] <- 'Glandular'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 18] <- 'Fibroblast'
polyp@meta.data$subset[polyp@meta.data$RNA_snn_res.0.8 == 19] <- 'PlasmaCell'

pdf("polyp_tsne_annot.pdf")
TSNEPlot(polyp, pt.size = 0.5, group.by = 'subset')
dev.off()

## plot umap with annotation
polyp1<-polyp
new.cluster.ids <- c('Basal','Apical','Glandular','Apical','Fibroblast','Endothelial','PlasmaCell','Basal',
                     'PlasmaCell','Apical','Myeloid','TCell','PlasmaCell','Ciliated','Glandular','Basal',
                     'MastCell','Glandular','Fibroblast','PlasmaCell')

names(new.cluster.ids) <- levels(polyp1)
polyp1 <- RenameIdents(polyp1, new.cluster.ids)

#Identify marker genes for cell types
polyp.markers.roc1= FindAllMarkers(polyp1,return.thresh = 0.65, only.pos = TRUE, test.use = "roc", verbose = T)
top.markers1<-polyp.markers.roc1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = myAUC)

write.table(polyp.markers.roc1, 'PolypALLTOT_markers_ct_roc.txt', sep = '\t', col.names= NA)

## umap plot
pdf("polyp_umap_annot.pdf")
DimPlot(polyp1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

## marker heatmap
pdf("polyp_marker_heatmap.pdf")
DoHeatmap(polyp1, features = top.markers1$gene) + NoLegend()
dev.off()

saveRDS(polyp,file = "polyp_cluster.RDS")
saveRDS(polyp1,file = "polyp_annot.RDS")

## plot heatmap of eqtm genes
load("../DNAm_data/eqtm_gene_list_2modules.Rdata")

keep1<-intersect(rownames(polyp1),DEG.gene_symbol1)
keep2<-intersect(rownames(polyp1),DEG.gene_symbol2)

pdf("heatmap_2modules_allcells.pdf",width = 16,height = 6)
DotPlot(polyp1,features = list("Module1"=keep1,"Module2"=keep2))+RotatedAxis()
dev.off()



