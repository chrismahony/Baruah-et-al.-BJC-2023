library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(gsfisher)
library(EnhancedVolcano)
options(bitmapType='cairo')


#preprocessing
AN009_counts<-Read10X("/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/count/AN009/outs/filtered_feature_bc_matrix")
AN014_counts<-Read10X("/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/count/AN014/outs/filtered_feature_bc_matrix")
AN017_counts<-Read10X("/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/count/AN017/outs/filtered_feature_bc_matrix")
AN018_counts<-Read10X("/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/count/AN018/outs/filtered_feature_bc_matrix")

AN009<- CreateSeuratObject(counts = AN009_counts, project = "AN009", min.cells = 0, min.features = 0)
AN014<- CreateSeuratObject(counts = AN014_counts, project = "AN014", min.cells = 0, min.features = 0)
AN017<- CreateSeuratObject(counts = AN017_counts, project = "AN017", min.cells = 0, min.features = 0)
AN018<- CreateSeuratObject(counts = AN018_counts, project = "AN018", min.cells = 0, min.features = 0)

AN009[["percent.mt"]] <- PercentageFeatureSet(AN009, pattern = "^MT-")
AN014[["percent.mt"]] <- PercentageFeatureSet(AN014, pattern = "^MT-")
AN017[["percent.mt"]] <- PercentageFeatureSet(AN017, pattern = "^MT-")
AN018[["percent.mt"]] <- PercentageFeatureSet(AN018, pattern = "^MT-")

VlnPlot(AN009, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN014, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN017, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN018, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AN009 <- subset(AN009, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
AN014 <- subset(AN014, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
AN017 <- subset(AN017, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
AN018 <- subset(AN018, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

VlnPlot(AN009, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN014, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN017, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AN018, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AN009 <- AN009 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

AN014 <- AN014 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

AN017 <- AN017 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

AN018 <- AN018 %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)


#sample AN009 not used due to very low cell numbers

list=c(AN014,AN017, AN018)
anchors <- FindIntegrationAnchors(object.list = list, reduction = "rpca", dims = 1:30)
aggr <- IntegrateData(anchorset = anchors, dims = 1:30)

aggr <- aggr %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

DimPlot(aggr, group.by = "orig.ident")

DimPlot(aggr, split.by = "orig.ident")

DimPlot(aggr, split.by = "orig.ident", cols=cols)

aggr <- FindNeighbors(object = aggr, reduction = "umap", dims = 1:2)


for (res in c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7)) {
   aggr <- FindClusters(aggr, resolution = res, algorithm = 3, graph.name = "integrated_nn")
                          
}

#library(clustree)
#clustree(aggr, assay="integrated_snn")

aggr <- FindClusters(aggr, resolution = c(0.02, 0.05, 0.07), algorithm = 3, graph.name = "integrated_nn")
aggr <- FindClusters(aggr, resolution = c(0.01), algorithm = 3, graph.name = "integrated_nn")


DimPlot(aggr, group.by = "integrated_nn_res.0.01", label = T)
Idents(aggr)<-'integrated_nn_res.0.01'
Markers<-FindAllMarkers(aggr, only.pos = T)
FeaturePlot(aggr, features = c("FBN1"))

table(aggr$integrated_nn_res.0.01)


##### marker gene expression - all cells

DimPlot(aggr, group.by = "integrated_nn_res.0.01", label = T)

DefaultAssay(aggr)<-'RNA'
FeaturePlot(aggr, features = "ALOX15")
FeaturePlot(aggr, features = "ALOX15B")
FeaturePlot(aggr, features = "CD14")
FeaturePlot(aggr, features = "INHBA")
FeaturePlot(aggr, features = "CSF1R")
FeaturePlot(aggr, features = "IL1B")
FeaturePlot(aggr, features = "VIM")
FeaturePlot(aggr, features = "THY1")
FeaturePlot(aggr, features = "PDPN")
FeaturePlot(aggr, features = "FAP")
FeaturePlot(aggr, features = "FN1")
FeaturePlot(aggr, features = "CD248")
FeaturePlot(aggr, features = "COL1A1")
FeaturePlot(aggr, features = "PI16")
FeaturePlot(aggr, features = "ACTA2")
FeaturePlot(aggr, features = "S100A4")
FeaturePlot(aggr, features = "IL1B")
FeaturePlot(aggr, features = "MMP14")
FeaturePlot(aggr, features = "LGAL1")
FeaturePlot(aggr, features = "LGAL9")
FeaturePlot(aggr, features = "VIM", max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(aggr, features = "CD68")
FeaturePlot(aggr, features = "CD163")
FeaturePlot(aggr, features = "CD86")
FeaturePlot(aggr, features = "CD3E")
FeaturePlot(aggr, features = "MT1M")
FeaturePlot(aggr, features = "SOX2")
FeaturePlot(aggr, features = "PMP2")
FeaturePlot(aggr, features = "S100B")
FeaturePlot(aggr, features = "DCN")
FeaturePlot(aggr, features = "LUM")
FeaturePlot(aggr, features = "CYP1B1")
FeaturePlot(aggr, features = "COL1A1")
FeaturePlot(aggr, features = "ACTA2")
FeaturePlot(aggr, features = "COL1A2")
FeaturePlot(aggr, features = "LUM")
FeaturePlot(aggr, features = "DCN")
FeaturePlot(aggr, features = "SOX2")
FeaturePlot(aggr, features = "PMP2")
FeaturePlot(aggr, features = "S100P")
table(aggr$integrated_nn_res.0.01)


DefaultAssay(aggr)<-'RNA'
DotPlot(aggr, features = c("MAFB", "CD163", "CD14", "APOD", "LUM", "COL1A1", "STAT4", "IL32", "CD3G", "CDH19", "ITGB8", "NCAM1", "SPARC"))

Markers_top<- Markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)

write.csv(Markers_top, "/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/Markers_top.csv")

#proportion plot
pt <- table(aggr$orig.ident, aggr$integrated_nn_res.0.01)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())


#macrophage analysis

macs <- aggr[,grepl("0", aggr$integrated_nn_res.0.01, ignore.case=TRUE)]

macs <- macs %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

DimPlot(macs, group.by = 'orig.ident')

library(harmony)
macs_safe<-macs
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(macs@reductions$umap@cell.embeddings),
  meta_data = macs@meta.data,
  vars_use  = "orig.ident",
  do_pca = FALSE
)
rownames(my_harmony_embeddings) <- rownames(macs@reductions$umap@cell.embeddings)
#store the harmony reduction as a custom dimensional reduction called 'harmony' in the default assay
macs[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(macs))
macs <- FindNeighbors(object = macs, reduction = "harmony", dims = 1:2)
macs <- FindClusters(object = macs, verbose = TRUE, algorithm = 1) # Louvain algorithm
macs <- RunUMAP(object = macs, reduction = "harmony", dims = 1:2)
DimPlot(macs_safe, group.by = "orig.ident", reduction = "umap")
DimPlot(macs_safe, split.by = "orig.ident")
DimPlot(macs, group.by = "orig.ident", reduction = "harmony")

#macrophage gene expression anlaysis

FeaturePlot(macs, features = 'CD163', reduction = 'harmony')

DefaultAssay(macs)<-'RNA'
FeaturePlot(macs, features = 'CD163', reduction = 'harmony')

FeaturePlot(macs, features = 'CD68', reduction = 'harmony')
FeaturePlot(macs, features = 'CD14', reduction = 'harmony')
FeaturePlot(macs, features = 'IL1B', reduction = 'harmony')
FeaturePlot(macs, features = 'INHBA', reduction = 'harmony')
FeaturePlot(macs, features = 'CSF1R', reduction = 'harmony')
FeaturePlot(macs, features = 'ALOX15', reduction = 'harmony')



FeaturePlot(macs, features = 'AREG', reduction = 'harmony')
FeaturePlot(macs, features = 'PLAUR', reduction = 'harmony')
FeaturePlot(macs, features = 'AFF3', reduction = 'harmony')

FeaturePlot(macs, features = 'PLCG2', reduction = 'harmony')
FeaturePlot(macs, features = 'NCKAP5', reduction = 'harmony')

FeaturePlot(macs, features = 'AUTS2', reduction = 'harmony')
FeaturePlot(macs, features = 'SPP1', reduction = 'harmony')
FeaturePlot(macs, features = 'SERPINE1', reduction = 'harmony')


FeaturePlot(macs, features = 'S100A4', reduction = 'harmony')
FeaturePlot(macs, features = 'TNF', reduction = 'harmony')
FeaturePlot(macs, features = 'IL10', reduction = 'harmony')
FeaturePlot(macs, features = 'IFNG', reduction = 'harmony')
FeaturePlot(macs, features = 'VEGFA', reduction = 'harmony')
FeaturePlot(macs, features = 'VEGFB', reduction = 'harmony')

FeaturePlot(macs, features = 'GFAP', reduction = 'harmony')
FeaturePlot(macs, features = 'MKI67', reduction = 'harmony')
FeaturePlot(macs, features = 'MERTK', reduction = 'harmony')
FeaturePlot(macs, features = 'CCL3', reduction = 'harmony')
FeaturePlot(macs, features = 'CCL4', reduction = 'harmony')
FeaturePlot(macs, features = 'TMEM119', reduction = 'harmony')
FeaturePlot(macs, features = 'P2RY12', reduction = 'harmony')


FeaturePlot(macs, features = 'TGFB', reduction = 'harmony')
FeaturePlot(macs, features = 'IBA1', reduction = 'harmony')
FeaturePlot(macs, features = 'CCL2', reduction = 'harmony')
FeaturePlot(macs, features = 'CCL5', reduction = 'harmony')
FeaturePlot(macs, features = 'TIE2', reduction = 'harmony')
FeaturePlot(macs, features = 'PDL1', reduction = 'harmony')
FeaturePlot(macs, features = 'SIGLEC15', reduction = 'harmony')
FeaturePlot(macs, features = 'TYRO3', reduction = 'harmony')
FeaturePlot(macs, features = 'AXL', reduction = 'harmony')

FeaturePlot(macs, features = 'IL8', reduction = 'harmony')
FeaturePlot(macs, features = 'SIRPA', reduction = 'harmony')

FeaturePlot(macs, features = 'EGFR', reduction = 'harmony')
#FeaturePlot(macs, features = 'RAGE', reduction = 'harmony')
FeaturePlot(macs, features = 'IL6', reduction = 'harmony')

FeaturePlot(macs, features = 'PLAUR', reduction = 'harmony')
FeaturePlot(macs, features = 'MIP1A', reduction = 'harmony')




FeaturePlot(macs, features = "SLC2A1", reduction = 'harmony')
FeaturePlot(macs, features = "PFKFB3", reduction = 'harmony')
FeaturePlot(macs, features = "SLC16A1", reduction = 'harmony')
FeaturePlot(macs, features = "SLC16A3", reduction = 'harmony')
FeaturePlot(macs, features = "S100P", reduction = 'harmony')


DotPlot(macs, features = c('PLAUR', "AXL"))

DefaultAssay(macs)<-"RNA"
DotPlot(macs, features = c('AURKB', "CDCA3", "ASPM", "HLA-DQA2", "HLA-DQB2", "MT1G", "ANKRD28"))+RotatedAxis()

#subcluster macs
macs <- FindNeighbors(object = macs, reduction = "harmony", dims = 1:2)
macs <- FindClusters(object = macs, verbose = TRUE, resolution = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1)) 
library(clustree)
clustree(macs)
cols <- ArchR::paletteDiscrete(macs@meta.data[, "integrated_snn_res.0.02"])

DimPlot(macs, group.by="integrated_snn_res.0.02", reduction = "harmony", cols=cols)

Idents(macs)<-'integrated_snn_res.0.02'
macs_markers<-FindAllMarkers(macs, only.pos = T)
macs_markers_top<- macs_markers %>%
    group_by(cluster) %>%
    top_n(n = 30, wt = avg_log2FC)
write.csv(macs_markers_top, '/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/macs_markers_top.csv')

acs_markers_top3<- macs_markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)


macs_markers_top10<- macs_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
macs<-ScaleData(macs)
DoHeatmap(macs, features = macs_markers_top10$gene)
DotPlot(macs, features = macs_markers_top3$gene)+RotatedAxis()


dotplot<-DotPlot(macs, features = unique(macs_markers_top10$gene))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot)


library(gsfisher)

getExpressedGenesFromSeuratObject <- function(seurat_object,
                                              clusters,
                                              min.pct=0.1)
{
  expressed <- c()
  for(cluster in clusters)
  {
    # get genes detected in the cluster
    cluster_cells <- names(seurat_object@active.ident[seurat_object@active.ident==cluster])
    clust_pcts <- apply(seurat_object@assays$RNA@data[,cluster_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_clust <- names(clust_pcts[clust_pcts>min.pct])
    
    # get genes detected in the other cells
    other_cells <- names(seurat_object@active.ident[seurat_object@active.ident!=cluster])
    other_pcts <- apply(seurat_object@assays$RNA@data[,other_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_other_cells <- names(other_pcts[other_pcts>min.pct])
    
    expressed <- c(expressed, detected_in_clust, detected_in_other_cells)
  }
  expressed <- unique(expressed)
  expressed
}

Idents(macs)<-'integrated_snn_res.0.02'
expressed_genes <- getExpressedGenesFromSeuratObject(macs,levels(macs@active.ident), min.pct=0.1)


#expressed_genes<-rownames(macs)
annotation_gs <- fetchAnnotation(species="hs", ensembl_version=NULL, ensembl_host=NULL)

index <- match(macs_markers$gene, annotation_gs$gene_name)
macs_markers$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- expressed_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

seurat_obj.res <- macs_markers
seurat_obj <- macs
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$ensembl),]
ensemblUni <- na.omit(ensemblUni)

go.results <- runGO.all(results=seurat_obj.res,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

write.csv(go.results.top , "/rds/projects/c/croftap-lab-data-2021/Paramita_SingleCell/go_results_top.csv")

#TF analysis

                        Idents(macs)<-'integrated_snn_res.0.02'
macs_f<-subset(macs, idents=c("0", "1", "2"))

macs_f <- macs_f %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)


library(OmnipathR)
library(dorothea)
library(decoupleR)
library(SCpubr)

network <- decoupleR::get_dorothea(organism = "human",
                                   levels = c("A", "B", "C"))


activities <- decoupleR::run_wmean(mat = as.matrix(macs_f@assays[['RNA']]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "mor",
                                   times = 100,
                                   minsize = 5)
Idents(macs)<-'integrated_snn_res.0.02'
out <- SCpubr::do_TFActivityPlot(sample = macs_f,
                                 activities = activities,
                                 )
p <- out$heatmaps$average_scores
p

Heatmap(ht, cluster_columns = T, cluster_rows = T)


#revirews comments
macs_test <- aggr[,grepl("0", aggr$integrated_nn_res.0.01, ignore.case=TRUE)]

macs_test <- macs_test %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

DimPlot(macs_test, group.by = 'orig.ident')

DimPlot(macs, group.by = 'orig.ident', reduction = "harmony")


Idents(aggr) <- 'integrated_nn_res.0.01'
DotPlot(aggr, features=c("CD163","CD14","CSFR1","INHBA","ALOX15","IL1B","VEGFA","TNF","IL10","S100A4","CCL3","CCL4"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


Idents(aggr) <- 'integrated_nn_res.0.01'
DotPlot(aggr, features=c("MRC1", "CD80", "AIF1"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


DotPlot(macs, features=c("MRC1", "CD80", "AIF1"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


DotPlot(macs, features=c("IL4", "IL13", "MRC1", "CCL1", "TNFS14", "IL10", "MERTK", "TGFB", "VEGFA"))+RotatedAxis()+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



FeaturePlot(macs, features=c("IL4", "IL13", "MRC1", "CCL1", "TNFS14", "IL10", "MERTK", "TGFB", "VEGFA"), reduction="harmony")

                        

