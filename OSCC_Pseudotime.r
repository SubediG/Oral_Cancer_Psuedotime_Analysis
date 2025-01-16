
setwd("/Users/lord/singleCell") 
list.files(all.files = TRUE, full.names = TRUE)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
library(harmony)
library(dplyr)
oral_cancer <- readRDS("oral_cancer.rds")
meta_data <- oral_cancer@meta.data



# Since original seurat dataset is misisng feature count and expression count, it was extrated by creating a seuratobject using expressin matrix
count_matrix <- oral_cancer@assays$RNA$counts
oral_cancer_seurat <- CreateSeuratObject(counts = count_matrix , project = "oral cancer" , min.cells = 3, min.features = 200)
meta_data2 <- oral_cancer_seurat@meta.data

#Adding feature count to raw dataset
oral_cancer$nFeature_RNA <- meta_data2$nFeature_RNA
oral_cancer$nCount_RNA <- meta_data2$nCount_RNA
oral_cancer$Mtpercent <- PercentageFeatureSet(oral_cancer, pattern = "^MT")
(oral_cancer@meta.data)


#Creating plot to identify required pre-processing
before_cutoff_plot <- VlnPlot(oral_cancer, features = c("nFeature_RNA", "nCount_RNA" , "Mtpercent"), ncol = 3)


# Creating a feature scatter plot to identify cutoffs and
FeatureScatter(oral_cancer, feature1 = "nFeature_RNA" , feature2 = "nCount_RNA") + geom_smooth(method = "lm")

# determining cutoffs based on the plot
oral_cancer <- subset(oral_cancer, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 16000)
after_scatter <- FeatureScatter(oral_cancer, feature1 = "nFeature_RNA" , feature2 = "nCount_RNA")+ geom_smooth(method = "lm")


# Identification of inerested cell types
unique(oral_cancer@meta.data$cell_type)
Idents(oral_cancer) <- oral_cancer$cell_type
osc_seurat <- subset(oral_cancer, idents = "oral mucosa squamous cell")
osc_metadata <- osc_seurat@meta.data


# Preparation for clustering and visualization
osc_seurat <- NormalizeData(osc_seurat)
osc_seurat <- FindVariableFeatures(osc_seurat)
osc_seurat <- ScaleData(osc_seurat)
osc_seurat <- RunPCA(osc_seurat)
ElbowPlot(osc_seurat)
osc_seurat <- FindNeighbors(osc_seurat, dims = 1:15)
osc_seurat <- FindClusters(osc_seurat, resolution = 0.6)
osc_seurat <- RunUMAP(osc_seurat, reduction = "pca", dims = 1:15)

# Plotting UMAP to visualize batch effects
p1 <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "development_stage_ontology_term_id", label = TRUE)
p3 <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "development_stage", label = TRUE)
p4 <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "library_preparation_batch", label = TRUE)
p1|p2
p3|p4


# Removing batch effect using harmony
osc_seurat <- RunHarmony(
  object = osc_seurat,
  group.by.vars = "library_preparation_batch",
  dims.use = 1:15 # Use the same dimensions you selected in FindNeighbors
)

# Re- running clustering using harmony
osc_seurat <- osc_seurat%>%
  RunUMAP(reduction = 'harmony',dims = 1:15) %>% 
  FindNeighbors(reduction = 'harmony' , dims = 1:15) %>% 
  FindClusters(resolution = 0.6)


p1_after <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "seurat_clusters", label = TRUE)
p2_after <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "development_stage_ontology_term_id", label = TRUE)
p3_after <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "development_stage", label = TRUE)
p4_after <- DimPlot(osc_seurat, reduction= 'umap' , group.by = "library_preparation_batch", label = TRUE)
p3|p3_after
p4|p4_after


# Creating a cds from the seurat
cds <- as.cell_data_set(osc_seurat)
cds
# To built trajectory i am using cluster that i created using seurat and embeding 
#it into cds assigning partition
new_partition <- c(rep(1, length(cds@colData@rownames)))
names(new_partition) <- cds@colData@rownames
new_partition <- as.factor(new_partition)
cds@clusters@listData$UMAP$partitions <- new_partition

# Assigning cluster information, which cell belong to which cluster
list_cluster <- osc_seurat@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assigning UMAP co-ordinates from seurat to cds
cds@int_colData@listData$reducedDims$UMAP <- osc_seurat@reductions$umap@cell.embeddings

#Learn trajectory
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds, 
           color_cells_by = "development_stage", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE, 
           label_leaves = FALSE,
           group_label_size = 5
) 


# From the cluster, it was identified that age 36 cluster will be the choice as a root node to order cell as most undifferentiated taking age factor in consideration.
cell_identifiers_seurat <- rownames(osc_seurat@meta.data)
root_ontology_id <- "36-year-old stage"
root_cells_identifiers <- cell_identifiers_seurat[osc_seurat$development_stage == root_ontology_id]


# Setting the root state in the trajectory using the identified root cells
cds <- order_cells(cds, root_cells = root_cells_identifiers)

plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE, 
           label_roots = TRUE,  # Mark root nodes
           label_leaves = FALSE, 
           show_trajectory_graph = TRUE
)


cds$MonoPseudo <- pseudotime(cds)
pseudotimedata <- as.data.frame(colData(cds))
ggplot(pseudotimedata, aes(MonoPseudo,reorder(development_stage, MonoPseudo, median), fill = development_stage)) + geom_boxplot()

# Identofying differentially expressed genes
def_genes<- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
def_genes %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head()

#plotting gene of interest to identify expression level in a cluster.
FeaturePlot(osc_seurat, features = c("ENSG00000178585" , "ENSG00000285646" , "ENSG00000117616", "ENSG00000175793"))

#Plotting features according to pseudotime expression in a cluster
osc_seurat$pseudotime <- pseudotime(cds)
FeaturePlot(osc_seurat, features = "pseudotime")
