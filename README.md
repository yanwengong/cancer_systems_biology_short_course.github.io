# Cell Ranger
## Cell Ranger Count
## Cell Ranger Aggregation 

# Seurat
Seurat is a R package for analyzing single-cell multi-omics data. 

## Software setup
To install Seurat v4.0+, R version 4.0 or greater is required. We also recommend installing Rstudio, which is an integrated development environment (IDE) designed for R. Rstudio makes R programming much easier.

Install **R**: https://www.r-project.org \\
Install **Rstudio**: https://www.rstudio.com/products/rstudio/

## Install R packages
The general ways to install a package in R is:
```
install.packages("package_name")
````
Then you can load the package by:
```
library(package_name)
````
Please note you only need to install package one time, and need to load package everytime before using it.

Please refer to **Seurat** webpage below for the detailed installation instruction.
https://satijalab.org/seurat/articles/install.html

We also recommend install **tidyverse**, which is a collection of R packages designed for data science. https://www.tidyverse.org

So here, we want to type the following code in Rstudio console to install and load the two of our pakcages.
```
install.packages("tidyverse")
install.packages("Seurat")
library(tidyverse)
library(Seurat)
````
## Prepare the source data
For this tutorial, we will analyze a dataset from "Alshetaiwi et al (2020). Defining the emergence of myeloid-derived suppressor cells in breast cancer using single-cell transcriptomics. Science immunology". It contains several immune cells profiled from tumor-bearing mice and healthy mice. We will learn how to identify different cell type in the dataset first. Then subset the data to focus on neutrophils and compare the neutrophils between two condition. 

Please first down the two expression matrix from the GEO links, and save the data at your desired location:\\
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131336 \\
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131337

## Prepare the directory for the project
It's a good habit to keep project related data in an organized folder. First, create the path in terminal, then assign the path for R as below.  
Note: "main_path" is designed to store the short course related data.
```
main_path <- "/scratch/yanweng/short_course/"
input_matrix_path <- file.path(main_path, "geo_download")
r_object_path <- file.path(main_path, "r_obj")
plot_path <- file.path(main_path, "plot")
intermediate_data_path <- file.path(main_path, "intermediate_data")
````

## Build the Seurat object
We will first read the two downloaded expression metrics into R, then create two Seurat object.
```
wt_count <- read.table(file.path(input_matrix_path, 
                                 "GSM4131336_Wild_Type_Expression_Matrix.txt"),
                       header = TRUE, row.names = 1)
wt <- CreateSeuratObject(counts = wt_count)

pymt_count <- read.table(file.path(input_matrix_path, 
                                 "GSM4131337_PYMT_Expression_Matrix.txt"),
                       header = TRUE, row.names = 1)
pymt <- CreateSeuratObject(counts = pymt_count)
```

*Question*: How to check the count matrix stored in the Seurat object?

## Combine the two Seurat objects

We will use the integration methods described in Stuart, Butler et al, 2019. 
```
##### standard integration #####

# normalize and identify variable features for each dataset independently
pymt_wt.list <- c(pymt, wt)
pymt_wt.list <- lapply(X = pymt_wt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pymt_wt.list)

#################################################
########### Perform integration #################
#################################################

anchors <- FindIntegrationAnchors(object.list = pymt_wt.list, anchor.features = features) # take super long
# this command creates an 'integrated' data assay
pymt_wt_integrated <- IntegrateData(anchorset = anchors)

# save the r object along the analysis process to avoid losing data
saveRDS(pymt_wt_integrated, file.path(r_object_path, "agg_pymt_wt_integrated.rds"))
```

## QC and selecting cells for further analysis

Before doing more analysis, we want to make sure the cells included are not doublets or low quality cells. There are a number of softwares for doublet identification, e.g. Scrublet, DoubletFinder, scDblFinder, etc. Please refer to this protocol by [Xi and Li] (https://www.sciencedirect.com/science/article/pii/S2405471220304592) for benchmarking these methods. They also include detailed code to run the softwares. 

In this tutorial, we will use commands that available in Seurat for QC.
* nFeature: the number of unique genes detected in each cell
* nCount: the total number of molecules detected within each cell
* percent.mt: The percentage of reads that map to the mitochondrial genome

```
pymt_wt_integrated[["percent.mt"]] <- PercentageFeatureSet(pymt_wt_integrated, pattern = "^mt-")
pdf(file.path(plot_path, "violin_qc.pdf"), width = 10, height = 3)
VlnPlot(pymt_wt_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
pymt_wt_integrated <- subset(pymt_wt_integrated, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
```
![QC](/docs/assets/violin_qc.png)

## Exam the inegrated results

We perform dimension reduction and clustering on the integrated dataset.

```
############################################################
########### Perform an integrated analysis #################
############################################################

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(pymt_wt_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pymt_wt_integrated <- ScaleData(pymt_wt_integrated, verbose = FALSE)
pymt_wt_integrated <- RunPCA(pymt_wt_integrated, npcs = 30, verbose = FALSE)
ElbowPlot(pymt_wt_integrated, ndims = 30) #16, 21, 30
pymt_wt_integrated <- RunUMAP(pymt_wt_integrated, reduction ="pca", dims= 1:20)
pymt_wt_integrated <- FindNeighbors(pymt_wt_integrated, reduction = "pca", dims = 1:20)
pymt_wt_integrated <- FindClusters(pymt_wt_integrated, resolution= 0.2)
# Visualization
p1 <- DimPlot(pymt_wt_integrated, reduction = "umap", group.by = "mouse_type")
p2 <- DimPlot(pymt_wt_integrated, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file.path(plot_path, "umap_mouse_type_pc20_res02.pdf"), width =5, height = 3)
p1
dev.off()

pdf(file.path(plot_path, "umap_cluster_pc20_re02.pdf"), width = 5, height = 3)
p2
dev.off()
```
<img src="/docs/assets/umap_mouse_type_cluster_pc20_re02.png" width="300">
<img src="/docs/assets/umap_mouse_type_cluster_pc20_re02.png" alt="umap_mouse_type_cluster_pc20_re02" width="500"/>
<img src="/docs/assets/umap_cluster_pc20_re02.png" alt="umap_cluster_pc20_re02" width="500"/>
## Exam Heatmap

We can plot heatmap to exam the overexpressed genes for each cluster. The heatmap can also help use decide on the clustering resolution. Here, we first identify the conserved markers (regardless of mouse type) for each cluster, save them, and generate the heatmap to display top markers for each cluster. 

```
##################################################################
########### Identify conserved cell type markers #################
##################################################################
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(pymt_wt_integrated) <- "RNA"

# find conserved marker
# they are differentially expressed compared to other clusters, 
# but have similar expression between the two groups (PyMT vs WT) you're actually comparing
for (i in c(Idents(pymt_wt_integrated))){
  name <- paste(as.character(i), "markers", sep = "_")
  print(name)
  ident.1 <- as.character(i)
  val <- FindConservedMarkers(pymt_wt_integrated, ident.1 = ident.1, 
                              grouping.var = "mouse_type", verbose = FALSE)
  file_name=paste("conserved_cluster", name, sep = "")
  write.csv(val, file.path(intermediate_data_path, file_name))
}

remove("conserve_marker_df")
for (i in unique(Idents(pymt_wt_integrated))){
  file_name <- paste("conserved_cluster", as.character(i), "_markers", sep = "")
  print(file_name)
  val <- read.csv(file.path(intermediate_data_path, file_name)) %>%
    mutate(cluster=i) %>% rename("gene"="X")
  if (exists("conserve_marker_df")){
    conserve_marker_df <- rbind(conserve_marker_df, val)
  } else{
    conserve_marker_df <- val
  }
  
}

conserve_marker_df<- conserve_marker_df %>% 
  rowwise() %>% 
  mutate(mean_logFC = mean(c(pymt_avg_logFC, wt_avg_logFC), na.rm = FALSE))%>%
  as.data.frame()
```

## Marker genes experssion on UMAP

We can show the known marker genes expression on UMAP.

```
pymt_wt_integrated <- ScaleData(pymt_wt_integrated, verbose = FALSE)
selected_markers <- c("Cd19", "Cd22", "Cd79a", "Cd3", "Cd4", "Cd8",
                      "Csf1r", "Ccr2", "Ly6g", "Cxcr2")
pdf(file.path(plot_path, "umap_marker_genes.pdf"), width =15, height = 6)
FeaturePlot(pymt_wt_integrated, features = selected_markers, ncol = 5)
dev.off()

```

<img src="/docs/assets/umap_marker_genes.png" alt="umap_marker_genes" width="500"/>
## Give cluster identities

By examing the heatmap and marker gene expressions, we can assign cell label to each cluster. Then, we save the integrated object for future use.

| Cluster ID        | Markers           | Cell Type  |
| ------------- |:-------------:| -----:|
| 0,2,3,5,6| Ly6g, Cxcr2      | neutrophil |
| 1,9,10   | Csf1r, Ccr2      | monoycte   |
| 4,7      | Cd19, Cd22, Cd79a| B cell     |
| 8        | Cd3e, Cd4, Cd8a. | T cell     |


```
new.cluster.ids <- c("neut", "mono", "neut", "neut", "B", "neut",
                     "neut", "B", "T", "mono", "mono")
names(new.cluster.ids) <- levels(pymt_wt_integrated)
pymt_wt_integrated <- RenameIdents(pymt_wt_integrated, new.cluster.ids)
pdf(file.path(plot_path, "umap_new_cluster.pdf"), width = 5, height = 3)
DimPlot(pymt_wt_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# save the integrated object
saveRDS(pymt_wt_integrated, file.path(r_object_path, 
                                      "agg_pymt_wt_integrated.rds"))
```

<img src="/docs/assets/umap_new_cluster.png" alt="umap_new_cluster" width="500"/>
## Subset and save neutrophils

```
neut <- subset(x = pymt_wt_integrated, idents = "neut")

saveRDS(neut, file.path(r_object_path, 
                                      "agg_pymt_wt_neut.rds"))
```

## Re-cluster on neutrophil

We subset and "zoom in" the neutrophil group. We want to identify the G-MDSC group, that expresses the known marker genes and enriched in the PyMT sample. 

```
DefaultAssay(neut) <- "integrated"

# Run the standard workflow for visualization and clustering
neut <- ScaleData(neut, verboseneutFALSE)
neut <- RunPCA(neut, npcs = 30, verbose = FALSE)
ElbowPlot(neut, ndims = 30) #16, 21, 30
neut <- RunUMAP(neut, reduction ="pca", dims= 1:20)
neut <- FindNeighbors(neut, reduction = "pca", dims = 1:20)
neut <- FindClusters(neut, resolution= 0.2)
# Visualization
p1 <- DimPlot(neut, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(neut, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file.path(plot_path, "umap_mouse_type_pc20_res02.pdf"), width =5, height = 3)
p1
dev.off()

pdf(file.path(plot_path, "umap_cluster_pc20_re02.pdf"), width = 5, height = 3)
p2
dev.off()
```

<img src="/docs/assets/umap_mouse_type_pc20_res02.png" alt="umap_mouse_type_pc20_res02" width="500"/>
<img src="/docs/assets/umap_cluster_pc20_re02.png" alt="umap_cluster_pc20_re02" width="500"/>
## Exam the known marker genes and compute gene signiture score

Similar as the analysis on the integrated object, we can plot the known G-MDSC markers on the UMAP. In addition, we can also compute the gene signiture score, which allows as to exam a larger list of genes by taking control genes into consideration. Here is the document for the function: https://satijalab.org/seurat/reference/addmodulescore

We are going to use the previous identified G-MDSC gene list. <paper heatmap here>

```
DefaultAssay(neut) <- "RNA"
selected_markers <- c("Cd84", "Ctsd", "Il1b", "Anxa1", "Arhgdib", "Prdx5")
pdf(file.path(plot_path, "umap_marker_genes.pdf"), width =12, height = 6)
FeaturePlot(neut, features = selected_markers, ncol = 3)
dev.off()

top10 <- c("Cd84", "Ctsd", "Arg2", "Pla2g7", "Il1b", "Clec4e", 
           "Il1f9", "Junb", "Wfdc17", "Clec4d", "BC100530") %>%
  list()

neut <- AddModuleScore(
  object = neut,
  features = top10,
  ctrl = 100,
  name = "top10_g_mdsc_sig"
)

pdf(paste(plot_path, "violinboxplot_gmdsc_score_top10.pdf", sep = "/"), 
    width = 5, height = 3)
VlnPlot(neut, features = "top10_g_mdsc_sig1", 
        pt.size=0, y.max=2, split.by = "orig.ident")+
  geom_boxplot(width=0.3, outlier.shape = NA) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")
dev.off()
```

<img src="/docs/assets/umap_marker_genes.png" alt="umap_marker_genes" width="500"/>
<img src="/docs/assets/violinboxplot_gmdsc_score_top10.png" alt="violinboxplot_gmdsc_score_top10" width="500"/>
        
From the UMAP and gene signiture scores, we can say that cluster 1 and 5 are likely to be G_MDSC cluster. 

## Mouse type fraction in each cluster
        
We want to check the mouse type fraction for each cluster, to confirm that cluter 1 and 5 are mainly PyMT cells.
        
```
## check fraction
data <- neut@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = n/sum(n))
#data$organ <- factor(data$organ, levels = c("BM", "spleen", "blood", "lung", "tumorMFP"))
pdf(file.path(plot_path, "Neut_cellpercent_dist_cluster_vsvs_mouse.pdf"), width = 6, height = 4)
data %>% ggplot(aes(x=orig.ident, y=percentage, fill=seurat_clusters)) +
  geom_bar(colour="black", stat="identity") + 
  facet_grid(.~ seurat_clusters) 
theme_bw() +
  ggtitle("Cell number percentage by cluster and mouse")
dev.off()
  
# last we save the neutrophil object for future analysis
saveRDS(neut, file.path(r_object_path, 
                        "agg_pymt_wt_neut.rds"))
```

<img src="/docs/assets/Neut_cellpercent_dist_cluster_vsvs_mouse.png" alt="Neut_cellpercent_dist_cluster_vsvs_mouse" width="500"/>
## Practice


* Identify the differentially expressed genes in G-MDSC compared with WT neutrophils.
* Subset the monocyte groups and identify the M-MDSC group.
  
