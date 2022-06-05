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
![Book logo](/docs/assets/violin_qc.png)
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
for (i in c("N0", "N1", "N2")){
  name <- paste(as.character(i), "markers", sep = "_")
  print(name)
  ident.1 <- as.character(i)
  val <- FindConservedMarkers(pymt_wt_integrated_major, ident.1 = ident.1, 
                              grouping.var = "mouse_type", verbose = FALSE)
  file_name=paste("conserved_cluster", name, sep = "")
  write.csv(val, file.path(intermediate_data_path, file_name))
}

remove("conserve_marker_df")
for (i in unique(Idents(pymt_wt_integrated_major))){
  file_name <- paste("conserved_cluster", as.character(i), "_markers", sep = "")
  print(file_name)
  val <- read.csv(file.path(intermediate_data_path, file_name)) %>%
    mutate(cluster=i) %>% rename("gene"="X")
  if (exists("conserve_marker_df")){
    # if (as.character(i) == "14"){
    #   new_col <- setdiff(colnames(conserve_marker_df), colnames(val))
    #   val[new_col] <- NA
    #   val <- val %>% select(colnames(conserve_marker_df))
    #   conserve_marker_df <- rbind(conserve_marker_df, val)
    # }
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

Italicized text is the *cat's meow*.


* This is the first list item.
* Here's the second list item.

    I need to add another paragraph below the second list item.

* And here's the third list item.


At the command prompt, type `nano`.

blow is a R code block
```
library(Seaurt)
install.packages()
````
code ends
