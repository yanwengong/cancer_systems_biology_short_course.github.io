# Basic Command Lines

We will learn some basic command lines to help us nevigate the system fast.

## Access terminal
To open the Mac command lines / terminal, we nevigate to the top bar and click "Go". In the drop-down window, click "Utilities". Then find the Terminal and double click to open it. 

To open a windows terminal :

* Click on Start and type cmd. Select Command Prompt from the list.

* Another option is to use Git bash ( can be downloaded here 
https://urldefense.com/v3/__https://gitforwindows.org/__;!!CzAuKJ42GuquVTTmVmPViYEvSg!Ja5jmSMmrc3R9baahKfqp8IkvfsrBq_oQtzN0ygfkw1keGS2ACMovvZAZT2gmZYTD4OSwt-rr8zDr2c$ )

## Basic command line commands
* pwd: get the current working directory/folder
* ls: list the files/subfolders under the current folder
* cd: get into one directory
* mkdir: make a new directory
* cp: make a copy
* vim: start to write a file
* sh: run shell script
* cat: dislay a file content 
* ...

## Practice

* Print and save the current working directory
* Create the following directories 
 >.. short_course
 >>   ... geo_download\
 >>   ... r_obj\
 >>   ... plot\
 >>   ... intermediate_data

# Cell Ranger
Cell Ranger is a set of analysis pipelines developed/maintained by 10x Genomics. Its functions include aligning reads, generating feature-barcode matrices, clustering and other analysis. Please see details at their [website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). 

We will give examples of their two functions here. You can write save the shell script as "file_name.sh" and run it in the terminal with "sh file_name.sh". 

## CellRanger Count
This function takes the fastq files as input and performs alignment, barcode count and UMI counting.

```
FASTQ_PATH=/share/crsp/lab/kkessenb/yanweng/mdsc/fastq/wt/HA_14WKS_WT_Blood_GEX_1_08262020/HA_14WKS_WT_Blood_GEX_1_08262020
RESULT_PATH="/share/crsp/lab/kkessenb/yanweng/mdsc/cellranger/wt"
cd $RESULT_PATH
cellranger count --id=blood \
                 --fastqs=$FASTQ_PATH \
                 --sample=HA_14WKS_WT_Blood_GEX_1_08262020 \
                 --transcriptome=/pub/yanweng/common_data/mm10/refdata-gex-mm10-2020-A
```
## CellRanger Aggregation 
This function aggregates outputs from multiple runs of cellranger count, normalizing those runs to the same sequencing depth and then recomputing the feature-barcode matrices and analysis on the combined data.

To aggregate samples, you need the sample id and the path to the molecule_info.h5 file produced by cellranger count. We will create a csv like below to store the information.

```
library_id,molecule_h5
spleen_pymt,/share/crsp/lab/kkessenb/yanweng/mdsc/cellranger/spleen2_re20210414/outs/molecule_info.h5
spleen_wt,/share/crsp/lab/kkessenb/yanweng/mdsc/cellranger/wt/lung/outs/molecule_info.h5
```
Then you run the following commands in the terminal.
```
cd /share/crsp/lab/kkessenb/yanweng/mdsc/cellranger
cellranger aggr --id=pymt_wt_agg_20210725 --csv=/pub/yanweng/kai/mdcs/script/pymt_wt_aggr_20210725.csv
```

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

We also recommend installing **tidyverse**, which is a collection of R packages designed for data science. https://www.tidyverse.org

So here, we want to type the following code in Rstudio console to install and load the two of our pakcages.
```
install.packages("tidyverse")
install.packages("Seurat")
library(tidyverse)
library(Seurat)
````

## Prepare the directory for the project
It's a good habit to keep project related data in an organized folder. First, create the path in terminal (by mkdir) or manually, then assign the path for R as below.  

To create the path in the terminal, open the terminal and type the following commands.

```
pwd 
# save this path somewhere # for me it's "/Users/yanwengong" (Mac user)
# for windows, it could be "C:/Users/yanwengong/Documents"
mkdir short_course
cd short_course
mkdir geo_download
mkdir r_obj
mkdir plot
mkdir intermediate_data
```

Then go back to Rstudio, and assign the created path as variables. 
```
# Note: you need to create one dir on your laptop for this course, and set it as main path here
main_path <- "/scratch/yanweng/short_course/" 
# Note: geo_download, r_obj, plot and intermediate_data are the subfolders under the previous main path
input_matrix_path <- file.path(main_path, "geo_download") # where you will store your downloaded GEO data
r_object_path <- file.path(main_path, "r_obj")
plot_path <- file.path(main_path, "plot")
intermediate_data_path <- file.path(main_path, "intermediate_data")
````

## Prepare the source data
For this tutorial, we will analyze a dataset from "Alshetaiwi et al (2020). Defining the emergence of myeloid-derived suppressor cells in breast cancer using single-cell transcriptomics. Science immunology". It contains several immune cells profiled from tumor-bearing mice and healthy mice. We will learn how to identify different cell type in the dataset first. Then subset the data to focus on neutrophils and compare the neutrophils between two condition. 

Please first down the two expression matrix from the GEO links, and save the data at your input_matrix_path:\\
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131336 \\
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131337

Then decompress the downloaded files by gunzip
```
gunzip GSM4131336_Wild_Type_Expression_Matrix.txt.gz
gunzip GSM4131337_PYMT_Expression_Matrix.txt.gz
```

Move the decompressed files into the previously defined input_matrix_path.

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

*Question*: How to check the count matrix stored in the Seurat object? What about the meta data?
```
GetAssayData(wt, slot="counts")
wt@meta.data
```

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
DefaultAssay(pymt_wt_integrated) <- "RNA"
pymt_wt_integrated[["percent.mt"]] <- PercentageFeatureSet(pymt_wt_integrated, pattern = "^mt-")
pdf(file.path(plot_path, "violin_qc.pdf"), width = 10, height = 3)
VlnPlot(pymt_wt_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
pymt_wt_integrated <- subset(pymt_wt_integrated, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
```
![QC](/docs/assets/violin_qc.png)

## Examine the inegrated results

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
# Note: if you do not have enought memory, type "memory.limit(size=500000)", then re-try
pymt_wt_integrated <- RunPCA(pymt_wt_integrated, npcs = 30, verbose = FALSE)
ElbowPlot(pymt_wt_integrated, ndims = 30) 
# Check Point: you can evaluate what is the best PC number to choose here, 
# and use it in the "dims" variable below![image](https://user-images.githubusercontent.com/11083640/173200532-2a631215-7327-4ba3-99b9-c49463a72c4d.png)
![image](https://user-images.githubusercontent.com/11083640/173200538-340dde12-217a-446e-951b-dc124052db8b.png)
![image](https://user-images.githubusercontent.com/11083640/173200540-ded19d7f-d0c0-4f13-a2b7-34c8185f0a32.png)

# dims= 1:20 means use the top 20 PCs
pymt_wt_integrated <- RunUMAP(pymt_wt_integrated, reduction ="pca", dims= 1:20)
pymt_wt_integrated <- FindNeighbors(pymt_wt_integrated, reduction = "pca", dims = 1:20)
pymt_wt_integrated <- FindClusters(pymt_wt_integrated, resolution= 0.2)
# Check Point: you can adjust resolution by changing the values in "dims". 
# For instance, if you want to decrease resolution, you can use "dims = 0.15"

# Visualization
p1 <- DimPlot(pymt_wt_integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pymt_wt_integrated, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file.path(plot_path, "umap_mouse_type_pc20_res02.pdf"), width =5, height = 3)
p1
dev.off()

pdf(file.path(plot_path, "umap_cluster_pc20_re02.pdf"), width = 5, height = 3)
p2
dev.off()
```
<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_mouse_type_cluster_pc20_re02.png" alt="umap_mouse_type_cluster_pc20_re02" width="400"/>
<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_cluster_pc20_re02.png" alt="umap_cluster_pc20_re02" width="400"/>
## Generate a Heatmap

We can plot heatmap to exam the overexpressed genes for each cluster. The heatmap can also help us view the top differentially expressed genes for each clusters so that we can assign cell type label. In addition, if we see two cluster share similar genes we can decide to merge them or change the clustering resolution. 


```
##################################################################
########### Identify cell type markers ###########################
##################################################################
# For performing differential expression after integration, we switch back to the original data
DefaultAssay(pymt_wt_integrated) <- "RNA"
Idents(pymt_wt_integrated) <- "seurat_clusters"
pymt_wt_integrated <- ScaleData(pymt_wt_integrated, verbose = FALSE)
markers <- FindAllMarkers(pymt_wt_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(markers, file.path(intermediate_data_path, "nonconserve_markers_integration_seurat_cluster.csv"))
top5 <- markers %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
pdf(file.path(plot_path, "heatmap_markers_integration_cellType_seurat_cluster.pdf"), width = 15, height = 10)
DoHeatmap(pymt_wt_integrated, 
          features = as.character(top5$gene), 
          label = FALSE,
          raster=FALSE) + NoLegend() 
dev.off()

```
![heatmap_markers_integration_cellType_seurat_cluster](/docs/assets/heatmap_markers_integration_cellType_seurat_cluster.png)
## Marker genes expression on UMAP

We can show the known marker genes expression on UMAP.

```
selected_markers <- c("Cd19", "Cd22", "Cd79a", "Cd3e", "Cd4", "Cd8a",
                      "Csf1r", "Ccr2", "Ly6g", "Cxcr2")
pdf(file.path(plot_path, "umap_marker_genes.pdf"), width =15, height = 6)
FeaturePlot(pymt_wt_integrated, features = selected_markers, ncol = 5)
dev.off()

```

<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_marker_genes.png" alt="umap_marker_genes">

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
# plot the UMAP
pdf(file.path(plot_path, "umap_new_cluster.pdf"), width = 5, height = 3)
DimPlot(pymt_wt_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# plot the heatmap 
markers <- FindAllMarkers(pymt_wt_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(intermediate_data_path, "markers_integration_cellType.csv"))
top10 <- markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
pdf(file.path(plot_path, "heatmap_markers_integration_cellTyper.pdf"), width = 15, height = 10)
DoHeatmap(pymt_wt_integrated, 
          features = as.character(top5$gene), 
          label = FALSE,
          raster=FALSE) + NoLegend() 
dev.off()
# save the integrated object
saveRDS(pymt_wt_integrated, file.path(r_object_path, 
                                      "agg_pymt_wt_integrated.rds"))
```
![heatmap_markers_integration_cellTyper](/docs/assets/heatmap_markers_integration_cellTyper.png)
<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_new_cluster.png" alt="umap_new_cluster" width="400"/>
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
neut <- ScaleData(neut, verbos = FALSE)
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

<img src="/cancer_systems_biology_short_course.github.io/docs/assets/neut_umap_mouse_type_pc20_res02.png" alt="neut_umap_mouse_type_pc20_res02" width="400"/>
<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_cluster_pc20_re02.png" alt="umap_cluster_pc20_re02" width="400"/>
## Examine the known marker genes and compute gene signature score

Similar as the analysis on the integrated object, we can plot the known G-MDSC markers on the UMAP. In addition, we can also compute the gene signature score, which allows as to exam a larger list of genes by taking control genes into consideration. Here is the document for the function: https://satijalab.org/seurat/reference/addmodulescore

We are going to use the previous identified G-MDSC gene list. <paper heatmap here>

```
DefaultAssay(neut) <- "RNA"
neut <- ScaleData(neut, verbos = FALSE)
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

<img src="/cancer_systems_biology_short_course.github.io/docs/assets/umap_marker_genes_neut.png" alt="umap_marker_genes_neut">
<img src="/cancer_systems_biology_short_course.github.io/docs/assets/violinboxplot_gmdsc_score_top10.png" alt="violinboxplot_gmdsc_score_top10" width="400"/>
        
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

<img src="/cancer_systems_biology_short_course.github.io/docs/assets/Neut_cellpercent_dist_cluster_vsvs_mouse.png" alt="Neut_cellpercent_dist_cluster_vsvs_mouse" width="400"/>
## Practice


* Identify the differentially expressed genes in G-MDSC compared with WT neutrophils.
  - Assign cluster1.pymt and cluster5.pymt as gMDSC, the other pymt clusters as neutrophils
  - Assign all cluster in pymt as neutrophils
  - Identify markers by FindAllMarkers
  - Plot heatmap
* Subset the monocyte groups and identify the M-MDSC group.
  - Similar analysis as above but on the monocyte groups
  
