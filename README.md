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

Before doing more analysis, we want to ...



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
