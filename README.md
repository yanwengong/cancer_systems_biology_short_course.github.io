# Cell Ranger
## Cell Ranger Count
## Cell Ranger Aggregation 

# Seurat
Seurat is a R package for analyzing single-cell multi-omics data. 

## Software setup
To install Seurat v4.0+, R version 4.0 or greater is required. We also recommend installing Rstudio, which is an integrated development environment (IDE) designed for R. Rstudio makes R programming much easier.

Install **R**: https://www.r-project.org
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

Please first down the two expression matrix from the GEO links, and save the data at your desired location:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131336
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4131337

## Prepare the directory for the project
It's a good habit to keep project related data in an organized folder. First, create the path in terminal, then assign the path for R as below.  
```
main_path <- "/scratch/yanweng/short_course/"
input_matrix_path <- file.path(main_path, "geo_download")
r_object_path <- file.path(main_path, "r_obj")
plot_path <- file.path(main_path, "plot")
intermediate_data_path <- file.path(main_path, "intermediate_data")
````

## Build the R object



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
