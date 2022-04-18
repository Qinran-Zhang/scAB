
# scAB: Multiresolution dissection of phenotype-associated cell states


## Capabilities
scAB is an R toolkit for simultaneously dissecting multiresolution phenotype-associated cell states and their predictive signatures. The multiresolution cell states include both coarse- and fine-grain cell states, which respectively correspond to all phenotype-associated cells and distinct subsets of phenotype-associated cells with distinct clinical significance. By integrating single-cell and bulk RNA-seq data with phenotype information, scAB exhibits superior performance in resolving multiresolution cell states and molecular signatures for clinical prognosis.


## Installation

scAB R package can be easily installed from Github using devtools:  

```
devtools::install_github("Qinran-Zhang/scAB")
```


### Installation of other dependencies
- Install [Seurat](https://github.com/satijalab/seurat) using `install.packages('Seurat')`. 
- Install [diptest](https://cran.r-project.org/web/packages/diptest/index.html) using `install.packages('diptest')`.
- Install [multimode](https://cran.r-project.org/web/packages/multimode/index.html) using `install.packages('multimode')`.
- Install [preprocessCore](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html) using `if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("preprocessCore")`.




