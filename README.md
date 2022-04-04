# scEFSC: Accurate Single-cell RNA-seq Data Analysis via Ensemble Consensus Clustering Based on Multiple Feature Selections
We developed a single-cell consensus clustering algorithm based on ensemble feature selection (scEFSC) for scRNA-seq data analysis in an ensemble manner. The input of the proposed scEFSC was the scRNA-seq data with rows as genes and columns as samples. To begin data pre-processing was performed using a log2 transformation to normalize the data and then genes detected in the normalized data in less than 2% of the cells were removed to filter out the low-level expressed genes from the scRNA-seq data. The overall framework of our proposed scEFSC is summarized in the figure as follow:
![Image text](https://raw.githubusercontent.com/Conan-Bian/scEFSC/main/img/scEFSC.png) 

As depicted in this figure, scEFSC consists of three important phases. In phase A, inspired by scDHA, we first employed a non-negative kernel autoencoder to pre-select 5000 genes to remove genes that were insignificant. After that, we propose multiple unsupervised feature selections including Low Variance, Laplacian Score, SPEC, and MCFS to remove genes that do not contribute significantly to the analysis of the scRNA-seq data. We then fed the derived feature subsets into the clustering algorithms. Among them, three of four feature selections are extensions of spectral model previously used for scRNA-seq data analysis. Indeed, we can observe that Low Variance is based on statistics. Laplacian Score and SPEC are based on similarity and MCFS is based on sparse learning. Unlike feature extraction methods, these feature selection methods do not change the original representation of the data and are considered to provide better readability and interpretability. 

In phase B, we applied several different scRNA-seq clustering algorithms to cluster the feature subsets obtained by the multiple feature selection models. Various scRNA-seq clustering methods exist to run in our model. These methods are based on different underlying mathematical formulations as described above, including SC3, CIDR, monocle, pcaReduce, Rphenograph, Seurat, SHARP, SINCERA and RaceID. For each feature subset derived from the different feature selection models, we applied the stated clustering methods to generate cluster labels to finally yield a set of individual cluster labels. In addition, to enhance the diversity of the individual cluster labels in the set, the pairwise Adjusted Rand Index (ARI) was employed to measure the similarity between any two individual clustering labels and then remove the method having similarity with the lowest variance. In phase C, a weighted-ensemble clustering method called wMetaC was used to obtain the final clustering result of the individual cluster labels.

# How to install:
- The package can be installed from this repository.
- install.packages("scEFSC_0.1.0.tar.gz",repos = NULL)
- If necessary, install miniconda: `reticulate::install_miniconda(force = T)`
- Install scikit-feature in python,  please visit https://jundongl.github.io/scikit-feature/index.html
- For more information about installation of keras, please visit https://keras.rstudio.com/

# How to use the package for new data 
To use our package for new data, the package includes these functions:  
- scEFSC: main function, doing clustering. The input is a matrix with rows as genes and columns as samples.
- In order to run scEFSC, you need to have feature_selection.py in your run directory. You can find feature_selection.py in the files.
- More detail about parameters for each function could be found in the manual.

# Example
library("scEFSC")

library("SingleCellExperiment")

#read data

d <- readRDS("data/yan.rds")

data <- assay(d)

label <- as.numeric(factor(d$cell_type1))

n <- length(unique(label))

#run scEFSC

scEFSC_labels <- scEFSC(data, n, normalize = F, dim1 = 5000,dim2 = 2000)

#scEFSC_labels is the list of labels obtained by scEFSC.
