# scEFSC
The scEFSC software package can perform clustering on single-cell RNA sequencing data. 

# How to install:
- The package can be installed from this repository.
- Install.packages("scEFSC_0.1.0.tar.gz",repos = NULL)
- If necessary, install miniconda: `reticulate::install_miniconda(force = T)`
- Install scikit-feature in python,  please visit https://jundongl.github.io/scikit-feature/index.html
- Install tensorflow and keras in python using: `keras::install_keras(tensorflow = "1.10.0")`
- For more information about installation of keras, please visit https://keras.rstudio.com/

# How to use the package for new data 
To use our package for new data, the package includes these functions:  
- scEFSC: main function, doing clustering. The input is a matrix with rows as genes and columns as samples.
- In order to run scEFSC, you need to have feature_selection.py in your run directory. You can find feature_selection.py in the files.
- More detail about parameters for each function could be found in the manual.
