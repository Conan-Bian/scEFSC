gene.filtering <- function(data.list, original_dim, batch_size, ncores.ind, wdecay, seed)
{
  or <- list()
  cl <- parallel::makeCluster(3, outfile = "/dev/null")
  registerDoParallel(cl, cores = 3)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
    library(tensorflow)
  })
  or <- foreach(i = seq(3)) %dopar%
    {
      if (is.null(seed))
      {
        config <- list()
        config$intra_op_parallelism_threads <- ncores.ind
        config$inter_op_parallelism_threads <- ncores.ind
        session_conf <- do.call(tf$ConfigProto, config)
        sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
        keras::k_set_session(session = sess)
      } else {
        set.seed((seed+i))
        use_session_with_seed((seed+i))
      }
      
      data.tmp <- data.list[[i]]
      batch_size <-round(nrow(data.tmp)/50)
      
      x <- keras::layer_input(shape = c(original_dim))
      h <- keras::layer_dense(x, 50, kernel_constraint = keras::constraint_nonneg())
      x_decoded_mean <- keras::layer_dense(h, original_dim)
      vae <- keras::keras_model(x, x_decoded_mean)
      magrittr::`%>%`(vae,keras::compile(optimizer =   tf$contrib$opt$AdamWOptimizer(wdecay,1e-3), loss = 'mse'))
      
      
      his <- magrittr::`%>%`(vae,keras::fit(
        data.tmp, data.tmp,
        shuffle = TRUE,
        epochs = 10,
        batch_size = batch_size,
        verbose = 0
      )) 
      
      W <- keras::get_weights(keras::get_layer(vae, index = 2))[[1]]
      Wsd <- matrixStats::rowSds(W)
      Wsd[is.na(Wsd)] <- 0
      Wsd <- (Wsd-min(Wsd))/(max(Wsd)-min(Wsd))
      Wsd
    }
  
  parallel::stopCluster(cl)
  or <- matrixStats::rowMeans2(as.matrix(data.frame(or)))
  or
}

scDHA_FS <- function(data = data, k = NULL, K = 3, n = 5000, ncores = 15L, gen_fil = T,sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-6
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > n))
  n <- ifelse(gen_fil, min(n, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 15
  
  if (gen_fil) {
    data.list <- lapply(seq(3), function(i) {
      if(!is.null(seed)) set.seed((seed+i))
      if(nrow(data) > 2000)
      {
        ind <- sample.int(nrow(data), 2000, replace = F, prob = sample.prob)
      } else {
        ind <- seq(nrow(data))
      }
      data.tmp <- as.matrix(data[ind,])
      data.tmp
    })
    
    or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)
    
    keep <- which(or > quantile(or,(1-min(n,original_dim)/original_dim)))
    
    da <- data[,keep]
    original_dim_reduce <- ncol(da)
    or.da <- da
  } else {
    da <- data
    original_dim_reduce <- ncol(da)
    or.da <- da
  }
  return(or.da)
}

Seurat_cluster <- function(data){
  rownames(data) <- paste0("GENE",c(1:dim(data)[1]))
  colnames(data) <- paste0("CELL",c(1:dim(data)[2]))
  seurat_obj <- CreateSeuratObject(data, project = "SEURAT")
  all.genes <- rownames(data)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj,features = rownames(data))
  seurat_obj <- FindNeighbors(seurat_obj,dims = 1:10)
  seurat_obj <- FindClusters(object = seurat_obj)
  return(as.numeric(Idents(seurat_obj)))
}

CIDR_cluster <- function(data, n){
  sc_cidr <- scDataConstructor(data)
  sc_cidr <- determineDropoutCandidates(sc_cidr)
  sc_cidr <- wThreshold(sc_cidr)
  sc_cidr <- scDissim(sc_cidr)
  sc_cidr <- scPCA(sc_cidr,plotPC = FALSE)
  sc_cidr <- nPC(sc_cidr)
  sc_cidr <- scCluster(sc_cidr, nCluster = n)
  return(sc_cidr@clusters)
}

monocle_cluster <- function(data, n){
  cds <- newCellDataSet(data)
  cds <- estimateSizeFactors(cds)
  cds <- reduceDimension(cds, max_components = 3, num_dim = 10, reduction_method = 'tSNE', verbose = T, perplexity =5
                         ,norm_method = "none",check_duplicates=FALSE)
  cds <- clusterCells(cds, num_clusters = n+1)
  labels <- as.numeric(levels(cds$Cluster))[cds$Cluster]
  return(labels)
}

RaceID_cluster <- function(data,n){
  data <- as.data.frame(data)
  rownames(data) <- 1:dim(data)[1]
  colnames(data) <- 1:dim(data)[2]
  sc <- SCseq(data)
  sc <- filterdata(sc, mintotal=1, minexpr=5, minnumber=1)
  sc <- compdist(sc,metric="pearson")
  sc <- clustexp(sc,sat = FALSE,clustnr=20,bootnr=50, cln=n,rseed=17000)
  labels <- sc@cluster$kpart
  return(labels)
}

scEFSC <- function(data, n, normalize = F, dim1 = 5000,dim2 = 2000){
  library("reticulate")
  library("scDHA")
  library("doParallel")
  library("foreach")
  library("keras")
  library("tensorflow")
  library("SC3")
  library("SingleCellExperiment")
  library("SIMLR")
  library("Seurat")
  library("SHARP")
  library("cidr")
  library("SINCERA")
  library("mclust")
  library("clues")
  library("pcaReduce")
  library("SIMLR")
  library("caret")
  library("Rphenograph")
  library("monocle")
  library("RaceID")
  
  source_python("feature_selection.py")
  
  #normalization
  if(normalize == F){
    data <- log2(data+1)
    data <- data[apply(data,1, function(x) sum(x>1) > floor(ncol(data)/50)),]
  }
  data <- t(data)
  #feature selection 1
  if(dim(data)[2] > dim1)
  {
    data <- scDHA_FS(data = data,k = n,n = dim1)
    print("Feature selection 1 is performed")
  }else
  {
    print("Feature selection 1 is not performed")
  }
  #feature selection 2
  if(dim(data)[2] > dim2)
  {
    data1 <- fs_lap_score(data,as.integer(dim2))
    data2 <- fs_low_variance(data,as.integer(dim2))
    data3 <- fs_SPEC(data,as.integer(dim2))
    data4 <- fs_MCFS(data,as.integer(dim2),n)
    print("Feature selection 2 is performed")
  }else
  {
    print("Feature selection 2 is not performed")
  }
  
  data <- t(data)
  data1 <- t(data1)
  data2 <- t(data2)
  data3 <- t(data3)
  data4 <- t(data4)
  
  #SC3
  sce <- SingleCellExperiment(assays = list(counts = data,logcounts = log2(data+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data))
  labels_SC3 <- as.numeric(dat[[1]])
  
  sce <- SingleCellExperiment(assays = list(counts = data1,logcounts = log2(data1+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data1)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data1))
  labels_SC3_1 <- as.numeric(dat[[1]])
  
  sce <- SingleCellExperiment(assays = list(counts = data2,logcounts = log2(data2+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data2)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data2))
  labels_SC3_2 <- as.numeric(dat[[1]])
  
  sce <- SingleCellExperiment(assays = list(counts = data3,logcounts = log2(data3+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data3)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data3))
  labels_SC3_3 <- as.numeric(dat[[1]])
  
  sce <- SingleCellExperiment(assays = list(counts = data4,logcounts = log2(data4+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data4)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data4))
  labels_SC3_4 <- as.numeric(dat[[1]])
  rm(sce,dat)
  print("SC3 is performed")
  
  #Seurat
  labels_Seurat <- Seurat_cluster(data)
  
  labels_Seurat_1 <- Seurat_cluster(data1)
  
  labels_Seurat_2 <- Seurat_cluster(data2)
  
  labels_Seurat_3 <- Seurat_cluster(data3)
  
  labels_Seurat_4 <- Seurat_cluster(data4)
  print("Seurat is performed")
  
  #SHARP
  res <- SHARP(data, prep = TRUE, n.cores = 5)
  labels_SHARP <- as.numeric(as.factor(res$pred_clusters)) 
  
  res <- SHARP(data1, prep = TRUE, n.cores = 5)
  labels_SHARP_1 <- as.numeric(as.factor(res$pred_clusters))
  
  res <- SHARP(data2, prep = TRUE, n.cores = 5)
  labels_SHARP_2 <- as.numeric(as.factor(res$pred_clusters))
  
  res <- SHARP(data3, prep = TRUE, n.cores = 5)
  labels_SHARP_3 <- as.numeric(as.factor(res$pred_clusters))
  
  res <- SHARP(data4, prep = TRUE, n.cores = 5)
  labels_SHARP_4 <- as.numeric(as.factor(res$pred_clusters))
  rm(res)
  print("SHARP is performed")
  
  #cidr
  labels_cidr <- CIDR_cluster(data, n)
  
  labels_cidr_1 <- CIDR_cluster(data1, n)
  
  labels_cidr_2 <- CIDR_cluster(data2, n)
  
  labels_cidr_3 <- CIDR_cluster(data3, n)
  
  labels_cidr_4 <- CIDR_cluster(data4, n)
  print("cidr is performed")
  
  #SINCERA
  dat <- apply(data, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA <- cutree(hc, k = n)
  
  dat <- apply(data1, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA_1 <- cutree(hc, k = n)
  
  dat <- apply(data2, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA_2 <- cutree(hc, k = n)
  
  dat <- apply(data3, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA_3 <- cutree(hc, k = n)
  
  dat <- apply(data4, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA_4 <- cutree(hc, k = n)
  print("SINCERA is performed")
  
  #pcaReduce
  pca <- PCAreduce(t(data),1,n-1,'M')
  labels_pcaReduce <- pca[[1]][,1]
  
  pca <- PCAreduce(t(data1),1,n-1,'M')
  labels_pcaReduce_1 <- pca[[1]][,1]
  
  pca <- PCAreduce(t(data2),1,n-1,'M')
  labels_pcaReduce_2 <- pca[[1]][,1]
  
  pca <- PCAreduce(t(data3),1,n-1,'M')
  labels_pcaReduce_3 <- pca[[1]][,1]
  
  pca <- PCAreduce(t(data4),1,n-1,'M')
  labels_pcaReduce_4 <- pca[[1]][,1]
  print("pcaReduce is performed")
  
  #Rphenograph
  Rphenograph_out <- Rphenograph(t(data),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph <- as.numeric(membership(Rphenograph_out[[2]]))
  
  Rphenograph_out <- Rphenograph(t(data1),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph_1 <- as.numeric(membership(Rphenograph_out[[2]]))
  
  Rphenograph_out <- Rphenograph(t(data2),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph_2 <- as.numeric(membership(Rphenograph_out[[2]]))
  
  Rphenograph_out <- Rphenograph(t(data3),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph_3 <- as.numeric(membership(Rphenograph_out[[2]]))
  
  Rphenograph_out <- Rphenograph(t(data4),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph_4 <- as.numeric(membership(Rphenograph_out[[2]]))
  print("Rphenograph is performed")
  
  #monocle
  labels_monocle <- monocle_cluster(data, n)
  
  labels_monocle_1 <- monocle_cluster(data1, n)
  
  labels_monocle_2 <- monocle_cluster(data2, n)
  
  labels_monocle_3 <- monocle_cluster(data3, n)
  
  labels_monocle_4 <- monocle_cluster(data4, n)
  print("monocle is performed")
  
  #RaceID
  labels_RaceID <- RaceID_cluster(data,n)
  
  labels_RaceID_1 <- RaceID_cluster(data1,n)
  
  labels_RaceID_2 <- RaceID_cluster(data2,n)
  
  labels_RaceID_3 <- RaceID_cluster(data3,n)
  
  labels_RaceID_4 <- RaceID_cluster(data4,n)
  print("RaceID is performed")
  
  #remove
  cluster_results <- rbind(labels_SC3, labels_SC3_1, labels_SC3_2, labels_SC3_3, labels_SC3_4, 
                           labels_Seurat, labels_Seurat_1, labels_Seurat_2, labels_Seurat_3, labels_Seurat_4, 
                           labels_SHARP, labels_SHARP_1, labels_SHARP_2, labels_SHARP_3, labels_SHARP_4,
                           labels_cidr, labels_cidr_1, labels_cidr_2, labels_cidr_3, labels_cidr_4,
                           labels_SINCERA, labels_SINCERA_1, labels_SINCERA_2, labels_SINCERA_3, labels_SINCERA_4,
                           labels_pcaReduce, labels_pcaReduce_1, labels_pcaReduce_2, labels_pcaReduce_3, labels_pcaReduce_4,
                           labels_Rphenograph, labels_Rphenograph_1, labels_Rphenograph_2, labels_Rphenograph_3, labels_Rphenograph_4,
                           labels_monocle, labels_monocle_1, labels_monocle_2, labels_monocle_3, labels_monocle_4,
                           labels_RaceID, labels_RaceID_1, labels_RaceID_2, labels_RaceID_2, labels_RaceID_3, labels_RaceID_4)
  ARI=matrix(0,45,45)
  for(i in 1:45){
    for(j in 1:45){
      ARI[i,j] <- mclust::adjustedRandIndex(unlist(cluster_results[i,]), unlist(cluster_results[j,]))
    }
  }
  m1 <- order(apply(ARI,1,var),decreasing=FALSE)[1:10]
  cluster_results <- cluster_results[-m1,]
  print("remove is performed")
  
  #wMetaC
  labels <- SHARP::wMetaC(t(cluster_results), enN.cluster = n)
  labels <- as.numeric(labels$finalC)
  print("wMetaC is performed")
  
  return(labels)
}






