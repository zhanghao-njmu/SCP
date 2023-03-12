#' Run Non-negative matrix factorization(NMF)
#'
#' Run a NMF dimensionality reduction.
#'
#' @param object An object
#' @param ... Arguments passed to \link[RcppML]{nmf}.
#'
#' @return Returns Seurat object with the NMF calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunNMF
#' @export RunNMF
#'
RunNMF <- function(object, ...) {
  UseMethod(generic = "RunNMF", object = object)
}

#' @param reduction.name dimensional reduction name, 'nmf' by default
#'
#' @param object
#' @param assay
#' @param slot
#' @param features
#' @param nbes
#' @param nmf.method
#' @param tol
#' @param maxit
#' @param rev.nmf
#' @param verbose
#' @param ndims.print
#' @param nfeatures.print
#' @param reduction.key
#' @param seed.use
#' @param ...
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand VariableFeatures DefaultAssay GetAssay
#' @export
#' @method RunNMF Seurat
#'
RunNMF.Seurat <- function(object, assay = NULL, slot = "data", features = NULL, nbes = 50,
                          nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                          ndims.print = 1:5, nfeatures.print = 30, reduction.name = "nmf",
                          reduction.key = "BE_", verbose = TRUE, seed.use = 42,
                          ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunNMF(
    object = assay.data,
    assay = assay,
    slot = slot,
    features = features,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param features Features to compute nmf on. If features=NULL, nmf will be run
#' using the variable features for the Assay.
#'
#' @param object
#' @param assay
#' @param slot
#' @param nbes
#' @param nmf.method
#' @param tol
#' @param maxit
#' @param rev.nmf
#' @param verbose
#' @param ndims.print
#' @param nfeatures.print
#' @param reduction.key
#' @param seed.use
#' @param ...
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
#' @method RunNMF Assay
#'
RunNMF.Assay <- function(object, assay = NULL, slot = "data", features = NULL, nbes = 50,
                         nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                         ndims.print = 1:5, nfeatures.print = 30,
                         reduction.key = "BE_", verbose = TRUE, seed.use = 42,
                         ...) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(object = object, slot = slot)
  features.var <- apply(
    X = data.use[features, ], MARGIN = 1,
    FUN = var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunNMF(
    object = data.use,
    assay = assay,
    slot = slot,
    nbes = nbes,
    nmf.method = nmf.method,
    tol = tol,
    maxit = maxit,
    rev.nmf = rev.nmf,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param assay Name of Assay nmf is being run on.
#'
#' @param slot Name of slot nmf is being run on.
#' @param nbes Total Number of BEs("basis experiment", aka "metagene") to compute and store (50 by default).
#' @param rev.nmf By default computes the nmf on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param verbose Print the top genes associated with high/low loadings for
#' the BEs.
#' @param ndims.print BEs to print genes for.
#' @param nfeatures.print Number of genes to print for each BE.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "BE_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param object
#' @param nmf.method
#' @param tol
#' @param maxit
#' @param ...
#'
#' @importFrom utils capture.output
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @export
#' @method RunNMF default
RunNMF.default <- function(object, assay = NULL, slot = "data", nbes = 50,
                           nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                           ndims.print = 1:5, nfeatures.print = 30, reduction.key = "BE_",
                           verbose = TRUE, seed.use = 42, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.nmf) {
    object <- t(object)
  }
  nbes <- min(nbes, nrow(x = object) - 1)
  if (nmf.method == "RcppML") {
    check_R("RcppML")
    require("Matrix", quietly = TRUE)
    nmf.results <- RcppML::nmf(
      A = t(object), k = nbes, tol = tol, maxit = maxit,
      seed = seed.use, verbose = verbose, ...
    )
    cell.embeddings <- nmf.results$w
    feature.loadings <- t(nmf.results$h)
    d <- nmf.results$d
    iter <- nmf.results$iter
    # mse <- mse(t(object), cell.embeddings, nmf.results$d, t(feature.loadings))
  }
  if (nmf.method == "NMF") {
    check_R("NMF")
    nmf.results <- NMF::nmf(x = as.matrix(t(object)), rank = nbes, seed = seed.use)
    cell.embeddings <- nmf.results@fit@W
    feature.loadings <- t(nmf.results@fit@H)
    d <- iter <- tol <- NULL
    # mse <- mse(t(object), cell.embeddings, NULL, t(feature.loadings))
  }

  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nbes)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, d = d, tol = tol, iter = iter)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = "\n"))
  }
  return(reduction.data)
}

#' Run multi-dimensional scaling(MDS)
#'
#' Run a MDS dimensionality reduction.
#'
#' @param object An object
#' @param ...
#'
#' @return Returns Seurat object with the MDS calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunMDS
#' @export RunMDS
#'
RunMDS <- function(object, ...) {
  UseMethod(generic = "RunMDS", object = object)
}

#' @param reduction.name dimensional reduction name, 'mds' by default
#'
#' @param object
#' @param assay
#' @param slot
#' @param features
#' @param nmds
#' @param dist.method
#' @param mds.method
#' @param rev.mds
#' @param verbose
#' @param ndims.print
#' @param nfeatures.print
#' @param reduction.key
#' @param seed.use
#' @param ...
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand VariableFeatures DefaultAssay GetAssay
#' @export
#' @method RunMDS Seurat
#'
RunMDS.Seurat <- function(object, assay = NULL, slot = "data",
                          features = NULL, nmds = 50, dist.method = "euclidean", mds.method = "cmdscale",
                          rev.mds = FALSE, ndims.print = 1:5, nfeatures.print = 30,
                          reduction.name = "mds", reduction.key = "MDS_",
                          verbose = TRUE, seed.use = 42, ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunMDS(
    object = assay.data,
    assay = assay,
    slot = slot,
    features = features,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param features Features to compute mds on. If features=NULL, mds will be run
#' using the variable features for the Assay.
#'
#' @param object
#' @param assay
#' @param slot
#' @param nmds
#' @param dist.method
#' @param mds.method
#' @param rev.mds
#' @param verbose
#' @param ndims.print
#' @param nfeatures.print
#' @param reduction.key
#' @param seed.use
#' @param ...
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
#' @method RunMDS Assay
#'
RunMDS.Assay <- function(object, assay = NULL, slot = "data",
                         features = NULL, nmds = 50, dist.method = "euclidean",
                         mds.method = "cmdscale", rev.mds = FALSE,
                         ndims.print = 1:5, nfeatures.print = 30,
                         reduction.key = "MDS_", verbose = TRUE, seed.use = 42, ...) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(object = object, slot = slot)
  features.var <- apply(
    X = data.use[features, ], MARGIN = 1,
    FUN = var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunMDS(
    object = data.use,
    assay = assay,
    slot = slot,
    nmds = nmds,
    dist.method = dist.method,
    mds.method = mds.method,
    rev.mds = rev.mds,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @param assay Name of Assay mds is being run on.
#'
#' @param slot Name of slot mds is being run on.
#' @param nmds Total Number of MDS to compute and store (50 by default).
#' @param rev.mds By default computes the mds on the cell x gene matrix. Setting
#' to true will compute it on gene x cell matrix.
#' @param verbose Print the top genes associated with high/low loadings for
#' the MDS.
#' @param ndims.print MDS to print genes for.
#' @param nfeatures.print Number of genes to print for each MDS.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "MDS_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param object
#' @param dist.method
#' @param mds.method
#' @param ...
#'
#' @importFrom parallelDist parDist
#' @importFrom utils capture.output
#' @importFrom Matrix t
#' @importFrom stats cmdscale
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#'
#' @rdname RunMDS
#' @concept dimensional_reduction
#' @export
#' @method RunMDS default
RunMDS.default <- function(object, assay = NULL, slot = "data",
                           nmds = 50, dist.method = "euclidean", mds.method = "cmdscale",
                           rev.mds = FALSE, ndims.print = 1:5, nfeatures.print = 30,
                           reduction.key = "MDS_", seed.use = 42, verbose = TRUE, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.mds) {
    object <- t(object)
  }
  nmds <- min(nmds, nrow(x = object) - 1)
  x <- t(as.matrix(object))
  dist.method.keep <- dist.method
  if (dist.method %in% c("pearson", "spearman")) {
    if (dist.method == "spearman") {
      x <- t(apply(x, 1, rank))
    }
    x <- t(apply(x, 1, function(x) x - mean(x)))
    dist.method <- "cosine"
  }
  cell.dist <- parDist(x = x, method = dist.method)
  if (mds.method == "cmdscale") {
    mds.results <- cmdscale(cell.dist, k = nmds, eig = TRUE)
  }
  if (mds.method == "isoMDS") {
    check_R("MASS")
    mds.results <- MASS::isoMDS(cell.dist, k = nmds)
  }
  if (mds.method == "sammon") {
    check_R("MASS")
    mds.results <- MASS::sammon(cell.dist, k = nmds)
  }
  cell.embeddings <- mds.results$points

  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:nmds)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dist.method = dist.method.keep, mds.method = mds.method)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = "\n"))
  }
  return(reduction.data)
}

#' Run GLMPCA
#'
#' @param object An object
#' @param ... Extra parameters passed to \code{\link[glmpca]{glmpca}}
#'
#' @return A Seurat object containing the output of GLMPCA stored as a DimReduc object.
#'
#' @export
#'
#' @rdname RunGLMPCA
#' @export RunGLMPCA
#'
RunGLMPCA <- function(object, ...) {
  UseMethod(generic = "RunGLMPCA", object = object)
}

#' @param object A Seurat object
#' @param L The number of dimensions to return (defaults to 5)
#' @param assay Assay to use, defaults to the default assay
#' @param features A list of features to use when performing GLM-PCA. If null, defaults to variable features.
#' @param reduction.name Name to store resulting DimReduc object as. Defaults to glmpca
#' @param reduction.key Key for resulting DimReduc. Defaults to GLMPC_
#' @param ... Extra parameters passed to \code{\link[glmpca]{glmpca}}
#' @param ...
#'
#' @rdname RunGLMPCA
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
#' @method RunGLMPCA Seurat
#'
RunGLMPCA.Seurat <- function(object,
                             L = 5,
                             fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                             assay = NULL,
                             slot = "counts",
                             features = NULL,
                             reduction.name = "glmpca",
                             reduction.key = "GLMPC_",
                             verbose = TRUE,
                             ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunGLMPCA(
    object = assay.data,
    assay = assay,
    slot = slot,
    L = L,
    fam = fam,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#' @param L
#' @param fam
#' @param assay
#' @param slot
#' @param features
#' @param reduction.name
#' @param reduction.key
#' @param verbose
#' @param ...
#'
#' @rdname RunGLMPCA
#' @concept dimensional_reduction
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
#' @method RunGLMPCA Assay
RunGLMPCA.Assay <- function(object,
                            L = 5,
                            fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                            assay = NULL,
                            slot = "counts",
                            features = NULL,
                            reduction.name = "glmpca",
                            reduction.key = "GLMPC_",
                            verbose = TRUE,
                            ...) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(object = object, slot = slot)
  features.var <- apply(
    X = data.use[features, ], MARGIN = 1,
    FUN = var
  )
  features.keep <- features[features.var > 0]
  data.use <- data.use[features.keep, ]
  reduction.data <- RunGLMPCA(
    object = data.use,
    assay = assay,
    slot = slot,
    L = L,
    fam = fam,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @param object
#' @param L
#' @param fam
#' @param assay
#' @param slot
#' @param features
#' @param reduction.name
#' @param reduction.key
#' @param verbose
#' @param ...
#'
#' @rdname RunGLMPCA
#' @concept dimensional_reduction
#' @importFrom Seurat DefaultAssay DefaultAssay<- CreateDimReducObject Tool<- LogSeuratCommand
#' @export
#' @method RunGLMPCA default
RunGLMPCA.default <- function(object,
                              L = 5,
                              fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                              assay = NULL,
                              slot = "counts",
                              features = NULL,
                              reduction.key = "GLMPC_",
                              verbose = TRUE,
                              ...) {
  check_R("glmpca")
  if (inherits(object, "dgCMatrix")) {
    object <- as_matrix(object)
  }
  glmpca_results <- glmpca::glmpca(Y = object, L = L, fam = fam, ...)
  glmpca_dimnames <- paste0(reduction.key, 1:L)
  factors <- as.matrix(glmpca_results$factors)
  loadings <- as.matrix(glmpca_results$loadings)
  colnames(x = factors) <- glmpca_dimnames
  colnames(x = loadings) <- glmpca_dimnames
  factors_l2_norm <- sqrt(colSums(factors^2))
  # strip S3 class "glmpca" to enable it to pass validObject()
  class(glmpca_results) <- NULL
  # save memory by removing factors and loadings since they are stored separately
  glmpca_results$factors <- glmpca_results$loadings <- NULL
  reduction <- CreateDimReducObject(
    embeddings = factors,
    key = reduction.key,
    loadings = loadings,
    stdev = factors_l2_norm,
    assay = assay,
    global = TRUE,
    misc = list(slot = slot, glmpca_results = glmpca_results)
  )
  return(reduction)
}

#' Run DiffusionMap(DM)
#'
#' Run a DM dimensionality reduction.
#'
#' @param object An object
#' @param ... Arguments passed to \link[destiny]{DiffusionMap}.
#'
#' @return Returns Seurat object with the DM calculation stored in the reductions slot.
#'
#' @export
#'
#' @rdname RunDM
#' @export RunDM
#'
RunDM <- function(object, ...) {
  UseMethod(generic = "RunDM", object = object)
}

#' @param object
#'
#' @param reduction
#' @param dims
#' @param assay
#' @param slot
#' @param features
#' @param ndcs
#' @param sigma
#' @param k
#' @param dist.method
#' @param reduction.key
#' @param reduction.name dimensional reduction name, 'dm' by default
#' @param seed.use
#' @param verbose
#' @param ...
#' @param graph
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
#' @method RunDM Seurat
#'
RunDM.Seurat <- function(object,
                         reduction = "pca", dims = 1:30, features = NULL, assay = NULL, slot = "data",
                         graph = NULL, ndcs = 2, sigma = "local", k = 30, dist.method = "euclidean",
                         reduction.name = "dm", reduction.key = "DM_", seed.use = 42, verbose = TRUE, ...) {
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < ndcs) {
      stop("Please provide as many or more features than ndcs: ",
        length(x = features), " features provided, ",
        ndcs, " Diffusion components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < ndcs) {
      stop("Please provide as many or more dims than ndcs: ",
        length(x = dims), " dims provided, ", ndcs,
        " DiffusionMap components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
  }
  reduction.data <- RunDM(
    object = data.use,
    assay = assay,
    slot = slot,
    ndcs = ndcs,
    sigma = sigma,
    k = k,
    dist.method = dist.method,
    reduction.key = reduction.key,
    seed.use = seed.use,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' @param object
#'
#' @param ndcs
#' @param sigma
#' @param k
#' @param dist.method
#' @param assay
#' @param slot
#' @param reduction.key
#' @param seed.use
#' @param verbose
#' @param ...
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @importFrom parallelDist parDist
#' @importFrom utils capture.output
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#' @export
#' @method RunDM matrix
#'
RunDM.matrix <- function(object, ndcs = 2, sigma = "local", k = 30, dist.method = "euclidean",
                         assay = NULL, slot = "data",
                         reduction.key = "DM_", seed.use = 42, verbose = TRUE, ...) {
  check_R("destiny")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ndcs <- min(ndcs, nrow(x = object) - 1)
  x <- as.matrix(object)
  dm.results <- destiny::DiffusionMap(data = x, n_eigs = ndcs, sigma = sigma, k = k, distance = dist.method, verbose = verbose)
  # set.seed(11);plot(dm.results)
  # knn <- get_knn(NULL, dists, k, distance, knn_params,verbose)
  # gene.relevance <- gene_relevance(coords=dm.results@eigenvectors,verbose=verbose)
  cell.embeddings <- dm.results@eigenvectors
  rownames(x = cell.embeddings) <- rownames(object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dm.results = dm.results)
  )
  return(reduction)
  # if (dist.method %in% c("pearson", "spearman")) {
  #   if (dist.method == "spearman") {
  #     x <- t(apply(x, 1, rank))
  #   }
  #   x <- t(apply(x, 1, function(x) x - mean(x)))
  #   dist.method <- "cosine"
  # }
  # cell.dist <- parDist(x = x, method = dist.method)
  # reduction <- RunDM(
  #   object = cell.dist,
  #   exprs = x,
  #   ndcs = ndcs,
  #   sigma = sigma,
  #   k = k,
  #   assay = assay,
  #   slot = slot,
  #   reduction.key = reduction.key,
  #   seed.use = seed.use,
  #   verbose = verbose,
  #   ...
  # )
  # return(reduction)
}

#' @param ndcs Total Number of DM to compute and store (2 by default).
#'
#' @param assay Name of Assay dm is being run on.
#' @param slot Name of slot dm is being run on.
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. "DM_" by default.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param verbose Show a progressbar and other progress information
#' @param object
#' @param exprs
#' @param sigma
#' @param k
#' @param ...
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @importFrom parallelDist parDist
#' @importFrom utils capture.output
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#' @export
RunDM.dist <- function(object,
                       ndcs = 2, sigma = "local", k = 30, slot = "data",
                       assay = NULL, reduction.key = "DM_", verbose = TRUE, seed.use = 42, ...) {
  check_R("destiny")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  # message("Number of nearest neighbors: ", k)
  dm.results <- destiny::DiffusionMap(distance = object, n_eigs = ndcs, sigma = sigma, k = k, verbose = verbose)
  set.seed(11)
  plot(dm.results)
  # knn <- get_knn(NULL, dists, k, distance, knn_params,verbose)
  # gene.relevance <- gene_relevance(coords=dm.results@eigenvectors,verbose=verbose)
  cell.embeddings <- dm.results@eigenvectors
  rownames(x = cell.embeddings) <- attr(object, "Labels")
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dm.results = dm.results)
  )
  return(reduction)
}

#' Run UMAP
#'
#' @param object
#' @param ...
#'
#' @export
#'
#' @rdname RunUMAP2
#' @export RunUMAP2
#'
RunUMAP2 <- function(object, ...) {
  UseMethod(generic = "RunUMAP2", object = object)
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param graph
#' @param assay
#' @param neighbor
#' @param slot
#' @param umap.method
#' @param reduction.model
#' @param return.model
#' @param n.neighbors
#' @param n.components
#' @param metric
#' @param n.epochs
#' @param learning.rate
#' @param min.dist
#' @param spread
#' @param set.op.mix.ratio
#' @param local.connectivity
#' @param repulsion.strength
#' @param negative.sample.rate
#' @param a
#' @param b
#' @param seed.use
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#'
#' @rdname RunUMAP2
#' @method RunUMAP2 Seurat
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand
#' @export
RunUMAP2.Seurat <- function(object,
                            reduction = "pca", dims = NULL, features = NULL, neighbor = NULL, graph = NULL,
                            assay = DefaultAssay(object = object), slot = "data",
                            umap.method = "uwot", reduction.model = NULL,
                            return.model = FALSE, n.neighbors = 30L, n.components = 2L,
                            metric = "cosine", n.epochs = 200L, spread = 1, min.dist = 0.3,
                            set.op.mix.ratio = 1, local.connectivity = 1L, negative.sample.rate = 5L,
                            a = NULL, b = NULL, learning.rate = 1, repulsion.strength = 1,
                            seed.use = 42L, verbose = TRUE,
                            reduction.name = "umap", reduction.key = "UMAP_",
                            ...) {
  if (sum(c(is.null(x = dims), is.null(x = features), is.null(neighbor), is.null(x = graph))) == 4) {
    stop("Please specify only one of the following arguments: dims, features, neighbor or graph")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n.components) {
      stop(
        "Please provide as many or more features than n.components: ",
        length(x = features),
        " features provided, ",
        n.components,
        " UMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- as.matrix(Embeddings(object[[reduction]])[, dims])
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n.components) {
      stop(
        "Please provide as many or more dims than n.components: ",
        length(x = dims),
        " dims provided, ",
        n.components,
        " UMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = neighbor)) {
    if (!inherits(x = object[[neighbor]], what = "Neighbor")) {
      stop(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[neighbor]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[neighbor]]
  } else if (!is.null(x = graph)) {
    if (!inherits(x = object[[graph]], what = "Graph")) {
      stop(
        "Please specify a Graph object name, ",
        "instead of the name of a ",
        class(object[[graph]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[graph]]
  } else {
    stop("Please specify one of dims, features, neighbor, or graph")
  }
  object[[reduction.name]] <- RunUMAP2(
    object = data.use, assay = assay,
    umap.method = umap.method, reduction.model = reduction.model,
    return.model = return.model, n.neighbors = n.neighbors, n.components = n.components,
    metric = metric, n.epochs = n.epochs, spread = spread, min.dist = min.dist,
    set.op.mix.ratio = set.op.mix.ratio, local.connectivity = local.connectivity, negative.sample.rate = negative.sample.rate,
    a = a, b = b, learning.rate = learning.rate, repulsion.strength = repulsion.strength,
    seed.use = seed.use, verbose = verbose, reduction.key = reduction.key
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#'
#' @param assay
#' @param umap.method
#' @param reduction.model
#' @param return.model
#' @param n.neighbors
#' @param n.components
#' @param metric
#' @param n.epochs
#' @param learning.rate
#' @param min.dist
#' @param spread
#' @param set.op.mix.ratio
#' @param local.connectivity
#' @param repulsion.strength
#' @param negative.sample.rate
#' @param a
#' @param b
#' @param seed.use
#' @param verbose
#' @param reduction.key
#' @param ...
#'
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#'
#' @rdname RunUMAP2
#' @method RunUMAP2 default
#' @concept dimensional_reduction
#' @importFrom SeuratObject Indices Distances as.sparse
#' @importFrom Matrix sparseMatrix
#' @export
RunUMAP2.default <- function(object, assay = NULL,
                             umap.method = "uwot", reduction.model = NULL,
                             return.model = FALSE, n.neighbors = 30L, n.components = 2L,
                             metric = "cosine", n.epochs = 200L, spread = 1, min.dist = 0.3,
                             set.op.mix.ratio = 1, local.connectivity = 1L, negative.sample.rate = 5L,
                             a = NULL, b = NULL, learning.rate = 1, repulsion.strength = 1,
                             seed.use = 42L, verbose = TRUE, reduction.key = "UMAP_", ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (return.model) {
    if (verbose) {
      message("UMAP will return its model")
    }
  }
  if (!is.null(x = reduction.model)) {
    if (verbose) {
      message("Running UMAP projection")
    }
    if (is.null(x = reduction.model) || !inherits(x = reduction.model, what = "DimReduc")) {
      stop("If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.",
        call. = FALSE
      )
    }
    model <- Misc(object = reduction.model, slot = "model")
    if (length(x = model) == 0) {
      stop("The provided reduction.model does not have a model stored.",
        call. = FALSE
      )
    }
    umap.method <- ifelse("layout" %in% names(model), "naive-predict", "uwot-predict")
  }
  n.epochs <- as.integer(n.epochs)
  n.neighbors <- as.integer(n.neighbors)
  n.components <- as.integer(n.components)
  local.connectivity <- as.integer(local.connectivity)
  negative.sample.rate <- as.integer(negative.sample.rate)

  if (inherits(x = object, what = "Neighbor")) {
    object <- list(idx = Indices(object), dist = Distances(object))
  }

  if (umap.method == "naive") {
    umap.config <- umap::umap.defaults
    umap.config$n_neighbors <- n.neighbors
    umap.config$n_components <- n.components
    umap.config$metric <- metric
    umap.config$n_epochs <- ifelse(is.null(n.epochs), 200, n.epochs)
    umap.config$spread <- spread
    umap.config$min_dist <- min.dist
    umap.config$set_op_mix_ratio <- set.op.mix.ratio
    umap.config$local_connectivity <- local.connectivity
    umap.config$a <- ifelse(is.null(a), NA, a)
    umap.config$b <- ifelse(is.null(b), NA, b)
    umap.config$gamma <- repulsion.strength
    umap.config$alpha <- learning.rate
    umap.config$negative_sample_rate <- negative.sample.rate
    umap.config$random_state <- seed.use
    umap.config$transform_state <- seed.use
    umap.config$verbose <- verbose
    if (is.na(umap.config$a) || is.na(umap.config$b)) {
      umap.config[c("a", "b")] <- umap:::find.ab.params(umap.config$spread, umap.config$min_dist)
      umap.config$min_dist <- umap::umap.defaults$min_dist
    }

    if (inherits(x = object, what = "dist")) {
      knn <- umap:::knn.from.dist(d = object, k = n.neighbors)
      out <- umap::umap(
        d = matrix(nrow = nrow(object)),
        config = umap.config, knn = knn
      )
      embeddings <- out$layout
      rownames(x = embeddings) <- attr(object, "Labels")
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "list")) {
      # object <- srt@neighbors$Seurat.nn
      # object <- list(idx = Indices(object), dist = Distances(object))

      # j <- as.numeric(t(object$idx))
      # i <- ((1:length(j)) - 1) %/% ncol(object$idx) + 1
      # graph <- as(object = sparseMatrix(
      #   i = i, j = j, x = as.numeric(t(object$dist)),
      #   dims = c(nrow(x = object$idx), nrow(x = object$idx))
      # ), Class = "Graph")
      # rownames(x = graph) <- rownames(object$idx)
      # colnames(x = graph) <- rownames(object$idx)
      # object <- graph
      knn <- umap::umap.knn(indexes = object[["idx"]], distances = object[["dist"]])
      out <- umap::umap(
        d = matrix(nrow = nrow(object[["idx"]])),
        config = umap.config, knn = knn
      )
      embeddings <- out$layout
      rownames(x = embeddings) <- rownames(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      if (!inherits(object, "dgCMatrix")) {
        object <- as.sparse(object[1:nrow(object), ])
      }
      diag(object) <- 0
      if (ncol(object) > 10000) {
        obs_sample <- sample(1:ncol(object), size = 10000)
      } else {
        obs_sample <- 1:ncol(object)
      }
      if (!isSymmetric(as.matrix(object[obs_sample, obs_sample]))) {
        stop("Graph must be a symmetric matrix.")
      }

      coo <- matrix(ncol = 3, nrow = length(object@x))
      coo[, 1] <- rep(1:ncol(object), diff(object@p))
      coo[, 2] <- object@i + 1
      coo[, 3] <- object@x
      colnames(coo) <- c("from", "to", "value")
      coo <- coo[order(coo[, 1], coo[, 2]), ]
      graph <- list(coo = coo, names = rownames(object), n.elements = nrow(object))
      class(graph) <- "coo"

      # diag(object) <- 0
      # if (inherits(object, what = "dgCMatrix")) {
      #   object <- as_matrix(object)
      # } else {
      #   object <- as.matrix(object)
      # }
      # if (!isSymmetric(object)) {
      #   stop("Graph must be a symmetric matrix.")
      # }
      # graph <- umap:::coo(object)
      # # coo <- function(x) {
      # #   if (!is(x, "Matrix")) {
      # #     stop("x must be a square matrix\n")
      # #   }
      # #   if (nrow(x) != ncol(x)) {
      # #     stop("x must be a square matrix\n")
      # #   }
      # #   nx <- nrow(x)
      # #   coo <- matrix(0, ncol = 3, nrow = nx * nx)
      # #   coo[, 1] <- rep(seq_len(nx), nx)
      # #   coo[, 2] <- rep(seq_len(nx), each = nx)
      # #   coo[, 3] <- as.vector(x)
      # #   colnames(coo) <- c("from", "to", "value")
      # #   coo <- coo[order(coo[, 1], coo[, 2]), ]
      # #   make.coo(coo, rownames(x), nrow(x))
      # # }
      # # make.coo <- function(x, names, n.elements) {
      # #   x <- x[, 1:3, drop = FALSE]
      # #   colnames(x) <- c("from", "to", "value")
      # #   result <- list(coo = x, names = names, n.elements = n.elements)
      # #   class(result) <- "coo"
      # #   result
      # # }

      umap.config$init <- "spectral"
      initial <- umap:::make.initial.embedding(V = graph$n.elements, config = umap.config, g = graph)
      embeddings <- umap:::naive.simplicial.set.embedding(g = graph, embedding = initial, config = umap.config)
      embeddings <- umap:::center.embedding(embeddings)
      rownames(x = embeddings) <- rownames(x = object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        warning("return.model does not support 'Graph' input.", immediate. = TRUE)
      }
      return(reduction)
    }
    if (inherits(x = object, what = "matrix") || inherits(x = object, what = "Matrix")) {
      # object <- srt@reductions$Seuratpca@cell.embeddings
      out <- umap::umap(d = object, config = umap.config, method = "naive")
      embeddings <- out$layout
      rownames(x = embeddings) <- rownames(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
  }

  if (umap.method == "uwot") {
    if (inherits(x = object, what = "dist")) {
      embeddings <- uwot::umap(
        X = object, n_neighbors = n.neighbors, n_threads = 1, n_components = n.components,
        metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
        min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a, b = b, verbose = verbose,
        ret_model = FALSE
      )
      rownames(x = embeddings) <- attr(object, "Labels")
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        warning("return.model does not support 'dist' input.", immediate. = TRUE)
      }
      return(reduction)
    }
    if (inherits(x = object, what = "list")) {
      # object <- srt@neighbors$Seurat.nn
      # object <- list(idx = Indices(object), dist = Distances(object))
      out <- uwot::umap(
        X = NULL, nn_method = object, n_threads = 1, n_components = n.components,
        metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
        min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a, b = b, verbose = verbose,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      if (!inherits(object, "dgCMatrix")) {
        object <- as.sparse(object[1:nrow(object), ])
      }
      diag(object) <- 0
      if (ncol(object) > 10000) {
        obs_sample <- sample(1:ncol(object), size = 10000)
      } else {
        obs_sample <- 1:ncol(object)
      }
      if (!isSymmetric(as.matrix(object[obs_sample, obs_sample]))) {
        stop("Graph must be a symmetric matrix.")
      }
      val <- split(object@x, rep(1:ncol(object), diff(object@p)))
      pos <- split(object@i + 1, rep(1:ncol(object), diff(object@p)))
      idx <- t(mapply(function(x, y) {
        out <- y[head(order(x, decreasing = TRUE), n.neighbors)]
        length(out) <- n.neighbors
        return(out)
      }, x = val, y = pos))
      connectivity <- t(mapply(function(x, y) {
        out <- y[head(order(x, decreasing = TRUE), n.neighbors)]
        length(out) <- n.neighbors
        out[is.na(out)] <- 0
        return(out)
      }, x = val, y = val))
      idx[is.na(idx)] <- sample(1:nrow(object), size = sum(connectivity == 0), replace = TRUE)
      nn <- list(idx = idx, dist = max(connectivity) - connectivity + min(diff(range(connectivity)), 1) / 1e50)
      # idx <- t(as.matrix(apply(object, 2, function(x) order(x, decreasing = TRUE)[1:n.neighbors])))
      # connectivity <- t(as.matrix(apply(object, 2, function(x) x[order(x, decreasing = TRUE)[1:n.neighbors]])))
      out <- uwot::umap(
        X = NULL, nn_method = nn, n_threads = 1, n_components = n.components,
        metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
        min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a, b = b, verbose = verbose,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }
      # object <- srt@graphs$BBKNN
      # object@x <- max(object@x) - object@x + 1e-10
      # min_neighbors <- min(rowSums(object > 0))
      # embeddings <- uwot::umap(
      #   X = object*100, n_neighbors = min(min_neighbors,n.neighbors), n_threads = 1, n_components = n.components,
      #   metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
      #   min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
      #   local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
      #   negative_sample_rate = negative.sample.rate,
      #   a = a, b = b, verbose = verbose,
      #   ret_model = FALSE
      # )
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
    if (inherits(x = object, what = "matrix") || inherits(x = object, what = "Matrix")) {
      # object <- srt@reductions$Seuratpca@cell.embeddings
      out <- uwot::umap(
        X = object, n_neighbors = n.neighbors, n_threads = 1, n_components = n.components,
        metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
        min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate,
        a = a, b = b, verbose = verbose,
        ret_model = return.model
      )
      if (return.model) {
        embeddings <- out$embedding
      } else {
        embeddings <- out
      }
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      if (return.model) {
        out$nn_index <- NULL
        Misc(reduction, slot = "model") <- out
      }
      return(reduction)
    }
  }

  if (umap.method == "naive-predict") {
    if (inherits(x = object, what = "matrix") || inherits(x = object, what = "Matrix")) {
      class(model) <- "umap"
      if (any(!colnames(model$data) %in% colnames(object))) {
        stop("query data must contain the same features with the model:\n", paste(head(colnames(model$data), 10), collapse = ","), " ......")
      }
      embeddings <- predict(model, object[, colnames(model$data)])
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    } else {
      stop("naive umap model only support 'matrix' input.")
    }
  }

  if (umap.method == "uwot-predict") {
    if (inherits(x = object, what = "list")) {
      if (ncol(object[["idx"]]) != model$n_neighbors) {
        warning(
          "Number of neighbors between query and reference is not equal to the number of neighbros within reference"
        )
        model$n_neighbors <- ncol(object[["idx"]])
      }
      if (is.null(model$num_precomputed_nns) || model$num_precomputed_nns == 0) {
        model$num_precomputed_nns <- 1
      }
      embeddings <- uwot::umap_transform(
        X = NULL, nn_method = object, model = model, n_epochs = n.epochs,
        n_threads = 1, verbose = verbose
      )
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
    if (inherits(x = object, what = "Graph")) {
      match_k <- t(as.matrix(apply(object, 2, function(x) order(x, decreasing = TRUE)[1:n.neighbors])))
      match_k_connectivity <- t(as.matrix(apply(object, 2, function(x) x[order(x, decreasing = TRUE)[1:n.neighbors]])))
      object <- list(idx = match_k, dist = max(match_k_connectivity) - match_k_connectivity)
      if (is.null(model$num_precomputed_nns) || model$num_precomputed_nns == 0) {
        model$num_precomputed_nns <- 1
      }
      embeddings <- uwot::umap_transform(
        X = NULL, nn_method = object, model = model, n_epochs = n.epochs,
        n_threads = 1, verbose = verbose
      )
      rownames(x = embeddings) <- row.names(object[["idx"]])
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
    if (inherits(x = object, what = "matrix") || inherits(x = object, what = "Matrix")) {
      embeddings <- uwot::umap_transform(
        X = object, model = model, n_epochs = n.epochs,
        n_threads = 1, verbose = verbose
      )
      rownames(x = embeddings) <- row.names(object)
      colnames(x = embeddings) <- paste0(reduction.key, 1:n.components)
      reduction <- CreateDimReducObject(
        embeddings = embeddings,
        key = reduction.key,
        assay = assay,
        global = TRUE
      )
      return(reduction)
    }
  }
}


#' Run PaCMAP
#'
#' @param object
#' @param ...
#'
#' @export
#'
#' @rdname RunPaCMAP
#' @export RunPaCMAP
#'
RunPaCMAP <- function(object, ...) {
  UseMethod(generic = "RunPaCMAP", object = object)
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param assay
#' @param slot
#' @param n.neighbors
#' @param seed.use
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#' @param MN_ratio
#' @param FP_ratio
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param distance_method
#' @param lr
#' @param num_iters
#' @param apply_pca
#' @param init
#' @param n_components
#'
#' @rdname RunPaCMAP
#' @method RunPaCMAP Seurat
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand
#' @export
RunPaCMAP.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                             assay = DefaultAssay(object = object), slot = "data",
                             n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                             pair_neighbors = NULL, pair_MN = NULL, pair_FP = NULL, distance_method = "euclidean",
                             lr = 1, num_iters = 450L, apply_pca = TRUE, init = "random",
                             reduction.name = "pacmap", reduction.key = "PaCMAP_",
                             verbose = TRUE, seed.use = 11L,
                             ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) < 1) {
    stop("Please specify only one of the following arguments: dims, features, or graph")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " PaCMAP components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " PaCMAP components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims or features")
  }
  object[[reduction.name]] <- RunPaCMAP(
    object = data.use, assay = assay,
    n_components = n_components, n.neighbors = n.neighbors, MN_ratio = MN_ratio, FP_ratio = FP_ratio,
    pair_neighbors = pair_neighbors, pair_MN = pair_MN, pair_FP = pair_FP, distance_method = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca, init = init,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#'
#' @param n.neighbors
#' @param n_components
#' @param seed.use
#' @param verbose
#' @param reduction.key
#' @param ...
#' @param MN_ratio
#' @param FP_ratio
#' @param pair_neighbors
#' @param pair_MN
#' @param pair_FP
#' @param distance_method
#' @param lr
#' @param num_iters
#' @param apply_pca
#' @param init
#' @param assay
#'
#' @rdname RunPaCMAP
#' @concept dimensional_reduction
#' @importFrom Seurat CreateDimReducObject
#' @importFrom reticulate import
#' @export
#' @method RunPaCMAP default
RunPaCMAP.default <- function(object, assay = NULL,
                              n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                              pair_neighbors = NULL, pair_MN = NULL, pair_FP = NULL, distance_method = "euclidean",
                              lr = 1, num_iters = 450L, apply_pca = TRUE, init = "random",
                              reduction.key = "PaCMAP_", verbose = TRUE, seed.use = 11L,
                              ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_Python("pacmap")
  pacmap <- import("pacmap")

  operator <- pacmap$PaCMAP(
    n_components = as.integer(n_components), n_neighbors = n.neighbors, MN_ratio = MN_ratio, FP_ratio = FP_ratio,
    pair_neighbors = pair_neighbors, pair_MN = pair_MN, pair_FP = pair_FP, distance = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca,
    verbose = verbose, random_state = as.integer(seed.use)
  )
  embedding <- operator$fit_transform(object, init = init)

  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key, assay = assay, global = TRUE
  )
  return(reduction)
}


#' Run PHATE
#'
#' @param object
#' @param ...
#'
#' @export
#'
#' @rdname RunPHATE
#' @export RunPHATE
#'
RunPHATE <- function(object, ...) {
  UseMethod(generic = "RunPHATE", object = object)
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param assay
#' @param slot
#' @param seed.use
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param n_components
#' @param knn
#' @param decay
#' @param n_landmark
#' @param t
#' @param gamma
#' @param n_pca
#' @param knn_dist
#' @param knn_max
#' @param mds
#' @param mds_dist
#' @param mds_solver
#' @param t_max
#' @param n_jobs
#' @param ...
#'
#' @rdname RunPHATE
#' @concept dimensional_reduction
#' @export
#' @importFrom Seurat LogSeuratCommand
#' @method RunPHATE Seurat
RunPHATE.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                            assay = DefaultAssay(object = object), slot = "data",
                            n_components = 2, knn = 5, decay = 40, n_landmark = 2000, t = "auto", gamma = 1,
                            n_pca = 100, knn_dist = "euclidean", knn_max = NULL, n_jobs = 1,
                            mds = "metric", mds_dist = "euclidean", mds_solver = "sgd", t_max = 100,
                            do_cluster = FALSE, n_clusters = "auto", max_clusters = 100,
                            reduction.name = "phate", reduction.key = "PHATE_",
                            verbose = TRUE, seed.use = 11L,
                            ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " PHATE components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " PHATE components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims or features")
  }
  object[[reduction.name]] <- RunPHATE(
    object = data.use, assay = assay,
    n_components = n_components, knn = knn, decay = decay, n_landmark = n_landmark, t = t, gamma = gamma,
    n_pca = n_pca, knn_dist = knn_dist, knn_max = knn_max, n_jobs = n_jobs,
    mds = mds, mds_dist = mds_dist, mds_solver = mds_solver, t_max = t_max,
    do_cluster = do_cluster, n_clusters = n_clusters, max_clusters = max_clusters,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#'
#' @param seed.use
#' @param verbose
#' @param reduction.key
#' @param ...
#' @param n_components
#' @param knn
#' @param decay
#' @param n_landmark
#' @param t
#' @param gamma
#' @param n_pca
#' @param knn_dist
#' @param knn_max
#' @param mds
#' @param mds_dist
#' @param mds_solver
#' @param t_max
#' @param n_jobs
#' @param assay
#' @param n_clusters
#' @param max_clusters
#'
#' @importFrom reticulate import
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#'
#' @rdname RunPHATE
#' @concept dimensional_reduction
#' @importFrom reticulate import
#' @export
#' @method RunPHATE default
RunPHATE.default <- function(object, assay = NULL,
                             n_components = 2, knn = 5, decay = 40, n_landmark = 2000, t = "auto", gamma = 1,
                             n_pca = 100, knn_dist = "euclidean", knn_max = NULL, n_jobs = 1,
                             mds = "metric", mds_dist = "euclidean", mds_solver = "sgd", t_max = 100,
                             do_cluster = FALSE, n_clusters = "auto", max_clusters = 100,
                             reduction.key = "PHATE_", verbose = TRUE, seed.use = 11L,
                             ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_Python("phate")
  phate <- import("phate")

  if (is.numeric(knn_max) && length(knn_max) > 0) {
    knn_max <- as.integer(knn_max)
  } else {
    knn_max <- NULL
  }
  operator <- phate$PHATE(
    n_components = as.integer(n_components),
    knn = as.integer(knn),
    decay = as.integer(decay),
    n_landmark = as.integer(n_landmark),
    t = as.character(t),
    gamma = as.numeric(gamma),
    n_pca = as.integer(n_pca),
    knn_dist = as.character(knn_dist),
    knn_max = knn_max,
    mds = as.character(mds),
    mds_dist = as.character(mds_dist),
    mds_solver = as.character(mds_solver),
    n_jobs = as.integer(n_jobs),
    random_state = as.integer(seed.use),
    verbose = as.integer(verbose)
  )
  embedding <- operator$fit_transform(object, t_max = as.integer(t_max))
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key, assay = assay, global = TRUE
  )
  if (isTRUE(do_cluster)) {
    if (is.numeric(n_clusters)) {
      n_clusters <- as.integer(n_clusters)
    }
    if (is.numeric(max_clusters)) {
      max_clusters <- as.integer(max_clusters)
    }
    clusters <- phate$cluster$kmeans(operator, n_clusters = n_clusters, max_clusters = max_clusters, random_state = as.integer(seed.use))
    clusters <- clusters + 1
    clusters <- factor(clusters, levels = sort(unique(clusters)))
    names(clusters) <- rownames(embedding)
    Misc(reduction, slot = "clusters") <- clusters
  }
  return(reduction)
}


#' Run TriMap
#'
#' @param object
#' @param ...
#'
#' @export
#'
#' @rdname RunTriMap
#' @export RunTriMap
#'
RunTriMap <- function(object, ...) {
  UseMethod(generic = "RunTriMap", object = object)
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param assay
#' @param slot
#' @param seed.use
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#' @param n_inliers
#' @param n_outliers
#' @param n_random
#' @param distance
#' @param lr
#' @param n_iters
#' @param triplets
#' @param weights
#' @param use_dist_matrix
#' @param knn_tuple
#' @param apply_pca
#' @param opt_method
#' @param n_components
#'
#' @rdname RunTriMap
#' @concept dimensional_reduction
#' @export
#' @importFrom Seurat LogSeuratCommand
#' @method RunTriMap Seurat
RunTriMap.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                             assay = DefaultAssay(object = object), slot = "data",
                             n_components = 2, n_inliers = 12, n_outliers = 4, n_random = 3, distance_method = "euclidean",
                             lr = 0.1, n_iters = 400, triplets = NULL, weights = NULL, use_dist_matrix = FALSE, knn_tuple = NULL,
                             apply_pca = TRUE, opt_method = "dbd",
                             reduction.name = "trimap", reduction.key = "TriMap_",
                             verbose = TRUE, seed.use = 11L,
                             ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " TriMap components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " TriMap components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
  }
  object[[reduction.name]] <- RunTriMap(
    object = data.use, assay = assay,
    n_components = n_components, n_inliers = n_inliers, n_outliers = n_outliers, n_random = n_random, distance_method = distance_method,
    lr = lr, n_iters = n_iters, triplets = triplets, weights = weights, use_dist_matrix = use_dist_matrix, knn_tuple = knn_tuple,
    apply_pca = apply_pca, opt_method = opt_method,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#'
#' @param seed.use
#' @param verbose
#' @param reduction.key
#' @param ...
#' @param assay
#' @param n_inliers
#' @param n_outliers
#' @param n_random
#' @param distance
#' @param lr
#' @param n_iters
#' @param triplets
#' @param weights
#' @param use_dist_matrix
#' @param knn_tuple
#' @param apply_pca
#' @param opt_method
#' @param n_components
#'
#' @importFrom reticulate import
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#'
#' @rdname RunTriMap
#' @concept dimensional_reduction
#' @importFrom reticulate import
#' @export
#' @method RunTriMap default
RunTriMap.default <- function(object, assay = NULL,
                              n_components = 2, n_inliers = 12, n_outliers = 4, n_random = 3, distance_method = "euclidean",
                              lr = 0.1, n_iters = 400, triplets = NULL, weights = NULL, use_dist_matrix = FALSE, knn_tuple = NULL,
                              apply_pca = TRUE, opt_method = "dbd",
                              reduction.key = "TriMap_", verbose = TRUE, seed.use = 11L,
                              ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_Python("trimap")
  trimap <- import("trimap")

  if (inherits(x = object, what = "dist")) {
    object <- as.matrix(object)
    use_dist_matrix <- TRUE
  }
  operator <- trimap$TRIMAP(
    n_dims = as.integer(n_components), n_inliers = as.integer(n_inliers), n_outliers = as.integer(n_outliers), n_random = as.integer(n_random), distance = distance_method,
    lr = lr, n_iters = as.integer(n_iters), triplets = triplets, weights = weights, use_dist_matrix = use_dist_matrix, knn_tuple = knn_tuple,
    apply_pca = apply_pca, opt_method = opt_method, verbose = verbose
  )
  embedding <- operator$fit_transform(object)
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  if (inherits(x = object, what = "dist")) {
    rownames(x = embedding) <- attr(object, "Labels")
  } else {
    rownames(x = embedding) <- rownames(object)
  }
  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key, assay = assay, global = TRUE
  )
  return(reduction)
}


#' Run LargeVis
#'
#' @param object
#' @param ...
#'
#' @rdname RunLargeVis
#' @export RunLargeVis
#'
RunLargeVis <- function(object, ...) {
  UseMethod(generic = "RunLargeVis", object = object)
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param assay
#' @param slot
#' @param seed.use
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#' @param init
#' @param perplexity
#' @param n_neighbors
#' @param metric
#' @param n_epochs
#' @param learning_rate
#' @param scale
#' @param init_sdev
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param nn_method
#' @param n_trees
#' @param search_k
#' @param n_threads
#' @param n_sgd_threads
#' @param grain_size
#' @param kernel
#' @param pca
#' @param pca_center
#' @param pcg_rand
#' @param fast_sgd
#' @param batch
#' @param opt_args
#' @param epoch_callback
#' @param pca_method
#' @param n_components
#'
#' @rdname RunLargeVis
#' @method RunLargeVis Seurat
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand
#' @export
RunLargeVis.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                               assay = DefaultAssay(object = object), slot = "data",
                               perplexity = 50, n_neighbors = perplexity * 3, n_components = 2, metric = "euclidean",
                               n_epochs = -1, learning_rate = 1, scale = "maxabs", init = "lvrandom", init_sdev = NULL,
                               repulsion_strength = 7, negative_sample_rate = 5, nn_method = NULL, n_trees = 50,
                               search_k = 2 * n_neighbors * n_trees, n_threads = NULL, n_sgd_threads = 0, grain_size = 1,
                               kernel = "gauss", pca = NULL, pca_center = TRUE, pcg_rand = TRUE, fast_sgd = FALSE,
                               batch = FALSE, opt_args = NULL, epoch_callback = NULL, pca_method = NULL,
                               reduction.name = "largevis", reduction.key = "LargeVis_",
                               verbose = TRUE, seed.use = 11L,
                               ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 3) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
    if (ncol(x = data.use) < n_components) {
      stop(
        "Please provide as many or more features than n_components: ",
        length(x = features),
        " features provided, ",
        n_components,
        " LargeVis components requested",
        call. = FALSE
      )
    }
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    if (length(x = dims) < n_components) {
      stop(
        "Please provide as many or more dims than n_components: ",
        length(x = dims),
        " dims provided, ",
        n_components,
        " LargeVis components requested",
        call. = FALSE
      )
    }
  } else {
    stop("Please specify one of dims, features")
  }
  object[[reduction.name]] <- RunLargeVis(
    object = data.use, assay = assay,
    perplexity = perplexity, n_neighbors = n_neighbors, n_components = n_components, metric = metric,
    n_epochs = n_epochs, learning_rate = learning_rate, scale = scale, init = init, init_sdev = init_sdev,
    repulsion_strength = repulsion_strength, negative_sample_rate = negative_sample_rate, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k, n_threads = n_threads, n_sgd_threads = n_sgd_threads, grain_size = grain_size,
    kernel = kernel, pca = pca, pca_center = pca_center, pcg_rand = pcg_rand, fast_sgd = fast_sgd,
    batch = batch, opt_args = opt_args, epoch_callback = epoch_callback, pca_method = pca_method,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param object
#'
#' @param seed.use
#' @param verbose
#' @param reduction.key
#' @param ...
#' @param init
#' @param assay
#' @param perplexity
#' @param n_neighbors
#' @param metric
#' @param n_epochs
#' @param learning_rate
#' @param scale
#' @param init_sdev
#' @param repulsion_strength
#' @param negative_sample_rate
#' @param nn_method
#' @param n_trees
#' @param search_k
#' @param n_threads
#' @param n_sgd_threads
#' @param grain_size
#' @param kernel
#' @param pca
#' @param pca_center
#' @param pcg_rand
#' @param fast_sgd
#' @param batch
#' @param opt_args
#' @param epoch_callback
#' @param pca_method
#' @param n_components
#'
#' @rdname RunLargeVis
#' @concept dimensional_reduction
#' @importFrom Seurat CreateDimReducObject
#' @export
#' @method RunLargeVis default
RunLargeVis.default <- function(object, assay = NULL,
                                perplexity = 50, n_neighbors = perplexity * 3, n_components = 2, metric = "euclidean",
                                n_epochs = -1, learning_rate = 1, scale = "maxabs", init = "lvrandom", init_sdev = NULL,
                                repulsion_strength = 7, negative_sample_rate = 5, nn_method = NULL, n_trees = 50,
                                search_k = 2 * n_neighbors * n_trees, n_threads = NULL, n_sgd_threads = 0, grain_size = 1,
                                kernel = "gauss", pca = NULL, pca_center = TRUE, pcg_rand = TRUE, fast_sgd = FALSE,
                                batch = FALSE, opt_args = NULL, epoch_callback = NULL, pca_method = NULL,
                                reduction.key = "LargeVis_", verbose = TRUE, seed.use = 11L,
                                ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  embedding <- uwot::lvish(
    X = object, perplexity = perplexity, n_neighbors = n_neighbors, n_components = n_components, metric = metric,
    n_epochs = n_epochs, learning_rate = learning_rate, scale = scale, init = init, init_sdev = init_sdev,
    repulsion_strength = repulsion_strength, negative_sample_rate = negative_sample_rate, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k, n_threads = n_threads, n_sgd_threads = n_sgd_threads, grain_size = grain_size,
    kernel = kernel, pca = pca, pca_center = pca_center, pcg_rand = pcg_rand, fast_sgd = fast_sgd,
    verbose = verbose, batch = batch, opt_args = opt_args, epoch_callback = epoch_callback, pca_method = pca_method
  )
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  if (inherits(x = object, what = "dist")) {
    rownames(x = embedding) <- attr(object, "Labels")
  } else {
    rownames(x = embedding) <- rownames(object)
  }
  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key, assay = assay, global = TRUE
  )
  return(reduction)
}

# RunIsomap <- function(object, ...) {
#   UseMethod(generic = "RunIsomap", object = object)
# }
#
# RunIsomap.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL, distance = NULL,
#                              assay = DefaultAssay(object = object), slot = "data",
#                              n_components = 2, epsilon = NULL, k = 30, dist.method = "euclidean",
#                              reduction.name = "Isomap", reduction.key = "Isomap_",
#                              verbose = TRUE, seed.use = 11L,
#                              ...) {
#   if (sum(c(is.null(x = dims), is.null(x = features), is.null(distance))) == 3) {
#     stop("Please specify only one of the following arguments: dims, features or distance")
#   }
#   if (!is.null(x = features)) {
#     data.use <- as.matrix(x = t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
#     if (ncol(x = data.use) < n_components) {
#       stop(
#         "Please provide as many or more features than n_components: ",
#         length(x = features),
#         " features provided, ",
#         n_components,
#         " Isomap components requested",
#         call. = FALSE
#       )
#     }
#   } else if (!is.null(x = dims)) {
#     data.use <- Embeddings(object[[reduction]])[, dims]
#     assay <- DefaultAssay(object = object[[reduction]])
#     if (length(x = dims) < n_components) {
#       stop(
#         "Please provide as many or more dims than n_components: ",
#         length(x = dims),
#         " dims provided, ",
#         n_components,
#         " Isomap components requested",
#         call. = FALSE
#       )
#     }
#   } else if (!is.null(x = distance)) {
#     if (!inherits(x = object[[distance]], what = "Graph")) {
#       stop(
#         "Please specify a Graph object name(save the distance matrix), ",
#         "instead of the name of a ",
#         class(object[[distance]]),
#         " object",
#         call. = FALSE
#       )
#     }
#     data.use <- as.dist(as.matrix(object[[distance]]))
#   } else {
#     stop("Please specify one of dims, features or distance")
#   }
#   object[[reduction.name]] <- RunIsomap(
#     object = data.use, assay = assay,
#     n_components = n_components, epsilon = epsilon, k = k, dist.method = dist.method,
#     reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
#   )
#   object <- LogSeuratCommand(object = object)
#   return(object)
# }
#
# RunIsomap.default <- function(object, assay = NULL,
#                               n_components = 2, epsilon = NULL, k = 30, dist.method = "euclidean",
#                               reduction.key = "Isomap_", verbose = TRUE, seed.use = 11L,
#                               ...) {
#   if (!is.null(x = seed.use)) {
#     set.seed(seed = seed.use)
#   }
#   check_R("vegan")
#   if (inherits(x = object, what = "matrix") || inherits(x = object, what = "Matrix")) {
#     x <- as.matrix(object)
#     if (dist.method %in% c("pearson", "spearman")) {
#       if (dist.method == "spearman") {
#         x <- t(apply(x, 1, rank))
#       }
#       x <- t(apply(x, 1, function(x) x - mean(x)))
#       dist.method <- "cosine"
#     }
#     cell.dist <- parDist(x = x, method = dist.method)
#   } else if (inherits(x = object, what = "dist")) {
#     cell.dist <- object
#   } else {
#     stop("object must be a matrix or dist object")
#   }
#
#   params <- list(
#     dist = cell.dist,
#     ndim = n_components
#   )
#   for (nm in c("epsilon", "k")) {
#     params[[nm]] <- get(nm)
#   }
#   out <- invoke(.fn = vegan::isomap, .args = params)
#   embeddings <- out$points
#   colnames(x = embeddings) <- paste0(reduction.key, seq_len(ncol(x = embeddings)))
#   rownames(x = embeddings) <- attr(cell.dist, "Labels")
#   reduction <- CreateDimReducObject(
#     embeddings = embeddings,
#     key = reduction.key, assay = assay, global = TRUE
#   )
#   return(reduction)
# }

#' Harmony single cell integration
#'
#' Run Harmony algorithm with Seurat pipelines.
#'
#' @param object An Seurat object
#'
#' @rdname RunHarmony2
#' @export RunHarmony2
#'
RunHarmony2 <- function(object, ...) {
  UseMethod(generic = "RunHarmony2", object = object)
}

#' @param object
#'
#' @param group.by.vars
#' @param reduction
#' @param dims.use
#' @param theta
#' @param lambda
#' @param sigma
#' @param nclust
#' @param tau
#' @param block.size
#' @param max.iter.harmony
#' @param max.iter.cluster
#' @param epsilon.cluster
#' @param epsilon.harmony
#' @param plot_convergence
#' @param verbose
#' @param reference_values
#' @param reduction.save
#' @param assay.use
#' @param project.dim
#' @param ...
#'
#' @rdname RunHarmony2
#' @method RunHarmony2 Seurat
#' @importFrom Seurat Embeddings RunPCA FetchData CreateDimReducObject ProjectDim LogSeuratCommand
#' @export
RunHarmony2.Seurat <- function(object,
                               group.by.vars,
                               reduction = "pca",
                               dims.use = NULL,
                               theta = NULL,
                               lambda = NULL,
                               sigma = 0.1,
                               nclust = NULL,
                               tau = 0,
                               block.size = 0.05,
                               max.iter.harmony = 10,
                               max.iter.cluster = 20,
                               epsilon.cluster = 1e-5,
                               epsilon.harmony = 1e-4,
                               plot_convergence = FALSE,
                               verbose = TRUE,
                               reference_values = NULL,
                               reduction.save = "Harmony",
                               assay.use = "RNA",
                               project.dim = TRUE,
                               ...) {
  check_R("immunogenomics/harmony")
  available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
  if (!(reduction %in% available.dimreduc)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  embedding <- embedding[, dims.use, drop = FALSE]
  metavars_df <- FetchData(object, group.by.vars)

  harmonyObject <- harmony::HarmonyMatrix(
    data_mat = embedding, meta_data = metavars_df, vars_use = group.by.vars,
    do_pca = FALSE, npcs = 0, theta = theta, lambda = lambda, sigma = sigma,
    nclust = nclust, tau = tau, block.size = block.size,
    max.iter.harmony = max.iter.harmony, max.iter.cluster = max.iter.cluster,
    epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony,
    plot_convergence = plot_convergence, return_object = TRUE, verbose = verbose,
    reference_values = reference_values
  )

  harmonyEmbed <- t(as.matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  harmonyClusters <- t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(embedding)
  colnames(harmonyClusters) <- paste0("R", seq_len(ncol(harmonyClusters)))

  suppressWarnings({
    object[[reduction.save]] <- CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use,
      key = reduction.save,
      misc = list(
        R = harmonyClusters,
        reduction_use = reduction,
        reduction_dims = dims.use
      )
    )
  })

  if (project.dim) {
    object <- ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}
