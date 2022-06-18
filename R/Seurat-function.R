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
#' @importFrom RcppML nmf mse
#' @importFrom utils capture.output
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#'
#' @rdname RunNMF
#' @concept dimensional_reduction
#' @export
#'
RunNMF.default <- function(object,
                           assay = NULL,
                           slot = "data",
                           nbes = 50,
                           nmf.method = "RcppML",
                           tol = 1e-5,
                           maxit = 100,
                           rev.nmf = FALSE,
                           verbose = TRUE,
                           ndims.print = 1:5,
                           nfeatures.print = 30,
                           reduction.key = "BE_",
                           seed.use = 42,
                           ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.nmf) {
    object <- t(object)
  }
  nbes <- min(nbes, nrow(x = object) - 1)
  if (nmf.method == "RcppML") {
    nmf.results <- nmf(
      A = t(object), k = nbes, tol = tol, maxit = maxit,
      seed = seed.use, verbose = verbose, ...
    )
    cell.embeddings <- nmf.results$w
    feature.loadings <- t(nmf.results$h)
    d <- nmf.results$d
    iter <- nmf.results$iter
    mse <- mse(t(object), cell.embeddings, nmf.results$d, t(feature.loadings))
  }
  if (nmf.method == "NMF") {
    check_R("NMF")
    nmf.results <- NMF::nmf(x = as.matrix(t(object)), rank = nbes, seed = seed.use)
    cell.embeddings <- nmf.results@fit@W
    feature.loadings <- t(nmf.results@fit@H)
    d <- iter <- tol <- NULL
    mse <- mse(t(object), cell.embeddings, NULL, t(feature.loadings))
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
    misc = list(slot = slot, mse = mse, d = d, tol = tol, iter = iter)
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
RunNMF.Assay <- function(object,
                         assay = NULL,
                         slot = "data",
                         features = NULL,
                         nbes = 50,
                         nmf.method = "RcppML",
                         tol = 1e-5,
                         maxit = 100,
                         rev.nmf = FALSE,
                         verbose = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.key = "BE_",
                         seed.use = 42,
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
RunNMF.Seurat <- function(object,
                          assay = NULL,
                          slot = "data",
                          features = NULL,
                          nbes = 50,
                          nmf.method = "RcppML",
                          tol = 1e-5,
                          maxit = 100,
                          rev.nmf = FALSE,
                          verbose = TRUE,
                          ndims.print = 1:5,
                          nfeatures.print = 30,
                          reduction.name = "nmf",
                          reduction.key = "BE_",
                          seed.use = 42,
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
#'
RunMDS.default <- function(object,
                           assay = NULL,
                           slot = "scale.data",
                           nmds = 50,
                           dist.method = "cosine",
                           mds.method = "cmdscale",
                           rev.mds = FALSE,
                           ndims.print = 1:5,
                           nfeatures.print = 30,
                           reduction.key = "MDS_",
                           seed.use = 42,
                           verbose = TRUE,
                           ...) {
  message("Use ", assay, "/", slot, " to measure '", dist.method, "' distance.")
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
    mds.results2 <- MASS::isoMDS(cell.dist, k = nmds)
  }
  if (mds.method == "sammon") {
    check_R("MASS")
    mds.results3 <- MASS::sammon(cell.dist, k = nmds)
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
RunMDS.Assay <- function(object,
                         assay = NULL,
                         slot = "scale.data",
                         features = NULL,
                         nmds = 50,
                         dist.method = "cosine",
                         mds.method = "cmdscale",
                         rev.mds = FALSE,
                         verbose = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.key = "MDS_",
                         seed.use = 42,
                         ...) {
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
RunMDS.Seurat <- function(object,
                          assay = NULL,
                          slot = "scale.data",
                          features = NULL,
                          nmds = 50,
                          dist.method = "cosine",
                          mds.method = "cmdscale",
                          rev.mds = FALSE,
                          verbose = TRUE,
                          ndims.print = 1:5,
                          nfeatures.print = 30,
                          reduction.name = "mds",
                          reduction.key = "MDS_",
                          seed.use = 42,
                          ...) {
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
#' @importFrom parallelDist parDist
#' @importFrom utils capture.output
#' @importFrom Seurat CreateDimReducObject
#' @importFrom rlang "%||%"
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @export
#'
RunDM.dist <- function(object,
                       exprs,
                       ndcs = 2,
                       sigma = "local",
                       k = NULL,
                       slot = "scale.data",
                       assay = NULL,
                       reduction.key = "DM_",
                       seed.use = 42,
                       verbose = TRUE,
                       ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (is.null(k)) {
    k <- destiny::find_dm_k(destiny:::dataset_n_observations(data = NULL, distances = object) - 1L)
  }
  message("Number of nearest neighbors: ", k)
  dm.results <- destiny::DiffusionMap(
    distance = object, sigma = sigma, k = k,
    n_eigs = ndcs, verbose = verbose, ...
  )
  # gene.relevance <- gene_relevance(coords=dm.results@eigenvectors, exprs=exprs,verbose=verbose)
  cell.embeddings <- dm.results@eigenvectors
  rownames(x = cell.embeddings) <- attr(object, "Labels")
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, dm.results = dm.results)
  )
  return(reduction.data)
}


#' @param object
#'
#' @param ndcs
#' @param sigma
#' @param k
#' @param knn_method
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
#' @importFrom Seurat CreateDimReducObject
#' @export
#' @method RunDM matrix
#'
RunDM.matrix <- function(object,
                         ndcs = 2,
                         sigma = "local",
                         k = NULL,
                         knn_method = "find_knn",
                         dist.method = "euclidean",
                         assay = NULL,
                         slot = "scale.data",
                         reduction.key = "DM_",
                         seed.use = 42,
                         verbose = TRUE,
                         ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ndcs <- min(ndcs, nrow(x = object) - 1)
  x <- as.matrix(object)
  if (knn_method == "parDist") {
    if (dist.method %in% c("pearson", "spearman")) {
      if (dist.method == "spearman") {
        x <- t(apply(x, 1, rank))
      }
      x <- t(apply(x, 1, function(x) x - mean(x)))
      dist.method <- "cosine"
    }
    cell.dist <- parDist(x = x, method = dist.method)
    reduction.data <- RunDM(
      object = cell.dist,
      exprs = x,
      ndcs = ndcs,
      sigma = sigma,
      k = k,
      assay = assay,
      slot = slot,
      reduction.key = reduction.key,
      seed.use = seed.use,
      verbose = verbose,
      ...
    )
  } else if (knn_method == "find_knn") {
    if (is.null(k)) {
      k <- destiny::find_dm_k(destiny:::dataset_n_observations(data = x) - 1L)
    }
    dm.results <- destiny::DiffusionMap(
      data = x, sigma = sigma, k = k, distance = dist.method,
      n_eigs = ndcs, verbose = verbose, ...
    )
    # gene.relevance <- gene_relevance(coords=dm.results@eigenvectors, exprs=exp,verbose=verbose)
    cell.embeddings <- dm.results@eigenvectors
    colnames(x = cell.embeddings) <- paste0(reduction.key, 1:ndcs)
    reduction.data <- CreateDimReducObject(
      embeddings = cell.embeddings,
      assay = assay,
      key = reduction.key,
      misc = list(slot = slot, dm.results = dm.results)
    )
  } else {
    stop("knn_method must be one of 'find_knn' and 'parDist'")
  }

  return(reduction.data)
}

#' @param reduction.name dimensional reduction name, 'dm' by default
#'
#' @param object
#' @param reduction
#' @param dims
#' @param assay
#' @param slot
#' @param features
#' @param ndcs
#' @param sigma
#' @param k
#' @param knn_method
#' @param dist.method
#' @param reduction.key
#' @param seed.use
#' @param verbose
#' @param ...
#'
#' @rdname RunDM
#' @concept dimensional_reduction
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
#' @method RunDM Seurat
#'
RunDM.Seurat <- function(object,
                         reduction = "pca",
                         dims = 1:30,
                         assay = NULL,
                         slot = "scale.data",
                         features = NULL,
                         ndcs = 2,
                         sigma = "local",
                         k = NULL,
                         knn_method = "find_knn",
                         dist.method = "euclidean",
                         reduction.name = "dm",
                         reduction.key = "DM_",
                         seed.use = 42,
                         verbose = TRUE,
                         ...) {
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
    knn_method = knn_method,
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
#' @param nn.name
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
#' @param uwot.sgd
#' @param seed.use
#' @param metric.kwds
#' @param angular.rp.forest
#' @param densmap
#' @param dens.lambda
#' @param dens.frac
#' @param dens.var.shift
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#'
#' @importFrom reticulate py_module_available py_set_seed import
#' @importFrom uwot umap_transform
#' @importFrom future nbrOfWorkers
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#' @importFrom umap umap.defaults
#'
#' @rdname RunUMAP2
#' @concept dimensional_reduction
#' @method RunUMAP2 default
#' @export
#'
RunUMAP2.default <- function(object, dims = NULL, reduction = "pca", features = NULL,
                             graph = NULL, assay = DefaultAssay(object = object), nn.name = NULL,
                             slot = "data", umap.method = "uwot", reduction.model = NULL,
                             return.model = FALSE, n.neighbors = 30L, n.components = 2L,
                             metric = "cosine", n.epochs = NULL, learning.rate = 1, min.dist = 0.3,
                             spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L,
                             repulsion.strength = 1, negative.sample.rate = 5L, a = NULL,
                             b = NULL, uwot.sgd = FALSE, seed.use = 42L, metric.kwds = NULL,
                             angular.rp.forest = FALSE, densmap = FALSE, dens.lambda = 2,
                             dens.frac = 0.3, dens.var.shift = 0.1, verbose = TRUE, reduction.name = "umap",
                             reduction.key = "UMAP_", ...) {
  Seurat:::CheckDots(...)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (umap.method != "umap-learn" && getOption(
    "Seurat.warn.umap.uwot",
    TRUE
  )) {
    warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric",
      "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'",
      "\nThis message will be shown once per session",
      call. = FALSE, immediate. = TRUE
    )
    options(Seurat.warn.umap.uwot = FALSE)
  }
  if (umap.method == "uwot-learn") {
    warning("'uwot-learn' is deprecated. Set umap.method = 'uwot' and return.model = TRUE")
    umap.method <- "uwot"
    return.model <- TRUE
  }
  if (densmap && umap.method != "umap-learn") {
    warning("densmap is only supported by umap-learn method. Method is changed to 'umap-learn'")
    umap.method <- "umap-learn"
  }
  if (return.model) {
    if (verbose) {
      message("UMAP will return its model")
    }
    if (!umap.method %in% c("uwot", "naive")) {
      stop("'umap.method' must be one of 'uwot' or 'naive' when 'return.model' is TRUE")
    }
  }
  if (inherits(x = object, what = "Neighbor")) {
    object <- list(idx = Indices(object), dist = Distances(object))
  }
  if (!is.null(x = reduction.model)) {
    if (verbose) {
      message("Running UMAP projection")
    }
    umap.method <- "uwot-predict"
  }
  umap.output <- switch(EXPR = umap.method,
    `umap-learn` = {
      if (!py_module_available(module = "umap")) {
        stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
      }
      if (!py_module_available(module = "sklearn")) {
        stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
      }
      if (!is.null(x = seed.use)) {
        py_set_seed(seed = seed.use)
      }
      if (typeof(x = n.epochs) == "double") {
        n.epochs <- as.integer(x = n.epochs)
      }
      umap_import <- import(module = "umap", delay_load = TRUE)
      sklearn <- import("sklearn", delay_load = TRUE)
      if (densmap && numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) <
        numeric_version(x = "0.5.0")) {
        stop("densmap is only supported by versions >= 0.5.0 of umap-learn. Upgrade umap-learn (e.g. pip install --upgrade umap-learn).")
      }
      random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
      umap.args <- list(
        n_neighbors = as.integer(x = n.neighbors),
        n_components = as.integer(x = n.components), metric = metric,
        n_epochs = n.epochs, learning_rate = learning.rate,
        min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
        local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
        negative_sample_rate = negative.sample.rate, random_state = random.state,
        a = a, b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest,
        verbose = verbose
      )
      if (numeric_version(x = umap_import$pkg_resources$get_distribution("umap-learn")$version) >=
        numeric_version(x = "0.5.0")) {
        umap.args <- c(umap.args, list(
          densmap = densmap,
          dens_lambda = dens.lambda, dens_frac = dens.frac,
          dens_var_shift = dens.var.shift, output_dens = FALSE
        ))
      }
      umap <- do.call(what = umap_import$UMAP, args = umap.args)
      umap$fit_transform(as.matrix(x = object))
    },
    uwot = {
      if (is.list(x = object)) {
        uwot::umap(
          X = NULL, nn_method = object, n_threads = nbrOfWorkers(),
          n_components = as.integer(x = n.components),
          metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
          min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose,
          ret_model = return.model
        )
      } else {
        uwot::umap(
          X = object, n_threads = nbrOfWorkers(), n_neighbors = as.integer(x = n.neighbors),
          n_components = as.integer(x = n.components),
          metric = metric, n_epochs = n.epochs, learning_rate = learning.rate,
          min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio,
          local_connectivity = local.connectivity, repulsion_strength = repulsion.strength,
          negative_sample_rate = negative.sample.rate,
          a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose,
          ret_model = return.model
        )
      }
    },
    `uwot-predict` = {
      if (metric == "correlation") {
        warning("UWOT does not implement the correlation metric, using cosine instead",
          call. = FALSE, immediate. = TRUE
        )
        metric <- "cosine"
      }
      if (is.null(x = reduction.model) || !inherits(
        x = reduction.model,
        what = "DimReduc"
      )) {
        stop("If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.",
          call. = FALSE
        )
      }
      model <- Misc(object = reduction.model, slot = "model")
      if (length(x = model) == 0) {
        stop("The provided reduction.model does not have a model stored. Please try running umot-learn on the object first",
          call. = FALSE
        )
      }
      if (is.list(x = object)) {
        if (ncol(object$idx) != model$n_neighbors) {
          warning(
            "Number of neighbors between query and reference ",
            "is not equal to the number of neighbros within reference"
          )
          model$n_neighbors <- ncol(object$idx)
        }
        umap_transform(
          X = NULL, nn_method = object, model = model,
          n_threads = nbrOfWorkers(), n_epochs = n.epochs,
          verbose = verbose
        )
      } else {
        umap_transform(
          X = object, model = model, n_threads = nbrOfWorkers(),
          n_epochs = n.epochs, verbose = verbose
        )
      }
    },
    naive = {
      umap.config <- umap.defaults
      umap.config$n_neighbors <- min(n.neighbors, ncol(object))
      umap.config$n_components <- min(n.components, ncol(object))
      umap.config$metric <- metric
      umap.config$n_epochs <- ifelse(is.null(n.epochs), 200, n.epochs)
      umap.config$min_dist <- min.dist
      umap.config$set_op_mix_ratio <- set.op.mix.ratio
      umap.config$local_connectivity <- local.connectivity
      umap.config$a <- ifelse(is.null(a), NA, a)
      umap.config$b <- ifelse(is.null(b), NA, b)
      umap.config$random_state <- seed.use
      umap.config$transform_state <- seed.use
      umap.config$verbose <- verbose
      umap::umap(object, config = umap.config)
    },
    stop("Unknown umap method: ", umap.method, call. = FALSE)
  )
  if (isTRUE(return.model)) {
    if (umap.method != "naive") {
      umap.output$nn_index <- NULL
      umap.model <- umap.output
      umap.output <- umap.output$embedding
    } else {
      umap.model <- umap.output
      umap.output <- umap.output$layout
    }
  }

  colnames(x = umap.output) <- paste0(reduction.key, 1:ncol(x = umap.output))
  if (inherits(x = object, what = "dist") && !is.null(attr(x = object, "Labels"))) {
    rownames(x = umap.output) <- attr(x = object, "Labels")
  } else if (is.list(x = object)) {
    rownames(x = umap.output) <- rownames(x = object$idx)
  } else {
    rownames(x = umap.output) <- rownames(x = object)
  }
  umap.reduction <- CreateDimReducObject(
    embeddings = umap.output,
    key = reduction.key, assay = assay, global = TRUE
  )
  if (return.model) {
    Misc(umap.reduction, slot = "model") <- umap.model
  }
  return(umap.reduction)
}

#' @param object
#'
#' @param assay
#' @param umap.method
#' @param n.components
#' @param metric
#' @param n.epochs
#' @param learning.rate
#' @param min.dist
#' @param spread
#' @param repulsion.strength
#' @param negative.sample.rate
#' @param a
#' @param b
#' @param uwot.sgd
#' @param seed.use
#' @param metric.kwds
#' @param densmap
#' @param densmap.kwds
#' @param verbose
#' @param reduction.key
#' @param ...
#'
#' @importFrom reticulate py_module_available import
#' @importFrom Seurat CreateDimReducObject
#'
#' @rdname RunUMAP2
#' @concept dimensional_reduction
#' @method RunUMAP2 Graph
#' @export
#'
RunUMAP2.Graph <- function(object, assay = NULL,
                           umap.method = "uwot",
                           n.components = 2L,
                           metric = "correlation",
                           n.epochs = 0L,
                           learning.rate = 1,
                           min.dist = 0.3,
                           spread = 1,
                           repulsion.strength = 1,
                           negative.sample.rate = 5L,
                           a = NULL,
                           b = NULL,
                           uwot.sgd = FALSE,
                           seed.use = 42L,
                           metric.kwds = NULL,
                           densmap = FALSE,
                           densmap.kwds = NULL,
                           verbose = TRUE,
                           reduction.key = "UMAP_",
                           ...) {
  Seurat:::CheckDots(...)
  if (umap.method != "umap-learn") {
    warning(
      "Running UMAP on Graph objects is only supported using the umap-learn method",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!py_module_available(module = "umap")) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!py_module_available(module = "numpy")) {
    stop("Cannot find numpy, please install through pip (e.g. pip install numpy).")
  }
  if (!py_module_available(module = "sklearn")) {
    stop("Cannot find sklearn, please install through pip (e.g. pip install scikit-learn).")
  }
  if (!py_module_available(module = "scipy")) {
    stop("Cannot find scipy, please install through pip (e.g. pip install scipy).")
  }
  np <- import("numpy", delay_load = TRUE)
  sp <- import("scipy", delay_load = TRUE)
  sklearn <- import("sklearn", delay_load = TRUE)
  umap <- import("umap", delay_load = TRUE)
  diag(x = object) <- 0
  data <- object
  object <- sp$sparse$coo_matrix(arg1 = object)
  ab.params <- umap$umap_$find_ab_params(spread = spread, min_dist = min.dist)
  a <- a %||% ab.params[[1]]
  b <- b %||% ab.params[[2]]
  n.epochs <- n.epochs %||% 0L
  random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
  umap.args <- list(
    data = data,
    graph = object,
    n_components = n.components,
    initial_alpha = learning.rate,
    a = a,
    b = b,
    gamma = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    n_epochs = as.integer(x = n.epochs),
    random_state = random.state,
    init = "spectral",
    metric = metric,
    metric_kwds = metric.kwds,
    verbose = verbose
  )
  if (numeric_version(x = umap$pkg_resources$get_distribution("umap-learn")$version) >=
    numeric_version(x = "0.5.0")) {
    umap.args <- c(umap.args, list(
      densmap = densmap,
      densmap_kwds = densmap.kwds,
      output_dens = FALSE
    ))
  }
  embeddings <- do.call(what = umap$umap_$simplicial_set_embedding, args = umap.args)
  if (length(x = embeddings) == 2) {
    embeddings <- embeddings[[1]]
  }
  rownames(x = embeddings) <- colnames(x = data)
  colnames(x = embeddings) <- paste0("UMAP_", 1:n.components)
  # center the embeddings on zero
  embeddings <- scale(x = embeddings, scale = FALSE)
  umap <- CreateDimReducObject(
    embeddings = embeddings,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(umap)
}

#' @param object
#'
#' @param reduction.model
#' @param ...
#'
#' @rdname RunUMAP2
#' @concept dimensional_reduction
#' @method RunUMAP2 Neighbor
#' @export
#'
RunUMAP2.Neighbor <- function(object,
                              reduction.model,
                              ...) {
  neighborlist <- list(
    "idx" = Indices(object),
    "dist" = Distances(object)
  )
  RunUMAP2(
    object = neighborlist,
    reduction.model = reduction.model,
    ...
  )
}

#' @param object
#'
#' @param dims
#' @param reduction
#' @param features
#' @param graph
#' @param assay
#' @param nn.name
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
#' @param uwot.sgd
#' @param seed.use
#' @param metric.kwds
#' @param angular.rp.forest
#' @param densmap
#' @param dens.lambda
#' @param dens.frac
#' @param dens.var.shift
#' @param verbose
#' @param reduction.name
#' @param reduction.key
#' @param ...
#'
#' @rdname RunUMAP2
#' @concept dimensional_reduction
#' @export
#' @importFrom Seurat LogSeuratCommand
#' @method RunUMAP2 Seurat
RunUMAP2.Seurat <- function(object, dims = NULL, reduction = "pca", features = NULL, graph = NULL,
                            assay = DefaultAssay(object = object), nn.name = NULL, slot = "data",
                            umap.method = "uwot", reduction.model = NULL, return.model = FALSE,
                            n.neighbors = 30L, n.components = 2L, metric = "cosine", n.epochs = NULL,
                            learning.rate = 1, min.dist = 0.3, spread = 1, set.op.mix.ratio = 1,
                            local.connectivity = 1L, repulsion.strength = 1, negative.sample.rate = 5L,
                            a = NULL, b = NULL, uwot.sgd = FALSE, seed.use = 42L, metric.kwds = NULL,
                            angular.rp.forest = FALSE, densmap = FALSE, dens.lambda = 2, dens.frac = 0.3,
                            dens.var.shift = 0.1, verbose = TRUE,
                            reduction.name = "umap", reduction.key = "UMAP_",
                            ...) {
  Seurat:::CheckDots(...)
  if (sum(c(is.null(x = dims), is.null(x = features), is.null(x = graph))) < 2) {
    stop("Please specify only one of the following arguments: dims, features, or graph")
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
    data.use <- Embeddings(object[[reduction]])[, dims]
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
  } else if (!is.null(x = nn.name)) {
    if (!inherits(x = object[[nn.name]], what = "Neighbor")) {
      stop(
        "Please specify a Neighbor object name, ",
        "instead of the name of a ",
        class(object[[nn.name]]),
        " object",
        call. = FALSE
      )
    }
    data.use <- object[[nn.name]]
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
    stop("Please specify one of dims, features, or graph")
  }
  object[[reduction.name]] <- RunUMAP2(
    object = data.use,
    reduction.model = reduction.model,
    return.model = return.model,
    assay = assay,
    umap.method = umap.method,
    n.neighbors = n.neighbors,
    n.components = n.components,
    metric = metric,
    n.epochs = n.epochs,
    learning.rate = learning.rate,
    min.dist = min.dist,
    spread = spread,
    set.op.mix.ratio = set.op.mix.ratio,
    local.connectivity = local.connectivity,
    repulsion.strength = repulsion.strength,
    negative.sample.rate = negative.sample.rate,
    a = a,
    b = b,
    uwot.sgd = uwot.sgd,
    seed.use = seed.use,
    metric.kwds = metric.kwds,
    angular.rp.forest = angular.rp.forest,
    densmap = densmap,
    dens.lambda = dens.lambda,
    dens.frac = dens.frac,
    dens.var.shift = dens.var.shift,
    reduction.key = reduction.key,
    verbose = verbose
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' Harmony single cell integration
#'
#' Run Harmony algorithm with Seurat pipelines.
#'
#' @param object An Seurat object
#' @inheritParams harmony::RunHarmony
#'
#'
#' @export
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
#' @importFrom Seurat Embeddings RunPCA FetchData CreateDimReducObject ProjectDim LogSeuratCommand
#' @importFrom harmony HarmonyMatrix
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
  if (reduction == "pca") {
    tryCatch(
      embedding <- Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) {
          message("Harmony needs PCA. Trying to run PCA now.")
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object,
            assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Embeddings(object, reduction = reduction)
  }
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
  metavars_df <- FetchData(object, group.by.vars)

  harmonyObject <- HarmonyMatrix(
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
      misc = list(R = harmonyClusters)
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
