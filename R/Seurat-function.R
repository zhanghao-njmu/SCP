#' Run NMF (non-negative matrix factorization)
#'
#' @param object An object. This can be a Seurat object, an Assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis. Default is NULL.
#' @param slot A character string specifying the slot name to be used for the analysis. Default is "data".
#' @param features A character vector specifying the features to be used for the analysis. Default is NULL, which uses all variable features.
#' @param nbes An integer specifying the number of basis vectors (components) to be computed. Default is 50.
#' @param nmf.method A character string specifying the NMF algorithm to be used. Currently supported values are "RcppML" and "NMF". Default is "RcppML".
#' @param tol A numeric value specifying the tolerance for convergence (only applicable when nmf.method is "RcppML"). Default is 1e-5.
#' @param maxit An integer specifying the maximum number of iterations for convergence (only applicable when nmf.method is "RcppML"). Default is 100.
#' @param rev.nmf A logical value indicating whether to perform reverse NMF (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param ndims.print An integer vector specifying the dimensions (number of basis vectors) to print in the output. Default is 1:5.
#' @param nfeatures.print An integer specifying the number of features to print in the output. Default is 30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "nmf".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "BE_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments passed to RcppML::\link[RcppML]{nmf} or NMF::\link[NMF]{nmf} function.
#'
#' @examples
#' pancreas_sub <- RunNMF(object = pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "nmf")
#'
#' @rdname RunNMF
#' @export
#'
RunNMF <- function(object, ...) {
  UseMethod(generic = "RunNMF", object = object)
}

#' @rdname RunNMF
#' @method RunNMF Seurat
#' @importFrom Seurat LogSeuratCommand VariableFeatures DefaultAssay GetAssay
#' @export
RunNMF.Seurat <- function(object, assay = NULL, slot = "data", features = NULL, nbes = 50,
                          nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                          ndims.print = 1:5, nfeatures.print = 30,
                          reduction.name = "nmf", reduction.key = "BE_",
                          verbose = TRUE, seed.use = 11, ...) {
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

#' @rdname RunNMF
#' @method RunNMF Assay
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
RunNMF.Assay <- function(object, assay = NULL, slot = "data", features = NULL, nbes = 50,
                         nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                         ndims.print = 1:5, nfeatures.print = 30,
                         reduction.key = "BE_", verbose = TRUE, seed.use = 11,
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

#' @rdname RunNMF
#' @method RunNMF default
#' @importFrom utils capture.output
#' @importFrom Matrix t
#' @importFrom Seurat CreateDimReducObject
#' @export
RunNMF.default <- function(object, assay = NULL, slot = "data", nbes = 50,
                           nmf.method = "RcppML", tol = 1e-5, maxit = 100, rev.nmf = FALSE,
                           ndims.print = 1:5, nfeatures.print = 30,
                           reduction.key = "BE_", verbose = TRUE, seed.use = 11, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.nmf) {
    object <- t(object)
  }
  nbes <- min(nbes, nrow(x = object) - 1)
  if (nmf.method == "RcppML") {
    check_R("zdebruine/RcppML")
    options("RcppML.verbose" = FALSE)
    options("RcppML.threads" = 0)
    if (!"package:Matrix" %in% search()) {
      attachNamespace("Matrix")
    }
    nmf.results <- RcppML::nmf(
      t(object),
      k = nbes, tol = tol, maxit = maxit,
      verbose = verbose, ...
    )
    cell.embeddings <- nmf.results$w
    feature.loadings <- t(nmf.results$h)
  }
  if (nmf.method == "NMF") {
    check_R("NMF")
    seed <- NMF::seed
    nmf.results <- NMF::nmf(x = as_matrix(t(object)), rank = nbes)
    cell.embeddings <- nmf.results@fit@W
    feature.loadings <- t(nmf.results@fit@H)
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
    misc = list(slot = slot, nmf.results = nmf.results)
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

#' Run MDS (multi-dimensional scaling)
#'
#' @param object An object. This can be a Seurat object, an assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis. Default is NULL.
#' @param slot A character string specifying the slot name to be used for the analysis. Default is "data".
#' @param features A character vector specifying the features to be used for the analysis. Default is NULL, which uses all variable features.
#' @param nmds An integer specifying the number of dimensions to be computed. Default is 50.
#' @param dist.method A character string specifying the distance metric to be used. Currently supported values are "euclidean", "chisquared","kullback", "jeffreys", "jensen", "manhattan", "maximum", "canberra", "minkowski", and "hamming". Default is "euclidean".
#' @param mds.method A character string specifying the MDS algorithm to be used. Currently supported values are "cmdscale", "isoMDS", and "sammon". Default is "cmdscale".
#' @param rev.mds A logical value indicating whether to perform reverse MDS (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "mds".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "MDS_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the stats::\link[stats]{cmdscale}, MASS::\link[MASS]{isoMDS} or MASS::\link[MASS]{sammon} function.
#'
#' @examples
#' pancreas_sub <- RunMDS(object = pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "mds")
#'
#' @rdname RunMDS
#' @export
#'
RunMDS <- function(object, ...) {
  UseMethod(generic = "RunMDS", object = object)
}

#' @rdname RunMDS
#' @method RunMDS Seurat
#' @importFrom Seurat LogSeuratCommand VariableFeatures DefaultAssay GetAssay
#' @export
RunMDS.Seurat <- function(object, assay = NULL, slot = "data",
                          features = NULL, nmds = 50, dist.method = "euclidean", mds.method = "cmdscale",
                          rev.mds = FALSE,
                          reduction.name = "mds", reduction.key = "MDS_",
                          verbose = TRUE, seed.use = 11, ...) {
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
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunMDS
#' @method RunMDS Assay
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
RunMDS.Assay <- function(object, assay = NULL, slot = "data",
                         features = NULL, nmds = 50, dist.method = "euclidean", mds.method = "cmdscale",
                         rev.mds = FALSE,
                         reduction.key = "MDS_", verbose = TRUE, seed.use = 11, ...) {
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
    reduction.key = reduction.key,
    seed.use = seed.use,
    ...
  )
  return(reduction.data)
}

#' @rdname RunMDS
#' @method RunMDS default
#' @importFrom proxyC dist
#' @importFrom utils capture.output
#' @importFrom Matrix t
#' @importFrom stats cmdscale as.dist
#' @importFrom Seurat CreateDimReducObject
#' @export
RunMDS.default <- function(object, assay = NULL, slot = "data",
                           nmds = 50, dist.method = "euclidean", mds.method = "cmdscale",
                           rev.mds = FALSE,
                           reduction.key = "MDS_", verbose = TRUE, seed.use = 11, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (rev.mds) {
    object <- t(object)
  }
  nmds <- min(nmds, nrow(x = object) - 1)
  x <- t(as_matrix(object))
  cell.dist <- as.dist(dist(x = x, method = dist.method))
  if (mds.method == "cmdscale") {
    mds.results <- cmdscale(cell.dist, k = nmds, eig = TRUE, ...)
  }
  if (mds.method == "isoMDS") {
    check_R("MASS")
    mds.results <- MASS::isoMDS(cell.dist, k = nmds, ...)
  }
  if (mds.method == "sammon") {
    check_R("MASS")
    mds.results <- MASS::sammon(cell.dist, k = nmds, ...)
  }
  cell.embeddings <- mds.results$points

  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:nmds)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    assay = assay,
    key = reduction.key,
    misc = list(slot = slot, mds.results = mds.results)
  )
  return(reduction.data)
}

#' Run GLMPCA (generalized version of principal components analysis)
#'
#' @param object An object. This can be a Seurat object, an assay object, or a matrix-like object.
#' @param assay A character string specifying the assay to be used for the analysis. Default is NULL.
#' @param slot A character string specifying the slot name to be used for the analysis. Default is "counts".
#' @param features A character vector specifying the features to be used for the analysis. Default is NULL, which uses all variable features.
#' @param L An integer specifying the number of components to be computed. Default is 5.
#' @param fam A character string specifying the family of the generalized linear model to be used. Currently supported values are "poi", "nb", "nb2", "binom", "mult", and "bern". Default is "poi".
#' @param rev.gmlpca A logical value indicating whether to perform reverse GLMPCA (i.e., transpose the input matrix) before running the analysis. Default is FALSE.
#' @param ndims.print An integer vector specifying the dimensions (number of components) to print in the output. Default is 1:5.
#' @param nfeatures.print An integer specifying the number of features to print in the output. Default is 30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "glmpca".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "GLMPC_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the \link[glmpca]{glmpca} function.
#'
#' @examples
#' pancreas_sub <- RunGLMPCA(object = pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "glmpca")
#'
#' @rdname RunGLMPCA
#' @export
RunGLMPCA <- function(object, ...) {
  UseMethod(generic = "RunGLMPCA", object = object)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA Seurat
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
RunGLMPCA.Seurat <- function(object, assay = NULL, slot = "counts",
                             features = NULL, L = 5, fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                             rev.gmlpca = FALSE, ndims.print = 1:5, nfeatures.print = 30,
                             reduction.name = "glmpca", reduction.key = "GLMPC_",
                             verbose = TRUE, seed.use = 11, ...) {
  features <- features %||% VariableFeatures(object = object)
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunGLMPCA(
    object = assay.data,
    assay = assay,
    slot = slot,
    L = L,
    fam = fam,
    rev.gmlpca = rev.gmlpca,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA Assay
#' @importFrom stats var
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
RunGLMPCA.Assay <- function(object, assay = NULL, slot = "counts",
                            features = NULL, L = 5, fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                            rev.gmlpca = FALSE, ndims.print = 1:5, nfeatures.print = 30,
                            reduction.key = "GLMPC_", verbose = TRUE, seed.use = 11, ...) {
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
    rev.gmlpca = rev.gmlpca,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @rdname RunGLMPCA
#' @method RunGLMPCA default
#' @importFrom Seurat DefaultAssay DefaultAssay<- CreateDimReducObject Tool<- LogSeuratCommand
#' @export
RunGLMPCA.default <- function(object, assay = NULL, slot = "counts",
                              features = NULL, L = 5, fam = c("poi", "nb", "nb2", "binom", "mult", "bern"),
                              rev.gmlpca = FALSE, ndims.print = 1:5, nfeatures.print = 30,
                              reduction.key = "GLMPC_", verbose = TRUE, seed.use = 11, ...) {
  check_R("glmpca")
  if (inherits(object, "dgCMatrix")) {
    object <- as_matrix(object)
  }
  fam <- match.arg(fam)
  glmpca_results <- glmpca::glmpca(Y = object, L = L, fam = fam, ...)
  glmpca_dimnames <- paste0(reduction.key, seq_len(L))
  factors <- as_matrix(glmpca_results$factors)
  loadings <- as_matrix(glmpca_results$loadings)
  colnames(x = factors) <- glmpca_dimnames
  colnames(x = loadings) <- glmpca_dimnames
  factors_l2_norm <- sqrt(colSums(factors^2))
  class(glmpca_results) <- NULL
  glmpca_results$factors <- glmpca_results$loadings <- NULL
  reduction.data <- CreateDimReducObject(
    embeddings = factors,
    key = reduction.key,
    loadings = loadings,
    stdev = factors_l2_norm,
    assay = assay,
    global = TRUE,
    misc = list(slot = slot, glmpca.results = glmpca_results)
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

#' Run DM (diffusion map)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used.
#' @param dims An integer vector specifying the dimensions to be used. Default is 1:30.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param ndcs An integer specifying the number of diffusion components (dimensions) to be computed. Default is 2.
#' @param sigma A character string specifying the diffusion scale parameter of the Gaussian kernel. Currently supported values are "local" (default) and "global".
#' @param k An integer specifying the number of nearest neighbors to be used for the construction of the graph. Default is 30.
#' @param dist.method A character string specifying the distance metric to be used for the construction of the knn graph. Currently supported values are "euclidean" and "cosine". Default is "euclidean".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "dm".
#' @param reduction.key A character string specifying the prefix for the column names of the basis vectors. Default is "DM_".
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param ... Additional arguments to be passed to the \link[destiny]{DiffusionMap} function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunDM(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "dm")
#'
#' @rdname RunDM
#' @export
#'
RunDM <- function(object, ...) {
  UseMethod(generic = "RunDM", object = object)
}

#' @rdname RunDM
#' @method RunDM Seurat
#' @importFrom Seurat LogSeuratCommand DefaultAssay GetAssayData Embeddings
#' @export
RunDM.Seurat <- function(object,
                         reduction = "pca", dims = 1:30,
                         features = NULL, assay = NULL, slot = "data",
                         ndcs = 2, sigma = "local", k = 30, dist.method = "euclidean",
                         reduction.name = "dm", reduction.key = "DM_",
                         verbose = TRUE, seed.use = 11, ...) {
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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

#' @rdname RunDM
#' @method RunDM default
#' @importFrom utils capture.output
#' @importFrom Seurat CreateDimReducObject
#' @export
RunDM.default <- function(object, assay = NULL, slot = "data",
                          ndcs = 2, sigma = "local", k = 30, dist.method = "euclidean",
                          reduction.key = "DM_", verbose = TRUE, seed.use = 11, ...) {
  check_R("destiny")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  dm.results <- destiny::DiffusionMap(data = as_matrix(object), n_eigs = ndcs, sigma = sigma, k = k, distance = dist.method, verbose = verbose, ...)

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
}

#' Run UMAP (Uniform Manifold Approximation and Projection)
#'
#' @param object An object. This can be a Seurat object, a matrix-like object, a Neighbor object, or a Graph object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param neighbor A character string specifying the name of the Neighbor object to be used. Default is NULL.
#' @param graph A character string specifying the name of the Graph object to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param umap.method A character string specifying the UMAP method to be used. Options are "naive" and uwot". Default is "uwot".
#' @param reduction.model A DimReduc object containing a pre-trained UMAP model. Default is NULL.
#' @param return.model A logical value indicating whether to return the UMAP model. Default is FALSE.
#' @param n.neighbors An integer specifying the number of nearest neighbors to be used. Default is 30.
#' @param n.components An integer specifying the number of UMAP components. Default is 2.
#' @param metric A character string specifying the metric or a function to be used for distance calculations. When using a string, available metrics are: euclidean, manhattan. Other available generalized metrics are: cosine, pearson, pearson2. Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal. When using metric.function as a function, the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns.  Default is "cosine".
#' @param n.epochs An integer specifying the number of iterations performed during layout optimization for UMAP. Default is 200.
#' @param spread A numeric value specifying the spread parameter for UMAP, used during automatic estimation of a/b parameters. Default is 1.
#' @param min.dist A numeric value specifying the minimum distance between UMAP embeddings, determines how close points appear in the final layout. Default is 0.3.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets. Both fuzzy set operations use the product t-norm. The value of this parameter should be between 0.0 and 1.0; a value of 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy intersection.
#' @param local.connectivity An integer specifying the local connectivity, used during construction of fuzzy simplicial set. Default is 1.
#' @param negative.sample.rate An integer specifying the negative sample rate for UMAP optimization. Determines how many non-neighbor points are used per point and per iteration during layout optimization. Default is 5.
#' @param a A numeric value specifying the parameter a for UMAP optimization. Contributes to gradient calculations during layout optimization. When left at NA, a suitable value will be estimated automatically. Default is NULL.
#' @param b A numeric value specifying the parameter b for UMAP optimization. Contributes to gradient calculations during layout optimization. When left at NA, a suitable value will be estimated automatically. Default is NULL.
#' @param learning.rate A numeric value specifying the initial value of "learning rate" of layout optimization. Default is 1.
#' @param repulsion.strength A numeric value determines, together with alpha, the learning rate of layout optimization. Default is 1.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "umap".
#' @param reduction.key A character string specifying the prefix for the column names of the UMAP embeddings. Default is "UMAP_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Unused argument.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunUMAP2(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "umap")
#'
#' @rdname RunUMAP2
#' @export
RunUMAP2 <- function(object, ...) {
  UseMethod(generic = "RunUMAP2", object = object)
}

#' @rdname RunUMAP2
#' @method RunUMAP2 Seurat
#' @importFrom Seurat LogSeuratCommand
#' @export
RunUMAP2.Seurat <- function(object,
                            reduction = "pca", dims = NULL, features = NULL, neighbor = NULL, graph = NULL,
                            assay = NULL, slot = "data",
                            umap.method = "uwot", reduction.model = NULL, n_threads = NULL,
                            return.model = FALSE, n.neighbors = 30L, n.components = 2L,
                            metric = "cosine", n.epochs = 200L, spread = 1, min.dist = 0.3,
                            set.op.mix.ratio = 1, local.connectivity = 1L, negative.sample.rate = 5L,
                            a = NULL, b = NULL, learning.rate = 1, repulsion.strength = 1,
                            reduction.name = "umap", reduction.key = "UMAP_",
                            verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features), is.null(neighbor), is.null(x = graph))) == 4) {
    stop("Please specify only one of the following arguments: dims, features, neighbor or graph")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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
    data.use <- as_matrix(Embeddings(object[[reduction]])[, dims])
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

#' @rdname RunUMAP2
#' @method RunUMAP2 default
#' @importFrom SeuratObject Indices Distances as.sparse CreateDimReducObject Misc<- Misc
#' @importFrom Matrix sparseMatrix
#' @export
RunUMAP2.default <- function(object, assay = NULL,
                             umap.method = "uwot", reduction.model = NULL, n_threads = NULL,
                             return.model = FALSE, n.neighbors = 30L, n.components = 2L,
                             metric = "cosine", n.epochs = 200L, spread = 1, min.dist = 0.3,
                             set.op.mix.ratio = 1, local.connectivity = 1L, negative.sample.rate = 5L,
                             a = NULL, b = NULL, learning.rate = 1, repulsion.strength = 1,
                             reduction.key = "UMAP_", verbose = TRUE, seed.use = 11L, ...) {
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
      if (!isSymmetric(as_matrix(object[obs_sample, obs_sample]))) {
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
      #   object <- as_matrix(object)
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
        X = object, n_neighbors = n.neighbors, n_threads = n_threads, n_components = n.components,
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
        X = NULL, nn_method = object, n_threads = n_threads, n_components = n.components,
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
      if (!isSymmetric(as_matrix(object[obs_sample, obs_sample]))) {
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
      idx[is.na(idx)] <- sample(1:nrow(object), size = sum(is.na(idx)), replace = TRUE)
      nn <- list(idx = idx, dist = max(connectivity) - connectivity + min(diff(range(connectivity)), 1) / 1e50)
      # idx <- t(as_matrix(apply(object, 2, function(x) order(x, decreasing = TRUE)[1:n.neighbors])))
      # connectivity <- t(as_matrix(apply(object, 2, function(x) x[order(x, decreasing = TRUE)[1:n.neighbors]])))
      out <- uwot::umap(
        X = NULL, nn_method = nn, n_threads = n_threads, n_components = n.components,
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
        X = object, n_neighbors = n.neighbors, n_threads = n_threads, n_components = n.components,
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
        n_threads = n_threads, verbose = verbose
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
      match_k <- t(as_matrix(apply(object, 2, function(x) order(x, decreasing = TRUE)[1:n.neighbors])))
      match_k_connectivity <- t(as_matrix(apply(object, 2, function(x) x[order(x, decreasing = TRUE)[1:n.neighbors]])))
      object <- list(idx = match_k, dist = max(match_k_connectivity) - match_k_connectivity)
      if (is.null(model$num_precomputed_nns) || model$num_precomputed_nns == 0) {
        model$num_precomputed_nns <- 1
      }
      embeddings <- uwot::umap_transform(
        X = NULL, nn_method = object, model = model, n_epochs = n.epochs,
        n_threads = n_threads, verbose = verbose
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
        n_threads = n_threads, verbose = verbose
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

#' Run PaCMAP (Pairwise Controlled Manifold Approximation)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param n_components An integer specifying the number of PaCMAP components. Default is 2.
#' @param n.neighbors An integer specifying the number of neighbors considered in the k-Nearest Neighbor graph. Default to 10 for dataset whose sample size is smaller than 10000. For large dataset whose sample size (n) is larger than 10000, the default value is: 10 + 15 * (log10(n) - 4).
#' @param MN_ratio A numeric value specifying the ratio of the ratio of the number of mid-near pairs to the number of neighbors. Default is 0.5.
#' @param FP_ratio A numeric value specifying the ratio of the ratio of the number of further pairs to the number of neighbors. Default is 2.
#' @param distance_method A character string specifying the distance metric to be used. Default is "euclidean".
#' @param lr A numeric value specifying the learning rate of the AdaGrad optimizer. Default is 1.
#' @param num_iters An integer specifying the number of iterations for PaCMAP optimization. Default is 450.
#' @param apply_pca A logical value indicating whether pacmap should apply PCA to the data before constructing the k-Nearest Neighbor graph. Using PCA to preprocess the data can largely accelerate the DR process without losing too much accuracy. Notice that this option does not affect the initialization of the optimization process. Default is TRUE.
#' @param init A character string specifying the initialization of the lower dimensional embedding. One of "pca" or "random". Default is "random".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "pacmap".
#' @param reduction.key A character string specifying the prefix for the column names of the PaCMAP embeddings. Default is "PaCMAP_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the pacmap.PaCMAP function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunPaCMAP(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "pacmap")
#'
#' @rdname RunPaCMAP
#' @export
RunPaCMAP <- function(object, ...) {
  UseMethod(generic = "RunPaCMAP", object = object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP Seurat
#' @importFrom Seurat LogSeuratCommand
#' @export
RunPaCMAP.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                             assay = NULL, slot = "data",
                             n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                             distance_method = "euclidean",
                             lr = 1, num_iters = 450L, apply_pca = TRUE, init = "random",
                             reduction.name = "pacmap", reduction.key = "PaCMAP_",
                             verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) < 1) {
    stop("Please specify only one of the following arguments: dims, features, or graph")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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
    distance_method = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca, init = init,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunPaCMAP
#' @method RunPaCMAP default
#' @importFrom Seurat CreateDimReducObject
#' @importFrom reticulate import
#' @export
RunPaCMAP.default <- function(object, assay = NULL,
                              n_components = 2, n.neighbors = NULL, MN_ratio = 0.5, FP_ratio = 2,
                              distance_method = "euclidean",
                              lr = 1, num_iters = 450L, apply_pca = TRUE, init = "random",
                              reduction.key = "PaCMAP_", verbose = TRUE, seed.use = 11L, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_Python("pacmap")
  pacmap <- import("pacmap")

  operator <- pacmap$PaCMAP(
    n_components = as.integer(n_components), n_neighbors = n.neighbors, MN_ratio = MN_ratio, FP_ratio = FP_ratio,
    distance = distance_method,
    lr = lr, num_iters = num_iters, apply_pca = apply_pca,
    verbose = verbose, random_state = as.integer(seed.use)
  )
  embedding <- operator$fit_transform(object, init = init)

  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  return(reduction)
}

#' Run PHATE (Potential of Heat-diffusion for Affinity-based Trajectory Embedding)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param n_components An integer specifying the number of PHATE components. Default is 2.
#' @param knn An integer specifying the number of nearest neighbors on which to build kernel. Default is 5.
#' @param decay An integer specifying the sets decay rate of kernel tails. Default is 40.
#' @param n_landmark An integer specifying the number of landmarks to use in fast PHATE. Default is 2000.
#' @param t A character string specifying the power to which the diffusion operator is powered. This sets the level of diffusion. If ‘auto’, t is selected according to the knee point in the Von Neumann Entropy of the diffusion operator. Default is "auto".
#' @param gamma A numeric value specifying the informational distance constant between -1 and 1. gamma=1 gives the PHATE log potential, gamma=0 gives a square root potential. Default is 1.
#' @param n_pca An integer specifying the number of principal components to use for calculating neighborhoods. For extremely large datasets, using n_pca < 20 allows neighborhoods to be calculated in roughly log(n_samples) time. Default is 100.
#' @param knn_dist A character string specifying the distance metric for k-nearest neighbors. Recommended values: "euclidean, "cosine, "precomputed". Default is "euclidean".
#' @param knn_max An integer specifying the maximum number of neighbors for which alpha decaying kernel is computed for each point. For very large datasets, setting knn_max to a small multiple of knn can speed up computation significantly. Default is NULL.
#' @param t_max An integer specifying the maximum \code{t} to test. Default is 100.
#' @param do_cluster A logical value indicating whether to perform clustering on the PHATE embeddings. Default is FALSE.
#' @param n_clusters An integer specifying the number of clusters to be identified. Default is "auto".
#' @param max_clusters An integer specifying the maximum number of clusters to test. Default is 100.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "phate".
#' @param reduction.key A character string specifying the prefix for the column names of the PHATE embeddings. Default is "PHATE_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the phate.PHATE function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunPHATE(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "phate")
#'
#' @rdname RunPHATE
#' @export
RunPHATE <- function(object, ...) {
  UseMethod(generic = "RunPHATE", object = object)
}

#' @rdname RunPHATE
#' @method RunPHATE Seurat
#' @importFrom Seurat LogSeuratCommand
#' @export
RunPHATE.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                            assay = NULL, slot = "data",
                            n_components = 2, knn = 5, decay = 40, n_landmark = 2000, t = "auto", gamma = 1,
                            n_pca = 100, knn_dist = "euclidean", knn_max = NULL, t_max = 100,
                            do_cluster = FALSE, n_clusters = "auto", max_clusters = 100,
                            reduction.name = "phate", reduction.key = "PHATE_",
                            verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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
    n_pca = n_pca, knn_dist = knn_dist, knn_max = knn_max, t_max = t_max,
    do_cluster = do_cluster, n_clusters = n_clusters, max_clusters = max_clusters,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' @rdname RunPHATE
#' @method RunPHATE default
#' @importFrom reticulate import
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#' @export
RunPHATE.default <- function(object, assay = NULL,
                             n_components = 2, knn = 5, decay = 40, n_landmark = 2000, t = "auto", gamma = 1,
                             n_pca = 100, knn_dist = "euclidean", knn_max = NULL, t_max = 100,
                             do_cluster = FALSE, n_clusters = "auto", max_clusters = 100,
                             reduction.key = "PHATE_", verbose = TRUE, seed.use = 11L, ...) {
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
    random_state = as.integer(seed.use),
    verbose = as.integer(verbose),
    ...
  )
  embedding <- operator$fit_transform(object, t_max = as.integer(t_max))
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)

  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key,
    assay = assay,
    global = TRUE
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

#' Run TriMap (Large-scale Dimensionality Reduction Using Triplets)
#'
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param n_components An integer specifying the number of TriMap components. Default is 2.
#' @param n_inliers An integer specifying the number of nearest neighbors for forming the nearest neighbor triplets. Default is 12.
#' @param n_outliers An integer specifying the number of outliers for forming the nearest neighbor triplets. Default is 4.
#' @param n_random An integer specifying the number of random triplets per point. Default is 3.
#' @param distance_method A character string specifying the distance metric for TriMap. Options are: "euclidean", "manhattan", "angular", "cosine", "hamming". Default is "euclidean".
#' @param lr A numeric value specifying the learning rate for TriMap. Default is 0.1.
#' @param n_iters An integer specifying the number of iterations for TriMap. Default is 400.
#' @param apply_pca A logical value indicating whether to apply PCA before the nearest-neighbor calculation. Default is TRUE.
#' @param opt_method A character string specifying the optimization method for TriMap. Options are: "dbd", "sd", "momentum". Default is "dbd".
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "trimap".
#' @param reduction.key A character string specifying the prefix for the column names of the TriMap embeddings. Default is "TriMap_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the trimap.TRIMAP function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunTriMap(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "trimap")
#'
#' @rdname RunTriMap
#' @export
RunTriMap <- function(object, ...) {
  UseMethod(generic = "RunTriMap", object = object)
}

#' @rdname RunTriMap
#' @method RunTriMap Seurat
#' @importFrom Seurat LogSeuratCommand
#' @export
RunTriMap.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                             assay = NULL, slot = "data",
                             n_components = 2, n_inliers = 12, n_outliers = 4, n_random = 3, distance_method = "euclidean",
                             lr = 0.1, n_iters = 400,
                             apply_pca = TRUE, opt_method = "dbd",
                             reduction.name = "trimap", reduction.key = "TriMap_",
                             verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 2) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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
    lr = lr, n_iters = n_iters,
    apply_pca = apply_pca, opt_method = opt_method,
    reduction.key = reduction.key, verbose = verbose, seed.use = seed.use
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunTriMap
#' @method RunTriMap default
#' @importFrom reticulate import
#' @importFrom Seurat CreateDimReducObject Misc<- Misc
#' @export
RunTriMap.default <- function(object, assay = NULL,
                              n_components = 2, n_inliers = 12, n_outliers = 4, n_random = 3, distance_method = "euclidean",
                              lr = 0.1, n_iters = 400,
                              apply_pca = TRUE, opt_method = "dbd",
                              reduction.key = "TriMap_", verbose = TRUE, seed.use = 11L, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  check_Python("trimap")
  trimap <- import("trimap")

  operator <- trimap$TRIMAP(
    n_dims = as.integer(n_components), n_inliers = as.integer(n_inliers), n_outliers = as.integer(n_outliers), n_random = as.integer(n_random), distance = distance_method,
    lr = lr, n_iters = as.integer(n_iters),
    apply_pca = apply_pca, opt_method = opt_method, verbose = verbose, ...
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

#' Run LargeVis (Dimensionality Reduction with a LargeVis-like method)
#'
#' @inheritParams uwot::lvish
#' @param object An object. This can be a Seurat object or a matrix-like object.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param n_components An integer specifying the number of LargeVis components. Default is 2.
#' @param pca_method Method to carry out any PCA dimensionality reduction when the pca parameter is specified. Allowed values are: "irlba", "rsvd", "bigstatsr", "svd", "auto"(the default. Uses "irlba", unless more than 50 case "svd" is used.)
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "largevis".
#' @param reduction.key A character string specifying the prefix for the column names of the LargeVis embeddings. Default is "LargeVis_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the \link[uwot]{lvish} function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunLargeVis(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "largevis")
#'
#' @rdname RunLargeVis
#' @export
#'
RunLargeVis <- function(object, ...) {
  UseMethod(generic = "RunLargeVis", object = object)
}

#' @rdname RunLargeVis
#' @method RunLargeVis Seurat
#' @importFrom Seurat LogSeuratCommand
#' @export
RunLargeVis.Seurat <- function(object, reduction = "pca", dims = NULL, features = NULL,
                               assay = NULL, slot = "data",
                               perplexity = 50, n_neighbors = perplexity * 3, n_components = 2, metric = "euclidean",
                               n_epochs = -1, learning_rate = 1, scale = "maxabs", init = "lvrandom", init_sdev = NULL,
                               repulsion_strength = 7, negative_sample_rate = 5, nn_method = NULL, n_trees = 50,
                               search_k = 2 * n_neighbors * n_trees, n_threads = NULL, n_sgd_threads = 0, grain_size = 1,
                               kernel = "gauss", pca = NULL, pca_center = TRUE, pcg_rand = TRUE, fast_sgd = FALSE,
                               batch = FALSE, opt_args = NULL, epoch_callback = NULL, pca_method = NULL,
                               reduction.name = "largevis", reduction.key = "LargeVis_",
                               verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features))) == 3) {
    stop("Please specify only one of the following arguments: dims, features")
  }
  if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- as_matrix(t(x = GetAssayData(object = object, slot = slot, assay = assay)[features, , drop = FALSE]))
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

#' @rdname RunLargeVis
#' @method RunLargeVis default
#' @importFrom Seurat CreateDimReducObject
#' @export
RunLargeVis.default <- function(object, assay = NULL,
                                perplexity = 50, n_neighbors = perplexity * 3, n_components = 2, metric = "euclidean",
                                n_epochs = -1, learning_rate = 1, scale = "maxabs", init = "lvrandom", init_sdev = NULL,
                                repulsion_strength = 7, negative_sample_rate = 5, nn_method = NULL, n_trees = 50,
                                search_k = 2 * n_neighbors * n_trees, n_threads = NULL, n_sgd_threads = 0, grain_size = 1,
                                kernel = "gauss", pca = NULL, pca_center = TRUE, pcg_rand = TRUE, fast_sgd = FALSE,
                                batch = FALSE, opt_args = NULL, epoch_callback = NULL, pca_method = NULL,
                                reduction.key = "LargeVis_", verbose = TRUE, seed.use = 11L, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  embedding <- uwot::lvish(
    X = object, perplexity = perplexity, n_neighbors = n_neighbors, n_components = n_components, metric = metric,
    n_epochs = n_epochs, learning_rate = learning_rate, scale = scale, init = init, init_sdev = init_sdev,
    repulsion_strength = repulsion_strength, negative_sample_rate = negative_sample_rate, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k, n_threads = n_threads, n_sgd_threads = n_sgd_threads, grain_size = grain_size,
    kernel = kernel, pca = pca, pca_center = pca_center, pcg_rand = pcg_rand, fast_sgd = fast_sgd,
    verbose = verbose, batch = batch, opt_args = opt_args, epoch_callback = epoch_callback, pca_method = pca_method,
    ...
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

#' Run Force-Directed Layout (Fruchterman-Reingold algorithm)
#'
#' @param object An object. This can be a Seurat object, a Neighbor object, or a Graph object.
#' @param reduction A character string specifying the reduction to be used. Default is NULL.
#' @param dims An integer vector specifying the dimensions to be used. Default is NULL.
#' @param features A character vector specifying the features to be used. Default is NULL.
#' @param assay A character string specifying the assay to be used. Default is NULL.
#' @param slot A character string specifying the slot name to be used. Default is "data".
#' @param graph A character string specifying the name of the Graph object to be used. Default is NULL.
#' @param neighbor A character string specifying the name of the Neighbor object to be used. Default is NULL.
#' @param k.param An integer specifying the number of nearest neighbors to consider. Default is 20.
#' @param ndim An integer specifying the number of dimensions for the force-directed layout. Default is 2.
#' @param niter An integer specifying the number of iterations for the force-directed layout. Default is 500.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "fr".
#' @param reduction.key A character string specifying the prefix for the column names of the force-directed layout embeddings. Default is "FR_".
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the \link[igraph]{layout_with_fr} function.
#'
#' @examples
#' pancreas_sub <- Seurat::FindVariableFeatures(pancreas_sub)
#' pancreas_sub <- RunFR(object = pancreas_sub, features = Seurat::VariableFeatures(pancreas_sub))
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "fr")
#'
#' @rdname RunFR
#' @export
RunFR <- function(object, ...) {
  UseMethod(generic = "RunFR", object = object)
}

#' @rdname RunFR
#' @method RunFR Seurat
#' @importFrom Seurat LogSeuratCommand DefaultAssay
#' @export
RunFR.Seurat <- function(object, reduction = NULL, dims = NULL, features = NULL,
                         assay = NULL, slot = "data",
                         graph = NULL, neighbor = NULL,
                         k.param = 20, ndim = 2, niter = 500,
                         reduction.name = "FR", reduction.key = "FR_",
                         verbose = TRUE, seed.use = 11L, ...) {
  if (sum(c(is.null(x = dims), is.null(x = features), is.null(neighbor), is.null(x = graph))) == 4) {
    stop("Please specify only one of the following arguments: dims, features, neighbor or graph")
  }
  if (!is.null(x = graph)) {
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
  } else if (!is.null(x = features)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- t(GetAssayData(object = object, slot = slot, assay = assay)[features, ])
    data.use <- FindNeighbors(data.use, k.param = k.param)[["snn"]]
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object = object[[reduction]])
    if (max(dims) > ncol(x = data.use)) {
      stop("More dimensions specified in dims than have been computed")
    }
    data.use <- data.use[, dims]
    data.use <- FindNeighbors(data.use, k.param = k.param)[["snn"]]
  } else {
    stop("Please specify one of dims, features, neighbor, or graph")
  }
  object[[reduction.name]] <- RunFR(
    object = data.use, ndim = ndim, niter = niter,
    seed.use = seed.use, verbose = verbose, reduction.key = reduction.key,
    ...
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @rdname RunFR
#' @method RunFR default
#' @importFrom Seurat CreateDimReducObject as.Graph
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr
#' @export
RunFR.default <- function(object, ndim = 2, niter = 500,
                          reduction.key = "FR_", verbose = TRUE, seed.use = 11L, ...) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (inherits(object, "Neighbor")) {
    object <- as.Graph(object)
  }
  g <- graph_from_adjacency_matrix(object, weighted = TRUE)
  embedding <- layout_with_fr(graph = g, dim = ndim, niter = niter, ...)
  colnames(x = embedding) <- paste0(reduction.key, seq_len(ncol(x = embedding)))
  rownames(x = embedding) <- rownames(object)
  reduction <- CreateDimReducObject(
    embeddings = embedding,
    key = reduction.key, assay = object@assay.used, global = TRUE
  )
  return(reduction)
}

#' Run Harmony algorithm
#'
#' This is a modified version of harmony::RunHarmony specifically designed for compatibility with RunSymphonyMap.
#'
#' @param object A Seurat object.
#' @param group.by.vars A character vector specifying the batch variable name.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims.use An integer vector specifying the dimensions to be used. Default is 1:30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "Harmony".
#' @param reduction.key A character string specifying the prefix for the column names of the Harmony embeddings. Default is "Harmony_".
#' @param project.dim A logical value indicating whether to project dimension reduction loadings. Default is TRUE.
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the \link[harmony]{RunHarmony} function.
#'
#' @examples
#' panc8_sub <- Standard_SCP(panc8_sub)
#' panc8_sub <- RunHarmony2(panc8_sub, group.by.vars = "tech", reduction = "Standardpca")
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = "Standardpca")
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = "Harmony")
#'
#' @rdname RunHarmony2
#' @export
RunHarmony2 <- function(object, ...) {
  UseMethod(generic = "RunHarmony2", object = object)
}

#' @rdname RunHarmony2
#' @method RunHarmony2 Seurat
#' @importFrom Seurat Embeddings RunPCA FetchData CreateDimReducObject ProjectDim LogSeuratCommand
#' @importFrom methods slot
#' @importFrom stats sd
#' @export
RunHarmony2.Seurat <- function(object, group.by.vars,
                               reduction = "pca", dims.use = 1:30,
                               project.dim = TRUE,
                               reduction.name = "Harmony", reduction.key = "Harmony_",
                               verbose = TRUE, seed.use = 11L, ...) {
  check_R("harmony@1.1.0")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }

  data.use <- Embeddings(object[[reduction]])
  if (max(dims.use) > ncol(data.use)) {
    stop("trying to use more dimensions than computed")
  }

  assay <- DefaultAssay(object = object[[reduction]])
  metavars_df <- FetchData(object, group.by.vars, cells = rownames(data.use))

  harmonyObject <- harmony::RunHarmony(
    data_mat = data.use[, dims.use, drop = FALSE],
    meta_data = metavars_df,
    vars_use = group.by.vars,
    verbose = verbose,
    return_object = TRUE,
    ...
  )

  harmonyEmbed <- t(as_matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(data.use)
  colnames(harmonyEmbed) <- paste0(reduction.name, "_", seq_len(ncol(harmonyEmbed)))

  harmonyClusters <- t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(data.use)
  colnames(harmonyClusters) <- paste0("R", seq_len(ncol(harmonyClusters)))

  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, sd)),
    assay = assay,
    key = reduction.key,
    misc = list(
      R = harmonyClusters,
      reduction_use = reduction,
      reduction_dims = dims.use
    )
  )

  if (project.dim) {
    object <- ProjectDim(
      object,
      reduction = reduction.name,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  object <- LogSeuratCommand(object = object)
  return(object)
}
