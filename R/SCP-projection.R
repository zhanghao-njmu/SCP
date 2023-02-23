# scArches
NULL

#' Single-cell reference mapping with KNN method
#'
#' @param srt_query
#' @param srt_ref
#' @param ref_umap
#' @param ref_group
#' @param features
#' @param nfeatures
#' @param query_reduction
#' @param ref_reduction
#' @param query_dims
#' @param ref_dims
#' @param projection_method
#' @param nn_method
#' @param k
#' @param distance_metric
#' @param vote_fun
#' @param query_assay
#' @param ref_assay
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "Seurat")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunKNNMap(srt_query = srt_query, srt_ref = srt_ref, ref_umap = "SeuratUMAP2D")
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom Seurat Reductions Embeddings FindVariableFeatures VariableFeatures GetAssayData FindNeighbors CreateDimReducObject DefaultAssay
#' @importFrom SeuratObject as.sparse
#' @importFrom Matrix t
#' @export
RunKNNMap <- function(srt_query, srt_ref, query_assay = NULL, ref_assay = NULL, ref_umap = NULL, ref_group = NULL,
                      features = NULL, nfeatures = 2000,
                      query_reduction = NULL, ref_reduction = NULL, query_dims = 1:30, ref_dims = 1:30,
                      projection_method = c("model", "knn"), nn_method = NULL, k = 30, distance_metric = "cosine", vote_fun = "mean") {
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
    ref_group <- "ref_group"
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(Reductions(srt_ref)[grep("umap", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_umap) == 0) {
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (projection_method == "model" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'")
  }
  simil_method <- c(
    "pearson", "spearman", "cosine", "correlation", "jaccard", "ejaccard", "dice", "edice",
    "hamman", "simple matching", "faith"
  )
  dist_method <- c(
    "euclidean", "chisquared", "kullback", "manhattan", "maximum", "canberra",
    "minkowski", "hamming"
  )
  if (!distance_metric %in% c(simil_method, dist_method)) {
    stop(distance_metric, " method is invalid.")
  }
  if (projection_method == "model") {
    model <- srt_ref[[ref_umap]]@misc$model
    if ("layout" %in% names(model)) {
      if (k != model$config$n_neighbors) {
        k <- model$config$n_neighbors
        message("Set k to ", k, " which is used in the umap model")
      }
    } else if ("embedding" %in% names(model)) {
      if (k != model$n_neighbors) {
        k <- model$n_neighbors
        message("Set k to ", k, " which is used in the umap model")
      }
    }
  }

  if (!is.null(query_reduction) && !is.null(ref_reduction)) {
    message("Use the reduction to calculate distance metric.")
    if (!is.null(query_dims) && !is.null(ref_dims) && length(query_dims) == length(ref_dims)) {
      query <- Embeddings(srt_query, reduction = query_reduction)[, query_dims]
      ref <- Embeddings(srt_ref, reduction = ref_reduction)[, ref_dims]
    } else {
      stop("query_dims and ref_dims must be provided with the same length.")
    }
  } else {
    message("Use the features to calculate distance metric.")
    status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
    message("Detected srt_query data type: ", status_query)
    status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
    message("Detected srt_ref data type: ", status_ref)
    if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
      warning("Data type is unknown or different between srt_query and srt_ref.", immediate. = TRUE)
    }
    if (length(features) == 0) {
      if (length(VariableFeatures(srt_ref, assay = ref_assay)) == 0) {
        srt_ref <- FindVariableFeatures(srt_ref, nfeatures = nfeatures, assay = ref_assay)
      }
      if (length(VariableFeatures(srt_query, assay = query_assay)) == 0) {
        srt_query <- FindVariableFeatures(srt_query, nfeatures = nfeatures, assay = query_assay)
      }
      features <- intersect(VariableFeatures(srt_query, assay = query_assay), VariableFeatures(srt_ref, assay = ref_assay))
    }
    features_common <- Reduce(intersect, list(features, rownames(srt_query[[query_assay]]), rownames(srt_ref[[ref_assay]])))
    message("Use ", length(features_common), " common features to calculate distance.")
    query <- t(GetAssayData(srt_query, slot = "data", assay = query_assay)[features_common, ])
    ref <- t(GetAssayData(srt_ref, slot = "data", assay = ref_assay)[features_common, ])
  }

  if (projection_method == "model" && "layout" %in% names(model) && is.null(ref_group)) {
    srt_query[["ref.umap"]] <- RunUMAP2(object = query, reduction.model = srt_ref[[ref_umap]], assay = query_assay)
    srt_query[["ref.umap"]]@misc[["reduction.model"]] <- ref_umap
    return(srt_query)
  }

  if (is.null(nn_method)) {
    if (as.numeric(nrow(query)) * as.numeric(nrow(ref)) < 1e9) {
      nn_method <- "raw"
    } else {
      nn_method <- "annoy"
    }
  }
  message("Use '", nn_method, "' method to find neighbors.")
  if (!nn_method %in% c("raw", "annoy", "rann")) {
    stop("nn_method must be one of raw, rann and annoy")
  }
  if (nn_method == "annoy" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when nn_method='annoy'")
  }

  if (nn_method %in% c("annoy", "rann")) {
    query.neighbor <- FindNeighbors(
      query = query, object = ref,
      k.param = k, nn.method = nn_method, annoy.metric = distance_metric,
      return.neighbor = TRUE
    )
    match_k <- query.neighbor@nn.idx
    rownames(match_k) <- rownames(query)
    match_k_cell <- apply(match_k, c(1, 2), function(x) rownames(ref)[x])
    knn_cells <- c(match_k_cell)
    match_k_distance <- query.neighbor@nn.dist
    rownames(match_k_distance) <- rownames(query)
    refumap_all <- srt_ref[[ref_umap]]@cell.embeddings[knn_cells, ]
    group <- rep(query.neighbor@cell.names, ncol(query.neighbor@nn.idx))
  } else {
    if (require("RcppParallel", quietly = TRUE)) {
      setThreadOptions()
    }
    if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
      if (distance_metric %in% c("pearson", "spearman")) {
        if (distance_metric == "spearman") {
          ref <- t(apply(ref, 1, rank))
          query <- t(apply(query, 1, rank))
        }
        distance_metric <- "correlation"
      }
      d <- 1 - proxyC::simil(
        x = as.sparse(ref),
        y = as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    } else if (distance_metric %in% dist_method) {
      d <- proxyC::dist(
        x = as.sparse(ref),
        y = as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    }
    if (k == 1) {
      match_k <- as.matrix(apply(d, 2, function(x) order(x, decreasing = FALSE)[1]))
      match_k_cell <- as.matrix(apply(d, 2, function(x) names(x)[order(x, decreasing = FALSE)[1]]))
      match_k_distance <- as.matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1]]))
    } else {
      match_k <- t(as.matrix(apply(d, 2, function(x) order(x, decreasing = FALSE)[1:k])))
      match_k_cell <- t(as.matrix(apply(d, 2, function(x) names(x)[order(x, decreasing = FALSE)[1:k]])))
      match_k_distance <- t(as.matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1:k]])))
    }
    knn_cells <- match_k_cell
    refumap_all <- srt_ref[[ref_umap]]@cell.embeddings[knn_cells, , drop = FALSE]
    group <- rep(colnames(d), k)
  }

  if (projection_method == "model") {
    if ("layout" %in% names(model)) {
      srt_query[["ref.umap"]] <- RunUMAP2(object = query, reduction.model = srt_ref[[ref_umap]], assay = query_assay)
      srt_query[["ref.umap"]]@misc[["reduction.model"]] <- ref_umap
    } else if ("embedding" %in% names(model)) {
      neighborlist <- list(idx = match_k, dist = match_k_distance)
      srt_query[["ref.umap"]] <- RunUMAP2(object = neighborlist, reduction.model = srt_ref[[ref_umap]], assay = query_assay)
      srt_query[["ref.umap"]]@misc[["reduction.model"]] <- ref_umap
    }
  } else {
    refumap <- aggregate(refumap_all, by = list(group), FUN = vote_fun)
    rownames(refumap) <- refumap[, 1]
    refumap[, 1] <- NULL
    colnames(refumap) <- paste0("refUMAP_", seq_len(ncol(refumap)))
    refumap <- as.matrix(refumap)
    srt_query[["ref.umap"]] <- CreateDimReducObject(embeddings = refumap, key = "refUMAP_", assay = query_assay, misc = list(reduction.model = ref_umap))
  }

  if (!is.null(ref_group)) {
    message("Predicting cell types based on ref_group.") ## slow
    level <- as.character(unique(srt_ref[["ref_group", drop = TRUE]]))
    if (k == 1) {
      match_best <- srt_ref[["ref_group", drop = TRUE]][match_k_cell[, 1]]
      names(match_best) <- names(match_k_cell[, 1])
    } else {
      rn <- rownames(match_k_cell)
      match_k_cell <- matrix(srt_ref[["ref_group", drop = TRUE]][match_k_cell],
        nrow = nrow(match_k_cell), ncol = ncol(match_k_cell)
      )
      rownames(match_k_cell) <- rn
      match_freq <- apply(match_k_cell, 1, table)
      if (!inherits(match_freq, "list")) {
        match_freq <- as.list(setNames(object = rep(k, nrow(match_k_cell)), rn))
        match_freq <- lapply(setNames(names(match_freq), names(match_freq)), function(x) setNames(k, match_k_cell[x, 1]))
      }
      match_prob <- lapply(match_freq, function(x) {
        x[level[!level %in% names(x)]] <- 0
        x <- x / sum(x)
        return(x)
      }) %>%
        bind_rows()
      match_prob <- as.matrix(match_prob)
      rownames(match_prob) <- names(match_freq)
      match_best <- apply(match_prob, 1, function(x) names(x)[order(x, decreasing = TRUE)][1])
    }
    srt_query[[paste0("predicted_", ref_group)]] <- match_best[colnames(srt_query)]
  }

  return(srt_query)
}

#' Single-cell reference mapping with PCA method
#'
#' @param srt_query
#' @param srt_ref
#' @param ref_pca
#' @param ref_dims
#' @param ref_umap
#' @param ref_group
#' @param projection_method
#' @param nn_method
#' @param k
#' @param distance_metric
#' @param vote_fun
#' @param query_assay
#' @param ref_assay
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "Seurat")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunPCAMap(srt_query = srt_query, srt_ref = srt_ref, ref_pca = "Seuratpca", ref_umap = "SeuratUMAP2D")
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom Seurat Reductions GetAssayData CreateDimReducObject ProjectUMAP
#' @export
RunPCAMap <- function(srt_query, srt_ref, query_assay = NULL, ref_assay = srt_ref[[ref_pca]]@assay.used,
                      ref_pca = NULL, ref_dims = 1:30, ref_umap = NULL, ref_group = NULL,
                      projection_method = c("model", "knn"), nn_method = NULL, k = 30, distance_metric = "cosine", vote_fun = "mean") {
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
  }

  if (is.null(ref_pca)) {
    ref_pca <- sort(Reductions(srt_ref)[grep("pca", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_pca) == 0) {
      stop("Cannot find PCA reduction in the srt_ref")
    } else {
      message("Set ref_pca to ", ref_pca)
    }
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(Reductions(srt_ref)[grep("umap", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_umap) == 0) {
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (projection_method == "model" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'")
  }

  pca.out <- srt_ref[[ref_pca]]
  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between srt_query and srt_ref.", immediate. = TRUE)
  }

  message("Run PCA projection")
  features <- rownames(pca.out@feature.loadings)
  center <- apply(GetAssayData(object = srt_ref, slot = "data", assay = ref_assay)[features, ], 1, mean)
  names(center) <- features
  sds <- apply(GetAssayData(object = srt_ref, slot = "data", assay = ref_assay)[features, ], 1, sd)
  names(sds) <- features
  rotation <- pca.out@feature.loadings

  features_common <- Reduce(intersect, list(features, rownames(srt_query[[query_assay]]), rownames(srt_ref[[ref_assay]])))
  message("Use ", length(features_common), " common features to calculate PC.")
  query_data <- t(GetAssayData(srt_query, slot = "data", assay = query_assay)[features_common, ])
  query_pca <- scale(query_data[, features_common], center[features_common], sds[features_common]) %*% rotation[features_common, ]
  # ggplot(data = as.data.frame(pca.out@cell.embeddings))+geom_point(aes(x=StandardPC_1,y=StandardPC_2 ))+geom_point(data = as.data.frame(query_pca),mapping = aes(x=StandardPC_1,y=StandardPC_2),color="red")
  srt_query[["ref.pca"]] <- CreateDimReducObject(embeddings = query_pca, key = pca.out@key, assay = query_assay)

  message("Run UMAP projection")
  srt_query <- RunKNNMap(
    srt_query = srt_query, query_assay = query_assay, srt_ref = srt_ref, ref_assay = ref_assay, ref_group = ref_group, ref_umap = ref_umap,
    query_reduction = "ref.pca", ref_reduction = ref_pca, query_dims = ref_dims, ref_dims = ref_dims,
    projection_method = projection_method, nn_method = nn_method, k = k, distance_metric = distance_metric, vote_fun = vote_fun
  )
  return(srt_query)
}

#' Single-cell reference mapping with Seurat method
#'
#' @param srt_query
#' @param srt_ref
#' @param ref_pca
#' @param ref_dims
#' @param ref_umap
#' @param ref_group
#' @param normalization.method
#' @param reduction_project_method
#' @param k.anchor
#' @param k.filter
#' @param k.score
#' @param k.weight
#' @param projection_method
#' @param nn_method
#' @param k
#' @param distance_metric
#' @param vote_fun
#' @param query_assay
#' @param ref_assay
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "Seurat")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunSeuratMap(srt_query = srt_query, srt_ref = srt_ref, ref_pca = "Seuratpca", ref_umap = "SeuratUMAP2D", k.weight = 50)
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom Seurat Reductions FindTransferAnchors TransferData IntegrateEmbeddings ProjectUMAP
#' @export
RunSeuratMap <- function(srt_query, srt_ref, query_assay = NULL, ref_assay = srt_ref[[ref_pca]]@assay.used,
                         ref_pca = NULL, ref_dims = 1:30, ref_umap = NULL, ref_group = NULL,
                         normalization.method = "LogNormalize", reduction_project_method = "pcaproject",
                         k.anchor = 5, k.filter = 200, k.score = 30, k.weight = 100,
                         projection_method = c("model", "knn"), nn_method = NULL, k = 30, distance_metric = "cosine", vote_fun = "mean") {
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  weight.reduction <- switch(reduction_project_method,
    "pcaproject" = "pcaproject",
    "lsiproject" = "lsiproject",
    "rpca" = "pcaproject",
    "cca" = "cca"
  )
  if (is.null(ref_pca)) {
    if (any(grepl("pca", Reductions(srt_ref), ignore.case = TRUE))) {
      ref_pca <- sort(Reductions(srt_ref)[grep("pca", Reductions(srt_ref))])[1]
    } else {
      cat("'ref_pca' is NUll and no pca reduction detected. Run Standard_SCP first.\n")
      srt_ref <- Standard_SCP(srt_ref)
      ref_pca <- "Standardpca"
    }
    cat("Set the ref_pca to '", ref_pca, "'\n", sep = "")
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(Reductions(srt_ref)[grep("umap", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_umap) == 0) {
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (projection_method == "model" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'")
  }

  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between srt_query and srt_ref.", immediate. = TRUE)
  }

  message("Run FindTransferAnchors")
  anchors <- FindTransferAnchors(
    query = srt_query, query.assay = query_assay, reference = srt_ref, normalization.method = normalization.method, dims = ref_dims,
    reference.reduction = ref_pca, reduction = reduction_project_method, reference.assay = ref_assay,
    k.anchor = k.anchor, k.filter = k.filter, k.score = k.score
  )
  if (reduction_project_method != "cca") {
    srt_query <- IntegrateEmbeddings(
      anchorset = anchors, reference = srt_ref, query = srt_query,
      reductions = reduction_project_method, new.reduction.name = "ref.pca",
      weight.reduction = weight.reduction, k.weight = k.weight
    )
  }

  message("Run UMAP projection")
  srt_query <- RunKNNMap(
    srt_query = srt_query, query_assay = query_assay, srt_ref = srt_ref, ref_assay = ref_assay, ref_group = ref_group, ref_umap = ref_umap,
    query_reduction = "ref.pca", ref_reduction = ref_pca, query_dims = ref_dims, ref_dims = ref_dims,
    projection_method = projection_method, nn_method = nn_method, k = k, distance_metric = distance_metric, vote_fun = vote_fun
  )

  return(srt_query)
}

#' Single-cell reference mapping with CSS method
#'
#' @param srt_query
#' @param srt_ref
#' @param ref_css
#' @param ref_umap
#' @param ref_group
#' @param projection_method
#' @param nn_method
#' @param k
#' @param distance_metric
#' @param vote_fun
#' @param query_assay
#' @param ref_assay
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "CSS")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunCSSMap(srt_query = srt_query, srt_ref = srt_ref, ref_css = "CSS", ref_umap = "CSSUMAP2D")
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom Seurat Reductions CreateDimReducObject ProjectUMAP
#' @export
RunCSSMap <- function(srt_query, srt_ref, query_assay = NULL, ref_assay = srt_ref[[ref_css]]@assay.used,
                      ref_css = NULL, ref_umap = NULL, ref_group = NULL,
                      projection_method = c("model", "knn"), nn_method = NULL, k = 30, distance_metric = "cosine", vote_fun = "mean") {
  check_R("quadbiolab/simspec")
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
    ref_group <- "ref_group"
  }
  if (is.null(ref_css)) {
    ref_css <- sort(Reductions(srt_ref)[grep("css", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_css) == 0) {
      stop("Cannot find CSS reduction in the srt_ref")
    } else {
      message("Set ref_css to ", ref_css)
      if (!"model" %in% names(srt_ref[[ref_css]]@misc) || !"sim2profiles" %in% names(srt_ref[[ref_css]]@misc$model)) {
        stop("CSS model is not in the reduction: ", ref_css)
      }
    }
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(Reductions(srt_ref)[grep("umap", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_umap) == 0) {
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (projection_method == "model" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'")
  }

  ref_assay <- srt_ref[[ref_css]]@assay.used
  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between srt_query and srt_ref.", immediate. = TRUE)
  }

  message("Run CSS projection")
  CSSmodel <- srt_ref[[ref_css]]@misc$model
  raw_assay <- DefaultAssay(srt_query)
  DefaultAssay(srt_query) <- query_assay
  srt_query <- simspec::css_project(object = srt_query, model = CSSmodel)
  DefaultAssay(srt_query) <- raw_assay

  message("Run UMAP projection")
  ref_dims <- seq_len(dim(srt_ref[[ref_css]])[2])
  srt_query <- RunKNNMap(
    srt_query = srt_query, query_assay = query_assay, srt_ref = srt_ref, ref_assay = ref_assay, ref_group = ref_group, ref_umap = ref_umap,
    query_reduction = "css_proj", ref_reduction = ref_css, query_dims = ref_dims, ref_dims = ref_dims,
    projection_method = projection_method, nn_method = nn_method, k = k, distance_metric = distance_metric, vote_fun = vote_fun
  )
  return(srt_query)
}

#' Single-cell reference mapping with Symphony method
#'
#' @param srt_query
#' @param srt_ref
#' @param ref_pca
#' @param ref_harmony
#' @param ref_umap
#' @param force
#' @param ref_group
#' @param projection_method
#' @param nn_method
#' @param k
#' @param distance_metric
#' @param vote_fun
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "Harmony")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunSymphonyMap(srt_query = srt_query, srt_ref = srt_ref, ref_pca = "Harmonypca", ref_harmony = "Harmony", ref_umap = "HarmonyUMAP2D")
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom Seurat Reductions GetAssayData DefaultAssay ProjectUMAP
#' @export
RunSymphonyMap <- function(srt_query, srt_ref, query_assay = NULL, ref_assay = srt_ref[[ref_pca]]@assay.used,
                           ref_pca = NULL, ref_harmony = NULL, ref_umap = NULL, ref_group = NULL,
                           projection_method = c("model", "knn"), nn_method = NULL, k = 30, distance_metric = "cosine", vote_fun = "mean") {
  check_R("immunogenomics/symphony")
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  if (!is.null(ref_group)) {
    if (length(ref_group) == ncol(srt_ref)) {
      srt_ref[["ref_group"]] <- ref_group
    } else if (length(ref_group) == 1) {
      if (!ref_group %in% colnames(srt_ref@meta.data)) {
        stop("ref_group must be one of the column names in the meta.data")
      } else {
        srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
      }
    } else {
      stop("Length of ref_group must be one or length of srt_ref.")
    }
    ref_group <- "ref_group"
  }
  if (is.null(ref_pca)) {
    ref_pca <- sort(Reductions(srt_ref)[grep("pca", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_pca) == 0) {
      stop("Cannot find PCA reduction in the srt_ref")
    } else {
      message("Set ref_pca to ", ref_pca)
    }
  }
  if (is.null(ref_harmony)) {
    ref_harmony <- sort(Reductions(srt_ref)[grep("harmony", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_harmony) == 0) {
      stop("Cannot find Harmony reduction in the srt_ref")
    } else {
      message("Set ref_harmony to ", ref_harmony)
    }
  }
  if (is.null(ref_umap)) {
    ref_umap <- sort(Reductions(srt_ref)[grep("umap", Reductions(srt_ref), ignore.case = TRUE)])[1]
    if (length(ref_umap) == 0) {
      stop("Cannot find UMAP reduction in the srt_ref")
    } else {
      message("Set ref_umap to ", ref_umap)
    }
  }
  ref_pca_dims <- srt_ref[[ref_harmony]]@misc$reduction_dims

  projection_method <- match.arg(projection_method)
  if (projection_method == "model" && !"model" %in% names(srt_ref[[ref_umap]]@misc)) {
    message("No UMAP model detected. Set the projection_method to 'knn'")
    projection_method <- "knn"
  }
  if (projection_method == "model" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when projection_method='model'")
  }

  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between srt_query and srt_ref.", immediate. = TRUE)
  }

  message("Build reference")
  ref <- buildReferenceFromSeurat(
    obj = srt_ref,
    assay = ref_assay,
    pca = ref_pca,
    pca_dims = ref_pca_dims,
    harmony = ref_harmony,
    umap = ref_umap
  )
  message("Run mapQuery")
  res <- mapQuery(
    exp_query = GetAssayData(srt_query, slot = "data", assay = query_assay),
    metadata_query = srt_query@meta.data,
    ref_obj = ref,
    vars = NULL,
    sigma = 0.1,
    verbose = TRUE
  )
  Z_pca_query <- res$Z_pca_query
  Zq_corr <- res$Zq_corr
  R_query <- res$R_query

  srt_query[["ref.pca"]] <- CreateDimReducObject(
    embeddings = t(Z_pca_query),
    loadings = ref$loadings,
    stdev = as.numeric(apply(Z_pca_query, 1, stats::sd)),
    assay = query_assay,
    key = "refpca_"
  )
  srt_query[["ref.harmony"]] <- CreateDimReducObject(
    embeddings = t(Zq_corr),
    stdev = as.numeric(apply(Zq_corr, 1, stats::sd)),
    assay = query_assay,
    key = "refharmony_",
    misc = list(R = R_query)
  )

  message("Run UMAP projection")
  ref_dims <- seq_len(dim(srt_ref[[ref_harmony]])[2])
  srt_query <- RunKNNMap(
    srt_query = srt_query, query_assay = query_assay, srt_ref = srt_ref, ref_assay = ref_assay, ref_group = ref_group, ref_umap = ref_umap,
    query_reduction = "ref.harmony", ref_reduction = ref_harmony, query_dims = ref_dims, ref_dims = ref_dims,
    projection_method = projection_method, nn_method = nn_method, k = k, distance_metric = distance_metric, vote_fun = vote_fun
  )

  return(srt_query)
}

#' @importFrom Matrix rowMeans
#' @importFrom Seurat Key Embeddings
#'
buildReferenceFromSeurat <- function(obj, assay = "RNA", pca = "pca", pca_dims = NULL, harmony = "harmony", umap = "umap") {
  if (!assay %in% c("RNA", "SCT")) {
    stop("Only supported assays are RNA or SCT.")
  }
  if (is.null(pca_dims)) {
    pca_dims <- seq_len(ncol(Embeddings(obj, pca)))
  }
  res <- list()
  ## TODO: check that these objects are all correctly initialized
  res$Z_corr <- t(Embeddings(obj, harmony))
  res$Z_orig <- t(Embeddings(obj, pca)[, pca_dims, drop = FALSE])
  message("Saved embeddings")

  res$R <- t(obj[[harmony]]@misc$R)
  message("Saved soft cluster assignments")

  if (assay == "RNA") {
    vargenes_means_sds <- data.frame(
      symbol = obj@assays[[assay]]@var.features,
      mean = rowMeans(obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ])
    )
    vargenes_means_sds$stddev <- symphony::rowSDs(
      A = obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ],
      row_means = vargenes_means_sds$mean
    )
  } else if (assay == "SCT") {
    vargenes_means_sds <- data.frame(
      symbol = obj@assays[[assay]]@var.features,
      mean = rowMeans(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ])
    )
    asdgc <- Matrix(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ], sparse = TRUE)
    vargenes_means_sds$stddev <- symphony::rowSDs(
      asdgc,
      vargenes_means_sds$mean
    )
  }

  res$vargenes_means_sds <- vargenes_means_sds
  message("Saved variable gene information for ", nrow(vargenes_means_sds), " genes.")

  res$loadings <- obj[[pca]]@feature.loadings[, pca_dims, drop = FALSE]
  message("Saved PCA loadings.")

  res$meta_data <- obj@meta.data
  message("Saved metadata.")

  ## Check UMAP
  if (is.null(obj[[umap]]@misc$model)) {
    error("uwot model not initialiazed in Seurat object. Please do RunUMAP with umap.method='uwot', return.model=TRUE first.")
  }
  res$umap <- obj[[umap]]@misc$model

  ## Build Reference!
  message("Calculate final L2 normalized reference centroids (Y_cos)")
  res$centroids <- t(symphony:::cosine_normalize_cpp(V = res$R %*% t(res$Z_corr), dim = 1))
  message("Calculate reference compression terms (Nr and C)")
  res$cache <- symphony:::compute_ref_cache(Rr = res$R, Zr = res$Z_corr)
  colnames(res$Z_orig) <- row.names(res$meta_data)
  rownames(res$Z_orig) <- paste0(Key(obj[[pca]]), seq_len(nrow(res$Z_corr)))
  colnames(res$Z_corr) <- row.names(res$meta_data)
  rownames(res$Z_corr) <- paste0(Key(obj[[harmony]]), seq_len(nrow(res$Z_corr)))
  message("Finished nicely.")
  return(res)
}


#' @importFrom Seurat Embeddings CreateSeuratObject CreateDimReducObject ProjectDim
mapQuery <- function(exp_query, metadata_query, ref_obj, vars = NULL,
                     sigma = 0.1, verbose = TRUE) {
  if (verbose) {
    message("Scaling and synchronizing query gene expression")
  }
  idx_shared_genes <- which(ref_obj$vargenes$symbol %in% rownames(exp_query))
  shared_genes <- ref_obj$vargenes$symbol[idx_shared_genes]
  if (verbose) {
    message("Found ", length(shared_genes), " reference variable genes in query dataset")
  }
  exp_query_scaled <- symphony::scaleDataWithStats(
    exp_query[shared_genes, ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes],
    1
  )
  exp_query_scaled_sync <- matrix(0,
    nrow = length(ref_obj$vargenes$symbol),
    ncol = ncol(exp_query)
  )
  exp_query_scaled_sync[idx_shared_genes, ] <- exp_query_scaled
  rownames(exp_query_scaled_sync) <- ref_obj$vargenes$symbol
  colnames(exp_query_scaled_sync) <- colnames(exp_query)
  if (verbose) {
    message("Project query cells using reference gene loadings")
  }
  Z_pca_query <- t(ref_obj$loadings) %*% exp_query_scaled_sync
  if (verbose) {
    message("Clustering query cells to reference centroids")
  }
  Z_pca_query_cos <- symphony:::cosine_normalize_cpp(V = Z_pca_query, dim = 2)
  R_query <- symphony:::soft_cluster(Y = ref_obj$centroids, Z = Z_pca_query_cos, sigma = sigma)
  if (verbose) {
    message("Correcting query batch effects")
  }
  if (!is.null(vars)) {
    design <- droplevels(metadata_query)[, vars] %>% as.data.frame()
    onehot <- design %>%
      purrr::map(function(.x) {
        if (length(unique(.x)) == 1) {
          rep(1, length(.x))
        } else {
          stats::model.matrix(~ 0 + .x)
        }
      }) %>%
      purrr::reduce(cbind)
    Xq <- cbind(1, intercept = onehot) %>% t()
  } else {
    Xq <- Matrix::Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))),
      sparse = TRUE
    )
  }
  Zq_corr <- symphony:::moe_correct_ref(
    Zq = as.matrix(Z_pca_query),
    Xq = as.matrix(Xq),
    Rq = as.matrix(R_query),
    Nr = as.matrix(ref_obj$cache[[1]]),
    RrZtr = as.matrix(ref_obj$cache[[2]])
  )
  colnames(Z_pca_query) <- row.names(metadata_query)
  rownames(Z_pca_query) <- paste0("PC_", seq_len(nrow(Zq_corr)))
  colnames(Zq_corr) <- row.names(metadata_query)
  rownames(Zq_corr) <- paste0("harmony_", seq_len(nrow(Zq_corr)))

  return(list(Z_pca_query = Z_pca_query, Zq_corr = Zq_corr, R_query = R_query))
}
