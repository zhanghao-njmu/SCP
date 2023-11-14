# StemID
NULL

#' RunKNNPredict
#'
#' This function performs KNN prediction to annotate cell types based on reference scRNA-seq or bulk RNA-seq data.
#'
#' @param srt_query An object of class Seurat to be annotated with cell types.
#' @param srt_ref An object of class Seurat storing the reference cells.
#' @param bulk_ref A cell atlas matrix, where cell types are represented by columns and genes are represented by rows, for example, SCP::ref_scHCL. Either `srt_ref` or `bulk_ref` must be provided.
#' @param query_group A character vector specifying the column name in the `srt_query` metadata that represents the cell grouping.
#' @param ref_group A character vector specifying the column name in the `srt_ref` metadata that represents the cell grouping.
#' @param query_assay A character vector specifying the assay to be used for the query data. Defaults to the default assay of the `srt_query` object.
#' @param ref_assay A character vector specifying the assay to be used for the reference data. Defaults to the default assay of the `srt_ref` object.
#' @param query_reduction A character vector specifying the dimensionality reduction method used for the query data. If NULL, the function will use the default reduction method specified in the `srt_query` object.
#' @param ref_reduction A character vector specifying the dimensionality reduction method used for the reference data. If NULL, the function will use the default reduction method specified in the `srt_ref` object.
#' @param query_dims A numeric vector specifying the dimensions to be used for the query data. Defaults to the first 30 dimensions.
#' @param ref_dims A numeric vector specifying the dimensions to be used for the reference data. Defaults to the first 30 dimensions.
#' @param query_collapsing A boolean value indicating whether the query data should be collapsed to group-level average expression values. If TRUE, the function will calculate the average expression values for each group in the query data and the annotation will be performed separately for each group. Otherwise it will use the raw expression values for each cell.
#' @param ref_collapsing A boolean value indicating whether the reference data should be collapsed to group-level average expression values. If TRUE, the function will calculate the average expression values for each group in the reference data and the annotation will be performed separately for each group. Otherwise it will use the raw expression values for each cell.
#' @param return_full_distance_matrix A boolean value indicating whether the full distance matrix should be returned. If TRUE, the function will return the distance matrix used for the KNN prediction, otherwise it will only return the annotated cell types.
#' @param features A character vector specifying the features (genes) to be used for the KNN prediction. If NULL, all the features in the query and reference data will be used.
#' @param features_type A character vector specifying the type of features to be used for the KNN prediction. Must be one of "HVF" (highly variable features) or "DE" (differentially expressed features). Defaults to "HVF".
#' @param feature_source A character vector specifying the source of the features to be used for the KNN prediction. Must be one of "both", "query", or "ref". Defaults to "both".
#' @param nfeatures An integer specifying the maximum number of features to be used for the KNN prediction. Defaults to 2000.
#' @param DEtest_param A list of parameters to be passed to the differential expression test function if `features_type` is set to "DE". Defaults to `list(max.cells.per.ident = 200, test.use = "wilcox")`.
#' @param DE_threshold Threshold used to filter the DE features. Default is "p_val < 0.05". If using "roc" test, \code{DE_threshold} should be needs to be reassigned. e.g. "power > 0.5".
#' @param nn_method A character vector specifying the method to be used for finding nearest neighbors. Must be one of "raw", "rann", or "annoy". Defaults to "raw".
#' @param distance_metric A character vector specifying the distance metric to be used for calculating similarity between cells. Must be one of "cosine", "euclidean", "manhattan", or "hamming". Defaults to "cosine".
#' @param k An integer specifying the number of nearest neighbors to be considered for the KNN prediction. Defaults to 30.
#' @param filter_lowfreq An integer specifying the threshold for filtering low-frequency cell types from the predicted results. Cell types with a frequency lower than `filter_lowfreq` will be labelled as "unreliable". Defaults to 0, which means no filtering will be performed.
#' @param prefix A character vector specifying the prefix to be added to the resulting annotations. Defaults to "KNNPredict".
#'
#' @examples
#' # Annotate cells using bulk RNA-seq data
#' data("pancreas_sub")
#' data("ref_scMCA")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' pancreas_sub <- RunKNNPredict(srt_query = pancreas_sub, bulk_ref = ref_scMCA)
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#'
#' # Removal of low credible cell types from the predicted results
#' pancreas_sub <- RunKNNPredict(srt_query = pancreas_sub, bulk_ref = ref_scMCA, filter_lowfreq = 30)
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#'
#' # Annotate clusters using bulk RNA-seq data
#' pancreas_sub <- RunKNNPredict(srt_query = pancreas_sub, query_group = "SubCellType", bulk_ref = ref_scMCA)
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#'
#' # Annotate using single cell RNA-seq data
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
#' panc8_sub <- RenameFeatures(panc8_sub, newnames = genenames)
#' panc8_sub <- check_srtMerge(panc8_sub, batch = "tech")[["srtMerge"]]
#'
#' pancreas_sub <- RunKNNPredict(srt_query = pancreas_sub, srt_ref = panc8_sub, ref_group = "celltype")
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "KNNPredict_simil")
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   ref_group = "celltype", ref_collapsing = FALSE
#' )
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "KNNPredict_prob")
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   query_group = "SubCellType", ref_group = "celltype"
#' )
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "KNNPredict_simil")
#'
#' # Annotate with DE gene instead of HVF
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   features_type = "DE", feature_source = "ref"
#' )
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "KNNPredict_simil")
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   query_group = "SubCellType", ref_group = "celltype",
#'   features_type = "DE", feature_source = "both"
#' )
#' CellDimPlot(pancreas_sub, group.by = "KNNPredict_classification", label = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "KNNPredict_simil")
#'
#' @importFrom methods as
#' @importFrom Matrix t colSums rowSums
#' @importFrom Seurat DefaultAssay GetAssayData FindVariableFeatures VariableFeatures AverageExpression FindNeighbors as.sparse
#' @importFrom proxyC simil dist
#' @export
#'
RunKNNPredict <- function(srt_query, srt_ref = NULL, bulk_ref = NULL,
                          query_group = NULL, ref_group = NULL,
                          query_assay = NULL, ref_assay = NULL,
                          query_reduction = NULL, ref_reduction = NULL,
                          query_dims = 1:30, ref_dims = 1:30,
                          query_collapsing = !is.null(query_group), ref_collapsing = TRUE,
                          return_full_distance_matrix = FALSE,
                          features = NULL, features_type = c("HVF", "DE"), feature_source = "both", nfeatures = 2000,
                          DEtest_param = list(max.cells.per.ident = 200, test.use = "wilcox"),
                          DE_threshold = "p_val_adj < 0.05",
                          nn_method = NULL, distance_metric = "cosine", k = 30,
                          filter_lowfreq = 0, prefix = "KNNPredict") {
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  features_type <- match.arg(features_type, choices = c("HVF", "DE"))
  if (is.null(query_reduction) + is.null(ref_reduction) == 1) {
    stop("query_reduction and ref_reduction must be both provided")
  }
  if (is.null(query_reduction) + is.null(ref_reduction) == 2) {
    use_reduction <- FALSE
  } else {
    use_reduction <- TRUE
  }
  if (!isTRUE(use_reduction)) {
    if (length(features) == 0) {
      if (features_type == "HVF" && feature_source %in% c("both", "query")) {
        if (length(VariableFeatures(srt_query, assay = query_assay)) == 0) {
          srt_query <- FindVariableFeatures(srt_query, nfeatures = nfeatures, assay = query_assay)
        }
        features_query <- VariableFeatures(srt_query, assay = query_assay)
      } else if (features_type == "DE" && feature_source %in% c("both", "query")) {
        if (is.null(query_group)) {
          stop("'query_group' must be provided when 'features_type' is 'DE' and 'feature_source' is 'both' or 'query'")
        } else {
          slot <- paste0("DEtest_", query_group)
          DEtest_param[["force"]] <- TRUE
          if (!slot %in% names(srt_query@tools) || length(grep(pattern = "AllMarkers", names(srt_query@tools[[slot]]))) == 0) {
            srt_query <- do.call(RunDEtest, c(list(srt = srt_query, group_by = query_group), DEtest_param))
          }
          if ("test.use" %in% names(DEtest_param)) {
            test.use <- DEtest_param[["test.use"]]
          } else {
            test.use <- "wilcox"
          }
          index <- grep(pattern = paste0("AllMarkers_", test.use), names(srt_query@tools[[slot]]))[1]
          de <- names(srt_query@tools[[slot]])[index]
          message("Use the DE features from ", de, " to calculate distance metric.")
          de_df <- srt_query@tools[[slot]][[de]]
          de_df <- de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), , drop = FALSE]
          rownames(de_df) <- seq_len(nrow(de_df))
          de_df <- de_df[order(de_df[["avg_log2FC"]], decreasing = TRUE), , drop = FALSE]
          de_top <- de_df[!duplicated(de_df[["gene"]]), , drop = FALSE]
          stat <- sort(table(de_top$group1))
          stat <- stat[stat > 0]
          mat <- matrix(FALSE, nrow = max(stat), ncol = length(stat))
          colnames(mat) <- names(stat)
          for (g in names(stat)) {
            mat[1:stat[g], g] <- TRUE
          }
          nfeatures <- sum(cumsum(rowSums(mat)) <= nfeatures)
          if (test.use == "roc") {
            features_query <- unlist(by(de_top, list(de_top[["group1"]]), function(x) {
              x <- x[order(x[["power"]], decreasing = TRUE), , drop = FALSE]
              head(x[["gene"]], nfeatures)
            }))
          } else {
            features_query <- unlist(by(de_top, list(de_top[["group1"]]), function(x) {
              x <- x[order(x[["p_val"]], decreasing = FALSE), , drop = FALSE]
              head(x[["gene"]], nfeatures)
            }))
          }
          message("DE features number of the query data: ", length(features_query))
        }
      } else {
        features_query <- rownames(srt_query[[query_assay]])
      }
    }
  }

  if (!is.null(bulk_ref)) {
    if (length(features) == 0) {
      features <- features_query
    }
    features_common <- intersect(features, rownames(bulk_ref))
    message("Use ", length(features_common), " features to calculate distance.")
    ref <- t(bulk_ref[features_common, , drop = FALSE])
  } else if (!is.null(srt_ref)) {
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
    } else {
      stop("ref_group must be provided.")
    }
    drop_cell <- colnames(srt_ref)[is.na(srt_ref[["ref_group", drop = TRUE]])]
    if (length(drop_cell) > 0) {
      message("Drop ", length(drop_cell), " cells with NA in the ref_group")
      srt_ref <- srt_ref[, setdiff(colnames(srt_ref), drop_cell)]
    }
    if (isTRUE(use_reduction)) {
      message("Use the reduction to calculate distance metric.")
      if (!is.null(query_dims) && !is.null(ref_dims) && length(query_dims) == length(ref_dims)) {
        query <- Embeddings(srt_query, reduction = query_reduction)[, query_dims]
        ref <- Embeddings(srt_ref, reduction = ref_reduction)[, ref_dims]
      } else {
        stop("query_dims and ref_dims must be provided with the same length.")
      }
    } else {
      if (length(features) == 0) {
        if (features_type == "HVF" && feature_source %in% c("both", "ref")) {
          message("Use the HVF to calculate distance metric.")
          if (length(VariableFeatures(srt_ref, assay = ref_assay)) == 0) {
            srt_ref <- FindVariableFeatures(srt_ref, nfeatures = nfeatures, assay = ref_assay)
          }
          features_ref <- VariableFeatures(srt_ref, assay = ref_assay)
        } else if (features_type == "DE" && feature_source %in% c("both", "ref")) {
          slot <- paste0("DEtest_", ref_group)
          DEtest_param[["force"]] <- TRUE
          if (!slot %in% names(srt_ref@tools) || length(grep(pattern = "AllMarkers", names(srt_ref@tools[[slot]]))) == 0) {
            srt_ref <- do.call(RunDEtest, c(list(srt = srt_ref, group_by = ref_group), DEtest_param))
          }
          if ("test.use" %in% names(DEtest_param)) {
            test.use <- DEtest_param[["test.use"]]
          } else {
            test.use <- "wilcox"
          }
          index <- grep(pattern = paste0("AllMarkers_", test.use), names(srt_ref@tools[[slot]]))[1]
          de <- names(srt_ref@tools[[slot]])[index]
          message("Use the DE features from ", de, " to calculate distance metric.")
          de_df <- srt_ref@tools[[slot]][[de]]
          de_df <- de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), , drop = FALSE]
          rownames(de_df) <- seq_len(nrow(de_df))
          de_df <- de_df[order(de_df[["avg_log2FC"]], decreasing = TRUE), , drop = FALSE]
          de_top <- de_df[!duplicated(de_df[["gene"]]), , drop = FALSE]
          stat <- sort(table(de_top$group1))
          stat <- stat[stat > 0]
          mat <- matrix(FALSE, nrow = max(stat), ncol = length(stat))
          colnames(mat) <- names(stat)
          for (g in names(stat)) {
            mat[1:stat[g], g] <- TRUE
          }
          nfeatures <- sum(cumsum(rowSums(mat)) <= nfeatures)
          if (test.use == "roc") {
            features_ref <- unlist(by(de_top, list(de_top[["group1"]]), function(x) {
              x <- x[order(x[["power"]], decreasing = TRUE), , drop = FALSE]
              head(x[["gene"]], nfeatures)
            }))
          } else {
            features_ref <- unlist(by(de_top, list(de_top[["group1"]]), function(x) {
              x <- x[order(x[["p_val"]], decreasing = FALSE), , drop = FALSE]
              head(x[["gene"]], nfeatures)
            }))
          }
          message("DE features number of the ref data: ", length(features_ref))
        } else {
          features_ref <- rownames(srt_ref[[ref_assay]])
        }
        features <- intersect(features_ref, features_query)
      }
      features_common <- Reduce(intersect, list(features, rownames(srt_query[[query_assay]]), rownames(srt_ref[[ref_assay]])))
      message("Use ", length(features_common), " features to calculate distance.")
      if (isTRUE(ref_collapsing)) {
        ref <- AverageExpression(object = srt_ref, features = features_common, slot = "data", assays = ref_assay, group.by = "ref_group", verbose = FALSE)[[1]]
        ref <- t(log1p(ref))
      } else {
        ref <- t(GetAssayData(srt_ref, slot = "data", assay = ref_assay)[features_common, ])
      }
    }
  } else {
    stop("srt_ref or bulk_ref must be provided at least one.")
  }

  if (!inherits(ref, "matrix")) {
    ref <- as_matrix(ref)
  }
  k <- min(c(k, nrow(ref)))

  if (!is.null(query_group)) {
    if (length(query_group) == ncol(srt_query)) {
      srt_query[["query_group"]] <- query_group
    } else if (length(query_group) == 1) {
      if (!query_group %in% colnames(srt_query@meta.data)) {
        stop("query_group must be one of the column names in the meta.data")
      } else {
        srt_query[["query_group"]] <- srt_query[[query_group]]
      }
    } else {
      stop("Length of query_group must be one or length of srt_query.")
    }
  }
  if (!isTRUE(use_reduction)) {
    query_assay <- query_assay %||% DefaultAssay(srt_query)
    if (isTRUE(query_collapsing)) {
      if (is.null(query_group)) {
        stop("query_group must be provided when query_collapsing is TRUE.")
      }
      query <- AverageExpression(object = srt_query, features = colnames(ref), slot = "data", assays = query_assay, group.by = "query_group", verbose = FALSE)[[1]]
      query <- t(log1p(query))
    } else {
      query <- t(GetAssayData(srt_query, slot = "data", assay = query_assay)[colnames(ref), , drop = FALSE])
    }
  }

  if (!isTRUE(use_reduction)) {
    status_dat <- check_DataType(data = query)
    message("Detected query data type: ", status_dat)
    status_ref <- check_DataType(data = ref)
    message("Detected reference data type: ", status_ref)
    if (status_ref != status_dat || any(status_dat == "unknown", status_ref == "unknown")) {
      warning("Data type is unknown or different between query and reference.", immediate. = TRUE)
    }
  }

  message("Calculate similarity...")
  # tst.expr <- matrix(data = 0, nrow = nrow(ref), ncol = dim(query)[2])
  # rownames(tst.expr) <- rownames(ref)
  # colnames(tst.expr) <- colnames(query)
  # features_common <- intersect(rownames(query), rownames(ref.expr))
  # tst.expr[features_common, ] <- as_matrix(query[features_common, ])
  #
  # library(RcppArmadillo)
  # library(RcppXPtrUtils)
  # library(parallelDist)
  #   query <- cbind(ref.expr, tst.expr)
  #   query <- apply(query, 2, function(x) x - mean(x))
  #   CosineCPP <- cppXPtr(
  #     "double customDist(const arma::mat &A, const arma::mat &B) {
  #   double similarity = arma::as_scalar(arma::dot(A, B))/(arma::as_scalar(arma::norm(A))*arma::as_scalar(arma::norm(B))) ;
  #   return(similarity);
  # }",
  #     depends = c("RcppArmadillo")
  #   )
  #   cors <- as_matrix(parDist(t(query), method = "custom", func = CosineCPP))
  #   cors <- cors[colnames(ref.expr), colnames(tst.expr)]
  #

  if (is.null(nn_method)) {
    if (as.numeric(nrow(query)) * as.numeric(nrow(ref)) >= 1e8) {
      nn_method <- "annoy"
    } else {
      nn_method <- "raw"
    }
  }
  message("Use '", nn_method, "' method to find neighbors.")
  if (!nn_method %in% c("raw", "annoy", "rann")) {
    stop("nn_method must be one of raw, rann and annoy")
  }
  if (nn_method == "annoy" && !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")) {
    stop("distance_metric must be one of euclidean, cosine, manhattan, and hamming when nn_method='annoy'")
  }
  if (isTRUE(return_full_distance_matrix) && nn_method != "raw") {
    warning("Distance matrix will not be returned besause nn_method is not 'raw'", immediate. = TRUE)
    return_full_distance_matrix <- FALSE
  }
  simil_method <- c(
    "cosine", "pearson", "spearman", "correlation", "jaccard", "ejaccard", "dice", "edice",
    "hamman", "simple matching", "faith"
  )
  dist_method <- c(
    "euclidean", "chisquared", "kullback", "manhattan", "maximum", "canberra",
    "minkowski", "hamming"
  )
  if (!distance_metric %in% c(simil_method, dist_method)) {
    stop(distance_metric, " method is invalid.")
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
    match_k_distance <- query.neighbor@nn.dist
    rownames(match_k_distance) <- rownames(query)
  } else {
    if (requireNamespace("RcppParallel", quietly = TRUE)) {
      RcppParallel::setThreadOptions()
    }
    if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
      if (distance_metric %in% c("pearson", "spearman")) {
        if (distance_metric == "spearman") {
          ref <- t(apply(ref, 1, rank))
          query <- t(apply(query, 1, rank))
        }
        distance_metric <- "correlation"
      }
      d <- 1 - simil(
        x = as.sparse(ref),
        y = as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    } else if (distance_metric %in% dist_method) {
      d <- dist(
        x = as.sparse(ref),
        y = as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    }
    if (k == 1) {
      match_k_cell <- as_matrix(apply(d, 2, function(x) names(x)[order(x, decreasing = FALSE)[1]]))
      match_k_distance <- as_matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1]]))
    } else {
      match_k_cell <- t(as_matrix(apply(d, 2, function(x) names(x)[order(x, decreasing = FALSE)[1:k]])))
      match_k_distance <- t(as_matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1:k]])))
    }
  }

  message("Predict cell type...")
  match_prob <- NULL
  if (!is.null(srt_ref) && (isFALSE(ref_collapsing) || isTRUE(use_reduction))) {
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
      match_prob <- do.call(rbind, lapply(match_freq, function(x) {
        x[level[!level %in% names(x)]] <- 0
        x <- x / sum(x)
        return(x)
      }))
      match_prob <- as_matrix(match_prob)
      rownames(match_prob) <- names(match_freq)
      match_best <- apply(match_prob, 1, function(x) names(x)[order(x, decreasing = TRUE)][1])
    }
  } else {
    match_best <- match_k_cell[, 1]
  }

  result <- list(
    features = features, nn_method = nn_method, distance_metric = distance_metric, k = k,
    other_params = list(
      query_group = query_group, query_reduction = query_reduction, query_assay = query_assay, query_dims = query_dims, query_collapsing = query_collapsing,
      ref_group = ref_group, ref_reduction = ref_reduction, ref_assay = ref_assay, ref_dims = ref_dims, ref_collapsing = ref_collapsing
    )
  )
  result[["match_best"]] <- match_best
  if (!is.null(match_prob)) {
    result[["match_prob"]] <- match_prob
  }
  result[["match_k_cell"]] <- match_k_cell
  result[["match_k_distance"]] <- match_k_distance
  if (isTRUE(return_full_distance_matrix)) {
    result[["distance_matrix"]] <- d[seq_len(nrow(d)), , drop = FALSE]
  }

  srt_query@tools[[paste0(prefix, "_classification")]] <- result

  if (isTRUE(query_collapsing)) {
    query_index <- as.character(srt_query[["query_group", drop = TRUE]])
  } else {
    query_index <- colnames(srt_query)
  }
  srt_query[[paste0(prefix, "_classification")]] <- match_best[query_index]
  if (!is.null(match_prob)) {
    srt_query[[paste0(prefix, "_prob")]] <- apply(match_prob, 1, max)[query_index]
  } else {
    distance <- match_k_distance[, 1]
    # srt_query[[paste0(prefix, "_score")]] <- ((max(distance) - distance) / diff(range(distance)))[query_index]
    if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
      srt_query[[paste0(prefix, "_simil")]] <- (1 - distance)[query_index]
    } else {
      srt_query[[paste0(prefix, "_dist")]] <- distance[query_index]
    }
  }

  if (is.numeric(filter_lowfreq) && filter_lowfreq > 0) {
    drop <- table(srt_query[[paste0(prefix, "_classification"), drop = TRUE]])
    drop <- names(drop)[drop <= filter_lowfreq]
    srt_query[[paste0(prefix, "_classification"), drop = TRUE]][srt_query[[paste0(prefix, "_classification"), drop = TRUE]] %in% drop] <- "unreliable"
  }

  return(srt_query)
}

#' Annotate single cells using scmap.
#'
#' @inheritParams RunKNNPredict
#' @param method The method to be used for scmap analysis. Can be any of "scmapCluster" or "scmapCell". The default value is "scmapCluster".
#' @param nfeatures The number of top features to be selected. The default value is 500.
#' @param threshold The threshold value on similarity to determine if a cell is assigned to a cluster. This should be a value between 0 and 1. The default value is 0.5.
#' @param k Number of clusters per group for k-means clustering when method is "scmapCell".
#'
#' @examples
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
#' panc8_sub <- RenameFeatures(panc8_sub, newnames = genenames)
#' panc8_sub <- check_srtMerge(panc8_sub, batch = "tech")[["srtMerge"]]
#'
#' # Annotation
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' pancreas_sub <- RunScmap(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   ref_group = "celltype", method = "scmapCluster"
#' )
#' CellDimPlot(pancreas_sub, group.by = "scmap_annotation")
#'
#' pancreas_sub <- RunScmap(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   ref_group = "celltype", method = "scmapCell"
#' )
#' CellDimPlot(pancreas_sub, group.by = "scmap_annotation")
#'
#' @importFrom Seurat GetAssayData
#' @export
RunScmap <- function(srt_query, srt_ref, ref_group = NULL, query_assay = "RNA", ref_assay = "RNA",
                     method = "scmapCluster", nfeatures = 500, threshold = 0.5, k = 10) {
  check_R("scmap")
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
  } else {
    stop("'ref_group' must be provided.")
  }

  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between query and ref.", immediate. = TRUE)
  }

  assays_query <- list(
    counts = GetAssayData(object = srt_query, assay = query_assay, slot = "counts"),
    logcounts = GetAssayData(object = srt_query, assay = query_assay, slot = "data")
  )
  sce_query <- as(SummarizedExperiment::SummarizedExperiment(assays = assays_query), Class = "SingleCellExperiment")
  SummarizedExperiment::rowData(sce_query)[["feature_symbol"]] <- rownames(sce_query)
  metadata_query <- srt_query[[]]
  SummarizedExperiment::colData(x = sce_query) <- S4Vectors::DataFrame(metadata_query)

  assays_ref <- list(
    counts = GetAssayData(object = srt_ref, assay = ref_assay, slot = "counts"),
    logcounts = GetAssayData(object = srt_ref, assay = ref_assay, slot = "data")
  )
  sce_ref <- as(SummarizedExperiment::SummarizedExperiment(assays = assays_ref), Class = "SingleCellExperiment")
  SummarizedExperiment::rowData(sce_ref)[["feature_symbol"]] <- rownames(sce_ref)
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(metadata_ref)

  message("Perform selectFeatures on the data...")
  sce_ref <- scmap::selectFeatures(sce_ref, n_features = nfeatures, suppress_plot = TRUE)
  features <- rownames(sce_ref)[SummarizedExperiment::rowData(sce_ref)$scmap_features]

  if (method == "scmapCluster") {
    message("Perform indexCluster on the data...")
    sce_ref <- scmap::indexCluster(sce_ref, cluster_col = ref_group)
    message("Perform scmapCluster on the data...")
    scmapCluster_results <- scmap::scmapCluster(
      projection = sce_query,
      index_list = list(S4Vectors::metadata(sce_ref)$scmap_cluster_index),
      threshold = threshold
    )
    if (!"scmap_cluster_labs" %in% names(scmapCluster_results)) {
      stop("scmap failed to run. Please check the warning message.")
    }
    srt_query@tools[["scmapCluster_results"]] <- scmapCluster_results
    srt_query$scmap_annotation <- scmapCluster_results$scmap_cluster_labs[, 1]
    srt_query$scmap_score <- scmapCluster_results$scmap_cluster_siml[, 1]
  } else if (method == "scmapCell") {
    message("Perform indexCell on the data...")
    sce_ref <- scmap::indexCell(sce_ref,
      M = ifelse(nfeatures <= 1000, nfeatures / 10, 100),
      k = sqrt(ncol(sce_ref))
    )
    message("Perform scmapCell on the data...")
    scmapCell_results <- scmap::scmapCell(
      projection = sce_query,
      index_list = list(S4Vectors::metadata(sce_ref)$scmap_cell_index),
      w = k
    )
    if (!"cells" %in% names(scmapCell_results[[1]])) {
      stop("scmap failed to run. Please check the warning message.")
    }
    srt_query@tools[["scmapCell_results"]] <- scmapCell_results[[1]]
    message("Perform scmapCell2Cluster on the data...")
    scmapCell_clusters <- scmap::scmapCell2Cluster(
      scmapCell_results = scmapCell_results,
      cluster_list = list(
        as.character(SummarizedExperiment::colData(sce_ref)[[ref_group]])
      ), w = k,
      threshold = threshold
    )
    srt_query@tools[["scmapCell_results"]][["scmapCell2Cluster"]] <- scmapCell_clusters
    srt_query$scmap_annotation <- scmapCell_clusters$scmap_cluster_labs[, 1]
    srt_query$scmap_score <- scmapCell_clusters$scmap_cluster_siml[, 1]
  }
  return(srt_query)
}

#' Annotate single cells using SingleR
#'
#' @inheritParams RunKNNPredict
#' @inheritParams SingleR::SingleR
#' @inheritParams SingleR::trainSingleR
#' @param genes "genes" parameter in \code{\link[SingleR]{SingleR}} function.
#' @param de.method "de.method" parameter in \code{\link[SingleR]{SingleR}} function.
#'
#' @examples
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
#' panc8_sub <- RenameFeatures(panc8_sub, newnames = genenames)
#' panc8_sub <- check_srtMerge(panc8_sub, batch = "tech")[["srtMerge"]]
#'
#' # Annotation
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' pancreas_sub <- RunSingleR(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   query_group = "Standardclusters", ref_group = "celltype",
#' )
#' CellDimPlot(pancreas_sub, group.by = "singler_annotation")
#'
#' pancreas_sub <- RunSingleR(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   query_group = NULL, ref_group = "celltype"
#' )
#' CellDimPlot(pancreas_sub, group.by = "singler_annotation")
#'
#' @importFrom Seurat GetAssayData
#' @export
RunSingleR <- function(srt_query, srt_ref, query_group = NULL, ref_group = NULL,
                       query_assay = "RNA", ref_assay = "RNA",
                       genes = "de", de.method = "wilcox", sd.thresh = 1, de.n = NULL,
                       aggr.ref = FALSE, aggr.args = list(),
                       quantile = 0.8, fine.tune = TRUE, tune.thresh = 0.05, prune = TRUE,
                       BPPARAM = BiocParallel::bpparam()) {
  check_R("SingleR")
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
  } else {
    stop("'ref_group' must be provided.")
  }
  if (!is.null(query_group)) {
    if (length(query_group) == ncol(srt_query)) {
      srt_query[["query_group"]] <- query_group
    } else if (length(query_group) == 1) {
      if (!query_group %in% colnames(srt_query@meta.data)) {
        stop("query_group must be one of the column names in the meta.data")
      } else {
        srt_query[["query_group"]] <- srt_query[[query_group]]
      }
    } else {
      stop("Length of query_group must be one or length of srt_query.")
    }
    query_group <- "query_group"
    method <- "SingleRCluster"
  } else {
    method <- "SingleRCell"
  }

  status_query <- check_DataType(data = GetAssayData(srt_query, slot = "data", assay = query_assay))
  message("Detected srt_query data type: ", status_query)
  status_ref <- check_DataType(data = GetAssayData(srt_ref, slot = "data", assay = ref_assay))
  message("Detected srt_ref data type: ", status_ref)
  if (status_ref != status_query || any(status_query == "unknown", status_ref == "unknown")) {
    warning("Data type is unknown or different between query and ref.", immediate. = TRUE)
  }

  assays_query <- list(
    counts = GetAssayData(object = srt_query, assay = query_assay, slot = "counts"),
    logcounts = GetAssayData(object = srt_query, assay = query_assay, slot = "data")
  )
  sce_query <- as(SummarizedExperiment::SummarizedExperiment(assays = assays_query), Class = "SingleCellExperiment")
  metadata_query <- srt_query[[]]
  SummarizedExperiment::colData(x = sce_query) <- S4Vectors::DataFrame(metadata_query)

  assays_ref <- list(
    counts = GetAssayData(object = srt_ref, assay = ref_assay, slot = "counts"),
    logcounts = GetAssayData(object = srt_ref, assay = ref_assay, slot = "data")
  )
  sce_ref <- as(SummarizedExperiment::SummarizedExperiment(assays = assays_ref), Class = "SingleCellExperiment")
  metadata_ref <- srt_ref[[]]
  SummarizedExperiment::colData(x = sce_ref) <- S4Vectors::DataFrame(metadata_ref)

  if (method == "SingleRCluster") {
    message("Perform ", method, " on the data...")
    SingleRCluster_results <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = factor(SummarizedExperiment::colData(sce_ref)[[ref_group]]),
      clusters = factor(SummarizedExperiment::colData(sce_query)[[query_group]]),
      de.method = de.method, genes = genes, sd.thresh = sd.thresh, de.n = de.n,
      aggr.ref = aggr.ref, aggr.args = aggr.args,
      quantile = quantile, fine.tune = fine.tune, tune.thresh = tune.thresh, prune = prune,
      BPPARAM = BPPARAM
    )
    names(SingleRCluster_results$labels) <- levels(factor(SummarizedExperiment::colData(sce_query)[[query_group]]))
    rownames(SingleRCluster_results$scores) <- levels(factor(SummarizedExperiment::colData(sce_query)[[query_group]]))
    if (isTRUE(prune)) {
      names(SingleRCluster_results$pruned.labels) <- levels(factor(SummarizedExperiment::colData(sce_query)[[query_group]]))
    }
    srt_query$singler_annotation <- if (isTRUE(prune)) {
      SingleRCluster_results$pruned.labels[as.character(srt_query$query_group)]
    } else {
      SingleRCluster_results$labels[as.character(srt_query$query_group)]
    }
    srt_query$singler_score <- sapply(as.character(unique(srt_query$query_group)), FUN = function(x) {
      if (isTRUE(prune)) {
        y <- SingleRCluster_results$pruned.labels[x]
      } else {
        y <- SingleRCluster_results$labels[x]
      }
      if (is.na(y)) {
        out <- NA
      } else {
        out <- SingleRCluster_results$scores[x, y]
      }
      return(out)
    })[as.character(srt_query$query_group)]
  } else if (method == "SingleRCell") {
    message("Perform ", method, " on the data...")
    SingleRCell_results <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = factor(SummarizedExperiment::colData(sce_ref)[[ref_group]]),
      clusters = NULL,
      de.method = de.method, genes = genes, sd.thresh = sd.thresh, de.n = de.n,
      aggr.ref = aggr.ref, aggr.args = aggr.args,
      quantile = quantile, fine.tune = fine.tune, tune.thresh = tune.thresh, prune = prune,
      BPPARAM = BPPARAM
    )
    srt_query$singler_annotation <- if (isTRUE(prune)) {
      SingleRCell_results$pruned.labels
    } else {
      SingleRCell_results$labels
    }
    srt_query$singler_score <- sapply(seq_len(ncol(srt_query)), FUN = function(x) {
      if (isTRUE(prune)) {
        y <- SingleRCell_results$pruned.labels[x]
      } else {
        y <- SingleRCell_results$labels[x]
      }
      if (is.na(y)) {
        out <- NA
      } else {
        out <- SingleRCell_results$scores[x, y]
      }
      return(out)
    })
  }
  return(srt_query)
}
