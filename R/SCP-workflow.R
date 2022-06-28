#' @param srt
#'
#' @param data
#' @param slot
#' @param assay
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
check_DataType <- function(srt, data = NULL, slot = "data", assay = NULL) {
  if (is.null(data)) {
    data <- GetAssayData(srt, slot = slot, assay = assay)
  }
  isfinite <- all(is.finite(range(data, na.rm = TRUE)))
  if ("dgCMatrix" %in% class(data)) {
    isfloat <- any(data@x %% 1 != 0, na.rm = TRUE)
  } else {
    isfloat <- any(data[, sample(seq_len(ncol(data)), min(ncol(data), 1000))] %% 1 != 0, na.rm = TRUE)
  }
  islog <- is.finite(expm1(x = max(data, na.rm = TRUE)))

  if (!isTRUE(isfinite)) {
    warning("Infinite values detected!", immediate. = TRUE)
    return("unknown")
  } else if (!isfloat) {
    return("raw_counts")
  } else if (isfloat && islog) {
    return("log_normalized_counts")
  } else if (isfloat && !islog) {
    return("raw_normalized_counts")
  }
}

#' @param srtList
#'
#' @param batch
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param vars_to_regress
#' @param seed
#' @param ...
#'
#' @importFrom Seurat SplitObject GetAssayData Assays NormalizeData FindVariableFeatures SCTransform SCTResults SelectIntegrationFeatures PrepSCTIntegration DefaultAssay DefaultAssay<- VariableFeatures VariableFeatures<-
#' @importFrom Matrix rowSums
#' @importFrom dplyr "%>%" arrange desc filter .data
#' @importFrom utils head
#' @export
#'
check_srtList <- function(srtList, batch = "orig.ident",
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                          vars_to_regress = NULL, seed = 11, ...) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  set.seed(seed)

  if (class(srtList) != "list" || any(sapply(srtList, class) != "Seurat")) {
    stop("'srtList' is not a list of Seurat object.")
  }
  if (!normalization_method %in% c("logCPM", "SCT")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT'")
  }
  if (normalization_method %in% c("SCT")) {
    check_R("glmGamPoi")
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global','separate'")
  }
  if (any(sapply(srtList, ncol) < 10)) {
    stop(paste0("Seurat objects in srtList contain less than 10 cells. srtList index: ", which(sapply(srtList, ncol) < 10)))
  }

  genelist <- lapply(srtList, function(x) {
    sort(rownames(GetAssayData(x, slot = "counts")))
  })
  if (length(unique(genelist)) != 1) {
    warning("'srtList' have different feature names! Will subset the common features for downstream analysis!", immediate. = TRUE)
    cf <- lapply(srtList, rownames) %>% Reduce(intersect, .)
    for (i in seq_along(srtList)) {
      srtList[[i]] <- subset(srtList[[i]], features = cf)
    }
  }

  celllist <- unlist(lapply(srtList, colnames))
  if (length(celllist) != length(unique(celllist))) {
    stop("'srtList' have duplicated cell names!")
  }

  if (length(batch) != 1 && length(batch) != length(srtList)) {
    stop("'batch' must be a vector to specify the batch column in Seurat object or a vector of the same length of the srtList!")
  }
  if (length(batch) == length(srtList)) {
    srtList_tmp <- list()
    for (bat in unique(batch)) {
      srtList_tmp[[bat]] <- Reduce(merge, srtList[batch == bat])
    }
    srtList <- srtList_tmp
  } else {
    if (!all(sapply(srtList, function(x) {
      batch %in% colnames(x@meta.data)
    }))) {
      stop(paste0("batch column('", batch, "') was not found in one or more object of the srtList!"))
    }
    for (i in seq_along(srtList)) {
      u <- unique(srtList[[i]][[batch, drop = TRUE]])
      if (length(u) > 1) {
        x <- SplitObject(srtList[[i]], split.by = batch)
        srtList[[i]] <- NULL
        srtList <- c(srtList, x)
      }
    }
  }

  status_i <- c()
  for (i in seq_along(srtList)) {
    if (!"RNA" %in% Assays(srtList[[i]])) {
      stop(paste("srtList", i, "does not contain 'RNA' assay."))
    }
    DefaultAssay(srtList[[i]]) <- "RNA"
    if (isTRUE(do_normalization)) {
      cat("Perform NormalizeData(logCPM) on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
      srtList[[i]] <- suppressWarnings(NormalizeData(object = srtList[[i]], normalization.method = "LogNormalize", verbose = FALSE))
    } else if (is.null(do_normalization)) {
      status <- check_DataType(srtList[[i]], slot = "data")
      if (length(status_i) > 0 && status != status_i) {
        stop("The normalization method for the data slot of srt ", i, " is detected as '", status, "', which is inconsistent with the previous.\nYou can manually set the 'do_normalization' parameter to TRUE or FALSE.")
      }
      if (status == "raw_counts") {
        cat("Data ", i, "/", length(srtList), " of the srtList is raw counts. Perform NormalizeData(logCPM) on the data ...\n", sep = "")
        srtList[[i]] <- suppressWarnings(NormalizeData(object = srtList[[i]], normalization.method = "LogNormalize", verbose = FALSE))
      }
      if (status == "log_normalized_counts") {
        cat("Data ", i, "/", length(srtList), " of the srtList has been log-normalized.\n", sep = "")
      }
      if (status == "raw_normalized_counts") {
        cat("Data ", i, "/", length(srtList), " of the srtList is normalized without log transformation. Perform log1p on the data...\n", sep = "")
        srtList[[i]][["RNA"]]@data <- log1p(srtList[[i]][["RNA"]]@data)
      }
      if (status == "unknown") {
        stop("Can not determine whether data ", i, " is log-normalized")
      }
    }
    if (is.null(HVF)) {
      if (isTRUE(do_HVF_finding) || is.null(do_HVF_finding) || length(VariableFeatures(srtList[[i]])) == 0) {
        cat("Perform FindVariableFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- FindVariableFeatures(srtList[[i]], nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE)
      }
    }

    if (normalization_method %in% c("SCT")) {
      if (isTRUE(do_normalization) || isTRUE(do_HVF_finding) || !"SCT" %in% Assays(srtList[[i]])) {
        cat("Perform SCTransform on the data", i, "of the srtList...\n")
        srtList[[i]] <- SCTransform(
          object = srtList[[i]],
          variable.features.n = nHVF,
          vars.to.regress = vars_to_regress,
          assay = "RNA",
          method = "glmGamPoi",
          verbose = FALSE
        )
      } else {
        DefaultAssay(srtList[[i]]) <- "SCT"
      }
      if (!"residual_variance" %in% colnames(srtList[[i]]@assays$SCT@meta.features)) {
        if (length(srtList[[i]]@assays$SCT@SCTModel.list) > 1) {
          index <- which(sapply(srtList[[i]]@assays$SCT@SCTModel.list, function(x) nrow(x@cell.attributes) == ncol(srtList[[i]])))
        } else {
          index <- 1
        }
        model <- srtList[[i]]@assays$SCT@SCTModel.list[[index]]
        feature.attr <- SCTResults(object = model, slot = "feature.attributes")
        nfeatures <- min(nHVF, nrow(x = feature.attr))
        top.features <- rownames(x = feature.attr)[head(order(feature.attr$residual_variance,
          decreasing = TRUE
        ), n = nfeatures)]
        VariableFeatures(object = srtList[[i]]) <- top.features
        srtList[[i]]@assays$SCT@meta.features <- feature.attr
      } else {
        VariableFeatures(srtList[[i]]) <- srtList[[i]]@assays$SCT@meta.features %>%
          arrange(desc(.data[["residual_variance"]])) %>%
          rownames() %>%
          head(n = nHVF)
      }
    }
  }

  if (is.null(HVF)) {
    if (HVF_source == "global") {
      cat("Perform global HVF calculation on the merged datasets from the srtList...\n")
      srtMerge <- Reduce(merge, srtList)
      srtMerge <- FindVariableFeatures(srtMerge, nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE)
      HVF <- VariableFeatures(srtMerge)
    }
    if (HVF_source == "separate") {
      cat("Use the separate HVF from the existed HVF in srtList...\n")
      # HVF_merge <- unlist(lapply(srtList, VariableFeatures))
      # HVF <- names(sort(table(HVF_merge), decreasing = TRUE))[1:nHVF]
      HVF <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
    }
  } else {
    cf <- Reduce(intersect, lapply(srtList, function(x) {
      rownames(GetAssayData(x, slot = "counts"))
    }))
    HVF <- HVF[HVF %in% cf]
  }

  hvf_sum <- lapply(srtList, function(x) {
    colSums(GetAssayData(x, slot = "counts")[HVF, ])
  })
  cell_all <- unlist(unname(hvf_sum))
  cell_abnormal <- names(cell_all)[cell_all == 0]
  if (length(cell_abnormal) > 0) {
    warning("Some cells do not express any of the highly variable features: ", paste(cell_abnormal, collapse = ","), immediate. = TRUE)
  }

  if (normalization_method %in% c("SCT")) {
    srtList <- PrepSCTIntegration(object.list = srtList, anchor.features = HVF, assay = "SCT", verbose = FALSE)
  }
  cat(paste0("[", Sys.time(), "]", " Finished checking.\n"))

  return(list(
    srtList = srtList,
    HVF = HVF
  ))
}

#' @param srtMerge
#'
#' @param batch
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param vars_to_regress
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData SplitObject SetAssayData VariableFeatures VariableFeatures<-
#' @export
check_srtMerge <- function(srtMerge, batch = "orig.ident",
                           do_normalization = NULL, normalization_method = "logCPM",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                           vars_to_regress = NULL, seed = 11, ...) {
  if (class(srtMerge) != "Seurat") {
    stop("'srtMerge' is not a Seurat object.")
  }
  if (length(batch) != 1) {
    stop("'batch' must be a vector to specify the batch column in srtMerge object!")
  }
  if (!batch %in% colnames(srtMerge@meta.data)) {
    stop(paste0("No batch column('", batch, "') found in the srtMerge object!"))
  }
  if (!is.factor(srtMerge[[batch, drop = TRUE]])) {
    srtMerge[[batch, drop = TRUE]] <- factor(srtMerge[[batch, drop = TRUE]],
      levels = unique(srtMerge[[batch, drop = TRUE]])
    )
  }
  srtMerge_raw <- srtMerge

  cat(paste0("[", Sys.time(), "]", " Spliting srtMerge into srtList... ...\n"))
  srtList <- SplitObject(object = srtMerge_raw, split.by = batch)

  checked <- check_srtList(
    srtList = srtList, batch = batch,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srtList <- checked[["srtList"]]
  HVF <- checked[["HVF"]]
  srtMerge <- Reduce(merge, srtList)
  VariableFeatures(srtMerge) <- HVF

  srtMerge <- SrtAppend(
    srt_raw = srtMerge, srt_append = srtMerge_raw, pattern = "",
    slots = "reductions", overwrite = TRUE, verbose = FALSE
  )
  return(list(
    srtMerge = srtMerge,
    srtList = srtList,
    HVF = HVF
  ))
}

#' @param srt
#'
#' @param HVF
#' @param do_normalization
#'
#' @importFrom Seurat DefaultAssay DefaultAssay<- GetAssayData NormalizeData VariableFeatures VariableFeatures<-
check_final <- function(srt, HVF, do_normalization) {
  raw_DefaultAssay <- DefaultAssay(object = srt)
  DefaultAssay(object = srt) <- "RNA"
  if (isTRUE(do_normalization)) {
    cat("Perform NormalizeData(logCPM) on the data [Final check]...\n", sep = "")
    srt <- suppressWarnings(NormalizeData(object = srt, normalization.method = "LogNormalize", verbose = FALSE))
  } else if (is.null(do_normalization)) {
    status <- check_DataType(srt, slot = "data")
    if (status == "raw_counts") {
      cat("Perform NormalizeData(logCPM) on the data [Final check]...\n")
      srt <- suppressWarnings(NormalizeData(object = srt, normalization.method = "LogNormalize", verbose = FALSE))
    }
    if (status == "raw_normalized_counts") {
      cat("Data is normalized without log transformation. Perform log1p on the data [Final check]...\n")
      srt[["RNA"]]@data <- log1p(srt[["RNA"]]@data)
    }
    if (status == "unknown") {
      warning("Unable to determine if data is log-normalized [Final check].", immediate. = TRUE)
    }
  }
  if (length(VariableFeatures(srt)) == 0) {
    HVF <- HVF[HVF %in% rownames(GetAssayData(srt, slot = "counts"))]
    VariableFeatures(srt) <- HVF
  }
  DefaultAssay(object = srt) <- raw_DefaultAssay
  return(srt)
}

#' @param srt
#'
#' @param newnames
#' @param oldnames
#' @param assays
#'
#' @importFrom Seurat Assays GetAssay
#' @importFrom methods slot
#' @export
RenameFeatures <- function(srt, newnames = NULL, oldnames = NULL, assays = NULL) {
  assays <- assays[assays %in% Assays(srt)] %||% Assays(srt)
  if (is.null(oldnames)) {
    oldnames <- rownames(srt)
  }
  if (!identical(length(newnames), length(oldnames))) {
    stop("'newnames' must be the length of features in the srt or the length of oldnames.")
  }
  names(newnames) <- oldnames
  for (assay in assays) {
    message("Rename features for the assay: ", assay)
    assay_obj <- GetAssay(srt, assay)
    for (d in c("meta.features", "scale.data", "counts", "data")) {
      index <- which(rownames(slot(assay_obj, d)) %in% names(newnames))
      rownames(slot(assay_obj, d))[index] <- newnames[rownames(slot(assay_obj, d))[index]]
    }
    if (length(slot(assay_obj, "var.features")) > 0) {
      index <- which(slot(assay_obj, "var.features") %in% names(newnames))
      slot(assay_obj, "var.features")[index] <- newnames[slot(assay_obj, "var.features")[index]]
    }
    srt[[assay]] <- assay_obj
  }
  return(srt)
}

#' @param srt
#'
#' @param features
#' @param reorder_by
#' @param slot
#' @param assay
#' @param log
#' @param distance_metric
#' @param reorder_FUN
#'
#' @importFrom Seurat VariableFeatures DefaultAssay DefaultAssay<- AverageExpression Idents<-
#' @importFrom stats hclust reorder as.dendrogram as.dist
#' @importFrom proxyC dist simil
#' @export
SrtReorder <- function(srt, features = NULL, reorder_by = NULL, slot = "data", assay = NULL, log = TRUE,
                       distance_metric = "euclidean", reorder_FUN = "mean") {
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }
  if (is.null(features)) {
    features <- VariableFeatures(srt)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  if (is.null(reorder_by)) {
    srt$ident <- Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  if (length(unique(srt[[reorder_by, drop = TRUE]])) == 1) {
    warning("Only one cluter found. No need to reorder.", immediate. = TRUE)
    return(srt)
  }
  simil_method <- c(
    "cosine", "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman",
    "simple matching", "faith"
  )
  dist_method <- c(
    "euclidean", "chisquared", "kullback", "manhattan", "maximum", "canberra",
    "minkowski", "hamming"
  )
  if (!distance_metric %in% c(simil_method, dist_method, "pearson", "spearman")) {
    stop(distance_metric, " method is invalid.")
  }

  data.avg <- AverageExpression(object = srt, features = features, slot = slot, assays = assay, group.by = "ident", verbose = FALSE)[[1]]
  if (isTRUE(log)) {
    data.avg <- log1p(data.avg)
  }
  mat <- t(x = data.avg[features, ])
  if (!"matrix" %in% class(mat)) {
    mat <- as.matrix(mat)
  }
  mat <- as(mat, "dgCMatrix")

  if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 - simil(as(mat, "dgCMatrix"), method = distance_metric)
  } else if (distance_metric %in% dist_method) {
    d <- dist(as(mat, "dgCMatrix"), method = distance_metric)
  }
  data.dist <- as.dist(d)
  hc <- hclust(d = data.dist)
  dd <- as.dendrogram(hc)
  dd_ordered <- reorder(dd, wts = colMeans(data.avg[features, , drop = FALSE]), agglo.FUN = reorder_FUN)
  ident_new <- unname(setNames(object = seq_along(labels(dd_ordered)), nm = labels(dd_ordered))[as.character(srt$ident)])
  ident_new <- factor(ident_new, levels = seq_along(labels(dd_ordered)))
  Idents(srt) <- srt$ident <- ident_new
  return(srt)
}

#' @param srt_raw
#'
#' @param srt_append
#' @param slots
#' @param pattern
#' @param overwrite
#' @param verbose
#'
#' @importFrom methods slotNames slot slot<-
#' @export
SrtAppend <- function(srt_raw, srt_append,
                      slots = slotNames(srt_append), pattern = NULL, overwrite = FALSE,
                      verbose = TRUE) {
  if (class(srt_raw) != "Seurat" || class(srt_append) != "Seurat") {
    stop("'srt_raw' or 'srt_append' is not a Seurat object.")
  }
  pattern <- pattern %||% ""
  for (slot_nm in slotNames(srt_append)) {
    if (!slot_nm %in% slots) {
      if (isTRUE(verbose)) {
        message("Slot ", slot_nm, " is not appended.")
      }
      next
    }
    if (identical(slot_nm, "active.ident") && isTRUE(overwrite)) {
      slot(srt_raw, name = "active.ident") <- slot(srt_append, name = "active.ident")
      next
    }
    for (info in names(slot(srt_append, name = slot_nm))) {
      if (is.null(info)) {
        if (length(slot(srt_append, name = slot_nm)) > 0 && isTRUE(overwrite)) {
          slot(srt_raw, name = slot_nm) <- slot(srt_append, name = slot_nm)
        }
        next
      }
      if (!grepl(pattern = pattern, x = info)) {
        if (isTRUE(verbose)) {
          message(info, " in slot ", slot_nm, " is not appended.")
        }
        next
      }
      if (!info %in% names(slot(srt_raw, name = slot_nm)) || isTRUE(overwrite)) {
        if (slot_nm %in% c("assays", "graphs", "neighbors", "reductions", "images")) {
          if (identical(slot_nm, "graphs")) {
            srt_raw@graphs[[info]] <- srt_append[[info]]
          } else {
            srt_raw[[info]] <- srt_append[[info]]
          }
        } else if (identical(slot_nm, "meta.data")) {
          srt_raw@meta.data[, info] <- NULL
          srt_raw@meta.data[[info]] <- srt_append@meta.data[colnames(srt_raw), info]
        } else {
          slot(srt_raw, name = slot_nm)[[info]] <- slot(srt_append, name = slot_nm)[[info]]
        }
      }
    }
  }
  return(srt_raw)
}

#' @param srt
#'
#' @param prefix
#' @param features
#' @param assay
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction
#' @param reduction_use
#' @param dims_use
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param verbose
#' @param seed
#'
#' @importFrom Seurat Embeddings RunPCA RunICA RunTSNE Reductions DefaultAssay DefaultAssay<- Key Key<-
#' @export
RunDimReduction <- function(srt, prefix = NULL, features = NULL, assay = NULL, split.by = NULL,
                            liner_reduction = NULL, liner_reduction_dims = 100, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                            nonliner_reduction = NULL, reduction_use = NULL, dims_use = NULL, graph_use = NULL, nonliner_reduction_dims = 2, nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                            verbose = TRUE, seed = 11) {
  set.seed(seed)
  if (!is.null(split.by)) {
    srt_sp <- SplitObject(srt, split.by = split.by)
  } else {
    srt_sp <- srt
  }
  if (!is.null(liner_reduction)) {
    if (!liner_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca")) {
      stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
    }
    if (is.null(liner_reduction_dims)) {
      liner_reduction_dims <- 100
    }
  }
  if (!is.null(nonliner_reduction)) {
    if (!isTRUE(nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
      stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
    }
    if (is.null(features) && is.null(reduction_use) && is.null(graph_use)) {
      stop("'features', 'reduction_use' or 'graph_use' must be provided when running nonliner reduction.")
    }
    if (is.null(prefix)) {
      prefix <- reduction_use
    }
  }
  srt_out <- lapply(c(srt_sp), function(srt) {
    assay_default <- DefaultAssay(srt)
    DefaultAssay(srt) <- assay %||% DefaultAssay(srt)
    if (!is.null(liner_reduction)) {
      if (!isTRUE(force_liner_reduction)) {
        if (liner_reduction %in% Reductions(srt)) {
          if (srt[[liner_reduction]]@assay.used == assay) {
            message("liner_reduction(", liner_reduction, ") is already existed. Skip calculation.")
            reduc <- srt[[liner_reduction]]
            Key(reduc) <- paste0(prefix, liner_reduction, "_")
            srt[[paste0(prefix, liner_reduction)]] <- reduc
            srt@misc[["Default_reduction"]] <- paste0(prefix, liner_reduction)
            DefaultAssay(srt) <- assay_default
            return(srt)
          } else {
            message("assay.used is ", srt[[liner_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the liner reduction")
          }
        }
      }
      if (is.null(features) || length(features) == 0) {
        message("No features provided. Use variable features.")
        if (length(VariableFeatures(srt)) == 0) {
          srt <- FindVariableFeatures(srt, verbose = FALSE)
        }
        features <- VariableFeatures(srt)
      }
      key_use <- switch(liner_reduction,
        "pca" = "PC_",
        "ica" = "IC_",
        "nmf" = "BE_",
        "mds" = "MDS_",
        "glmpca" = "GLMPC_"
      )
      if (liner_reduction == "pca") {
        srt <- RunPCA(
          object = srt, features = features, npcs = liner_reduction_dims,
          reduction.name = paste0(prefix, liner_reduction),
          reduction.key = paste0(prefix, key_use),
          verbose = verbose, seed.use = seed
        )
        pca.out <- srt[[paste0(prefix, liner_reduction)]]
        center <- rowMeans(GetAssayData(object = srt, slot = "scale.data")[features, ])
        model <- list(sdev = pca.out@stdev, rotation = pca.out@feature.loadings, center = center, scale = FALSE, x = pca.out@cell.embeddings)
        class(model) <- "prcomp"
        srt@reductions[[paste0(prefix, liner_reduction)]]@misc[["model"]] <- model
      }
      if (liner_reduction == "glmpca") {
        check_R(c("R.utils", "zhanghao-njmu/seurat-wrappers", "glmpca"), pkg_names = c("R.utils", "SeuratWrappers", "glmpca"))
        srt <- SeuratWrappers::RunGLMPCA(
          object = srt, features = features, L = liner_reduction_dims,
          reduction.name = paste0(prefix, liner_reduction),
          reduction.key = paste0(prefix, key_use),
          verbose = verbose
        )
      }
      if (liner_reduction == "ica") {
        srt <- RunICA(
          object = srt, features = features, nics = liner_reduction_dims,
          reduction.name = paste0(prefix, liner_reduction),
          reduction.key = paste0(prefix, key_use),
          verbose = verbose, seed.use = seed
        )
      }
      if (liner_reduction == "nmf") {
        srt <- RunNMF(
          object = srt, features = features, nbes = liner_reduction_dims,
          reduction.name = paste0(prefix, liner_reduction),
          reduction.key = paste0(prefix, key_use),
          verbose = verbose, seed.use = seed
        )
      }
      if (liner_reduction == "mds") {
        srt <- RunMDS(
          object = srt, features = features, nmds = liner_reduction_dims,
          dist.method = liner_reduction_distance,
          reduction.name = paste0(prefix, liner_reduction),
          reduction.key = paste0(prefix, key_use),
          verbose = verbose, seed.use = seed
        )
      }
      if (liner_reduction %in% c("glmpca", "nmf")) {
        dims_estimate <- 1:liner_reduction_dims
      }

      dim_est <- tryCatch(expr = {
        intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(prefix, liner_reduction)), k = 20, iterations = 100)[["dim.est"]]
      }, error = function(e) {
        message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
        return(NA)
      })
      if (!is.na(dim_est)) {
        if (ceiling(dim_est) < 10) {
          dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, liner_reduction))), 10))
        } else {
          dims_estimate <- seq_len(ceiling(dim_est))
        }
      } else {
        dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, liner_reduction))), 30))
        message("Set the dims_estimate to ", paste0(range(dims_estimate), collapse = ":"), " for '", liner_reduction, "'")
      }
      srt@reductions[[paste0(prefix, liner_reduction)]]@misc[["dims_estimate"]] <- dims_estimate
      srt@misc[["Default_reduction"]] <- paste0(prefix, liner_reduction)
      DefaultAssay(srt) <- assay_default
    } else if (!is.null(nonliner_reduction)) {
      if (!isTRUE(force_nonliner_reduction)) {
        if (nonliner_reduction %in% Reductions(srt)) {
          if (srt[[nonliner_reduction]]@assay.used == assay) {
            message("nonliner_reduction(", nonliner_reduction, ") is already existed. Skip calculation.")
            reduc <- srt[[nonliner_reduction]]
            Key(reduc) <- paste0(prefix, nonliner_reduction, "_")
            srt[[paste0(prefix, nonliner_reduction)]] <- reduc
            srt@misc[["Default_reduction"]] <- paste0(prefix, nonliner_reduction)
            DefaultAssay(srt) <- assay_default
            return(srt)
          } else {
            message("assay.used is ", srt[[nonliner_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the liner reduction")
          }
        }
      }
      if (nonliner_reduction == "umap") {
        srt <- suppressWarnings(RunUMAP2(
          object = srt, umap.method = "uwot", n.components = nonliner_reduction_dims,
          features = features, reduction = reduction_use, dims = dims_use, graph = graph_use,
          metric = nonliner_reduction_distance,
          reduction.name = paste0(prefix, "UMAP", nonliner_reduction_dims, "D"),
          reduction.key = paste0(prefix, "UMAP", nonliner_reduction_dims, "D_"),
          return.model = TRUE, verbose = verbose, seed.use = seed
        ))
      }
      if (nonliner_reduction == "umap-naive") {
        srt <- suppressWarnings(RunUMAP2(
          object = srt, umap.method = "naive", n.components = nonliner_reduction_dims,
          features = features, reduction = reduction_use, dims = dims_use, graph = graph_use,
          metric = nonliner_reduction_distance,
          reduction.name = paste0(prefix, "UMAP", nonliner_reduction_dims, "D"),
          reduction.key = paste0(prefix, "UMAP", nonliner_reduction_dims, "D_"),
          return.model = TRUE, verbose = verbose, seed.use = seed
        ))
        nonliner_reduction <- "umap"
      }
      if (nonliner_reduction == "umap-learn") {
        srt <- suppressWarnings(RunUMAP2(
          object = srt, umap.method = "umap-learn", n.components = as.integer(nonliner_reduction_dims),
          features = features, reduction = reduction_use, dims = dims_use, graph = graph_use,
          metric = nonliner_reduction_distance,
          reduction.name = paste0(prefix, "UMAP", nonliner_reduction_dims, "D"),
          reduction.key = paste0(prefix, "UMAP", nonliner_reduction_dims, "D_"),
          return.model = FALSE, verbose = verbose, seed.use = seed
        ))
        nonliner_reduction <- "umap"
      }
      if (nonliner_reduction == "tsne") {
        srt <- RunTSNE(
          object = srt, tsne.method = "Rtsne", dim.embed = nonliner_reduction_dims,
          features = features, reduction = reduction_use, dims = dims_use,
          reduction.name = paste0(prefix, "TSNE", nonliner_reduction_dims, "D"),
          reduction.key = paste0(prefix, "TSNE", nonliner_reduction_dims, "D_"),
          num_threads = 0, verbose = verbose, seed.use = seed
        )
      }
      if (nonliner_reduction == "dm") {
        srt <- RunDM(
          object = srt, ndcs = nonliner_reduction_dims,
          features = features, reduction = reduction_use, dims = dims_use,
          dist.method = nonliner_reduction_distance,
          reduction.name = paste0(prefix, "DM", nonliner_reduction_dims, "D"),
          reduction.key = paste0(prefix, "DM", nonliner_reduction_dims, "D_"),
          verbose = verbose, seed.use = seed
        )
      }
      srt@reductions[[paste0(prefix, toupper(nonliner_reduction), nonliner_reduction_dims, "D")]]@misc[["dims_use"]] <- dims_use
      srt@reductions[[paste0(prefix, toupper(nonliner_reduction), nonliner_reduction_dims, "D")]]@misc[["reduction_use"]] <- reduction_use
      srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nonliner_reduction))
      DefaultAssay(srt) <- assay_default
    }
    return(srt)
  })
  if (length(srt_out) == 1) {
    return(srt_out[[1]])
  } else {
    return(srt_out)
  }
}

#' Uncorrected_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_distance
#' @param liner_reduction_dims_use
#' @param force_liner_reduction
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData SetAssayData VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Uncorrected_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                                  do_normalization = NULL, normalization_method = "logCPM",
                                  do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                                  do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                  liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_distance = "euclidean", liner_reduction_dims_use = NULL, force_liner_reduction = FALSE,
                                  nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                                  do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                  seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate", "edge"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate','edge','date','dca','sauice','scVI'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  cat("Perform integration(Uncorrected) on the data...\n")
  srtIntegrated <- Standard_SCP(
    srt = srtMerge, prefix = "Uncorrected",
    do_normalization = do_normalization, normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
    do_scaling = do_scaling, vars_to_regress = vars_to_regress, regression_model = regression_model,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, liner_reduction_distance = liner_reduction_distance, liner_reduction_dims_use = liner_reduction_dims_use, force_liner_reduction = force_liner_reduction,
    nonliner_reduction = nonliner_reduction, nonliner_reduction_dims = nonliner_reduction_dims, nonliner_reduction_distance = nonliner_reduction_distance, force_nonliner_reduction = force_nonliner_reduction,
    do_cluster_finding = do_cluster_finding, cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = cluster_reorder,
    seed = seed
  )
  srtIntegrated[[paste0("Uncorrectedclusters")]] <- srtIntegrated[[paste0("Uncorrected", liner_reduction, "clusters")]]
  srtIntegrated[[paste0("Uncorrected", liner_reduction, "clusters")]] <- NULL
  for (nr in nonliner_reduction) {
    for (n in nonliner_reduction_dims) {
      reduc <- srtIntegrated@reductions[[paste0("Uncorrected", liner_reduction, toupper(nr), n, "D")]]
      Key(reduc) <- paste0("Uncorrected", toupper(nr), n, "D_")
      srtIntegrated@reductions[[paste0("Uncorrected", liner_reduction, toupper(nr), n, "D")]] <- NULL
      srtIntegrated@reductions[[paste0("Uncorrected", toupper(nr), n, "D")]] <- reduc
    }
    srtIntegrated@misc[["Default_reduction"]] <- paste0("Uncorrected", toupper(nr))
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Uncorrected_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "Uncorrected|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Seurat_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData FindIntegrationAnchors IntegrateData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents
#' @importFrom dplyr "%>%"
#' @export
Seurat_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                             do_normalization = NULL, normalization_method = "logCPM",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                             liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                             nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                             do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                             FindIntegrationAnchors_params = list(), IntegrateData_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    message("Using 'rpca' workflow...")
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
        cat("Perform ScaleData on the data ", i, " ...\n", sep = "")
        srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
      }
      cat("Perform linear dimension reduction (pca) on the data ", i, " ...\n", sep = "")
      srt <- RunDimReduction(
        srt = srt, prefix = "", features = HVF, liner_reduction_distance = liner_reduction_distance,
        liner_reduction = "pca", liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
        verbose = FALSE, seed = seed
      )
      srtList[[i]] <- srt
    }
  }

  cat("Perform FindIntegrationAnchors on the data...\n")
  params1 <- list(
    object.list = srtList,
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize",
      "SCT" = "SCT"
    ),
    anchor.features = HVF,
    verbose = FALSE
  )
  for (nm in names(FindIntegrationAnchors_params)) {
    params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
  }
  srt_anchors <- do.call(FindIntegrationAnchors, args = params1)

  cat("Perform integration(Seurat) on the data...\n")
  params2 <- list(
    anchorset = srt_anchors,
    new.assay.name = "Seurat",
    normalization.method = switch(normalization_method,
      "logCPM" = "LogNormalize",
      "SCT" = "SCT"
    ),
    features.to.integrate = HVF,
    verbose = FALSE
  )
  for (nm in names(IntegrateData_params)) {
    params2[[nm]] <- IntegrateData_params[[nm]]
  }
  srtIntegrated <- do.call(IntegrateData, args = params2)

  DefaultAssay(srtIntegrated) <- "Seurat"
  VariableFeatures(srtIntegrated[["Seurat"]]) <- HVF

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data...\n")
    srtIntegrated <- ScaleData(object = srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Seurat", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtIntegrated@reductions[[paste0("Seurat", liner_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Seurat", liner_reduction), dims = liner_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("Seurat", liner_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Seurat", liner_reduction, "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Seurat",
        reduction_use = paste0("Seurat", liner_reduction), dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF


  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "Seurat|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' scVI_integrate
#'
#' @param srtMerge
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param force_nonliner_reduction
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param SCVI_params
#' @param seed
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
scVI_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                           do_normalization = NULL, normalization_method = "logCPM",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                           nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                           do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                           SCVI_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  check_Python("scvi-tools")
  scvi <- import("scvi")
  scipy <- import("scipy")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  adata <- srt_to_adata(srtMerge[HVF, ])
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])
  scvi$model$SCVI$setup_anndata(adata, batch_key = batch)

  params <- list(
    adata = adata
  )
  for (nm in names(SCVI_params)) {
    params[[nm]] <- SCVI_params[[nm]]
  }
  model <- do.call(scvi$model$SCVI, args = params)
  model$train()

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  corrected <- t(as.matrix(model$get_normalized_expression()))
  srtIntegrated[["scVIcorrected"]] <- CreateAssayObject(counts = corrected)
  VariableFeatures(srtIntegrated[["scVIcorrected"]]) <- HVF
  latent <- as.matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("scvi_", seq_len(ncol(latent)))
  srtIntegrated[["scVI"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(srtIntegrated))
  scVI_dims_use <- seq_len(ncol(latent))
  srtIntegrated[["scVI"]]@misc[["dims_estimate"]] <- scVI_dims_use

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "scVI", dims = scVI_dims_use, force.recalc = TRUE, graph.name = paste0("scVI", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("scVI", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["scVIclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "scVI",
        reduction_use = "scVI", dims_use = scVI_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["scVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "scVI|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' MNN_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param fastMNN_dims_use
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
MNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                          nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          mnnCorrect_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  srtIntegrated <- Reduce(merge, srtList)
  VariableFeatures(srtIntegrated) <- HVF

  sceList <- lapply(srtList, function(x) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = GetAssayData(x, slot = "data", assay = DefaultAssay(x))[HVF, ]))
    sce@assays@data$logcounts <- as.matrix(sce@assays@data$logcounts)
    return(sce)
  })

  cat("Perform integration(MNN) on the data...\n")
  params <- list(
    sceList,
    cos.norm.out = FALSE
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- do.call(batchelor::mnnCorrect, args = params)

  srtIntegrated[["MNNcorrected"]] <- CreateAssayObject(counts = out@assays@data$corrected)
  VariableFeatures(srtIntegrated[["MNNcorrected"]]) <- HVF
  DefaultAssay(srtIntegrated) <- "MNNcorrected"

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data...\n")
    srtIntegrated <- ScaleData(object = srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "MNN", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtIntegrated@reductions[[paste0("MNN", liner_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("MNN", liner_reduction), dims = liner_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("MNN", liner_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("MNN", liner_reduction, "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["MNNclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "MNN",
        reduction_use = paste0("MNN", liner_reduction), dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["MNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "MNN|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' fastMNN_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param fastMNN_dims_use
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
fastMNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                              fastMNN_dims_use = NULL,
                              nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              fastMNN_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  srtIntegrated <- Reduce(merge, srtList)
  VariableFeatures(srtIntegrated) <- HVF

  sceList <- lapply(srtList, function(x) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = GetAssayData(x, slot = "data", assay = DefaultAssay(x))[HVF, ]))
    sce@assays@data$logcounts <- as.matrix(sce@assays@data$logcounts)
    return(sce)
  })

  cat("Perform integration(fastMNN) on the data...\n")
  params <- list(
    sceList
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- do.call(batchelor::fastMNN, args = params)

  srtIntegrated[["fastMNNcorrected"]] <- CreateAssayObject(counts = as.matrix(out@assays@data$reconstructed))
  VariableFeatures(srtIntegrated[["fastMNNcorrected"]]) <- HVF
  reduction <- out@int_colData$reducedDims$corrected
  colnames(reduction) <- paste0("fastMNN_", seq_len(ncol(reduction)))
  srtIntegrated[["fastMNN"]] <- CreateDimReducObject(embeddings = reduction, key = "fastMNN_", assay = DefaultAssay(srtIntegrated))

  if (is.null(fastMNN_dims_use)) {
    dim_est <- tryCatch(expr = {
      intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "fastMNN"), k = 20, iterations = 100)[["dim.est"]]
    }, error = function(e) {
      message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
      return(NA)
    })
    if (!is.na(dim_est)) {
      if (ceiling(dim_est) < 10) {
        fastMNN_dims_use <- seq_len(min(ncol(Embeddings(srtIntegrated, reduction = "fastMNN")), 10))
      } else {
        fastMNN_dims_use <- seq_len(ceiling(dim_est))
      }
    } else {
      fastMNN_dims_use <- seq_len(min(ncol(Embeddings(srtIntegrated, reduction = "fastMNN")), 30))
      message("Set the dims_estimate to ", fastMNN_dims_use, " for 'fastMNN'")
    }
  }
  srtIntegrated@reductions[["fastMNN"]]@misc[["dims_estimate"]] <- fastMNN_dims_use

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "fastMNN", dims = fastMNN_dims_use, force.recalc = TRUE, graph.name = paste0("fastMNN", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("fastMNN", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["fastMNNclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "fastMNN",
        reduction_use = "fastMNN", dims_use = fastMNN_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["fastMNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "fastMNN|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Harmony_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param Harmony_dims_use
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Harmony_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                              liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                              Harmony_dims_use = NULL,
                              nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              RunHarmony_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_R("harmony")

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "Harmony", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("Harmony", liner_reduction)]]@misc[["dims_estimate"]]
  }

  cat("Perform integration(Harmony) on the data...\n")
  params <- list(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = DefaultAssay(srtMerge),
    reduction = paste0("Harmony", liner_reduction),
    dims.use = liner_reduction_dims_use,
    reduction.save = "Harmony",
    max.iter.harmony = 100,
    verbose = FALSE
  )
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srtIntegrated <- do.call(RunHarmony2, args = params)

  if (is.null(Harmony_dims_use)) {
    dim_est <- tryCatch(expr = {
      intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srtIntegrated, reduction = "Harmony"), k = 20, iterations = 100)[["dim.est"]]
    }, error = function(e) {
      message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
      return(NA)
    })
    if (!is.na(dim_est)) {
      if (ceiling(dim_est) < 10) {
        Harmony_dims_use <- seq_len(min(ncol(Embeddings(srtIntegrated, reduction = "Harmony")), 10))
      } else {
        Harmony_dims_use <- seq_len(ceiling(dim_est))
      }
    } else {
      Harmony_dims_use <- seq_len(min(ncol(Embeddings(srtIntegrated, reduction = "Harmony")), 30))
      message("Set the dims_estimate to ", Harmony_dims_use, " for 'Harmony'")
    }
  }
  srtIntegrated[["Harmony"]]@misc[["dims_estimate"]] <- Harmony_dims_use

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Harmony", dims = Harmony_dims_use, force.recalc = TRUE, graph.name = paste0("Harmony", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Harmony", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["Harmonyclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Harmony",
        reduction_use = "Harmony", dims_use = Harmony_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Harmony_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "Harmony|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Scanorama_integrate
#'
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @importFrom Matrix t
#' @importFrom dplyr "%>%"
#' @importFrom reticulate import
#' @importFrom stats sd
#' @export
Scanorama_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                                do_normalization = NULL, normalization_method = "logCPM",
                                do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                                do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                                nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                                do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                scanorama_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_Python("scanorama", envname = "SCP")
  scanorama <- import("scanorama")

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  cat("Perform integration(Scanorama) on the data...\n")
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srtList)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = srtList[[i]], slot = "data", assay = DefaultAssay(srtList[[1]]))))
    genelist[[i]] <- rownames(srtList[[i]])
  }
  params <- list(
    datasets_full = assaylist,
    genes_list = genelist,
    return_dimred = TRUE,
    return_dense = TRUE,
    union = FALSE
  )
  for (nm in names(scanorama_params)) {
    params[[nm]] <- scanorama_params[[nm]]
  }
  integrated.corrected.data <- do.call(scanorama$correct, args = params)

  cor_value <- t(do.call(rbind, integrated.corrected.data[[2]]))
  rownames(cor_value) <- integrated.corrected.data[[3]]
  colnames(cor_value) <- unlist(sapply(assaylist, rownames))
  dim_reduction <- do.call(rbind, integrated.corrected.data[[1]])
  rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
  colnames(dim_reduction) <- paste0("scanorama_", seq_len(ncol(dim_reduction)))
  stdevs <- apply(dim_reduction, MARGIN = 2, FUN = sd)

  srtIntegrated <- Reduce(merge, srtList)
  srtIntegrated[["Scanoramacorrected"]] <- CreateAssayObject(data = cor_value)
  srtIntegrated[["Scanorama"]] <- CreateDimReducObject(embeddings = dim_reduction, stdev = stdevs, key = "scanorama_", assay = DefaultAssay(srtIntegrated))
  VariableFeatures(srtIntegrated[["Scanoramacorrected"]]) <- HVF
  DefaultAssay(srtIntegrated) <- "Scanoramacorrected"

  cat("Perform ScaleData on the data...\n")
  srtIntegrated <- ScaleData(srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Scanorama", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtIntegrated@reductions[[paste0("Scanorama", liner_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Scanorama", liner_reduction), dims = liner_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("Scanorama", liner_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("Scanorama", liner_reduction, "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["Scanoramaclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "Scanorama",
        reduction_use = paste0("Scanorama", liner_reduction), dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Scanorama_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "Scanorama|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' BBKNN_integrate
#'
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- as.Graph Embeddings FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom Matrix t
#' @importFrom dplyr "%>%"
#' @importFrom reticulate import
#' @export
BBKNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                            nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            bbknn_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_Python("bbknn", envname = "SCP")
  bbknn <- import("bbknn")

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "BBKNN", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("BBKNN", liner_reduction)]]@misc[["dims_estimate"]]
  }

  cat("Perform integration(BBKNN) on the data...\n")
  pca <- Embeddings(srtMerge, reduction = paste0("BBKNN", liner_reduction))[, liner_reduction_dims_use]
  params <- list(
    pca = pca,
    batch_list = srtMerge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- do.call(bbknn$bbknn_matrix, args = params)

  bbknn_graph <- as(as(bem[[2]], "CsparseMatrix"), "dgCMatrix")
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- colnames(srtMerge)
  bbknn_graph <- as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- DefaultAssay(srtMerge)
  srtMerge@graphs[["BBKNN"]] <- bbknn_graph
  srtIntegrated <- srtMerge

  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "BBKNN", resolution = cluster_resolution, algorithm = cluster_algorithm_index, verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["BBKNNclusters"]] <- Idents(srtIntegrated)
  }

  cat("Perform nonlinear dimension reduction (umap) on the data...\n")
  message("BBKNN only support umap-learn method.")
  check_Python("umap-learn", envname = "SCP")
  for (n in nonliner_reduction_dims) {
    srtIntegrated <- RunDimReduction(
      srt = srtIntegrated, prefix = "BBKNN",
      graph_use = "BBKNN",
      nonliner_reduction = "umap-learn", nonliner_reduction_dims = n,
      nonliner_reduction_distance = nonliner_reduction_distance,
      force_nonliner_reduction = force_nonliner_reduction,
      verbose = FALSE, seed = seed
    )
  }
  srtIntegrated@misc[["Default_reduction"]] <- "BBKNNUMAP"
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["BBKNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "BBKNN|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' CSS_integrate
#'
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param CSS_param
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
CSS_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                          nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          CSS_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(liner_reduction) > 1) {
    warning("Only the first method in the 'liner_reduction' will be used.", immediate. = TRUE)
    liner_reduction <- liner_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_R("quadbiolab/simspec")
  suppressPackageStartupMessages(require("qlcMatrix"))

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data...\n")
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform linear dimension reduction (", liner_reduction, ") on the data...\n", sep = "")
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "CSS", features = HVF, liner_reduction_distance = liner_reduction_distance,
    liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- srtMerge@reductions[[paste0("CSS", liner_reduction)]]@misc[["dims_estimate"]]
  }

  cat("Perform integration(CSS) on the data...\n")
  params <- list(
    object = srtMerge,
    use_dr = paste0("CSS", liner_reduction),
    dims_use = liner_reduction_dims_use,
    var_genes = HVF,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  for (nm in names(CSS_params)) {
    params[[nm]] <- CSS_params[[nm]]
  }
  srtIntegrated <- do.call(simspec::cluster_sim_spectrum, args = params)

  CSS_dims_use <- seq_len(ncol(Embeddings(srtIntegrated, reduction = "CSS")))
  srtIntegrated@reductions[["CSS"]]@misc[["dims_estimate"]] <- CSS_dims_use
  if (any(is.na(srtIntegrated@reductions[["CSS"]]@cell.embeddings))) {
    stop("NA detected in the CSS embeddings. You can try to use a lower resolution value in the CSS_param.")
  }

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "CSS", dims = CSS_dims_use, force.recalc = TRUE, graph.name = paste0("CSS", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("CSS", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["CSSclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "CSS",
        reduction_use = "CSS", dims_use = CSS_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["CSS_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "CSS|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' LIGER_integrate
#' @param srtMerge
#'
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
LIGER_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            RunOptimizeALS_params = list(), RunQuantileNorm_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_R(c("R.utils", "zhanghao-njmu/seurat-wrappers", "rliger"), pkg_names = c("R.utils", "SeuratWrappers", "rliger"))
  suppressPackageStartupMessages(library("SeuratWrappers"))

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat("Perform ScaleData on the data(do not center for LIGER)...\n")
    srtMerge_prep <- ScaleData(object = srtMerge, features = HVF, split.by = batch, do.center = FALSE, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat("Perform integration(LIGER) on the data...\n")
  params1 <- list(
    object = srtMerge_prep,
    k = 20,
    lambda = 5,
    split.by = batch,
    verbose = FALSE
  )
  for (nm in names(RunOptimizeALS_params)) {
    params1[[nm]] <- RunOptimizeALS_params[[nm]]
  }
  srtMerge_prep <- do.call(RunOptimizeALS, args = params1)

  params2 <- list(
    object = srtMerge_prep,
    reduction.name = "LIGER",
    reduction.key = "LIGER_",
    split.by = batch,
    verbose = FALSE
  )
  for (nm in names(RunQuantileNorm_params)) {
    params2[[nm]] <- RunQuantileNorm_params[[nm]]
  }
  srtIntegrated <- do.call(RunQuantileNorm, args = params2)

  srtMerge_prep <- NULL
  LIGER_dims_use <- seq_len(ncol(Embeddings(srtIntegrated, reduction = "LIGER")))
  srtIntegrated@reductions[["LIGER"]]@misc[["dims_estimate"]] <- LIGER_dims_use

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "LIGER", dims = LIGER_dims_use, force.recalc = TRUE, graph.name = paste0("LIGER", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("LIGER", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["LIGERclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "LIGER",
        reduction_use = "LIGER", dims_use = LIGER_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["LIGER_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "LIGER|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Conos_integrate
#'
#' @param srtMerge
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param liner_reduction_dims
#' @param liner_reduction_dims_use
#' @param liner_reduction_distance
#' @param force_liner_reduction
#' @param nonliner_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @return
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Conos_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            liner_reduction_dims = 100, liner_reduction_dims_use = NULL, liner_reduction_distance = "euclidean", force_liner_reduction = FALSE,
                            nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            buildGraph_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("conos")
  require("conos", quietly = TRUE)
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  srtIntegrated <- Reduce(merge, srtList)
  VariableFeatures(srtIntegrated) <- HVF

  liner_reduction <- "pca"
  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
      cat("Perform ScaleData on the data ", i, " ...\n", sep = "")
      srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    cat("Perform linear dimension reduction (", liner_reduction, ") on the data ", i, " ...\n", sep = "")
    srt <- RunDimReduction(
      srt = srt, prefix = "", features = HVF, liner_reduction_distance = liner_reduction_distance,
      liner_reduction = liner_reduction, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
      verbose = FALSE, seed = seed
    )
    srtList[[i]] <- srt
  }

  if (is.null(liner_reduction_dims_use)) {
    liner_reduction_dims_use <- 1:max(unlist(lapply(srtList, function(srt) srt@reductions[[liner_reduction]]@misc[["dims_estimate"]])))
  }

  cat("Perform integration(Conos) on the data...\n")
  srtList_con <- Conos$new(srtList)
  params <- list(
    space = "PCA",
    ncomps = max(liner_reduction_dims_use)
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  do.call(srtList_con[["buildGraph"]], args = params)

  conos_graph <- as.matrix(srtList_con$graph)
  conos_graph <- as.Graph(conos_graph)
  conos_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["Conos"]] <- conos_graph

  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "Conos", resolution = cluster_resolution, algorithm = cluster_algorithm_index, verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["Conosclusters"]] <- Idents(srtIntegrated)
  }

  cat("Perform nonlinear dimension reduction (umap) on the data...\n")
  message("Conos only support umap-learn method.")
  check_Python("umap-learn", envname = "SCP")
  for (n in nonliner_reduction_dims) {
    srtIntegrated <- RunDimReduction(
      srt = srtIntegrated, prefix = "Conos",
      graph_use = "Conos",
      nonliner_reduction = "umap-learn", nonliner_reduction_dims = n,
      nonliner_reduction_distance = nonliner_reduction_distance,
      force_nonliner_reduction = force_nonliner_reduction,
      verbose = FALSE, seed = seed
    )
  }
  srtIntegrated@misc[["Default_reduction"]] <- "ConosUMAP"
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Conos_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "Conos|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' ZINBWaVE_integrate
#'
#' @param srtMerge
#' @param batch
#' @param append
#' @param srtList
#' @param do_normalization
#' @param normalization_method
#' @param do_HVF_finding
#' @param HVF_source
#' @param HVF_method
#' @param nHVF
#' @param HVF
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param nonliner_reduction
#' @param nonliner_reduction_dims
#' @param nonliner_reduction_distance
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param seed
#' @param ...
#'
#' @return
#' @importFrom Seurat CreateSeuratObject as.SingleCellExperiment GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
ZINBWaVE_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                               do_normalization = NULL, normalization_method = "logCPM",
                               do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                               do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                               nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine",
                               do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                               zinbwave_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  check_R("zinbwave")

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- lapply(srtList, colnames) %>%
      unlist() %>%
      unique() %>%
      sort()
    cell2 <- colnames(srtMerge) %>%
      unique() %>%
      sort()
    if (!identical(cell1, cell2)) {
      stop("srtList and srtMerge have different cells.")
    }
  }
  if (!is.null(srtMerge)) {
    srtMerge_raw <- srtMerge
  } else {
    srtMerge_raw <- NULL
  }
  if (!is.null(srtList)) {
    checked <- check_srtList(
      srtList = srtList, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding, do_scaling = do_scaling,
      normalization_method = normalization_method,
      HVF_source = HVF_source, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  status <- check_DataType(srtMerge, slot = "counts")
  if (status != "raw_counts") {
    stop("All Seurat objects must have the raw counts!")
  }

  sce <- as.SingleCellExperiment(CreateSeuratObject(
    counts = GetAssayData(srtMerge, slot = "counts", assay = DefaultAssay(srtMerge))[HVF, ],
    meta.data = srtMerge@meta.data[, batch, drop = FALSE]
  ))
  sce@assays@data$counts <- as.matrix(sce@assays@data$counts)

  cat("Perform integration(ZINBWaVE) on the data...\n")
  params <- list(
    Y = sce,
    X = paste0("~", batch),
    which_genes = HVF,
    K = 50
  )
  if (ncol(sce) > 10000) {
    params <- c(params, list(prop_fit = 10000 / ncol(sce)))
  }
  for (nm in names(zinbwave_params)) {
    params[[nm]] <- zinbwave_params[[nm]]
  }
  if (ncol(sce) <= 10000) {
    sce_zinbwave <- do.call(zinbwave::zinbwave, args = params)
  } else {
    state <- 1
    while (state != 0) {
      tryCatch(expr = {
        sce_zinbwave <- do.call(zinbwave::zinbsurf, args = params)
        state <- 0
      }, error = function(error) {
        message(error)
        cat("Resampling the cells for zinbsurf ......\n")
        state <- state + 1
        if (state >= 3) {
          stop("Integration failed.\n")
        }
      })
    }
  }

  reduction <- sce_zinbwave@int_colData$reducedDims$zinbwave
  colnames(reduction) <- paste0("ZINBWaVE_", seq_len(ncol(reduction)))
  srtMerge[["ZINBWaVE"]] <- CreateDimReducObject(embeddings = reduction, key = "ZINBWaVE_", assay = DefaultAssay(srtMerge))
  srtIntegrated <- srtMerge
  srtMerge <- sce_zinbwave <- sce <- NULL

  ZINBWaVE_dims_use <- seq_len(ncol(Embeddings(srtIntegrated, reduction = "ZINBWaVE")))
  srtIntegrated@reductions[["ZINBWaVE"]]@misc[["dims_estimate"]] <- ZINBWaVE_dims_use

  srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "ZINBWaVE", dims = ZINBWaVE_dims_use, force.recalc = TRUE, graph.name = paste0("ZINBWaVE", "_", c("KNN", "SNN")), verbose = FALSE)
  if (isTRUE(do_cluster_finding)) {
    cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
    srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0("ZINBWaVE", "_SNN"), verbose = FALSE)
    if (isTRUE(cluster_reorder)) {
      cat("Reorder clusters...\n", sep = "")
      srtIntegrated <- SrtReorder(srtIntegrated, slot = "data", features = HVF, reorder_by = "seurat_clusters")
    }
    srtIntegrated[["seurat_clusters"]] <- NULL
    srtIntegrated[["ZINBWaVEclusters"]] <- Idents(srtIntegrated)
  }

  for (nr in nonliner_reduction) {
    cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
    for (n in nonliner_reduction_dims) {
      srtIntegrated <- RunDimReduction(
        srt = srtIntegrated, prefix = "ZINBWaVE",
        reduction_use = "ZINBWaVE", dims_use = liner_reduction_dims_use,
        nonliner_reduction = nr, nonliner_reduction_dims = n,
        nonliner_reduction_distance = nonliner_reduction_distance,
        force_nonliner_reduction = force_nonliner_reduction,
        verbose = FALSE, seed = seed
      )
    }
  }
  DefaultAssay(srtIntegrated) <- "RNA"
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["ZINBWaVE_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = "ZINBWaVE|RNA|Default_reduction", overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}


#' Standard_SCP
#'
#' Single cell pipeline for the single dataset.
#'
#' @inheritParams Integration_SCP
#' @param srt A \code{Seurat} object.
#' @param prefix The prefix used to name the result.
#'
#' @return A \code{Seurat} object containing the result.
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- Standard_SCP(srt = pancreas1k)
#' ClassDimPlot(pancreas1k, group.by = "CellType")
#' @importFrom Seurat Assays GetAssayData NormalizeData SCTransform SCTResults ScaleData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%" filter arrange desc
#' @importFrom Matrix rowSums
#' @export
Standard_SCP <- function(srt, prefix = "Standard",
                         do_normalization = NULL, normalization_method = "logCPM",
                         do_HVF_finding = TRUE, HVF_method = "vst", nHVF = 2000, HVF = NULL,
                         do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                         liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_distance = "euclidean", liner_reduction_dims_use = NULL, force_liner_reduction = FALSE,
                         nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                         do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                         seed = 11) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  reduc_test <- c(reduc_test, Reductions(srt))
  if (!liner_reduction %in% reduc_test) {
    stop("'liner_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonliner_reduction %in% c("umap", "umap-naive", "umap-learn", "tsne", "dm", "fdg", "isomap", "dbmap", "phate"))) {
    stop("'nonliner_reduction' must be one of 'umap', 'tsne', 'dm', 'fdg', 'isomap', 'dbmap', 'phate'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "slm" = 3,
    "leiden" = 4
  )

  time_start <- Sys.time()
  cat(paste0("[", time_start, "] ", "Start Standard_SCP\n"))
  set.seed(seed)

  checked <- check_srtList(
    srtList = list(srt), batch = "",
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = "separate", nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srt <- checked$srtList[[1]]
  HVF <- checked$HVF

  liner_reduction_dims <- min(liner_reduction_dims, ncol(srt) - 1)
  if (!is.null(liner_reduction_dims_use) && max(liner_reduction_dims_use) > liner_reduction_dims) {
    liner_reduction_dims <- max(liner_reduction_dims_use)
  }
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
    if (normalization_method != "SCT") {
      cat("Perform ScaleData on the data...\n")
      srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
  }

  for (lr in liner_reduction) {
    cat("Perform linear dimension reduction (", lr, ") on the data...\n", sep = "")
    srt <- RunDimReduction(
      srt = srt, prefix = prefix, features = HVF, liner_reduction_distance = liner_reduction_distance,
      liner_reduction = lr, liner_reduction_dims = liner_reduction_dims, force_liner_reduction = force_liner_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(liner_reduction_dims_use)) {
      liner_reduction_dims_use <- srt@reductions[[paste0(prefix, lr)]]@misc[["dims_estimate"]]
    }

    srt <- FindNeighbors(object = srt, reduction = paste0(prefix, lr), dims = liner_reduction_dims_use, force.recalc = TRUE, graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")), verbose = FALSE)
    if (isTRUE(do_cluster_finding)) {
      cat("Perform FindClusters (", cluster_algorithm, ") on the data...\n", sep = "")
      srt <- FindClusters(object = srt, resolution = cluster_resolution, algorithm = cluster_algorithm_index, graph.name = paste0(prefix, lr, "_SNN"), verbose = FALSE)
      if (isTRUE(cluster_reorder)) {
        cat("Reorder clusters...\n", sep = "")
        srt <- SrtReorder(srt, features = HVF, reorder_by = "seurat_clusters", slot = "data")
      }
      srt[["seurat_clusters"]] <- NULL
      srt[[paste0(prefix, lr, "clusters")]] <- Idents(srt)
    }

    for (nr in nonliner_reduction) {
      cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
      for (n in nonliner_reduction_dims) {
        srt <- RunDimReduction(
          srt = srt, prefix = NULL,
          reduction_use = paste0(prefix, lr), dims_use = liner_reduction_dims_use,
          nonliner_reduction = nr, nonliner_reduction_dims = n,
          nonliner_reduction_distance = nonliner_reduction_distance,
          verbose = FALSE, seed = seed
        )
      }
    }
  }
  VariableFeatures(srt) <- srt@misc[["Standard_HVF"]] <- HVF

  srt[[paste0(prefix, "clusters")]] <- srt[[paste0(prefix, liner_reduction[1], "clusters")]]
  for (nr in nonliner_reduction) {
    for (n in nonliner_reduction_dims) {
      reduc <- srt@reductions[[paste0(prefix, liner_reduction[1], toupper(nr), n, "D")]]
      srt@reductions[[paste0(prefix, toupper(nr), n, "D")]] <- reduc
    }
    srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nr))
  }

  DefaultAssay(srt) <- "RNA"

  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "Standard_SCP done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

  return(srt)
}

#' Integration_SCP
#'
#' Single cell pipeline for the integration of multiple datasets.
#'
#' @param srtList A list of \code{Seurat} object.
#' @param srtMerge A merged \code{Seurat} object with batch information.
#' @param batch Metadata column name containing the batch information.
#' @param append Whether append results into the \code{srtMerge}. Only valid when srtMerge is provided.
#' @param integration_method Integration method. Can be one of "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ZINBWaVE".
#' @param do_normalization Whether to normalize the data. If NULL, will automatically determine.
#' @param normalization_method Normalization method.Can be one of "logCPM", "SCT".
#' @param do_HVF_finding Whether to find the high variable features(HVF). If NULL, will automatically determine.
#' @param HVF_source Source of the HVF. Can be one of "separate" and "global".
#' @param nHVF HVF number to use.
#' @param HVF Custom high variable features.
#' @param do_scaling Whether to scale the data. If NULL, will automatically determine.
#' @param vars_to_regress Variables to regress out.
#' @param regression_model Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param liner_reduction Liner reduction method name. Can be one of "pca", "ica", "nmf", "mds", "glmpca".
#' @param liner_reduction_dims Dimensions to calculate when performing liner reduction.
#' @param liner_reduction_distance A distance method if liner reduction need.
#' @param liner_reduction_dims_use Which dimensions to use.
#' @param nonliner_reduction Non-liner reduction method name. Can be one of "umap", "umap-naive", "umap-learn", "tsne", "dm".
#' @param nonliner_reduction_dims Dimensions to calculate when performing non-liner reduction.
#' @param nonliner_reduction_distance A distance method if non-liner reduction need.
#' @param cluster_algorithm Algorithm for modularity optimization when finding clusters.
#' @param cluster_resolution Cluster resolution parameter.
#' @param cluster_reorder Whether to reorder the cluster names using hierarchical clustering.
#' @param seed Set a random seed.
#' @param HVF_method
#' @param force_liner_reduction
#' @param do_cluster_finding
#' @param ...
#'
#' @return A \code{Seurat} object containing the result.
#'
#' @examples
#' if (!require("SeuratData", quietly = TRUE)) {
#'   devtools::install_github("zhanghao-njmu/seurat-data")
#' }
#' library(SeuratData)
#' library(cowplot)
#' suppressWarnings(InstallData("panc8"))
#' data("panc8")
#' cell_sub <- unlist(lapply(split(colnames(panc8), panc8$tech), function(x) sample(x, size = 500)))
#' panc8 <- subset(panc8, cells = cell_sub)
#' plist <- list()
#' for (method in c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ZINBWaVE")) {
#'   panc8 <- Integration_SCP(panc8, batch = "tech", integration_method = method)
#'   plist[[method]] <- ClassDimPlot(panc8, group.by = c("tech", "celltype"), theme_use = "theme_blank")
#'   print(plist[[method]])
#' }
#' @export
Integration_SCP <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            integration_method = "Uncorrected",
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            liner_reduction = "pca", liner_reduction_dims = 100, liner_reduction_distance = "euclidean", liner_reduction_dims_use = NULL, force_liner_reduction = FALSE,
                            nonliner_reduction = "umap", nonliner_reduction_dims = c(2, 3), nonliner_reduction_distance = "cosine", force_nonliner_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            seed = 11, ...) {
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("Must be provided with one of the 'srtList' and 'srtMerge'")
  }
  if (length(integration_method) == 1 && integration_method %in% c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ZINBWaVE")) {
    args1 <- mget(names(formals()))
    args2 <- as.list(match.call())
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("integration_method", "...")]

    time_start <- Sys.time()
    cat(paste0("[", time_start, "] ", paste0("Start ", integration_method, "_integrate"), "\n"))
    srtIntegrated <- do.call(
      what = paste0(integration_method, "_integrate"),
      args = args1[names(args1) %in% formalArgs(paste0(integration_method, "_integrate"))]
    )
    time_end <- Sys.time()
    cat(paste0("[", time_end, "] ", paste0(integration_method, "_integrate done\n")))
    cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

    return(srtIntegrated)
  } else {
    stop(paste(integration_method, "is not a suppoted integration method!"))
  }
}
