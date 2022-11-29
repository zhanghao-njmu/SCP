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
  if (inherits(data, "dgCMatrix")) {
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
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                          vars_to_regress = NULL, seed = 11, ...) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  set.seed(seed)

  if (!inherits(srtList, "list") || any(sapply(srtList, function(x) !inherits(x, "Seurat")))) {
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
      if (isTRUE(HVF_intersect)) {
        HVF_sort <- sort(table(unlist(lapply(srtList, VariableFeatures))), decreasing = TRUE)
        HVF_sort <- HVF_sort[HVF_sort >= HVF_min_intersection]
        HVF <- names(head(HVF_sort, nHVF))
      } else {
        HVF <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
      }
      if (length(HVF) == 0) {
        stop("No HVF available.")
      }
    }
  } else {
    cf <- Reduce(intersect, lapply(srtList, function(x) {
      rownames(GetAssayData(x, slot = "counts"))
    }))
    HVF <- HVF[HVF %in% cf]
  }
  message("Number of available HVF: ", length(HVF))

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
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                           vars_to_regress = NULL, seed = 11, ...) {
  if (!inherits(srtMerge, "Seurat")) {
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
    HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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

#' Attempt to recover raw counts from the normalized matrix.
#'
#' @param srt A Seurat object.
#' @param assay Name of assay to recover counts.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#'  However, due to decimal point preservation during normalization, the calculated nCount is usually a floating point number close to the integer.
#'  The tolerance is its difference from the integer. Default is 0.1
#' @param verbose
#' @examples
#' data("pancreas_sub")
#' raw_counts <- pancreas_sub@assays$RNA@counts
#'
#' # Normalized the data
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#'
#' # Now replace counts with the log-normalized data matrix
#' pancreas_sub@assays$RNA@counts <- pancreas_sub@assays$RNA@data
#'
#' # Recover the counts and compare with the raw counts matrix
#' pancreas_sub <- RecoverCounts(pancreas_sub)
#' identical(raw_counts, pancreas_sub@assays$RNA@counts)
#'
#' @importFrom Seurat GetAssayData SetAssayData
#' @importFrom SeuratObject as.sparse
#' @export
RecoverCounts <- function(srt, assay = "RNA", min_count = c(1, 2, 3), tolerance = 0.1, verbose = TRUE) {
  counts <- GetAssayData(srt, assay = assay, slot = "counts")
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as.sparse(counts[1:nrow(counts), ])
  }
  status <- check_DataType(data = counts)
  if (status == "raw_counts") {
    if (isTRUE(verbose)) {
      message("The data is already raw counts.")
    }
    return(srt)
  }
  if (status == "log_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be log-normalized.")
    }
    counts <- expm1(counts)
  }
  if (status == "raw_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be normalized without log transformation.")
    }
  }
  sf <- unique(round(colSums(counts)))
  if (isTRUE(verbose)) {
    message("The presumed scale factor: ", paste0(head(sf, 10), collapse = ", "))
  }
  if (length(sf) == 1) {
    counts <- counts / sf
    elements <- split(counts@x, rep(1:ncol(counts), diff(counts@p)))
    min_norm <- sapply(elements, min)
    nCount <- NULL
    for (m in min_count) {
      if (is.null(nCount)) {
        presumed_nCount <- m / min_norm
        diff_value <- abs(presumed_nCount - round(presumed_nCount))
        if (max(diff_value) < tolerance) {
          nCount <- round(presumed_nCount)
        }
      }
    }
    if (is.null(nCount)) {
      warning("The presumed nCount of some cells is not valid: ", paste0(head(colnames(counts)[diff_value < tolerance], 10), collapse = ","), ", ...", immediate. = TRUE)
      return(srt)
    }
    counts@x <- round(counts@x * rep(nCount, diff(counts@p)))
    srt <- SetAssayData(srt, new.data = counts, assay = assay, slot = "counts")
    srt[[paste0("nCount_", assay)]] <- nCount
  } else {
    warning("Scale factor is not unique. No changes to be made.", immediate. = TRUE)
  }
  return(srt)
}

#' Rename features for the Seurat object
#'
#' @param srt
#'
#' @param newnames
#' @param assays
#'
#' @examples
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(stringr::str_to_title(rownames(panc8_sub)))
#' panc8_rename <- RenameFeatures(panc8_sub, newnames = genenames)
#' head(rownames(panc8_rename))
#'
#' @importFrom Seurat Assays GetAssay
#' @importFrom methods slot
#' @export
RenameFeatures <- function(srt, newnames = NULL, assays = NULL) {
  assays <- assays[assays %in% Assays(srt)] %||% Assays(srt)
  if (is.null(names(newnames))) {
    if (!identical(length(newnames), nrow(srt))) {
      stop("'newnames' must be named or the length of features in the srt.")
    }
    if (length(unique(sapply(pancreas_sub@assays, nrow))) > 1) {
      stop("Assays in the srt object have different number of features. Please use a named vectors.")
    }
    names(newnames) <- rownames(srt)
  }
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

#' Rename clusters for the Seurat object
#'
#' @param srt
#'
#' @param newnames
#' @param assays
#' @examples
#' data("pancreas_sub")
#' levels(pancreas_sub@meta.data[["SubCellType"]])
#'
#' # Rename all clusters
#' pancreas_sub <- RenameClusters(pancreas_sub, group.by = "SubCellType", newnames = letters[1:8], createnew = "newclusters")
#' pancreas_sub@meta.data[["newclusters"]] <- factor(pancreas_sub@meta.data[["newclusters"]], levels = letters[1:8])
#' ClassDimPlot(pancreas_sub, "newclusters")
#'
#' # Rename specified clusters
#' pancreas_sub <- RenameClusters(pancreas_sub, group.by = "SubCellType", newnames = c("Alpha" = "a", "Beta" = "b"), createnew = "newclusters")
#' ClassDimPlot(pancreas_sub, "newclusters")
#'
#' # Merge and rename clusters
#' pancreas_sub <- RenameClusters(pancreas_sub, group.by = "SubCellType", nameslist = list("EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")), createnew = "newclusters")
#' ClassDimPlot(pancreas_sub, "newclusters")
#'
#' @importFrom stats setNames
#' @export
RenameClusters <- function(srt, group.by, newnames = NULL, nameslist = list(), createnew = NULL) {
  if (missing(group.by)) {
    stop("group.by must be provided")
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    stop(paste0(group.by, " is not in the meta.data of srt object."))
  }
  if (length(nameslist) > 0) {
    names_assign <- setNames(rep(names(nameslist), sapply(nameslist, length)), nm = unlist(nameslist))
  } else {
    if (is.null(names(newnames))) {
      if (!is.factor(srt@meta.data[[group.by]])) {
        stop("'newnames' must be named when srt@meta.data[[group.by]] is not a factor")
      }
      if (!identical(length(newnames), length(unique(srt@meta.data[[group.by]])))) {
        stop("'newnames' must be named or the length of ", length(unique(srt@meta.data[[group.by]])))
      }
      names(newnames) <- levels(srt@meta.data[[group.by]])
    }
    names_assign <- newnames
  }
  if (is.factor(srt@meta.data[[group.by]])) {
    levels <- levels(srt@meta.data[[group.by]])
  } else {
    levels <- NULL
  }
  if (is.null(createnew)) {
    createnew <- group.by
  }
  index <- which(as.character(srt@meta.data[[group.by]]) %in% names(names_assign))
  srt@meta.data[[createnew]] <- as.character(srt@meta.data[[group.by]])
  srt@meta.data[[createnew]][index] <- names_assign[srt@meta.data[[createnew]][index]]
  if (!is.null(levels)) {
    levels[levels %in% names(names_assign)] <- names_assign[levels[levels %in% names(names_assign)]]
    srt@meta.data[[createnew]] <- factor(srt@meta.data[[createnew]], levels = unique(levels))
  }
  return(srt)
}


#' Reorder idents by the gene expression
#'
#' @param srt
#' @param features
#' @param reorder_by
#' @param slot
#' @param assay
#' @param log
#' @param distance_metric
#' @param reorder_FUN
#'
#' @importFrom Seurat VariableFeatures DefaultAssay DefaultAssay<- AverageExpression Idents<-
#' @importFrom SeuratObject as.sparse
#' @importFrom stats hclust reorder as.dendrogram as.dist
#' @importFrom Matrix t colMeans
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

  data.avg <- AverageExpression(object = srt, features = features, slot = slot, assays = assay, group.by = "ident", verbose = FALSE)[[1]][features, ]
  if (isTRUE(log)) {
    data.avg <- log1p(data.avg)
  }
  mat <- t(x = data.avg[features, ])
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as.sparse(mat[1:nrow(mat), ])
  }

  if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 - proxyC::simil(as.sparse(mat[1:nrow(mat), ]), method = distance_metric)
  } else if (distance_metric %in% dist_method) {
    d <- proxyC::dist(as.sparse(mat[1:nrow(mat), ]), method = distance_metric)
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

#' Append a Seurat object to another
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
  if (!inherits(srt_raw, "Seurat") || !inherits(srt_append, "Seurat")) {
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

#' RunDimReduction
#'
#' @param srt
#' @param prefix
#' @param features
#' @param assay
#' @param linear_reduction "pca", "ica", "nmf", "mds", "glmpca"
#' @param linear_reduction_dims
#' @param force_linear_reduction
#' @param nonlinear_reduction "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"
#' @param reduction_use
#' @param reduction_dims
#' @param nonlinear_reduction_dims
#' @param verbose
#' @param seed
#' @param slot
#' @param linear_reduction_params
#' @param neighbor_use
#' @param graph_use
#' @param distance_use
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#'
#' @importFrom Seurat Embeddings RunPCA RunICA RunTSNE Reductions DefaultAssay DefaultAssay<- Key Key<-
#' @export
RunDimReduction <- function(srt, prefix = "", features = NULL, assay = NULL, slot = "data",
                            linear_reduction = NULL, linear_reduction_dims = 100,
                            linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = NULL, nonlinear_reduction_dims = 2,
                            reduction_use = NULL, reduction_dims = NULL, neighbor_use = NULL,
                            graph_use = NULL, distance_use = NULL,
                            nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            verbose = TRUE, seed = 11) {
  set.seed(seed)
  if (!is.null(linear_reduction)) {
    if (any(!linear_reduction %in% c("pca", "ica", "nmf", "mds", "glmpca", Reductions(srt))) || length(linear_reduction) > 1) {
      stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
    }
    if (is.null(linear_reduction_dims)) {
      linear_reduction_dims <- 100
    }
  }
  if (!is.null(nonlinear_reduction)) {
    if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", Reductions(srt))) || length(nonlinear_reduction) > 1) {
      stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
    }
    if (is.null(features) && is.null(reduction_use) && is.null(graph_use)) {
      stop("'features', 'reduction_use' or 'graph_use' must be provided when running nonlinear reduction.")
    }
  }
  assay <- assay %||% DefaultAssay(srt)
  if (!is.null(linear_reduction)) {
    if (!isTRUE(force_linear_reduction)) {
      if (linear_reduction %in% Reductions(srt)) {
        if (srt[[linear_reduction]]@assay.used == assay) {
          message("linear_reduction(", linear_reduction, ") is already existed. Skip calculation.")
          reduc <- srt[[linear_reduction]]
          Key(reduc) <- paste0(prefix, linear_reduction, "_")
          srt[[paste0(prefix, linear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
          return(srt)
        } else {
          message("assay.used is ", srt[[linear_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the linear reduction")
        }
      }
    }
    if (is.null(features) || length(features) == 0) {
      message("No features provided. Use variable features.")
      if (length(VariableFeatures(srt)) == 0) {
        srt <- FindVariableFeatures(srt, assay = assay, verbose = FALSE)
      }
      features <- VariableFeatures(srt, assay = assay)
    }
    fun_use <- switch(linear_reduction,
      "pca" = "RunPCA",
      "ica" = "RunICA",
      "nmf" = "RunNMF",
      "mds" = "RunMDS",
      "glmpca" = "RunGLMPCA"
    )
    key_use <- switch(linear_reduction,
      "pca" = "PC_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    components_nm <- switch(linear_reduction,
      "pca" = "npcs",
      "ica" = "nics",
      "nmf" = "nbes",
      "mds" = "nmds",
      "glmpca" = "L"
    )
    params <- list(
      object = srt, assay = assay, slot = slot,
      features = features, components_nm = linear_reduction_dims,
      reduction.name = paste0(prefix, linear_reduction),
      reduction.key = paste0(prefix, key_use),
      verbose = verbose, seed.use = seed
    )
    if (fun_use == "RunICA") {
      params <- params[!names(params) %in% "slot"]
    }
    if (fun_use == "RunGLMPCA") {
      params[["slot"]] <- "counts"
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(linear_reduction_params)) {
      params[[nm]] <- linear_reduction_params[[nm]]
    }
    srt <- invoke(.fn = fun_use, .args = params)

    if (is.null(rownames(srt[[paste0(prefix, linear_reduction)]]))) {
      rownames(srt[[paste0(prefix, linear_reduction)]]@cell.embeddings) <- colnames(srt)
    }
    if (linear_reduction == "pca") {
      pca.out <- srt[[paste0(prefix, linear_reduction)]]
      center <- rowMeans(GetAssayData(object = srt, slot = "scale.data")[features, ])
      model <- list(sdev = pca.out@stdev, rotation = pca.out@feature.loadings, center = center, scale = FALSE, x = pca.out@cell.embeddings)
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["model"]] <- model
    }
    if (linear_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:linear_reduction_dims
    } else {
      dim_est <- tryCatch(expr = {
        intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(prefix, linear_reduction)), k = 20, iterations = 100)[["dim.est"]]
      }, error = function(e) {
        message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
        return(NA)
      })
      if (!is.na(dim_est)) {
        if (ceiling(dim_est) < 10) {
          dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 10))
        } else {
          dims_estimate <- seq_len(ceiling(dim_est))
        }
      } else {
        dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 30))
      }
    }
    message("dims_estimate is ", paste0(range(dims_estimate), collapse = ":"), " for '", linear_reduction, "'")
    srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["dims_estimate"]] <- dims_estimate
    srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
  } else if (!is.null(nonlinear_reduction)) {
    if (!isTRUE(force_nonlinear_reduction)) {
      if (nonlinear_reduction %in% Reductions(srt)) {
        if (srt[[nonlinear_reduction]]@assay.used == assay) {
          message("nonlinear_reduction(", nonlinear_reduction, ") is already existed. Skip calculation.")
          reduc <- srt[[nonlinear_reduction]]
          Key(reduc) <- paste0(prefix, nonlinear_reduction, "_")
          srt[[paste0(prefix, nonlinear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction)
          return(srt)
        } else {
          message("assay.used is ", srt[[nonlinear_reduction]]@assay.used, ", which is not the same as the ", assay, " specified. Recalculate the linear reduction")
        }
      }
    }
    if (!is.null(neighbor_use) && !nonlinear_reduction %in% c("umap", "umap-naive")) {
      stop("'neighbor_use' only support 'umap' or 'umap-naive' method")
    }
    if (!is.null(graph_use) && !nonlinear_reduction %in% c("umap", "umap-naive")) {
      stop("'graph_use' only support 'umap' or 'umap-naive' method")
    }
    if (!is.null(distance_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "trimap", "largevis")) {
      stop("'distance_use' only support 'umap', 'umap-naive', 'trimap', 'largevis' method")
    }
    fun_use <- switch(nonlinear_reduction,
      "umap" = "RunUMAP2",
      "umap-naive" = "RunUMAP2",
      "tsne" = "RunTSNE",
      "dm" = "RunDM",
      "phate" = "RunPHATE",
      "pacmap" = "RunPaCMAP",
      "trimap" = "RunTriMap",
      "largevis" = "RunLargeVis"
    )
    components_nm <- switch(nonlinear_reduction,
      "umap" = "n.components",
      "umap-naive" = "n.components",
      "tsne" = "dim.embed",
      "dm" = "ndcs",
      "phate" = "n_components",
      "pacmap" = "n_components",
      "trimap" = "n_components",
      "largevis" = "n_components"
    )
    other_params <- switch(nonlinear_reduction,
      "umap" = list(umap.method = "uwot", return.model = TRUE),
      "umap-naive" = list(umap.method = "naive", return.model = TRUE),
      "tsne" = list(tsne.method = "Rtsne", num_threads = 0, check_duplicates = FALSE),
      "dm" = list(),
      "phate" = list(),
      "pacmap" = list(),
      "trimap" = list(),
      "largevis" = list()
    )
    nonlinear_reduction_sim <- toupper(gsub(pattern = "-.*", replacement = "", x = nonlinear_reduction))
    params <- list(
      object = srt, assay = assay, slot = slot, components_nm = nonlinear_reduction_dims,
      features = features, reduction = reduction_use, dims = reduction_dims,
      reduction.name = paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D"),
      reduction.key = paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D_"),
      verbose = verbose, seed.use = seed
    )
    if (!is.null(neighbor_use)) {
      params[["neighbor"]] <- neighbor_use
    }
    if (!is.null(graph_use)) {
      params[["graph"]] <- graph_use
    }
    if (!is.null(distance_use)) {
      params[["distance"]] <- distance_use
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(other_params)) {
      params[[nm]] <- other_params[[nm]]
    }
    for (nm in names(nonlinear_reduction_params)) {
      params[[nm]] <- nonlinear_reduction_params[[nm]]
    }
    srt <- invoke(.fn = fun_use, .args = params)

    srt@reductions[[paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D")]]@misc[["reduction_dims"]] <- reduction_dims
    srt@reductions[[paste0(prefix, nonlinear_reduction_sim, nonlinear_reduction_dims, "D")]]@misc[["reduction_use"]] <- reduction_use
    srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction_sim)
  }
  return(srt)
}

#' Uncorrected_integrate
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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param seed
#'
#' @importFrom Seurat GetAssayData SetAssayData VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Uncorrected_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                                  do_normalization = NULL, normalization_method = "logCPM",
                                  do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                                  do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                  linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                                  nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                  do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                  seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "edge"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis','edge','date','dca','sauice','scVI'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Uncorrected) on the data...\n"))
  srtIntegrated <- Standard_SCP(
    srt = srtMerge, prefix = "Uncorrected",
    do_normalization = do_normalization, normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
    do_scaling = do_scaling, vars_to_regress = vars_to_regress, regression_model = regression_model,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_dims_use = linear_reduction_dims_use, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    nonlinear_reduction = nonlinear_reduction, nonlinear_reduction_dims = nonlinear_reduction_dims, nonlinear_reduction_params = nonlinear_reduction_params, force_nonlinear_reduction = force_nonlinear_reduction,
    do_cluster_finding = do_cluster_finding, cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = cluster_reorder,
    seed = seed
  )
  srtIntegrated[[paste0("Uncorrectedclusters")]] <- srtIntegrated[[paste0("Uncorrected", linear_reduction, "clusters")]]
  srtIntegrated[[paste0("Uncorrected", linear_reduction, "clusters")]] <- NULL
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      reduc <- srtIntegrated@reductions[[paste0("Uncorrected", linear_reduction, toupper(nr), n, "D")]]
      Key(reduc) <- paste0("Uncorrected", toupper(nr), n, "D_")
      srtIntegrated@reductions[[paste0("Uncorrected", linear_reduction, toupper(nr), n, "D")]] <- NULL
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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param FindIntegrationAnchors_params
#' @param IntegrateData_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData FindIntegrationAnchors IntegrateData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents
#' @importFrom dplyr "%>%"
#' @export
Seurat_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                             do_normalization = NULL, normalization_method = "logCPM",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                             linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                             FindIntegrationAnchors_params = list(), IntegrateData_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    cat(paste0("[", Sys.time(), "]", " Use 'rpca' workflow...\n"))
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
        cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
        srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
      }
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (pca) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt = srt, prefix = "", features = HVF,
        linear_reduction = "pca", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srtList[[i]] <- srt
    }
  }
  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
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
  srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)

  cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
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
  srtIntegrated <- invoke(.fn = IntegrateData, .args = params2)

  DefaultAssay(srtIntegrated) <- "Seurat"
  VariableFeatures(srtIntegrated[["Seurat"]]) <- HVF

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(object = srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Seurat", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("Seurat", linear_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Seurat", linear_reduction), dims = linear_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("Seurat", linear_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("Seurat", linear_reduction, "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["Seuratclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "Seurat",
            reduction_use = paste0("Seurat", linear_reduction), reduction_dims = linear_reduction_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param force_nonlinear_reduction
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param SCVI_params
#' @param seed
#' @param nonlinear_reduction_params
#' @param num_threads
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
scVI_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                           do_normalization = NULL, normalization_method = "logCPM",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                           nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                           do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                           SCVI_params = list(), num_threads = 8, seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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
  if (.Platform$OS.type == "windows" && !exist_pkg(pkg = "scvi-tools", envname = "SCP")) {
    suppressWarnings(system2(command = reticulate::virtualenv_python("SCP"), args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver", stdout = TRUE))
  }

  check_Python("scvi-tools")
  scvi <- import("scvi")
  scipy <- import("scipy")
  set.seed(seed)

  scvi$settings$num_threads <- as.integer(num_threads)

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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  adata <- srt_to_adata(srtMerge[HVF, ], assay_X = DefaultAssay(srtMerge), assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])
  scvi$model$SCVI$setup_anndata(adata, batch_key = batch)

  params <- list(
    adata = adata
  )
  for (nm in names(SCVI_params)) {
    params[[nm]] <- SCVI_params[[nm]]
  }
  model <- invoke(.fn = scvi$model$SCVI, .args = params)
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

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "scVI", dims = scVI_dims_use, force.recalc = TRUE, graph.name = paste0("scVI", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("scVI", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["scVIclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "scVI",
            reduction_use = "scVI", reduction_dims = scVI_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param do_scaling
#' @param vars_to_regress
#' @param regression_model
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param linear_reduction_params
#' @param force_linear_reduction
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param mnnCorrect_params
#' @param seed
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
MNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          mnnCorrect_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(MNN) on the data...\n"))
  params <- list(
    sceList,
    cos.norm.out = FALSE,
    BPPARAM = BiocParallel::bpparam()
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- invoke(.fn = batchelor::mnnCorrect, .args = params)

  srtIntegrated[["MNNcorrected"]] <- CreateAssayObject(counts = out@assays@data$corrected)
  VariableFeatures(srtIntegrated[["MNNcorrected"]]) <- HVF
  DefaultAssay(srtIntegrated) <- "MNNcorrected"

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(object = srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "MNN", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("MNN", linear_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("MNN", linear_reduction), dims = linear_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("MNN", linear_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("MNN", linear_reduction, "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["MNNclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "MNN",
            reduction_use = paste0("MNN", linear_reduction), reduction_dims = linear_reduction_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param fastMNN_dims_use
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param fastMNN_params
#' @param seed
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
fastMNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                              fastMNN_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              fastMNN_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(fastMNN) on the data...\n"))
  params <- list(
    sceList,
    BPPARAM = BiocParallel::bpparam()
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke(.fn = batchelor::fastMNN, .args = params)

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

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "fastMNN", dims = fastMNN_dims_use, force.recalc = TRUE, graph.name = paste0("fastMNN", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("fastMNN", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["fastMNNclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "fastMNN",
            reduction_use = "fastMNN", reduction_dims = fastMNN_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param Harmony_dims_use
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param RunHarmony_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Harmony_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                              do_normalization = NULL, normalization_method = "logCPM",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                              linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                              Harmony_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              RunHarmony_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "Harmony", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("Harmony", linear_reduction)]]@misc[["dims_estimate"]]
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Harmony) on the data...\n"))
  params <- list(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = DefaultAssay(srtMerge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srtIntegrated <- invoke(.fn = RunHarmony2, .args = params)

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

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Harmony", dims = Harmony_dims_use, force.recalc = TRUE, graph.name = paste0("Harmony", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("Harmony", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["Harmonyclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "Harmony",
            reduction_use = "Harmony", reduction_dims = Harmony_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            force_nonlinear_reduction = force_nonlinear_reduction,
            nonlinear_reduction_params = nonlinear_reduction_params,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param scanorama_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @importFrom Matrix t
#' @importFrom dplyr "%>%"
#' @importFrom reticulate import
#' @importFrom stats sd
#' @export
Scanorama_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                                do_normalization = NULL, normalization_method = "logCPM",
                                do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                                do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                                nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                scanorama_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Scanorama) on the data...\n"))
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
    verbose = FALSE
  )
  for (nm in names(scanorama_params)) {
    params[[nm]] <- scanorama_params[[nm]]
  }
  integrated.corrected.data <- invoke(.fn = scanorama$correct, .args = params)

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

  cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
  srtIntegrated <- ScaleData(srtIntegrated, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "Scanorama", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtIntegrated@reductions[[paste0("Scanorama", linear_reduction)]]@misc[["dims_estimate"]]
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("Scanorama", linear_reduction), dims = linear_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("Scanorama", linear_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("Scanorama", linear_reduction, "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["Scanoramaclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "Scanorama",
            reduction_use = paste0("Scanorama", linear_reduction), reduction_dims = linear_reduction_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param bbknn_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- as.Graph Embeddings FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom Matrix t
#' @importFrom dplyr "%>%"
#' @importFrom reticulate import
#' @export
BBKNN_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            bbknn_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "BBKNN", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("BBKNN", linear_reduction)]]@misc[["dims_estimate"]]
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(BBKNN) on the data...\n"))
  emb <- Embeddings(srtMerge, reduction = paste0("BBKNN", linear_reduction))[, linear_reduction_dims_use]
  params <- list(
    pca = emb,
    batch_list = srtMerge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke(.fn = bbknn$bbknn_matrix, .args = params)

  bbknn_graph <- as.sparse(bem[[2]][1:nrow(bem[[2]]), ])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- DefaultAssay(srtMerge)
  srtMerge@graphs[["BBKNN"]] <- bbknn_graph
  srtIntegrated <- srtMerge

  srtIntegrated <- tryCatch(
    {
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "BBKNN", resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["BBKNNclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "BBKNN",
            graph_use = "BBKNN",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param linear_reduction
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param CSS_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
CSS_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                          do_normalization = NULL, normalization_method = "logCPM",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          CSS_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "CSS", features = HVF,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- srtMerge@reductions[[paste0("CSS", linear_reduction)]]@misc[["dims_estimate"]]
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(CSS) on the data...\n"))
  params <- list(
    object = srtMerge,
    use_dr = paste0("CSS", linear_reduction),
    dims_use = linear_reduction_dims_use,
    var_genes = HVF,
    label_tag = batch,
    reduction.name = "CSS",
    reduction.key = "CSS_",
    verbose = FALSE
  )
  for (nm in names(CSS_params)) {
    params[[nm]] <- CSS_params[[nm]]
  }
  srtIntegrated <- invoke(.fn = simspec::cluster_sim_spectrum, .args = params)

  CSS_dims_use <- seq_len(ncol(Embeddings(srtIntegrated, reduction = "CSS")))
  srtIntegrated@reductions[["CSS"]]@misc[["dims_estimate"]] <- CSS_dims_use
  if (any(is.na(srtIntegrated@reductions[["CSS"]]@cell.embeddings))) {
    stop("NA detected in the CSS embeddings. You can try to use a lower resolution value in the CSS_param.")
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "CSS", dims = CSS_dims_use, force.recalc = TRUE, graph.name = paste0("CSS", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("CSS", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["CSSclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "CSS",
            reduction_use = "CSS", reduction_dims = CSS_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param optimizeALS_params
#' @param quantilenorm_params
#' @param seed
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
LIGER_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            optimizeALS_params = list(), quantilenorm_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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
  check_R("rliger")

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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data"))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data(do not center for LIGER)...\n"))
    srtMerge <- ScaleData(object = srtMerge, features = HVF, split.by = batch, do.center = FALSE, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  split.cells <- split(x = colnames(x = srtMerge), f = srtMerge[[batch]])
  scale.data <- lapply(X = split.cells, FUN = function(x) {
    return(t(x = GetAssayData(object = srtMerge, slot = "scale.data")[, x]))
  })

  cat(paste0("[", Sys.time(), "]", " Perform integration(LIGER) on the data...\n"))
  params1 <- list(
    object = scale.data,
    k = 20,
    verbose = FALSE
  )
  for (nm in names(optimizeALS_params)) {
    params1[[nm]] <- optimizeALS_params[[nm]]
  }
  out1 <- invoke(.fn = rliger::optimizeALS, .args = params1)
  cat("\n")
  colnames(x = out1$W) <- colnames(scale.data[[1]])
  reduction1 <- do.call(what = "rbind", args = out1$H)
  colnames(reduction1) <- paste0("riNMF_", seq_len(ncol(reduction1)))
  loadings1 <- t(x = out1$W)
  rownames(loadings1) <- colnames(scale.data[[1]])
  colnames(loadings1) <- paste0("riNMF_", seq_len(ncol(loadings1)))
  srtMerge[["iNMF_raw"]] <- CreateDimReducObject(
    embeddings = reduction1,
    loadings = loadings1,
    assay = DefaultAssay(srtMerge),
    key = "riNMF_"
  )

  embeddings <- sapply(
    X = SplitObject(object = srtMerge, split.by = batch),
    FUN = function(x) {
      return(Embeddings(object = x[["iNMF_raw"]]))
    }, simplify = FALSE, USE.NAMES = TRUE
  )
  num.samples <- vapply(X = embeddings, FUN = nrow, FUN.VALUE = integer(length = 1L))
  ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  params2 <- list(
    object = embeddings,
    ref_dataset = ref_dataset
  )
  for (nm in names(quantilenorm_params)) {
    params2[[nm]] <- quantilenorm_params[[nm]]
  }
  out2 <- invoke(.fn = rliger::quantile_norm, .args = params2)
  srtMerge[["LIGER"]] <- CreateDimReducObject(
    embeddings = out2$H.norm,
    assay = DefaultAssay(srtMerge),
    key = "LIGER_"
  )
  srtIntegrated <- srtMerge
  srtMerge <- NULL

  LIGER_dims_use <- seq_len(ncol(Embeddings(srtIntegrated, reduction = "LIGER")))
  srtIntegrated@reductions[["LIGER"]]@misc[["dims_estimate"]] <- LIGER_dims_use

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "LIGER", dims = LIGER_dims_use, force.recalc = TRUE, graph.name = paste0("LIGER", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("LIGER", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["LIGERclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "LIGER",
            reduction_use = "LIGER", reduction_dims = LIGER_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param linear_reduction_dims
#' @param linear_reduction_dims_use
#' @param force_linear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param linear_reduction
#' @param linear_reduction_params
#' @param nonlinear_reduction
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param buildGraph_params
#' @param num_threads
#' @param seed
#'
#' @return
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
#' @export
Conos_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            buildGraph_params = list(), num_threads = 2, seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive'.")
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
  }
  srtIntegrated <- Reduce(merge, srtList)
  VariableFeatures(srtIntegrated) <- HVF

  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
      srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data ", i, " ...\n"))
    srt <- RunDimReduction(
      srt = srt, prefix = "", features = HVF,
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    srt[["pca"]] <- srt[[linear_reduction]]
    srtList[[i]] <- srt
  }
  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:max(unlist(lapply(srtList, function(srt) srt@reductions[[linear_reduction]]@misc[["dims_estimate"]])))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Conos) on the data...\n"))
  srtList_con <- conos::Conos$new(srtList, n.cores = num_threads)
  params <- list(
    ncomps = max(linear_reduction_dims_use),
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invoke(.fn = srtList_con[["buildGraph"]], .args = params)

  conos_graph <- as.matrix(srtList_con$graph)
  conos_graph <- as.Graph(conos_graph)
  conos_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["Conos"]] <- conos_graph

  srtIntegrated <- tryCatch(
    {
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, graph.name = "Conos", resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["Conosclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat("Perform nonlinear dimension reduction (", nr, ") on the data...\n", sep = "")
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "Conos",
            graph_use = "Conos",
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' @param nonlinear_reduction
#' @param nonlinear_reduction_dims
#' @param do_cluster_finding
#' @param cluster_algorithm
#' @param cluster_resolution
#' @param cluster_reorder
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param zinbwave_params
#' @param seed
#'
#' @return
#' @importFrom Seurat CreateSeuratObject as.SingleCellExperiment GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%"
ZINBWaVE_integrate <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                               do_normalization = NULL, normalization_method = "logCPM",
                               do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                               do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                               nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                               do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                               zinbwave_params = list(), seed = 11) {
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg", envname = "SCP")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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
      HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_intersect = HVF_intersect, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
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

  cat(paste0("[", Sys.time(), "]", " Perform integration(ZINBWaVE) on the data...\n"))
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
    sce_zinbwave <- invoke(.fn = zinbwave::zinbwave, .args = params)
  } else {
    state <- 1
    while (state != 0) {
      tryCatch(expr = {
        sce_zinbwave <- invoke(zinbwave::zinbsurf, .args = params)
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

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "ZINBWaVE", dims = ZINBWaVE_dims_use, force.recalc = TRUE, graph.name = paste0("ZINBWaVE", "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("ZINBWaVE", "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, slot = "data", features = HVF, reorder_by = "seurat_clusters")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["ZINBWaVEclusters"]] <- Idents(srtIntegrated)
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing FindClusters. Skip this step...")
      return(srtIntegrated)
    }
  )

  srtIntegrated <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
        for (n in nonlinear_reduction_dims) {
          srtIntegrated <- RunDimReduction(
            srt = srtIntegrated, prefix = "ZINBWaVE",
            reduction_use = "ZINBWaVE", reduction_dims = linear_reduction_dims_use,
            nonlinear_reduction = nr, nonlinear_reduction_dims = n,
            nonlinear_reduction_params = nonlinear_reduction_params,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = FALSE, seed = seed
          )
        }
      }
      srtIntegrated
    },
    error = function(error) {
      message(error)
      message("Error when performing nonlinear dimension reduction. Skip this step...")
      return(srtIntegrated)
    }
  )

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
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(srt = pancreas_sub)
#' ClassDimPlot(pancreas_sub, group.by = "CellType")
#'
#' # Use a combination of different linear or non-linear dimension reduction methods
#' pancreas_sub <- Standard_SCP(
#'   srt = pancreas_sub,
#'   linear_reduction = c("pca", "ica", "nmf", "mds", "glmpca"),
#'   nonlinear_reduction = c("umap", "tsne", "phate", "pacmap", "trimap")
#' )
#' names(pancreas_sub@reductions)
#'
#' @importFrom Seurat Assays GetAssayData NormalizeData SCTransform SCTResults ScaleData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom dplyr "%>%" filter arrange desc
#' @importFrom Matrix rowSums
#' @export
Standard_SCP <- function(srt, prefix = "Standard",
                         do_normalization = NULL, normalization_method = "logCPM",
                         do_HVF_finding = TRUE, HVF_method = "vst", nHVF = 2000, HVF = NULL,
                         do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                         linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                         nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                         do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                         seed = 11) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  if (!normalization_method %in% c("logCPM", "SCT", "Linnorm", "Scran", "Scone", "DESeq2")) {
    stop("'normalization_method' must be one of: 'logCPM','SCT','Linnorm','Scran','Scone','DESeq2'")
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  reduc_test <- c(reduc_test, Reductions(srt))
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
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

  linear_reduction_dims <- min(linear_reduction_dims, ncol(srt) - 1)
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data"))))) {
    if (normalization_method != "SCT") {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srt <- ScaleData(object = srt, features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
  }

  for (lr in linear_reduction) {
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", lr, ") on the data...\n"))
    srt <- RunDimReduction(
      srt = srt, prefix = prefix, features = HVF,
      linear_reduction = lr, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- srt@reductions[[paste0(prefix, lr)]]@misc[["dims_estimate"]]
    }

    srt <- tryCatch(
      {
        srt <- FindNeighbors(object = srt, reduction = paste0(prefix, lr), dims = linear_reduction_dims_use, force.recalc = TRUE, graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")), verbose = FALSE)
        if (isTRUE(do_cluster_finding)) {
          cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
          srt <- FindClusters(object = srt, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0(prefix, lr, "_SNN"), verbose = FALSE)
          if (isTRUE(cluster_reorder)) {
            cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
            srt <- SrtReorder(srt, features = HVF, reorder_by = "seurat_clusters", slot = "data")
          }
          srt[["seurat_clusters"]] <- NULL
          srt[[paste0(prefix, lr, "clusters")]] <- Idents(srt)
        }
        srt
      },
      error = function(error) {
        message(error)
        message("Error when performing FindClusters. Skip this step...")
        return(srt)
      }
    )

    srt <- tryCatch(
      {
        for (nr in nonlinear_reduction) {
          cat(paste0("[", Sys.time(), "]", " Perform nonlinear dimension reduction (", nr, ") on the data...\n"))
          for (n in nonlinear_reduction_dims) {
            srt <- RunDimReduction(
              srt = srt, prefix = paste0(prefix, lr),
              reduction_use = paste0(prefix, lr), reduction_dims = linear_reduction_dims_use,
              nonlinear_reduction = nr, nonlinear_reduction_dims = n,
              nonlinear_reduction_params = nonlinear_reduction_params,
              force_nonlinear_reduction = force_nonlinear_reduction,
              verbose = FALSE, seed = seed
            )
          }
        }
        srt
      },
      error = function(error) {
        message(error)
        message("Error when performing nonlinear dimension reduction. Skip this step...")
        return(srt)
      }
    )
  }

  if (paste0(prefix, linear_reduction[1], "clusters") %in% colnames(srt@meta.data)) {
    srt[[paste0(prefix, "clusters")]] <- srt[[paste0(prefix, linear_reduction[1], "clusters")]]
  }
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      if (paste0(prefix, linear_reduction[1], toupper(nr), n, "D") %in% names(srt@reductions)) {
        reduc <- srt@reductions[[paste0(prefix, linear_reduction[1], toupper(nr), n, "D")]]
        srt@reductions[[paste0(prefix, toupper(nr), n, "D")]] <- reduc
      }
    }
    srt@misc[["Default_reduction"]] <- paste0(prefix, toupper(nr))
  }

  DefaultAssay(srt) <- "RNA"
  VariableFeatures(srt) <- srt@misc[["Standard_HVF"]] <- HVF

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
#' @param integration_method Integration method. Can be one of "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos".
#' @param do_normalization Whether to normalize the data. If NULL, will automatically determine.
#' @param normalization_method Normalization method.Can be one of "logCPM", "SCT".
#' @param do_HVF_finding Whether to find the high variable features(HVF). If NULL, will automatically determine.
#' @param HVF_source Source of the HVF. Can be one of "separate" and "global".
#' @param nHVF HVF number to use.
#' @param HVF Custom high variable features.
#' @param do_scaling Whether to scale the data. If NULL, will automatically determine.
#' @param vars_to_regress Variables to regress out.
#' @param regression_model Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are "linear" (default), "poisson", and "negbinom".
#' @param linear_reduction linear reduction method name. Can be one of "pca", "ica", "nmf", "mds", "glmpca".
#' @param linear_reduction_dims Dimensions to calculate when performing linear reduction.
#' @param linear_reduction_dims_use Which dimensions to use when performing the nonlinear reduction.
#' @param nonlinear_reduction Non-linear reduction method name. Can be one of "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis".
#' @param nonlinear_reduction_dims Dimensions to calculate when performing non-linear reduction.
#' @param cluster_algorithm Algorithm for modularity optimization when finding clusters. Can be one of "louvain", "slm", "leiden".
#' @param cluster_resolution Cluster resolution parameter.
#' @param cluster_reorder Whether to reorder the cluster names using hierarchical clustering.
#' @param seed Set a random seed.
#' @param HVF_method
#' @param force_linear_reduction
#' @param do_cluster_finding
#' @param ...
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#'
#' @return A \code{Seurat} object containing the result.
#'
#' @examples
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Uncorrected"
#' )
#' ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Seurat"
#' )
#' ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Seurat",
#'   HVF_intersect = TRUE
#' )
#' ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Seurat",
#'   HVF_intersect = TRUE, HVF_min_intersection = length(unique(panc8_sub$tech))
#' )
#' ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' \dontrun{
#' for (method in c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos")) {
#'   panc8_sub <- Integration_SCP(panc8_sub,
#'     batch = "tech", integration_method = method,
#'     nonlinear_reduction = "umap"
#'   )
#'   print(ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = paste0(method, "UMAP2D"), theme_use = "theme_blank"))
#' }
#'
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Seurat",
#'   nonlinear_reduction = c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")
#' )
#' for (reduc in c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")) {
#'   print(ClassDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = paste0("Seurat", reduc, "2D"), theme_use = "theme_blank"))
#' }
#' }
#'
#' @export
Integration_SCP <- function(srtMerge = NULL, batch = "orig.ident", append = TRUE, srtList = NULL,
                            integration_method = "Uncorrected",
                            do_normalization = NULL, normalization_method = "logCPM",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_intersect = FALSE, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 100, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            seed = 11, ...) {
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("Must be provided with one of the 'srtList' and 'srtMerge'")
  }
  if (length(integration_method) == 1 && integration_method %in% c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos")) {
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
