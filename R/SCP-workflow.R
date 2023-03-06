#' Check and report the type of data
#'
#' @param srt A Seurat object
#' @param data Use a data instead of getting it from a Seurat object.
#' @param slot Specific slot to get data from.
#' @param assay Specific assay to get data from.
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
check_DataType <- function(srt, data = NULL, slot = "data", assay = NULL) {
  if (is.null(data)) {
    assay <- assay %||% DefaultAssay(srt)
    data <- GetAssayData(srt, slot = slot, assay = assay)
  }
  isfinite <- all(is.finite(range(data, na.rm = TRUE)))
  if (inherits(data, "dgCMatrix")) {
    isfloat <- any(data@x %% 1 != 0, na.rm = TRUE)
  } else {
    isfloat <- any(data[, sample(seq_len(ncol(data)), min(ncol(data), 1000))] %% 1 != 0, na.rm = TRUE)
  }
  islog <- is.finite(expm1(x = max(data, na.rm = TRUE)))
  isnegative <- any(data < 0)

  if (!isTRUE(isfinite)) {
    warning("Infinite values detected!", immediate. = TRUE)
    return("unknown")
  } else if (isTRUE(isnegative)) {
    warning("Negative values detected!", immediate. = TRUE)
    return("unknown")
  } else {
    if (!isfloat) {
      return("raw_counts")
    } else if (isfloat && islog) {
      return("log_normalized_counts")
    } else if (isfloat && !islog) {
      if (isFALSE(isnegative)) {
        return("raw_normalized_counts")
      } else {
        return("unknown")
      }
    }
  }
}

#' Check and preprocess a list of seurat objects
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat SplitObject GetAssayData Assays NormalizeData FindVariableFeatures SCTransform SCTResults SelectIntegrationFeatures PrepSCTIntegration DefaultAssay DefaultAssay<- VariableFeatures VariableFeatures<-
#' @importFrom Signac RunTFIDF
#' @importFrom Matrix rowSums
#' @importFrom utils head
#' @export
#'
check_srtList <- function(srtList, batch = NULL, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                          nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          vars_to_regress = NULL, seed = 11, ...) {
  cat(paste0("[", Sys.time(), "]", " Checking srtList... ...\n"))
  set.seed(seed)

  if (!inherits(srtList, "list") || any(sapply(srtList, function(x) !inherits(x, "Seurat")))) {
    stop("'srtList' is not a list of Seurat object.")
  }
  if (!normalization_method %in% c("LogNormalize", "SCT", "TFIDF")) {
    stop("'normalization_method' must be one of: 'LogNormalize', 'SCT', 'TFIDF'")
  }
  if (normalization_method %in% c("SCT")) {
    check_R("glmGamPoi")
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global', 'separate'")
  }
  if (any(sapply(srtList, ncol) < 3)) {
    stop(paste0(
      "Seurat objects in srtList contain less than 3 cells. srtList index: ",
      paste0(which(sapply(srtList, ncol) < 3), collapse = ",")
    ))
  }

  if (is.null(assay)) {
    default_assay <- unique(sapply(srtList, DefaultAssay))
    if (length(default_assay) != 1) {
      stop("The default assay name of the Seurat object in the srtlist is inconsistent.")
    } else {
      assay <- default_assay
    }
  }

  assay_type <- unique(sapply(srtList, function(srt) class(srt[[assay]])))
  if (length(assay_type) != 1) {
    stop("The assay type of the Seurat object in the srtlist is inconsistent.")
  } else {
    if (assay_type == "Assay") {
      type <- "RNA"
    } else if (assay_type == "ChromatinAssay") {
      type <- "Chromatin"
    }
  }

  features_list <- lapply(srtList, function(srt) {
    sort(rownames(srt[[assay]]))
  })
  if (length(unique(features_list)) != 1) {
    if (type == "Chromatin") {
      warning("The peaks in assay ", assay, " is different between batches. Creating a common set of peaks and may take a long time...")
      srtMerge <- Reduce(merge, srtList)
      srtList <- SplitObject(object = srtMerge, split.by = batch)
    }
    cf <- Reduce(intersect, lapply(srtList, function(srt) rownames(srt[[assay]])))
    warning("'srtList' have different feature names! Will subset the common features(", length(cf), ") for downstream analysis!", immediate. = TRUE)
    for (i in seq_along(srtList)) {
      srtList[[i]][[assay]] <- subset(srtList[[i]][[assay]], features = cf)
    }
  }

  celllist <- unlist(lapply(srtList, colnames))
  if (length(celllist) != length(unique(celllist))) {
    stop("'srtList' have duplicated cell names!")
  }

  if (length(batch) != 1 && length(batch) != length(srtList)) {
    stop("'batch' must be a character to specify the batch column in the meta.data or a vector of the same length of the srtList!")
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
        srtList[[i]] <- character(0)
        srtList <- c(srtList, x)
      }
    }
    srtList <- srtList[sapply(srtList, length) > 0]
    srtList_batch <- sapply(srtList, function(x) unique(x[[batch, drop = TRUE]]))
    batch_to_merge <- names(which(table(srtList_batch) > 1))
    if (length(batch_to_merge) > 0) {
      for (b in batch_to_merge) {
        index <- which(srtList_batch == b)
        srtList_tmp <- Reduce(merge, srtList[index])
        for (i in index) {
          srtList[[i]] <- character(0)
        }
        srtList <- c(srtList, srtList_tmp)
      }
    }
    srtList <- srtList[sapply(srtList, length) > 0]
  }

  for (i in seq_along(srtList)) {
    if (!assay %in% Assays(srtList[[i]])) {
      stop(paste0("srtList ", i, " does not contain '", assay, "' assay."))
    }
    DefaultAssay(srtList[[i]]) <- assay
    if (isTRUE(do_normalization)) {
      if (normalization_method == "LogNormalize") {
        cat("Perform NormalizeData(LogNormalize) on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- suppressWarnings(NormalizeData(object = srtList[[i]], assay = assay, normalization.method = "LogNormalize", verbose = FALSE))
      }
      if (normalization_method == "TFIDF") {
        cat("Perform RunTFIDF on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- suppressWarnings(RunTFIDF(object = srtList[[i]], assay = assay, verbose = FALSE))
      }
    } else if (is.null(do_normalization)) {
      status <- check_DataType(srtList[[i]], slot = "data", assay = assay)
      if (status == "log_normalized_counts") {
        cat("Data ", i, "/", length(srtList), " of the srtList has been log-normalized.\n", sep = "")
      }
      if (status %in% c("raw_counts", "raw_normalized_counts")) {
        if (normalization_method == "LogNormalize") {
          cat("Data ", i, "/", length(srtList), " of the srtList is ", status, ". Perform NormalizeData(LogNormalize) on the data ...\n", sep = "")
          srtList[[i]] <- suppressWarnings(NormalizeData(object = srtList[[i]], assay = assay, normalization.method = "LogNormalize", verbose = FALSE))
        }
        if (normalization_method == "TFIDF") {
          cat("Data ", i, "/", length(srtList), " of the srtList is ", status, ". Perform RunTFIDF on the data ...\n", sep = "")
          srtList[[i]] <- suppressWarnings(RunTFIDF(object = srtList[[i]], assay = assay, verbose = FALSE))
        }
      }
      if (status == "unknown") {
        warning("Can not determine whether data ", i, " is log-normalized...\n", immediate. = TRUE)
      }
    }
    if (is.null(HVF)) {
      if (isTRUE(do_HVF_finding) || is.null(do_HVF_finding) || length(VariableFeatures(srtList[[i]], assay = assay)) == 0) {
        # if (type == "RNA") {
        cat("Perform FindVariableFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        srtList[[i]] <- suppressWarnings(FindVariableFeatures(srtList[[i]], assay = assay, nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE))
        # }
        # if (type == "Chromatin") {
        #   cat("Perform FindTopFeatures on the data ", i, "/", length(srtList), " of the srtList...\n", sep = "")
        #   srtList[[i]] <- suppressWarnings(FindTopFeatures(srtList[[i]], assay = assay, min.cutoff = HVF_min_cutoff, verbose = FALSE))
        # }
      }
    }

    if (normalization_method %in% c("SCT") && type == "RNA") {
      if (isTRUE(do_normalization) || isTRUE(do_HVF_finding) || !"SCT" %in% Assays(srtList[[i]])) {
        cat("Perform SCTransform on the data", i, "of the srtList...\n")
        srtList[[i]] <- SCTransform(
          object = srtList[[i]],
          variable.features.n = nHVF,
          vars.to.regress = vars_to_regress,
          assay = assay,
          method = "glmGamPoi",
          new.assay.name = "SCT",
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
      } else {
        feature.attr <- srtList[[i]]@assays$SCT@meta.features
      }
      nfeatures <- min(nHVF, nrow(x = feature.attr))
      top.features <- rownames(x = feature.attr)[head(order(feature.attr$residual_variance, decreasing = TRUE), n = nfeatures)]
      VariableFeatures(srtList[[i]], assay = DefaultAssay(srtList[[i]])) <- top.features
      srtList[[i]]@assays$SCT@meta.features <- feature.attr
    }
  }

  if (is.null(HVF)) {
    if (HVF_source == "global") {
      cat("Use the global HVF from merged dataset...\n")
      srtMerge <- Reduce(merge, srtList)
      # if (type == "RNA") {
      srtMerge <- FindVariableFeatures(srtMerge, assay = DefaultAssay(srtMerge), nfeatures = nHVF, selection.method = HVF_method, verbose = FALSE)
      # }
      # if (type == "Chromatin") {
      #   srtMerge <- FindTopFeatures(srtMerge, assay = DefaultAssay(srtMerge), min.cutoff = HVF_min_cutoff, verbose = FALSE)
      # }
      HVF <- VariableFeatures(srtMerge)
    }
    if (HVF_source == "separate") {
      cat("Use the separate HVF from srtList...\n")
      # if (type == "RNA") {
      HVF <- SelectIntegrationFeatures(object.list = srtList, nfeatures = nHVF, verbose = FALSE)
      HVF_sort <- sort(table(unlist(lapply(srtList, VariableFeatures))), decreasing = TRUE)
      HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      HVF <- intersect(HVF, names(HVF_filter))
      # }
      # if (type == "Chromatin") {
      #   nHVF <- min(sapply(srtList, function(srt) length(VariableFeatures(srt))))
      #   HVF_sort <- sort(table(unlist(lapply(srtList, VariableFeatures))), decreasing = TRUE)
      #   HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      #   HVF <- names(head(HVF_filter, nHVF))
      # }
      if (length(HVF) == 0) {
        stop("No HVF available.")
      }
    }
  } else {
    cf <- Reduce(intersect, lapply(srtList, function(srt) {
      rownames(GetAssayData(srt, slot = "counts", assay = DefaultAssay(srt)))
    }))
    HVF <- HVF[HVF %in% cf]
  }
  message("Number of available HVF: ", length(HVF))

  hvf_sum <- lapply(srtList, function(srt) {
    colSums(GetAssayData(srt, slot = "counts", assay = DefaultAssay(srt))[HVF, , drop = FALSE])
  })
  cell_all <- unlist(unname(hvf_sum))
  cell_abnormal <- names(cell_all)[cell_all == 0]
  if (length(cell_abnormal) > 0) {
    warning("Some cells do not express any of the highly variable features: ", paste(cell_abnormal, collapse = ","), immediate. = TRUE)
  }

  if (normalization_method == "SCT" && type == "RNA") {
    srtList <- PrepSCTIntegration(object.list = srtList, anchor.features = HVF, assay = "SCT", verbose = FALSE)
  }
  cat(paste0("[", Sys.time(), "]", " Finished checking.\n"))

  return(list(
    srtList = srtList,
    HVF = HVF,
    assay = assay,
    type = type
  ))
}

#' Check and preprocess a merged seurat object
#' @inheritParams Integration_SCP
#' @importFrom Seurat GetAssayData SplitObject SetAssayData VariableFeatures VariableFeatures<-
#' @export
check_srtMerge <- function(srtMerge, batch = NULL, assay = NULL,
                           do_normalization = NULL, normalization_method = "LogNormalize",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst",
                           nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                           vars_to_regress = NULL, seed = 11, ...) {
  if (!inherits(srtMerge, "Seurat")) {
    stop("'srtMerge' is not a Seurat object.")
  }
  if (length(batch) != 1) {
    stop("'batch' must be provided to specify the batch column in the meta.data")
  }
  if (!batch %in% colnames(srtMerge@meta.data)) {
    stop(paste0("No batch column('", batch, "') found in the meta.data"))
  }
  if (!is.factor(srtMerge[[batch, drop = TRUE]])) {
    srtMerge[[batch, drop = TRUE]] <- factor(srtMerge[[batch, drop = TRUE]],
      levels = unique(srtMerge[[batch, drop = TRUE]])
    )
  }
  assay <- assay %||% DefaultAssay(srtMerge)
  srtMerge_raw <- srtMerge

  cat(paste0("[", Sys.time(), "]", " Spliting srtMerge into srtList by column ", batch, "... ...\n"))
  srtList <- SplitObject(object = srtMerge_raw, split.by = batch)

  checked <- check_srtList(
    srtList = srtList, batch = batch, assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source, HVF_method = HVF_method, nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srtList <- checked[["srtList"]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  srtMerge <- Reduce(merge, srtList)

  srtMerge <- SrtAppend(
    srt_raw = srtMerge, srt_append = srtMerge_raw, pattern = "",
    slots = "reductions", overwrite = TRUE, verbose = FALSE
  )
  if (normalization_method == "SCT" && type == "RNA") {
    DefaultAssay(srtMerge) <- "SCT"
  } else {
    DefaultAssay(srtMerge) <- assay
  }
  VariableFeatures(srtMerge) <- HVF

  return(list(
    srtMerge = srtMerge,
    srtList = srtList,
    HVF = HVF,
    assay = assay,
    type = type
  ))
}

#' Attempt to recover raw counts from the normalized matrix.
#'
#' @param srt A Seurat object.
#' @param assay Name of assay to recover counts.
#' @param trans The transformation function to applied when data is presumed to be log-normalized.
#' @param min_count Minimum UMI count of genes.
#' @param tolerance When recovering the raw counts, the nCount of each cell is theoretically calculated as an integer.
#'  However, due to decimal point preservation during normalization, the calculated nCount is usually a floating point number close to the integer.
#'  The tolerance is its difference from the integer. Default is 0.1
#' @param sf Set the scaling factor manually.
#' @param verbose Show messages.
#'
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
#' @importFrom Seurat GetAssayData SetAssayData
#' @importFrom SeuratObject as.sparse
#' @export
RecoverCounts <- function(srt, assay = NULL, trans = c("expm1", "exp", "none"), min_count = c(1, 2, 3), tolerance = 0.1, sf = NULL, verbose = TRUE) {
  assay <- assay %||% DefaultAssay(srt)
  counts <- GetAssayData(srt, assay = assay, slot = "counts")
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as.sparse(counts[1:nrow(counts), , drop = FALSE])
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
    trans <- match.arg(trans)
    if (trans %in% c("expm1", "exp")) {
      if (isTRUE(verbose)) {
        message("Perform ", trans, " on the raw data.")
      }
      counts <- do.call(trans, list(counts))
    }
  }
  if (status == "raw_normalized_counts") {
    if (isTRUE(verbose)) {
      message("The data is presumed to be normalized without log transformation.")
    }
  }
  if (is.null(sf)) {
    sf <- unique(round(colSums(counts)))
    if (isTRUE(verbose)) {
      message("The presumed scale factor: ", paste0(head(sf, 10), collapse = ", "))
    }
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
#' @param srt A Seurat object.
#' @param newnames A vector with the same length of features in Seurat object, or characters named with old features.
#' @param assays Assays to rename.
#'
#' @examples
#' data("panc8_sub")
#' head(rownames(panc8_sub))
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
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
    if (length(unique(sapply(pancreas_sub@assays[assays], nrow))) > 1) {
      stop("Assays in the srt object have different number of features. Please use a named vectors.")
    }
    names(newnames) <- rownames(srt[[assays[1]]])
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
#' @param srt A Seurat object.
#' @param group.by The old group used to rename cells.
#' @param nameslist A named list of new cluster value.
#' @param name The name of the new cluster stored in the Seurat object.
#' @param keep_levels If the old group is a factor, keep the order of the levels.
#'
#' @examples
#' data("pancreas_sub")
#' levels(pancreas_sub@meta.data[["SubCellType"]])
#'
#' # Rename all clusters
#' pancreas_sub <- RenameClusters(pancreas_sub, group.by = "SubCellType", nameslist = letters[1:8])
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Rename specified clusters
#' pancreas_sub <- RenameClusters(pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list("a" = "Alpha", "b" = "Beta")
#' )
#' CellDimPlot(pancreas_sub, "newclusters")
#'
#' # Merge and rename clusters
#' pancreas_sub <- RenameClusters(pancreas_sub,
#'   group.by = "SubCellType",
#'   nameslist = list("EndocrineClusters" = c("Alpha", "Beta", "Epsilon", "Delta")),
#'   name = "Merged", keep_levels = TRUE
#' )
#' CellDimPlot(pancreas_sub, "Merged")
#'
#' @importFrom stats setNames
#' @export
RenameClusters <- function(srt, group.by, nameslist = list(), name = "newclusters", keep_levels = FALSE) {
  if (missing(group.by)) {
    stop("group.by must be provided")
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    stop(paste0(group.by, " is not in the meta.data of srt object."))
  }
  if (length(nameslist) > 0 && is.null(names(nameslist))) {
    names(nameslist) <- levels(srt@meta.data[[group.by]])
  }
  if (is.list(nameslist) && length(nameslist) > 0) {
    names_assign <- setNames(rep(names(nameslist), sapply(nameslist, length)), nm = unlist(nameslist))
  } else {
    if (is.null(names(nameslist))) {
      if (!is.factor(srt@meta.data[[group.by]])) {
        stop("'nameslist' must be named when srt@meta.data[[group.by]] is not a factor")
      }
      if (!identical(length(nameslist), length(unique(srt@meta.data[[group.by]])))) {
        stop("'nameslist' must be named or the length of ", length(unique(srt@meta.data[[group.by]])))
      }
      names(nameslist) <- levels(srt@meta.data[[group.by]])
    }
    names_assign <- nameslist
  }
  if (all(!names(names_assign) %in% srt@meta.data[[group.by]])) {
    stop("No group name mapped.")
  }
  if (is.factor(srt@meta.data[[group.by]])) {
    levels <- levels(srt@meta.data[[group.by]])
  } else {
    levels <- NULL
  }
  index <- which(as.character(srt@meta.data[[group.by]]) %in% names(names_assign))
  srt@meta.data[[name]] <- as.character(srt@meta.data[[group.by]])
  srt@meta.data[[name]][index] <- names_assign[srt@meta.data[[name]][index]]
  if (!is.null(levels)) {
    levels[levels %in% names(names_assign)] <- names_assign[levels[levels %in% names(names_assign)]]
    if (isFALSE(keep_levels)) {
      levels <- unique(c(names_assign, levels))
    } else {
      levels <- unique(levels)
    }
    srt@meta.data[[name]] <- factor(srt@meta.data[[name]], levels = levels)
  }
  return(srt)
}

#' Reorder idents by the gene expression
#'
#' @param srt A Seurat object.
#' @param features Features used to reorder idents.
#' @param reorder_by Reorder groups instead of idents.
#' @param slot Specific slot to get data from.
#' @param assay Specific assay to get data from.
#' @param log Whether log1p transformation needs to be applied. Default is \code{TRUE}.
#' @param distance_metric Metric to compute distance. Default is "euclidean".
#'
#' @importFrom Seurat VariableFeatures DefaultAssay DefaultAssay<- AverageExpression Idents<-
#' @importFrom SeuratObject as.sparse
#' @importFrom stats hclust reorder as.dendrogram as.dist
#' @importFrom Matrix t colMeans
#' @export
SrtReorder <- function(srt, features = NULL, reorder_by = NULL, slot = "data", assay = NULL, log = TRUE,
                       distance_metric = "euclidean") {
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  features <- intersect(x = features, y = rownames(x = srt))
  if (is.null(reorder_by)) {
    srt$ident <- Idents(srt)
  } else {
    srt$ident <- srt[[reorder_by, drop = TRUE]]
  }
  if (length(unique(srt[[reorder_by, drop = TRUE]])) == 1) {
    warning("Only one cluter found.", immediate. = TRUE)
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

  data.avg <- AverageExpression(object = srt, features = features, slot = slot, assays = assay, group.by = "ident", verbose = FALSE)[[1]][features, , drop = FALSE]
  if (isTRUE(log)) {
    data.avg <- log1p(data.avg)
  }
  mat <- t(x = data.avg[features, , drop = FALSE])
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as.sparse(mat[1:nrow(mat), , drop = FALSE])
  }

  if (distance_metric %in% c(simil_method, "pearson", "spearman")) {
    if (distance_metric %in% c("pearson", "spearman")) {
      if (distance_metric == "spearman") {
        mat <- t(apply(mat, 1, rank))
      }
      distance_metric <- "correlation"
    }
    d <- 1 - proxyC::simil(as.sparse(mat[1:nrow(mat), , drop = FALSE]), method = distance_metric)
  } else if (distance_metric %in% dist_method) {
    d <- proxyC::dist(as.sparse(mat[1:nrow(mat), , drop = FALSE]), method = distance_metric)
  }
  data.dist <- as.dist(d)
  hc <- hclust(d = data.dist)
  dd <- as.dendrogram(hc)
  dd_ordered <- reorder(dd, wts = colMeans(data.avg[features, , drop = FALSE]), agglo.FUN = base::mean)
  ident_new <- unname(setNames(object = seq_along(labels(dd_ordered)), nm = labels(dd_ordered))[as.character(srt$ident)])
  ident_new <- factor(ident_new, levels = seq_along(labels(dd_ordered)))
  Idents(srt) <- srt$ident <- ident_new
  return(srt)
}

#' Append a Seurat object to another
#'
#' @param srt_raw A Seurat object to be appended.
#' @param srt_append New Seurat object to append.
#' @param slots slots names.
#' @param pattern A character string containing a regular expression. All data with matching names will be considered for appending.
#' @param overwrite Whether to overwrite.
#' @param verbose Show messages.
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
          } else if (identical(slot_nm, "assays")) {
            if (!info %in% Assays(srt_raw)) {
              srt_raw[[info]] <- srt_append[[info]]
            } else {
              srt_raw[[info]]@counts <- srt_append[[info]]@counts
              srt_raw[[info]]@data <- srt_append[[info]]@data
              srt_raw[[info]]@key <- srt_append[[info]]@key
              srt_raw[[info]]@var.features <- srt_append[[info]]@var.features
              srt_raw[[info]]@misc <- srt_append[[info]]@misc
              srt_raw[[info]]@meta.features <- cbind(srt_raw[[info]]@meta.features, srt_append[[info]]@meta.features[
                rownames(srt_raw[[info]]@meta.features),
                setdiff(colnames(srt_append[[info]]@meta.features), colnames(srt_raw[[info]]@meta.features))
              ])
            }
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

#' Run dimensionality reduction
#'
#' @param srt A Seurat object.
#' @param prefix The prefix used to name the result.
#' @param features Use features expression data to run linear or nonlinear dimensionality reduction.
#' @param assay Specific assay to get data from.
#' @param slot Specific slot to get data from.
#' @param linear_reduction Method of linear dimensionality reduction. Options are "pca", "ica", "nmf", "mds", "glmpca".
#' @param linear_reduction_dims Total number of dimensions to compute and store for \code{linear_reduction}.
#' @param linear_reduction_params Other parameters passed to the \code{linear_reduction} method.
#' @param force_linear_reduction Whether force to do linear dimensionality reduction.
#' @param nonlinear_reduction Method of nonlinear dimensionality reduction. Options are "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"
#' @param nonlinear_reduction_dims Total number of dimensions to compute and store for \code{nonlinear_reduction}.
#' @param reduction_use Which dimensional reduction to use as input for \code{nonlinear_reduction}.
#' @param reduction_dims Which dimensions to use as input for \code{nonlinear_reduction}, used only if \code{features} is \code{NULL}.
#' @param neighbor_use Name of neighbor to use for the \code{nonlinear_reduction}.
#' @param graph_use Name of graph to use for the \code{nonlinear_reduction}.
#' @param nonlinear_reduction_params  Other parameters passed to the \code{nonlinear_reduction} method.
#' @param force_nonlinear_reduction Whether force to do nonlinear dimensionality reduction.
#' @param verbose Show messages.
#' @param seed Set a seed.
#'
#' @importFrom Seurat Embeddings RunPCA RunICA RunTSNE Reductions DefaultAssay DefaultAssay<- Key Key<-
#' @importFrom Signac RunSVD
#' @export
RunDimReduction <- function(srt, prefix = "", features = NULL, assay = NULL, slot = "data",
                            linear_reduction = NULL, linear_reduction_dims = 50,
                            linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = NULL, nonlinear_reduction_dims = 2,
                            reduction_use = NULL, reduction_dims = NULL,
                            graph_use = NULL, neighbor_use = NULL,
                            nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            verbose = TRUE, seed = 11) {
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    type <- "Chromatin"
  } else {
    type <- "RNA"
  }
  linear_reduction_dims <- min(linear_reduction_dims, nrow(srt[[assay]]) - 1, ncol(srt[[assay]]) - 1, na.rm = TRUE)
  nonlinear_reduction_dims <- min(nonlinear_reduction_dims, nrow(srt[[assay]]) - 1, ncol(srt[[assay]]) - 1, na.rm = TRUE)
  if (!is.null(linear_reduction)) {
    if (any(!linear_reduction %in% c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt))) || length(linear_reduction) > 1) {
      stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
    }
  }
  if (!is.null(nonlinear_reduction)) {
    if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", Reductions(srt))) || length(nonlinear_reduction) > 1) {
      stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
    }
    if (is.null(features) && is.null(reduction_use) && is.null(neighbor_use) && is.null(graph_use)) {
      stop("'features', 'reduction_use', 'graph_use', or 'neighbor_use' must be provided when running non-linear dimensionality reduction.")
    }
    if (!is.null(features)) {
      message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Features(length:", length(features), ") as input")
    }
    if (!is.null(reduction_use)) {
      message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Reduction(", reduction_use, ", dims_range:", min(reduction_dims), "-", max(reduction_dims), ") as input")
    }
    if (!is.null(neighbor_use)) {
      message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Neighbors(", neighbor_use, ") as input")
    }
    if (!is.null(graph_use)) {
      message("Non-linear dimensionality reduction(", nonlinear_reduction, ") using Graphs(", graph_use, ") as input")
    }
  }
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
      if (length(VariableFeatures(srt, assay = assay)) == 0) {
        srt <- FindVariableFeatures(srt, assay = assay, verbose = FALSE)
      }
      features <- VariableFeatures(srt, assay = assay)
    }
    fun_use <- switch(linear_reduction,
      "pca" = "RunPCA",
      "svd" = "RunSVD",
      "ica" = "RunICA",
      "nmf" = "RunNMF",
      "mds" = "RunMDS",
      "glmpca" = "RunGLMPCA"
    )
    key_use <- switch(linear_reduction,
      "pca" = "PC_",
      "svd" = "LSI_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    components_nm <- switch(linear_reduction,
      "pca" = "npcs",
      "svd" = "n",
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
    if (fun_use %in% c("RunSVD", "RunICA")) {
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
      center <- rowMeans(GetAssayData(object = srt, slot = "scale.data", assay = assay)[features, , drop = FALSE])
      model <- list(sdev = pca.out@stdev, rotation = pca.out@feature.loadings, center = center, scale = FALSE, x = pca.out@cell.embeddings)
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["model"]] <- model
    }
    if (linear_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:linear_reduction_dims
    } else {
      dim_est <- tryCatch(expr = {
        min(
          intrinsicDimension::maxLikGlobalDimEst(data = Embeddings(srt, reduction = paste0(prefix, linear_reduction)), k = 20)[["dim.est"]],
          ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction)))
        )
      }, error = function(e) {
        message("Can not estimate intrinsic dimensions with maxLikGlobalDimEst.")
        return(NA)
      })
      if (!is.na(dim_est)) {
        dims_estimate <- seq_len(max(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 10), ceiling(dim_est)))
      } else {
        dims_estimate <- seq_len(min(ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))), 30))
      }
    }
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
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData SetAssayData VariableFeatures VariableFeatures<-
#' @export
Uncorrected_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                                  do_normalization = NULL, normalization_method = "LogNormalize",
                                  do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                                  do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                  linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                                  nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                  do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                  seed = 11) {
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
    check_Python("leidenalg")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Uncorrected) on the data...\n"))
  srtIntegrated <- Standard_SCP(
    srt = srtMerge, prefix = "Uncorrected", assay = DefaultAssay(srtMerge),
    do_normalization = do_normalization, normalization_method = normalization_method,
    do_HVF_finding = do_HVF_finding, nHVF = nHVF, HVF = HVF, HVF_method = HVF_method,
    do_scaling = do_scaling, vars_to_regress = vars_to_regress, regression_model = regression_model,
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_dims_use = linear_reduction_dims_use, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    nonlinear_reduction = nonlinear_reduction, nonlinear_reduction_dims = nonlinear_reduction_dims, nonlinear_reduction_params = nonlinear_reduction_params, force_nonlinear_reduction = force_nonlinear_reduction,
    do_cluster_finding = do_cluster_finding, cluster_algorithm = cluster_algorithm, cluster_resolution = cluster_resolution, cluster_reorder = cluster_reorder,
    seed = seed
  )
  if ("Uncorrectedclusters" %in% colnames(srtIntegrated@meta.data)) {
    srtIntegrated[[paste0("Uncorrectedclusters")]] <- srtIntegrated[[paste0("Uncorrected", linear_reduction, "clusters")]]
    srtIntegrated[[paste0("Uncorrected", linear_reduction, "clusters")]] <- NULL
  }
  for (nr in nonlinear_reduction) {
    for (n in nonlinear_reduction_dims) {
      if (paste0("Uncorrected", linear_reduction, toupper(nr), n, "D") %in% Reductions(srtIntegrated)) {
        reduc <- srtIntegrated@reductions[[paste0("Uncorrected", linear_reduction, toupper(nr), n, "D")]]
        Key(reduc) <- paste0("Uncorrected", toupper(nr), n, "D_")
        srtIntegrated@reductions[[paste0("Uncorrected", linear_reduction, toupper(nr), n, "D")]] <- NULL
        srtIntegrated@reductions[[paste0("Uncorrected", toupper(nr), n, "D")]] <- reduc
      }
    }
    srtIntegrated@misc[["Default_reduction"]] <- paste0("Uncorrected", toupper(nr))
  }
  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Uncorrected_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Uncorrected|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Seurat_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData FindIntegrationAnchors IntegrateData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents
#' @importFrom Signac RunTFIDF
#' @export
Seurat_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                             linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                             FindIntegrationAnchors_params = list(), IntegrateData_params = list(), IntegrateEmbeddings_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "svd", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  set.seed(seed)
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    if (normalization_method == "TFIDF") {
      srtMerge <- Reduce(merge, srtList)
      VariableFeatures(srtMerge) <- HVF
    }
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'rlsi' integration workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
    FindIntegrationAnchors_params[["reduction"]] <- "rlsi"
    if (is.null(FindIntegrationAnchors_params[["dims"]])) {
      FindIntegrationAnchors_params[["dims"]] <- 2:min(linear_reduction_dims, 30)
    }
    srtMerge <- suppressWarnings(RunTFIDF(object = srtMerge, assay = DefaultAssay(srtMerge), verbose = FALSE))
    srtMerge <- RunDimReduction(
      srt = srtMerge, prefix = "", features = HVF, assay = DefaultAssay(srtMerge),
      linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    srtMerge[["lsi"]] <- srtMerge[["svd"]]
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (svd) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt = srt, prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "svd", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srt[["lsi"]] <- srt[["svd"]]
      srtList[[i]] <- srt
    }
  }

  if (isTRUE(FindIntegrationAnchors_params[["reduction"]] == "rpca")) {
    cat(paste0("[", Sys.time(), "]", " Use 'rpca' integration workflow...\n"))
    for (i in seq_along(srtList)) {
      srt <- srtList[[i]]
      if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data", assay = DefaultAssay(srt)))))) {
        cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
        srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
      }
      cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (pca) on the data ", i, " ...\n"))
      srt <- RunDimReduction(
        srt = srt, prefix = "", features = HVF, assay = DefaultAssay(srt),
        linear_reduction = "pca", linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
        verbose = FALSE, seed = seed
      )
      srtList[[i]] <- srt
    }
  }

  if (is.null(names(srtList))) {
    names(srtList) <- paste0("srt_", seq_along(srtList))
  }

  if (normalization_method %in% c("LogNormalize", "SCT")) {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = normalization_method,
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
      normalization.method = normalization_method,
      features.to.integrate = HVF,
      verbose = FALSE
    )
    for (nm in names(IntegrateData_params)) {
      params2[[nm]] <- IntegrateData_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateData, .args = params2)

    DefaultAssay(srtIntegrated) <- "Seurat"
    VariableFeatures(srtIntegrated[["Seurat"]]) <- HVF

    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srtIntegrated <- ScaleData(object = srtIntegrated, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }

    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
    srtIntegrated <- RunDimReduction(
      srt = srtIntegrated, prefix = "Seurat", features = HVF, assay = DefaultAssay(srtIntegrated),
      linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 1:linear_reduction_dims
    }
  } else if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " Perform FindIntegrationAnchors on the data...\n"))
    params1 <- list(
      object.list = srtList,
      normalization.method = "LogNormalize",
      anchor.features = HVF,
      reduction = "rlsi",
      verbose = FALSE
    )
    for (nm in names(FindIntegrationAnchors_params)) {
      params1[[nm]] <- FindIntegrationAnchors_params[[nm]]
    }
    suppressWarnings({
      srt_anchors <- invoke(.fn = FindIntegrationAnchors, .args = params1)
    })

    cat(paste0("[", Sys.time(), "]", " Perform integration(Seurat) on the data...\n"))
    params2 <- list(
      anchorset = srt_anchors,
      reductions = srtMerge[["lsi"]],
      new.reduction.name = "Seuratlsi",
      verbose = FALSE
    )
    for (nm in names(IntegrateEmbeddings_params)) {
      params2[[nm]] <- IntegrateEmbeddings_params[[nm]]
    }
    srtIntegrated <- invoke(.fn = IntegrateEmbeddings, .args = params2)

    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use <- 2:linear_reduction_dims
    }
    linear_reduction <- "lsi"
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Seurat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Seurat|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' scVI_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom reticulate import
#' @export
scVI_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                           do_normalization = NULL, normalization_method = "LogNormalize",
                           do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                           scVI_dims_use = NULL,
                           nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                           do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                           model = "SCVI", SCVI_params = list(), PEAKVI_params = list(), num_threads = 8, seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  if (.Platform$OS.type == "windows" && !exist_Python_pkgs(packages = "scvi-tools")) {
    suppressWarnings(system2(command = conda_python(), args = "-m pip install jax[cpu]===0.3.20 -f https://whls.blob.core.windows.net/unstable/index.html --use-deprecated legacy-resolver", stdout = TRUE))
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
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  adata <- srt_to_adata(srtMerge, features = HVF, assay_X = DefaultAssay(srtMerge), assay_layers = NULL, verbose = FALSE)
  adata[["X"]] <- scipy$sparse$csr_matrix(adata[["X"]])

  if (model == "SCVI") {
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
    DefaultAssay(srtIntegrated) <- "scVIcorrected"
    VariableFeatures(srtIntegrated[["scVIcorrected"]]) <- HVF
  } else if (model == "PEAKVI") {
    message("Assay is ChromatinAssay. Using PeakVI workflow.")
    scvi$model$PEAKVI$setup_anndata(adata, batch_key = batch)
    params <- list(
      adata = adata
    )
    for (nm in names(PEAKVI_params)) {
      params[[nm]] <- PEAKVI_params[[nm]]
    }
    model <- invoke(.fn = scvi$model$PEAKVI, .args = params)
    model$train()
    srtIntegrated <- srtMerge
    srtMerge <- NULL
  }

  latent <- as.matrix(model$get_latent_representation())
  rownames(latent) <- colnames(srtIntegrated)
  colnames(latent) <- paste0("scvi_", seq_len(ncol(latent)))
  srtIntegrated[["scVI"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(srtIntegrated))
  if (is.null(scVI_dims_use)) {
    scVI_dims_use <- 1:ncol(srtIntegrated[["scVI"]]@cell.embeddings)
  }

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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["scVI_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|scVI|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' MNN_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
MNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          mnnCorrect_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  sceList <- lapply(srtList, function(srt) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = GetAssayData(srt, slot = "data", assay = DefaultAssay(srt))[HVF, , drop = FALSE]))
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- as_matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(MNN) on the data...\n"))
  params <- list(
    sceList,
    cos.norm.out = FALSE
  )
  for (nm in names(mnnCorrect_params)) {
    params[[nm]] <- mnnCorrect_params[[nm]]
  }
  out <- invoke(.fn = batchelor::mnnCorrect, .args = params)

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["MNNcorrected"]] <- CreateAssayObject(counts = out@assays@data$corrected)
  VariableFeatures(srtIntegrated[["MNNcorrected"]]) <- HVF
  DefaultAssay(srtIntegrated) <- "MNNcorrected"

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(object = srtIntegrated, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "MNN", features = HVF, assay = DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["MNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|MNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' fastMNN_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat CreateSeuratObject GetAssayData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
fastMNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                              do_normalization = NULL, normalization_method = "LogNormalize",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                              fastMNN_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              fastMNN_params = list(), seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("batchelor")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  sceList <- lapply(srtList, function(srt) {
    sce <- as.SingleCellExperiment(CreateSeuratObject(counts = GetAssayData(srt, slot = "data", assay = DefaultAssay(srt))[HVF, , drop = FALSE]))
    if (inherits(sce@assays@data$logcounts, "dgCMatrix")) {
      sce@assays@data$logcounts <- as_matrix(sce@assays@data$logcounts)
    }
    return(sce)
  })
  if (is.null(names(sceList))) {
    names(sceList) <- paste0("sce_", seq_along(sceList))
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(fastMNN) on the data...\n"))
  params <- list(
    sceList
  )
  for (nm in names(fastMNN_params)) {
    params[[nm]] <- fastMNN_params[[nm]]
  }
  out <- invoke(.fn = batchelor::fastMNN, .args = params)

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["fastMNNcorrected"]] <- CreateAssayObject(counts = as.matrix(out@assays@data$reconstructed))
  DefaultAssay(srtIntegrated) <- "fastMNNcorrected"
  VariableFeatures(srtIntegrated[["fastMNNcorrected"]]) <- HVF
  reduction <- out@int_colData$reducedDims$corrected
  colnames(reduction) <- paste0("fastMNN_", seq_len(ncol(reduction)))
  srtIntegrated[["fastMNN"]] <- CreateDimReducObject(embeddings = reduction, key = "fastMNN_", assay = "fastMNNcorrected")

  if (is.null(fastMNN_dims_use)) {
    fastMNN_dims_use <- 1:ncol(srtIntegrated[["fastMNN"]]@cell.embeddings)
  }

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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["fastMNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|fastMNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Harmony_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
Harmony_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                              do_normalization = NULL, normalization_method = "LogNormalize",
                              do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                              do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                              linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                              Harmony_dims_use = NULL,
                              nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                              do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                              RunHarmony_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("harmony")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "Harmony", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Harmony) on the data...\n"))
  message("Harmony integration using Reduction(", paste0("Harmony", linear_reduction), ", dims_range:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
  params <- list(
    object = srtMerge,
    group.by.vars = batch,
    assay.use = DefaultAssay(srtMerge),
    reduction = paste0("Harmony", linear_reduction),
    dims.use = linear_reduction_dims_use,
    reduction.save = "Harmony",
    verbose = FALSE
  )
  if (nrow(GetAssayData(srtMerge, slot = "scale.data", assay = DefaultAssay(srtMerge))) == 0) {
    params[["project.dim"]] <- FALSE
  }
  for (nm in names(RunHarmony_params)) {
    params[[nm]] <- RunHarmony_params[[nm]]
  }
  srtIntegrated <- invoke(.fn = RunHarmony2, .args = params)

  if (is.null(Harmony_dims_use)) {
    Harmony_dims_use <- 1:ncol(srtIntegrated[["Harmony"]]@cell.embeddings)
  }

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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Harmony_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Harmony|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Scanorama_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @importFrom Matrix t
#' @importFrom reticulate import
#' @importFrom stats sd
#' @export
Scanorama_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                                do_normalization = NULL, normalization_method = "LogNormalize",
                                do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                                do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                                Scanorama_dims_use = NULL,
                                linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                                nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                                do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                                return_corrected = FALSE, scanorama_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_Python("scanorama")
  scanorama <- import("scanorama")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    srtList <- SplitObject(object = srtMerge, split.by = batch)
    checked <- check_srtList(
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  srtIntegrated <- Reduce(merge, srtList)

  cat(paste0("[", Sys.time(), "]", " Perform integration(Scanorama) on the data...\n"))
  assaylist <- list()
  genelist <- list()
  for (i in seq_along(srtList)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(object = srtList[[i]], slot = "data", assay = DefaultAssay(srtList[[i]]))[HVF, , drop = FALSE]))
    genelist[[i]] <- HVF
  }
  if (isTRUE(return_corrected)) {
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
    corrected <- invoke(.fn = scanorama$correct, .args = params)

    cor_value <- t(do.call(rbind, corrected[[2]]))
    rownames(cor_value) <- corrected[[3]]
    colnames(cor_value) <- unlist(sapply(assaylist, rownames))
    srtIntegrated[["Scanoramacorrected"]] <- CreateAssayObject(data = cor_value)
    VariableFeatures(srtIntegrated[["Scanoramacorrected"]]) <- HVF

    dim_reduction <- do.call(rbind, corrected[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0("Scanorama_", seq_len(ncol(dim_reduction)))
  } else {
    params <- list(
      datasets_full = assaylist,
      genes_list = genelist,
      verbose = FALSE
    )
    for (nm in names(scanorama_params)) {
      params[[nm]] <- scanorama_params[[nm]]
    }
    integrated <- invoke(.fn = scanorama$integrate, .args = params)

    dim_reduction <- do.call(rbind, integrated[[1]])
    rownames(dim_reduction) <- unlist(sapply(assaylist, rownames))
    colnames(dim_reduction) <- paste0("Scanorama_", seq_len(ncol(dim_reduction)))
  }
  srtIntegrated[["Scanorama"]] <- CreateDimReducObject(embeddings = dim_reduction, key = "Scanorama_", assay = DefaultAssay(srtIntegrated))

  if (is.null(Scanorama_dims_use)) {
    Scanorama_dims_use <- 1:ncol(srtIntegrated[["Scanorama"]]@cell.embeddings)
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = "Scanorama", dims = Scanorama_dims_use, force.recalc = TRUE, graph.name = paste0("Scanorama_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = "Scanorama_SNN", verbose = FALSE)
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
            reduction_use = "Scanorama", reduction_dims = Scanorama_dims_use,
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Scanorama_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Scanorama|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' BBKNN_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- as.Graph Embeddings FindClusters Idents VariableFeatures VariableFeatures<- as.sparse
#' @importFrom Matrix t
#' @importFrom reticulate import
#' @export
BBKNN_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            bbknn_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_Python("bbknn")
  bbknn <- import("bbknn")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "BBKNN", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(BBKNN) on the data...\n"))
  message("BBKNN integration using Reduction(", paste0("BBKNN", linear_reduction), ", dims_range:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
  emb <- Embeddings(srtMerge, reduction = paste0("BBKNN", linear_reduction))[, linear_reduction_dims_use, drop = FALSE]
  params <- list(
    pca = emb,
    batch_list = srtMerge[[batch, drop = TRUE]]
  )
  for (nm in names(bbknn_params)) {
    params[[nm]] <- bbknn_params[[nm]]
  }
  bem <- invoke(.fn = bbknn$bbknn_matrix, .args = params)
  n.neighbors <- bem[[3]]$n_neighbors
  srtIntegrated <- srtMerge

  bbknn_graph <- as.sparse(bem[[2]][1:nrow(bem[[2]]), , drop = FALSE])
  rownames(bbknn_graph) <- colnames(bbknn_graph) <- rownames(emb)
  bbknn_graph <- as.Graph(bbknn_graph)
  bbknn_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["BBKNN"]] <- bbknn_graph

  bbknn_dist <- t(as.sparse(bem[[1]][1:nrow(bem[[1]]), , drop = FALSE]))
  rownames(bbknn_dist) <- colnames(bbknn_dist) <- rownames(emb)
  val <- split(bbknn_dist@x, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  pos <- split(bbknn_dist@i + 1, rep(1:ncol(bbknn_dist), diff(bbknn_dist@p)))
  idx <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = F), n.neighbors)]
    length(out) <- n.neighbors
    return(out)
  }, x = val, y = pos))
  idx[is.na(idx)] <- seq_len(nrow(idx))
  dist <- t(mapply(function(x, y) {
    out <- y[head(order(x, decreasing = F), n.neighbors)]
    length(out) <- n.neighbors
    out[is.na(out)] <- 0
    return(out)
  }, x = val, y = val))
  srtIntegrated[["BBKNN_neighbors"]] <- new(Class = "Neighbor", nn.idx = idx, nn.dist = dist, alg.info = list(), cell.names = rownames(emb))
  nonlinear_reduction_params[["n.neighbors"]] <- n.neighbors

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
            srt = srtIntegrated, prefix = "BBKNN", neighbor_use = "BBKNN_neighbors",
            nonlinear_reduction = nonlinear_reduction, nonlinear_reduction_dims = n,
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["BBKNN_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|BBKNN|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' CSS_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
CSS_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                          do_normalization = NULL, normalization_method = "LogNormalize",
                          do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                          do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                          linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                          CSS_dims_use = NULL,
                          nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                          do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                          CSS_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("quadbiolab/simspec")
  suppressPackageStartupMessages(requireNamespace("qlcMatrix"))
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtMerge, slot = "scale.data", assay = DefaultAssay(srtMerge)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtMerge <- ScaleData(object = srtMerge, assay = DefaultAssay(srtMerge), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtMerge <- RunDimReduction(
    srt = srtMerge, prefix = "CSS", features = HVF, assay = DefaultAssay(srtMerge),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(CSS) on the data...\n"))
  message("CSS integration using Reduction(", paste0("CSS", linear_reduction), ", dims_range:", min(linear_reduction_dims_use), "-", max(linear_reduction_dims_use), ") as input")
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

  if (any(is.na(srtIntegrated@reductions[["CSS"]]@cell.embeddings))) {
    stop("NA detected in the CSS embeddings. You can try to use a lower resolution value in the CSS_param.")
  }
  if (is.null(CSS_dims_use)) {
    CSS_dims_use <- 1:ncol(srtIntegrated[["CSS"]]@cell.embeddings)
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["CSS_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|CSS|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' LIGER_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @export
LIGER_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            LIGER_dims_use = NULL,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            optimizeALS_params = list(), quantilenorm_params = list(), seed = 11) {
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("rliger")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  scale.data <- list()
  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data", assay = DefaultAssay(srt)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, do.center = FALSE, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    scale.data[[i]] <- t(x = GetAssayData(object = srt, slot = "scale.data", assay = DefaultAssay(srt)))
  }

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
  if (is.null(LIGER_dims_use)) {
    LIGER_dims_use <- 1:ncol(srtIntegrated[["LIGER"]]@cell.embeddings)
  }

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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["LIGER_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|LIGER|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Conos_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- Embeddings FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom igraph as_adjacency_matrix
#' @export
Conos_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            buildGraph_params = list(), num_threads = 2, seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'umap-naive'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("conos")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }
  srtIntegrated <- srtMerge
  srtMerge <- NULL

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  for (i in seq_along(srtList)) {
    srt <- srtList[[i]]
    if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data", assay = DefaultAssay(srt)))))) {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data ", i, " ...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data ", i, " ...\n"))
    srt <- RunDimReduction(
      srt = srt, prefix = "", features = HVF, assay = DefaultAssay(srt),
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
    maxdims <- linear_reduction_dims
  } else {
    maxdims <- max(linear_reduction_dims_use)
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Conos) on the data...\n"))
  message("Conos using ", linear_reduction, "(dims_max:", maxdims, ") as input")
  srtList_con <- conos::Conos$new(srtList, n.cores = num_threads)
  params <- list(
    ncomps = maxdims,
    k = 15,
    verbose = FALSE
  )
  for (nm in names(buildGraph_params)) {
    params[[nm]] <- buildGraph_params[[nm]]
  }
  invoke(.fn = srtList_con[["buildGraph"]], .args = params)
  conos_graph <- as_adjacency_matrix(srtList_con$graph, type = "both", attr = "weight", names = TRUE, sparse = TRUE)
  conos_graph <- as.Graph(conos_graph)
  conos_graph@assay.used <- DefaultAssay(srtIntegrated)
  srtIntegrated@graphs[["Conos"]] <- conos_graph
  nonlinear_reduction_params[["n.neighbors"]] <- params[["k"]]

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
            srt = srtIntegrated, prefix = "Conos", graph_use = "Conos",
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["Conos_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|Conos|Default_reduction"), overwrite = TRUE, verbose = FALSE)
    return(srtMerge_raw)
  } else {
    return(srtIntegrated)
  }
}

#' Combat_integrate
#'
#' @inheritParams Integration_SCP
#'
#' @importFrom Seurat GetAssayData ScaleData SetAssayData DefaultAssay DefaultAssay<- SplitObject CreateAssayObject CreateDimReducObject Embeddings FindNeighbors FindClusters Idents
#' @export
ComBat_integrate <- function(srtMerge = NULL, batch = NULL, append = TRUE, srtList = NULL, assay = NULL,
                             do_normalization = NULL, normalization_method = "LogNormalize",
                             do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                             do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                             linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                             nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                             do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                             ComBat_params = list(), seed = 11) {
  if (length(linear_reduction) > 1) {
    warning("Only the first method in the 'linear_reduction' will be used.", immediate. = TRUE)
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srtMerge)) {
    reduc_test <- c(reduc_test, Reductions(srtMerge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  check_R("sva")
  set.seed(seed)

  if (is.null(srtList) && is.null(srtMerge)) {
    stop("srtList and srtMerge were all empty.")
  }
  if (!is.null(srtList) && !is.null(srtMerge)) {
    cell1 <- sort(unique(unlist(lapply(srtList, colnames))))
    cell2 <- sort(unique(colnames(srtMerge)))
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
      srtList = srtList, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtList <- checked[["srtList"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srtMerge <- Reduce(merge, srtList)
    VariableFeatures(srtMerge) <- HVF
  }
  if (is.null(srtList) && !is.null(srtMerge)) {
    checked <- check_srtMerge(
      srtMerge = srtMerge, batch = batch, assay = assay,
      do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source, HVF_method = HVF_method,
      nHVF = nHVF, HVF_min_intersection = HVF_min_intersection, HVF = HVF,
      vars_to_regress = vars_to_regress, seed = seed
    )
    srtMerge <- checked[["srtMerge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  cat(paste0("[", Sys.time(), "]", " Perform integration(Combat) on the data...\n"))
  dat <- GetAssayData(srtMerge, slot = "data", assay = DefaultAssay(srtMerge))[HVF, , drop = FALSE]
  batch <- srtMerge[[batch, drop = TRUE]]
  params <- list(
    dat = dat,
    batch = batch
  )
  for (nm in names(ComBat_params)) {
    params[[nm]] <- ComBat_params[[nm]]
  }
  corrected <- suppressMessages(invoke(.fn = sva::ComBat, .args = params))

  srtIntegrated <- srtMerge
  srtMerge <- NULL
  srtIntegrated[["ComBatcorrected"]] <- CreateAssayObject(data = corrected)
  DefaultAssay(srtIntegrated) <- "ComBatcorrected"
  VariableFeatures(srtIntegrated[["ComBatcorrected"]]) <- HVF

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srtIntegrated, slot = "scale.data", assay = DefaultAssay(srtIntegrated)))))) {
    cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
    srtIntegrated <- ScaleData(srtIntegrated, assay = DefaultAssay(srtIntegrated), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
  }

  cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", linear_reduction, ") on the data...\n"))
  srtIntegrated <- RunDimReduction(
    srt = srtIntegrated, prefix = "ComBat", features = HVF, assay = DefaultAssay(srtIntegrated),
    linear_reduction = linear_reduction, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
    verbose = FALSE, seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- 1:linear_reduction_dims
    if (normalization_method == "TFIDF") {
      linear_reduction_dims_use <- 2:max(linear_reduction_dims_use)
    }
  }

  srtIntegrated <- tryCatch(
    {
      srtIntegrated <- FindNeighbors(object = srtIntegrated, reduction = paste0("ComBat", linear_reduction), dims = linear_reduction_dims_use, force.recalc = TRUE, graph.name = paste0("ComBat", linear_reduction, "_", c("KNN", "SNN")), verbose = FALSE)
      if (isTRUE(do_cluster_finding)) {
        cat(paste0("[", Sys.time(), "]", " Perform FindClusters (", cluster_algorithm, ") on the data...\n"))
        srtIntegrated <- FindClusters(object = srtIntegrated, resolution = cluster_resolution, algorithm = cluster_algorithm_index, method = "igraph", graph.name = paste0("ComBat", linear_reduction, "_SNN"), verbose = FALSE)
        if (isTRUE(cluster_reorder)) {
          cat(paste0("[", Sys.time(), "]", " Reorder clusters...\n"))
          srtIntegrated <- SrtReorder(srtIntegrated, features = HVF, reorder_by = "seurat_clusters", slot = "data")
        }
        srtIntegrated[["seurat_clusters"]] <- NULL
        srtIntegrated[["ComBatclusters"]] <- Idents(srtIntegrated)
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
            srt = srtIntegrated, prefix = "ComBat",
            reduction_use = paste0("ComBat", linear_reduction), reduction_dims = linear_reduction_dims_use,
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

  DefaultAssay(srtIntegrated) <- assay
  VariableFeatures(srtIntegrated) <- srtIntegrated@misc[["ComBat_HVF"]] <- HVF

  if (isTRUE(append) && !is.null(srtMerge_raw)) {
    srtMerge_raw <- SrtAppend(srt_raw = srtMerge_raw, srt_append = srtIntegrated, pattern = paste0(assay, "|ComBat|Default_reduction"), overwrite = TRUE, verbose = FALSE)
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
#' @return A \code{Seurat} object.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(srt = pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType")
#'
#' # Use a combination of different linear or non-linear dimension reduction methods
#' linear_reductions <- c("pca", "ica", "nmf", "mds", "glmpca")
#' pancreas_sub <- Standard_SCP(
#'   srt = pancreas_sub,
#'   linear_reduction = linear_reductions,
#'   nonlinear_reduction = "umap"
#' )
#' plist1 <- lapply(linear_reductions, function(lr) CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = paste0("Standard", lr, "UMAP2D"), theme = "theme_blank"))
#' patchwork::wrap_plots(plotlist = plist1)
#'
#' nonlinear_reductions <- c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")
#' pancreas_sub <- Standard_SCP(
#'   srt = pancreas_sub,
#'   linear_reduction = "pca",
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' plist2 <- lapply(nonlinear_reductions, function(nr) CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = paste0("Standardpca", toupper(nr), "2D"), theme = "theme_blank"))
#' patchwork::wrap_plots(plotlist = plist2)
#'
#' @importFrom Seurat Assays GetAssayData NormalizeData SCTransform SCTResults ScaleData SetAssayData DefaultAssay DefaultAssay<- FindNeighbors FindClusters Idents VariableFeatures VariableFeatures<-
#' @importFrom Matrix rowSums
#' @export
Standard_SCP <- function(srt, prefix = "Standard", assay = NULL,
                         do_normalization = NULL, normalization_method = "LogNormalize",
                         do_HVF_finding = TRUE, HVF_method = "vst", nHVF = 2000, HVF = NULL,
                         do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                         linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                         nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                         do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                         seed = 11) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  if (any(!linear_reduction %in% c("pca", "svd", "ica", "nmf", "mds", "glmpca", Reductions(srt)))) {
    stop("'linear_reduction' must be one of 'pca', 'svd', 'ica', 'nmf', 'mds', 'glmpca'.")
  }
  if (!is.null(linear_reduction_dims_use) && max(linear_reduction_dims_use) > linear_reduction_dims) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }
  if (any(!nonlinear_reduction %in% c("umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"))) {
    stop("'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis'.")
  }
  if (!cluster_algorithm %in% c("louvain", "slm", "leiden")) {
    stop("'cluster_algorithm' must be one of 'louvain', 'slm', 'leiden'.")
  }
  if (cluster_algorithm == "leiden") {
    check_Python("leidenalg")
  }
  cluster_algorithm_index <- switch(tolower(cluster_algorithm),
    "louvain" = 1,
    "louvain_refined" = 2,
    "slm" = 3,
    "leiden" = 4
  )

  time_start <- Sys.time()
  set.seed(seed)

  cat(paste0("[", time_start, "] ", "Start Standard_SCP\n"))

  checked <- check_srtList(
    srtList = list(srt), batch = "", assay = assay,
    do_normalization = do_normalization, do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = "separate", HVF_method = HVF_method, nHVF = nHVF, HVF = HVF,
    vars_to_regress = vars_to_regress, seed = seed
  )
  srt <- checked[["srtList"]][[1]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]

  if (normalization_method == "TFIDF") {
    cat(paste0("[", Sys.time(), "]", " normalization_method is 'TFIDF'. Use 'lsi' workflow...\n"))
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  if (isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% rownames(GetAssayData(srt, slot = "scale.data", assay = DefaultAssay(srt)))))) {
    if (normalization_method != "SCT") {
      cat(paste0("[", Sys.time(), "]", " Perform ScaleData on the data...\n"))
      srt <- ScaleData(object = srt, assay = DefaultAssay(srt), features = HVF, vars.to.regress = vars_to_regress, model.use = regression_model, verbose = FALSE)
    }
  }

  for (lr in linear_reduction) {
    cat(paste0("[", Sys.time(), "]", " Perform linear dimension reduction (", lr, ") on the data...\n"))
    srt <- RunDimReduction(
      srt = srt, prefix = prefix, features = HVF, assay = DefaultAssay(srt),
      linear_reduction = lr, linear_reduction_dims = linear_reduction_dims, linear_reduction_params = linear_reduction_params, force_linear_reduction = force_linear_reduction,
      verbose = FALSE, seed = seed
    )
    if (is.null(linear_reduction_dims_use)) {
      linear_reduction_dims_use_current <- 1:linear_reduction_dims
      if (normalization_method == "TFIDF") {
        linear_reduction_dims_use_current <- 2:max(linear_reduction_dims_use_current)
      }
    } else {
      linear_reduction_dims_use_current <- linear_reduction_dims_use
    }

    srt <- tryCatch(
      {
        srt <- FindNeighbors(object = srt, reduction = paste0(prefix, lr), dims = linear_reduction_dims_use_current, force.recalc = TRUE, graph.name = paste0(prefix, lr, "_", c("KNN", "SNN")), verbose = FALSE)
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
              reduction_use = paste0(prefix, lr), reduction_dims = linear_reduction_dims_use_current,
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

  DefaultAssay(srt) <- assay
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
#' @param srtMerge A merged \code{Seurat} object with batch information.
#' @param batch Metadata column name containing the batch information.
#' @param append Whether append results into the \code{srtMerge}. Only valid when srtMerge is provided.
#' @param srtList A list of \code{Seurat} object.
#' @param assay
#' @param integration_method Integration method. Can be one of "Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat".
#' @param do_normalization Whether to normalize the data. If NULL, will automatically determine.
#' @param normalization_method Normalization method.Can be one of "LogNormalize", "SCT".
#' @param do_HVF_finding Whether to find the high variable features(HVF). If NULL, will automatically determine.
#' @param nHVF HVF number to use.
#' @param HVF_source Source of the HVF. Can be one of "separate" and "global".
#' @param HVF_method
#' @param HVF Custom high variable features.
#' @param HVF_min_intersection
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
#' @param force_linear_reduction
#' @param do_cluster_finding
#' @param linear_reduction_params
#' @param nonlinear_reduction_params
#' @param force_nonlinear_reduction
#' @param ...
#'
#' @return A \code{Seurat} object.
#'
#' @examples
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech", integration_method = "Uncorrected"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' \dontrun{
#' integration_methods <- c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat")
#' for (method in integration_methods) {
#'   panc8_sub <- Integration_SCP(
#'     srtMerge = panc8_sub, batch = "tech",
#'     integration_method = method, nonlinear_reduction = "umap"
#'   )
#'   print(CellDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = paste0(method, "UMAP2D"), theme_use = "theme_blank"))
#' }
#'
#' nonlinear_reductions <- c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")
#' panc8_sub <- Integration_SCP(
#'   srtMerge = panc8_sub, batch = "tech", integration_method = "Seurat",
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' for (nr in nonlinear_reductions) {
#'   print(CellDimPlot(panc8_sub, group.by = c("tech", "celltype"), reduction = paste0("Seurat", nr, "2D"), theme_use = "theme_blank"))
#' }
#' }
#'
#' @export
Integration_SCP <- function(srtMerge = NULL, batch, append = TRUE, srtList = NULL, assay = NULL,
                            integration_method = "Uncorrected",
                            do_normalization = NULL, normalization_method = "LogNormalize",
                            do_HVF_finding = TRUE, HVF_source = "separate", HVF_method = "vst", nHVF = 2000, HVF_min_intersection = 1, HVF = NULL,
                            do_scaling = TRUE, vars_to_regress = NULL, regression_model = "linear",
                            linear_reduction = "pca", linear_reduction_dims = 50, linear_reduction_dims_use = NULL, linear_reduction_params = list(), force_linear_reduction = FALSE,
                            nonlinear_reduction = "umap", nonlinear_reduction_dims = c(2, 3), nonlinear_reduction_params = list(), force_nonlinear_reduction = TRUE,
                            do_cluster_finding = TRUE, cluster_algorithm = "louvain", cluster_resolution = 0.6, cluster_reorder = TRUE,
                            seed = 11, ...) {
  if (is.null(srtList) && is.null(srtMerge)) {
    stop("Neither 'srtList' nor 'srtMerge' was found.")
  }
  if (length(integration_method) == 1 && integration_method %in% c("Uncorrected", "Seurat", "scVI", "MNN", "fastMNN", "Harmony", "Scanorama", "BBKNN", "CSS", "LIGER", "Conos", "ComBat")) {
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
    stop(paste(integration_method, "is not a supported integration method!"))
  }
}
