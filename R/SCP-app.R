#' CreateDataFile
#'
#' Creates a data file in HDF5 format from a Seurat object.
#'
#' @param srt The Seurat object.
#' @param DataFile Path to the output data file. If not provided, the file will be named "Data.hdf5" in the current directory.
#' @param name Name of the dataset. If not provided, the name will default to the Seurat object's project name.
#' @param assays Character vector specifying the assays to include in the data file. Default is "RNA".
#' @param slots Character vector specifying the slots to include in the data file. Default is "data".
#' @param compression_level Compression level for the HDF5 dataset. Default is 6.
#' @param overwrite Logical value indicating whether to overwrite existing data in the data file. Default is FALSE.
#'
#' @seealso \code{\link{CreateMetaFile}} \code{\link{PrepareSCExplorer}} \code{\link{FetchH5}} \code{\link{RunSCExplorer}}
#'
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratObject as.sparse
#' @importFrom rhdf5 h5ls h5createFile h5createGroup h5delete h5write
#' @importFrom HDF5Array writeTENxMatrix
#' @export
CreateDataFile <- function(srt, DataFile, name = NULL, assays = "RNA", slots = "data", compression_level = 6, overwrite = FALSE) {
  if (missing(DataFile) || is.null(DataFile)) {
    DataFile <- "Data.hdf5"
  }
  if (!file.exists(DataFile)) {
    h5createFile(DataFile)
  }
  if (is.null(name)) {
    name <- srt@project.name
    message("Set the dataset name to ", name)
  }
  if (substr(name, 1, 1) != "/") {
    name <- paste0("/", name)
  }
  if (!name %in% h5ls(DataFile)$group) {
    h5createGroup(file = DataFile, group = name)
  }

  message("Write the expression matrix to hdf5 file: ", DataFile)
  for (assay in assays) {
    for (slot in slots) {
      data <- t(GetAssayData(srt, slot = slot, assay = assay))
      if (isTRUE(overwrite)) {
        try(h5delete(file = DataFile, name = paste0(name, "/", assay, "/", slot)), silent = TRUE)
      }
      if (paste0(name, "/", assay, "/", slot) %in% h5ls(DataFile)$group) {
        message("Group ", paste0(name, "/", assay, "/", slot), " already exists in the ", DataFile)
      } else {
        if (!paste0(name, "/", assay) %in% h5ls(DataFile)$group) {
          h5createGroup(file = DataFile, group = paste0(name, "/", assay))
        }
        if (!inherits(data, "dgCMatrix")) {
          data <- as.sparse(data[1:nrow(data), ])
        }
        writeTENxMatrix(x = data, filepath = DataFile, group = paste0(name, "/", assay, "/", slot), level = compression_level)
      }
    }
  }
  if (isTRUE(overwrite)) {
    try(h5delete(file = DataFile, name = paste0(name, "/Default_assay")), silent = TRUE)
  }
  if (paste0(name, "/Default_assay") %in% h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/Default_assay"), " already exists in the ", DataFile)
  } else {
    h5write(obj = DefaultAssay(srt), file = DataFile, name = paste0(name, "/Default_assay"), level = compression_level)
  }
  if (isTRUE(overwrite)) {
    try(h5delete(file = DataFile, name = paste0(name, "/cells")), silent = TRUE)
  }
  if (paste0(name, "/cells") %in% h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/cells"), " already exists in the ", DataFile)
  } else {
    h5write(obj = colnames(srt), file = DataFile, name = paste0(name, "/cells"), level = compression_level)
  }
  if (isTRUE(overwrite)) {
    try(h5delete(file = DataFile, name = paste0(name, "/features")), silent = TRUE)
  }
  if (paste0(name, "/features") %in% h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/features"), " already exists in the ", DataFile)
  } else {
    h5write(obj = unique(unlist(lapply(srt@assays, rownames))), file = DataFile, name = paste0(name, "/features"), level = compression_level)
  }
  return(invisible(NULL))
}

#' CreateMetaFile
#'
#' Creates a meta file in HDF5 format from a Seurat object.
#'
#' @param srt The Seurat object.
#' @param MetaFile Path to the output meta file. If not provided, the file will be named "Meta.hdf5" in the current directory.
#' @param name Name of the dataset. If not provided, the name will default to the Seurat object's project name.
#' @param write_tools A logical value indicating whether to write the tools information to the meta file. Default is FALSE.
#' @param write_misc A logical value indicating whether to write the miscellaneous information to the meta file. Default is FALSE.
#' @param ignore_nlevel The number of levels above which a metadata field will be ignored. Default is 100.
#' @param compression_level The level of compression for the meta file. Default is 6.
#' @param overwrite A logical value indicating whether to overwrite existing metadata and reductions in the meta file. Default is FALSE.
#'
#'  @seealso \code{\link{CreateDataFile}} \code{\link{PrepareSCExplorer}} \code{\link{FetchH5}} \code{\link{RunSCExplorer}}
#'
#' @importFrom Seurat Key Reductions Embeddings
#' @importFrom Matrix sparseMatrix
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset h5delete h5ls h5write
#' @export
CreateMetaFile <- function(srt, MetaFile, name = NULL, write_tools = FALSE, write_misc = FALSE, ignore_nlevel = 100, compression_level = 6, overwrite = FALSE) {
  if (missing(MetaFile) || is.null(MetaFile)) {
    MetaFile <- "Meta.hdf5"
  }
  if (!file.exists(MetaFile)) {
    h5createFile(MetaFile)
  }
  if (is.null(name)) {
    name <- srt@project.name
    message("Set the dataset name to ", name)
  }
  if (substr(name, 1, 1) != "/") {
    name <- paste0("/", name)
  }
  if (!name %in% h5ls(MetaFile)$group) {
    h5createGroup(file = MetaFile, group = name)
  }

  message("Write the meta information to hdf5 file: ", MetaFile)
  if (!paste0(name, "/metadata") %in% h5ls(MetaFile)$group) {
    h5createGroup(file = MetaFile, group = paste0(name, "/metadata"))
  }
  meta_asfeatures <- c()
  meta_asgroups <- c()
  for (var in colnames(srt[[]])) {
    meta <- srt[[var, drop = TRUE]][colnames(srt)]
    if (inherits(meta, "factor")) {
      levels <- levels(meta)
      meta <- as.character(meta)
      attr(meta, "levels") <- levels
      write.attributes <- TRUE
    } else if (inherits(meta, "logical")) {
      levels <- c("TRUE", "FALSE")
      meta <- as.character(meta)
      attr(meta, "levels") <- levels
      write.attributes <- TRUE
    } else {
      write.attributes <- FALSE
    }
    if (is.numeric(meta)) {
      meta <- as.double(meta)
      meta_asfeatures <- c(meta_asfeatures, var)
    } else {
      if (length(unique(meta)) > ignore_nlevel) {
        warning("The number of categories in ", var, " is greater than ", ignore_nlevel, ". It will be ignored.", immediate. = TRUE)
      } else {
        meta_asgroups <- c(meta_asgroups, var)
      }
    }
    if (isTRUE(overwrite)) {
      try(h5delete(file = MetaFile, name = paste0(name, "/metadata/", var)), silent = TRUE)
    }
    if (paste0(name, "/metadata/", var) %in% paste0(h5ls(MetaFile)$group, "/", h5ls(MetaFile)$name)) {
      message("Group ", paste0(name, "/metadata/", var), " already exists in the ", MetaFile)
    } else {
      if (all(is.na(meta))) {
        warning("All of values in ", var, " is NA. It will be ignored.", immediate. = TRUE)
      } else {
        h5write(obj = meta, file = MetaFile, name = paste0(name, "/metadata/", var), write.attributes = write.attributes, level = compression_level)
      }
    }
  }

  if (isTRUE(overwrite)) {
    try(h5delete(file = MetaFile, name = paste0(name, "/metadata.stat")), silent = TRUE)
  }
  if (paste0(name, "/metadata.stat") %in% h5ls(MetaFile)$group) {
    message("Group ", paste0(name, "/metadata.stat"), " already exists in the ", MetaFile)
  } else {
    h5createGroup(file = MetaFile, group = paste0(name, "/metadata.stat"))
    h5write(obj = meta_asfeatures, file = MetaFile, name = paste0(name, "/metadata.stat/asfeatures"), level = compression_level)
    h5write(obj = meta_asgroups, file = MetaFile, name = paste0(name, "/metadata.stat/asgroups"), level = compression_level)
  }

  if (!paste0(name, "/reductions") %in% h5ls(MetaFile)$group) {
    h5createGroup(file = MetaFile, group = paste0(name, "/reductions"))
  }
  for (reduction in Reductions(srt)) {
    emb <- Embeddings(srt, reduction)[colnames(srt), ]
    attr(emb, "key") <- as.character(Key(srt[[reduction]]))
    if (isTRUE(overwrite)) {
      try(h5delete(file = MetaFile, name = paste0(name, "/reductions/", reduction)), silent = TRUE)
    }
    if (paste0(name, "/reductions/", reduction) %in% paste0(h5ls(MetaFile)$group, "/", h5ls(MetaFile)$name)) {
      message("Group ", paste0(name, "/reductions/", reduction), " already exists in the ", MetaFile)
    } else {
      h5createDataset(
        file = MetaFile, dataset = paste0(name, "/reductions/", reduction),
        dims = dim(emb), chunk = c(nrow(emb), min(10, ncol(emb))), level = compression_level
      )
      suppressWarnings(h5write(obj = emb, file = MetaFile, name = paste0(name, "/reductions/", reduction), write.attributes = TRUE, level = compression_level))
    }
  }
  if (isTRUE(overwrite)) {
    try(h5delete(file = MetaFile, name = paste0(name, "/reductions.stat")), silent = TRUE)
  }
  if (paste0(name, "/reductions.stat") %in% h5ls(MetaFile)$group) {
    message("Group ", paste0(name, "/reductions.stat"), " already exists in the ", MetaFile)
  } else {
    h5createGroup(file = MetaFile, group = paste0(name, "/reductions.stat"))
    h5write(obj = DefaultReduction(srt), file = MetaFile, name = paste0(name, "/reductions.stat/Default_reduction"), level = compression_level)
  }

  if (isTRUE(write_misc)) {
    if (isTRUE(overwrite)) {
      try(h5delete(file = MetaFile, name = paste0(name, "/misc")), silent = TRUE)
    }
    if (paste0(name, "/misc") %in% h5ls(MetaFile)$group) {
      message("Group ", paste0(name, "/misc"), " already exists in the ", MetaFile)
    } else {
      h5write(obj = srt@misc, file = MetaFile, name = paste0(name, "/misc"), level = compression_level)
    }
  }
  if (isTRUE(write_tools)) {
    if (isTRUE(overwrite)) {
      try(h5delete(file = MetaFile, name = paste0(name, "/tools")), silent = TRUE)
    }
    if (paste0(name, "/tools") %in% h5ls(MetaFile)$group) {
      message("Group ", paste0(name, "/tools"), " already exists in the ", MetaFile)
    } else {
      h5write(obj = srt@tools, file = MetaFile, name = paste0(name, "/tools"), level = compression_level)
    }
  }

  return(invisible(NULL))
}

#' Prepare Seurat objects for the SCExplorer
#'
#' This function prepares one or multiple Seurat objects for the SCExplorer app. It takes a Seurat object or a list of Seurat objects as input and outputs two hdf5 files: one for the data and one for the metadata.
#'
#' @inheritParams CreateDataFile
#' @inheritParams CreateMetaFile
#' @param object A Seurat object or a list of Seurat objects.
#' @param base_dir The base directory where the SCExplorer hdf5 files will be written. Default is "SCExplorer".
#'
#' @seealso \code{\link{CreateDataFile}} \code{\link{CreateMetaFile}} \code{\link{FetchH5}} \code{\link{RunSCExplorer}}
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
#' }
#' @importFrom Seurat Reductions Assays DefaultAssay
#' @export
PrepareSCExplorer <- function(object,
                              base_dir = "SCExplorer", DataFile = "Data.hdf5", MetaFile = "Meta.hdf5",
                              assays = "RNA", slots = c("counts", "data"),
                              ignore_nlevel = 100, write_tools = FALSE, write_misc = FALSE,
                              compression_level = 6, overwrite = FALSE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  if (!dir.exists(base_dir)) {
    message("Create SCExplorer base directory: ", base_dir)
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  DataFile_full <- paste0(base_dir, "/", DataFile)
  MetaFile_full <- paste0(base_dir, "/", MetaFile)

  if (!is.list(object)) {
    object <- list(object)
  }
  if (any(sapply(object, function(x) !inherits(x, "Seurat")))) {
    stop("'object' must be one Seurat object or a list of Seurat object.")
  }
  if (length(names(object)) > 0 && length(names(object)) != length(object)) {
    stop("The object is named, but the name length is not equal to the number of elements.")
  }
  if (length(names(object)) == 0) {
    names(object) <- make.names(sapply(object, function(x) x@project.name), unique = TRUE)
    message("Set the project name of each seurat object to their dataset name")
  }

  for (i in seq_along(object)) {
    nm <- names(object)[i]
    srt <- object[[nm]]
    message("Prepare data for object: ", nm)
    if (length(Reductions(srt)) == 0) {
      stop("No reduction found in the Seurat object ", i)
    }
    if (!any(assays %in% Assays(srt))) {
      warning("Assay:", assays[!assays %in% Assays(srt)], " is not in the Seurat object ", i, immediate. = TRUE)
      assays <- assays[assays %in% Assays(srt)]
      if (length(assays) == 0) {
        warning("No assays found in the Seurat object ", i, ". Use the default assay to create data file.")
        assays <- DefaultAssay(srt)
      }
    }
    if (!"orig.ident" %in% colnames(srt[[]])) {
      srt[["orig.ident"]] <- "SeuratObject"
    }
    CreateDataFile(
      srt = srt, DataFile = DataFile_full, name = nm, assays = assays, slots = slots,
      compression_level = compression_level, overwrite = overwrite
    )
    CreateMetaFile(
      srt = srt, MetaFile = MetaFile_full, name = nm, ignore_nlevel = ignore_nlevel,
      write_tools = write_tools, write_misc = write_misc,
      compression_level = compression_level, overwrite = overwrite
    )
  }
  return(invisible(NULL))
}

#' Fetch data from the hdf5 file
#'
#' This function fetches data from an hdf5 file. It can fetch gene expression data, metadata, and reduction data from the specified file and returns a Seurat object.
#'
#' @param DataFile The path to the hdf5 file containing the data.
#' @param MetaFile The path to the hdf5 file containing the metadata.
#' @param name The name of the dataset in the hdf5 file. If not specified, the function will attempt to find the shared group name in both files.
#' @param features The names of the genes or features to fetch. If specified, only these features will be fetched.
#' @param slot The slot for the counts in the hdf5 file. If not specified, the first slot will be used.
#' @param assay The name of the assay to use. If not specified, the default assay in the hdf5 file will be used.
#' @param metanames The names of the metadata columns to fetch.
#' @param reduction The name of the reduction to fetch.
#'
#' @return A Seurat object with the fetched data.
#'
#' @seealso \code{\link{CreateDataFile}} \code{\link{CreateMetaFile}} \code{\link{PrepareSCExplorer}} \code{\link{RunSCExplorer}}
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
#' srt <- FetchH5(
#'   DataFile = "./SCExplorer/Data.hdf5",
#'   MetaFile = "./SCExplorer/Meta.hdf5",
#'   features = c("Ins1", "Ghrl"),
#'   metanames = c("SubCellType", "Phase"),
#'   reduction = "UMAP"
#' )
#' CellDimPlot(srt, group.by = c("SubCellType", "Phase"), reduction = "UMAP")
#' FeatureDimPlot(srt, features = c("Ins1", "Ghrl"), reduction = "UMAP")
#' }
#' @importFrom Seurat CreateAssayObject CreateDimReducObject
#' @importFrom Matrix sparseMatrix
#' @importFrom rhdf5 h5ls h5read h5readAttributes
#' @importFrom HDF5Array TENxMatrix
#' @export
FetchH5 <- function(DataFile, MetaFile, name = NULL,
                    features = NULL, slot = NULL, assay = NULL,
                    metanames = NULL,
                    reduction = NULL) {
  if (missing(DataFile) || missing(MetaFile)) {
    stop("'DataFile', 'MetaFile' must be provided.")
  }
  data_group <- h5ls(DataFile)$group
  meta_group <- h5ls(MetaFile)$group

  if (is.null(name)) {
    group <- intersect(data_group, meta_group)
    group <- group[group != "/"]
    group <- sapply(group, function(x) substr(x, 2, nchar(x)))
    if (length(group) == 0) {
      stop("Can not find the shared group name in the DataFile and the MetaFile. They may not correspond to the same project.")
    } else if (length(group) > 1) {
      message(length(group), " possible dataset names were found: ", paste(group, collapse = ", "), "\nUse the first one.")
    }
    name <- group[1]
  }
  if (substr(name, 1, 1) != "/") {
    name <- paste0("/", name)
  }

  all_features <- h5read(DataFile, name = paste0(name, "/features"))
  all_cells <- h5read(DataFile, name = paste0(name, "/cells"))
  meta_struc <- h5ls(MetaFile)
  meta_features_name <- h5read(MetaFile, name = paste0(name, "/metadata.stat/asfeatures"))
  meta_groups_name <- h5read(MetaFile, name = paste0(name, "/metadata.stat/asgroups"))
  reduction_name <- meta_struc[meta_struc$group == paste0(name, "/reductions"), "name"]

  if (!is.null(features) && any(!features %in% c(all_features, meta_features_name))) {
    warning("Can not find the features: ", paste0(features[!features %in% c(all_features, meta_features_name)], collapse = ","), immediate. = TRUE)
  }
  gene_features <- features[features %in% c(all_features)]
  meta_features <- features[features %in% c(meta_features_name)]

  if (!is.null(metanames) && any(!metanames %in% c(meta_features_name, meta_groups_name))) {
    warning("Can not find the meta information: ", paste0(metanames[!metanames %in% c(meta_features_name, meta_groups_name)], collapse = ","), immediate. = TRUE)
  }
  metanames <- metanames[metanames %in% c(meta_features_name, meta_groups_name)]

  if (length(gene_features) == 0 && length(meta_features) == 0 && length(metanames) == 0) {
    stop("No features or meta information found.")
  }

  if (length(gene_features) > 0) {
    if (is.null(assay)) {
      assay <- as.character(h5read(DataFile, name = paste0(name, "/Default_assay")))
    }
    if (is.null(slot)) {
      slots <- as.character(unique(na.omit(sapply(strsplit(data_group[grep(name, data_group)], "/"), function(x) x[4]))))
      slot <- ifelse("counts" %in% slots, "counts", slots[1])
    }
    if (!paste0(name, "/", assay, "/", slot) %in% h5ls(DataFile)[["group"]]) {
      stop("There is no ", paste0(name, "/", assay, "/", slot), " in DataFile, please write it first using the PrepareSCExplorer function")
    }
    data <- TENxMatrix(filepath = DataFile, group = paste0(name, "/", assay, "/", slot))
    gene_features <- gene_features[gene_features %in% colnames(data)]
  }

  if (length(gene_features) > 0) {
    counts <- t(as(data[, gene_features, drop = FALSE], "dgCMatrix")) # matrix,sparseMatrix,dgCMatrix,dgRMatrix
    AssayObject <- CreateAssayObject(counts = counts)
    srt_tmp <- CreateSeuratObject2(assay = assay %||% "RNA", counts = AssayObject)
  } else {
    counts <- matrix(data = 0, ncol = length(all_cells), dimnames = list("empty", all_cells))
    AssayObject <- CreateAssayObject(counts = counts)
    srt_tmp <- CreateSeuratObject2(assay = assay %||% "RNA", counts = AssayObject)
  }

  if (length(c(metanames, meta_features)) > 0) {
    for (i in unique(c(metanames, meta_features))) {
      meta <- h5read(MetaFile, name = paste0(name, "/metadata/", i))
      if (is.array(meta)) {
        if (is.integer(meta)) {
          meta <- as.integer(meta)
        } else if (is.numeric(meta)) {
          meta <- as.numeric(meta)
        } else if (is.character(meta)) {
          meta <- as.character(meta)
        }
      }
      meta_attr <- h5readAttributes(MetaFile, name = paste0(name, "/metadata/", i))
      if ("levels" %in% names(meta_attr)) {
        meta <- factor(meta, levels = meta_attr$levels)
      }
      srt_tmp@meta.data[, i] <- meta
    }
  }

  if (!is.null(reduction)) {
    reduction <- reduction_name[agrep(reduction, reduction_name)]
    if (length(reduction) == 0) {
      stop("Can not find the reduction: ", as.character(reduction))
    }
    for (i in reduction) {
      reduction <- as.matrix(h5read(MetaFile, name = paste0(name, "/reductions/", i)))
      reduction_attr <- h5readAttributes(MetaFile, name = paste0(name, "/reductions/", i))
      rownames(reduction) <- all_cells
      colnames(reduction) <- paste0(as.character(reduction_attr$key), seq_len(ncol(reduction)))
      srt_tmp@reductions[[i]] <- CreateDimReducObject(
        embeddings = reduction,
        key = as.character(reduction_attr$key),
        assay = assay %||% DefaultAssay(srt_tmp)
      )
    }
  }
  return(srt_tmp)
}

#' @importFrom Seurat Key AddMetaData
#' @importFrom utils packageVersion
#' @importFrom methods new
CreateSeuratObject2 <- function(counts, project = "SeuratProject", assay = "RNA", names.field = 1, names.delim = "_", meta.data = NULL, idents = "SeuratObject", ...) {
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(
      x = rownames(x = meta.data),
      y = colnames(x = counts)
    ))) {
      warning("Some cells in meta.data not present in provided counts matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), y = colnames(x = counts)), , drop = FALSE]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = counts))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  if (!length(x = Key(object = counts)) || !nchar(x = Key(object = counts))) {
    Key(object = counts) <- tolower(x = assay)
  }
  assay.list <- list(counts)
  names(x = assay.list) <- assay
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name",
      call. = FALSE, immediate. = TRUE
    )
  }
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = counts))
  }
  names(x = idents) <- colnames(x = counts)
  object <- new(
    Class = "Seurat", assays = assay.list, meta.data = data.frame(row.names = colnames(x = counts)),
    active.assay = assay, active.ident = idents, project.name = project,
    version = packageVersion(pkg = "SeuratObject")
  )
  object@meta.data[["orig.ident"]] <- idents
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

#' RunSCExplorer
#'
#' @param base_dir A string. The base directory of the SCExplorer app. Default is "SCExplorer".
#' @param DataFile A string. The name of the HDF5 file that stores data matrices for each dataset. Default is "Data.hdf5".
#' @param MetaFile A string. The name of the HDF5 file that stores metadata for each dataset. Default is "Meta.hdf5".
#' @param title A string. The title of the SCExplorer app. Default is "SCExplorer".
#' @param initial_dataset A string. The initial dataset to be loaded into the app. Default is NULL.
#' @param initial_reduction A string. The initial dimensional reduction method to be loaded into the app. Default is NULL.
#' @param initial_group A string. The initial variable to group cells in the app. Default is NULL.
#' @param initial_feature A string. The initial feature to be loaded into the app. Default is NULL.
#' @param initial_assay A string. The initial assay to be loaded into the app. Default is NULL.
#' @param initial_slot A string. The initial slot to be loaded into the app. Default is NULL.
#' @param initial_label A string. Whether to add labels in the initial plot. Default is FALSE.
#' @param initial_cell_palette A string. The initial color palette for cells. Default is "Paired".
#' @param initial_feature_palette A string. The initial color palette for features. Default is "Spectral".
#' @param initial_theme A string. The initial theme for plots. Default is "theme_scp".
#' @param initial_size A numeric. The initial size of plots. Default is 4.
#' @param initial_ncol A numeric. The initial number of columns for arranging plots. Default is 3.
#' @param initial_arrange A logical. Whether to use "Row" as the initial arrangement. Default is TRUE.
#' @param initial_raster A logical. Whether to perform rasterization in the initial plot. By default, it is set to automatic, meaning it will be TRUE if the number of cells in the initial datasets exceeds 100,000.
#' @param session_workers A numeric. The number of workers for concurrent execution in an asynchronous programming session. Default is 2.
#' @param plotting_workers A numeric. The number of threads per worker for parallel plotting. Default is 8.
#' @param create_script A logical. Whether to create the SCExplorer app script. Default is TRUE.
#' @param style_script A logical. Whether to style the SCExplorer app script. Default is TRUE.
#' @param overwrite A logical. Whether to overwrite existing files. Default is FALSE.
#' @param return_app A logical. Whether to return the SCExplorer app. Default is TRUE.
#'
#' @seealso \code{\link{CreateDataFile}} \code{\link{CreateMetaFile}} \code{\link{PrepareSCExplorer}} \code{\link{FetchH5}}
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(srtMerge = panc8_sub, batch = "tech", integration_method = "Seurat")
#'
#' PrepareSCExplorer(list(mouse_pancreas = pancreas_sub, human_pancreas = panc8_sub), base_dir = "./SCExplorer", overwrite = TRUE)
#'
#' # Create the app.R script
#' app <- RunSCExplorer(
#'   base_dir = "./SCExplorer",
#'   initial_dataset = "mouse_pancreas",
#'   initial_group = "CellType",
#'   initial_feature = "Neurog3",
#'   session_workers = 2,
#'   overwrite = TRUE
#' )
#' list.files("./SCExplorer") # This directory can be used as site directory for Shiny Server.
#'
#' # Run shiny app
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#' # Note: If SCP installed in the isolated environment using renv, you need to add `renv::activate(project = "path/to/SCP_env")` to the app.R script.
#'
#' ####################################################################################################################
#' # You can deploy the app on the self-hosted shiny server(https://www.rstudio.com/products/shiny/shiny-server/).
#' # Or deploy the app on the website(https://www.shinyapps.io) for free:
#'
#' ### step1: set the repository URL for Bioconductor packages and update them to the latest version
#' # options(repos = BiocManager::repositories())
#' # BiocManager::install(ask = FALSE)
#'
#' ### step2: install "rsconnect" package and authorize your account
#' # install.packages("rsconnect")
#' # library(rsconnect)
#' # setAccountInfo(name = "<NAME>", token = "<TOKEN>", secret = "<SECRET>")
#'
#' ### step3: deploy the app
#' # deployApp("./SCExplorer")
#' }
#'
#' @importFrom utils packageVersion
#' @export
RunSCExplorer <- function(base_dir = "SCExplorer",
                          DataFile = "Data.hdf5",
                          MetaFile = "Meta.hdf5",
                          title = "SCExplorer",
                          initial_dataset = NULL,
                          initial_reduction = NULL,
                          initial_group = NULL,
                          initial_feature = NULL,
                          initial_assay = NULL,
                          initial_slot = NULL,
                          initial_label = FALSE,
                          initial_cell_palette = "Paired",
                          initial_feature_palette = "Spectral",
                          initial_theme = "theme_scp",
                          initial_size = 4,
                          initial_ncol = 3,
                          initial_arrange = NULL,
                          initial_raster = NULL,
                          session_workers = 2,
                          plotting_workers = 8,
                          create_script = TRUE,
                          style_script = require("styler", quietly = TRUE),
                          overwrite = FALSE,
                          return_app = TRUE) {
  check_R(c("rhdf5", "HDF5Array", "shiny@1.6.0", "ggplot2", "ragg", "htmlwidgets", "plotly", "bslib", "future", "promises", "BiocParallel"))
  DataFile_full <- paste0(base_dir, "/", DataFile)
  MetaFile_full <- paste0(base_dir, "/", MetaFile)
  if (!file.exists(DataFile_full) || !file.exists(MetaFile_full)) {
    stop("Please create the DataFile and MetaFile using PrepareSCExplorer function first!")
  }

  main_code <- '
if (!file.exists("Rplots.pdf")) {
  file.create("Rplots.pdf")
}
data_group <- rhdf5::h5ls(DataFile)$group
meta_group <- rhdf5::h5ls(MetaFile)$group
group <- intersect(data_group, meta_group)
group <- group[group != "/"]
group <- as.character(sapply(group, function(x) substr(x, 2, nchar(x))))
if (length(group) == 0) {
  stop("Can not find the shared group names in the DataFile and the MetaFile. They may not correspond to the same project.")
}
if (is.null(initial_dataset)) {
  initial_dataset <- group[1]
}
if (substr(initial_dataset, 1, 1) == "/") {
  initial_dataset <- substr(initial_dataset, 2, nchar(initial_dataset))
}
if (!initial_dataset %in% group) {
  stop("Dataset ", group, " is not in the DataFile and the MetaFile")
}

assays <- unique(na.omit(sapply(strsplit(data_group[grep(initial_dataset, data_group)], "/"), function(x) x[3])))
slots <- unique(na.omit(sapply(strsplit(data_group[grep(initial_dataset, data_group)], "/"), function(x) x[4])))
if (is.null(initial_assay)) {
  initial_assay <- as.character(rhdf5::h5read(DataFile, name = paste0("/", initial_dataset, "/Default_assay")))
}
if (is.null(initial_slot)) {
  initial_slot <- ifelse("data" %in% slots, "data", slots[1])
}
if (!initial_assay %in% assays) {
  stop("initial_assay is not in the dataset ", initial_dataset, " in the DataFile")
}
if (!initial_slot %in% slots) {
  stop("initial_slot is not in the dataset ", initial_slot, " in the DataFile")
}
all_cells <- rhdf5::h5read(DataFile, name = paste0("/", initial_dataset, "/cells"))

data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", initial_dataset, "/", initial_assay, "/", initial_slot))
all_features <- colnames(data)

meta_struc <- rhdf5::h5ls(MetaFile)
meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", initial_dataset, "/metadata.stat/asfeatures"))
meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", initial_dataset, "/metadata.stat/asgroups"))
reduction_name <- meta_struc[meta_struc$group == paste0("/", initial_dataset, "/reductions"), "name"]
default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", initial_dataset, "/reductions.stat/Default_reduction")))

if (is.null(initial_reduction)) {
  initial_reduction <- default_reduction
}
if (is.null(initial_group)) {
  initial_group <- "orig.ident"
}
if (is.null(initial_feature)) {
  initial_feature <- meta_features_name[1]
}
if (is.null(initial_arrange)) {
  initial_arrange <- TRUE
}
if (is.null(initial_raster)) {
  initial_raster <- length(all_cells) > 1e5
}

palette_list <- SCP::palette_list
theme_list <- list(
  SCP = c("theme_scp", "theme_blank"),
  ggplot2 = c("theme_classic", "theme_linedraw", "theme_minimal", "theme_void", "theme_grey", "theme_dark", "theme_light")
)
themes <- setNames(rep(names(theme_list), sapply(theme_list, length)), unlist(theme_list))
panel_raster <- FALSE

ui <- fluidPage(
  theme = page_theme,
  navbarPage(
    title = title,
    # 1. Cell dimensional reduction plot ----------------------------------------------------------------------
    tabPanel(
      title = "Cell dimensional reduction plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset1",
            label = "Select a dataset",
            choices = group,
            selected = initial_dataset
          ),
          selectInput(
            inputId = "reduction1",
            label = "Select a reduction",
            choices = reduction_name,
            selected = initial_reduction
          ),
          selectInput(
            inputId = "group1",
            label = "Select a variable to group cells",
            choices = meta_groups_name,
            selected = initial_group
          ),
          selectInput(
            inputId = "split1",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          selectInput(
            inputId = "palette1",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_cell_palette
          ),
          selectInput(
            inputId = "theme1",
            label = "Select a theme",
            choices = names(themes),
            selected = initial_theme
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "label1",
                label = "Label",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = initial_label,
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "raster1",
                label = "Raster",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = initial_raster,
                inline = TRUE
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "pt_size1",
                label = "Point size",
                value = 1,
                min = 0.1,
                max = 10,
                step = 0.5,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "size1",
                label = "Panel size",
                value = initial_size,
                min = 1,
                max = 10,
                step = 0.1,
                width = "150px"
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "ncol1",
                label = "Number of columns",
                value = initial_ncol,
                min = 1,
                max = 100,
                step = 1,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "arrange1",
                label = "Arrange by",
                choices = c("Row" = TRUE, "Column" = FALSE),
                selected = initial_arrange
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              actionButton(inputId = "submit1", label = "Submit", icon("play"), class = "btn-info", width = "150px")
            ),
            column(
              width = 6, align = "center",
              downloadButton(outputId = "download1", label = "Download", class = "btn-warning", width = "150px")
            )
          )
        ),
        mainPanel(
          width = 9,
          fluidPage(
            tabsetPanel(
              tabPanel(
                title = "2D plot",
                column(
                  width = 12, offset = 0, style = "padding:0px;margin:0%",
                  div(
                    style = "overflow-x: auto;",
                    uiOutput("plot1")
                  )
                )
              ),
              tabPanel(
                title = "3D plot",
                column(
                  width = 12, offset = 0, style = "padding:0px;margin:0%",
                  div(
                    style = "overflow-x: auto;",
                    plotly::plotlyOutput("plot1_3d", height = "100%", width = "100%")
                  )
                )
              )
            )
          )
        )
      )
    ),
    # 2. Feature dimensional reduction plot ----------------------------------------------------------------------
    tabPanel(
      title = "Feature dimensional reduction plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset2",
            label = "Select a dataset",
            choices = group,
            selected = initial_dataset
          ),
          selectInput(
            inputId = "reduction2",
            label = "Select a reduction",
            choices = reduction_name,
            selected = initial_reduction
          ),
          selectInput(
            inputId = "split2",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = "assays2",
                label = "Select a assay",
                choices = assays,
                selected = initial_assay
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = "slots2",
                label = "Select a slot",
                choices = slots,
                selected = initial_slot
              )
            ),
          ),
          selectizeInput(
            inputId = "features2",
            label = "Select features",
            choices = NULL,
            selected = initial_feature,
            multiple = TRUE,
            options = list(maxOptions = 20, maxItems = 20)
          ),
          textAreaInput(
            inputId = "feature_area2",
            label = "Input features",
            height = "100px",
            placeholder = paste(sample(all_features, 4), collapse = "\\n")
          ),
          selectInput(
            inputId = "palette2",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_feature_palette
          ),
          selectInput(
            inputId = "theme2",
            label = "Select a theme",
            choices = names(themes),
            selected = initial_theme
          ),
          fluidRow(
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "coExp2",
                label = "Co-expression",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              ),
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "scale2",
                label = "Color scale",
                choices = list("feature" = "feature", "all" = "all"),
                selected = "feature",
                inline = TRUE
              ),
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "raster2",
                label = "Raster",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = initial_raster,
                inline = TRUE
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "pt_size2",
                label = "Point size",
                value = 1,
                min = 0.1,
                max = 10,
                step = 0.5,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "size2",
                label = "Panel size",
                value = initial_size,
                min = 1,
                max = 10,
                step = 0.1,
                width = "150px"
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "ncol2",
                label = "Number of columns",
                value = initial_ncol,
                min = 1,
                max = 100,
                step = 1,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "arrange2",
                label = "Arrange by",
                choices = c("Row" = TRUE, "Column" = FALSE),
                selected = initial_arrange
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              actionButton(inputId = "submit2", label = "Submit", icon("play"), class = "btn-info", width = "150px")
            ),
            column(
              width = 6, align = "center",
              downloadButton(outputId = "download2", label = "Download", class = "btn-warning", width = "150px")
            )
          )
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel(
              title = "2D plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  uiOutput("plot2")
                )
              )
            ),
            tabPanel(
              title = "3D plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  plotly::plotlyOutput("plot2_3d", height = "100%", width = "100%")
                )
              )
            )
          )
        )
      )
    ),
    # 3. Cell statistical plot ----------------------------------------------------------------------
    tabPanel(
      title = "Cell statistical plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset3",
            label = "Select a dataset",
            choices = group,
            selected = initial_dataset
          ),
          selectInput(
            inputId = "group3",
            label = "Select a variable to group cells",
            choices = c("None", meta_groups_name),
            selected = initial_group
          ),
          selectInput(
            inputId = "groupuse3",
            label = "Group use",
            choices = "All",
            selected = "All",
            multiple = TRUE
          ),
          selectInput(
            inputId = "split3",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          selectInput(
            inputId = "stat3",
            label = "Select a variable to be counted",
            choices = meta_groups_name,
            selected = initial_group
          ),
          fluidRow(
            column(
              width = 4,
              selectInput(
                inputId = "plottype3",
                label = "Plot type",
                choices = c("bar", "rose", "ring", "pie", "trend", "area", "dot"),
                selected = "bar"
              )
            ),
            column(
              width = 4,
              selectInput(
                inputId = "stattype3",
                label = "Measurement",
                choices = c("percent", "count"),
                selected = "stack"
              )
            ),
            column(
              width = 4,
              selectInput(
                inputId = "position3",
                label = "Position",
                choices = c("stack", "dodge"),
                selected = "stack"
              )
            ),
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "label3",
                label = "Label",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = initial_label,
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "flip3",
                label = "Flip",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
          ),
          selectInput(
            inputId = "palette3",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_cell_palette
          ),
          selectInput(
            inputId = "theme3",
            label = "Select a theme",
            choices = names(themes),
            selected = initial_theme
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "aspectratio3",
                label = "Aspect ratio",
                choices = c("auto", "custom"),
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              conditionalPanel(
                condition = "input.aspectratio3 == \'custom\'",
                numericInput(
                  inputId = "aspectratio_value3",
                  label = NULL,
                  value = 1,
                  min = 0,
                  max = 100,
                  step = 0.1,
                  width = "150px"
                )
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "labelsize3",
                label = "Label size",
                value = 3.5,
                min = 1,
                max = 10,
                step = 0.1,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "size3",
                label = "Panel size",
                value = initial_size,
                min = 1,
                max = 10,
                step = 0.1,
                width = "150px"
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "ncol3",
                label = "Number of columns",
                value = initial_ncol,
                min = 1,
                max = 100,
                step = 1,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "arrange3",
                label = "Arrange by",
                choices = c("Row" = TRUE, "Column" = FALSE),
                selected = initial_arrange
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              actionButton(inputId = "submit3", label = "Submit", icon("play"), class = "btn-info", width = "150px")
            ),
            column(
              width = 6, align = "center",
              downloadButton(outputId = "download3", label = "Download", class = "btn-warning", width = "150px")
            )
          )
        ),
        mainPanel(
          width = 9,
          fluidPage(
            tabsetPanel(
              tabPanel(
                title = "Statistical plot",
                column(
                  width = 12, offset = 0, style = "padding:0px;margin:0%",
                  div(
                    style = "overflow-x: auto;",
                    uiOutput("plot3")
                  )
                )
              )
            )
          )
        )
      )
    ),
    # 4. Feature statistical plot ----------------------------------------------------------------------
    tabPanel(
      title = "Feature statistical plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset4",
            label = "Select a dataset",
            choices = group,
            selected = initial_dataset
          ),
          selectInput(
            inputId = "group4",
            label = "Select a variable to group cells",
            choices = c("None", meta_groups_name),
            selected = initial_group
          ),
          selectInput(
            inputId = "groupuse4",
            label = "Group use",
            choices = "All",
            selected = "All",
            multiple = TRUE
          ),
          selectInput(
            inputId = "split4",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = "assays4",
                label = "Select a assay",
                choices = assays,
                selected = initial_assay
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = "slots4",
                label = "Select a slot",
                choices = slots,
                selected = initial_slot
              )
            ),
          ),
          selectInput(
            inputId = "plottype4",
            label = "Select a plot type",
            choices = c("violin", "box", "bar", "dot", "col"),
            selected = "violin"
          ),
          selectizeInput(
            inputId = "features4",
            label = "Select features",
            choices = NULL,
            selected = initial_feature,
            multiple = TRUE,
            options = list(maxOptions = 20, maxItems = 20)
          ),
          textAreaInput(
            inputId = "feature_area4",
            label = "Input features",
            height = "100px",
            placeholder = paste(sample(all_features, 4), collapse = "\\n")
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "plotby4",
                label = "Plot by",
                choices = c("group", "feature"),
                selected = "group",
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "fillby4",
                label = "Fill by",
                choices = c("group", "feature", "expression"),
                selected = "group",
                inline = TRUE
              )
            )
          ),
          fluidRow(
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "coExp4",
                label = "Co-expression",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "stack4",
                label = "Stack",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "flip4",
                label = "Flip",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            )
          ),
          fluidRow(
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "addbox4",
                label = "Add box",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "addpoint4",
                label = "Add point",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 4, align = "center",
              radioButtons(
                inputId = "addtrend4",
                label = "Add trend",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            )
          ),
          selectInput(
            inputId = "palette4",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_cell_palette
          ),
          selectInput(
            inputId = "theme4",
            label = "Select a theme",
            choices = names(themes),
            selected = initial_theme
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "aspectratio4",
                label = "Aspect ratio",
                choices = c("auto", "custom"),
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              conditionalPanel(
                condition = "input.aspectratio4 == \'custom\'",
                numericInput(
                  inputId = "aspectratio_value4",
                  label = NULL,
                  value = 1,
                  min = 0,
                  max = 100,
                  step = 0.1,
                  width = "150px"
                )
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "sameylims4",
                label = "Same y-axis",
                choices = c("Yes" = TRUE, "No" = FALSE),
                selected = FALSE,
                inline = TRUE
              )
            ),
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "size4",
                label = "Panel size",
                value = initial_size,
                min = 1,
                max = 10,
                step = 0.1,
                width = "150px"
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              numericInput(
                inputId = "ncol4",
                label = "Number of columns",
                value = initial_ncol,
                min = 1,
                max = 100,
                step = 1,
                width = "150px"
              )
            ),
            column(
              width = 6, align = "center",
              radioButtons(
                inputId = "arrange4",
                label = "Arrange by",
                choices = c("Row" = TRUE, "Column" = FALSE),
                selected = initial_arrange
              )
            )
          ),
          fluidRow(
            column(
              width = 6, align = "center",
              actionButton(inputId = "submit4", label = "Submit", icon("play"), class = "btn-info", width = "150px")
            ),
            column(
              width = 6, align = "center",
              downloadButton(outputId = "download4", label = "Download", class = "btn-warning", width = "150px")
            )
          )
        ),
        mainPanel(
          width = 9,
          fluidPage(
            tabsetPanel(
              tabPanel(
                title = "Statistical plot",
                column(
                  width = 12, offset = 0, style = "padding:0px;margin:0%",
                  div(
                    style = "overflow-x: auto;",
                    uiOutput("plot4")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Initial  ----------------------------------------------------------------
  promisedData <- reactiveValues()

  get_attr <- function(x, attr, verbose = FALSE) {
    width <- attr(x, "size")$width
    height <- attr(x, "size")$height
    units <- attr(x, "size")$units
    dpi <- attr(x, "dpi")
    if (verbose) {
      message(paste("width:", width, "height:", height, "units:", units, "dpi:", dpi))
    }
    if (attr == "width") {
      return(width)
    }
    if (attr == "height") {
      return(height)
    }
    if (attr == "units") {
      return(units)
    }
    if (attr == "dpi") {
      return(dpi)
    }
  }

  # change dataset  ----------------------------------------------------------------
  observe({
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset1, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/reductions.stat/Default_reduction")))
    all_cells <- rhdf5::h5read(DataFile, name = paste0("/", input$dataset1, "/cells"))
    updateSelectInput(session, "reduction1", choices = reduction_name, selected = intersect(c(initial_reduction, default_reduction), reduction_name)[1])
    updateSelectInput(session, "group1", choices = meta_groups_name, selected = intersect(c(initial_group, "orig.ident"), meta_groups_name)[1])
    updateSelectInput(session, "split1", choices = c("None", meta_groups_name), selected = "None")
    updateRadioButtons(session, "raster1", choices = c("Yes" = TRUE, "No" = FALSE), selected = initial_raster)
  }) %>% bindEvent(input$dataset1, ignoreNULL = TRUE, ignoreInit = FALSE)

  observe({
    assays <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[3])))
    slots <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[4])))
    default_assay <- as.character(rhdf5::h5read(DataFile, name = paste0("/", input$dataset2, "/Default_assay")))
    default_slot <- ifelse("data" %in% slots, "data", slots[1])
    assay <- intersect(c(initial_assay, default_assay), assays)[1]
    slot <- intersect(c(initial_slot, default_slot), slots)[1]
    updateSelectInput(session, "assays2", choices = assays, selected = assay)
    updateSelectInput(session, "slots2", choices = slots, selected = slot)
    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", input$dataset2, "/", assay, "/", slot))
    all_features <- colnames(data)
    all_cells <- rhdf5::h5read(DataFile, name = paste0("/", input$dataset2, "/cells"))
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/metadata.stat/asfeatures"))
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset2, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/reductions.stat/Default_reduction")))
    updateSelectInput(session, "reduction2", choices = reduction_name, selected = intersect(c(initial_reduction, default_reduction), reduction_name)[1])
    updateSelectizeInput(session, "features2",
      choices = c(meta_features_name, all_features), selected = intersect(c(initial_feature, meta_features_name[1]), c(all_features, meta_features_name))[1],
      options = list(maxOptions = 20, maxItems = 20), server = TRUE
    )
    updateSelectInput(session, "split2", choices = c("None", meta_groups_name), selected = "None")
    updateRadioButtons(session, "raster2", choices = c("Yes" = TRUE, "No" = FALSE), selected = initial_raster)
  }) %>% bindEvent(input$dataset2, ignoreNULL = TRUE, ignoreInit = FALSE)

  observe({
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset3, "/metadata.stat/asgroups"))
    updateSelectInput(session, "stat3", choices = meta_groups_name, selected = "orig.ident")
    updateSelectInput(session, "group3", choices = meta_groups_name, selected = intersect(c(initial_group, "orig.ident"), meta_groups_name)[1])
    updateSelectInput(session, "split3", choices = c("None", meta_groups_name), selected = "None")
  }) %>% bindEvent(input$dataset3, ignoreNULL = TRUE, ignoreInit = FALSE)

  observe({
    assays <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset4, data_group)], "/"), function(x) x[3])))
    slots <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset4, data_group)], "/"), function(x) x[4])))
    default_assay <- as.character(rhdf5::h5read(DataFile, name = paste0("/", input$dataset4, "/Default_assay")))
    default_slot <- ifelse("data" %in% slots, "data", slots[1])
    assay <- intersect(c(initial_assay, default_assay), assays)[1]
    slot <- intersect(c(initial_slot, default_slot), slots)[1]
    updateSelectInput(session, "assays4", choices = assays, selected = assay)
    updateSelectInput(session, "slots4", choices = slots, selected = slot)
    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", input$dataset4, "/", assay, "/", slot))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset4, "/metadata.stat/asfeatures"))
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset4, "/metadata.stat/asgroups"))
    updateSelectizeInput(session, "features4",
      choices = c(meta_features_name, all_features), selected = intersect(c(initial_feature, meta_features_name[1]), c(all_features, meta_features_name))[1],
      options = list(maxOptions = 20, maxItems = 20), server = TRUE
    )
    updateSelectInput(session, "split4", choices = c("None", meta_groups_name), selected = "None")
    updateSelectInput(session, "group4", choices = meta_groups_name, selected = intersect(c(initial_group, "orig.ident"), meta_groups_name)[1])
  }) %>% bindEvent(input$dataset4, ignoreNULL = TRUE, ignoreInit = FALSE)

  observe({
    if (input$group3 != "None") {
      groups <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset3, "/metadata/", input$group3))
      grouplevels <- rhdf5::h5readAttributes(MetaFile, name = paste0("/", input$dataset3, "/metadata/", input$group3))
      if ("levels" %in% names(grouplevels)) {
        groups <- factor(groups, levels = grouplevels$levels)
      }
      updateSelectInput(session, "groupuse3", choices = c("All", if (is.null(levels(groups))) unique(groups) else levels(groups)), selected = "All")
    } else {
      updateSelectInput(session, "groupuse3", choices = "All", selected = "All")
    }
  }) %>% bindEvent(input$group3, ignoreNULL = TRUE, ignoreInit = FALSE)

  observe({
    if (input$group4 != "None") {
      groups <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset4, "/metadata/", input$group4))
      grouplevels <- rhdf5::h5readAttributes(MetaFile, name = paste0("/", input$dataset4, "/metadata/", input$group4))
      if ("levels" %in% names(grouplevels)) {
        groups <- factor(groups, levels = grouplevels$levels)
      }
      updateSelectInput(session, "groupuse4", choices = c("All", if (is.null(levels(groups))) unique(groups) else levels(groups)), selected = "All")
    } else {
      updateSelectInput(session, "groupuse4", choices = "All", selected = "All")
    }
  }) %>% bindEvent(input$group4, ignoreNULL = TRUE, ignoreInit = FALSE)

  # submit1  ----------------------------------------------------------------
  r1 <- reactive({
    dataset1 <- input$dataset1
    reduction1 <- input$reduction1
    group1 <- input$group1
    if (input$split1 == "None") {
      split1 <- NULL
    } else {
      split1 <- input$split1
    }
    label1 <- input$label1
    raster1 <- input$raster1
    palette1 <- input$palette1
    theme1 <- input$theme1
    size1 <- input$size1
    pt_size1 <- input$pt_size1
    ncol1 <- input$ncol1
    byrow1 <- input$arrange1

    # lapply(grep("1$",names(input),value = TRUE), function(x)print(paste0(x,":",input[[x]])))

    promisedData[["p1_dim"]] <- NULL
    promisedData[["p1_3d"]] <- NULL
    promises::future_promise(
      {
        # print("******************************** New task ********************************")
        # print(">>> fetch data:")
        # print(system.time(
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset1,
          metanames = unique(c(group1, split1)), reduction = reduction1
        )
        # ))

        theme1 <- get(theme1, envir = asNamespace(themes[theme1]))

        # print(">>> plot:")
        # print(system.time(
        p1_dim <- SCP::CellDimPlot(srt_tmp,
          group.by = group1, split.by = split1, reduction = reduction1, raster = raster1, pt.size = pt_size1,
          label = label1, palette = palette1, theme_use = theme1,
          ncol = ncol1, byrow = byrow1, force = TRUE
        )
        # ))

        # print(">>> panel_fix:")
        # print(system.time(
        p1_dim <- SCP::panel_fix(SCP::slim_data(p1_dim), height = size1, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        # ))
        attr(p1_dim, "dpi") <- 300
        plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
        if (isTRUE(plot3d)) {
          p1_3d <- SCP::CellDimPlot3D(srt_tmp, group.by = group1, pt.size = pt_size1 * 2, reduction = reduction1, palette = palette1, force = TRUE)
        } else {
          p1_3d <- NULL
        }
        return(list(p1_dim, p1_3d))
      },
      seed = TRUE
    )
  }) %>%
    bindCache(
      input$dataset1, input$reduction1, input$group1, input$split1,
      input$palette1, input$theme1, input$label1, input$raster1,
      input$pt_size1, input$size1, input$ncol1, input$arrange1
    ) %>%
    bindEvent(input$submit1, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    prog <- Progress$new(min = 1, max = 10)
    prog$set(value = 3, message = "Fetch data...", detail = "[Cell dimensional reduction plot]")
    r1()$then(function(x) {
      promisedData[["p1_dim"]] <- x[[1]]
      promisedData[["p1_3d"]] <- x[[2]]
      width <- get_attr(x[[1]], "width")
      height <- get_attr(x[[1]], "height")
      dpi <- get_attr(x[[1]], "dpi")

      prog$set(value = 8, message = "Render plot...", detail = "[Cell dimensional reduction plot]")
      # print("renderPlot:")
      # print(system.time(
      output$plot1 <- renderUI({
        renderPlot(
          {
            x[[1]]
          },
          width = width * 96,
          height = height * 96,
          res = 96
        )
      })
      # ))

      # print("renderPlotly:")
      # print(system.time(
      output$plot1_3d <- plotly::renderPlotly({
        x[[2]]
      })
      # ))
    }) %>%
      finally(~ {
        prog$set(value = 10, message = "Done.", detail = "[Cell dimensional reduction plot]")
        prog$close()
      })
  }) %>% bindEvent(input$submit1, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$download1 <- downloadHandler(
    filename = function() {
      paste0("CellDimPlot-", format(Sys.time(), "%Y%m%d%H%M%S"), ".zip")
    },
    content = function(file) {
      width <- get_attr(promisedData[["p1_dim"]], "width")
      height <- get_attr(promisedData[["p1_dim"]], "height")
      dpi <- get_attr(promisedData[["p1_dim"]], "dpi")

      temp1 <- tempfile(pattern = "CellDimPlot-", fileext = ".png")
      ggplot2::ggsave(filename = temp1, plot = promisedData[["p1_dim"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      temp2 <- tempfile(pattern = "CellDimPlot-", fileext = ".pdf")
      ggplot2::ggsave(filename = temp2, plot = promisedData[["p1_dim"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      if (!is.null(promisedData[["p1_3d"]])) {
        temp3 <- tempfile(pattern = "CellDimPlot3D-", fileext = ".html")
        htmlwidgets::saveWidget(
          widget = plotly::as_widget(promisedData[["p1_3d"]]),
          file = temp3
        )
        unlink(gsub("\\\\.html", "_files", temp3), recursive = TRUE)
      } else {
        temp3 <- NULL
      }

      zip(zipfile = file, flags = "-jq", files = c(temp1, temp2, temp3))
    },
    contentType = "application/zip"
  )

  # submit2  ----------------------------------------------------------------
  r2 <- reactive({
    dataset2 <- input$dataset2
    reduction2 <- input$reduction2
    if (input$split2 == "None") {
      split2 <- NULL
    } else {
      split2 <- input$split2
    }
    assays2 <- input$assays2
    slots2 <- input$slots2
    features2 <- input$features2
    feature_area2 <- input$feature_area2
    coExp2 <- input$coExp2
    scale2 <- input$scale2
    raster2 <- input$raster2
    palette2 <- input$palette2
    theme2 <- input$theme2
    size2 <- input$size2
    pt_size2 <- input$pt_size2
    ncol2 <- input$ncol2
    byrow2 <- input$arrange2

    # lapply(grep("2$",names(input),value = TRUE), function(x)print(paste0(x,":",input[[x]])))

    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", dataset2, "/", assays2, "/", slots2))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", dataset2, "/metadata.stat/asfeatures"))

    feature_area2 <- gsub(x = unlist(strsplit(feature_area2, "(\\r)|(\\n)", perl = TRUE)), pattern = " ", replacement = "")
    features2 <- c(as.character(features2), as.character(feature_area2))
    features2 <- unique(features2[features2 %in% c(all_features, meta_features_name)])
    if (length(features2) == 0) {
      features2 <- intersect(c(initial_feature, meta_features_name[1]), c(all_features, meta_features_name))[1]
    }

    promisedData[["p2_dim"]] <- NULL
    promisedData[["p2_3d"]] <- NULL
    promises::future_promise(
      {
        # print("******************************** New task ********************************")
        # print(">>> fetch data:")
        # print(system.time(
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset2,
          features = features2, slot = slots2, assay = assays2,
          metanames = split2, reduction = reduction2
        )
        # ))

        theme2 <- get(theme2, envir = asNamespace(themes[theme2]))

        # print(">>> plot:")
        # print(system.time(
        p2_dim <- SCP::FeatureDimPlot(
          srt = srt_tmp, features = features2, split.by = split2, reduction = reduction2, slot = "data", raster = raster2, pt.size = pt_size2,
          calculate_coexp = coExp2, keep_scale = scale2, palette = palette2, theme_use = theme2,
          ncol = ncol2, byrow = byrow2, force = TRUE
        )
        # ))

        # print(">>> panel_fix:")
        # print(system.time(
        p2_dim <- SCP::panel_fix(SCP::slim_data(p2_dim), height = size2, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        # ))
        attr(p2_dim, "dpi") <- 300
        plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
        if (isTRUE(plot3d)) {
          p2_3d <- SCP::FeatureDimPlot3D(
            srt = srt_tmp, features = features2, reduction = reduction2, pt.size = pt_size2 * 2,
            calculate_coexp = coExp2, force = TRUE
          )
        } else {
          p2_3d <- NULL
        }
        return(list(p2_dim, p2_3d))
      },
      seed = TRUE
    )
  }) %>%
    bindCache(
      input$dataset2, input$reduction2, input$split2, input$assays2, input$slots2,
      input$features2, input$feature_area2,
      input$palette2, input$theme2, input$coExp2, input$scale2, input$raster2,
      input$pt_size2, input$size2, input$ncol2, input$arrange2
    ) %>%
    bindEvent(input$submit2, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$download2 <- downloadHandler(
    filename = function() {
      paste0("FeatureDimPlot-", format(Sys.time(), "%Y%m%d%H%M%S"), ".zip")
    },
    content = function(file) {
      width <- get_attr(promisedData[["p2_dim"]], "width")
      height <- get_attr(promisedData[["p2_dim"]], "height")
      dpi <- get_attr(promisedData[["p2_dim"]], "dpi")

      temp1 <- tempfile(pattern = "FeatureDimPlot-", fileext = ".png")
      ggplot2::ggsave(filename = temp1, plot = promisedData[["p2_dim"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      temp2 <- tempfile(pattern = "FeatureDimPlot-", fileext = ".pdf")
      ggplot2::ggsave(filename = temp2, plot = promisedData[["p2_dim"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      if (!is.null(promisedData[["p2_3d"]])) {
        temp3 <- tempfile(pattern = "FeatureDimPlot3D-", fileext = ".html")
        htmlwidgets::saveWidget(
          widget = plotly::as_widget(promisedData[["p2_3d"]]),
          file = temp3
        )
        unlink(gsub("\\\\.html", "_files", temp3), recursive = TRUE)
      } else {
        temp3 <- NULL
      }

      zip(zipfile = file, flags = "-jq", files = c(temp1, temp2, temp3))
    },
    contentType = "application/zip"
  )

  observe({
    prog <- Progress$new(min = 1, max = 10)
    prog$set(value = 3, message = "Fetch data...", detail = "[Feature dimensional reduction plot]")
    r2()$then(function(x) {
      promisedData[["p2_dim"]] <- x[[1]]
      promisedData[["p2_3d"]] <- x[[2]]
      width <- get_attr(x[[1]], "width")
      height <- get_attr(x[[1]], "height")
      dpi <- get_attr(x[[1]], "dpi")

      prog$set(value = 8, message = "Render plot...", detail = "[Feature dimensional reduction plot]")
      # print("renderPlot:")
      # print(system.time(
      output$plot2 <- renderUI({
        renderPlot(
          {
            x[[1]]
          },
          width = width * 96,
          height = height * 96,
          res = 96
        )
      })
      # ))

      # print("renderPlotly:")
      # print(system.time(
      output$plot2_3d <- plotly::renderPlotly({
        x[[2]]
      })
      # ))
    }) %>%
      finally(~ {
        prog$set(value = 10, message = "Done.", detail = "[Feature dimensional reduction plot]")
        prog$close()
      })
  }) %>% bindEvent(input$submit2, ignoreNULL = FALSE, ignoreInit = FALSE)

  # submit3  ----------------------------------------------------------------
  r3 <- reactive({
    dataset3 <- input$dataset3
    plottype3 <- input$plottype3
    stattype3 <- input$stattype3
    position3 <- input$position3
    stat3 <- input$stat3
    group3 <- input$group3
    groupuse3 <- input$groupuse3
    if (input$group3 == "None") {
      group3 <- NULL
    } else {
      group3 <- input$group3
    }
    if (input$split3 == "None") {
      split3 <- NULL
    } else {
      split3 <- input$split3
    }
    label3 <- input$label3
    flip3 <- input$flip3
    palette3 <- input$palette3
    theme3 <- input$theme3
    labelsize3 <- input$labelsize3
    size3 <- input$size3
    ncol3 <- input$ncol3
    byrow3 <- input$arrange3
    if (input$aspectratio3 == "auto") {
      aspect.ratio <- NULL
    } else {
      aspect.ratio <- input$aspectratio_value3
    }

    # lapply(grep("3$",names(input),value = TRUE), function(x)print(paste0(x,":",input[[x]])))

    promisedData[["p3"]] <- NULL
    promises::future_promise(
      {
        # print("******************************** New task ********************************")
        # print(">>> fetch data:")
        # print(system.time(
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset3,
          metanames = unique(c(stat3, group3, split3))
        )
        # ))

        theme3 <- get(theme3, envir = asNamespace(themes[theme3]))

        if (!is.null(group3)) {
          if ("All" %in% groupuse3) {
            groupuse3 <- unique(srt_tmp[[group3, drop = TRUE]])
          }
          cells <- colnames(srt_tmp)[srt_tmp@meta.data[[group3]] %in% groupuse3]
        } else {
          cells <- colnames(srt_tmp)
        }

        if (is.null(aspect.ratio)) {
          aspect.ratio <- ifelse(is.null(group3), 5, 5 / max(length(unique(srt_tmp@meta.data[cells, group3])), 1))
        }

        # print(">>> plot:")
        # print(system.time(
        p3 <- SCP::CellStatPlot(
          srt = srt_tmp, stat.by = stat3, group.by = group3, split.by = split3, cells = cells,
          plot_type = plottype3, stat_type = stattype3, position = position3,
          label = label3, label.size = labelsize3, flip = flip3, palette = palette3, theme_use = theme3,
          aspect.ratio = as.numeric(aspect.ratio), # must be class of numeric instead of integer
          ncol = ncol3, byrow = byrow3, force = TRUE
        )
        # ))

        # print(">>> panel_fix:")
        # print(system.time(
        if (flip3) {
          p3 <- SCP::panel_fix(SCP::slim_data(p3), width = size3, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        } else {
          p3 <- SCP::panel_fix(SCP::slim_data(p3), height = size3, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        }
        # ))
        attr(p3, "dpi") <- 300
        return(p3)
      },
      seed = TRUE
    )
  }) %>%
    bindCache(
      input$dataset3, input$group3, input$split3, input$stat3,
      input$plottype3, input$stattype3, input$position3,
      input$label3, input$flip3, input$palette3, input$theme3,
      input$aspectratio3, input$aspectratio_value3,
      input$labelsize3, input$size3, input$ncol3, input$arrange3
    ) %>%
    bindEvent(input$submit3, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    prog <- Progress$new(min = 1, max = 10)
    prog$set(value = 3, message = "Fetch data...", detail = "[Cell statistical plot]")
    r3()$then(function(x) {
      promisedData[["p3"]] <- x
      width <- get_attr(x, "width")
      height <- get_attr(x, "height")
      dpi <- get_attr(x, "dpi")

      prog$set(value = 8, message = "Render plot...", detail = "[Cell statistical plot]")
      # print("renderPlot:")
      # print(system.time(
      output$plot3 <- renderUI({
        renderPlot(
          {
            x
          },
          width = width * 96,
          height = height * 96,
          res = 96
        )
      })
      # ))
    }) %>%
      finally(~ {
        prog$set(value = 10, message = "Done.", detail = "[Cell statistical plot]")
        prog$close()
      })
  }) %>% bindEvent(input$submit3, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$download3 <- downloadHandler(
    filename = function() {
      paste0("CellStatPlot-", format(Sys.time(), "%Y%m%d%H%M%S"), ".zip")
    },
    content = function(file) {
      width <- get_attr(promisedData[["p3"]], "width")
      height <- get_attr(promisedData[["p3"]], "height")
      dpi <- get_attr(promisedData[["p3"]], "dpi")

      temp1 <- tempfile(pattern = "CellStatPlot-", fileext = ".png")
      ggplot2::ggsave(filename = temp1, plot = promisedData[["p3"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      temp2 <- tempfile(pattern = "CellStatPlot-", fileext = ".pdf")
      ggplot2::ggsave(filename = temp2, plot = promisedData[["p3"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      zip(zipfile = file, flags = "-jq", files = c(temp1, temp2))
    },
    contentType = "application/zip"
  )

  # submit4  ----------------------------------------------------------------
  r4 <- reactive({
    dataset4 <- input$dataset4
    group4 <- input$group4
    groupuse4 <- input$groupuse4
    if (input$group4 == "None") {
      group4 <- NULL
    } else {
      group4 <- input$group4
    }
    if (input$split4 == "None") {
      split4 <- NULL
    } else {
      split4 <- input$split4
    }
    assays4 <- input$assays4
    slots4 <- input$slots4
    plottype4 <- input$plottype4
    features4 <- input$features4
    feature_area4 <- input$feature_area4
    plotby4 <- input$plotby4
    fillby4 <- input$fillby4
    coExp4 <- input$coExp4
    stack4 <- input$stack4
    flip4 <- input$flip4
    addbox4 <- input$addbox4
    addpoint4 <- input$addpoint4
    addtrend4 <- input$addtrend4
    palette4 <- input$palette4
    theme4 <- input$theme4
    sameylims4 <- input$sameylims4
    size4 <- input$size4
    ncol4 <- input$ncol4
    byrow4 <- input$arrange4
    if (input$aspectratio4 == "auto") {
      aspect.ratio <- NULL
    } else {
      aspect.ratio <- input$aspectratio_value4
    }

    # lapply(grep("4$",names(input),value = TRUE), function(x)print(paste0(x,":",input[[x]])))

    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", dataset4, "/", assays4, "/", slots4))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", dataset4, "/metadata.stat/asfeatures"))

    feature_area4 <- gsub(x = unlist(strsplit(feature_area4, "(\\r)|(\\n)", perl = TRUE)), pattern = " ", replacement = "")
    features4 <- c(as.character(features4), as.character(feature_area4))
    features4 <- unique(features4[features4 %in% c(all_features, meta_features_name)])
    if (length(features4) == 0) {
      features4 <- intersect(c(initial_feature, meta_features_name[1]), c(all_features, meta_features_name))[1]
    }

    promisedData[["p4"]] <- NULL
    promises::future_promise(
      {
        # print("******************************** New task ********************************")
        # print(">>> fetch data:")
        # print(system.time(
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset4,
          features = features4, slot = slots4, assay = assays4,
          metanames = unique(c(group4, split4))
        )
        # ))

        theme4 <- get(theme4, envir = asNamespace(themes[theme4]))

        if (!is.null(group4)) {
          if ("All" %in% groupuse4) {
            groupuse4 <- unique(srt_tmp[[group4, drop = TRUE]])
          }
          cells <- colnames(srt_tmp)[srt_tmp@meta.data[[group4]] %in% groupuse4]
        } else {
          cells <- colnames(srt_tmp)
        }

        if (is.null(aspect.ratio)) {
          aspect.ratio <- ifelse(is.null(group4), 5, 5 / max(length(unique(srt_tmp@meta.data[cells, group4])), 1))
        }

        # print(">>> plot:")
        # print(system.time(
        p4 <- SCP::FeatureStatPlot(
          srt = srt_tmp, stat.by = features4, group.by = group4, split.by = split4, cells = cells, slot = "data", plot_type = plottype4,
          calculate_coexp = coExp4, stack = stack4, flip = flip4,
          add_box = addbox4, add_point = addpoint4, add_trend = addtrend4,
          plot.by = plotby4, fill.by = fillby4, palette = palette4, theme_use = theme4, same.y.lims = sameylims4,
          aspect.ratio = as.numeric(aspect.ratio), # must be class of numeric instead of integer
          ncol = ncol4, byrow = byrow4, force = TRUE
        )
        # ))

        # print(">>> panel_fix:")
        # print(system.time(
        if (flip4) {
          p4 <- SCP::panel_fix(SCP::slim_data(p4), width = size4, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        } else {
          p4 <- SCP::panel_fix(SCP::slim_data(p4), height = size4, units = "in", raster = panel_raster, BPPARAM = BPPARAM, verbose = FALSE)
        }

        # ))
        attr(p4, "dpi") <- 300
        return(p4)
      },
      seed = TRUE
    )
  }) %>%
    bindCache(
      input$dataset4, input$group4, input$groupuse4, input$split4,
      input$assays4, input$slots4, input$plottype4,
      input$features4, input$feature_area4,
      input$coExp4, input$stack4, input$flip4,
      input$addbox4, input$addpoint4, input$addtrend4,
      input$fillby4, input$palette4, input$theme4,
      input$aspectratio4, input$aspectratio_value4,
      input$sameylims4, input$size4, input$ncol4, input$arrange4
    ) %>%
    bindEvent(input$submit4, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    prog <- Progress$new(min = 1, max = 10)
    prog$set(value = 3, message = "Fetch data...", detail = "[Feature statistical plot]")
    r4()$then(function(x) {
      promisedData[["p4"]] <- x
      width <- get_attr(x, "width")
      height <- get_attr(x, "height")
      dpi <- get_attr(x, "dpi")

      prog$set(value = 8, message = "Render plot...", detail = "[Feature statistical plot]")
      # print("renderPlot:")
      # print(system.time(
      output$plot4 <- renderUI({
        renderPlot(
          {
            x
          },
          width = width * 96,
          height = height * 96,
          res = 96
        )
      })
      # ))
    }) %>%
      finally(~ {
        prog$set(value = 10, message = "Done.", detail = "[Feature statistical plot]")
        prog$close()
      })
  }) %>% bindEvent(input$submit4, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$download4 <- downloadHandler(
    filename = function() {
      paste0("FeatureStatPlot-", format(Sys.time(), "%Y%m%d%H%M%S"), ".zip")
    },
    content = function(file) {
      width <- get_attr(promisedData[["p4"]], "width")
      height <- get_attr(promisedData[["p4"]], "height")
      dpi <- get_attr(promisedData[["p4"]], "dpi")

      temp1 <- tempfile(pattern = "FeatureStatPlot-", fileext = ".png")
      ggplot2::ggsave(filename = temp1, plot = promisedData[["p4"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      temp2 <- tempfile(pattern = "FeatureStatPlot-", fileext = ".pdf")
      ggplot2::ggsave(filename = temp2, plot = promisedData[["p4"]], width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)

      zip(zipfile = file, flags = "-jq", files = c(temp1, temp2))
    },
    contentType = "application/zip"
  )
}
  '

  main_code <- readLines(textConnection(main_code))
  args <- mget(names(formals()))
  args <- args[!names(args) %in% c("base_dir", "create_script", "style_script", "overwrite", "return_app")]
  args_code <- NULL
  for (varnm in names(args)) {
    args_code <- c(args_code, paste0(varnm, "=", deparse(args[[varnm]])))
  }
  app_code <- c(
    "# !/usr/bin/env Rscript",
    "if (!requireNamespace('SCP', quietly = TRUE)) {
      if (!requireNamespace('devtools', quietly = TRUE)) {install.packages('devtools')}
      devtools::install_github('zhanghao-njmu/SCP')
    }",
    "options(SCP_virtualenv_init = FALSE)",
    paste0("app_SCP_version <- package_version('", as.character(packageVersion("SCP")), "')"),
    paste0("if (utils::packageVersion('SCP') < app_SCP_version) {
      stop(paste0('SCExplorer requires SCP >= ", as.character(packageVersion("SCP")), "'))
    }"),
    "SCP::check_R(c('rhdf5', 'HDF5Array', 'shiny@1.6.0', 'ggplot2', 'ragg', 'htmlwidgets', 'plotly', 'bslib', 'future', 'promises', 'BiocParallel'))",
    "library(shiny)",
    "library(bslib)",
    "library(future)",
    "library(promises)",
    "library(BiocParallel)",
    "library(ggplot2)",
    "library(rlang)",
    args_code,
    "plan(multisession, workers = session_workers)",
    "if (.Platform$OS.type == 'windows') {
      BPPARAM = SerialParam()
    } else {
      BPPARAM = MulticoreParam(workers = plotting_workers)
    }",
    "page_theme <- bs_theme(bootswatch = 'zephyr')",
    main_code,
    "shinyApp(ui = ui, server = server)"
  )
  temp <- tempfile("SCExplorer")
  writeLines(app_code, temp)
  if (isTRUE(create_script)) {
    app_file <- paste0(base_dir, "/app.R")
    if (!file.exists(app_file) || isTRUE(overwrite)) {
      message("Create the SCExplorer app script: ", app_file)
      suppressWarnings(file.remove(app_file))
      file.copy(from = temp, to = app_file, overwrite = TRUE)
      if (isTRUE(style_script)) {
        message("Styling the script...")
        invisible(capture.output(styler::style_file(app_file)))
      }
    } else {
      message("app.R already exists. You may regenerate it with 'overwrite = TRUE'.")
    }
  }
  unlink(temp)

  if (isTRUE(return_app)) {
    app <- shiny::shinyAppDir(base_dir)
    return(app)
  } else {
    return(invisible(NULL))
  }
}
