#' CreateDataFile
#'
#' @param srt
#'
#' @param DataFile
#' @param name
#' @param assays
#' @param slots
#' @param compression_level
#' @param overwrite
#'
#' @importFrom Seurat GetAssayData
#' @importFrom SeuratObject as.sparse
#' @export
CreateDataFile <- function(srt, DataFile, name = NULL, assays = "RNA", slots = "data", compression_level = 6, overwrite = FALSE) {
  check_R("HDF5Array", "rhdf5")
  if (missing(DataFile) || is.null(DataFile)) {
    DataFile <- "Data.hdf5"
  }
  if (!file.exists(DataFile)) {
    rhdf5::h5createFile(DataFile)
  }
  if (is.null(name)) {
    name <- srt@project.name
    message("Set the dataset name to ", name)
  }
  if (substr(name, 1, 1) != "/") {
    name <- paste0("/", name)
  }
  if (!name %in% rhdf5::h5ls(DataFile)$group) {
    rhdf5::h5createGroup(file = DataFile, group = name)
  }

  message("Write the expression matrix to hdf5 file: ", DataFile)
  for (assay in assays) {
    for (slot in slots) {
      data <- t(GetAssayData(srt, slot = slot, assay = assay))
      if (isTRUE(overwrite)) {
        try(rhdf5::h5delete(file = DataFile, name = paste0(name, "/", assay, "/", slot)), silent = TRUE)
      }
      if (paste0(name, "/", assay, "/", slot) %in% rhdf5::h5ls(DataFile)$group) {
        message("Group ", paste0(name, "/", assay, "/", slot), " already exists in the ", DataFile)
      } else {
        if (!paste0(name, "/", assay) %in% rhdf5::h5ls(DataFile)$group) {
          rhdf5::h5createGroup(file = DataFile, group = paste0(name, "/", assay))
        }
        if (!inherits(data, "dgCMatrix")) {
          data <- as.sparse(data[1:nrow(data), ])
        }
        HDF5Array::writeTENxMatrix(x = data, filepath = DataFile, group = paste0(name, "/", assay, "/", slot), level = compression_level)
      }
    }
  }
  if (isTRUE(overwrite)) {
    try(rhdf5::h5delete(file = DataFile, name = paste0(name, "/Default_assay")), silent = TRUE)
  }
  if (paste0(name, "/Default_assay") %in% rhdf5::h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/Default_assay"), " already exists in the ", DataFile)
  } else {
    rhdf5::h5write(obj = DefaultAssay(srt), file = DataFile, name = paste0(name, "/Default_assay"), level = compression_level)
  }
  if (isTRUE(overwrite)) {
    try(rhdf5::h5delete(file = DataFile, name = paste0(name, "/cells")), silent = TRUE)
  }
  if (paste0(name, "/cells") %in% rhdf5::h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/cells"), " already exists in the ", DataFile)
  } else {
    rhdf5::h5write(obj = colnames(srt), file = DataFile, name = paste0(name, "/cells"), level = compression_level)
  }
  if (isTRUE(overwrite)) {
    try(rhdf5::h5delete(file = DataFile, name = paste0(name, "/features")), silent = TRUE)
  }
  if (paste0(name, "/features") %in% rhdf5::h5ls(DataFile)$group) {
    message("Group ", paste0(name, "/features"), " already exists in the ", DataFile)
  } else {
    rhdf5::h5write(obj = unique(unlist(lapply(srt@assays, rownames))), file = DataFile, name = paste0(name, "/features"), level = compression_level)
  }
  return(invisible(NULL))
}

#' CreateMetaFile
#'
#' @param srt
#'
#' @param MetaFile
#' @param name
#' @param write_tools
#' @param write_misc
#' @param ignore_nlevel
#' @param compression_level
#' @param overwrite
#'
#' @importFrom Seurat Key Reductions Embeddings
#' @importFrom Matrix sparseMatrix
#' @export
CreateMetaFile <- function(srt, MetaFile, name = NULL, write_tools = FALSE, write_misc = FALSE, ignore_nlevel = 100, compression_level = 6, overwrite = FALSE) {
  check_R("rhdf5")
  if (missing(MetaFile) || is.null(MetaFile)) {
    MetaFile <- "Meta.hdf5"
  }
  if (!file.exists(MetaFile)) {
    rhdf5::h5createFile(MetaFile)
  }
  if (is.null(name)) {
    name <- srt@project.name
    message("Set the dataset name to ", name)
  }
  if (substr(name, 1, 1) != "/") {
    name <- paste0("/", name)
  }
  if (!name %in% rhdf5::h5ls(MetaFile)$group) {
    rhdf5::h5createGroup(file = MetaFile, group = name)
  }

  message("Write the meta information to hdf5 file: ", MetaFile)
  if (!paste0(name, "/metadata") %in% rhdf5::h5ls(MetaFile)$group) {
    rhdf5::h5createGroup(file = MetaFile, group = paste0(name, "/metadata"))
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
    if (inherits(meta, "numeric")) {
      meta <- as.double(meta)
      meta_asfeatures <- c(meta_asfeatures, var)
    }
    if (!is.numeric(meta)) {
      if (length(unique(meta)) > ignore_nlevel) {
        warning("The number of categories in ", var, " is greater than ", ignore_nlevel, ", it will be ignored.", immediate. = TRUE)
      } else {
        meta_asgroups <- c(meta_asgroups, var)
      }
    }
    if (isTRUE(overwrite)) {
      try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/metadata/", var)), silent = TRUE)
    }
    if (paste0(name, "/metadata/", var) %in% paste0(rhdf5::h5ls(MetaFile)$group, "/", rhdf5::h5ls(MetaFile)$name)) {
      message("Group ", paste0(name, "/metadata/", var), " already exists in the ", MetaFile)
    } else {
      rhdf5::h5write(obj = meta, file = MetaFile, name = paste0(name, "/metadata/", var), write.attributes = write.attributes, level = compression_level)
    }
  }

  if (isTRUE(overwrite)) {
    try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/metadata.stat")), silent = TRUE)
  }
  if (paste0(name, "/metadata.stat") %in% rhdf5::h5ls(MetaFile)$group) {
    message("Group ", paste0(name, "/metadata.stat"), " already exists in the ", MetaFile)
  } else {
    rhdf5::h5createGroup(file = MetaFile, group = paste0(name, "/metadata.stat"))
    rhdf5::h5write(obj = meta_asfeatures, file = MetaFile, name = paste0(name, "/metadata.stat/asfeatures"), level = compression_level)
    rhdf5::h5write(obj = meta_asgroups, file = MetaFile, name = paste0(name, "/metadata.stat/asgroups"), level = compression_level)
  }

  if (!paste0(name, "/reductions") %in% rhdf5::h5ls(MetaFile)$group) {
    rhdf5::h5createGroup(file = MetaFile, group = paste0(name, "/reductions"))
  }
  for (reduction in Reductions(srt)) {
    emb <- Embeddings(srt, reduction)[colnames(srt), ]
    attr(emb, "key") <- as.character(Key(srt[[reduction]]))
    if (isTRUE(overwrite)) {
      try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/reductions/", reduction)), silent = TRUE)
    }
    if (paste0(name, "/reductions/", reduction) %in% paste0(rhdf5::h5ls(MetaFile)$group, "/", rhdf5::h5ls(MetaFile)$name)) {
      message("Group ", paste0(name, "/reductions/", reduction), " already exists in the ", MetaFile)
    } else {
      rhdf5::h5createDataset(
        file = MetaFile, dataset = paste0(name, "/reductions/", reduction),
        dims = dim(emb), chunk = c(nrow(emb), min(10, ncol(emb))), level = compression_level
      )
      suppressWarnings(rhdf5::h5write(obj = emb, file = MetaFile, name = paste0(name, "/reductions/", reduction), write.attributes = TRUE, level = compression_level))
    }
  }
  if (isTRUE(overwrite)) {
    try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/reductions.stat")), silent = TRUE)
  }
  if (paste0(name, "/reductions.stat") %in% rhdf5::h5ls(MetaFile)$group) {
    message("Group ", paste0(name, "/reductions.stat"), " already exists in the ", MetaFile)
  } else {
    rhdf5::h5createGroup(file = MetaFile, group = paste0(name, "/reductions.stat"))
    rhdf5::h5write(obj = DefaultReduction(srt), file = MetaFile, name = paste0(name, "/reductions.stat/Default_reduction"), level = compression_level)
  }

  if (isTRUE(write_misc)) {
    if (isTRUE(overwrite)) {
      try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/misc")), silent = TRUE)
    }
    if (paste0(name, "/misc") %in% rhdf5::h5ls(MetaFile)$group) {
      message("Group ", paste0(name, "/misc"), " already exists in the ", MetaFile)
    } else {
      rhdf5::h5write(obj = srt@misc, file = MetaFile, name = paste0(name, "/misc"), level = compression_level)
    }
  }
  if (isTRUE(write_tools)) {
    if (isTRUE(overwrite)) {
      try(rhdf5::h5delete(file = MetaFile, name = paste0(name, "/tools")), silent = TRUE)
    }
    if (paste0(name, "/tools") %in% rhdf5::h5ls(MetaFile)$group) {
      message("Group ", paste0(name, "/tools"), " already exists in the ", MetaFile)
    } else {
      rhdf5::h5write(obj = srt@tools, file = MetaFile, name = paste0(name, "/tools"), level = compression_level)
    }
  }

  return(invisible(NULL))
}

#' Prepare Seurat object data for the SCExplorer into hdf5 file
#'
#' @param object
#' @param base_dir
#' @param DataFile
#' @param MetaFile
#' @param assays
#' @param slots
#' @param ignore_nlevel
#' @param write_tools
#' @param write_misc
#' @param compression_level
#' @param overwrite
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
  base_dir <- R.utils::getAbsolutePath(base_dir)
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
#' @param DataFile
#' @param MetaFile
#' @param name
#' @param features
#' @param slot
#' @param assay
#' @param metanames
#' @param reduction
#'
#' @examples
#' \dontrun{
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' PrepareSCExplorer(pancreas_sub, base_dir = "./SCExplorer")
#' srt <- FetchH5(DataFile = "./SCExplorer/Data.hdf5", MetaFile = "./SCExplorer/Meta.hdf5", features = c("Ins1", "Ghrl"), metanames = c("SubCellType", "Phase"), reduction = "UMAP")
#' CellDimPlot(srt, group.by = c("SubCellType", "Phase"), reduction = "UMAP")
#' FeatureDimPlot(srt, features = c("Ins1", "Ghrl"), reduction = "UMAP")
#' }
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject
#' @importFrom Matrix sparseMatrix
#' @importFrom rlang %||%
#' @export
FetchH5 <- function(DataFile, MetaFile, name = NULL,
                    features = NULL, slot = NULL, assay = NULL,
                    metanames = NULL,
                    reduction = NULL) {
  check_R(c("HDF5Array", "rhdf5"))
  if (missing(DataFile) || missing(MetaFile)) {
    stop("'DataFile', 'MetaFile' must be provided.")
  }
  data_group <- rhdf5::h5ls(DataFile)$group
  meta_group <- rhdf5::h5ls(MetaFile)$group

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

  all_features <- rhdf5::h5read(DataFile, name = paste0(name, "/features"))
  all_cells <- rhdf5::h5read(DataFile, name = paste0(name, "/cells"))
  meta_struc <- rhdf5::h5ls(MetaFile)
  meta_features_name <- rhdf5::h5read(MetaFile, name = paste0(name, "/metadata.stat/asfeatures"))
  meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0(name, "/metadata.stat/asgroups"))
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
      assay <- as.character(rhdf5::h5read(DataFile, name = paste0(name, "/Default_assay")))
    }
    if (is.null(slot)) {
      slots <- as.character(unique(na.omit(sapply(strsplit(data_group[grep(name, data_group)], "/"), function(x) x[4]))))
      slot <- ifelse("counts" %in% slots, "counts", slots[1])
    }
    if (!paste0(name, "/", assay, "/", slot) %in% rhdf5::h5ls(DataFile)[["group"]]) {
      stop("There is no ", paste0(name, "/", assay, "/", slot), " in DataFile, please write it first using the PrepareSCExplorer function")
    }
    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0(name, "/", assay, "/", slot))
    gene_features <- gene_features[gene_features %in% colnames(data)]
  }

  if (length(gene_features) > 0) {
    srt_tmp <- CreateSeuratObject(assay = assay %||% "RNA", counts = as(t(data[, gene_features, drop = FALSE]), "matrix")) # sparseMatrix
  } else {
    srt_tmp <- CreateSeuratObject(assay = assay %||% "RNA", counts = matrix(data = 0, ncol = length(all_cells), dimnames = list("empty", all_cells)))
  }

  if (length(c(metanames, meta_features)) > 0) {
    for (i in unique(c(metanames, meta_features))) {
      meta <- rhdf5::h5read(MetaFile, name = paste0(name, "/metadata/", i))
      if (is.array(meta)) {
        if (is.integer(meta)) {
          meta <- as.integer(meta)
        } else if (is.numeric(meta)) {
          meta <- as.numeric(meta)
        } else if (is.character(meta)) {
          meta <- as.character(meta)
        }
      }
      meta_attr <- rhdf5::h5readAttributes(MetaFile, name = paste0(name, "/metadata/", i))
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
      reduction <- as.matrix(rhdf5::h5read(MetaFile, name = paste0(name, "/reductions/", i)))
      reduction_attr <- rhdf5::h5readAttributes(MetaFile, name = paste0(name, "/reductions/", i))
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

#' RunSCExplorer
#'
#' @param base_dir
#' @param DataFile
#' @param MetaFile
#' @param title
#' @param initial_dataset
#' @param initial_reduction
#' @param initial_group
#' @param initial_feature
#' @param initial_assay
#' @param initial_slot
#' @param initial_label
#' @param initial_cell_palette
#' @param initial_feature_palette
#' @param initial_theme
#' @param initial_theme
#' @param initial_coExp
#' @param initial_size
#' @param initial_ncol
#' @param initial_arrange
#' @param initial_panel_dpi
#' @param initial_dpi
#' @param create_script
#' @param style_script
#' @param overwrite
#' @param return_app
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
#' app <- RunSCExplorer(base_dir = "./SCExplorer", overwrite = TRUE)
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
#' @importFrom rlang %||%
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
                          initial_label = "No",
                          initial_cell_palette = "Paired",
                          initial_feature_palette = "Spectral",
                          initial_theme = "theme_scp",
                          initial_coExp = "No",
                          initial_size = 4,
                          initial_ncol = 3,
                          initial_arrange = "Row",
                          initial_dpi = 100,
                          initial_raster = FALSE,
                          create_script = TRUE,
                          style_script = require("styler", quietly = TRUE),
                          overwrite = FALSE,
                          return_app = TRUE) {
  check_R(c("HDF5Array", "rhdf5", "shiny@1.6.0", "ragg", "bslib", "future", "promises"))
  DataFile_full <- paste0(base_dir, "/", DataFile)
  MetaFile_full <- paste0(base_dir, "/", MetaFile)
  if (!file.exists(DataFile_full) || !file.exists(MetaFile_full)) {
    stop("Please create the DataFile and MetaFile using PrepareSCExplorer function first!")
  }

  main_code <- '
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

palette_list <- SCP:::palette_list

ui <- fluidPage(
  theme = page_theme,
  navbarPage(
    title = title,
    # 1. Cell dimensional reduction plot ----------------------------------------------------------------------
    tabPanel(
      "Cell dimensional reduction plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset1",
            label = "Select a dataset",
            choices = group,
            selected = substr(initial_dataset, 2, nchar(initial_dataset))
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
          radioButtons(
            inputId = "label1",
            label = "Whether to add labels",
            choices = c("Yes", "No"),
            selected = initial_label,
            inline = TRUE
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
            choices = c("theme_scp", "theme_blank"),
            selected = initial_theme
          ),
          selectInput(
            inputId = "raster1",
            label = "Raster",
            choices = c("TRUE", "FALSE"),
            selected = initial_raster
          ),
          numericInput(
            inputId = "size1",
            label = "Panel size",
            value = initial_size,
            min = 1,
            max = 10,
            step = 0.1,
            width = "150px"
          ),
          numericInput(
            inputId = "plot_dpi1",
            label = "Resolution of the plot",
            value = initial_dpi,
            min = 50,
            max = 1000,
            step = 50,
            width = "150px"
          ),
          numericInput(
            inputId = "ncol1",
            label = "Number of columns",
            value = initial_ncol,
            min = 1,
            max = 100,
            step = 1,
            width = "150px"
          ),
          radioButtons(
            inputId = "arrange1",
            label = "Arrange by",
            choices = c("Row", "Column"),
            selected = initial_arrange,
            inline = TRUE
          ),
          actionButton(inputId = "submit1", label = "Submit", icon("play")),
          helpText("Click the submit button to update the plot displayed in the right panel.")
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
      "Feature dimensional reduction plot",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset2",
            label = "Select a dataset",
            choices = group,
            selected = substr(initial_dataset, 2, nchar(initial_dataset))
          ),
          selectInput(
            inputId = "reduction2",
            label = "Select a reduction",
            choices = reduction_name,
            selected = initial_reduction
          ),
          selectInput(
            inputId = "assays2",
            label = "Select a assay",
            choices = assays,
            selected = initial_assay
          ),
          selectInput(
            inputId = "slots2",
            label = "Select a slot",
            choices = slots,
            selected = initial_slot
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
            height = "200px",
            placeholder = paste(sample(all_features, 4), collapse = "\n")
          ),
          selectInput(
            inputId = "split2",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          radioButtons(
            inputId = "coExp2",
            label = "Calculate co-expression?",
            choices = c("Yes", "No"),
            selected = initial_coExp,
            inline = TRUE
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
            choices = c("theme_scp", "theme_blank"),
            selected = initial_theme
          ),
          selectInput(
            inputId = "raster2",
            label = "Raster",
            choices = c("TRUE", "FALSE"),
            selected = initial_raster
          ),
          numericInput(
            inputId = "size2",
            label = "Panel size",
            value = initial_size,
            min = 1,
            max = 10,
            step = 0.1,
            width = "150px"
          ),
          numericInput(
            inputId = "plot_dpi2",
            label = "Resolution of the plot",
            value = initial_dpi,
            min = 50,
            max = 1000,
            step = 50,
            width = "150px"
          ),
          numericInput(
            inputId = "ncol2",
            label = "Number of columns",
            value = initial_ncol,
            min = 1,
            max = 100,
            step = 1,
            width = "150px"
          ),
          radioButtons(
            inputId = "arrange2",
            label = "Arrange by",
            choices = c("Row", "Column"),
            selected = initial_arrange,
            inline = TRUE
          ),
          actionButton(inputId = "submit2", label = "Submit", icon("play")),
          helpText("Click the submit button to update the plot displayed in the right panel.")
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
    # 3. Statistical plot of cells ----------------------------------------------------------------------
    tabPanel(
      "Statistical plot of cells",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset3",
            label = "Select a dataset",
            choices = group,
            selected = substr(initial_dataset, 2, nchar(initial_dataset))
          ),
          selectInput(
            inputId = "plottype3",
            label = "Select a plot type",
            choices = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", "upset"),
            selected = "bar"
          ),
          selectInput(
            inputId = "stattype3",
            label = "Select a value type",
            choices = c("percent", "count"),
            selected = "stack"
          ),
          selectInput(
            inputId = "position3",
            label = "Select a position",
            choices = c("stack", "dodge"),
            selected = "stack"
          ),
          selectizeInput(
            inputId = "stat3",
            label = "Select a variable to be counted",
            choices = meta_groups_name,
            selected = initial_group,
            multiple = TRUE,
            options = list(maxOptions = 20, maxItems = 7)
          ),
          selectInput(
            inputId = "group3",
            label = "Select a variable to group cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          selectInput(
            inputId = "split3",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          radioButtons(
            inputId = "label3",
            label = "Whether to add labels",
            choices = c("Yes", "No"),
            selected = initial_label,
            inline = TRUE
          ),
          selectInput(
            inputId = "palette3",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_cell_palette
          ),
          numericInput(
            inputId = "size3",
            label = "Panel size",
            value = initial_size,
            min = 1,
            max = 10,
            step = 0.1,
            width = "150px"
          ),
          numericInput(
            inputId = "plot_dpi3",
            label = "Resolution of the plot",
            value = initial_dpi,
            min = 50,
            max = 1000,
            step = 50,
            width = "150px"
          ),
          numericInput(
            inputId = "ncol3",
            label = "Number of columns",
            value = initial_ncol,
            min = 1,
            max = 100,
            step = 1,
            width = "150px"
          ),
          radioButtons(
            inputId = "arrange3",
            label = "Arrange by",
            choices = c("Row", "Column"),
            selected = initial_arrange,
            inline = TRUE
          ),
          actionButton(inputId = "submit3", label = "Submit", icon("play")),
          helpText("Click the submit button to update the plot displayed in the right panel.")
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
    # 4. Statistical plot of features ----------------------------------------------------------------------
    tabPanel(
      "Statistical plot of features",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "dataset4",
            label = "Select a dataset",
            choices = group,
            selected = substr(initial_dataset, 2, nchar(initial_dataset))
          ),
          selectInput(
            inputId = "assays4",
            label = "Select a assay",
            choices = assays,
            selected = initial_assay
          ),
          selectInput(
            inputId = "slots4",
            label = "Select a slot",
            choices = slots,
            selected = initial_slot
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
            height = "200px",
            placeholder = paste(sample(all_features, 4), collapse = "\n")
          ),
          selectInput(
            inputId = "group4",
            label = "Select a variable to group cells",
            choices = c("None", meta_groups_name),
            selected = initial_group
          ),
          selectInput(
            inputId = "split4",
            label = "Select a variable to split cells",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          radioButtons(
            inputId = "coExp4",
            label = "Calculate co-expression?",
            choices = c("Yes", "No"),
            selected = initial_coExp,
            inline = TRUE
          ),
          selectInput(
            inputId = "palette4",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_cell_palette
          ),
          numericInput(
            inputId = "size4",
            label = "Panel size",
            value = initial_size,
            min = 1,
            max = 10,
            step = 0.1,
            width = "150px"
          ),
          numericInput(
            inputId = "plot_dpi4",
            label = "Resolution of the plot",
            value = initial_dpi,
            min = 50,
            max = 1000,
            step = 50,
            width = "150px"
          ),
          numericInput(
            inputId = "ncol4",
            label = "Number of columns",
            value = initial_ncol,
            min = 1,
            max = 100,
            step = 1,
            width = "150px"
          ),
          radioButtons(
            inputId = "arrange4",
            label = "Arrange by",
            choices = c("Row", "Column"),
            selected = initial_arrange,
            inline = TRUE
          ),
          actionButton(inputId = "submit4", label = "Submit", icon("play")),
          helpText("Click the submit button to update the plot displayed in the right panel.")
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
    dpi <- attr(x, "dpi")
    if (verbose) {
      message(paste("width:", width, "height:", height, "dpi:", dpi))
    }
    if (attr == "dpi") {
      return(dpi)
    }
    if (attr == "width") {
      return(width * dpi)
    }
    if (attr == "height") {
      return(height * dpi)
    }
  }

  updateSelectizeInput(session, "features2", choices = c(meta_features_name, all_features), selected = initial_feature, server = TRUE)
  updateSelectizeInput(session, "features4", choices = c(meta_features_name, all_features), selected = initial_feature, server = TRUE)

  # change dataset  ----------------------------------------------------------------
  observe({
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset1, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/reductions.stat/Default_reduction")))
    updateSelectInput(session, "reduction1", choices = reduction_name, selected = default_reduction)
    updateSelectInput(session, "group1", choices = meta_groups_name, selected = "orig.ident")
    updateSelectInput(session, "split1", choices = c("None", meta_groups_name), selected = "None")
  }) %>% bindEvent(input$dataset1, ignoreNULL = TRUE, ignoreInit = TRUE)

  observe({
    assays <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[3])))
    slots <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[4])))
    default_assay <- as.character(rhdf5::h5read(DataFile, name = paste0("/", input$dataset2, "/Default_assay")))
    default_slot <- ifelse("data" %in% slots, "data", slots[1])
    updateSelectInput(session, "assays2", choices = assays, selected = default_assay)
    updateSelectInput(session, "slots2", choices = slots, selected = default_slot)
    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", input$dataset2, "/", default_assay, "/", default_slot))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/metadata.stat/asfeatures"))
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset2, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/reductions.stat/Default_reduction")))
    updateSelectInput(session, "reduction2", choices = reduction_name, selected = default_reduction)
    updateSelectizeInput(session, "features2",
      choices = c(meta_features_name, all_features), selected = meta_features_name[1],
      options = list(maxOptions = 20, maxItems = 20), server = TRUE
    )
    updateSelectInput(session, "split2", choices = c("None", meta_groups_name), selected = "None")
    updateSelectInput(session, "group2", choices = meta_groups_name, selected = "orig.ident")
  }) %>% bindEvent(input$dataset2, ignoreNULL = TRUE, ignoreInit = TRUE)

  observe({
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset3, "/metadata.stat/asgroups"))
    updateSelectizeInput(session, "stat3", choices = meta_groups_name, selected = "orig.ident", options = list(maxOptions = 20, maxItems = 7))
    updateSelectInput(session, "group3", choices = meta_groups_name, selected = "orig.ident")
    updateSelectInput(session, "split3", choices = c("None", meta_groups_name), selected = "None")
  }) %>% bindEvent(input$dataset3, ignoreNULL = TRUE, ignoreInit = TRUE)

  observe({
    assays <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset4, data_group)], "/"), function(x) x[3])))
    slots <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset4, data_group)], "/"), function(x) x[4])))
    default_assay <- as.character(rhdf5::h5read(DataFile, name = paste0("/", input$dataset4, "/Default_assay")))
    default_slot <- ifelse("data" %in% slots, "data", slots[1])
    updateSelectInput(session, "assays2", choices = assays, selected = default_assay)
    updateSelectInput(session, "slots2", choices = slots, selected = default_slot)
    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", input$dataset4, "/", default_assay, "/", default_slot))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset4, "/metadata.stat/asfeatures"))
    meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset4, "/metadata.stat/asgroups"))
    updateSelectizeInput(session, "features2",
      choices = c(meta_features_name, all_features), selected = meta_features_name[1],
      options = list(maxOptions = 20, maxItems = 20), server = TRUE
    )
    updateSelectInput(session, "split2", choices = c("None", meta_groups_name), selected = "None")
    updateSelectInput(session, "group2", choices = meta_groups_name, selected = "orig.ident")
  }) %>% bindEvent(input$dataset4, ignoreNULL = TRUE, ignoreInit = TRUE)

  # submit1  ----------------------------------------------------------------
  r1 <- reactive({
    dataset1 <- input$dataset1
    group1 <- input$group1
    if (is.null(group1)) {
      message("input group1 is null")
      group1 <- "orig.ident"
    }
    if (input$split1 == "None") {
      split1 <- NULL
    } else {
      split1 <- input$split1
    }
    reduction1 <- input$reduction1
    label1 <- input$label1
    raster1 <- input$raster1 == "TRUE"
    palette1 <- input$palette1
    theme1 <- input$theme1
    size1 <- input$size1
    plot_dpi1 <- input$plot_dpi1
    ncol1 <- input$ncol1
    arrange1 <- input$arrange1

    # message("dataset1:", dataset1)
    # message("group1:", group1)
    # message("split1:", split1)
    # message("reduction1:", reduction1)
    # message("label1:", label1)
    # message("raster1:", raster1)
    # message("palette1:", palette1)
    # message("theme1:", theme1)
    # message("size1:", size1)
    # message("plot_dpi1:", plot_dpi1)
    # message("ncol1:", ncol1)
    # message("arrange1:", arrange1)

    # message(paste0("run r1: ", group1))

    promisedData[["p1_dim"]] <- NULL
    promisedData[["p1_3d"]] <- NULL
    promises::future_promise(
      {
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset1,
          metanames = unique(c(group1, split1)), reduction = reduction1
        )
        p1_dim <- SCP::CellDimPlot(srt_tmp,
          group.by = group1, split.by = split1, reduction = reduction1, raster = raster1,
          label = ifelse(label1 == "Yes", TRUE, FALSE), palette = palette1, theme_use = theme1,
          ncol = ncol1, byrow = ifelse(arrange1 == "Row", TRUE, FALSE), force = TRUE
        )
        p1_dim <- SCP::panel_fix(p1_dim, height = size1, raster = FALSE, verbose = FALSE)
        attr(p1_dim, "dpi") <- plot_dpi1
        plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
        if (isTRUE(plot3d)) {
          p1_3d <- SCP::CellDimPlot3D(srt_tmp, group.by = group1, reduction = reduction1, palette = palette1, force = TRUE)
        } else {
          p1_3d <- NULL
        }
        return(list(p1_dim, p1_3d))
      },
      seed = TRUE
    )
  }) %>%
    bindCache(input$dataset1, input$group1, input$split1, input$reduction1, input$label1, input$raster1, input$palette1, input$theme1, input$size1, input$plot_dpi1, input$ncol1, input$arrange1) %>%
    bindEvent(input$submit1, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    r1()$then(function(x) {
      promisedData[["p1_dim"]] <- x[[1]]
      promisedData[["p1_3d"]] <- x[[2]]
    })
  }) %>% bindEvent(input$submit1, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$plot1 <- renderUI({
    renderPlot(
      {
        req(promisedData[["p1_dim"]])
      },
      width = get_attr(req(promisedData[["p1_dim"]]), "width"),
      height = get_attr(req(promisedData[["p1_dim"]]), "height"),
      res = get_attr(req(promisedData[["p1_dim"]]), "dpi")
    )
  })

  output$plot1_3d <- plotly::renderPlotly({
    req(promisedData[["p1_3d"]])
  })


  # submit2  ----------------------------------------------------------------
  r2 <- reactive({
    dataset2 <- input$dataset2
    assays2 <- input$assays2
    slots2 <- input$slots2
    features2 <- input$features2
    feature_area2 <- input$feature_area2
    coExp2 <- input$coExp2
    reduction2 <- input$reduction2
    raster2 <- input$raster2 == "TRUE"
    if (input$split2 == "None") {
      split2 <- NULL
    } else {
      split2 <- input$split2
    }
    palette2 <- input$palette2
    theme2 <- input$theme2
    size2 <- input$size2
    plot_dpi2 <- input$plot_dpi2
    ncol2 <- input$ncol2
    arrange2 <- input$arrange2

    # message("dataset2:", dataset2)
    # message("assays2:", assays2)
    # message("slots2:", slots2)
    # message("features2:", features2)
    # message("coExp2:", coExp2)
    # message("split2:", split2)
    # message("reduction2:", reduction2)
    # message("raster2:", raster2)
    # message("palette2:", palette2)
    # message("theme2:", theme2)
    # message("size2:", size2)
    # message("plot_dpi2:", plot_dpi2)
    # message("ncol2:", ncol2)
    # message("arrange2:", arrange2)

    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", dataset2, "/", assays2, "/", slots2))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", dataset2, "/metadata.stat/asfeatures"))

    if (is.null(features2)) {
      message("input feature2 is null")
      features2 <- meta_features_name[1]
    }
    feature_area <- gsub(x = unlist(strsplit(feature_area2, "(\\r)|(\\n)", perl = TRUE)), pattern = " ", replacement = "")
    features2 <- c(as.character(features2), as.character(feature_area))
    features2 <- unique(features2[features2 %in% c(all_features, meta_features_name)])

    promisedData[["p2_dim"]] <- NULL
    promisedData[["p2_3d"]] <- NULL
    promises::future_promise(
      {
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset2,
          features = features2, slot = slots2, assay = assays2,
          metanames = split2, reduction = reduction2
        )
        p2_dim <- SCP::FeatureDimPlot(
          srt = srt_tmp, features = features2, split.by = split2, reduction = reduction2, slot = "data", raster = raster2,
          calculate_coexp = ifelse(coExp2 == "Yes", TRUE, FALSE), palette = palette2, theme_use = theme2,
          ncol = ncol2, byrow = ifelse(arrange2 == "Row", TRUE, FALSE)
        )
        p2_dim <- SCP::panel_fix(p2_dim, height = size2, raster = FALSE, verbose = FALSE)
        attr(p2_dim, "dpi") <- plot_dpi2
        plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
        if (isTRUE(plot3d)) {
          p2_3d <- SCP::FeatureDimPlot3D(
            srt = srt_tmp, features = features2, reduction = reduction2,
            calculate_coexp = ifelse(coExp2 == "Yes", TRUE, FALSE)
          )
        } else {
          p2_3d <- NULL
        }
        return(list(p2_dim, p2_3d))
      },
      seed = TRUE
    )
  }) %>%
    bindCache(input$dataset2, input$assays2, input$slots2, input$features2, input$feature_area2, input$coExp2, input$reduction2, input$raster2, input$split2, input$palette2, input$theme2, input$size2, input$plot_dpi2, input$ncol2, input$arrange2) %>%
    bindEvent(input$submit2, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    r2()$then(function(x) {
      promisedData[["p2_dim"]] <- x[[1]]
      promisedData[["p2_3d"]] <- x[[2]]
    })
  }) %>% bindEvent(input$submit2, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$plot2 <- renderUI({
    renderPlot(
      {
        req(promisedData[["p2_dim"]])
      },
      width = get_attr(req(promisedData[["p2_dim"]]), "width"),
      height = get_attr(req(promisedData[["p2_dim"]]), "height"),
      res = get_attr(req(promisedData[["p2_dim"]]), "dpi")
    )
  })

  output$plot2_3d <- plotly::renderPlotly({
    req(promisedData[["p2_3d"]])
  })


  # submit3  ----------------------------------------------------------------
  r3 <- reactive({
    dataset3 <- input$dataset3
    plottype3 <- input$plottype3
    stattype3 <- input$stattype3
    position3 <- input$position3
    stat3 <- input$stat3
    group3 <- input$group3
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
    palette3 <- input$palette3
    size3 <- input$size3
    plot_dpi3 <- input$plot_dpi3
    ncol3 <- input$ncol3
    arrange3 <- input$arrange3

    # message("dataset3:", dataset3)
    # message("plottype3:", plottype3)
    # message("stattype3:", stattype3)
    # message("position3:", position3)
    # message("stat3:", stat3)
    # message("group3:", group3)
    # message("split3:", split3)
    # message("label3:", label3)
    # message("palette3:", palette3)
    # message("size3:", size3)
    # message("plot_dpi3:", plot_dpi3)
    # message("ncol3:", ncol3)
    # message("arrange3:", arrange3)

    promisedData[["p3"]] <- NULL
    promises::future_promise(
      {
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset3,
          metanames = unique(c(stat3, group3, split3))
        )
        p3 <- SCP::CellStatPlot(
          srt = srt_tmp, stat.by = stat3, group.by = group3, split.by = split3,
          plot_type = plottype3, stat_type = stattype3, position = position3,
          label = ifelse(label3 == "Yes", TRUE, FALSE), palette = palette3,
          aspect.ratio = 8 / max(length(unique(srt_tmp[[group3, drop = TRUE]])), 1),
          ncol = ncol3, byrow = ifelse(arrange3 == "Row", TRUE, FALSE), force = TRUE
        )
        p3 <- SCP::panel_fix(p3, height = size3, raster = FALSE, verbose = FALSE)
        attr(p3, "dpi") <- plot_dpi3
        return(p3)
      },
      seed = TRUE
    )
  }) %>%
    bindCache(input$dataset3, input$plottype3, input$stattype3, input$position3, input$stat3, input$group3, input$split3, input$label3, input$palette3, input$size3, input$plot_dpi3, input$ncol3, input$arrange3) %>%
    bindEvent(input$submit3, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    r3()$then(function(x) {
      promisedData[["p3"]] <- x
    })
  }) %>% bindEvent(input$submit3, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$plot3 <- renderUI({
    renderPlot(
      {
        req(promisedData[["p3"]])
      },
      width = get_attr(req(promisedData[["p3"]]), "width"),
      height = get_attr(req(promisedData[["p3"]]), "height"),
      res = get_attr(req(promisedData[["p3"]]), "dpi")
    )
  })

  # submit4  ----------------------------------------------------------------
  r4 <- reactive({
    dataset4 <- input$dataset4
    assays4 <- input$assays4
    slots4 <- input$slots4
    plottype4 <- input$plottype4
    features4 <- input$features4
    feature_area4 <- input$feature_area4
    coExp4 <- input$coExp4
    group4 <- input$group4
    if (input$split4 == "None") {
      split4 <- NULL
    } else {
      split4 <- input$split4
    }
    palette4 <- input$palette4
    size4 <- input$size4
    plot_dpi4 <- input$plot_dpi4
    ncol4 <- input$ncol4
    arrange4 <- input$arrange4

    data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", dataset4, "/", assays4, "/", slots4))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", dataset4, "/metadata.stat/asfeatures"))

    if (is.null(features4)) {
      message("input feature4 is null")
      features4 <- meta_features_name[1]
    }
    feature_area <- gsub(x = unlist(strsplit(feature_area4, "(\\r)|(\\n)", perl = TRUE)), pattern = " ", replacement = "")
    features4 <- c(as.character(features4), as.character(feature_area))
    features4 <- unique(features4[features4 %in% c(all_features, meta_features_name)])

    promisedData[["p4"]] <- NULL
    promises::future_promise(
      {
        srt_tmp <- SCP::FetchH5(
          DataFile = DataFile, MetaFile = MetaFile, name = dataset4,
          features = features4, slot = slots4, assay = assays4,
          metanames = unique(c(group4, split4))
        )
        p4 <- SCP::FeatureStatPlot(
          srt = srt_tmp, stat.by = features4, group.by = group4, split.by = split4, plot_type = plottype4,
          calculate_coexp = ifelse(coExp4 == "Yes", TRUE, FALSE), palette = palette4,
          aspect.ratio = 8 / max(length(unique(srt_tmp[[group4, drop = TRUE]])), 1),
          ncol = ncol4, byrow = ifelse(arrange4 == "Row", TRUE, FALSE), force = TRUE
        )
        p4 <- SCP::panel_fix(p4, height = size4, raster = FALSE, verbose = FALSE)
        attr(p4, "dpi") <- plot_dpi4
        return(p4)
      },
      seed = TRUE
    )
  }) %>%
    bindCache(input$dataset4, input$assays4, input$slots4, input$plottype4, input$features4, input$feature_area4, input$coExp4, input$group4, input$split4, input$palette4, input$size4, input$plot_dpi4, input$ncol4, input$arrange4) %>%
    bindEvent(input$submit4, ignoreNULL = FALSE, ignoreInit = FALSE)

  observe({
    r4()$then(function(x) {
      promisedData[["p4"]] <- x
    })
  }) %>% bindEvent(input$submit4, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$plot4 <- renderUI({
    renderPlot(
      {
        req(promisedData[["p4"]])
      },
      width = get_attr(req(promisedData[["p4"]]), "width"),
      height = get_attr(req(promisedData[["p4"]]), "height"),
      res = get_attr(req(promisedData[["p4"]]), "dpi")
    )
  })
}
  '

  main_code <- readLines(textConnection(main_code))
  main_code <- gsub("\\\\r", "\\\\\\\r", main_code)
  main_code <- gsub("\\\\n", "\\\\\\\n", main_code)
  args <- mget(names(formals()))
  args <- args[!names(args) %in% c("base_dir", "create_script", "style_script", "overwrite", "return_app")]
  for (varnm in names(args)) {
    main_code <- c(paste0(varnm, "=", deparse(args[[varnm]])), main_code)
  }
  app_code <- c(
    "# !/usr/bin/env Rscript",
    "if (!requireNamespace('SCP', quietly = TRUE)) {
      if (!requireNamespace('devtools', quietly = TRUE)) {install.packages('devtools')}
      devtools::install_github('zhanghao-njmu/SCP')
    }",
    "options(SCP_virtualenv_init = FALSE)",
    paste0("app_SCP_version <- package_version('", as.character(utils::packageVersion("SCP")), "')"),
    paste0("if (utils::packageVersion('SCP') < app_SCP_version) {
      stop(paste0('SCExplorer requires SCP >= ", as.character(utils::packageVersion("SCP")), "'))
    }"),
    "SCP::check_R(c('HDF5Array', 'rhdf5', 'shiny@1.6.0', 'ragg', 'bslib', 'future', 'promises'))",
    "library(shiny)",
    "library(bslib)",
    "library(future)",
    "library(promises)",
    "plan(multisession, workers = min(availableCores() - 1, 64))",
    "page_theme <- bs_theme(bootswatch = 'zephyr')",
    main_code,
    "shinyApp(ui = ui, server = server)"
  )
  temp <- tempfile("SCExplorer")
  writeLines(app_code, temp)
  wd <- getwd()
  setwd(base_dir)
  source(temp)
  setwd(wd)
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
