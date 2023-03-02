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
  check_R("HDF5Array")
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
      slot <- ifelse("data" %in% slots, "data", slots[1])
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
#' @param initial_palette1
#' @param initial_palette2
#' @param initial_theme1
#' @param initial_theme2
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
                          initial_palette1 = "Paired",
                          initial_palette2 = "Spectral",
                          initial_theme1 = "theme_scp",
                          initial_theme2 = "theme_scp",
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
  check_R(c("shiny", "shinycssloaders", "ragg"))
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

  initial_srt_tmp <- SCP::FetchH5(
    DataFile = DataFile, MetaFile = MetaFile, name = initial_dataset,
    features = initial_feature, slot = initial_slot, assay = initial_assay,
    metanames = initial_group, reduction = initial_reduction
  )

  initial_p1_dim <- SCP::CellDimPlot(
    srt = initial_srt_tmp, group.by = initial_group, reduction = initial_reduction, raster = initial_raster,
    label = ifelse(initial_label == "Yes", TRUE, FALSE), palette = initial_palette1, theme_use = initial_theme1,
    ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE), force = TRUE
  )
  initial_p1_dim <- SCP::panel_fix(initial_p1_dim, height = initial_size, raster = FALSE, verbose = FALSE)
  initial_p2_dim <- SCP::FeatureDimPlot(
    srt = initial_srt_tmp, features = initial_feature, reduction = initial_reduction, slot = "data", raster = initial_raster,
    calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE), palette = initial_palette2, theme_use = initial_theme2,
    ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE)
  )
  initial_p2_dim <- SCP::panel_fix(initial_p2_dim, height = initial_size, raster = FALSE, verbose = FALSE)
  initial_p2_vln <- SCP::FeatureStatPlot(
    srt = initial_srt_tmp, stat.by = initial_feature, group.by = initial_group,
    calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE), palette = initial_palette2,
    ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE), force = TRUE
  )
  initial_p2_vln <- SCP::panel_fix(initial_p2_vln, height = initial_size, raster = FALSE, verbose = FALSE)

  initial_plot3d <- max(sapply(names(initial_srt_tmp@reductions), function(r) dim(initial_srt_tmp[[r]])[2])) >= 3
  if (isTRUE(initial_plot3d)) {
    initial_p1_3d <- SCP::CellDimPlot3D(
      srt = initial_srt_tmp, group.by = initial_group, reduction = initial_reduction, palette = initial_palette1,
      force = TRUE
    )
    initial_p2_3d <- SCP::FeatureDimPlot3D(
      srt = initial_srt_tmp, features = initial_feature, reduction = initial_reduction,
      calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE)
    )
  } else {
    initial_p1_3d <- initial_p2_3d <- NULL
  }

  palette_list <- SCP:::palette_list
  # ui ----------------------------------------------------------------------
  ui <- fluidPage(
    navbarPage(
      title = title,
      tabPanel(
        "Cell group view",
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
            inputId = "class1",
            label = "Select a categorical variable",
            choices = meta_groups_name,
            selected = initial_group
          ),
          selectInput(
            inputId = "split1",
            label = "Select a variable for splitting",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          radioButtons(
            inputId = "label1",
            label = "Whether labeled",
            choices = c("Yes", "No"),
            selected = initial_label,
            inline = TRUE
          ),
          selectInput(
            inputId = "palette1",
            label = "Select a palette",
            choices = names(palette_list),
            selected = initial_palette1
          ),
          selectInput(
            inputId = "theme1",
            label = "Select a theme",
            choices = c("theme_scp", "theme_blank"),
            selected = initial_theme1
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
                    shinycssloaders::withSpinner(plotOutput("plot1", height = "100%", width = "100%"))
                  )
                )
              ),
              tabPanel(
                title = "3D plot",
                column(
                  width = 12, offset = 0, style = "padding:0px;margin:0%",
                  div(
                    style = "overflow-x: auto;",
                    shinycssloaders::withSpinner(plotly::plotlyOutput("plot1_3d", height = "100%", width = "100%"))
                  )
                )
              )
            )
          )
        )
      ),
      tabPanel(
        "Feature expression view",
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
            label = "Select genes",
            choices = NULL,
            selected = initial_feature,
            multiple = TRUE,
            options = list(maxOptions = 20, maxItems = 20)
          ),
          textAreaInput(
            inputId = "gene_area",
            label = "Input genes",
            height = "200px",
            placeholder = paste(sample(all_features, 4), collapse = "\n")
          ),
          selectInput(
            inputId = "split2",
            label = "Select a variable for splitting",
            choices = c("None", meta_groups_name),
            selected = "None"
          ),
          selectInput(
            inputId = "group2",
            label = "Select a grouping variable",
            choices = c("None", meta_groups_name),
            selected = initial_group
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
            selected = initial_palette2
          ),
          selectInput(
            inputId = "theme2",
            label = "Select a theme",
            choices = c("theme_scp", "theme_blank"),
            selected = initial_theme2
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
                  shinycssloaders::withSpinner(plotOutput("plot2", height = "100%", width = "100%"))
                )
              )
            ),
            tabPanel(
              title = "3D plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  shinycssloaders::withSpinner(plotly::plotlyOutput("plot2_3d", height = "100%", width = "100%"))
                )
              )
            ),
            tabPanel(
              title = "Violin plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  shinycssloaders::withSpinner(plotOutput("plot2_vln", height = "100%", width = "100%"))
                )
              )
            )
          )
        )
      )
    )
  )

  # server ------------------------------------------------------------------
  server <- function(input, output, session) {
    # initial  ----------------------------------------------------------------
    updateSelectizeInput(session, "features2", choices = c(meta_features_name, all_features), selected = initial_feature, server = TRUE)

    output$plot1 <- renderPlot(
      {
        initial_p1_dim
      },
      width = attr(initial_p1_dim, "size")$width * initial_dpi,
      height = attr(initial_p1_dim, "size")$height * initial_dpi,
      res = initial_dpi
    )

    output$plot1_3d <- plotly::renderPlotly(
      {
        initial_p1_3d
      }
    )

    output$plot2 <- renderPlot(
      {
        initial_p2_dim
      },
      width = attr(initial_p2_dim, "size")$width * initial_dpi,
      height = attr(initial_p2_dim, "size")$height * initial_dpi,
      res = initial_dpi
    )

    output$plot2_3d <- plotly::renderPlotly(
      {
        initial_p2_3d
      }
    )

    output$plot2_vln <- renderPlot(
      {
        initial_p2_vln
      },
      width = attr(initial_p2_vln, "size")$width * initial_dpi,
      height = attr(initial_p2_vln, "size")$height * initial_dpi,
      res = initial_dpi
    )

    # change dataset  ----------------------------------------------------------------
    observeEvent(input$dataset1, {
      meta_groups_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/metadata.stat/asgroups"))
      reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset1, "/reductions"), "name"]
      default_reduction <- as.character(rhdf5::h5read(MetaFile, name = paste0("/", input$dataset1, "/reductions.stat/Default_reduction")))
      updateSelectInput(session, "reduction1", choices = reduction_name, selected = default_reduction)
      updateSelectInput(session, "class1", choices = meta_groups_name, selected = "orig.ident")
      updateSelectInput(session, "split1", choices = c("None", meta_groups_name), selected = "None")
    })

    observeEvent(input$dataset2, {
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
    })

    # submit1  ----------------------------------------------------------------
    observeEvent(input$submit1, ignoreInit = FALSE, {
      if (input$split1 == "None") {
        split1 <- NULL
      } else {
        split1 <- input$split1
      }

      # message("DataFile:", DataFile)
      # message("MetaFile:", MetaFile)
      # message("input$dataset1:", input$dataset1)
      # message("initial_slot:", initial_slot)
      # message("c(input$class1, split1):", c(input$class1, split1))
      # message("input$reduction1:", input$reduction1)

      srt_tmp <- SCP::FetchH5(
        DataFile = DataFile, MetaFile = MetaFile, name = input$dataset1,
        features = NULL, slot = initial_slot, assay = NULL,
        metanames = unique(c(input$class1, split1)), reduction = input$reduction1
      )

      raster1 <- input$raster1 == "TRUE"
      p1_dim <- SCP::CellDimPlot(
        srt = srt_tmp, group.by = input$class1, split.by = split1, reduction = input$reduction1, raster = raster1,
        label = ifelse(input$label1 == "Yes", TRUE, FALSE), palette = input$palette1, theme_use = input$theme1,
        ncol = input$ncol1, byrow = ifelse(input$arrange1 == "Row", TRUE, FALSE), force = TRUE
      )
      p1_dim <- SCP::panel_fix(p1_dim, height = input$size1, raster = FALSE, verbose = FALSE)

      output$plot1 <- renderPlot(
        {
          p1_dim
        },
        width = attr(p1_dim, "size")$width * input$plot_dpi1,
        height = attr(p1_dim, "size")$height * input$plot_dpi1,
        res = input$plot_dpi1
      )

      plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
      if (isTRUE(plot3d)) {
        p1_3d <- SCP::CellDimPlot3D(
          srt = srt_tmp, group.by = input$class1, reduction = input$reduction1, palette = input$palette1,
          force = TRUE
        )
      } else {
        p1_3d <- NULL
      }

      output$plot1_3d <- plotly::renderPlotly(
        {
          p1_3d
        }
      )
    })

    # submit2  ----------------------------------------------------------------
    observeEvent(input$submit2, ignoreInit = FALSE, {
      data <- HDF5Array::TENxMatrix(filepath = DataFile, group = paste0("/", input$dataset2, "/", initial_assay, "/", initial_slot))
      all_features <- colnames(data)
      meta_features_name <- rhdf5::h5read(MetaFile, name = paste0("/", input$dataset2, "/metadata.stat/asfeatures"))

      if (input$split2 == "None") {
        split2 <- NULL
      } else {
        split2 <- input$split2
      }
      input_features <- input$features2
      if (is.null(input_features)) {
        message("input gene is null")
        input_features <- meta_features_name[1]
      }
      gene_area <- gsub(x = unlist(strsplit(input$gene_area, "(\\r)|(\\n)", perl = TRUE)), pattern = " ", replacement = "")
      input_features <- c(as.character(input_features), as.character(gene_area))
      input_features <- unique(input_features[input_features %in% c(all_features, meta_features_name)])

      # message("DataFile:", DataFile)
      # message("MetaFile:", MetaFile)
      # message("input$dataset2:", input$dataset2)
      # message("input_features:", paste0(input_features,collapse = ","))
      # message("initial_slot:", initial_slot)
      # message("initial_assay:", initial_assay)
      # message("c(input$group2, split2):", c(input$group2, split2))
      # message("input$reduction2:", input$reduction2)

      srt_tmp <- SCP::FetchH5(
        DataFile = DataFile, MetaFile = MetaFile, name = input$dataset2,
        features = input_features, slot = input$slots2, assay = input$assays2,
        metanames = unique(c(input$group2, split2)), reduction = input$reduction2
      )

      raster2 <- input$raster2 == "TRUE"
      p2_dim <- SCP::FeatureDimPlot(
        srt = srt_tmp, features = input_features, split.by = split2, reduction = input$reduction2, slot = "data", raster = raster2,
        calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE), palette = input$palette2, theme_use = input$theme2,
        ncol = input$ncol2, byrow = ifelse(input$arrange2 == "Row", TRUE, FALSE)
      )
      p2_dim <- SCP::panel_fix(p2_dim, height = input$size2, raster = FALSE, verbose = FALSE)

      output$plot2 <- renderPlot(
        {
          p2_dim
        },
        width = attr(p2_dim, "size")$width * input$plot_dpi2,
        height = attr(p2_dim, "size")$height * input$plot_dpi2,
        res = input$plot_dpi2
      )

      p2_vln <- SCP::FeatureStatPlot(
        srt = srt_tmp, stat.by = input_features, group.by = input$group2, split.by = split2,
        calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE), palette = input$palette2,
        ncol = input$ncol2, byrow = ifelse(input$arrange2 == "Row", TRUE, FALSE), force = TRUE
      )
      p2_vln <- SCP::panel_fix(p2_vln, height = input$size2, raster = FALSE, verbose = FALSE)

      output$plot2_vln <- renderPlot(
        {
          p2_vln
        },
        width = attr(p2_vln, "size")$width * input$plot_dpi2,
        height = attr(p2_vln, "size")$height * input$plot_dpi2,
        res = input$plot_dpi2
      )

      plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
      if (isTRUE(plot3d)) {
        p2_3d <- SCP::FeatureDimPlot3D(
          srt = srt_tmp, features = input_features, reduction = input$reduction2,
          calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE)
        )
      } else {
        p2_3d <- NULL
      }

      output$plot2_3d <- plotly::renderPlotly(
        {
          p2_3d
        }
      )
    })
  }
  '

  main_code <- readLines(textConnection(main_code))
  main_code <- gsub("\\\\r", "\\\\\\\r", main_code)
  main_code <- gsub("\\\\n", "\\\\\\\n", main_code)
  args <- mget(names(formals()))
  args <- args[!names(args) %in% c("base_dir", "return_app", "create_script", "style_script", "overwrite")]
  for (varnm in names(args)) {
    main_code <- c(paste0(varnm, "=", deparse(args[[varnm]])), main_code)
  }
  main_code <- c(
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
    "SCP::check_R(c('shiny', 'shinycssloaders'))",
    "library(shiny)",
    main_code
  )
  main_code <- c(main_code, "shinyApp(ui = ui, server = server)")
  temp <- tempfile("SCExplorer")
  writeLines(main_code, temp)
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
