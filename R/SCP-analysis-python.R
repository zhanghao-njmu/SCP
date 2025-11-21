srt_to_adata <- function(srt, features = NULL,
                         assay_X = "RNA", slot_X = "counts",
                         assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                         convert_tools = FALSE, convert_misc = FALSE, verbose = TRUE) {
  check_Python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }

  if (is.null(features)) {
    features <- rownames(srt[[assay_X]])
  }
  if (length(slot_layers) == 1) {
    slot_layers <- rep(slot_layers, length(assay_layers))
    names(slot_layers) <- assay_layers
  } else if (length(slot_layers) != length(assay_layers)) {
    stop("slot_layers must be one character or the same length of the assay_layers")
  }

  sc <- import("scanpy", convert = FALSE)
  np <- import("numpy", convert = FALSE)

  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }

  var <- srt[[assay_X]]@meta.features[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(as.character(var[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }
  if (length(VariableFeatures(srt, assay = assay_X) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in% VariableFeatures(srt, assay = assay_X)
  }

  # X <- t(as_matrix(slot(srt@assays[[assay_X]], slot_X)[features, , drop = FALSE]))
  X <- t(GetAssayData(srt, assay = assay_X, slot = slot_X)[features, , drop = FALSE])
  adata <- sc$AnnData(
    X = np_array(X, dtype = np$float32),
    obs = obs,
    var = cbind(data.frame(features = features), var)
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_X]) {
    if (assay %in% assay_layers) {
      layer <- t(GetAssayData(srt, assay = assay, slot = slot_layers[assay]))
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          stop(
            "The following features in the '", assay_X, "' assay can not be found in the '", assay, "' assay:\n  ",
            paste0(head(colnames(X)[!colnames(X) %in% colnames(layer)], 10), collapse = ","), "..."
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      if (isTRUE(verbose)) {
        message("Assay '", assay, "' is in the srt object but not converted.")
      }
    }
  }
  if (length(layer_list) > 0) {
    adata$layers <- layer_list
  }

  reduction_list <- list()
  for (reduction in names(srt@reductions)) {
    reduction_list[[paste0(reduction)]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  # varm_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   if (ncol(srt[[reduction]]@feature.loadings) > 0) {
  #     varm_list[[paste0(reduction, "_feature.loadings")]] <- srt[[reduction]]@feature.loadings
  #   }
  # }
  # if (length(varm_list) > 0) {
  #   adata$varm <- varm_list
  # }

  obsp_list <- list()
  for (graph in names(srt@graphs)) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in names(srt@neighbors)) {
    obsp_list[[neighbor]] <- srt[[neighbor]]
  }
  if (length(obsp_list) > 0) {
    adata$obsp <- obsp_list
  }

  uns_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   uns_list[[paste0(reduction,"_stdev")]] <- srt[[reduction]]@stdev
  #   uns_list[[paste0(reduction,"_misc")]] <- srt[[reduction]]@misc
  # }
  # for (image in names(srt@images)) {
  #   uns_list[[image]] <- srt[[image]]
  # }
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@misc[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'misc' slot is not converted.")
    }
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'tools' slot is not converted.")
    }
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}

#' Convert an anndata object to a seurat object using reticulate
#'
#' @param adata a connected python anndata object.
#'
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(adata = adata, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' srt <- adata_to_srt(adata)
#' srt
#'
#' ### Or convert a h5ad file to Seurat object
#' # library(reticulate)
#' # check_Python("scanpy")
#' # sc <- import("scanpy")
#' # adata <- sc$read_h5ad("pancreas.h5ad")
#' # srt <- adata_to_srt(adata)
#' # srt
#' @importFrom Seurat CreateSeuratObject CreateAssayObject CreateDimReducObject AddMetaData
#' @importFrom SeuratObject as.Graph as.sparse
#' @importFrom reticulate iterate
#' @importFrom Matrix t
#' @export
adata_to_srt <- function(adata) {
  if (!inherits(adata, "python.builtin.object")) {
    stop("'adata' is not a python.builtin.object.")
  }
  sc <- import("scanpy", convert = TRUE)
  np <- import("numpy", convert = TRUE)

  x <- t(py_to_r_auto(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- as.sparse(x[1:nrow(x), , drop = FALSE])
  }
  rownames(x) <- py_to_r_auto(adata$var_names$values)
  colnames(x) <- py_to_r_auto(adata$obs_names$values)

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(py_to_r_auto(adata$obs))
    colnames(metadata) <- make.names(colnames(metadata))
  }

  srt <- CreateSeuratObject(counts = x, meta.data = metadata)

  if (inherits(adata$layers, "python.builtin.object")) {
    keys <- iterate(adata$layers$keys())
  } else {
    keys <- names(adata$layers)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      layer <- py_to_r_auto(adata$layers[[k]])
      if (!inherits(layer, c("Matrix", "matrix"))) {
        stop(paste0("The object in '", k, "' layers is not a matrix: ", paste0(class(adata$layers[[k]]), collapse = ",")))
      }
      layer <- t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- as.sparse(layer[1:nrow(layer), , drop = FALSE])
      }
      rownames(layer) <- py_to_r_auto(adata$var_names$values)
      colnames(layer) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateAssayObject(counts = layer)
    }
  }

  if (inherits(adata$obsm, "python.builtin.object")) {
    keys <- iterate(adata$obsm$keys())
  } else {
    keys <- names(adata$obsm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsm <- tryCatch(py_to_r_auto(adata$obsm[[k]]), error = identity)
      if (inherits(obsm, "error")) {
        warning("'obsm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- as_matrix(obsm)
      }
      k <- gsub(pattern = "^X_", replacement = "", x = py_to_r_auto(k))
      colnames(obsm) <- paste0(k, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateDimReducObject(embeddings = obsm, assay = "RNA", key = paste0(gsub(pattern = "_", replacement = "", x = k), "_"))
    }
  }

  if (inherits(adata$obsp, "python.builtin.object")) {
    keys <- iterate(adata$obsp$keys())
  } else {
    keys <- names(adata$obsp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsp <- tryCatch(py_to_r_auto(adata$obsp[[k]]), error = identity)
      if (inherits(obsp, "error")) {
        warning("'obsp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- as.sparse(obsp[1:nrow(obsp), , drop = FALSE])
      }
      colnames(obsp) <- py_to_r_auto(adata$obs_names$values)
      rownames(obsp) <- py_to_r_auto(adata$obs_names$values)
      obsp <- as.Graph(obsp[seq_len(nrow(obsp)), , drop = FALSE])
      DefaultAssay(object = obsp) <- "RNA"
      srt[[py_to_r_auto(k)]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- AddMetaData(srt[["RNA"]], metadata = as.data.frame(py_to_r_auto(adata$var)))
  }

  if (inherits(adata$varm, "python.builtin.object")) {
    keys <- iterate(adata$varm$keys())
  } else {
    keys <- names(adata$varm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varm <- tryCatch(py_to_r_auto(adata$varm[[k]]), error = identity)
      if (inherits(varm, "error")) {
        warning("'varm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- as_matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r_auto(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.loadings"]][[py_to_r_auto(k)]] <- varm
    }
  }

  if (inherits(adata$varp, "python.builtin.object")) {
    keys <- iterate(adata$varp$keys())
  } else {
    keys <- names(adata$varp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varp <- tryCatch(py_to_r_auto(adata$varp[[k]]), error = identity)
      if (inherits(varp, "error")) {
        warning("'varp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- as_matrix(varp)
      }
      colnames(varp) <- py_to_r_auto(adata$var_names$values)
      rownames(varp) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.graphs"]][[py_to_r_auto(k)]] <- varp
    }
  }

  if (inherits(adata$uns, "python.builtin.object")) {
    keys <- iterate(adata$uns$keys())
  } else {
    keys <- names(adata$uns)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      uns <- tryCatch(py_to_r_auto(adata$uns[[k]]), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r_auto(k)]] <- uns
      } else {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
    }
  }
  return(srt)
}

maxDepth <- function(x, depth = 0) {
  if (is.list(x)) {
    return(max(unlist(lapply(x, maxDepth, depth + 1))))
  } else {
    return(depth)
  }
}

#' @importFrom reticulate py_to_r
check_python_element <- function(x, depth = maxDepth(x)) {
  if (depth == 0 || !is.list(x) || !inherits(x, "python.builtin.object")) {
    if (inherits(x, "python.builtin.object")) {
      x_r <- tryCatch(py_to_r(x), error = identity)
      if (inherits(x_r, "error")) {
        return(x)
      } else {
        return(x_r)
      }
    } else {
      return(x)
    }
  } else {
    raw_depth <- maxDepth(x)
    x <- lapply(x, function(element) {
      if (inherits(element, "python.builtin.object")) {
        element_r <- tryCatch(py_to_r(element), error = identity)
        if (inherits(element_r, "error")) {
          return(element)
        } else {
          return(element_r)
        }
      } else {
        return(element)
      }
    })
    cur_depth <- maxDepth(x)
    if (cur_depth > raw_depth) {
      depth <- depth + 1
    }
    x_checked <- lapply(x, check_python_element, depth - 1)
    return(x_checked)
  }
}

#' Run PAGA analysis
#'
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay_X Assay to convert as the main data matrix (X) in the anndata object.
#' @param slot_X Slot name for assay_X in the Seurat object.
#' @param assay_layers Assays to convert as layers in the anndata object.
#' @param slot_layers Slot names for the assay_layers in the Seurat object.
#' @param adata An anndata object.
#' @param group_by Variable to use for grouping cells in the Seurat object.
#' @param linear_reduction Linear reduction method to use, e.g., "PCA".
#' @param nonlinear_reduction Non-linear reduction method to use, e.g., "UMAP".
#' @param basis The basis to use for reduction, e.g., "UMAP".
#' @param n_pcs Number of principal components to use for linear reduction. Default is 30.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph. Default is 30.
#' @param use_rna_velocity Whether to use RNA velocity for PAGA analysis. Default is FALSE.
#' @param vkey The name of the RNA velocity data to use if \code{use_rna_velocity} is TRUE. Default is "stochastic".
#' @param embedded_with_PAGA Whether to embed data using PAGA layout. Default is FALSE.
#' @param paga_layout The layout for plotting PAGA graph. See See \href{https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.paga.html}{layout} param in scanpy.pl.paga function.
#' @param threshold The threshold for plotting PAGA graph. Edges for weights below this threshold will not draw.
#' @param point_size The point size for plotting.
#' @param infer_pseudotime Whether to infer pseudotime.
#' @param root_group The group to use as the root for pseudotime inference.
#' @param root_cell The cell to use as the root for pseudotime inference.
#' @param n_dcs TThe number of diffusion components to use for pseudotime inference.
#' @param n_branchings Number of branchings to detect.
#' @param min_group_size The minimum size of a group (as a fraction of the total number of cells) to consider it as a potential branching point.
#' @param palette The palette to use for coloring cells.
#' @param palcolor A vector of colors to use as the palette.
#' @param show_plot Whether to show the PAGA plot.
#' @param dpi The DPI (dots per inch) for saving the PAGA plot.
#' @param save Whether to save the PAGA plots.
#' @param dirpath The directory to save the PAGA plots.
#' @param fileprefix The file prefix to use for the PAGA plots.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object. Default is TRUE.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{PAGAPlot}} \code{\link{CellDimPlot}} \code{\link{RunSCVELO}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, assay_X = "RNA", group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "draw_graph_fr")
#' PAGAPlot(pancreas_sub, reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga)
#'
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP",
#'   embedded_with_PAGA = TRUE, infer_pseudotime = TRUE, root_group = "Ductal"
#' )
#' head(pancreas_sub[[]])
#' names(pancreas_sub@reductions)
#' FeatureDimPlot(pancreas_sub, features = "dpt_pseudotime", reduction = "PAGAUMAP2D")
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "PAGAUMAP2D", paga = pancreas_sub@misc$paga)
#'
#' @export
RunPAGA <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                    adata = NULL, group_by = NULL,
                    linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                    n_pcs = 30, n_neighbors = 30, use_rna_velocity = FALSE, vkey = "stochastic",
                    embedded_with_PAGA = FALSE, paga_layout = "fr", threshold = 0.1, point_size = 20,
                    infer_pseudotime = FALSE, root_group = NULL, root_cell = NULL, n_dcs = 10, n_branchings = 0, min_group_size = 0.01,
                    palette = "Paired", palcolor = NULL,
                    show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                    return_seurat = !is.null(srt)) {
  check_Python("scanpy")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$PAGA, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(paga)|(distances)|(connectivities)|(draw_graph)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run scVelo workflow
#'
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells. This function runs scVelo workflow on a Seurat object.
#'
#' @inheritParams RunPAGA
#' @param mode Velocity estimation model to use, e.g., "stochastic".
#' @param fitting_by Method used to fit gene velocities, e.g., "stochastic".
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' @param knn The number of nearest neighbors for magic.MAGIC.
#' @param t power to which the diffusion operator is powered for magic.MAGIC.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' @param n_pcs Number of principal components (PCs) used for velocity estimation.
#' @param n_neighbors Number of nearest neighbors used for velocity estimation.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines. When density = 2 (default), the domain is divided into a 60x60 grid, whereas density linearly scales this grid. Each cell in the grid can have, at most, one traversing streamline.
#' @param arrow_length Length of arrows.
#' @param arrow_size Size of arrows.
#' @param arrow_density Amount of velocities to show.
#' @param denoise Boolean flag indicating whether to denoise.
#' @param denoise_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param kinetics Boolean flag indicating whether to estimate RNA kinetics.
#' @param kinetics_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param calculate_velocity_genes Boolean flag indicating whether to calculate velocity genes.
#' @param top_n The number of top features to plot.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{VelocityPlot}} \code{\link{CellDimPlot}} \code{\link{RunPAGA}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, assay_X = "RNA", group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' head(pancreas_sub[[]])
#' names(pancreas_sub@assays)
#'
#' FeatureDimPlot(pancreas_sub, c("stochastic_length", "stochastic_confidence"))
#' FeatureDimPlot(pancreas_sub, "stochastic_pseudotime")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = NA, velocity = "stochastic")
#'
#' pancreas_sub <- Standard_SCP(pancreas_sub, normalization_method = "SCT", nonlinear_reduction = "tsne")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, assay_X = "SCT", group_by = "SubCellType", linear_reduction = "Standardpca", nonlinear_reduction = "StandardTSNE2D")
#'
#' @export
RunSCVELO <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                      adata = NULL, group_by = NULL,
                      linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                      mode = "stochastic", fitting_by = "stochastic",
                      magic_impute = FALSE, knn = 5, t = 2,
                      min_shared_counts = 30, n_pcs = 30, n_neighbors = 30,
                      stream_smooth = NULL, stream_density = 2,
                      arrow_length = 5, arrow_size = 5, arrow_density = 0.5,
                      denoise = FALSE, denoise_topn = 3, kinetics = FALSE, kinetics_topn = 100,
                      calculate_velocity_genes = FALSE, top_n = 6, n_jobs = 1,
                      palette = "Paired", palcolor = NULL,
                      show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                      return_seurat = !is.null(srt)) {
  check_Python("scvelo")
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$SCVELO, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = paste0("(velocity)|(distances)|(connectivities)|(Ms)|(Mu)|(", paste(mode, collapse = ")|("), ")|(paga)"), overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run Palantir analysis
#'
#' @inheritParams RunPAGA
#' @param dm_n_components The number of diffusion components to calculate.
#' @param dm_alpha Normalization parameter for the diffusion operator.
#' @param dm_n_eigs Number of eigen vectors to use.
#' @param early_group Name of the group to start Palantir analysis from.
#' @param early_cell Name of the cell to start Palantir analysis from.
#' @param terminal_groups Character vector specifying terminal groups for Palantir analysis.
#' @param terminal_cells Character vector specifying terminal cells for Palantir analysis.
#' @param num_waypoints Number of waypoints to be included.
#' @param scale_components Should the cell fate probabilities be scaled for each component independently?
#' @param use_early_cell_as_start Should the starting cell for each terminal group be set as early_cell?
#' @param adjust_early_cell Whether to adjust the early cell to the cell with the minimum pseudotime value.
#' @param adjust_terminal_cells hether to adjust the terminal cells to the cells with the maximum pseudotime value for each terminal group.
#' @param max_iterations Maximum number of iterations for pseudotime convergence.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPalantir(
#'   srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP",
#'   early_group = "Ductal", use_early_cell_as_start = TRUE,
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#' head(pancreas_sub[[]])
#' FeatureDimPlot(pancreas_sub, c("palantir_pseudotime", "palantir_diff_potential"))
#' FeatureDimPlot(pancreas_sub, paste0(c("Alpha", "Beta", "Delta", "Epsilon"), "_diff_potential"))
#'
#' @export
RunPalantir <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, group_by = NULL,
                        linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                        n_pcs = 30, n_neighbors = 30, dm_n_components = 10, dm_alpha = 0, dm_n_eigs = NULL,
                        early_group = NULL, early_cell = NULL, terminal_cells = NULL, terminal_groups = NULL,
                        num_waypoints = 1200, scale_components = TRUE, use_early_cell_as_start = TRUE,
                        adjust_early_cell = FALSE, adjust_terminal_cells = FALSE,
                        max_iterations = 25, n_jobs = 8,
                        point_size = 20, palette = "Paired", palcolor = NULL,
                        show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                        return_seurat = !is.null(srt)) {
  check_Python("palantir")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by) && any(!is.null(early_group), !is.null(terminal_groups))) {
    stop("'group_by' must be provided when early_group or terminal_groups provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
  }
  if (is.null(early_cell) && is.null(early_group)) {
    stop("'early_cell' or 'early_group' must be provided.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$Palantir, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(palantir)|(dm_kernel)|(_diff_potential)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run WOT analysis
#'
#' @inheritParams RunPAGA
#' @param time_field A character string specifying the column name in `adata.obs` or `srt@meta.data` that contains the time information.
#' @param growth_iters An integer specifying the number of growth iterations to perform during the OT Model computation. Default is 3.
#' @param tmap_out A character string specifying the path to store the computed transport maps.
#' @param time_from A numeric value specifying the starting time point for trajectory and fate analysis.
#' @param time_to A numeric value specifying the ending time point for trajectory and fate analysis. If not provided, only trajectory and fate analysis for the specified `time_from` will be performed.
#' @param get_coupling A logical value indicating whether to compute and store the coupling matrix between the specified `time_from` and `time_to`. Default is FALSE.
#' @param recalculate A logical value indicating whether to recalculate the transport maps even if they already exist at the specified `tmap_out` location. Default is FALSE.
#'
#' @seealso \code{\link{srt_to_adata}}
#' @export
RunWOT <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                   adata = NULL, group_by = NULL,
                   time_field = "Time", growth_iters = 3L, tmap_out = "tmaps/tmap_out",
                   time_from = NULL, time_to = NULL, get_coupling = FALSE, recalculate = FALSE,
                   palette = "Paired", palcolor = NULL,
                   show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                   return_seurat = !is.null(srt)) {
  check_Python("wot")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(time_field)) {
    stop("'time_field' must be provided.")
  }
  if (is.null(time_from)) {
    stop("'time_from' must be provided.")
  }
  if (isTRUE(get_coupling) && is.null(time_to)) {
    warning("The 'get_coupling' paramter is only valid when 'time_to' is specified.")
  }

  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$WOT, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(trajectory_)|(fates_)|(transition_)|(coupling_)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

RunCellRank <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, group_by = NULL, n_jobs = 1,
                        linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                        mode = "stochastic", fitting_by = "stochastic",
                        magic_impute = FALSE, knn = 5, t = 2,
                        min_shared_counts = 30, n_pcs = 30, n_neighbors = 30, approx = TRUE,
                        stream_smooth = NULL, stream_density = 2,
                        arrow_size = 5, arrow_length = 5, arrow_density = 0.5,
                        s_genes = NULL, g2m_genes = NULL, calculate_velocity_genes = FALSE,
                        denoise = FALSE, kinetics = FALSE, axis = "equal",
                        show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                        return_seurat = !is.null(srt)) {
  check_Python("cellrank")
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
  }
  mode <- as.list(mode)
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$CellRank, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      return(SrtAppend(srt_raw = srt, srt_append = srt_out))
    }
  } else {
    return(adata)
  }
}
