RunSlingshot <- function(srt, group.by, reduction = NULL, dims = NULL, start = NULL, end = NULL, prefix = NULL,
                         reverse = FALSE, align_start = FALSE, show_plot = TRUE, lineage_palette = "Dark2", seed = 11, ...) {
  if (missing(group.by)) {
    stop("group.by is missing")
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (is.null(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }
  srt_sub <- srt[, !is.na(srt[[group.by, drop = TRUE]])]
  if (is.null(dims)) {
    dims <- 1:2
  }

  set.seed(seed)
  sl <- slingshot(
    data = as.data.frame(srt_sub[[reduction]]@cell.embeddings[, dims]),
    clusterLabels = as.character(srt_sub[[group.by, drop = TRUE]]),
    start.clus = start, end.clus = end, ...
  )

  srt@tools[[paste("Slingshot", group.by, reduction, sep = "_")]] <- sl
  df <- as.data.frame(slingPseudotime(sl))
  colnames(df) <- paste0(prefix, colnames(df))
  if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
      df <- max(df, na.rm = TRUE) - df
    }
  }
  srt <- AddMetaData(srt, metadata = df)
  srt <- AddMetaData(srt, metadata = slingBranchID(sl), col.name = paste0(prefix, "BranchID"))

  if (isTRUE(show_plot)) {
    if (ncol(srt[[reduction]]@cell.embeddings) == 2 || ncol(srt[[reduction]]@cell.embeddings) > 3) {
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 2, type = "lineages", col = "black")
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 3, col = 1:length(slingshot::SlingshotDataSet(sl)@lineages))
      p <- CellDimPlot(srt, group.by = group.by, reduction = reduction, dims = c(1, 2), lineages = colnames(df))
      print(p)
    } else if (ncol(srt[[reduction]]@cell.embeddings) == 3) {
      p <- CellDimPlot3D(srt, group.by = group.by, reduction = reduction, lineages = colnames(df))
      print(p)
    }
  }
  return(srt)
}

#' Run Monocle2 analysis
#'
#' Runs the Monocle2 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param slot The slot in the Seurat object to use for analysis. Default is "counts".
#' @param expressionFamily The distribution family to use for modeling gene expression. Default is "negbinomial.size".
#' @param features A vector of gene names or indices specifying the features to use in the analysis. Defaults to NULL, in which case features were determined by \code{feature_type}.
#' @param feature_type The type of features to use in the analysis. Possible values are "HVF" for highly variable features
#'   or "Disp" for features selected based on dispersion. Default is "HVF".
#' @param disp_filter A string specifying the filter to use when \code{feature_type} is "Disp". Default is
#'   "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit".
#' @param max_components The maximum number of dimensions to use for dimensionality reduction. Default is 2.
#' @param reduction_method The dimensionality reduction method to use. Possible values are "DDRTree" and "UMAP". Default is "DDRTree".
#' @param norm_method The normalization method to use. Possible values are "log" and "none". Default is "log".
#' @param residualModelFormulaStr A model formula specifying the effects to subtract. Default is NULL.
#' @param pseudo_expr Amount to increase expression values before dimensionality reduction. Default is 1.
#' @param root_state The state to use as the root of the trajectory. If NULL, will prompt for user input.
#' @param seed An integer specifying the random seed to use. Default is 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   pancreas_sub <- RunMonocle2(srt = pancreas_sub)
#'   names(pancreas_sub@tools$Monocle2)
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#'
#'   pancreas_sub <- RunMonocle2(
#'     srt = pancreas_sub,
#'     feature_type = "Disp", disp_filter = "mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#' }
#' @importFrom Seurat DefaultAssay CreateDimReducObject GetAssayData VariableFeatures FindVariableFeatures
#' @importFrom SeuratObject as.sparse
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @importFrom utils select.list
#' @export
RunMonocle2 <- function(srt, assay = NULL, slot = "counts", expressionFamily = "negbinomial.size",
                        features = NULL, feature_type = "HVF", disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit",
                        max_components = 2, reduction_method = "DDRTree", norm_method = "log", residualModelFormulaStr = NULL, pseudo_expr = 1,
                        root_state = NULL, seed = 11) {
  set.seed(seed)
  check_R(c("monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM"))

  if (!"package:DDRTree" %in% search()) {
    attachNamespace("DDRTree")
  }

  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- as.sparse(GetAssayData(srt, assay = assay, slot = slot))
  p_data <- srt@meta.data
  f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
  pd <- new("AnnotatedDataFrame", data = p_data)
  fd <- new("AnnotatedDataFrame", data = f_data)
  cds <- monocle::newCellDataSet(expr_matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily = do.call(get(expressionFamily, envir = getNamespace("VGAM")), args = list())
  )
  if (any(c("negbinomial", "negbinomial.size") %in% expressionFamily)) {
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- BiocGenerics::estimateDispersions(cds)
  }
  if (is.null(features)) {
    if (feature_type == "HVF") {
      features <- VariableFeatures(srt, assay = assay)
      if (length(features) == 0) {
        features <- VariableFeatures(FindVariableFeatures(srt, assay = assay), assay = assay)
      }
    }
    if (feature_type == "Disp") {
      features <- subset(monocle::dispersionTable(cds), eval(rlang::parse_expr(disp_filter)))$gene_id
    }
  }
  message("features number: ", length(features))
  cds <- monocle::setOrderingFilter(cds, features)
  p <- monocle::plot_ordering_genes(cds)
  print(p)

  cds <- monocle::reduceDimension(
    cds = cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr
  )
  cds <- orderCells(cds)

  embeddings <- t(cds@reducedDimS)
  colnames(embeddings) <- paste0(cds@dim_reduce_type, "_", 1:ncol(embeddings))
  srt[[cds@dim_reduce_type]] <- CreateDimReducObject(embeddings = embeddings, key = paste0(cds@dim_reduce_type, "_"), assay = assay)
  srt[["Monocle2_State"]] <- cds[["State"]]
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimS))
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimK))
  }
  edge_df <- as_data_frame(cds@minSpanningTree)
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  trajectory <- geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
  p <- CellDimPlot(srt, group.by = "Monocle2_State", reduction = reduction_method, label = TRUE, force = TRUE) +
    trajectory
  print(p)
  if (is.null(root_state)) {
    root_state <- select.list(sort(unique(cds[["State"]])), title = "Select the root state to order cells:")
    if (root_state == "" || length(root_state) == 0) {
      root_state <- NULL
    }
  }
  cds <- orderCells(cds, root_state = root_state)
  srt[["Monocle2_State"]] <- cds[["State"]]
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(cds = cds, features = features, trajectory = trajectory)

  p1 <- CellDimPlot(srt, group.by = "Monocle2_State", reduction = reduction_method, label = TRUE, force = TRUE) + trajectory
  p2 <- FeatureDimPlot(srt, features = "Monocle2_Pseudotime", reduction = reduction_method) + trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}

orderCells <- function(cds, root_state = NULL, num_paths = NULL, reverse = NULL) {
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths <- 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- as_matrix(stats::dist(adjusted_S))
    cellPairwiseDistances(cds) <- as_matrix(stats::dist(adjusted_S))
    gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::minimum.spanning.tree(gp)
    monocle::minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- monocle:::pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, monocle::cellPairwiseDistances(cds), num_paths, FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    monocle::minSpanningTree(cds) <- igraph::as.undirected(order_list$cell_ordering_tree)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
    monocle::minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- igraph::as.undirected(order_list$cell_ordering_tree)
  } else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    K_old <- monocle::reducedDimK(cds)
    old_dp <- monocle::cellPairwiseDistances(cds)
    old_mst <- monocle::minSpanningTree(cds)
    old_A <- monocle::reducedDimA(cds)
    old_W <- monocle::reducedDimW(cds)
    cds <- project2MST(cds, monocle:::project_point_to_line_segment)
    monocle::minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(igraph::V(old_mst)$name == root_cell, arr.ind = TRUE)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- igraph::V(monocle::minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(igraph::degree(monocle::minSpanningTree(cds)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering_new_pseudotime[row.names(Biobase::pData(cds)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      cds[["State"]] <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
    cds@reducedDimK <- K_old
    cds@cellPairwiseDistances <- old_dp
    cds@minSpanningTree <- old_mst
    cds@reducedDimA <- old_A
    cds@reducedDimW <- old_W
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}
project2MST <- function(cds, Projection_Method) {
  dp_mst <- monocle::minSpanningTree(cds)
  Z <- monocle::reducedDimS(cds)
  Y <- monocle::reducedDimK(cds)
  cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as_matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(igraph::V(dp_mst)[suppressWarnings(nei(closest_vertex_names[i], mode = "all"))])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- monocle:::projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      if (!inherits(projection, "matrix")) {
        projection <- as_matrix(projection)
      }
      P[, i] <- projection[which(distance == min(distance))[1], ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- as_matrix(stats::dist(t(P)))
  min_dist <- min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  monocle::cellPairwiseDistances(cds) <- dp
  gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- igraph::minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}
extract_ddrtree_ordering <- function(cds, root_cell, verbose = TRUE) {
  dp <- monocle::cellPairwiseDistances(cds)
  dp_mst <- monocle::minSpanningTree(cds)
  curr_state <- 1
  res <- list(subtree = dp_mst, root = root_cell)
  states <- rep(1, ncol(dp))
  names(states) <- igraph::V(dp_mst)$name
  pseudotimes <- rep(0, ncol(dp))
  names(pseudotimes) <- igraph::V(dp_mst)$name
  parents <- rep(NA, ncol(dp))
  names(parents) <- igraph::V(dp_mst)$name
  mst_traversal <- igraph::graph.dfs(dp_mst,
    root = root_cell, mode = "all",
    unreachable = FALSE, father = TRUE
  )
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  for (i in 1:length(mst_traversal$order)) {
    curr_node <- mst_traversal$order[i]
    curr_node_name <- igraph::V(dp_mst)[curr_node]$name
    if (is.na(mst_traversal$father[curr_node]) == FALSE) {
      parent_node <- mst_traversal$father[curr_node]
      parent_node_name <- igraph::V(dp_mst)[parent_node]$name
      parent_node_pseudotime <- pseudotimes[parent_node_name]
      parent_node_state <- states[parent_node_name]
      curr_node_pseudotime <- parent_node_pseudotime +
        dp[curr_node_name, parent_node_name]
      if (igraph::degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
    } else {
      parent_node <- NA
      parent_node_name <- NA
      curr_node_pseudotime <- 0
    }
    curr_node_state <- curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  ordering_df <- data.frame(
    sample_name = names(states), cell_state = factor(states),
    pseudo_time = as.vector(pseudotimes), parent = parents
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

#' Run Monocle3 analysis
#'
#' Runs the Monocle3 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param slot The slot in the Seurat object to use for analysis. Default is "counts".
#' @param reduction The reduction used. Defaults to NULL, in which case the default reduction of the Seurat object is used.
#' @param clusters The cluster variable in the Seurat object to use for analysis. Defaults to NULL, in which case use Monocle clusters is used.
#' @param graph The name of the graph slot in the Seurat object to use for analysis. Defaults to NULL, in which case Monocle graph is used.
#' @param partition_qval The q-value threshold for partitioning cells. Defaults to 0.05.
#' @param k The number of nearest neighbors to consider for clustering. Defaults to 50.
#' @param cluster_method The clustering method to use. Defaults to "louvain".
#' @param num_iter The number of iterations for clustering. Defaults to 2.
#' @param resolution The resolution parameter for clustering. Defaults to NULL.
#' @param use_partition Whether to use partitions to learn disjoint graph in each partition. If not specified, user will be prompted for input. Defaults to NULL.
#' @param close_loop Whether to close loops in the graph. Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param root_cells The root cells to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param seed The random seed to use for reproducibility. Defaults to 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   # Use Monocle clusters to infer the trajectories
#'   pancreas_sub <- RunMonocle3(srt = pancreas_sub, reduction = "UMAP")
#'   names(pancreas_sub@tools$Monocle3)
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   milestones <- pancreas_sub@tools$Monocle3$milestones
#'
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory + milestones
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "UMAP", theme_use = "theme_blank") + trajectory
#'
#'   ## Select the lineage using monocle3::choose_graph_segments
#'   # cds <- pancreas_sub@tools$Monocle3$cds
#'   # cds_sub <- monocle3::choose_graph_segments(cds, starting_pr_node = NULL, ending_pr_nodes = NULL)
#'   # pancreas_sub$Lineages_1 <- NA
#'   # pancreas_sub$Lineages_1[colnames(cds_sub)] <- pancreas_sub$Monocle3_Pseudotime[colnames(cds_sub)]
#'   # CellDimPlot(pancreas_sub, group.by = "SubCellType", lineages = "Lineages_1", lineages_span = 0.1, theme_use = "theme_blank")
#'
#'   # Use Seurat clusters to infer the trajectories
#'   pancreas_sub <- Standard_SCP(pancreas_sub)
#'   CellDimPlot(pancreas_sub, group.by = c("Standardclusters", "CellType"), label = TRUE, theme_use = "theme_blank")
#'   pancreas_sub <- RunMonocle3(srt = pancreas_sub, clusters = "Standardclusters")
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "StandardUMAP2D", theme_use = "theme_blank") + trajectory
#'
#'   # Use custom graphs and cell clusters to infer the partitions and trajectories, respectively
#'   pancreas_sub <- Standard_SCP(pancreas_sub, cluster_resolution = 5)
#'   CellDimPlot(pancreas_sub, group.by = c("Standardclusters", "CellType"), label = TRUE)
#'   pancreas_sub <- RunMonocle3(
#'     srt = pancreas_sub,
#'     clusters = "Standardclusters", graph = "Standardpca_SNN"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "StandardUMAP2D", theme_use = "theme_blank") + trajectory
#' }
#' @importFrom SeuratObject as.sparse Embeddings Loadings Stdev
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color
#' @importFrom utils packageVersion select.list
#' @export
RunMonocle3 <- function(srt, assay = NULL, slot = "counts",
                        reduction = DefaultReduction(srt), clusters = NULL, graph = NULL, partition_qval = 0.05,
                        k = 50, cluster_method = "louvain", num_iter = 2, resolution = NULL,
                        use_partition = NULL, close_loop = TRUE,
                        root_pr_nodes = NULL, root_cells = NULL, seed = 11) {
  set.seed(seed)
  if (!requireNamespace("monocle3", quietly = TRUE) || packageVersion("monocle3") < package_version("1.2.0")) {
    check_R("cole-trapnell-lab/monocle3", force = TRUE)
  }
  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- as.sparse(GetAssayData(srt, assay = assay, slot = slot))
  p_data <- srt@meta.data
  f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
  cds <- monocle3::new_cell_data_set(
    expression_data = expr_matrix,
    cell_metadata = p_data,
    gene_metadata = f_data
  )
  if (!"Size_Factor" %in% colnames(cds@colData)) {
    size.factor <- paste0("nCount_", assay)
    if (size.factor %in% colnames(srt@meta.data)) {
      cds[["Size_Factor"]] <- cds[[size.factor, drop = TRUE]]
    }
  }
  reduction <- reduction %||% DefaultReduction(srt)
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- Embeddings(srt[[reduction]])
  loadings <- Loadings(object = srt[[reduction]])
  if (length(loadings) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["gene_loadings"]] <- loadings
  }
  stdev <- Stdev(object = srt[[reduction]])
  if (length(stdev) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["prop_var_expl"]] <- stdev
  }

  if (!is.null(clusters)) {
    if (!is.null(graph)) {
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = srt[[graph]],
        weighted = TRUE
      )
      cluster_result <- list(
        g = g,
        relations = NULL,
        distMatrix = "matrix",
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(as.factor(srt[[clusters, drop = TRUE]])),
          modularity = NA_real_
        )
      )
      if (length(unique(cluster_result$optim_res$membership)) > 1) {
        cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, cluster_result$optim_res, partition_qval)
        partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
        partitions <- as.factor(partitions)
      } else {
        partitions <- rep(1, ncol(srt))
      }
      names(partitions) <- colnames(cds)
      cds@clusters[["UMAP"]] <- list(
        cluster_result = cluster_result,
        partitions = partitions,
        clusters = as.factor(srt[[clusters, drop = TRUE]])
      )
      cds[["clusters"]] <- cds[[clusters]]
      cds <- monocle3:::add_citation(cds, "clusters")
      cds <- monocle3:::add_citation(cds, "partitions")
    } else {
      cds <- monocle3::cluster_cells(cds,
        reduction_method = "UMAP", partition_qval = partition_qval,
        k = k, cluster_method = cluster_method, num_iter = num_iter, resolution = resolution
      )
      cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters <- as.factor(srt[[clusters, drop = TRUE]])
    }
  } else {
    cds <- monocle3::cluster_cells(cds,
      reduction_method = "UMAP", partition_qval = partition_qval,
      k = k, cluster_method = cluster_method, num_iter = num_iter, resolution = resolution
    )
    cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
  }
  srt[["Monocle3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
  srt[["Monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions
  p1 <- CellDimPlot(srt, "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE)
  p2 <- CellDimPlot(srt, "Monocle3_clusters", reduction = reduction, label = FALSE, force = TRUE)
  print(p1)
  Sys.sleep(1)
  print(p2)
  if (is.null(use_partition)) {
    use_partition <- select.list(c(TRUE, FALSE), title = "Whether to use partitions to learn disjoint graph in each partition?")
    if (use_partition == "" || length(use_partition) == 0) {
      use_partition <- TRUE
    }
  }
  cds <- monocle3::learn_graph(cds = cds, use_partition = use_partition, close_loop = close_loop)

  reduced_dim_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  edge_df <- as_data_frame(cds@principal_graph[["UMAP"]])
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
  mst_leaf_nodes <- monocle3:::leaf_nodes(cds, "UMAP")
  mst_root_nodes <- monocle3:::root_nodes(cds, "UMAP")
  pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
  point_df <- data.frame(nodes = names(pps), x = reduced_dim_coords[pps, 1], y = reduced_dim_coords[pps, 2])
  point_df[, "is_branch"] <- names(pps) %in% names(mst_branch_nodes)
  trajectory <- list(
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
  )
  milestones <- list(
    geom_point(
      data = point_df[point_df[["is_branch"]] == FALSE, , drop = FALSE], aes(x = x, y = y),
      shape = 21, color = "white", fill = "black", size = 3, stroke = 1
    ),
    geom_point(
      data = point_df[point_df[["is_branch"]] == TRUE, , drop = FALSE], aes(x = x, y = y),
      shape = 21, color = "white", fill = "red", size = 3, stroke = 1
    ),
    new_scale_color(),
    geom_text_repel(
      data = point_df, aes(x = x, y = y, label = nodes, color = is_branch),
      fontface = "bold", min.segment.length = 0,
      point.size = 3, max.overlaps = 100,
      bg.color = "white", bg.r = 0.1, size = 3.5
    ),
    scale_color_manual(values = setNames(c("red", "black"), nm = c(TRUE, FALSE)))
  )
  p <- CellDimPlot(srt, group.by = "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE) +
    trajectory + milestones
  print(p)

  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    root_pr_nodes <- select.list(names(pps), title = "Select the root nodes to order cells, or leave blank for interactive selection:", multiple = TRUE)
    if (root_pr_nodes == "" || length(root_pr_nodes) == 0) {
      root_pr_nodes <- NULL
    }
  }
  cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes, root_cells = root_cells)
  pseudotime <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime[is.infinite(pseudotime)] <- NA
  srt[["Monocle3_Pseudotime"]] <- pseudotime
  srt@tools$Monocle3 <- list(cds = cds, trajectory = trajectory, milestones = milestones)

  p1 <- CellDimPlot(srt, group.by = "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE) +
    trajectory
  p2 <- FeatureDimPlot(srt, features = "Monocle3_Pseudotime", reduction = reduction) +
    theme(legend.position = "none") +
    trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}

#' RunDynamicFeatures
#'
#' Calculates dynamic features for lineages in a single-cell RNA-seq dataset.
#'
#' @param srt A Seurat object.
#' @param lineages A character vector specifying the lineage names for which dynamic features should be calculated.
#' @param features A character vector specifying the features (genes or metadata variables) for which dynamic features should be calculated. If NULL, n_candidates must be provided.
#' @param suffix A character vector specifying the suffix to append to the output slot names for each lineage. Defaults to the lineage names.
#' @param n_candidates An integer specifying the number of candidate features to select when features is NULL. Defaults to 1000.
#' @param minfreq An integer specifying the minimum frequency threshold for candidate features. Features with a frequency less than minfreq will be excluded. Defaults to 5.
#' @param family A character or character vector specifying the family of distributions to use for the generalized additive models (GAMs). If family is set to NULL, the appropriate family will be automatically determined based on the data. If length(family) is 1, the same family will be used for all features. Otherwise, family must have the same length as features.
#' @param slot A character vector specifying the slot in the Seurat object to use. Default is "counts".
#' @param assay A character vector specifying the assay in the Seurat object to use. Default is NULL.
#' @param libsize A numeric or numeric vector specifying the library size correction factors for each cell. If NULL, the library size correction factors will be calculated based on the expression matrix. If length(libsize) is 1, the same value will be used for all cells. Otherwise, libsize must have the same length as the number of cells in srt. Defaults to NULL.
#' @param BPPARAM A BiocParallelParam object specifying the parallelization parameters for computing dynamic features. Defaults to BiocParallel::bpparam().
#' @param seed An integer specifying the seed for random number generation. Defaults to 11.
#'
#' @return Returns the modified Seurat object with the calculated dynamic features stored in the tools slot.
#'
#' @seealso \code{\link{DynamicHeatmap}} \code{\link{DynamicPlot}} \code{\link{RunDynamicEnrichment}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#' names(pancreas_sub@tools$DynamicFeatures_Lineage1)
#' head(pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures)
#' ht <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 6, reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers
#' @export
RunDynamicFeatures <- function(srt, lineages, features = NULL, suffix = lineages,
                               n_candidates = 1000, minfreq = 5,
                               family = NULL,
                               slot = "counts", assay = NULL, libsize = NULL,
                               BPPARAM = BiocParallel::bpparam(), seed = 11) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  assay <- assay %||% DefaultAssay(srt)

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start RunDynamicFeatures"))
  message("Workers: ", bpworkers(BPPARAM))

  check_R("mgcv")
  meta <- c()
  gene <- c()
  if (!is.null(features)) {
    gene <- features[features %in% rownames(srt[[assay]])]
    meta <- features[features %in% colnames(srt@meta.data)]
    isnum <- sapply(srt@meta.data[, meta, drop = FALSE], is.numeric)
    if (!all(isnum)) {
      warning(paste0(meta[!isnum], collapse = ","), " is not numeric and will be dropped.", immediate. = TRUE)
      meta <- meta[isnum]
    }
    features <- c(gene, meta)
    if (length(features) == 0) {
      stop("No feature found in the srt object.")
    }
  }

  Y <- GetAssayData(srt, slot = slot, assay = assay)
  if (is.null(libsize)) {
    status <- check_DataType(srt, assay = assay, slot = "counts")
    if (status != "raw_counts") {
      Y_libsize <- setNames(rep(1, ncol(srt)), colnames(srt))
    } else {
      Y_libsize <- colSums(GetAssayData(srt, slot = "counts", assay = assay))
    }
  } else {
    if (length(libsize) == 1) {
      Y_libsize <- setNames(rep(libsize, ncol(srt)), colnames(srt))
    } else if (length(libsize) == ncol(srt)) {
      Y_libsize <- setNames(libsize, colnames(srt))
    } else {
      stop("libsize must be length of 1 or the number of cells.")
    }
  }

  if (length(meta) > 0) {
    Y <- rbind(Y, t(srt@meta.data[, meta, drop = FALSE]))
  }

  features_list <- c()
  srt_sub_list <- list()
  for (l in lineages) {
    srt_sub <- subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])])
    if (is.null(features)) {
      if (is.null(n_candidates)) {
        stop("'features' or 'n_candidates' must provided at least one.")
      }
      HVF <- VariableFeatures(FindVariableFeatures(srt_sub, nfeatures = n_candidates, assay = assay), assay = assay)
      HVF_counts <- srt_sub[[assay]]@counts[HVF, , drop = FALSE]
      HVF <- HVF[apply(HVF_counts, 1, function(x) {
        length(unique(x))
      }) >= minfreq]
      features_list[[l]] <- HVF
    } else {
      features_list[[l]] <- features
    }
    srt_sub_list[[l]] <- srt_sub
  }
  features <- unique(unlist(features_list))
  gene <- features[features %in% rownames(srt[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  message("Number of candidate features(union): ", length(features))

  if (slot == "counts") {
    gene_status <- status
  }
  gene_status <- status <- check_DataType(srt, assay = assay, slot = slot)
  meta_status <- sapply(meta, function(x) {
    check_DataType(data = srt[[x]])
  })
  if (is.null(family)) {
    family <- rep("gaussian", length(features))
    names(family) <- features
    family[names(meta_status)[meta_status == "raw_counts"]] <- "nb"
    if (gene_status == "raw_counts") {
      family[gene] <- "nb"
    }
  } else {
    if (length(family) == 1) {
      family <- rep(family, length(features))
      names(family) <- features
    }
    if (length(family) != length(features)) {
      stop("'family' must be one character or a vector of the same length as features.")
    }
  }

  for (i in seq_along(lineages)) {
    l <- lineages[i]
    srt_sub <- srt_sub_list[[l]]
    t <- srt_sub[[l, drop = TRUE]]
    t <- t[is.finite(t)]
    t_ordered <- t[order(t)]
    Y_ordered <- as_matrix(Y[features, names(t_ordered), drop = FALSE])
    l_libsize <- Y_libsize[names(t_ordered)]
    raw_matrix <- as_matrix(cbind(data.frame(pseudotime = t_ordered), t(Y_ordered)))

    # df <- data.frame(x = rowMeans(Y_ordered), y = MatrixGenerics::rowVars(Y_ordered))
    # p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
    #   geom_point() +
    #   geom_abline(slope = 1, intercept = 0, color = "red") +
    #   labs(title = l, subtitle = "Var VS Mean") +
    #   theme_scp()
    # print(p)

    message("Calculate dynamic features for ", l, "...")
    system.time({
      gam_out <- bplapply(seq_len(nrow(Y_ordered)), function(n, Y_ordered, t_ordered, l_libsize, family) {
        feature_nm <- rownames(Y_ordered)[n]
        family_current <- family[feature_nm]
        if (min(Y_ordered[feature_nm, ]) < 0 && family_current %in% c("nb", "poisson", "binomial")) {
          warning("Negative values detected. Replace family with 'gaussian' for the feature: ", feature_nm, immediate. = TRUE)
          family_use <- "gaussian"
        } else {
          family_use <- family_current
        }
        if (slot == "counts" && family_use != "gaussian" && !feature_nm %in% meta) {
          l_libsize <- l_libsize
        } else {
          l_libsize <- rep(median(Y_libsize), ncol(Y_ordered))
        }
        sizefactror <- median(Y_libsize) / l_libsize
        mod <- mgcv::gam(y ~ s(x, bs = "cs") + offset(log(l_libsize)),
          family = family_use,
          data = data.frame(y = Y_ordered[feature_nm, ], x = t_ordered, l_libsize = l_libsize)
        )
        pre <- predict(mod, type = "link", se.fit = TRUE)
        upr <- pre$fit + (2 * pre$se.fit)
        lwr <- pre$fit - (2 * pre$se.fit)
        upr <- mod$family$linkinv(upr)
        lwr <- mod$family$linkinv(lwr)
        res <- summary(mod)
        fitted <- fitted(mod)
        pvalue <- res$s.table[[4]]
        dev.expl <- res$dev.expl
        r.sq <- res$r.sq
        fitted.values <- fitted * sizefactror
        upr.values <- upr * sizefactror
        lwr.values <- lwr * sizefactror
        exp_ncells <- sum(Y_ordered[feature_nm, ] > min(Y_ordered[feature_nm, ]), na.rm = TRUE)
        peaktime <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
        valleytime <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])

        # ggplot(data = data.frame(
        #   x = t_ordered,
        #   raw = FetchData(srt_sub, vars = feature_nm, slot = "counts")[names(t_ordered), feature_nm, drop = TRUE],
        #   fitted = fitted.values,
        #   upr.values = upr.values,
        #   lwr.values = lwr.values,
        #   l_libsize = l_libsize
        # )) +
        #   geom_point(aes(x = x, y = raw), color = "black", size = 0.5) +
        #   geom_point(aes(x = x, y = fitted), color = "red", size = 0.5) +
        #   geom_path(aes(x = x, y = upr.values), color = "blue") +
        #   geom_path(aes(x = x, y = lwr.values), color = "green")

        # a <- data.frame(x = t_ordered, y = Y_ordered[feature_nm, ])
        # qplot(a$x, a$y)
        # length(unique(a$y) > 5)
        # a <- a[a$y > 0, ]
        # qplot(a$x, a$y)

        return(list(
          features = feature_nm, exp_ncells = exp_ncells,
          r.sq = r.sq, dev.expl = dev.expl,
          peaktime = peaktime, valleytime = valleytime,
          pvalue = pvalue, fitted.values = fitted.values,
          upr.values = upr.values, lwr.values = lwr.values
        ))
      }, BPPARAM = BPPARAM, Y_ordered = Y_ordered, t_ordered = t_ordered, l_libsize = l_libsize, family = family)
    })
    fitted_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["fitted.values"]]))
    colnames(fitted_matrix) <- rownames(Y_ordered)
    fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)

    upr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["upr.values"]]))
    colnames(upr_matrix) <- rownames(Y_ordered)
    upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)

    lwr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["lwr.values"]]))
    colnames(lwr_matrix) <- rownames(Y_ordered)
    lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)

    DynamicFeatures <- as.data.frame(do.call(rbind.data.frame, lapply(gam_out, function(x) x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")])))
    char_var <- c("features")
    numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
    DynamicFeatures[, char_var] <- lapply(DynamicFeatures[, char_var, drop = FALSE], as.character)
    DynamicFeatures[, numb_var] <- lapply(DynamicFeatures[, numb_var, drop = FALSE], as.numeric)
    rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
    DynamicFeatures[, "padjust"] <- p.adjust(DynamicFeatures[, "pvalue", drop = TRUE])

    res <- list(
      DynamicFeatures = DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = fitted_matrix,
      upr_matrix = upr_matrix,
      lwr_matrix = lwr_matrix,
      libsize = l_libsize,
      lineages = l,
      family = family
    )
    srt@tools[[paste0("DynamicFeatures_", suffix[i])]] <- res
  }

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "RunDynamicFeatures done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  return(srt)
}

#' RunDynamicEnrichment
#'
#' This function calculates gene-set scores from the specified database (\code{db}) for each lineage using the specified scoring method (\code{score_method}).
#' It then treats these scores as expression values and uses them as input to the RunDynamicFeatures function to identify dynamically enriched terms along the lineage.
#'
#' @inheritParams RunEnrichment
#' @inheritParams DynamicHeatmap
#' @inheritParams CellScoring
#' @param score_method The method to use for scoring. Can be "Seurat", "AUCell", or "UCell". Defaults to "Seurat".
#'
#' @seealso \code{\link{RunDynamicFeatures}} \code{\link{DynamicHeatmap}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = "Lineage1", n_candidates = 200)
#' ht1 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   cell_annotation = "SubCellType",
#'   n_split = 4
#' )
#' ht1$plot
#'
#' pancreas_sub <- RunDynamicEnrichment(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   score_method = "AUCell",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' ht2 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   assay = "GO_BP",
#'   lineages = "Lineage1_GO_BP",
#'   cell_annotation = "SubCellType",
#'   n_split = 4,
#'   split_method = "kmeans-peaktime"
#' )
#' ht2$plot
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @export
RunDynamicEnrichment <- function(srt, lineages,
                                 score_method = "AUCell",
                                 slot = "data", assay = NULL,
                                 min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05,
                                 IDtype = "symbol", species = "Homo_sapiens",
                                 db = "GO_BP", db_update = FALSE, db_version = "latest", convert_species = TRUE,
                                 Ensembl_version = 103, mirror = NULL,
                                 TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                                 BPPARAM = BiocParallel::bpparam(), seed = 11) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  assay <- assay %||% DefaultAssay(srt)

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start RunDynamicFeatures"))
  message("Workers: ", bpworkers(BPPARAM))

  feature_union <- c()
  cell_union <- c()
  dynamic <- list()
  for (l in lineages) {
    if (!paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      stop(l, " info not found in the srt object. Should perform RunDynamicFeatures first!")
    }
    DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
    DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > min_expcells & DynamicFeatures$r.sq > r.sq & DynamicFeatures$dev.expl > dev.expl & DynamicFeatures$padjust < padjust, , drop = FALSE]
    dynamic[[l]] <- DynamicFeatures
    feature_union <- c(feature_union, DynamicFeatures[, "features"])
    cell_union <- c(cell_union, rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]))
  }
  feature_union <- unique(feature_union)

  if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species, db = db, db_update = db_update, db_version = db_version,
      db_IDtypes = IDtype, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    db <- "custom"
    db_list <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
  }

  for (i in seq_along(db)) {
    term <- db[i]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", IDtype)]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], , drop = FALSE]

    term_use <- unique(TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% feature_union, "Term"])
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, "Term"] %in% term_use, , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% term_use, , drop = FALSE]
    rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% rownames(srt[[assay]]), , drop = FALSE]
    feature_list <- split(TERM2GENE_tmp[, IDtype], TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"])
    GSSize <- sapply(feature_list, length)
    feature_list <- feature_list[GSSize >= minGSSize & GSSize <= maxGSSize]

    srt <- CellScoring(
      srt = srt,
      features = feature_list,
      method = score_method,
      classification = FALSE,
      slot = slot,
      assay = assay,
      name = term,
      new_assay = TRUE,
      BPPARAM = BPPARAM
    )
    srt <- RunDynamicFeatures(
      srt = srt,
      lineages = lineages,
      features = rownames(srt[[term]]@counts),
      suffix = paste(lineages, term, sep = "_"),
      assay = term
    )
  }

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "RunDynamicEnrichment done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  return(srt)
}

#' @importFrom reticulate py_to_r
py_to_r_auto <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    x <- py_to_r(x)
  }
  return(x)
}

#' Convert a seurat object to an anndata object using reticulate
#'
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @param srt A Seurat object.
#' @param assay_X Assay to convert as the main data matrix (X) in the anndata object.
#' @param slot_X Slot name for assay_X in the Seurat object.
#' @param assay_layers Assays to convert as layers in the anndata object.
#' @param slot_layers Slot names for the assay_layers in the Seurat object.
#' @param convert_tools Logical indicating whether to convert the tool-specific data.
#' @param convert_misc Logical indicating whether to convert the miscellaneous data.
#' @param features Optional vector of features to include in the anndata object. Defaults to all features in assay_X.
#' @param verbose Logical indicating whether to print verbose messages during the conversion process.
#'
#' @return An anndata object.
#'
#' @return A \code{anndata} object.
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ### Or save as an h5ad file or a loom file
#' # adata$write_h5ad("pancreas_sub.h5ad")
#' # adata$write_loom("pancreas_sub.loom", write_obsm_varm = TRUE)
#'
#' @importFrom reticulate import np_array
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @export
