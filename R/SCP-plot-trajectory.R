LineagePlot <- function(srt, lineages, reduction = NULL, dims = c(1, 2), cells = NULL,
                        trim = c(0.01, 0.99), span = 0.75,
                        palette = "Dark2", palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                        linewidth = 1, line_bg = "white", line_bg_stroke = 0.5,
                        whiskers = FALSE, whiskers_linewidth = 0.5, whiskers_alpha = 0.5,
                        aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                        legend.position = "right", legend.direction = "vertical",
                        theme_use = "theme_scp", theme_args = list(),
                        return_layer = FALSE, seed = 11) {
  set.seed(seed)

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_lineages <- srt@meta.data[, unique(lineages), drop = FALSE]
  dat <- cbind(dat_dim, dat_lineages[row.names(dat_dim), , drop = FALSE])
  dat[, "cell"] <- rownames(dat)
  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
  }

  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  colors <- palette_scp(lineages, palette = palette, palcolor = palcolor)
  axes <- paste0(reduction_key, dims)
  fitted_list <- lapply(lineages, function(l) {
    trim_pass <- dat[[l]] > quantile(dat[[l]], trim[1], na.rm = TRUE) & dat[[l]] < quantile(dat[[l]], trim[2], na.rm = TRUE)
    na_pass <- !is.na(dat[[l]])
    index <- which(trim_pass & na_pass)
    index <- index[order(dat[index, l])]
    dat_sub <- dat[index, , drop = FALSE]
    # if (is.null(weights)) {
    weights_used <- rep(1, nrow(dat_sub))
    # } else {
    # weights_used <- dat_sub[[weights]]
    # }
    fitted <- lapply(axes, function(x) {
      loess(formula(paste(x, l, sep = "~")), weights = weights_used, data = dat_sub, span = span, degree = 2)$fitted
    })
    names(fitted) <- axes
    fitted[["index"]] <- index
    return(fitted)
  })
  names(fitted_list) <- lineages

  curve_layer <- lapply(lineages, function(l) {
    dat_smooth <- as.data.frame(fitted_list[[l]])
    colnames(dat_smooth) <- c(paste0("Axis_", 1:(ncol(dat_smooth) - 1)), "index")
    dat_smooth[, "Lineages"] <- factor(l, levels = lineages)
    dat_smooth <- unique(na.omit(dat_smooth))
    curve <- list()
    if (isTRUE(whiskers)) {
      dat_smooth[, "raw_Axis_1"] <- dat[dat_smooth[, "index"], axes[1]]
      dat_smooth[, "raw_Axis_2"] <- dat[dat_smooth[, "index"], axes[2]]
      curve <- c(curve, geom_segment(
        data = dat_smooth, mapping = aes(x = Axis_1, y = Axis_2, xend = raw_Axis_1, yend = raw_Axis_2, color = Lineages),
        linewidth = whiskers_linewidth, alpha = whiskers_alpha,
        show.legend = TRUE, inherit.aes = FALSE
      ))
    }
    curve <- c(
      curve,
      geom_path(
        data = dat_smooth, mapping = aes(x = Axis_1, y = Axis_2), color = line_bg,
        linewidth = linewidth + line_bg_stroke, arrow = lineages_arrow,
        show.legend = TRUE, inherit.aes = FALSE
      ),
      geom_path(
        data = dat_smooth, mapping = aes(x = Axis_1, y = Axis_2, color = Lineages),
        linewidth = linewidth, arrow = lineages_arrow,
        show.legend = TRUE, inherit.aes = FALSE
      )
    )
    return(curve)
  })
  curve_layer <- c(unlist(curve_layer), list(scale_color_manual(values = colors)))

  lab_layer <- list(labs(title = title, subtitle = subtitle, x = xlab, y = ylab))
  theme_layer <- list(do.call(theme_use, theme_args) +
    theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    ))

  if (isTRUE(return_layer)) {
    return(list(
      curve_layer = curve_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  } else {
    return(ggplot() +
      curve_layer +
      lab_layer +
      theme_layer)
  }
}


#' PAGA plot
#'
#' This function generates a PAGA plot based on the given Seurat object and PAGA result.
#'
#' @param srt A Seurat object containing a PAGA result.
#' @param paga The PAGA result from the Seurat object. Defaults to \code{srt\$misc\$paga}.
#' @param type The type of plot to generate. Possible values are "connectivities" (default) and "connectivities_tree".
#' @param reduction The type of reduction to use for the plot. Defaults to the default reduction in the Seurat object.
#' @param dims The dimensions of the reduction to use for the plot. Defaults to the first two dimensions.
#' @param cells The cells to include in the plot. Defaults to all cells.
#' @param show_transition A logical value indicating whether to display transitions between different cell states. Defaults to \code{FALSE}.
#' @param node_palette The color palette to use for node coloring. Defaults to "Paired".
#' @param node_palcolor A vector of colors to use for node coloring. Defaults to \code{NULL}.
#' @param node_size The size of the nodes in the plot. Defaults to 4.
#' @param node_alpha The transparency of the nodes in the plot. Defaults to 1.
#' @param node_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param node_highlight_color The color to use for highlighting the nodes. Defaults to "red".
#' @param label A logical value indicating whether to display labels for the nodes. Defaults to \code{FALSE}.
#' @param label.size The size of the labels. Defaults to 3.5.
#' @param label.fg The color of the label text. Defaults to "white".
#' @param label.bg The background color of the labels. Defaults to "black".
#' @param label.bg.r The transparency of the label background color. Defaults to 0.1.
#' @param label_insitu A logical value indicating whether to use in-situ labeling for the nodes. Defaults to \code{FALSE}.
#' @param label_repel A logical value indicating whether to use repel mode for labeling nodes. Defaults to \code{FALSE}.
#' @param label_repulsion The repulsion factor for repel mode. Defaults to 20.
#' @param label_point_size The size of the points in the labels. Defaults to 1.
#' @param label_point_color The color of the points in the labels. Defaults to "black".
#' @param label_segment_color The color of the lines connecting the points in the labels. Defaults to "black".
#' @param edge_threshold The threshold for including edges in the plot. Defaults to 0.01.
#' @param edge_line The type of line to use for the edges. Possible values are "straight" and "curved". Defaults to "straight".
#' @param edge_line_curvature The curvature factor for curved edges. Defaults to 0.3.
#' @param edge_line_angle The angle for curved edges. Defaults to 90.
#' @param edge_size The size range for the edges. Defaults to c(0.2, 1).
#' @param edge_color The color of the edges. Defaults to "grey40".
#' @param edge_alpha The transparency of the edges. Defaults to 0.5.
#' @param edge_shorten The amount to shorten the edges. Defaults to 0.
#' @param edge_offset The offset for curved edges. Defaults to 0.
#' @param edge_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param edge_highlight_color The color to use for highlighting the edges. Defaults to "red".
#' @param transition_threshold The threshold for including transitions in the plot. Defaults to 0.01.
#' @param transition_line The type of line to use for the transitions. Possible values are "straight" and "curved". Defaults to "straight".
#' @param transition_line_curvature The curvature factor for curved transitions. Defaults to 0.3.
#' @param transition_line_angle The angle for curved transitions. Defaults to 90.
#' @param transition_size The size range for the transitions. Defaults to c(0.2, 1).
#' @param transition_color The color of the transitions. Defaults to "black".
#' @param transition_alpha The transparency of the transitions. Defaults to 1.
#' @param transition_arrow_type The type of arrow to use for the transitions. Possible values are "closed", "open", and "triangle". Defaults to "closed".
#' @param transition_arrow_angle The angle of the arrow for transitions. Defaults to 20.
#' @param transition_arrow_length The length of the arrow for transitions. Defaults to unit(0.02, "npc").
#' @param transition_shorten The amount to shorten the transitions. Defaults to 0.05.
#' @param transition_offset The offset for curved transitions. Defaults to 0.
#' @param transition_highlight The group(s) to highlight in the plot. Defaults to \code{NULL}.
#' @param transition_highlight_color The color to use for highlighting the transitions. Defaults to "red".
#' @param aspect.ratio The aspect ratio of the plot. Defaults to 1.
#' @param title The title of the plot. Defaults to "PAGA".
#' @param subtitle The subtitle of the plot. Defaults to \code{NULL}.
#' @param xlab The label for the x-axis. Defaults to \code{NULL}.
#' @param ylab The label for the y-axis. Defaults to \code{NULL}.
#' @param legend.position The position of the legend. Possible values are "right", "left", "bottom", and "top". Defaults to "right".
#' @param legend.direction The direction of the legend. Possible values are "vertical" and "horizontal". Defaults to "vertical".
#' @param theme_use The name of the theme to use for the plot. Defaults to "theme_scp".
#' @param theme_args A list of arguments to pass to the theme function. Defaults to an empty list.
#' @param return_layer A logical value indicating whether to return the plot as a ggplot2 layer. Defaults to \code{FALSE}.
#'
#' @seealso \code{\link{RunPAGA}} \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' PAGAPlot(pancreas_sub)
#' PAGAPlot(pancreas_sub, type = "connectivities_tree")
#' PAGAPlot(pancreas_sub, reduction = "PCA")
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
#' PAGAPlot(pancreas_sub, edge_shorten = 0.05)
#' PAGAPlot(pancreas_sub, label = TRUE)
#' PAGAPlot(pancreas_sub, label = TRUE, label_insitu = TRUE)
#' PAGAPlot(pancreas_sub, label = TRUE, label_insitu = TRUE, label_repel = TRUE)
#' PAGAPlot(pancreas_sub, edge_line = "curved")
#' PAGAPlot(pancreas_sub, node_size = "GroupSize")
#' PAGAPlot(pancreas_sub, node_highlight = "Ductal")
#' PAGAPlot(pancreas_sub, edge_highlight = paste("Pre-endocrine", levels(pancreas_sub$SubCellType), sep = "-"))
#'
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' PAGAPlot(pancreas_sub, show_transition = TRUE)
#' PAGAPlot(pancreas_sub, show_transition = TRUE, transition_offset = 0.02)
#'
#' @importFrom Seurat Reductions Key Embeddings
#' @export
PAGAPlot <- function(srt, paga = srt@misc$paga, type = "connectivities",
                     reduction = NULL, dims = c(1, 2), cells = NULL, show_transition = FALSE,
                     node_palette = "Paired", node_palcolor = NULL, node_size = 4, node_alpha = 1,
                     node_highlight = NULL, node_highlight_color = "red",
                     label = FALSE, label.size = 3.5, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                     label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                     label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                     edge_threshold = 0.01, edge_line = c("straight", "curved"), edge_line_curvature = 0.3, edge_line_angle = 90,
                     edge_size = c(0.2, 1), edge_color = "grey40", edge_alpha = 0.5,
                     edge_shorten = 0, edge_offset = 0, edge_highlight = NULL, edge_highlight_color = "red",
                     transition_threshold = 0.01, transition_line = c("straight", "curved"), transition_line_curvature = 0.3, transition_line_angle = 90,
                     transition_size = c(0.2, 1), transition_color = "black", transition_alpha = 1,
                     transition_arrow_type = "closed", transition_arrow_angle = 20, transition_arrow_length = unit(0.02, "npc"),
                     transition_shorten = 0.05, transition_offset = 0, transition_highlight = NULL, transition_highlight_color = "red",
                     aspect.ratio = 1, title = "PAGA", subtitle = NULL, xlab = NULL, ylab = NULL,
                     legend.position = "right", legend.direction = "vertical",
                     theme_use = "theme_scp", theme_args = list(),
                     return_layer = FALSE) {
  if (is.null(paga)) {
    stop("Cannot find the paga result.")
  }
  if (type == "connectivities_tree") {
    use_triangular <- "both"
    edge_threshold <- 0
  } else {
    use_triangular <- "upper"
  }
  connectivities <- paga[[type]]
  transition <- paga[["transitions_confidence"]]
  groups <- paga[["groups"]]
  if (!is.factor(srt@meta.data[[groups]])) {
    srt@meta.data[[groups]] <- factor(srt@meta.data[[groups]])
  }
  if (nlevels(srt@meta.data[[groups]]) != nrow(connectivities)) {
    stop("nlevels in ", groups, " is not identical with the group in paga")
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- as.data.frame(srt@reductions[[reduction]]@cell.embeddings)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_dim <- dat_dim[, paste0(reduction_key, dims)]
  dat_dim[[groups]] <- srt@meta.data[rownames(dat_dim), groups]
  dat <- aggregate(dat_dim[, paste0(reduction_key, dims)], by = list(dat_dim[[groups]]), FUN = median)
  colnames(dat)[1] <- groups
  rownames(dat) <- dat[[groups]]
  dat[["GroupSize"]] <- as.numeric(table(dat_dim[[groups]])[rownames(dat)])
  colnames(connectivities) <- rownames(connectivities) <- rownames(dat)
  if (!is.null(transition)) {
    colnames(transition) <- rownames(transition) <- rownames(dat)
  }
  if (!isTRUE(show_transition)) {
    transition <- NULL
  } else if (isTRUE(show_transition) && is.null(transition)) {
    warning("transitions_confidence need to be calculated first.", immediate. = TRUE)
  }
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])

  if (!is.null(cells)) {
    dat <- dat[intersect(rownames(dat), cells), , drop = FALSE]
    connectivities <- connectivities[rownames(dat), rownames(dat)]
  }

  out <- GraphPlot(
    node = dat, edge = as_matrix(connectivities), node_coord = paste0(reduction_key, dims),
    node_group = groups, node_palette = node_palette, node_palcolor = node_palcolor, node_size = node_size, node_alpha = node_alpha,
    node_highlight = node_highlight, node_highlight_color = node_highlight_color,
    label = label, label.size = label.size, label.fg = label.fg, label.bg = label.bg, label.bg.r = label.bg.r,
    label_insitu = label_insitu, label_repel = label_repel, label_repulsion = label_repulsion,
    label_point_size = label_point_size, label_point_color = label_point_color, label_segment_color = label_segment_color,
    edge_threshold = edge_threshold, use_triangular = use_triangular,
    edge_line = edge_line, edge_line_curvature = edge_line_curvature, edge_line_angle = edge_line_angle, edge_size = edge_size, edge_color = edge_color, edge_alpha = edge_alpha,
    edge_shorten = edge_shorten, edge_offset = edge_offset,
    edge_highlight = edge_highlight, edge_highlight_color = edge_highlight_color,
    transition = transition, transition_threshold = transition_threshold, transition_line = transition_line, transition_line_curvature = transition_line_curvature, transition_line_angle = transition_line_angle,
    transition_color = transition_color, transition_size = transition_size, transition_alpha = transition_alpha,
    transition_arrow_type = transition_arrow_type, transition_arrow_angle = transition_arrow_angle, transition_arrow_length = transition_arrow_length,
    transition_shorten = transition_shorten, transition_offset = transition_offset,
    transition_highlight = transition_highlight, transition_highlight_color = transition_highlight_color,
    aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
    legend.position = legend.position, legend.direction = legend.direction,
    theme_use = theme_use, theme_args = theme_args,
    return_layer = return_layer
  )
  return(out)
}

#' GraphPlot
#'
#' A function to plot a graph with nodes and edges.
#'
#' @param node A data frame representing the nodes of the graph.
#' @param edge A matrix representing the edges of the graph.
#' @param transition A matrix representing the transitions between nodes.
#' @param node_coord A character vector specifying the names of the columns in \code{node} that represent the x and y coordinates.
#' @param node_group A character vector specifying the name of the column in \code{node} that represents the grouping of the nodes.
#' @param node_palette A character vector specifying the name of the color palette for node groups.
#' @param node_palcolor A character vector specifying the names of the colors for each node group.
#' @param node_size A numeric value or column name of \code{node} specifying the size of the nodes.
#' @param node_alpha A numeric value or column name of \code{node} specifying the transparency of the nodes.
#' @param node_highlight A character vector specifying the names of nodes to highlight.
#' @param node_highlight_color A character vector specifying the color for highlighting nodes.
#' @param label A logical value indicating whether to show labels for the nodes.
#' @param label.size A numeric value specifying the size of the labels.
#' @param label.fg A character vector specifying the foreground color of the labels.
#' @param label.bg A character vector specifying the background color of the labels.
#' @param label.bg.r A numeric value specifying the background color transparency of the labels.
#' @param label_insitu A logical value indicating whether to display the node group labels in situ or as numeric values.
#' @param label_repel A logical value indicating whether to use force-directed label repulsion.
#' @param label_repulsion A numeric value specifying the repulsion force for labels.
#' @param label_point_size A numeric value specifying the size of the label points.
#' @param label_point_color A character vector specifying the color of the label points.
#' @param label_segment_color A character vector specifying the color for the label segments.
#' @param edge_threshold A numeric value specifying the threshold for removing edges.
#' @param use_triangular A character vector specifying which part of the edge matrix to use (upper, lower, both).
#' @param edge_line A character vector specifying the type of line for edges (straight, curved).
#' @param edge_line_curvature A numeric value specifying the curvature of curved edges.
#' @param edge_line_angle A numeric value specifying the angle of curved edges.
#' @param edge_color A character vector specifying the color of the edges.
#' @param edge_size A numeric vector specifying the range of edge sizes.
#' @param edge_alpha A numeric value specifying the transparency of the edges.
#' @param edge_shorten A numeric value specifying the length of the edge shorten.
#' @param edge_offset A numeric value specifying the length of the edge offset.
#' @param edge_highlight A character vector specifying the names of edges to highlight.
#' @param edge_highlight_color A character vector specifying the color for highlighting edges.
#' @param transition_threshold A numeric value specifying the threshold for removing transitions.
#' @param transition_line A character vector specifying the type of line for transitions (straight, curved).
#' @param transition_line_curvature A numeric value specifying the curvature of curved transitions.
#' @param transition_line_angle A numeric value specifying the angle of curved transitions.
#' @param transition_color A character vector specifying the color of the transitions.
#' @param transition_size A numeric vector specifying the range of transition sizes.
#' @param transition_alpha A numeric value specifying the transparency of the transitions.
#' @param transition_arrow_type A character vector specifying the type of arrow for transitions (closed, open).
#' @param transition_arrow_angle A numeric value specifying the angle of the transition arrow.
#' @param transition_arrow_length A numeric value specifying the length of the transition arrow.
#' @param transition_shorten A numeric value specifying the length of the transition shorten.
#' @param transition_offset A numeric value specifying the length of the transition offset.
#' @param transition_highlight A character vector specifying the names of transitions to highlight.
#' @param transition_highlight_color A character vector specifying the color for highlighting transitions.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' @param title A character value specifying the title of the plot.
#' @param subtitle A character value specifying the subtitle of the plot.
#' @param xlab A character value specifying the label for the x-axis.
#' @param ylab A character value specifying the label for the y-axis.
#' @param legend.position A character value specifying the position of the legend.
#' @param legend.direction A character value specifying the direction of the legend.
#' @param theme_use A character value specifying the theme to use.
#' @param theme_args A list of arguments to be passed to the theme.
#' @param return_layer A logical value indicating whether to return the layers of the plot instead of the plot itself.
#'
#' @seealso \code{\link{CellDimPlot}}
#'
#' @importFrom ggplot2 scale_size_identity scale_size_continuous scale_size_discrete scale_alpha_identity scale_alpha_continuous scale_alpha_discrete geom_curve geom_segment geom_point scale_color_manual guide_legend guides labs aes scale_linewidth_continuous
#' @importFrom ggnewscale new_scale
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow
#' @importFrom reshape2 melt
#' @export
GraphPlot <- function(node, edge, transition = NULL,
                      node_coord = c("x", "y"), node_group = NULL, node_palette = "Paired", node_palcolor = NULL, node_size = 4, node_alpha = 1,
                      node_highlight = NULL, node_highlight_color = "red",
                      label = FALSE, label.size = 3.5, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                      label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                      label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                      edge_threshold = 0.01, use_triangular = c("upper", "lower", "both"), edge_line = c("straight", "curved"), edge_line_curvature = 0.3, edge_line_angle = 90,
                      edge_color = "grey40", edge_size = c(0.2, 1), edge_alpha = 0.5,
                      edge_shorten = 0, edge_offset = 0, edge_highlight = NULL, edge_highlight_color = "red",
                      transition_threshold = 0.01, transition_line = c("straight", "curved"), transition_line_curvature = 0.3, transition_line_angle = 90,
                      transition_color = "black", transition_size = c(0.2, 1), transition_alpha = 1,
                      transition_arrow_type = "closed", transition_arrow_angle = 20, transition_arrow_length = unit(0.02, "npc"),
                      transition_shorten = 0.05, transition_offset = 0, transition_highlight = NULL, transition_highlight_color = "red",
                      aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                      legend.position = "right", legend.direction = "vertical",
                      theme_use = "theme_scp", theme_args = list(),
                      return_layer = FALSE) {
  use_triangular <- match.arg(use_triangular)
  edge_line <- match.arg(edge_line)
  transition_line <- match.arg(transition_line)
  if (!is.data.frame(node)) {
    stop("'node' must be a data.frame object.")
  }
  if (!is.matrix(edge)) {
    stop("'edge' must be a matrix object.")
  }
  if (!identical(nrow(edge), ncol(edge))) {
    stop("nrow and ncol is not identical in edge matrix")
  }
  if (!identical(nrow(edge), nrow(node))) {
    stop("nrow is not identical between edge and node.")
  }
  if (!identical(rownames(edge), rownames(node))) {
    warning("rownames of node is not identical with edge matrix. They will correspond according to the order.", immediate. = TRUE)
    colnames(edge) <- rownames(edge) <- rownames(node) <- rownames(node) %||% colnames(edge) %||% rownames(edge)
  }
  if (!all(node_coord %in% colnames(node))) {
    stop("Cannot find the node_coord ", paste(node_coord[!node_coord %in% colnames(node)], collapse = ","), " in the node column")
  }
  if (!is.null(transition)) {
    if (!identical(nrow(transition), nrow(node))) {
      stop("nrow is not identical between transition and node.")
    }
    if (!identical(rownames(transition), rownames(node))) {
      warning("rownames of node is not identical with transition matrix. They will correspond according to the order.", immediate. = TRUE)
      colnames(transition) <- rownames(transition) <- rownames(node) <- rownames(node) %||% colnames(transition) %||% rownames(transition)
    }
  }
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  node <- as.data.frame(node)
  node[["x"]] <- node[[node_coord[1]]]
  node[["y"]] <- node[[node_coord[2]]]
  node[["node_name"]] <- rownames(node)
  node_group <- node_group %||% "node_name"
  node_size <- node_size %||% 5
  node_alpha <- node_alpha %||% 1
  if (!node_group %in% colnames(node)) {
    node[["node_group"]] <- node_group
  } else {
    node[["node_group"]] <- node[[node_group]]
  }
  if (!is.factor(node[["node_group"]])) {
    node[["node_group"]] <- factor(node[["node_group"]], levels = unique(node[["node_group"]]))
  }
  if (!node_size %in% colnames(node)) {
    if (!is.numeric(node_size)) {
      node_size <- 5
    }
    node[["node_size"]] <- node_size
    scale_size <- scale_size_identity()
  } else {
    node[["node_size"]] <- node[[node_size]]
    if (is.numeric(node[[node_size]])) {
      scale_size <- scale_size_continuous(name = node_size)
    } else {
      scale_size <- scale_size_discrete()
    }
  }
  if (!node_alpha %in% colnames(node)) {
    if (!is.numeric(node_alpha)) {
      node_alpha <- 1
    }
    node[["node_alpha"]] <- node_alpha
    scale_alpha <- scale_alpha_identity()
  } else {
    node[["node_alpha"]] <- node[[node_alpha]]
    if (is.numeric(node[[node_alpha]])) {
      scale_alpha <- scale_alpha_continuous()
    } else {
      scale_alpha <- scale_alpha_discrete()
    }
  }

  if (isTRUE(label) && !isTRUE(label_insitu)) {
    label_use <- paste0(1:nlevels(node[["node_group"]]), ": ", levels(node[["node_group"]]))
  } else {
    label_use <- levels(node[["node_group"]])
  }
  global_size <- sqrt(max(node[["x"]], na.rm = TRUE)^2 + max(node[["y"]], na.rm = TRUE)^2)

  edge[edge <= edge_threshold] <- NA
  if (use_triangular == "upper") {
    edge[lower.tri(edge)] <- NA
  } else if (use_triangular == "lower") {
    edge[upper.tri(edge)] <- NA
  }
  edge_df <- melt(edge, na.rm = TRUE, stringsAsFactors = FALSE)
  if (nrow(edge_df) == 0) {
    edge_layer <- NULL
  } else {
    colnames(edge_df) <- c("from", "to", "size")
    edge_df[, "from"] <- as.character(edge_df[, "from"])
    edge_df[, "to"] <- as.character(edge_df[, "to"])
    edge_df[, "size"] <- as.numeric(edge_df[, "size"])
    edge_df[, "x"] <- node[edge_df[, "from"], "x"]
    edge_df[, "y"] <- node[edge_df[, "from"], "y"]
    edge_df[, "xend"] <- node[edge_df[, "to"], "x"]
    edge_df[, "yend"] <- node[edge_df[, "to"], "y"]
    rownames(edge_df) <- edge_df[, "edge_name"] <- paste0(edge_df[, "from"], "-", edge_df[, "to"])
    edge_df <- segementsDf(edge_df, global_size * edge_shorten, global_size * edge_shorten, global_size * edge_offset)

    linetype <- ifelse(is.null(transition), 1, 2)
    if (edge_line == "straight") {
      edge_layer <- list(
        geom_segment(
          data = edge_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
          lineend = "round", linejoin = "mitre", linetype = linetype, color = edge_color, alpha = edge_alpha,
          inherit.aes = FALSE, show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[edge_df[["edge_name"]] %in% edge_highlight, , drop = FALSE]
        edge_layer <- c(
          edge_layer,
          list(
            geom_segment(
              data = edge_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
              lineend = "round", linejoin = "mitre", linetype = linetype, color = edge_highlight_color, alpha = 1,
              inherit.aes = FALSE, show.legend = FALSE
            )
          )
        )
      }
    } else {
      edge_layer <- list(
        geom_curve(
          data = edge_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
          curvature = edge_line_curvature, angle = edge_line_angle,
          lineend = "round", linetype = linetype, color = edge_color, alpha = edge_alpha,
          inherit.aes = FALSE, show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[edge_df[["edge_name"]] %in% edge_highlight, , drop = FALSE]
        edge_layer <- c(
          edge_layer,
          list(
            geom_curve(
              data = edge_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
              curvature = edge_line_curvature, angle = edge_line_angle,
              lineend = "round", linetype = linetype, color = edge_highlight_color, alpha = 1,
              inherit.aes = FALSE, show.legend = FALSE
            )
          )
        )
      }
    }
    edge_layer <- c(edge_layer, list(
      scale_linewidth_continuous(range = range(edge_size), guide = "none"),
      new_scale("linewidth")
    ))
  }

  if (!is.null(transition)) {
    trans2 <- trans1 <- as_matrix(transition)
    trans1[lower.tri(trans1)] <- 0
    trans2[upper.tri(trans2)] <- 0
    trans <- t(trans1) - trans2
    trans[abs(trans) <= transition_threshold] <- NA
    trans_df <- melt(trans, na.rm = TRUE, stringsAsFactors = FALSE)
    if (nrow(trans_df) == 0) {
      trans_layer <- NULL
    } else {
      trans_df <- as.data.frame(t(apply(trans_df, 1, function(x) {
        if (as.numeric(x[3]) < 0) {
          return(c(x[c(2, 1)], -as.numeric(x[3])))
        } else {
          return(x)
        }
      })))
      colnames(trans_df) <- c("from", "to", "size")
      trans_df[, "from"] <- as.character(trans_df[, "from"])
      trans_df[, "to"] <- as.character(trans_df[, "to"])
      trans_df[, "size"] <- as.numeric(trans_df[, "size"])
      trans_df[, "x"] <- node[trans_df[, "from"], "x"]
      trans_df[, "y"] <- node[trans_df[, "from"], "y"]
      trans_df[, "xend"] <- node[trans_df[, "to"], "x"]
      trans_df[, "yend"] <- node[trans_df[, "to"], "y"]
      rownames(trans_df) <- trans_df[, "trans_name"] <- paste0(trans_df[, "from"], "-", trans_df[, "to"])
      trans_df <- segementsDf(trans_df, global_size * transition_shorten, global_size * transition_shorten, global_size * transition_offset)

      if (transition_line == "straight") {
        trans_layer <- list(
          geom_segment(
            data = trans_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
            arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
            lineend = "round", linejoin = "mitre", color = transition_color, alpha = transition_alpha,
            inherit.aes = FALSE, show.legend = FALSE
          )
        )
        if (!is.null(transition_highlight)) {
          trans_df_highlight <- trans_df[trans_df[["trans_name"]] %in% transition_highlight, , drop = FALSE]
          trans_layer <- c(
            trans_layer,
            list(
              geom_segment(
                data = trans_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
                arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
                lineend = "round", linejoin = "mitre", color = transition_highlight_color, alpha = 1,
                inherit.aes = FALSE, show.legend = FALSE
              )
            )
          )
        }
      } else {
        trans_layer <- list(
          geom_curve(
            data = trans_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
            arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
            curvature = transition_line_curvature, angle = transition_line_angle,
            lineend = "round", color = transition_color, alpha = transition_alpha,
            inherit.aes = FALSE, show.legend = FALSE
          )
        )
        if (!is.null(edge_highlight)) {
          trans_df_highlight <- trans_df[trans_df[["trans_name"]] %in% transition_highlight, , drop = FALSE]
          trans_layer <- c(
            trans_layer,
            list(
              geom_curve(
                data = trans_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
                arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
                curvature = transition_line_curvature, angle = transition_line_angle,
                lineend = "round", color = transition_highlight_color, alpha = 1,
                inherit.aes = FALSE, show.legend = FALSE
              )
            )
          )
        }
      }
      trans_layer <- c(trans_layer, list(
        scale_linewidth_continuous(range = range(transition_size), guide = "none"),
        new_scale("linewidth")
      ))
    }
  } else {
    trans_layer <- NULL
  }

  node_layer <- list(
    geom_point(data = node, aes(x = x, y = y, size = node_size * 1.2), color = "black", show.legend = FALSE, inherit.aes = FALSE),
    geom_point(data = node, aes(x = x, y = y, size = node_size, color = node_group, alpha = node_alpha), inherit.aes = FALSE)
  )
  if (!is.null(node_highlight)) {
    node_highlight <- node[node[["node_name"]] %in% node_highlight, , drop = FALSE]
    node_layer <- c(
      node_layer,
      list(
        geom_point(data = node_highlight, aes(x = x, y = y, size = node_size * 1.3), color = node_highlight_color, show.legend = FALSE, inherit.aes = FALSE),
        geom_point(data = node_highlight, aes(x = x, y = y, size = node_size, color = node_group, alpha = node_alpha), inherit.aes = FALSE)
      )
    )
  }
  node_layer <- c(node_layer, list(
    scale_color_manual(
      name = node_group, values = palette_scp(node[["node_group"]], palette = node_palette, palcolor = node_palcolor), labels = label_use,
      guide = guide_legend(
        title.hjust = 0,
        order = 1,
        override.aes = list(size = 4, alpha = 1)
      )
    ),
    scale_size,
    scale_alpha
  ))

  if (isTRUE(label)) {
    if (isTRUE(label_insitu)) {
      node[, "label"] <- node[["node_group"]]
    } else {
      node[, "label"] <- as.numeric(node[["node_group"]])
    }
    if (isTRUE(label_repel)) {
      label_layer <- list(geom_point(
        data = node, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
        color = label_point_color, size = label_point_size, inherit.aes = FALSE
      ), geom_text_repel(
        data = node, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
        fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
        point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
        color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
      ))
    } else {
      label_layer <- list(geom_text_repel(
        data = node, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
        fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
        point.size = NA, max.overlaps = 100, force = 0,
        color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
      ))
    }
  } else {
    label_layer <- NULL
  }

  lab_layer <- list(labs(title = title, subtitle = subtitle, x = xlab, y = ylab))
  theme_layer <- list(do.call(theme_use, theme_args) +
    theme(
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction
    ))

  if (isTRUE(return_layer)) {
    return(list(
      edge_layer = edge_layer,
      trans_layer = trans_layer,
      node_layer = node_layer,
      label_layer = label_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  } else {
    return(ggplot() +
      edge_layer +
      trans_layer +
      node_layer +
      label_layer +
      lab_layer +
      theme_layer)
  }
}

#' Shorten and offset the segment
#'
#' This function takes a data frame representing segments in a plot and shortens
#' and offsets them based on the provided arguments.
#'
#' @param data A data frame containing the segments. It should have columns
#'             'x', 'y', 'xend', and 'yend' representing the start and end points
#'             of each segment.
#' @param shorten_start The amount to shorten the start of each segment by.
#' @param shorten_end The amount to shorten the end of each segment by.
#' @param offset The amount to offset each segment by.
#'
#' @return The modified data frame with the shortened and offset segments.
#'
#' @examples
#' library(ggplot2)
#' tempNodes <- data.frame("x" = c(10, 40), "y" = c(10, 30))
#' data <- data.frame("x" = c(10, 40), "y" = c(10, 30), "xend" = c(40, 10), "yend" = c(30, 10))
#' ggplot(tempNodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(data = data, aes(x = x, xend = xend, y = y, yend = yend))
#'
#' ggplot(tempNodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(
#'     data = segementsDf(data, shorten_start = 2, shorten_end = 3, offset = 1),
#'     aes(x = x, xend = xend, y = y, yend = yend)
#'   )
#' @export
segementsDf <- function(data, shorten_start, shorten_end, offset) {
  data$dx <- data$xend - data$x
  data$dy <- data$yend - data$y
  data$dist <- sqrt(data$dx^2 + data$dy^2)
  data$px <- data$dx / data$dist
  data$py <- data$dy / data$dist

  data$x <- data$x + data$px * shorten_start
  data$y <- data$y + data$py * shorten_start
  data$xend <- data$xend - data$px * shorten_end
  data$yend <- data$yend - data$py * shorten_end
  data$x <- data$x - data$py * offset
  data$xend <- data$xend - data$py * offset
  data$y <- data$y + data$px * offset
  data$yend <- data$yend + data$px * offset

  return(data)
}

#' Velocity Plot
#'
#' This function creates a velocity plot for a given Seurat object. The plot shows the velocity vectors of the cells in a specified reduction space.
#'
#' @param srt A Seurat object.
#' @param reduction Name of the reduction in the Seurat object to use for plotting.
#' @param dims Indices of the dimensions to use for plotting.
#' @param cells Cells to include in the plot. If NULL, all cells will be included.
#' @param velocity Name of the velocity to use for plotting.
#' @param plot_type Type of plot to create. Can be "raw", "grid", or "stream".
#' @param group_by Name of the column in the Seurat object metadata to group the cells by. Defaults to NULL.
#' @param group_palette Name of the palette to use for coloring the groups. Defaults to "Paired".
#' @param group_palcolor Colors to use for coloring the groups. Defaults to NULL.
#' @param n_neighbors Number of neighbors to include for the density estimation. Defaults to ceiling(ncol(srt@assays[[1]]) / 50).
#' @param density Propotion of cells to plot. Defaults to 1 (plot all cells).
#' @param smooth Smoothing parameter for density estimation. Defaults to 0.5.
#' @param scale Scaling factor for the velocity vectors. Defaults to 1.
#' @param min_mass Minimum mass value for the density-based cutoff. Defaults to 1.
#' @param cutoff_perc Percentile value for the density-based cutoff. Defaults to 5.
#' @param arrow_angle Angle of the arrowheads. Defaults to 20.
#' @param arrow_color Color of the arrowheads. Defaults to "black".
#' @param streamline_L Length of the streamlines. Defaults to 5.
#' @param streamline_minL Minimum length of the streamlines. Defaults to 1.
#' @param streamline_res Resolution of the streamlines. Defaults to 1.
#' @param streamline_n Number of streamlines to plot. Defaults to 15.
#' @param streamline_width Width of the streamlines. Defaults to c(0, 0.8).
#' @param streamline_alpha Alpha transparency of the streamlines. Defaults to 1.
#' @param streamline_color Color of the streamlines. Defaults to NULL.
#' @param streamline_palette Name of the palette to use for coloring the streamlines. Defaults to "RdYlBu".
#' @param streamline_palcolor Colors to use for coloring the streamlines. Defaults to NULL.
#' @param streamline_bg_color Background color of the streamlines. Defaults to "white".
#' @param streamline_bg_stroke Stroke width of the streamlines background. Defaults to 0.5.
#' @param aspect.ratio Aspect ratio of the plot. Defaults to 1.
#' @param title Title of the plot. Defaults to "Cell velocity".
#' @param subtitle Subtitle of the plot. Defaults to NULL.
#' @param xlab x-axis label. Defaults to NULL.
#' @param ylab y-axis label. Defaults to NULL.
#' @param legend.position Position of the legend. Defaults to "right".
#' @param legend.direction Direction of the legend. Defaults to "vertical".
#' @param theme_use Name of the theme to use for plotting. Defaults to "theme_scp".
#' @param theme_args List of theme arguments for customization. Defaults to list().
#' @param return_layer Whether to return the plot layers as a list. Defaults to FALSE.
#' @param seed Random seed for reproducibility. Defaults to 11.
#'
#' @seealso \code{\link{RunSCVELO}} \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' VelocityPlot(pancreas_sub, reduction = "UMAP")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", group_by = "SubCellType")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "grid")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream", streamline_color = "black")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream", streamline_color = "black", arrow_color = "red")
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom ggplot2 ggplot aes geom_segment scale_size scale_alpha scale_color_gradientn scale_color_manual guide_legend
#' @importFrom grid arrow unit
#' @importFrom reshape2 melt
#' @export
VelocityPlot <- function(srt, reduction, dims = c(1, 2), cells = NULL, velocity = "stochastic", plot_type = c("raw", "grid", "stream"),
                         group_by = NULL, group_palette = "Paired", group_palcolor = NULL,
                         n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50), density = 1, smooth = 0.5, scale = 1, min_mass = 1, cutoff_perc = 5,
                         arrow_angle = 20, arrow_color = "black",
                         streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
                         streamline_width = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
                         streamline_bg_color = "white", streamline_bg_stroke = 0.5,
                         aspect.ratio = 1, title = "Cell velocity", subtitle = NULL, xlab = NULL, ylab = NULL,
                         legend.position = "right", legend.direction = "vertical",
                         theme_use = "theme_scp", theme_args = list(),
                         return_layer = FALSE, seed = 11) {
  set.seed(seed)

  plot_type <- match.arg(plot_type)
  check_R("metR")

  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  V_reduction <- paste0(velocity, "_", reduction)
  if (!V_reduction %in% names(srt@reductions)) {
    stop("Cannot find the velocity embedding ", V_reduction, ".")
  }
  X_emb <- srt@reductions[[reduction]]@cell.embeddings[, dims]
  V_emb <- srt@reductions[[V_reduction]]@cell.embeddings[, dims]
  if (!is.null(cells)) {
    X_emb <- X_emb[intersect(rownames(X_emb), cells), , drop = FALSE]
    V_emb <- V_emb[intersect(rownames(V_emb), cells), , drop = FALSE]
  }

  reduction_key <- srt@reductions[[reduction]]@key
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  if (plot_type == "raw") {
    if (!is.null(density) && (density > 0 && density < 1)) {
      s <- ceiling(density * nrow(X_emb))
      ix_choice <- sample(seq_len(nrow(X_emb)), size = s, replace = FALSE)
      X_emb <- X_emb[ix_choice, ]
      V_emb <- V_emb[ix_choice, ]
    }
    if (!is.null(scale)) {
      V_emb <- V_emb * scale
    }
    df_field <- cbind.data.frame(X_emb, V_emb)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(max(df_field[["x"]], na.rm = TRUE)^2 + max(df_field[["y"]], na.rm = TRUE)^2)
    df_field[["length_perc"]] <- df_field[["length"]] / global_size

    if (!is.null(group_by)) {
      df_field[["group_by"]] <- srt@meta.data[rownames(df_field), group_by, drop = TRUE]
      velocity_layer <- list(
        geom_segment(
          data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v, color = group_by),
          arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle),
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        scale_color_manual(
          name = group_by, values = palette_scp(df_field[["group_by"]], palette = group_palette, palcolor = group_palcolor),
          guide = guide_legend(title.hjust = 0, order = 1, override.aes = list(linewidth = 2, alpha = 1))
        )
      )
    } else {
      velocity_layer <- list(
        geom_segment(
          data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v),
          color = arrow_color,
          arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle),
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        )
      )
    }
  }
  if (plot_type == "grid") {
    res <- compute_velocity_on_grid(X_emb, V_emb,
      density = density, smooth = smooth, n_neighbors = n_neighbors,
      min_mass = min_mass, scale = scale
    )
    X_grid <- res$X_grid
    V_grid <- res$V_grid

    df_field <- cbind.data.frame(X_grid, V_grid)
    colnames(df_field) <- c("x", "y", "u", "v")
    df_field[["length"]] <- sqrt(df_field[["u"]]^2 + df_field[["v"]]^2)
    global_size <- sqrt(max(df_field[["x"]], na.rm = TRUE)^2 + max(df_field[["y"]], na.rm = TRUE)^2)
    df_field[["length_perc"]] <- df_field[["length"]] / global_size
    velocity_layer <- list(
      geom_segment(
        data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v),
        color = arrow_color,
        arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle),
        lineend = "round", linejoin = "mitre", inherit.aes = FALSE
      )
    )
  }
  if (plot_type == "stream") {
    check_R("metR")
    res <- compute_velocity_on_grid(X_emb, V_emb,
      density = density, smooth = smooth, n_neighbors = n_neighbors,
      min_mass = min_mass, scale = 1, cutoff_perc = cutoff_perc,
      adjust_for_stream = TRUE
    )
    X_grid <- res$X_grid
    V_grid <- res$V_grid

    # if (!is.null(density) && (density > 0 & density < 1)) {
    #   s <- ceiling(density * ncol(X_grid))
    #   ix_choice <- sample(1:ncol(X_grid), size = s, replace = FALSE)
    #   X_grid <- X_grid[, ix_choice]
    #   V_grid <- V_grid[, ix_choice, ix_choice]
    # }

    df_field <- expand.grid(X_grid[1, ], X_grid[2, ])
    colnames(df_field) <- c("x", "y")
    u <- melt(t(V_grid[1, , ]))
    v <- melt(t(V_grid[2, , ]))
    df_field[, "u"] <- u$value
    df_field[, "v"] <- v$value
    df_field[is.na(df_field)] <- 0

    if (!is.null(streamline_color)) {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, size = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, size = max(streamline_width, na.rm = TRUE), color = streamline_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, linetype = 0, color = arrow_color,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        )
      )
    } else {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, size = max(streamline_width, na.rm = TRUE) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v, size = after_stat(step), color = sqrt(after_stat(dx)^2 + after_stat(dy)^2)),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, alpha = streamline_alpha,
          arrow = NULL, lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, linetype = 0, color = arrow_color,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        scale_color_gradientn(
          name = "Velocity", colors = palette_scp(palette = streamline_palette, palcolor = streamline_palcolor),
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
        ),
        scale_size(range = range(streamline_width), guide = "none")
      )
    }
  }

  lab_layer <- list(labs(title = title, subtitle = subtitle, x = xlab, y = ylab))
  theme_layer <- list(do.call(theme_use, theme_args) + theme(
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction
  ))

  if (isTRUE(return_layer)) {
    return(list(
      velocity_layer = velocity_layer,
      lab_layer = lab_layer,
      theme_layer = theme_layer
    ))
  } else {
    return(ggplot() +
      velocity_layer +
      lab_layer +
      lab_layer +
      theme_layer)
  }
}

#' Compute velocity on grid
#' The original python code is on https://github.com/theislab/scvelo/blob/master/scvelo/plotting/velocity_embedding_grid.py
#'
#' @param X_emb A matrix of dimension n_obs x n_dim specifying the embedding coordinates of the cells.
#' @param V_emb A matrix of dimension n_obs x n_dim specifying the velocity vectors of the cells.
#' @param density An optional numeric value specifying the density of the grid points along each dimension. Default is 1.
#' @param smooth An optional numeric value specifying the smoothing factor for the velocity vectors. Default is 0.5.
#' @param n_neighbors An optional numeric value specifying the number of nearest neighbors for each grid point. Default is ceiling(n_obs / 50).
#' @param min_mass An optional numeric value specifying the minimum mass required for a grid point to be considered. Default is 1.
#' @param scale An optional numeric value specifying the scaling factor for the velocity vectors. Default is 1.
#' @param adjust_for_stream A logical value indicating whether to adjust the velocity vectors for streamlines. Default is FALSE.
#' @param cutoff_perc An optional numeric value specifying the percentile cutoff for removing low-density grid points. Default is 5.
#'
#' @importFrom SeuratObject as.sparse
#' @importFrom proxyC dist
#' @export
compute_velocity_on_grid <- function(X_emb, V_emb,
                                     density = NULL, smooth = NULL, n_neighbors = NULL, min_mass = NULL,
                                     scale = 1, adjust_for_stream = FALSE, cutoff_perc = NULL) {
  n_obs <- nrow(X_emb)
  n_dim <- ncol(X_emb)

  density <- density %||% 1
  smooth <- smooth %||% 0.5
  n_neighbors <- n_neighbors %||% ceiling(n_obs / 50)
  min_mass <- min_mass %||% 1
  cutoff_perc <- cutoff_perc %||% 5

  grs <- list()
  for (dim_i in 1:n_dim) {
    m <- min(X_emb[, dim_i], na.rm = TRUE)
    M <- max(X_emb[, dim_i], na.rm = TRUE)
    # m <- m - 0.01 * abs(M - m)
    # M <- M + 0.01 * abs(M - m)
    gr <- seq(m, M, length.out = ceiling(50 * density))
    grs <- c(grs, list(gr))
  }
  X_grid <- as_matrix(expand.grid(grs))

  d <- dist(
    x = as.sparse(X_emb),
    y = as.sparse(X_grid),
    method = "euclidean",
    use_nan = TRUE
  )
  neighbors <- t(as_matrix(apply(d, 2, function(x) order(x, decreasing = FALSE)[1:n_neighbors])))
  dists <- t(as_matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1:n_neighbors]])))

  # ggplot() +
  #   annotate(geom = "point", x = X_grid[, 1], y = X_grid[, 2]) +
  #   annotate(geom = "point", x = X_grid[1, 1], y = X_grid[1, 2], color = "blue") +
  #   annotate(geom = "point", x = X_grid[neighbors[1, ], 1], y = X_grid[neighbors[1, ], 2], color = "red")

  weight <- dnorm(dists, sd = mean(sapply(grs, function(g) g[2] - g[1])) * smooth)
  p_mass <- p_mass_V <- rowSums(weight)
  p_mass_V[p_mass_V < 1] <- 1

  # qplot(dists[,1],weight[,1])
  # qplot(py$dists[,1],py$weight[,1])

  neighbors_emb <- array(V_emb[neighbors, seq_len(ncol(V_emb))],
    dim = c(dim(neighbors), dim(V_emb)[2])
  )
  V_grid <- apply((neighbors_emb * c(weight)), c(1, 3), sum)
  V_grid <- V_grid / p_mass_V

  # qplot(V_grid[,1],V_grid[,2])
  # qplot(py$V_grid[,1],py$V_grid[,2])

  if (isTRUE(adjust_for_stream)) {
    X_grid <- matrix(c(unique(X_grid[, 1]), unique(X_grid[, 2])), nrow = 2, byrow = TRUE)
    ns <- floor(sqrt(length(V_grid[, 1])))
    V_grid <- reticulate::array_reshape(t(V_grid), c(2, ns, ns))

    mass <- sqrt(apply(V_grid**2, c(2, 3), sum))
    min_mass <- 10**(min_mass - 6) # default min_mass = 1e-5
    min_mass[min_mass > max(mass, na.rm = TRUE) * 0.9] <- max(mass, na.rm = TRUE) * 0.9
    cutoff <- reticulate::array_reshape(mass, dim = c(ns, ns)) < min_mass

    length <- t(apply(apply(abs(neighbors_emb), c(1, 3), mean), 1, sum))
    length <- reticulate::array_reshape(length, dim = c(ns, ns))
    cutoff <- cutoff | length < quantile(length, cutoff_perc / 100)
    V_grid[1, , ][cutoff] <- NA
  } else {
    min_mass <- min_mass * quantile(p_mass, 0.99) / 100
    X_grid <- X_grid[p_mass > min_mass, ]
    V_grid <- V_grid[p_mass > min_mass, ]
    if (!is.null(scale)) {
      V_grid <- V_grid * scale
    }
  }
  return(list(X_grid = X_grid, V_grid = V_grid))
}

#' VolcanoPlot
#'
#' Generate a volcano plot based on differential expression analysis results.
#'
#' @param srt An object of class `SummarizedExperiment` containing the results of differential expression analysis.
#' @param group_by A character vector specifying the column in `srt` to group the samples by. Default is `NULL`.
#' @param test.use A character string specifying the type of statistical test to use. Default is "wilcox".
#' @param DE_threshold A character string specifying the threshold for differential expression. Default is "avg_log2FC > 0 & p_val_adj < 0.05".
#' @param x_metric A character string specifying the metric to use for the x-axis. Default is "diff_pct".
#' @param palette A character string specifying the color palette to use for the plot. Default is "RdBu".
#' @param palcolor A character string specifying the color for the palette. Default is `NULL`.
#' @param pt.size A numeric value specifying the size of the points. Default is 1.
#' @param pt.alpha A numeric value specifying the transparency of the points. Default is 1.
#' @param cols.highlight A character string specifying the color for highlighted points. Default is "black".
#' @param sizes.highlight A numeric value specifying the size of the highlighted points. Default is 1.
#' @param alpha.highlight A numeric value specifying the transparency of the highlighted points. Default is 1.
#' @param stroke.highlight A numeric value specifying the stroke width for the highlighted points. Default is 0.5.
#' @param nlabel An integer value specifying the number of labeled points per group. Default is 5.
#' @param features_label A character vector specifying the feature labels to plot. Default is `NULL`.
#' @param label.fg A character string specifying the color for the labels' foreground. Default is "black".
#' @param label.bg A character string specifying the color for the labels' background. Default is "white".
#' @param label.bg.r A numeric value specifying the radius of the rounding of the labels' background. Default is 0.1.
#' @param label.size A numeric value specifying the size of the labels. Default is 4.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot. Default is `NULL`.
#' @param xlab A character string specifying the x-axis label. Default is the value of `x_metric`.
#' @param ylab A character string specifying the y-axis label. Default is "-log10(p-adjust)".
#' @param theme_use A character string specifying the theme to use for the plot. Default is "theme_scp".
#' @param theme_args A list of theme arguments to pass to the `theme_use` function. Default is an empty list.
#' @param combine A logical value indicating whether to combine the plots for each group into a single plot. Default is `TRUE`.
#' @param nrow An integer value specifying the number of rows in the combined plot. Default is `NULL`.
#' @param ncol An integer value specifying the number of columns in the combined plot. Default is `NULL`.
#' @param byrow A logical value indicating whether to arrange the plots by row in the combined plot. Default is `TRUE`.
#'
#' @seealso \code{\link{RunDEtest}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType", only.pos = FALSE)
#' VolcanoPlot(pancreas_sub, group_by = "CellType")
#' VolcanoPlot(pancreas_sub, group_by = "CellType", DE_threshold = "abs(diff_pct) > 0.3 & p_val_adj < 0.05")
#' VolcanoPlot(pancreas_sub, group_by = "CellType", x_metric = "avg_log2FC")
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline labs scale_color_gradientn guide_colorbar facet_wrap position_jitter
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots
#' @export
