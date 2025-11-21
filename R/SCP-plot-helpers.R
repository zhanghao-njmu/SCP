drop_data <- function(p) {
  UseMethod(generic = "drop_data", object = p)
}

#' @param p A \code{ggplot} object.
#' @export
#' @rdname drop_data
#' @method drop_data ggplot
drop_data.ggplot <- function(p) {
  p <- ggplot2:::plot_clone(p)

  # fix the scales for x/y axis and 'fill', 'color', 'shape',...
  for (i in seq_along(p$scales$scales)) {
    if (inherits(p$scales$scales[[i]], "ScaleDiscrete")) {
      p$scales$scales[[i]][["drop"]] <- FALSE
    }
    if (inherits(p$scales$scales[[i]], "ScaleContinuous")) {
      limits <- p$scales$scales[[i]]$get_limits()
      if (p$scales$scales[[i]]$aesthetics[1] == "x") {
        p$coordinates$limits$x <- limits
      }
      if (p$scales$scales[[i]]$aesthetics[1] == "y") {
        p$coordinates$limits$y <- limits
      }
    }
  }

  vars <- get_vars(p)
  # drop main data
  if (length(p$data) > 0) {
    vars_modified <- names(which(sapply(p$data[, intersect(colnames(p$data), vars), drop = FALSE], class) == "character"))
    for (v in vars_modified) {
      p$data[[v]] <- as.factor(p$data[[v]])
    }
    p$data <- p$data[1, , drop = FALSE]
  }

  # drop layer data
  for (i in seq_along(p$layers)) {
    if (length(p$layers[[i]]$data) > 0) {
      vars_modified <- names(which(sapply(p$layers[[i]]$data[, intersect(colnames(p$layers[[i]]$data), vars), drop = FALSE], class) == "character"))
      for (v in vars_modified) {
        p$layers[[i]]$data[[v]] <- as.factor(p$layers[[i]]$data[[v]])
      }
      p$layers[[i]]$data <- p$layers[[i]]$data[1, , drop = FALSE]
    }
  }

  return(p)
}

#' @param p A \code{patchwork} object.
#' @export
#' @rdname drop_data
#' @method drop_data patchwork
drop_data.patchwork <- function(p) {
  for (i in seq_along(p$patches$plots)) {
    p$patches$plots[[i]] <- drop_data(p$patches$plots[[i]])
  }
  p <- drop_data.ggplot(p)
  return(p)
}

#' @export
#' @rdname drop_data
#' @method drop_data default
drop_data.default <- function(p) {
  return(p)
}

#' Drop unused data from the plot to reduce the object size
#'
#' @param p A \code{ggplot} object or a \code{patchwork} object.
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
#'   geom_point()
#' object.size(p)
#' colnames(p$data)
#'
#' p_slim <- slim_data(p)
#' object.size(p_slim)
#' colnames(p_slim$data)
#'
#' @export
slim_data <- function(p) {
  UseMethod(generic = "slim_data", object = p)
}

#' @param p A \code{ggplot} object.
#' @export
#' @rdname slim_data
#' @method slim_data ggplot
slim_data.ggplot <- function(p) {
  vars <- get_vars(p)
  if (length(vars) > 0) {
    p$data <- p$data[, intersect(colnames(p$data), vars), drop = FALSE]
    for (i in seq_along(p$layers)) {
      if (length(p$layers[[i]]$data) > 0) {
        p$layers[[i]]$data <- p$layers[[i]]$data[, intersect(colnames(p$layers[[i]]$data), vars), drop = FALSE]
      }
    }
  }
  return(p)
}

#' @param p A \code{patchwork} object.
#' @export
#' @rdname slim_data
#' @method slim_data patchwork
slim_data.patchwork <- function(p) {
  for (i in seq_along(p$patches$plots)) {
    p$patches$plots[[i]] <- slim_data(p$patches$plots[[i]])
  }
  p <- slim_data.ggplot(p)
  return(p)
}

#' @export
#' @method slim_data default
slim_data.default <- function(p) {
  return(p)
}


#' Get used vars in a ggplot object
#' @param p A \code{ggplot} object.
#' @param reverse If \code{TRUE} then will return unused vars.
#' @param verbose Whether to print messages.
#' @export
get_vars <- function(p, reverse, verbose = FALSE) {
  mappings <- c(
    as.character(p$mapping),
    unlist(lapply(p$layers, function(x) as.character(x$mapping))),
    unlist(lapply(p$layers, function(x) names(p$layers[[1]]$aes_params))),
    names(p$facet$params$facets), names(p$facet$params$rows), names(p$facet$params$cols)
  )
  vars <- unique(unlist(strsplit(gsub("[~\\[\\]\\\"\\(\\)]", " ", unique(mappings), perl = TRUE), " ")))
  vars_used <- intersect(unique(c(colnames(p$data), unlist(lapply(p$layers, function(x) colnames(x$data))))), vars)
  if (verbose) {
    message(
      "vars_used: ", paste0(vars_used, collapse = ","), "\n",
      "vars_notused: ", paste0(setdiff(names(p$data), vars), collapse = ",")
    )
  }
  return(vars_used)
}

#' Convert a color with arbitrary transparency to a fixed color
#'
#' This function takes a vector of colors and an alpha level and converts the colors
#' to fixed colors with the specified alpha level.
#'
#' @param colors Color vectors.
#' @param alpha Alpha level in [0,1]
#' @examples
#' colors <- c("red", "blue", "green")
#' adjcolors(colors, 0.5)
#' ggplot2::alpha(colors, 0.5)
#'
#' show_palettes(list(
#'   "raw" = colors,
#'   "adjcolors" = adjcolors(colors, 0.5),
#'   "ggplot2::alpha" = ggplot2::alpha(colors, 0.5)
#' ))
#'
#' @export
adjcolors <- function(colors, alpha) {
  color_df <- as.data.frame(col2rgb(colors) / 255)
  colors_out <- sapply(color_df, function(color) {
    color_rgb <- RGBA2RGB(list(color, alpha))
    return(rgb(color_rgb[1], color_rgb[2], color_rgb[3]))
  })
  return(colors_out)
}

#' Blend colors
#'
#' This function blends a list of colors using the specified blend mode.
#'
#' @param colors Color vectors.
#' @param mode Blend mode. One of "blend", "average", "screen", or "multiply".
#'
#' @seealso \code{\link{FeatureDimPlot}}
#'
#' @examples
#' blend <- c("red", "green", blendcolors(c("red", "green"), mode = "blend"))
#' average <- c("red", "green", blendcolors(c("red", "green"), mode = "average"))
#' screen <- c("red", "green", blendcolors(c("red", "green"), mode = "screen"))
#' multiply <- c("red", "green", blendcolors(c("red", "green"), mode = "multiply"))
#' show_palettes(list("blend" = blend, "average" = average, "screen" = screen, "multiply" = multiply))
#'
#' @export
blendcolors <- function(colors, mode = c("blend", "average", "screen", "multiply")) {
  mode <- match.arg(mode)
  colors <- colors[!is.na(colors)]
  if (length(colors) == 0) {
    return(NA)
  }
  if (length(colors) == 1) {
    return(colors)
  }
  rgb <- as.list(as.data.frame(col2rgb(colors) / 255))
  Clist <- lapply(rgb, function(x) {
    list(x, 1)
  })
  blend_color <- BlendRGBList(Clist, mode = mode)
  blend_color <- rgb(blend_color[1], blend_color[2], blend_color[3])
  return(blend_color)
}

RGBA2RGB <- function(RGBA, BackGround = c(1, 1, 1)) {
  A <- RGBA[[length(RGBA)]]
  RGB <- RGBA[[-length(RGBA)]] * A + BackGround * (1 - A)
  return(RGB)
}

Blend2Color <- function(C1, C2, mode = "blend") {
  c1 <- C1[[1]]
  c1a <- C1[[2]]
  c2 <- C2[[1]]
  c2a <- C2[[2]]
  A <- 1 - (1 - c1a) * (1 - c2a)
  if (A < 1.0e-6) {
    return(list(c(0, 0, 0), 1))
  }
  if (mode == "blend") {
    out <- (c1 * c1a + c2 * c2a * (1 - c1a)) / A
    A <- 1
  }
  if (mode == "average") {
    out <- (c1 + c2) / 2
    out[out > 1] <- 1
  }
  if (mode == "screen") {
    out <- 1 - (1 - c1) * (1 - c2)
  }
  if (mode == "multiply") {
    out <- c1 * c2
  }
  return(list(out, A))
}

BlendRGBList <- function(Clist, mode = "blend", RGB_BackGround = c(1, 1, 1)) {
  N <- length(Clist)
  ClistUse <- Clist
  while (N != 1) {
    temp <- ClistUse
    ClistUse <- list()
    for (C in temp[1:(length(temp) - 1)]) {
      c1 <- C[[1]]
      a1 <- C[[2]]
      c2 <- temp[[length(temp)]][[1]]
      a2 <- temp[[length(temp)]][[2]]
      ClistUse <- append(ClistUse, list(Blend2Color(C1 = list(c1, a1 * (1 - 1 / N)), C2 = list(c2, a2 * 1 / N), mode = mode)))
    }
    N <- length(ClistUse)
  }
  Result <- list(ClistUse[[1]][[1]], ClistUse[[1]][[2]])
  Result <- RGBA2RGB(Result, BackGround = RGB_BackGround)
  return(Result)
}

#' Visualize cell groups on a 2-dimensional reduction plot
#'
#' Plotting cell points on a reduced 2D plane and coloring according to the groups.
#'
#' @param srt A Seurat object.
#' @param group.by Name of one or more meta.data columns to group (color) cells by (for example, orig.ident).
#' @param reduction Which dimensionality reduction to use. If not specified, will use the reduction returned by \code{\link{DefaultReduction}}.
#' @param split.by Name of a column in meta.data column to split plot by.
#' @param palette Name of a color palette name collected in SCP. Default is "Paired".
#' @param palcolor Custom colors used to create a color palette.
#' @param bg_color Color value for background(NA) points.
#' @param pt.size Point size.
#' @param pt.alpha Point transparency.
#' @param cells.highlight A vector of cell names to highlight.
#' @param cols.highlight Color used to highlight the cells.
#' @param sizes.highlight Size of highlighted cell points.
#' @param alpha.highlight Transparency of highlighted cell points.
#' @param stroke.highlight Border width of highlighted cell points.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top").
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical")
#' @param combine Combine plots into a single \code{patchwork} object. If \code{FALSE}, return a list of ggplot objects.
#' @param nrow Number of rows in the combined plot.
#' @param ncol Number of columns in the combined plot.
#' @param byrow Logical value indicating if the plots should be arrange by row (default) or by column.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param show_na Whether to assign a color from the color palette to NA group. If \code{FALSE}, cell points with NA level will colored by \code{bg_color}.
#' @param show_stat Whether to show statistical information on the plot.
#' @param label Whether to label the cell groups.
#' @param label_insitu Whether to place the raw labels (group names) in the center of the cells with the corresponding group. Default is \code{FALSE}, which using numbers instead of raw labels.
#' @param label.size Size of labels.
#' @param label.fg Foreground color of label.
#' @param label.bg Background color of label.
#' @param label.bg.r Background ratio of label.
#' @param label_repel Logical value indicating whether the label is repel away from the center points.
#' @param label_repulsion Force of repulsion between overlapping text labels. Defaults to 20.
#' @param label_point_size Size of the center points.
#' @param label_point_color Color of the center points.
#' @param label_segment_color Color of the line segment for labels.
#' @param add_density Whether to add a density layer on the plot.
#' @param density_color Color of the density contours lines.
#' @param density_filled Whether to add filled contour bands instead of contour lines.
#' @param density_filled_palette Color palette used to fill contour bands.
#' @param density_filled_palcolor Custom colors used to fill contour bands.
#' @param lineages Lineages/pseudotime to add to the plot. If specified, curves will be fitted using \code{\link[stats]{loess}} method.
#' @param lineages_trim Trim the leading and the trailing data in the lineages.
#' @param lineages_span The parameter α which controls the degree of smoothing in \code{\link[stats]{loess}} method.
#' @param lineages_palette Color palette used for lineages.
#' @param lineages_palcolor Custom colors used for lineages.
#' @param lineages_arrow Set arrows of the lineages. See \code{\link[grid]{arrow}}.
#' @param lineages_linewidth Width of fitted curve lines for lineages.
#' @param lineages_line_bg Background color of curve lines for lineages.
#' @param lineages_line_bg_stroke Border width of curve lines background.
#' @param lineages_whiskers Whether to add whiskers for lineages.
#' @param lineages_whiskers_linewidth Width of whiskers for lineages.
#' @param lineages_whiskers_alpha Transparency of whiskers for lineages.
#' @param stat.by The name of a metadata column to stat.
#' @param stat_type Set stat types ("percent" or "count").
#' @param stat_plot_type Set the statistical plot type.
#' @param stat_plot_size Set the statistical plot size. Defaults to 0.1
#' @param stat_plot_palette Color palette used in statistical plot.
#' @param stat_palcolor Custom colors used in statistical plot
#' @param stat_plot_position Position adjustment in statistical plot.
#' @param stat_plot_alpha Transparency of the statistical plot.
#' @param stat_plot_label Whether to add labels in the statistical plot.
#' @param stat_plot_label_size Label size in the statistical plot.
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param edge_size Size of edges.
#' @param edge_alpha Transparency of edges.
#' @param edge_color Color of edges.
#' @param paga Specify the calculated paga results to add a PAGA graph layer to the plot.
#' @param paga_type PAGA plot type. "connectivities" or "connectivities_tree".
#' @param paga_node_size Size of the nodes in PAGA plot.
#' @param paga_edge_threshold Threshold of edge connectivities in PAGA plot.
#' @param paga_edge_size Size of edges in PAGA plot.
#' @param paga_edge_color Color of edges in PAGA plot.
#' @param paga_edge_alpha Transparency of edges in PAGA plot.
#' @param paga_show_transition Whether to show transitions between edges.
#' @param paga_transition_threshold Threshold of transition edges in PAGA plot.
#' @param paga_transition_size Size of transition edges in PAGA plot.
#' @param paga_transition_color Color of transition edges in PAGA plot.
#' @param paga_transition_alpha Transparency of transition edges in PAGA plot.
#' @param velocity Specify the calculated RNA velocity mode to add a velocity layer to the plot.
#' @param velocity_plot_type Set the velocity plot type.
#' @param velocity_n_neighbors Set the number of neighbors used in velocity plot.
#' @param velocity_density Set the density value used in velocity plot.
#' @param velocity_smooth Set the smooth value used in velocity plot.
#' @param velocity_scale Set the scale value used in velocity plot.
#' @param velocity_min_mass Set the min_mass value used in velocity plot.
#' @param velocity_cutoff_perc Set the cutoff_perc value used in velocity plot.
#' @param velocity_arrow_color Color of arrows in velocity plot.
#' @param velocity_arrow_angle Angle of arrows in velocity plot.
#' @param streamline_L Typical length of a streamline in x and y units
#' @param streamline_minL Minimum length of segments to show.
#' @param streamline_res Resolution parameter (higher numbers increases the resolution).
#' @param streamline_n Number of points to draw.
#' @param streamline_width Size of streamline.
#' @param streamline_alpha Transparency of streamline.
#' @param streamline_color Color of streamline.
#' @param streamline_palette Color palette used for streamline.
#' @param streamline_palcolor Custom colors used for streamline.
#' @param streamline_bg_color Background color of streamline.
#' @param streamline_bg_stroke Border width of streamline background.
#' @param hex Whether to chane the plot type from point to the hexagonal bin.
#' @param hex.count Whether show cell counts in each hexagonal bin.
#' @param hex.bins Number of hexagonal bins.
#' @param hex.binwidth Hexagonal bin width.
#' @param hex.linewidth Border width of hexagonal bins.
#' @param raster Convert points to raster format, default is NULL which automatically rasterizes if plotting more than 100,000 cells
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).
#' @param theme_use Theme used. Can be a character string or a theme function. For example, \code{"theme_blank"} or \code{ggplot2::theme_classic}.
#' @param aspect.ratio Aspect ratio of the panel.
#' @param title The text for the title.
#' @param subtitle The text for the subtitle for the plot which will be displayed below the title.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param force Whether to force drawing regardless of maximum levels in any cell group is greater than 100.
#' @param cells Subset cells to plot.
#' @param theme_args Other arguments passed to the \code{theme_use}.
#' @param seed Random seed set for reproducibility
#'
#' @seealso \code{\link{FeatureDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", theme_use = "theme_blank")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP") %>% panel_fix(height = 2, raster = TRUE, dpi = 30)
#'
#' # Highlight cells
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP",
#'   cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"]
#' )
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", split.by = "Phase", reduction = "UMAP",
#'   cells.highlight = TRUE, theme_use = "theme_blank", legend.position = "none"
#' )
#'
#' # Add group labels
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", label = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label.fg = "orange", label.bg = "red", label.size = 5
#' )
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE
#' )
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP",
#'   label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red"
#' )
#'
#' # Add various shape of marks
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_expand = unit(1, "mm"))
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_alpha = 0.3)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_linetype = 2)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_type = "ellipse")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_type = "rect")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_mark = TRUE, mark_type = "circle")
#'
#' # Add a density layer
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_density = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", add_density = TRUE, density_filled = TRUE)
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP",
#'   add_density = TRUE, density_filled = TRUE, density_filled_palette = "Blues",
#'   cells.highlight = TRUE
#' )
#'
#' # Add statistical charts
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", stat.by = "Phase")
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", stat.by = "Phase", stat_plot_type = "ring", stat_plot_label = TRUE, stat_plot_size = 0.15)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", stat.by = "Phase", stat_plot_type = "bar", stat_type = "count", stat_plot_position = "dodge")
#'
#' # Chane the plot type from point to the hexagonal bin
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", hex = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", hex = TRUE, hex.bins = 20)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", hex = TRUE, hex.count = FALSE)
#'
#' # Show neighbors graphs on the plot
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", graph = "Standardpca_SNN")
#' CellDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", graph = "Standardpca_SNN", edge_color = "grey80")
#'
#' # Show lineages on the plot based on the pseudotime
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", show_plot = FALSE)
#' FeatureDimPlot(pancreas_sub, features = paste0("Lineage", 1:3), reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3))
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_whiskers = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
#'
#' # Show PAGA results on the plot
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga, paga_type = "connectivities_tree")
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2,
#'   label = TRUE, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
#'   paga = pancreas_sub@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
#'   legend.position = "none", theme_use = "theme_blank"
#' )
#'
#' # Show RNA velocity results on the plot
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", mode = "stochastic", return_seurat = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga, paga_show_transition = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = NA, velocity = "stochastic")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "grid")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "grid", velocity_scale = 1.5)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "stream")
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2,
#'   label = TRUE, label_insitu = TRUE,
#'   velocity = "stochastic", velocity_plot_type = "stream", velocity_arrow_color = "yellow",
#'   velocity_density = 2, velocity_smooth = 1, streamline_n = 20, streamline_color = "black",
#'   legend.position = "none", theme_use = "theme_blank"
#' )
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom dplyr group_by "%>%" .data
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_manual scale_fill_manual guides guide_legend geom_hex geom_path theme_void annotation_custom scale_linewidth_continuous after_stat
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color new_scale_fill new_scale
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom ggforce geom_mark_hull geom_mark_ellipse geom_mark_circle geom_mark_rect
#' @importFrom grid arrow unit
#' @importFrom patchwork wrap_plots
#' @importFrom stats median loess aggregate
#' @importFrom utils askYesNo
#' @export
#'
