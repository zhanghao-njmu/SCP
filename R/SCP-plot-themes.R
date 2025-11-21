#' SCP theme
#'
#' The default theme for SCP plot function.
#'
#' @param aspect.ratio Aspect ratio of the panel.
#' @param base_size Base font size
#' @param ... Arguments passed to the \code{\link[ggplot2]{theme}}.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
#'   geom_point()
#' p + theme_scp()
#' @importFrom ggplot2 theme element_blank element_text element_rect margin
#' @export
theme_scp <- function(aspect.ratio = NULL, base_size = 12, ...) {
  text_size_scale <- base_size / 12
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(size = 12 * text_size_scale, color = "black"),
    plot.title = element_text(size = 14 * text_size_scale, colour = "black", vjust = 1),
    plot.subtitle = element_text(size = 13 * text_size_scale, hjust = 0, margin = margin(b = 3)),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_blank(),
    axis.title = element_text(size = 13 * text_size_scale, colour = "black"),
    axis.text = element_text(size = 12 * text_size_scale, colour = "black"),
    strip.text = element_text(size = 12.5 * text_size_scale, colour = "black", hjust = 0.5, margin = margin(3, 3, 3, 3)),
    strip.background = element_rect(fill = "transparent", linetype = 0),
    strip.switch.pad.grid = unit(-1, "pt"),
    strip.switch.pad.wrap = unit(-1, "pt"),
    strip.placement = "outside",
    legend.title = element_text(size = 12 * text_size_scale, colour = "black", hjust = 0),
    legend.text = element_text(size = 11 * text_size_scale, colour = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.key.size = unit(10, "pt"),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1)
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  return(out)
}

#' Blank theme
#'
#' This function creates a theme with all elements blank except for axis lines and labels.
#' It can optionally add coordinate axes in the plot.
#'
#' @param add_coord Whether to add coordinate arrows. Default is \code{TRUE}.
#' @param xlen_npc The length of the x-axis arrow in "npc".
#' @param ylen_npc The length of the y-axis arrow in "npc".
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param lab_size Label size.
#' @param ... Arguments passed to the \code{\link[ggplot2]{theme}}.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
#'   geom_point()
#' p + theme_blank()
#' p + theme_blank(xlab = "x-axis", ylab = "y-axis", lab_size = 16)
#' @importFrom ggplot2 theme element_blank margin annotation_custom coord_cartesian
#' @importFrom grid grobTree gList linesGrob textGrob arrow gpar
#' @export
theme_blank <- function(add_coord = TRUE, xlen_npc = 0.15, ylen_npc = 0.15, xlab = "", ylab = "", lab_size = 12, ...) {
  args1 <- list(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(10, "pt"),
    plot.margin = margin(lab_size + 2, lab_size + 2, lab_size + 2, lab_size + 2, unit = "points")
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  if (isTRUE(add_coord)) {
    g <- grobTree(gList(
      linesGrob(x = unit(c(0, xlen_npc), "npc"), y = unit(c(0, 0), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = xlab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3, hjust = 0, gp = gpar(fontsize = lab_size)),
      linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = ylab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3, hjust = 0, rot = 90, gp = gpar(fontsize = lab_size))
    ))
    return(list(
      list(annotation_custom(g)),
      list(theme_scp() + out),
      list(coord_cartesian(clip = "off"))
    ))
  } else {
    return(list(
      list(theme_scp() + out)
    ))
  }
}

#' Color palettes collected in SCP.
#'
#' @param x A vector of character/factor or numeric values. If missing, numeric values 1:n will be used as x.
#' @param n The number of colors to return for numeric values.
#' @param palette Palette name. All available palette names can be queried with \code{show_palettes()}.
#' @param palcolor Custom colors used to create a color palette.
#' @param type Type of \code{x}. Can be one of "auto", "discrete" or "continuous". The default is "auto", which automatically detects if \code{x} is a numeric value.
#' @param matched If \code{TRUE}, will return a color vector of the same length as \code{x}.
#' @param reverse Whether to invert the colors.
#' @param NA_keep Whether to keep the color assignment to NA in \code{x}.
#' @param NA_color Color assigned to NA if NA_keep is \code{TRUE}.
#'
#' @seealso \code{\link{show_palettes}} \code{\link{palette_list}}
#'
#' @examples
#' x <- c(1:3, NA, 3:5)
#' (pal1 <- palette_scp(x, palette = "Spectral"))
#' (pal2 <- palette_scp(x, palcolor = c("red", "white", "blue")))
#' (pal3 <- palette_scp(x, palette = "Spectral", n = 10))
#' (pal4 <- palette_scp(x, palette = "Spectral", n = 10, reverse = TRUE))
#' (pal5 <- palette_scp(x, palette = "Spectral", matched = TRUE))
#' (pal6 <- palette_scp(x, palette = "Spectral", matched = TRUE, NA_keep = TRUE))
#' (pal7 <- palette_scp(x, palette = "Paired", type = "discrete"))
#' show_palettes(list(pal1, pal2, pal3, pal4, pal5, pal6, pal7))
#'
#' all_palettes <- show_palettes(return_palettes = TRUE)
#' names(all_palettes)
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
#'
palette_scp <- function(x, n = 100, palette = "Paired", palcolor = NULL, type = "auto",
                        matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80") {
  palette_list <- SCP::palette_list
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }
  if (!palette %in% names(palette_list)) {
    stop("The palette name (", palette, ") is invalid! You can check the available palette names with 'show_palettes()'. Or pass palette colors via the 'palcolor' parameter.")
  }
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) {
    palcolor <- palette_list[[palette]]
  }
  if (!is.null(names(palcolor))) {
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }
  pal_n <- length(palcolor)

  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }

  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x))
    }
    n_x <- nlevels(x)
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- colorRampPalette(palcolor)(n_x)
    } else {
      color <- ifelse(rep(n_x, n_x) <= pal_n,
        palcolor[1:n_x],
        colorRampPalette(palcolor)(n_x)
      )
    }
    names(color) <- levels(x)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x) && all(!is.na(x))) {
      stop("'x' must be type of numeric when use continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, n))
    } else if (length(unique(na.omit(as.numeric(x)))) == 1) {
      values <- as.factor(rep(unique(na.omit(as.numeric(x))), n))
    } else {
      if (isTRUE(matched)) {
        values <- cut(x, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
      } else {
        values <- cut(1:100, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
      }
    }

    n_x <- nlevels(values)
    color <- ifelse(rep(n_x, n_x) <= pal_n,
      palcolor[1:n_x],
      colorRampPalette(palcolor)(n_x)
    )
    names(color) <- levels(values)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      if (all(is.na(x))) {
        color <- NA_color
      } else if (length(unique(na.omit(x))) == 1) {
        color <- color[as.character(unique(na.omit(x)))]
        color[is.na(color)] <- NA_color
      } else {
        color <- color[as.character(values)]
        color[is.na(color)] <- NA_color
      }
    }
  }

  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}

#' Show the color palettes
#'
#' This function displays color palettes using ggplot2.
#'
#' @param palettes A list of color palettes. If `NULL`, uses default palettes.
#' @param type A character vector specifying the type of palettes to include. Default is "discrete".
#' @param index A numeric vector specifying the indices of the palettes to include. Default is `NULL`.
#' @param palette_names A character vector specifying the names of the SCP palettes to include. Default is `NULL`.
#' @param return_names A logical value indicating whether to return the names of the selected palettes. Default is `TRUE`.
#' @param return_palettes A logical value indicating whether to return the colors of selected palettes. Default is `FALSE`.
#'
#' @seealso \code{\link{palette_scp}} \code{\link{palette_list}}
#'
#' @examples
#' show_palettes(palettes = list(c("red", "blue", "green"), c("yellow", "purple", "orange")))
#' all_palettes <- show_palettes(return_palettes = TRUE)
#' names(all_palettes)
#' all_palettes[["simspec"]]
#' show_palettes(index = 1:10)
#' show_palettes(type = "discrete", index = 1:10)
#' show_palettes(type = "continuous", index = 1:10)
#' show_palettes(palette_names = c("Paired", "nejm", "simspec", "Spectral", "jet"), return_palettes = TRUE)
#'
#' @importFrom ggplot2 ggplot geom_col scale_fill_manual scale_x_continuous element_blank
#' @export
show_palettes <- function(palettes = NULL, type = c("discrete", "continuous"), index = NULL, palette_names = NULL, return_names = TRUE, return_palettes = FALSE) {
  palette_list <- SCP::palette_list
  if (!is.null(palettes)) {
    palette_list <- palettes
  } else {
    palette_list <- palette_list[unlist(lapply(palette_list, function(x) isTRUE(attr(x, "type") %in% type)))]
  }
  index <- index[index %in% seq_along(palette_list)]
  if (!is.null(index)) {
    palette_list <- palette_list[index]
  }
  if (is.null(names(palette_list))) {
    names(palette_list) <- seq_along(palette_list)
  }
  if (is.null(palette_names)) {
    palette_names <- palette_names %||% names(palette_list)
  }
  if (any(!palette_names %in% names(palette_list))) {
    stop(paste("Can not find the palettes: ", paste0(palette_names[!palette_names %in% names(palette_list)], collapse = ",")))
  }
  palette_list <- palette_list[palette_names]

  df <- data.frame(palette = rep(names(palette_list), sapply(palette_list, length)), color = unlist(palette_list))
  df[["palette"]] <- factor(df[["palette"]], levels = rev(unique(df[["palette"]])))
  df[["color_order"]] <- factor(seq_len(nrow(df)), levels = seq_len(nrow(df)))
  df[["proportion"]] <- as.numeric(1 / table(df$palette)[df$palette])

  p <- ggplot(data = df, aes(y = .data[["palette"]], x = .data[["proportion"]], fill = .data[["color_order"]])) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = df[["color"]]) +
    scale_x_continuous(expand = c(0, 0), trans = "reverse") +
    theme_scp(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank()
    )
  print(p)

  if (isTRUE(return_palettes)) {
    return(palette_list)
  }
  if (isTRUE(return_names)) {
    return(palette_names)
  }
}

#' Set the panel width/height of a plot object to a fixed value.
#'
#' The ggplot object, when stored, can only specify the height and width of the entire plot, not the panel.
#' The latter is obviously more important to control the final result of a plot.
#' This function can set the panel width/height of plot to a fixed value and rasterize it.
#'
#' @param x A ggplot object, a grob object, or a combined plot made by patchwork or cowplot package.
#' @param panel_index Specify the panel to be fixed. If NULL, will fix all panels.
#' @param respect If a logical, this indicates whether row heights and column widths should respect each other.
#' @param width The desired width of the fixed panels.
#' @param height The desired height of the fixed panels.
#' @param margin The margin to add around each panel, in inches. The default is 1 inch.
#' @param padding The padding to add around each panel, in inches. The default is 0 inches.
#' @param units The units in which \code{height}, \code{width} and \code{margin} are given. Can be \code{mm}, \code{cm}, \code{in}, etc. See \code{\link[grid]{unit}}.
#' @param raster Whether to rasterize the panel.
#' @param dpi Plot resolution.
#' @param BPPARAM An \code{\link[BiocParallel]{BiocParallelParam}} instance determining the parallel back-end to be used during building the object made by patchwork package.
#' @param return_grob If \code{TRUE} then return a grob object instead of a wrapped \code{patchwork} object.
#' @param save NULL or the file name used to save the plot.
#' @param bg_color Plot background color.
#' @param verbose Whether to print messages.
#' @param ... Unused.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
#'   geom_point() +
#'   facet_wrap(~gear, nrow = 2)
#' # fix the size of panel
#' panel_fix(p, width = 5, height = 3, units = "cm")
#' # rasterize the panel
#' panel_fix(p, width = 5, height = 3, units = "cm", raster = TRUE, dpi = 30)
#'
#' # panel_fix will build and render the plot when the input is a ggplot object.
#' # so after panel_fix, the size of the object will be changed.
#' object.size(p)
#' object.size(panel_fix(p, width = 5, height = 3, units = "cm"))
#'
#' ## save the plot with appropriate size
#' # p_fix <- panel_fix(p, width = 5, height = 3, units = "cm")
#' # plot_size <- attr(p_fix, "size")
#' # ggsave(
#' #   filename = "p_fix.png", plot = p_fix,
#' #   units = plot_size$units, width = plot_size$width, height = plot_size$height
#' # )
#' ## or save the plot directly
#' # p_fix <- panel_fix(p, width = 5, height = 3, units = "cm", save = "p_fix.png")
#'
#' # fix the panel of the plot combined by patchwork
#' data("pancreas_sub")
#' p1 <- CellDimPlot(pancreas_sub, "Phase", aspect.ratio = 1) # ggplot object
#' p2 <- FeatureDimPlot(pancreas_sub, "Ins1", aspect.ratio = 0.5) # ggplot object
#' p <- p1 / p2 # plot is combined by patchwork
#' # fix the panel size for each plot, the width will be calculated automatically based on aspect.ratio
#' panel_fix(p, height = 1)
#'
#' # fix the panel of the plot combined by plot_grid
#' if (requireNamespace("cowplot", quietly = TRUE)) {
#'   p1 <- CellDimPlot(pancreas_sub, c("Phase", "SubCellType"), label = TRUE) # plot is combined by patchwork
#'   p2 <- FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg"), label = TRUE) # plot is combined by patchwork
#'   p <- cowplot::plot_grid(p1, p2, nrow = 2) # plot is combined by plot_grid
#'   # fix the size of panel for each plot
#'   panel_fix(p, height = 1)
#'   # rasterize the panel while keeping all labels and text in vector format
#'   panel_fix(p, height = 1, raster = TRUE, dpi = 30)
#' }
#'
#' # fix the panel of the heatmap
#' ht1 <- GroupHeatmap(pancreas_sub,
#'   features = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   group.by = c("CellType", "SubCellType"),
#'   show_row_names = TRUE
#' )
#' # the size of the heatmap is not fixed and can be resized by zooming the viewport.
#' ht1$plot
#' # fix the size of the heatmap according the current viewport
#' panel_fix(ht1$plot)
#' # rasterize the heatmap body
#' panel_fix(ht1$plot, raster = TRUE, dpi = 30)
#' # fix the size of overall heatmap including annotation and legend
#' panel_fix(ht1$plot, height = 4, width = 6)
#'
#' ht2 <- GroupHeatmap(pancreas_sub,
#'   features = pancreas_sub[["RNA"]]@var.features,
#'   group.by = "SubCellType",
#'   n_split = 5, nlabel = 20,
#'   db = "GO_BP", species = "Mus_musculus", anno_terms = TRUE,
#'   height = 4, width = 1 # Heatmap body size for two groups
#' )
#' # the size of the heatmap is already fixed
#' ht2$plot
#' # when no height/width is specified, panel_fix does not change the size of the heatmap.
#' panel_fix(ht2$plot)
#' # rasterize the heatmap body
#' panel_fix(ht2$plot, raster = TRUE, dpi = 30)
#' # however, gene labels on the left and enrichment annotations on the right cannot be adjusted
#' panel_fix(ht2$plot, height = 5, width = 10)
#'
#' @importFrom ggplot2 ggsave
#' @importFrom grid grob unit convertWidth convertHeight convertUnit is.unit unitType
#' @export
#'
panel_fix <- function(x = NULL, panel_index = NULL, respect = NULL,
                      width = NULL, height = NULL, margin = 1, padding = 0, units = "in",
                      raster = FALSE, dpi = 300, BPPARAM = BiocParallel::SerialParam(),
                      return_grob = FALSE, bg_color = "white", save = NULL, verbose = FALSE, ...) {
  if (!inherits(x, "gtable")) {
    tryCatch(
      {
        gtable <- as_gtable(x, BPPARAM = BPPARAM)
      },
      error = function(error) {
        stop(error, "\nCannot convert the x to a gtable object.")
      }
    )
  } else {
    gtable <- x
  }
  args <- as.list(match.call())[-1]
  depth <- args[["depth"]]
  if (is.null(depth)) {
    depth <- 1
  }

  if (is.null(panel_index)) {
    non_zero <- grep(pattern = "zeroGrob", vapply(gtable$grobs, as.character, character(1)), invert = TRUE)
    panel_index <- grep("panel|full", gtable[["layout"]][["name"]])
    panel_index <- intersect(panel_index, non_zero)
  }
  if (length(panel_index) == 0 && length(gtable$grobs) == 1) {
    panel_index <- 1
  }
  add_margin <- TRUE
  for (i in panel_index) {
    geom_index <- grep("GeomDrawGrob", names(gtable$grobs[[i]][["children"]]))
    if (length(geom_index) > 0) {
      if (isTRUE(verbose)) {
        message("panel ", i, " is detected as generated by plot_grid.")
      }
      for (j in geom_index) {
        subgrob <- gtable$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]]
        # print(subgrob$grobs[[1]][["children"]])
        if (length(subgrob$grobs[[1]][["children"]]) > 0 && all(sapply(subgrob$grobs[[1]][["children"]], function(x) inherits(x, "recordedGrob")))) {
          subgrob <- panel_fix_overall(subgrob$grobs[[1]][["children"]], width = width, height = height, margin = padding, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
        } else {
          subgrob <- panel_fix(subgrob, width = width, height = height, margin = padding, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
        }
        gtable$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]] <- subgrob
        # print(paste0("plot_width:",plot_width," plot_height:",plot_height))
      }
      sum_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$width)
      sum_height <- convertHeight(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$height)
      gtable <- panel_fix_overall(gtable, panel_index = i, width = sum_width, height = sum_height, margin = ifelse(depth == 1, margin, 0), units = units, raster = FALSE, return_grob = TRUE)
    } else if (gtable$grobs[[i]]$name == "layout" || inherits(x, "patchwork")) {
      if (isTRUE(verbose)) {
        message("panel ", i, " is detected as generated by patchwork.")
      }
      # if (i == panel_index[1] && length(panel_index) > 1 && isTRUE(verbose)) {
      #   message("More than 2 panels detected. panel_fix may not work as expected.")
      # }
      subgrob <- gtable$grobs[[i]]
      if (length(subgrob[["children"]]) > 0 && all(sapply(subgrob[["children"]], function(x) inherits(x, "recordedGrob")))) {
        subgrob <- panel_fix_overall(subgrob[["children"]], width = width, height = height, margin = 0, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
      } else {
        subgrob <- panel_fix(subgrob, width = width, height = height, margin = 0, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
      }
      gtable$grobs[[i]] <- subgrob
      layout <- gtable$layout
      layout[["rowranges"]] <- lapply(seq_len(nrow(layout)), function(n) layout$t[n]:layout$b[n])
      layout[["colranges"]] <- lapply(seq_len(nrow(layout)), function(n) layout$l[n]:layout$r[n])
      p_row <- c(layout$t[i], layout$b[i])
      p_col <- c(layout$l[i], layout$r[i])
      background_index <- grep("background", layout$name)
      background_index <- background_index[order(layout$z[background_index], decreasing = TRUE)]
      for (bgi in background_index) {
        if (all(p_row %in% layout[["rowranges"]][[bgi]]) && all(p_col %in% layout[["colranges"]][[bgi]])) {
          p_background_index <- bgi
          break
        }
      }
      gtable <- gtable_add_rows(gtable, heights = unit(padding, units), pos = layout$t[p_background_index] - 1)
      gtable <- gtable_add_rows(gtable, heights = unit(padding, units), pos = layout$b[p_background_index])
      gtable <- gtable_add_cols(gtable, widths = unit(padding, units), pos = layout$l[p_background_index] - 1)
      gtable <- gtable_add_cols(gtable, widths = unit(padding, units), pos = layout$r[p_background_index])
      sum_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE)
      sum_height <- convertHeight(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE)

      gtable <- panel_fix_overall(gtable, panel_index = i, width = sum_width, height = sum_height, margin = ifelse(depth == 1 & add_margin, margin, 0), units = units, raster = FALSE, respect = TRUE, return_grob = TRUE)
      if (depth == 1 & add_margin) {
        add_margin <- FALSE
      }
    } else {
      # print("fix the gtable")
      gtable <- panel_fix_overall(gtable, panel_index = i, width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
    }
  }

  if (!is.null(respect)) {
    gtable$respect <- respect
  }

  if (isTRUE(return_grob)) {
    return(gtable)
  } else {
    p <- wrap_plots(gtable) + theme(plot.background = element_rect(fill = bg_color, color = bg_color))
    if (units != "null") {
      plot_width <- convertWidth(sum(gtable[["widths"]]), unitTo = units, valueOnly = TRUE)
      plot_height <- convertHeight(sum(gtable[["heights"]]), unitTo = units, valueOnly = TRUE)
      attr(p, "size") <- list(width = plot_width, height = plot_height, units = units)
    }

    if (!is.null(save) && is.character(save) && nchar(save) > 0) {
      if (units == "null") {
        stop("units can not be 'null' if want to save the plot.")
      }
      filename <- normalizePath(save)
      if (isTRUE(verbose)) {
        message("Save the plot to the file: ", filename)
      }
      if (!dir.exists(dirname(filename))) {
        dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      }
      ggsave(
        plot = p, filename = filename, width = plot_width, height = plot_height, units = units, dpi = dpi, limitsize = FALSE
      )
    }
    return(p)
  }
}

#' @rdname panel_fix
#' @importFrom ggplot2 ggsave zeroGrob
#' @importFrom gtable gtable_add_padding
#' @importFrom grid grob unit unitType convertWidth convertHeight convertUnit viewport grid.draw rasterGrob grobTree addGrob
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @export
panel_fix_overall <- function(x, panel_index = NULL, respect = NULL,
                              width = NULL, height = NULL, margin = 1, units = "in",
                              raster = FALSE, dpi = 300, BPPARAM = BiocParallel::SerialParam(),
                              return_grob = FALSE, bg_color = "white", save = NULL, verbose = TRUE) {
  if (!inherits(x, "gtable")) {
    if (inherits(x, "gTree")) {
      x <- x[["children"]]
    }
    tryCatch(
      {
        gtable <- as_gtable(x, BPPARAM = BPPARAM)
      },
      error = function(error) {
        stop(error, "\nCannot convert the x to a gtable object")
      }
    )
  } else {
    gtable <- x
  }

  if (is.null(panel_index)) {
    non_zero <- grep(pattern = "zeroGrob", vapply(gtable$grobs, as.character, character(1)), invert = TRUE)
    panel_index <- grep("panel|full", gtable[["layout"]][["name"]])
    panel_index <- intersect(panel_index, non_zero)
  }
  if (length(panel_index) == 0 && length(gtable$grobs) == 1) {
    panel_index <- 1
  }
  if (!length(width) %in% c(0, 1, length(panel_index)) || !length(height) %in% c(0, 1, length(panel_index))) {
    stop("The length of 'width' and 'height' must be 1 or the length of panels.")
  }

  if (inherits(x, "gList")) {
    panel_index <- 1
    panel_index_h <- panel_index_w <- list(1)
    w_comp <- h_comp <- list(unit(1, "null"))
    w <- h <- list(unit(1, "null"))
  } else if (length(panel_index) > 0) {
    panel_index_w <- panel_index_h <- list()
    w_comp <- h_comp <- list()
    w <- h <- list()
    for (i in seq_along(panel_index)) {
      index <- panel_index[i]
      panel_index_h[[i]] <- sort(unique(c(gtable[["layout"]][["t"]][index], gtable[["layout"]][["b"]][index])))
      panel_index_w[[i]] <- sort(unique(c(gtable[["layout"]][["l"]][index], gtable[["layout"]][["r"]][index])))
      w_comp[[i]] <- gtable[["widths"]][seq(min(panel_index_w[[i]]), max(panel_index_w[[i]]))]
      h_comp[[i]] <- gtable[["heights"]][seq(min(panel_index_h[[i]]), max(panel_index_h[[i]]))]

      if (length(w_comp[[i]]) == 1) {
        w[[i]] <- w_comp[[i]]
      } else if (length(w_comp[[i]]) > 1 && any(unitType(w_comp[[i]]) == "null")) {
        w[[i]] <- unit(1, units = "null")
      } else {
        w[[i]] <- sum(w_comp[[i]])
      }
      if (length(h_comp[[i]]) == 1) {
        h[[i]] <- h_comp[[i]]
      } else if (length(h_comp[[i]]) > 1 && any(unitType(h_comp[[i]]) == "null")) {
        h[[i]] <- unit(1, units = "null")
      } else {
        h[[i]] <- sum(h_comp[[i]])
      }
    }
  } else {
    stop("No panel detected.")
  }

  if (units != "null") {
    raw_w <- sapply(w, function(x) convertWidth(x, unitTo = units, valueOnly = TRUE))
    raw_h <- sapply(h, function(x) convertHeight(x, unitTo = units, valueOnly = TRUE))
    for (i in seq_along(w)) {
      if (unitType(w[[i]]) == "null" || convertUnit(w[[i]], unitTo = "cm", valueOnly = TRUE) < 1e-10) {
        raw_w[i] <- 0
      }
    }
    for (i in seq_along(h)) {
      if (unitType(h[[i]]) == "null" || convertUnit(h[[i]], unitTo = "cm", valueOnly = TRUE) < 1e-10) {
        raw_h[i] <- 0
      }
    }
    if (isTRUE(gtable$respect)) {
      raw_aspect <- sapply(h, as.vector) / sapply(w, as.vector)
    } else {
      if (all(raw_w != 0) && all(raw_h != 0)) {
        raw_aspect <- raw_h / raw_w
      } else {
        raw_aspect <- convertHeight(unit(1, "npc"), "cm", valueOnly = TRUE) / convertWidth(unit(1, "npc"), "cm", valueOnly = TRUE)
      }
    }

    if (is.null(width) && is.null(height)) {
      width <- raw_w
      height <- raw_h
      if (all(width == 0) && all(height == 0)) {
        width <- convertWidth(unit(1, "npc"), units, valueOnly = TRUE)
        height <- convertHeight(unit(1, "npc"), units, valueOnly = TRUE)
        if (isTRUE(gtable$respect)) {
          if (raw_aspect <= 1) {
            height <- width * raw_aspect
          } else {
            width <- height / raw_aspect
          }
        }
      }
    }

    for (i in seq_along(raw_aspect)) {
      if (is.finite(raw_aspect[i]) && raw_aspect[i] != 0) {
        if (is.null(width[i]) || is.na(width[i]) || width[i] == 0) {
          width[i] <- height[i] / raw_aspect[i]
        }
        if (is.null(height[i]) || is.na(height[i]) || height[i] == 0) {
          height[i] <- width[i] * raw_aspect[i]
        }
      }
    }

    for (i in seq_along(width)) {
      if (inherits(width[i], "unit")) {
        width[i] <- convertWidth(width[i], unitTo = units, valueOnly = TRUE)
      }
    }
    for (i in seq_along(height)) {
      if (inherits(height[i], "unit")) {
        height[i] <- convertHeight(height[i], unitTo = units, valueOnly = TRUE)
      }
    }
  }

  if (length(width) == 1) {
    width <- rep(width, length(panel_index))
  }
  if (length(height) == 1) {
    height <- rep(height, length(panel_index))
  }
  for (i in seq_along(panel_index)) {
    if (!is.null(width)) {
      width_unit <- width[i] / length(w_comp[[i]])
      gtable[["widths"]][seq(min(panel_index_w[[i]]), max(panel_index_w[[i]]))] <- rep(unit(width_unit, units = units), length(w_comp[[i]]))
    }
    if (!is.null(height)) {
      height_unit <- height[i] / length(h_comp[[i]])
      gtable[["heights"]][seq(min(panel_index_h[[i]]), max(panel_index_h[[i]]))] <- rep(unit(height_unit, units = units), length(h_comp[[i]]))
    }
  }
  gtable <- gtable_add_padding(gtable, padding = unit(margin, units = units))

  # print(paste0("width:", width))
  # print(paste0("height:", height))

  if (isTRUE(raster)) {
    check_R(c("png", "ragg"))
    for (i in seq_along(panel_index)) {
      index <- panel_index[i]
      g <- g_new <- gtable$grobs[[index]]
      vp <- g$vp
      childrenOrder <- g$childrenOrder
      if (is.null(g$vp)) {
        g$vp <- viewport()
      }
      # child_list <- NULL
      for (j in seq_along(g[["children"]])) {
        child <- g[["children"]][[j]]
        child_nm <- names(g[["children"]])[j]
        if (!is.null(child$vp) ||
          any(grepl("(text)|(label)", child_nm)) ||
          any(grepl("(text)|(segments)|(legend)", class(child$list[[1]])))) {
          zero <- zeroGrob()
          zero$name <- g[["children"]][[j]]$name
          g[["children"]][[j]] <- zero
          # child_list[[child_nm]] <- child
          # print(j)
        } else if (inherits(child$list[[1]], "grob") || is.null(child$list[[1]])) {
          g_new[["children"]][[j]] <- zeroGrob()
          # print(j)
        }
      }
      temp <- tempfile(fileext = "png")
      ragg::agg_png(temp, width = width[i], height = height[i], bg = "transparent", res = dpi, units = units)
      grid.draw(g)
      dev.off()
      g_ras <- rasterGrob(png::readPNG(temp, native = TRUE))
      unlink(temp)
      # g <- do.call(grobTree, c(list(g_ras), child_list))
      g <- addGrob(g_new, g_ras)
      g$vp <- vp
      g$childrenOrder <- c(g_ras$name, childrenOrder)
      gtable$grobs[[index]] <- g

      # grid.draw(gtable)
    }
  }

  if (!is.null(respect)) {
    gtable$respect <- respect
  }

  if (isTRUE(return_grob)) {
    return(gtable)
  } else {
    p <- wrap_plots(gtable) + theme(plot.background = element_rect(fill = bg_color, color = bg_color))
    if (units != "null") {
      plot_width <- convertWidth(sum(gtable[["widths"]]), unitTo = units, valueOnly = TRUE)
      plot_height <- convertHeight(sum(gtable[["heights"]]), unitTo = units, valueOnly = TRUE)
      attr(p, "size") <- list(width = plot_width, height = plot_height, units = units)
    }

    if (!is.null(save) && is.character(save) && nchar(save) > 0) {
      if (units == "null") {
        stop("units can not be 'null' if want to save the plot.")
      }
      filename <- normalizePath(save)
      if (isTRUE(verbose)) {
        message("Save the plot to the file: ", filename)
      }
      if (!dir.exists(dirname(filename))) {
        dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      }
      ggsave(
        plot = p, filename = filename, width = plot_width, height = plot_height, units = units,
        dpi = dpi, limitsize = FALSE
      )
    }
    return(p)
  }
}

#' Drop all data in the plot (only one observation is kept)
#'
#' @param p A \code{ggplot} object or a \code{patchwork} object.
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars, aes(x = mpg, y = wt, colour = cyl)) +
#'   geom_point() +
#'   scale_x_continuous(limits = c(10, 30)) +
#'   scale_y_continuous(limits = c(1, 6)) +
#'   theme_scp()
#' object.size(p)
#'
#' p_drop <- drop_data(p)
#' object.size(p_drop)
#'
#' p / p_drop
#'
#' @export
