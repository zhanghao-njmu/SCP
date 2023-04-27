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
#'
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
    legend.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1),
    complete = TRUE
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
#' @param add_coord Whether to add coordinate arrows. Default is \code{TRUE}.
#' @param xlen_npc The length of the x-axis arrow in "npc".
#' @param ylen_npc The length of the y-axis arrow in "npc".
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param lab_size Label size.
#' @param ... Arguments passed to the \code{\link[ggplot2]{theme}}.
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
  if (isTRUE(add_coord) && isTRUE(xlab != "")) {
    x_space <- lab_size + 2
  } else {
    x_space <- 0
  }
  if (isTRUE(add_coord) && isTRUE(ylab != "")) {
    y_space <- lab_size + 2
  } else {
    y_space <- 0
  }
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
    plot.margin = margin(0, 0, x_space, y_space, unit = "points"),
    complete = FALSE
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
#' @seealso
#' \code{\link{show_palettes}}
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
#' \dontrun{
#' if (interactive()) {
#'   check_R(c("stringr", "RColorBrewer", "ggsci", "Redmonder", "rcartocolor", "nord", "viridis", "pals", "oompaBase", "dichromat", "jcolors"))
#'   library(stringr)
#'   library(RColorBrewer)
#'   library(ggsci)
#'   library(Redmonder)
#'   library(rcartocolor)
#'   library(nord)
#'   library(viridis)
#'   library(pals)
#'   library(dichromat)
#'   library(jcolors)
#'   brewer.pal.info <- RColorBrewer::brewer.pal.info
#'   ggsci_db <- ggsci:::ggsci_db
#'   redmonder.pal.info <- Redmonder::redmonder.pal.info
#'   metacartocolors <- rcartocolor::metacartocolors
#'   rownames(metacartocolors) <- metacartocolors$Name
#'   nord_palettes <- nord::nord_palettes
#'   viridis_names <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
#'   viridis_palettes <- lapply(setNames(viridis_names, viridis_names), function(x) viridis::viridis(100, option = x))
#'   ocean_names <- names(pals:::syspals)[grep("ocean", names(pals:::syspals))]
#'   ocean_palettes <- pals:::syspals[ocean_names]
#'   dichromat_palettes <- dichromat::colorschemes
#'   jcolors_names <- paste0("jcolors-", c("default", "pal2", "pal3", "pal4", "pal5", "pal6", "pal7", "pal8", "pal9", "pal10", "pal11", "pal12", "rainbow"))
#'   custom_names <- c("jet", "simspec")
#'   custom_palettes <- list(
#'     oompaBase::jetColors(N = 100),
#'     c("#c22b86", "#f769a1", "#fcc5c1", "#253777", "#1d92c0", "#9ec9e1", "#015b33", "#42aa5e", "#d9f0a2", "#E66F00", "#f18c28", "#FFBB61")
#'   )
#'   names(custom_palettes) <- custom_names
#'
#'   palette_list <- list()
#'   all_colors <- c(
#'     rownames(brewer.pal.info), names(ggsci_db), rownames(redmonder.pal.info),
#'     rownames(metacartocolors), names(nord_palettes), names(viridis_palettes),
#'     ocean_names, names(dichromat_palettes), jcolors_names,
#'     custom_names
#'   )
#'   for (pal in all_colors) {
#'     if (!pal %in% all_colors) {
#'       stop(paste0("Invalid pal Must be one of ", paste0(all_colors, collapse = ",")))
#'     }
#'     if (pal %in% rownames(brewer.pal.info)) {
#'       pal_n <- brewer.pal.info[pal, "maxcolors"]
#'       pal_category <- brewer.pal.info[pal, "category"]
#'       if (pal_category == "div") {
#'         palcolor <- rev(brewer.pal(name = pal, n = pal_n))
#'       } else {
#'         if (pal == "Paired") {
#'           palcolor <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
#'         } else {
#'           palcolor <- brewer.pal(name = pal, n = pal_n)
#'         }
#'       }
#'       if (pal_category == "qual") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% names(ggsci_db)) {
#'       if (pal %in% c("d3", "uchicago", "material")) {
#'         for (subpal in names(ggsci_db[[pal]])) {
#'           palcolor <- ggsci_db[[pal]][[subpal]]
#'           if (pal == "material") {
#'             attr(palcolor, "type") <- "continuous"
#'           } else {
#'             attr(palcolor, "type") <- "discrete"
#'           }
#'           palette_list[[paste0(pal, "-", subpal)]] <- palcolor
#'         }
#'         next
#'       } else {
#'         palcolor <- ggsci_db[[pal]][[1]]
#'         if (pal == "gsea") {
#'           attr(palcolor, "type") <- "continuous"
#'         } else {
#'           attr(palcolor, "type") <- "discrete"
#'         }
#'       }
#'     } else if (pal %in% rownames(redmonder.pal.info)) {
#'       pal_n <- redmonder.pal.info[pal, "maxcolors"]
#'       pal_category <- redmonder.pal.info[pal, "category"]
#'       if (pal_category == "div") {
#'         palcolor <- rev(redmonder.pal(name = pal, n = pal_n))
#'       } else {
#'         palcolor <- redmonder.pal(name = pal, n = pal_n)
#'       }
#'       if (pal_category == "qual") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% rownames(metacartocolors)) {
#'       pal_n <- metacartocolors[pal, "Max_n"]
#'       palcolor <- carto_pal(name = pal, n = pal_n)
#'       if (pal_category == "qualitative") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% names(nord_palettes)) {
#'       palcolor <- nord_palettes[[pal]]
#'       attr(palcolor, "type") <- "discrete"
#'     } else if (pal %in% names(viridis_palettes)) {
#'       palcolor <- viridis_palettes[[pal]]
#'       attr(palcolor, "type") <- "continuous"
#'     } else if (pal %in% names(ocean_palettes)) {
#'       palcolor <- ocean_palettes[[pal]]
#'       attr(palcolor, "type") <- "continuous"
#'     } else if (pal %in% names(dichromat_palettes)) {
#'       palcolor <- dichromat_palettes[[pal]]
#'       if (pal %in% c("Categorical.12", "SteppedSequential.5")) {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% jcolors_names) {
#'       palcolor <- jcolors(palette = gsub("jcolors-", "", pal))
#'       if (pal %in% paste0("jcolors-", c("pal10", "pal11", "pal12", "rainbow"))) {
#'         attr(palcolor, "type") <- "continuous"
#'       } else {
#'         attr(palcolor, "type") <- "discrete"
#'       }
#'     } else if (pal %in% custom_names) {
#'       palcolor <- custom_palettes[[pal]]
#'       if (pal %in% c("jet")) {
#'         attr(palcolor, "type") <- "continuous"
#'       } else {
#'         attr(palcolor, "type") <- "discrete"
#'       }
#'     }
#'     palette_list[[pal]] <- palcolor
#'   }
#'   # usethis::use_data(palette_list, internal = TRUE)
#' }
#' }
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
#'
palette_scp <- function(x, n = 100, palette = "Paired", palcolor = NULL, type = "auto",
                        matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80") {
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
#' @param palettes A list of custom color palettes to be shown.
#' @param type Specifies the type of color palettes collected in SCP.
#' @param index The index of the palette in SCP palette list. If \code{NULL}, show all palettes collected in SCP.
#' @param palette_names Specifies name of the color palettes collected in SCP.
#' @param return_names Whether to return the palette names.
#' @param return_palettes Whether to return the palettes.
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
#' @importFrom rlang %||%
#' @importFrom graphics par text rect
#' @export
show_palettes <- function(palettes = NULL, type = c("discrete", "continuous"), index = NULL, palette_names = NULL, return_names = TRUE, return_palettes = FALSE) {
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

  par(mar = c(0, 0, 0, 0) + 0.1)
  plot(0, 0,
    type = "n", axes = FALSE, bty = "n", xlab = "", ylab = "",
    xlim = c(0, 1), ylim = c(-length(palette_list) - 1, -1)
  )
  for (i in seq_len(length(palette_list))) {
    colors_len <- length(palette_list[[i]])
    breaks <- seq(from = 0, to = 1, length = colors_len + 1)
    text(0, -i, palette_names[i], pos = 4)
    rect(
      xleft = breaks[1:colors_len], xright = breaks[1:colors_len + 1],
      ytop = -0.15 - i, ybottom = -0.8 - i,
      col = palette_list[[i]], border = NA
    )
  }
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
#' @param width The width of the panel.
#' @param height The height of the panel.
#' @param margin Margins around the plot.
#' @param units The units in which \code{height}, \code{width} and \code{margin} are given. Can be \code{mm}, \code{cm}, \code{in}, etc. See \code{\link[grid]{unit}}.
#' @param raster Whether to rasterize the panel.
#' @param dpi Plot resolution.
#' @param BPPARAM An \code{\link[BiocParallel]{BiocParallelParam}} instance determining the parallel back-end to be used during building the object made by patchwork package.
#' @param return_grob If \code{TRUE} then return a grob object instead of a wrapped \code{patchwork} object.
#' @param save NULL or the file name used to save the plot.
#' @param bg_color Plot background color.
#' @param verbose Whether to print messages.
#' @param device
#' @param ...
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
#' library(cowplot)
#' p1 <- CellDimPlot(pancreas_sub, c("Phase", "SubCellType"), label = TRUE) # plot is combined by patchwork
#' p2 <- FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg"), label = TRUE) # plot is combined by patchwork
#' p <- plot_grid(p1, p2, nrow = 2) # plot is combined by plot_grid
#' # fix the size of panel for each plot
#' panel_fix(p, height = 1)
#' # rasterize the panel while keeping all labels and text in vector format
#' panel_fix(p, height = 1, raster = TRUE, dpi = 30)
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
                      width = NULL, height = NULL, margin = 0, units = "in",
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
          subgrob <- panel_fix_single(subgrob$grobs[[1]][["children"]], width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
        } else {
          subgrob <- panel_fix(subgrob, width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
        }
        gtable$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]] <- subgrob
        # print(paste0("plot_width:",plot_width," plot_height:",plot_height))
      }
      sum_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$width)
      sum_height <- convertHeight(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE) / as.numeric(gtable$grobs[[i]][["children"]][[j]]$vp$height)
      gtable <- panel_fix_single(gtable, panel_index = i, width = sum_width, height = sum_height, margin = ifelse(depth == 1, margin, 0), units = units, raster = FALSE, return_grob = TRUE)
    } else if (gtable$grobs[[i]]$name == "layout" || inherits(x, "patchwork")) {
      if (isTRUE(verbose)) {
        message("panel ", i, " is detected as generated by patchwork.")
      }
      # if (i == panel_index[1] && length(panel_index) > 1 && isTRUE(verbose)) {
      #   message("More than 2 panels detected. panel_fix may not work as expected.")
      # }
      subgrob <- gtable$grobs[[i]]
      if (length(subgrob[["children"]]) > 0 && all(sapply(subgrob[["children"]], function(x) inherits(x, "recordedGrob")))) {
        subgrob <- panel_fix_single(subgrob[["children"]], width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
      } else {
        subgrob <- panel_fix(subgrob, width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
      }
      gtable$grobs[[i]] <- subgrob
      sum_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE)
      sum_height <- convertHeight(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE)
      gtable <- panel_fix_single(gtable, panel_index = i, width = sum_width, height = sum_height, margin = ifelse(depth == 1, margin, 0), units = units, raster = FALSE, respect = TRUE, return_grob = TRUE)
    } else {
      # print("fix the gtable")
      gtable <- panel_fix_single(gtable, panel_index = i, width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi, return_grob = TRUE)
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
      filename <- R.utils::getAbsolutePath(save)
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

#' Fix the width/height of panels in a single plot object.
#'
#' @inheritParams panel_fix
#'
#' @importFrom ggplot2 ggsave zeroGrob
#' @importFrom gtable gtable_add_padding
#' @importFrom grid grob unit unitType convertWidth convertHeight convertUnit viewport grid.draw rasterGrob grobTree addGrob
#' @importFrom patchwork wrap_plots
#' @export
panel_fix_single <- function(x, panel_index = NULL, respect = NULL,
                             width = NULL, height = NULL, margin = 0, units = "in",
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
    check_R(c("png", "ragg", "grDevices"))
    for (i in seq_along(panel_index)) {
      index <- panel_index[i]
      do_sub <- TRUE
      depth <- 1
      while (isTRUE(do_sub) && depth <= 10) {
        g <- gtable$grobs[[index]]
        if (!inherits(gtable$grobs[[index]], "gtable")) {
          do_sub <- FALSE
          g_new <- g
        }
        depth <- depth + 1
      }
      # g <- g_bk
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
      grDevices::dev.off()
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
      filename <- R.utils::getAbsolutePath(save)
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
  p$data <- p$data[, intersect(colnames(p$data), vars), drop = FALSE]
  for (i in seq_along(p$layers)) {
    if (length(p$layers[[i]]$data) > 0) {
      p$layers[[i]]$data <- p$layers[[i]]$data[, intersect(colnames(p$layers[[i]]$data), vars), drop = FALSE]
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
#' @param colors Color vectors.
#' @param alpha Alpha level in [0,1]
#' @examples
#' colors <- c("red", "blue", "green")
#' adjcolors(colors, 0.5)
#'
#' library(scales)
#' alpha(colors, 0.5)
#'
#' show_palettes(list(colors, adjcolors(colors, 0.5), alpha(colors, 0.5)))
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
#' @param colors Color vectors.
#' @param mode Blend mode.
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
#' @param stat_plot_palette Color palette used in statistical plot
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
#' @param velocity_arrow_flank Flank of arrows in velocity plot.
#' @param streamline_L Typical length of a streamline in x and y units
#' @param streamline_minL Minimum length of segments to show.
#' @param streamline_res Resolution parameter (higher numbers increases the resolution).
#' @param streamline_n Number of points to draw.
#' @param streamlinewidth Size of streamline.
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
#' @return A single ggplot/patchwork object if combine = TRUE; otherwise, a list of ggplot objects
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", theme_use = "theme_blank", show_stat = FALSE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP") %>% panel_fix(height = 2, raster = TRUE, dpi = 30)
#'
#' # Label and highlight cell points
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP", label = TRUE, label_insitu = TRUE,
#'   cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"]
#' )
#' CellDimPlot(pancreas_sub,
#'   group.by = "SubCellType", split.by = "Phase", reduction = "UMAP", show_stat = FALSE,
#'   cells.highlight = TRUE, theme_use = "theme_blank", legend.position = "none"
#' )
#'
#' # Add a density layer
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", label = TRUE, add_density = TRUE)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", label = TRUE, add_density = TRUE, density_filled = TRUE)
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
#'   show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
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
#'   show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
#' )
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom dplyr group_by "%>%" .data
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_manual scale_fill_manual guides guide_legend geom_hex geom_path theme_void annotation_custom scale_linewidth_continuous after_stat
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color new_scale_fill new_scale
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom patchwork wrap_plots
#' @importFrom stats median loess aggregate
#' @importFrom utils askYesNo
#' @importFrom rlang %||%
#' @export
#'
CellDimPlot <- function(srt, group.by, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL,
                        show_na = FALSE, show_stat = TRUE,
                        pt.size = NULL, pt.alpha = 1, palette = "Paired", palcolor = NULL, bg_color = "grey80",
                        label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                        label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                        label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                        cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                        add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_palcolor = NULL,
                        lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75,
                        lineages_palette = "Dark2", lineages_palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                        lineages_linewidth = 1, lineages_line_bg = "white", lineages_line_bg_stroke = 0.5,
                        lineages_whiskers = FALSE, lineages_whiskers_linewidth = 0.5, lineages_whiskers_alpha = 0.5,
                        stat.by = NULL, stat_type = "percent", stat_plot_type = "pie", stat_plot_position = c("stack", "dodge"), stat_plot_size = 0.1,
                        stat_plot_palette = "Set1", stat_palcolor = NULL, stat_plot_alpha = 1, stat_plot_label = FALSE, stat_plot_label_size = 3,
                        graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
                        paga = NULL, paga_type = "connectivities", paga_node_size = 4,
                        paga_edge_threshold = 0.01, paga_edge_size = c(0.2, 1), paga_edge_color = "grey40", paga_edge_alpha = 0.5,
                        paga_transition_threshold = 0.01, paga_transition_size = c(0.2, 1), paga_transition_color = "black", paga_transition_alpha = 1, paga_show_transition = FALSE,
                        velocity = NULL, velocity_plot_type = "raw", velocity_n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50),
                        velocity_density = 1, velocity_smooth = 0.5, velocity_scale = 1, velocity_min_mass = 1, velocity_cutoff_perc = 5,
                        velocity_arrow_color = "black", velocity_arrow_angle = 20, velocity_arrow_flank = 0.8,
                        streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
                        streamlinewidth = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
                        streamline_bg_color = "white", streamline_bg_stroke = 0.5,
                        hex = FALSE, hex.linewidth = 0.5, hex.count = TRUE, hex.bins = 50, hex.binwidth = NULL,
                        raster = NULL, raster.dpi = c(512, 512),
                        aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                        legend.position = "right", legend.direction = "vertical",
                        theme_use = "theme_scp", theme_args = list(),
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in unique(c(group.by, split.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
    if (isTRUE(show_na) && any(is.na(srt@meta.data[[i]]))) {
      raw_levels <- unique(c(levels(srt@meta.data[[i]]), "NA"))
      srt@meta.data[[i]] <- as.character(srt@meta.data[[i]])
      srt@meta.data[[i]][is.na(srt@meta.data[[i]])] <- "NA"
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = raw_levels)
    }
  }
  for (l in lineages) {
    if (!l %in% colnames(srt@meta.data)) {
      stop(paste0(l, " is not in the meta.data of srt object."))
    }
  }
  if (!is.null(graph) && !graph %in% names(srt@graphs)) {
    stop("Graph ", graph, " is not exist in the srt object.")
  }
  if (!is.null(graph)) {
    graph <- srt@graphs[[graph]]
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }

  dat_meta <- srt@meta.data[, unique(c(group.by, split.by)), drop = FALSE]
  nlev <- sapply(dat_meta, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 100 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_use <- cbind(dat_dim, dat_meta[row.names(dat_dim), , drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) > 1e5)
  if (isTRUE(raster)) {
    check_R("exaexa/scattermore")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector")
    }
  }
  if (!is.null(stat.by)) {
    subplots <- CellStatPlot(srt,
      cells = cells,
      stat.by = stat.by, group.by = group.by, split.by = split.by,
      stat_type = stat_type, plot_type = stat_plot_type, position = stat_plot_position,
      palette = stat_plot_palette, palcolor = stat_palcolor, alpha = stat_plot_alpha,
      label = stat_plot_label, label.size = stat_plot_label_size,
      legend.position = "bottom", legend.direction = legend.direction,
      theme_use = theme_void, theme_args = theme_args,
      individual = TRUE, combine = FALSE
    )
  }
  if (!is.null(lineages)) {
    lineages_layers <- LineagePlot(srt,
      cells = cells,
      lineages = lineages, reduction = reduction, dims = dims,
      trim = lineages_trim, span = lineages_span,
      palette = lineages_palette, palcolor = lineages_palcolor, lineages_arrow = lineages_arrow,
      linewidth = lineages_linewidth, line_bg = lineages_line_bg, line_bg_stroke = lineages_line_bg_stroke,
      whiskers = lineages_whiskers, whiskers_linewidth = lineages_whiskers_linewidth, whiskers_alpha = lineages_whiskers_alpha,
      aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = "bottom", legend.direction = legend.direction,
      theme_use = theme_void, theme_args = theme_args,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(paga)) {
    if (split.by != "All_cells") {
      stop("paga can only plot on the non-split data")
    }
    paga_layers <- PAGAPlot(srt,
      cells = cells,
      paga = paga, type = paga_type, reduction = reduction, dims = dims,
      node_palette = palette, node_palcolor = palcolor, node_size = paga_node_size,
      edge_threshold = paga_edge_threshold, edge_size = paga_edge_size, edge_color = paga_edge_color, edge_alpha = paga_edge_alpha,
      transition_threshold = paga_transition_threshold, transition_size = paga_transition_size, transition_color = paga_transition_color, transition_alpha = paga_transition_alpha, show_transition = paga_show_transition,
      aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = "bottom", legend.direction = legend.direction,
      theme_use = theme_void, theme_args = theme_args,
      return_layer = TRUE
    )
    paga_layers <- paga_layers[!names(paga_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(velocity)) {
    if (split.by != "All_cells") {
      stop("velocity can only plot on the non-split data")
    }
    velocity_layers <- VelocityPlot(srt,
      cells = cells,
      reduction = reduction, dims = dims, velocity = velocity, plot_type = velocity_plot_type, group_by = group.by, group_palette = palette, group_palcolor = palcolor,
      n_neighbors = velocity_n_neighbors, density = velocity_density, smooth = velocity_smooth, scale = velocity_scale, min_mass = velocity_min_mass, cutoff_perc = velocity_cutoff_perc,
      arrow_color = velocity_arrow_color, arrow_angle = velocity_arrow_angle, arrow_flank = velocity_arrow_flank,
      streamline_L = streamline_L, streamline_minL = streamline_minL, streamline_res = streamline_res, streamline_n = streamline_n,
      streamlinewidth = streamlinewidth, streamline_alpha = streamline_alpha, streamline_color = streamline_color, streamline_palette = streamline_palette, streamline_palcolor = streamline_palcolor,
      streamline_bg_color = streamline_bg_color, streamline_bg_stroke = streamline_bg_stroke,
      aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = "bottom", legend.direction = legend.direction,
      theme_use = theme_void, theme_args = theme_args,
      return_layer = TRUE
    )
    velocity_layers <- velocity_layers[!names(velocity_layers) %in% c("lab_layer", "theme_layer")]
  }
  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  comb <- expand.grid(split = levels(dat_use[[split.by]]), group = group.by, stringsAsFactors = FALSE)
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["group"]])
  plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
    g <- comb[i, "group"]
    s <- comb[i, "split"]
    colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor, NA_keep = TRUE)
    dat <- dat_use
    cells_mask <- dat[[split.by]] != s
    # dat[[split.by]] <- NULL
    dat[[g]][cells_mask] <- NA
    legend_list <- list()
    labels_tb <- table(dat[[g]])
    labels_tb <- labels_tb[labels_tb != 0]
    cells.highlight_use <- cells.highlight
    if (isTRUE(cells.highlight_use)) {
      cells.highlight_use <- rownames(dat)[!is.na(dat[[g]])]
    }
    if (isTRUE(label_insitu)) {
      if (isTRUE(show_stat)) {
        label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
      } else {
        label_use <- paste0(names(labels_tb))
      }
    } else {
      if (isTRUE(label)) {
        if (isTRUE(show_stat)) {
          label_use <- paste0(seq_along(labels_tb), ": ", names(labels_tb), "(", labels_tb, ")")
        } else {
          label_use <- paste0(seq_along(labels_tb), ": ", names(labels_tb))
        }
      } else {
        if (isTRUE(show_stat)) {
          label_use <- paste0(names(labels_tb), "(", labels_tb, ")")
        } else {
          label_use <- paste0(names(labels_tb))
        }
      }
    }
    var_nm <- setNames(
      object = c("x", "y", "group.by"),
      nm = c(paste0(reduction_key, dims[1]), paste0(reduction_key, dims[2]), g)
    )
    colnames(dat)[colnames(dat) %in% names(var_nm)] <- var_nm[colnames(dat)[colnames(dat) %in% names(var_nm)]]
    dat[, "split.by"] <- s
    dat <- dat[order(dat[, "group.by"], decreasing = FALSE, na.last = FALSE), , drop = FALSE]
    naindex <- which(is.na(dat[, "group.by"]))
    naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
    dat <- dat[c(1:naindex, sample((min(naindex + 1, nrow(dat))):nrow(dat))), , drop = FALSE]
    if (isTRUE(show_stat)) {
      subtitle_use <- subtitle %||% paste0(s, " nCells:", sum(!is.na(dat[["group.by"]])))
    } else {
      subtitle_use <- subtitle
    }
    if (isTRUE(add_density)) {
      if (isTRUE(density_filled)) {
        filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_palcolor)
        density <- list(
          stat_density_2d(
            geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = after_stat(density)),
            contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
          ),
          scale_fill_gradientn(name = "Density", colours = filled_color),
          new_scale_fill()
        )
      } else {
        density <- geom_density_2d(aes(x = .data[["x"]], y = .data[["y"]]),
          color = density_color, inherit.aes = FALSE, show.legend = FALSE
        )
      }
    } else {
      density <- NULL
    }
    if (!is.null(graph)) {
      net_mat <- as.matrix(x = graph)[rownames(dat), rownames(dat)]
      net_mat[net_mat == 0] <- NA
      net_mat[upper.tri(net_mat)] <- NA
      net_df <- reshape2::melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
      net_df[, "value"] <- as.numeric(net_df[, "value"])
      net_df[, "Var1"] <- as.character(net_df[, "Var1"])
      net_df[, "Var2"] <- as.character(net_df[, "Var2"])
      net_df[, "x"] <- dat[net_df[, "Var1"], "x"]
      net_df[, "y"] <- dat[net_df[, "Var1"], "y"]
      net_df[, "xend"] <- dat[net_df[, "Var2"], "x"]
      net_df[, "yend"] <- dat[net_df[, "Var2"], "y"]
      net <- list(
        geom_segment(
          data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = value),
          color = edge_color, alpha = edge_alpha, show.legend = FALSE
        ),
        scale_linewidth_continuous(range = edge_size)
      )
    } else {
      net <- NULL
    }

    p <- ggplot(dat) +
      net +
      density +
      labs(title = title, subtitle = subtitle_use, x = xlab, y = ylab) +
      scale_x_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE), max(dat_dim[, paste0(reduction_key, dims[1])], na.rm = TRUE))) +
      scale_y_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE), max(dat_dim[, paste0(reduction_key, dims[2])], na.rm = TRUE))) +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction
      )
    if (split.by != "All_cells") {
      p <- p + facet_grid(. ~ split.by)
    }

    if (isTRUE(raster)) {
      p <- p + scattermore::geom_scattermore(
        data = dat[is.na(dat[, "group.by"]), , drop = FALSE],
        mapping = aes(x = .data[["x"]], y = .data[["y"]]), color = bg_color,
        pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
      ) + scattermore::geom_scattermore(
        data = dat[!is.na(dat[, "group.by"]), , drop = FALSE],
        mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
        pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
      )
    } else if (isTRUE(hex)) {
      check_R("hexbin")
      if (isTRUE(hex.count)) {
        p <- p + geom_hex(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["group.by"]], color = .data[["group.by"]], alpha = after_stat(count)),
          linewidth = hex.linewidth, bins = hex.bins, binwidth = hex.binwidth
        )
      } else {
        p <- p + geom_hex(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["group.by"]], color = .data[["group.by"]]),
          linewidth = hex.linewidth, bins = hex.bins, binwidth = hex.binwidth
        )
      }
    } else {
      p <- p + geom_point(
        mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
        size = pt.size, alpha = pt.alpha
      )
    }
    if (!is.null(cells.highlight_use) && !isTRUE(hex)) {
      cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
      if (nrow(cell_df) > 0) {
        if (isTRUE(raster)) {
          p <- p + scattermore::geom_scattermore(
            data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
            pointsize = floor(sizes.highlight) + stroke.highlight, alpha = alpha.highlight, pixels = raster.dpi
          ) +
            scattermore::geom_scattermore(
              data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
              pointsize = floor(sizes.highlight), alpha = alpha.highlight, pixels = raster.dpi
            )
        } else {
          p <- p + geom_point(
            data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
            size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
          ) +
            geom_point(
              data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
              size = sizes.highlight, alpha = alpha.highlight
            )
        }
      }
    }
    p <- p + scale_color_manual(
      name = paste0(g, ":"),
      values = colors[names(labels_tb)],
      labels = label_use,
      na.value = bg_color,
      guide = guide_legend(
        title.hjust = 0,
        keywidth = 0,
        keyheight = 0,
        default.unit = "inch",
        order = 1,
        override.aes = list(size = 4.5, alpha = 1)
      )
    ) + scale_fill_manual(
      name = paste0(g, ":"),
      values = colors[names(labels_tb)],
      labels = label_use,
      na.value = bg_color,
      guide = guide_legend(
        title.hjust = 0,
        keywidth = 0,
        keyheight = 0,
        default.unit = "inch",
        order = 1
      )
    )
    p_base <- p

    if (!is.null(stat.by)) {
      coor_df <- aggregate(p$data[, c("x", "y")], by = list(p$data[["group.by"]]), FUN = median)
      colnames(coor_df)[1] <- "group"
      x_range <- diff(layer_scales(p)$x$range$range)
      y_range <- diff(layer_scales(p)$y$range$range)
      stat_plot <- subplots[paste0(g, ":", levels(dat[, "group.by"]), ":", s)]
      names(stat_plot) <- levels(dat[, "group.by"])

      stat_plot_list <- list()
      for (i in seq_len(nrow(coor_df))) {
        stat_plot_list[[i]] <- annotation_custom(as_grob(stat_plot[[coor_df[i, "group"]]] + theme_void() + theme(legend.position = "none")),
          xmin = coor_df[i, "x"] - x_range * stat_plot_size / 2, ymin = coor_df[i, "y"] - y_range * stat_plot_size / 2,
          xmax = coor_df[i, "x"] + x_range * stat_plot_size / 2, ymax = coor_df[i, "y"] + y_range * stat_plot_size / 2
        )
      }
      p <- p + stat_plot_list
      legend_list[["stat.by"]] <- get_legend(stat_plot[[coor_df[i, "group"]]] + theme(legend.position = "bottom"))
    }
    if (!is.null(lineages)) {
      lineages_layers <- c(list(new_scale_color()), lineages_layers)
      suppressMessages({
        legend_list[["lineages"]] <- get_legend(ggplot() +
          lineages_layers +
          theme_scp(
            legend.position = "bottom",
            legend.direction = legend.direction
          ))
      })
      p <- suppressWarnings({
        p + lineages_layers + theme(legend.position = "none")
      })
      if (is.null(legend_list[["lineages"]])) {
        legend_list["lineages"] <- list(NULL)
      }
    }
    if (!is.null(paga)) {
      paga_layers <- c(list(new_scale_color()), paga_layers)
      if (g != paga$groups) {
        suppressMessages({
          legend_list[["paga"]] <- get_legend(ggplot() +
            paga_layers +
            theme_scp(
              legend.position = "bottom",
              legend.direction = legend.direction
            ))
        })
      }
      p <- suppressWarnings({
        p + paga_layers + theme(legend.position = "none")
      })
      if (is.null(legend_list[["paga"]])) {
        legend_list["paga"] <- list(NULL)
      }
    }
    if (!is.null(velocity)) {
      velocity_layers <- c(list(new_scale_color()), list(new_scale("size")), velocity_layers)
      if (velocity_plot_type != "raw") {
        suppressMessages({
          legend_list[["velocity"]] <- get_legend(ggplot() +
            velocity_layers +
            theme_scp(
              legend.position = "bottom",
              legend.direction = legend.direction
            ))
        })
      }
      p <- suppressWarnings({
        p + velocity_layers + theme(legend.position = "none")
      })
      if (is.null(legend_list[["velocity"]])) {
        legend_list["velocity"] <- list(NULL)
      }
    }
    if (isTRUE(label)) {
      label_df <- aggregate(p$data[, c("x", "y")], by = list(p$data[["group.by"]]), FUN = median)
      colnames(label_df)[1] <- "label"
      label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
      if (!isTRUE(label_insitu)) {
        label_df[, "label"] <- seq_len(nrow(label_df))
      }
      if (isTRUE(label_repel)) {
        p <- p + geom_point(
          data = label_df, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
          color = label_point_color, size = label_point_size
        ) + geom_text_repel(
          data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
          fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
          point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
        )
      } else {
        p <- p + geom_text_repel(
          data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
          fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
          point.size = NA, max.overlaps = 100,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
        )
      }
    }
    if (length(legend_list) > 0) {
      legend_list <- legend_list[!sapply(legend_list, is.null)]
      legend_base <- get_legend(p_base + theme_scp(
        legend.position = "bottom",
        legend.direction = legend.direction
      ))
      if (legend.direction == "vertical") {
        legend <- do.call(cbind, c(list(base = legend_base), legend_list))
      } else {
        legend <- do.call(rbind, c(list(base = legend_base), legend_list))
      }
      gtable <- as_grob(p + theme(legend.position = "none"))
      gtable <- add_grob(gtable, legend, legend.position)
      p <- wrap_plots(gtable)
    }
    return(p)
  })

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' Visualize feature values on a 2-dimensional reduction plot
#'
#' Plotting cell points on a reduced 2D plane and coloring according to the values of the features.
#'
#' @param srt A Seurat object.
#' @param features Vector of features to plot. Features can be gene names in Assay or names of numeric columns in meta.data.
#' @param reduction Which dimensionality reduction to use. If not specified, will use the reduction returned by \code{\link{DefaultReduction}}.
#' @param split.by Name of a column in meta.data to split plot by.
#' @param palette Name of a color palette name collected in SCP.
#' @param palcolor Custom colors used to create a color palette.
#' @param pt.size Point size for plotting.
#' @param pt.alpha Point transparency.
#' @param keep_scale How to handle the color scale across multiple plots. Options are:
#' \itemize{
#'   \item{NULL (no scaling):}{ Each individual plot is scaled to the maximum expression value of the feature in the condition provided to 'split.by'. Be aware setting NULL will result in color scales that are not comparable between plots.}
#'   \item{"feature" (default; by row/feature scaling):}{ The plots for each individual feature are scaled to the maximum expression of the feature across the conditions provided to 'split.by'.}
#'   \item{"all" (universal scaling):}{ The plots for all features and conditions are scaled to the maximum expression value for the feature with the highest overall expression.}
#' }
#' @param cells.highlight A vector of cell names to highlight.
#' @param cols.highlight Color used to highlight the cells.
#' @param sizes.highlight Size of highlighted cells.
#' @param alpha.highlight Transparency of highlighted cell points.
#' @param stroke.highlight Border width of highlighted cell points.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top").
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical")
#' @param combine Combine plots into a single \code{patchwork} object. If \code{FALSE}, return a list of ggplot objects.
#' @param nrow Number of rows in the combined plot.
#' @param ncol Number of columns in the combined plot.
#' @param byrow Logical value indicating if the plots should be arrange by row (default) or by column.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions.
#' @param slot Which slot to pull expression data from? Default is \code{data}.
#' @param assay Which assay to pull expression data from. If \code{NULL}, will use the assay returned by \code{\link{DefaultAssay}}.
#' @param show_stat Whether to show statistical information on the plot.
#' @param calculate_coexp Whether to calculate the co-expression value (geometric mean) of the features.
#' @param compare_features Whether to show the values of multiple features on a single plot.
#' @param color_blend_mode Blend mode to use when \code{compare_features = TRUE}
#' @param bg_cutoff Background cutoff. Points with feature values lower than the cutoff will be considered as background and will be colored with \code{bg_color}.
#' @param bg_color Color value for background points.
#' @param lower_quantile,upper_quantile,lower_cutoff,upper_cutoff Vector of minimum and maximum cutoff values or quantile values for each feature.
#' @param add_density Whether to add a density layer on the plot.
#' @param density_color Color of the density contours lines.
#' @param density_filled Whether to add filled contour bands instead of contour lines.
#' @param density_filled_palette Color palette used to fill contour bands.
#' @param density_filled_palcolor Custom colors used to fill contour bands.
#' @param label Whether the feature name is labeled in the center of the location of cells wieh high expression.
#' @param label.size Size of labels.
#' @param label.fg Foreground color of label.
#' @param label.bg Background color of label.
#' @param label.bg.r Background ratio of label.
#' @param label_insitu Whether the labels is feature names instead of numbers. Valid only when \code{compare_features = TRUE}.
#' @param label_repel Logical value indicating whether the label is repel away from the center location.
#' @param label_repulsion Force of repulsion between overlapping text labels. Defaults to 20.
#' @param label_point_size Size of the center points.
#' @param label_point_color Color of the center points
#' @param label_segment_color Color of the line segment for labels.
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
#' @param graph Specify the graph name to add edges between cell neighbors to the plot.
#' @param edge_size Size of edges.
#' @param edge_alpha Transparency of edges.
#' @param edge_color Color of edges.
#' @param hex Whether to chane the plot type from point to the hexagonal bin.
#' @param hex.bins Number of hexagonal bins.
#' @param hex.binwidth Hexagonal bin width.
#' @param hex.color Border color of hexagonal bins.
#' @param hex.linewidth Border width of hexagonal bins.
#' @param raster Convert points to raster format, default is NULL which automatically rasterizes if plotting more than 100,000 cells
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).
#' @param theme_use Theme used. Can be a character string or a theme function. For example, \code{"theme_blank"} or \code{ggplot2::theme_classic}.
#' @param aspect.ratio Aspect ratio of the panel.
#' @param title The text for the title.
#' @param subtitle The text for the subtitle for the plot which will be displayed below the title.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param force Whether to force drawing regardless of the number of features greater than 100.
#' @param cells Subset cells to plot.
#' @param theme_args Other arguments passed to the \code{theme_use}.
#' @param seed Random seed set for reproducibility
#'
#' @return A single ggplot object if combine = TRUE; otherwise, a list of ggplot objects
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP")
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", bg_cutoff = -Inf)
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", theme_use = "theme_blank", show_stat = FALSE)
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP") %>% panel_fix(height = 2, raster = TRUE, dpi = 30)
#'
#' # Label and highlight cell points
#' FeatureDimPlot(pancreas_sub,
#'   features = "Rbp4", reduction = "UMAP", label = TRUE,
#'   cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Delta"]
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = "Rbp4", split.by = "Phase", reduction = "UMAP", show_stat = FALSE,
#'   cells.highlight = TRUE, theme_use = "theme_blank", legend.position = "none"
#' )
#'
#' # Add a density layer
#' FeatureDimPlot(pancreas_sub, features = "Rbp4", reduction = "UMAP", label = TRUE, add_density = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "Rbp4", reduction = "UMAP", label = TRUE, add_density = TRUE, density_filled = TRUE)
#'
#' # Chane the plot type from point to the hexagonal bin
#' FeatureDimPlot(pancreas_sub, features = "Rbp4", reduction = "UMAP", hex = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "Rbp4", reduction = "UMAP", hex = TRUE, hex.bins = 20)
#'
#' # Show lineages on the plot based on the pseudotime
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' FeatureDimPlot(pancreas_sub, features = "Lineage3", reduction = "UMAP", lineages = "Lineage3")
#' FeatureDimPlot(pancreas_sub, features = "Lineage3", reduction = "UMAP", lineages = "Lineage3", lineages_whiskers = TRUE)
#' FeatureDimPlot(pancreas_sub, features = "Lineage3", reduction = "UMAP", lineages = "Lineage3", lineages_span = 0.1)
#'
#' FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg", "Sst", "Ghrl"), reduction = "UMAP")
#' FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg", "Sst", "Ghrl"), reduction = "UMAP", lower_quantile = 0, upper_quantile = 0.8)
#' FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg", "Sst", "Ghrl"), reduction = "UMAP", lower_cutoff = 1, upper_cutoff = 4)
#' FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg", "Sst", "Ghrl"), reduction = "UMAP", bg_cutoff = 2, lower_cutoff = 2, upper_cutoff = 4)
#' FeatureDimPlot(pancreas_sub, c("Ins1", "Gcg", "Sst", "Ghrl"), reduction = "UMAP", keep_scale = "all")
#' FeatureDimPlot(pancreas_sub, c("Sst", "Ghrl"), split.by = "Phase", reduction = "UMAP", keep_scale = "feature")
#'
#' FeatureDimPlot(pancreas_sub,
#'   features = c("Ins1", "Gcg", "Sst", "Ghrl"), pt.size = 1,
#'   compare_features = TRUE, color_blend_mode = "blend",
#'   label = TRUE, label_insitu = TRUE
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "blend", title = "blend",
#'   label = TRUE, label_insitu = TRUE
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "average", title = "average",
#'   label = TRUE, label_insitu = TRUE
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "screen", title = "screen",
#'   label = TRUE, label_insitu = TRUE
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "multiply", title = "multiply",
#'   label = TRUE, label_insitu = TRUE
#' )
#'
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom dplyr group_by "%>%" .data
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colourbar geom_hex stat_summary_hex geom_path scale_linewidth_continuous after_stat
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom gtable gtable_add_cols
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @importFrom methods slot
#' @export
#'
FeatureDimPlot <- function(srt, features, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL, slot = "data", assay = NULL,
                           show_stat = TRUE,
                           palette = ifelse(isTRUE(compare_features), "Set1", "Spectral"), palcolor = NULL,
                           pt.size = NULL, pt.alpha = 1, bg_cutoff = 0, bg_color = "grey80",
                           keep_scale = NULL, lower_quantile = 0, upper_quantile = 0.99, lower_cutoff = NULL, upper_cutoff = NULL,
                           add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_palcolor = NULL,
                           cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                           calculate_coexp = FALSE, compare_features = FALSE, color_blend_mode = c("blend", "average", "screen", "multiply"),
                           label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                           label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                           label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                           lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75,
                           lineages_palette = "Dark2", lineages_palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                           lineages_linewidth = 1, lineages_line_bg = "white", lineages_line_bg_stroke = 0.5,
                           lineages_whiskers = FALSE, lineages_whiskers_linewidth = 0.5, lineages_whiskers_alpha = 0.5,
                           graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
                           hex = FALSE, hex.linewidth = 0.5, hex.color = "grey90", hex.bins = 50, hex.binwidth = NULL,
                           raster = NULL, raster.dpi = c(512, 512),
                           aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  color_blend_mode <- match.arg(color_blend_mode)
  require("ggrepel", quietly = TRUE)

  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in c(split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
  }
  for (l in lineages) {
    if (!l %in% colnames(srt@meta.data)) {
      stop(paste0(l, " is not in the meta.data of srt object."))
    }
  }
  if (!is.null(graph) && !graph %in% names(srt@graphs)) {
    stop("Graph ", graph, " is not exist in the srt object.")
  }
  if (!is.null(graph)) {
    graph <- srt@graphs[[graph]]
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }

  features <- unique(features)
  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names:", paste0(intersect(features_gene, features_meta), collapse = ","))
  }

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression", immediate. = TRUE)
    }
    status <- check_DataType(srt, slot = slot, assay = assay)
    message("Data type detected in ", slot, " slot: ", status)
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      srt@meta.data[["CoExp"]] <- apply(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE], 2, function(x) exp(mean(log(x))))
    } else if (status == "log_normalized_counts") {
      srt@meta.data[["CoExp"]] <- apply(expm1(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE]), 2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene) > 0) {
    if (all(rownames(srt@assays[[assay]]) %in% features_gene)) {
      dat_gene <- t(slot(srt@assays[[assay]], slot))
    } else {
      dat_gene <- t(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE])
    }
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as.matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  reduction_key <- srt@reductions[[reduction]]@key
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_sp <- srt@meta.data[, split.by, drop = FALSE]
  dat_use <- cbind(dat_dim, dat_sp[row.names(dat_dim), , drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) > 1e5)
  if (isTRUE(raster)) {
    check_R("exaexa/scattermore")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector")
    }
  }

  if (!is.null(lineages)) {
    lineages_layers <- LineagePlot(srt,
      lineages = lineages, reduction = reduction, dims = dims,
      trim = lineages_trim, span = lineages_span,
      palette = lineages_palette, palcolor = lineages_palcolor, lineages_arrow = lineages_arrow,
      linewidth = lineages_linewidth, line_bg = lineages_line_bg, line_bg_stroke = lineages_line_bg_stroke,
      whiskers = lineages_whiskers, whiskers_linewidth = lineages_whiskers_linewidth, whiskers_alpha = lineages_whiskers_alpha,
      aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = legend.position, legend.direction = legend.direction,
      theme_use = theme_void, theme_args = theme_args,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
  }

  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (identical(theme_use, "theme_blank")) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  if (isTRUE(compare_features) && length(features) > 1) {
    dat_all <- cbind(dat_use, dat_exp[row.names(dat_use), features, drop = FALSE])
    dat_split <- split.data.frame(dat_all, dat_all[[split.by]])
    plist <- lapply(levels(dat_sp[[split.by]]), function(s) {
      dat <- dat_split[[ifelse(split.by == "All_cells", 1, s)]][, , drop = FALSE]
      for (f in features) {
        dat[, f][dat[, f] <= bg_cutoff] <- NA
        if (any(is.infinite(dat[, f]))) {
          dat[, f][which(dat[, f] == max(dat[, f], na.rm = TRUE))] <- max(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
          dat[, f][which(dat[, f] == min(dat[, f], na.rm = TRUE))] <- min(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
        }
      }
      var_nm <- setNames(
        object = c("x", "y"),
        nm = c(paste0(reduction_key, dims[1]), paste0(reduction_key, dims[2]))
      )
      colnames(dat)[colnames(dat) %in% names(var_nm)] <- var_nm[colnames(dat)[colnames(dat) %in% names(var_nm)]]
      dat[, "split.by"] <- s
      dat[, "features"] <- paste(features, collapse = "|")
      subtitle_use <- subtitle %||% s
      colors <- palette_scp(features, type = "discrete", palette = palette, palcolor = palcolor)
      colors_list <- list()
      value_list <- list()
      pal_list <- list()
      temp_geom <- list()
      legend_list <- list()
      for (i in seq_along(colors)) {
        colors_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = TRUE, matched = TRUE, palcolor = c(adjcolors(colors[i], 0.1), colors[i]))
        pal_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = FALSE, matched = FALSE, palcolor = c(adjcolors(colors[i], 0.1), colors[i]))
        value_list[[i]] <- seq(min(dat[, names(colors)[i]], na.rm = TRUE), max(dat[, names(colors)[i]], na.rm = TRUE), length.out = 100)
        temp_geom[[i]] <- list(
          geom_point(data = dat, mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[[names(colors)[i]]])),
          scale_color_gradientn(
            colours = pal_list[[i]],
            values = rescale(value_list[[i]]), na.value = bg_color,
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0)
          ),
          new_scale_color()
        )
        legend_list[[i]] <- get_legend(ggplot(dat, aes(x = .data[["x"]], y = .data[["y"]])) +
          temp_geom[[i]] +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction
          ))
      }
      for (j in seq_len(nrow(dat))) {
        dat[j, "color_blend"] <- blendcolors(sapply(colors_list, function(x) x[j]), mode = color_blend_mode)
      }
      dat["color_value"] <- colSums(col2rgb(dat[, "color_blend"]))
      dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_value"] <- NA
      dat <- dat[order(dat[, "color_value"], decreasing = TRUE, na.last = FALSE), , drop = FALSE]
      dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_blend"] <- bg_color
      cells.highlight_use <- cells.highlight
      if (isTRUE(cells.highlight_use)) {
        cells.highlight_use <- rownames(dat)[dat[["color_blend"]] != bg_color]
      }
      if (!is.null(graph)) {
        net_mat <- as.matrix(x = graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- reshape2::melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
        net_df[, "value"] <- as.numeric(net_df[, "value"])
        net_df[, "Var1"] <- as.character(net_df[, "Var1"])
        net_df[, "Var2"] <- as.character(net_df[, "Var2"])
        net_df[, "x"] <- dat[net_df[, "Var1"], "x"]
        net_df[, "y"] <- dat[net_df[, "Var1"], "y"]
        net_df[, "xend"] <- dat[net_df[, "Var2"], "x"]
        net_df[, "yend"] <- dat[net_df[, "Var2"], "y"]
        net <- list(
          geom_segment(
            data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = value),
            color = edge_color, alpha = edge_alpha, show.legend = FALSE
          ),
          scale_linewidth_continuous(range = edge_size)
        )
      } else {
        net <- NULL
      }
      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_palcolor)
          density <- list(
            stat_density_2d(
              geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = after_stat(density)),
              contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
            ),
            scale_fill_gradientn(name = "Density", colours = filled_color),
            new_scale_fill()
          )
        } else {
          density <- geom_density_2d(aes(x = .data[["x"]], y = .data[["y"]]),
            color = density_color,
            inherit.aes = FALSE
          )
        }
      } else {
        density <- NULL
      }
      p <- ggplot(dat) +
        net +
        density +
        labs(title = title, subtitle = s, x = xlab, y = ylab) +
        scale_x_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[1])], na.rm = TRUE), max(dat_use[, paste0(reduction_key, dims[1])], na.rm = TRUE))) +
        scale_y_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE), max(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE)))
      if (split.by == "All_cells") {
        p <- p + facet_grid(. ~ features)
      } else {
        p <- p + facet_grid(split.by ~ features)
      }
      p <- p + do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = "none",
          legend.direction = legend.direction
        )
      if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(
          data = dat[dat[, "color_blend"] == bg_color, , drop = FALSE],
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
          pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
        ) +
          scattermore::geom_scattermore(
            data = dat[dat[, "color_blend"] != bg_color, , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
            pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
          ) +
          scale_color_identity() +
          new_scale_color()
      } else {
        p <- p + geom_point(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
          size = pt.size, alpha = pt.alpha
        ) +
          scale_color_identity() +
          new_scale_color()
      }

      if (!is.null(cells.highlight_use)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
        if (nrow(cell_df) > 0) {
          if (isTRUE(raster)) {
            p <- p + scattermore::geom_scattermore(
              data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
              pointsize = floor(sizes.highlight) + stroke.highlight, alpha = alpha.highlight, pixels = raster.dpi
            ) +
              scattermore::geom_scattermore(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
                pointsize = floor(sizes.highlight), alpha = alpha.highlight, pixels = raster.dpi
              ) +
              scale_color_identity() +
              new_scale_color()
          } else {
            p <- p +
              geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
                size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
              ) +
              geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
                size = sizes.highlight, alpha = alpha.highlight
              ) +
              scale_color_identity() +
              new_scale_color()
          }
        }
      }

      legend2 <- NULL
      if (isTRUE(label)) {
        label_df <- reshape2::melt(p$data, measure.vars = features)
        label_df <- label_df %>%
          group_by(variable) %>%
          filter(value >= quantile(value, 0.95, na.rm = TRUE) & value <= quantile(value, 0.99, na.rm = TRUE)) %>%
          reframe(x = median(.data[["x"]]), y = median(.data[["y"]])) %>%
          as.data.frame()
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
        label_df[, "rank"] <- seq_len(nrow(label_df))
        if (isTRUE(label_insitu)) {
          if (isTRUE(label_repel)) {
            p <- p + geom_point(
              data = label_df, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
              color = label_point_color, size = label_point_size
            ) + geom_text_repel(
              data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]], color = .data[["label"]]),
              fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
              point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
              color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, show.legend = FALSE
            )
          } else {
            p <- p + geom_text_repel(
              data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]], color = .data[["label"]]),
              fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
              point.size = NA, max.overlaps = 100,
              color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, show.legend = FALSE
            )
          }
          p <- p + scale_color_manual(
            name = "Label:",
            values = adjcolors(colors[label_df$label], 0.5),
            labels = label_df$label,
            na.value = bg_color
          )
        } else {
          if (isTRUE(label_repel)) {
            p <- p + geom_point(
              data = label_df, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
              color = "black", size = pt.size + 1
            ) +
              geom_text_repel(
                data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["rank"]], color = .data[["label"]]),
                fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
                point.size = pt.size + 1, max.overlaps = 100, force = label_repulsion,
                bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, key_glyph = "point"
              )
          } else {
            p <- p + geom_text_repel(
              data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["rank"]], color = .data[["label"]]),
              fontface = "bold", min.segment.length = 0, segment.colour = label_segment_color,
              point.size = NA, max.overlaps = 100,
              bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, key_glyph = "point"
            )
          }
          p <- p + scale_color_manual(
            name = "Label:",
            values = adjcolors(colors[label_df$label], 0.5),
            labels = paste(label_df$rank, label_df$label, sep = ": "),
            na.value = bg_color
          ) +
            guides(colour = guide_legend(override.aes = list(color = colors[label_df$label]), order = 1)) +
            theme(legend.position = "none")
          legend2 <- get_legend(p +
            do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              legend.position = legend.position,
              legend.direction = legend.direction
            ))
        }
      }

      legend_nrow <- min(ceiling(sqrt(length(legend_list))), 3)
      total <- length(legend_list)
      leg_list <- list()
      n <- 1
      for (i in 1:total) {
        if (i == 1 || is.null(leg)) {
          leg <- legend_list[[i]]
        } else {
          leg <- cbind(leg, legend_list[[i]])
        }
        if (i %% legend_nrow == 0) {
          leg_list[[n]] <- leg
          leg <- NULL
          n <- n + 1
        }
        if (i %% legend_nrow != 0 && i == total) {
          ncol_insert <- dim(leg_list[[n - 1]])[2] - dim(leg)[2]
          for (col_insert in 1:ncol_insert) {
            leg <- gtable_add_cols(leg, sum(leg_list[[n - 1]]$widths) / ncol_insert, -1)
          }
          leg_list[[n]] <- leg
        }
      }
      legend <- do.call(rbind, leg_list)
      if (!is.null(lineages)) {
        lineages_layers <- c(list(new_scale_color()), lineages_layers)
        suppressMessages({
          legend_curve <- get_legend(ggplot() +
            lineages_layers +
            theme_scp())
        })
        legend <- add_grob(legend, legend_curve, "top")
        p <- suppressMessages({
          p + lineages_layers + theme(legend.position = "none")
        })
      }

      gtable <- as_grob(p)
      gtable <- add_grob(gtable, legend, legend.position)
      if (!is.null(legend2)) {
        gtable <- add_grob(gtable, legend2, legend.position)
      }
      p <- wrap_plots(gtable)
      return(p)
    })
    names(plist) <- paste0(levels(dat_sp[[split.by]]), ":", paste0(features, collapse = "|"))
  } else {
    comb <- expand.grid(split = levels(dat_sp[[split.by]]), feature = features, stringsAsFactors = FALSE)
    rownames(comb) <- paste0(comb[["split"]], ":", comb[["feature"]])
    dat_all <- cbind(dat_use, dat_exp[row.names(dat_use), features, drop = FALSE])
    dat_split <- split.data.frame(dat_all, dat_all[[split.by]])
    colors <- palette_scp(type = "continuous", palette = palette, palcolor = palcolor)
    plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
      f <- comb[i, "feature"]
      s <- comb[i, "split"]
      dat <- dat_split[[ifelse(split.by == "All_cells", 1, s)]][, c(colnames(dat_use), f), drop = FALSE]
      dat[, f][dat[, f] <= bg_cutoff] <- NA
      if (any(is.infinite(dat[, f]))) {
        dat[, f][dat[, f] == max(dat[, f], na.rm = TRUE)] <- max(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
        dat[, f][dat[, f] == min(dat[, f], na.rm = TRUE)] <- min(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
      }
      var_nm <- setNames(
        object = c("x", "y", "value"),
        nm = c(paste0(reduction_key, dims[1]), paste0(reduction_key, dims[2]), f)
      )
      colnames(dat)[colnames(dat) %in% names(var_nm)] <- var_nm[colnames(dat)[colnames(dat) %in% names(var_nm)]]
      dat <- dat[order(dat[, "value"], method = "radix", decreasing = FALSE, na.last = FALSE), , drop = FALSE]
      dat[, "features"] <- f
      cells.highlight_use <- cells.highlight
      if (isTRUE(cells.highlight_use)) {
        cells.highlight_use <- rownames(dat)[!is.na(dat[["value"]])]
      }
      legend_list <- list()
      if (isTRUE(show_stat)) {
        subtitle_use <- subtitle %||% paste0(s, " nPos:", sum(dat[["value"]] > 0, na.rm = TRUE), ", ", round(sum(dat[["value"]] > 0, na.rm = TRUE) / nrow(dat) * 100, 2), "%")
      } else {
        subtitle_use <- subtitle
      }
      if (all(is.na(dat[["value"]]))) {
        colors_value <- rep(0, 100)
      } else {
        if (is.null(keep_scale)) {
          colors_value <- seq(lower_cutoff %||% quantile(dat[["value"]], lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(dat[["value"]], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
        } else {
          if (keep_scale == "feature") {
            colors_value <- seq(lower_cutoff %||% quantile(dat_exp[, f], lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(dat_exp[, f], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
          }
          if (keep_scale == "all") {
            all_values <- as.matrix(dat_exp[, features])
            colors_value <- seq(lower_cutoff %||% quantile(all_values, lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(all_values, upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
          }
        }
      }
      dat[which(dat[, "value"] > max(colors_value, na.rm = TRUE)), "value"] <- max(colors_value, na.rm = TRUE)
      dat[which(dat[, "value"] < min(colors_value, na.rm = TRUE)), "value"] <- min(colors_value, na.rm = TRUE)
      if (!is.null(graph)) {
        net_mat <- as.matrix(x = graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- reshape2::melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
        net_df[, "value"] <- as.numeric(net_df[, "value"])
        net_df[, "Var1"] <- as.character(net_df[, "Var1"])
        net_df[, "Var2"] <- as.character(net_df[, "Var2"])
        net_df[, "x"] <- dat[net_df[, "Var1"], "x"]
        net_df[, "y"] <- dat[net_df[, "Var1"], "y"]
        net_df[, "xend"] <- dat[net_df[, "Var2"], "x"]
        net_df[, "yend"] <- dat[net_df[, "Var2"], "y"]
        net <- list(
          geom_segment(
            data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, linewidth = value),
            color = edge_color, alpha = edge_alpha, show.legend = FALSE
          ),
          scale_linewidth_continuous(range = edge_size)
        )
      } else {
        net <- NULL
      }
      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_palcolor)
          density <- list(
            stat_density_2d(
              geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = after_stat(density)),
              contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
            ),
            scale_fill_gradientn(name = "Density", colours = filled_color),
            new_scale_fill()
          )
        } else {
          density <- geom_density_2d(aes(x = .data[["x"]], y = .data[["y"]]),
            color = density_color,
            inherit.aes = FALSE
          )
        }
      } else {
        density <- NULL
      }
      p <- ggplot(dat) +
        net +
        density +
        labs(title = title, subtitle = subtitle_use, x = xlab, y = ylab) +
        scale_x_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[1])], na.rm = TRUE), max(dat_use[, paste0(reduction_key, dims[1])], na.rm = TRUE))) +
        scale_y_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE), max(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE))) +
        guides(color = guide_colourbar(
          barwidth = 0.9,
          barheight = 4,
          frame.colour = "black",
          ticks.colour = "black"
        )) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(
          data = dat[is.na(dat[, "value"]), , drop = FALSE],
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
          pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
        ) +
          scattermore::geom_scattermore(
            data = dat[!is.na(dat[, "value"]), , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
            pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
          )
      } else if (isTRUE(hex)) {
        check_R("hexbin")
        dat_na <- dat[is.na(dat[["value"]]), , drop = FALSE]
        dat_hex <- dat[!is.na(dat[["value"]]), , drop = FALSE]
        if (nrow(dat_na) > 0) {
          p <- p + geom_hex(
            data = dat[is.na(dat[["value"]]), , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]]),
            fill = bg_color, color = hex.color,
            linewidth = hex.linewidth, bins = hex.bins, binwidth = hex.binwidth
          )
        }
        if (nrow(dat_hex) > 0) {
          p <- p + stat_summary_hex(
            data = dat_hex,
            mapping = aes(x = .data[["x"]], y = .data[["y"]], z = .data[["value"]]),
            color = hex.color, linewidth = hex.linewidth, bins = hex.bins, binwidth = hex.binwidth
          )
          if (all(is.na(dat[["value"]]))) {
            p <- p + scale_fill_gradient(
              name = "", na.value = bg_color
            )
          } else {
            p <- p + scale_fill_gradientn(
              name = "", colours = colors, values = rescale(colors_value), limits = range(colors_value), na.value = bg_color
            )
          }
          p <- p + new_scale_fill()
        }
      } else {
        p <- p + geom_point(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
          size = pt.size, alpha = pt.alpha
        )
      }
      if (!is.null(cells.highlight_use) && !isTRUE(hex)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
        if (nrow(cell_df) > 0) {
          if (isTRUE(raster)) {
            p <- p + scattermore::geom_scattermore(
              data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
              pointsize = floor(sizes.highlight) + stroke.highlight, alpha = alpha.highlight, pixels = raster.dpi
            ) +
              scattermore::geom_scattermore(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
                pointsize = floor(sizes.highlight), alpha = alpha.highlight, pixels = raster.dpi
              )
          } else {
            p <- p +
              geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
                size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
              ) +
              geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
                size = sizes.highlight, alpha = alpha.highlight
              )
          }
        }
      }
      if (nrow(dat) > 0) {
        if (split.by == "All_cells") {
          p <- p + facet_grid(. ~ features)
        } else {
          p <- p + facet_grid(formula(paste0(split.by, "~features")))
        }
      }
      if (all(is.na(dat[["value"]]))) {
        p <- p + scale_colour_gradient(
          name = "", na.value = bg_color, aesthetics = c("color")
        )
      } else {
        p <- p + scale_color_gradientn(
          name = "", colours = colors, values = rescale(colors_value), limits = range(colors_value), na.value = bg_color, aesthetics = c("color")
        )
      }
      p <- p + guides(
        color = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)
      )
      p_base <- p

      if (!is.null(lineages)) {
        lineages_layers <- c(list(new_scale_color()), lineages_layers)
        suppressMessages({
          legend_list[["lineages"]] <- get_legend(ggplot() +
            lineages_layers +
            theme_scp(
              legend.position = "bottom",
              legend.direction = legend.direction
            ))
        })
        p <- suppressWarnings({
          p + lineages_layers + theme(legend.position = "none")
        })
        if (is.null(legend_list[["lineages"]])) {
          legend_list["lineages"] <- list(NULL)
        }
      }
      if (isTRUE(label)) {
        label_df <- p$data %>%
          filter(value >= quantile(value, 0.95, na.rm = TRUE) & value <= quantile(value, 0.99, na.rm = TRUE)) %>%
          reframe(x = median(.data[["x"]]), y = median(.data[["y"]])) %>%
          as.data.frame()
        label_df[, "label"] <- f
        label_df[, "rank"] <- seq_len(nrow(label_df))
        if (isTRUE(label_repel)) {
          p <- p + annotate(
            geom = "point", x = label_df[["x"]], y = label_df[["y"]],
            color = "black", size = pt.size + 1
          ) + annotate(
            geom = "text_repel", x = label_df[["x"]], y = label_df[["y"]], label = label_df[["label"]],
            fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
            point.size = pt.size + 1, max.overlaps = 100, force = label_repulsion,
            color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size
          )
        } else {
          p <- p + annotate(
            geom = "text_repel", x = label_df[["x"]], y = label_df[["y"]], label = label_df[["label"]],
            fontface = "bold",
            point.size = NA, max.overlaps = 100,
            color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size
          )
        }
      }
      if (length(legend_list) > 0) {
        legend_list <- legend_list[!sapply(legend_list, is.null)]
        legend_base <- get_legend(p_base)
        if (legend.direction == "vertical") {
          legend <- do.call(cbind, c(list(base = legend_base), legend_list))
        } else {
          legend <- do.call(rbind, c(list(base = legend_base), legend_list))
        }
        gtable <- as_grob(p + theme(legend.position = "none"))
        gtable <- add_grob(gtable, legend, legend.position)
        p <- wrap_plots(gtable)
      }
      return(p)
    })
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' 3D-Dimensional reduction plot for cell classification visualization.
#'
#' Plotting cell points on a reduced 3D space and coloring according to the groups of the cells.
#'
#' @inheritParams CellDimPlot
#' @param shape.highlight Shape of the cell to highlight. See \href{https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol}{scattergl-marker-symbol}
#' @param width Width in pixels, defaults to automatic sizing.
#' @param height Height in pixels, defaults to automatic sizing.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' CellDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#'
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#' CellDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D", lineages = "Lineage1")
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom utils askYesNo
#' @export
CellDimPlot3D <- function(srt, group.by, reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                          palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = 1.5,
                          cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                          lineages = NULL, lineage_palette = "Dark2", span = 0.75, arrow_reverse = FALSE,
                          width = NULL, height = NULL, save = NULL, force = FALSE) {
  check_R("plotly")
  bg_color <- col2hex(bg_color)
  cols.highlight <- col2hex(cols.highlight)

  for (i in c(group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
  }
  for (l in lineages) {
    if (!l %in% colnames(srt@meta.data)) {
      stop(paste0(l, " is not in the meta.data of srt object."))
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt, min_dim = 3)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction, min_dim = 3)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (ncol(srt@reductions[[reduction]]@cell.embeddings) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  reduction_key <- srt@reductions[[reduction]]@key
  if (is.null(axis_labs) || length(axis_labs) != 3) {
    xlab <- paste0(reduction_key, dims[1])
    ylab <- paste0(reduction_key, dims[2])
    zlab <- paste0(reduction_key, dims[3])
  } else {
    xlab <- axis_labs[1]
    ylab <- axis_labs[2]
    zlab <- axis_labs[3]
  }
  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    check_R("htmlwidgets")
    if (!grepl(".html$", save)) {
      stop("'save' must be a string with .html as a suffix.")
    }
  }

  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_use <- cbind(dat_dim[colnames(srt@assays[[1]]), , drop = FALSE], srt@meta.data[colnames(srt@assays[[1]]), , drop = FALSE])
  nlev <- sapply(dat_use[, group.by, drop = FALSE], nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 100 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  if (!is.null(lineages)) {
    dat_lineages <- srt@meta.data[, unique(lineages), drop = FALSE]
    dat_use <- cbind(dat_use, dat_lineages[row.names(dat_use), , drop = FALSE])
  }
  if (!is.factor(dat_use[[group.by]])) {
    dat_use[[group.by]] <- factor(dat_use[[group.by]], levels = unique(dat_use[[group.by]]))
  }
  dat_use[["group.by"]] <- dat_use[[group.by]]
  if (any(is.na(dat_use[[group.by]]))) {
    n <- as.character(dat_use[[group.by]])
    n[is.na(n)] <- "NA"
    dat_use[[group.by]] <- factor(n, levels = c(levels(dat_use[[group.by]]), "NA"))
  }

  dat_use[["color"]] <- dat_use[[group.by]]
  colors <- palette_scp(dat_use[["group.by"]], palette = palette, palcolor = palcolor, NA_color = bg_color, NA_keep = TRUE)

  dat_use[[paste0(reduction_key, dims[1], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[1])]]
  dat_use[[paste0(reduction_key, dims[2], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[2])]]
  dat_use[[paste0(reduction_key, dims[3], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[3])]]
  cells.highlight_use <- cells.highlight
  if (isTRUE(cells.highlight_use)) {
    cells.highlight_use <- rownames(dat_use)[dat_use[[group.by]] != "NA"]
  }
  if (!is.null(cells.highlight_use)) {
    cells.highlight_use <- cells.highlight_use[cells.highlight_use %in% rownames(dat_use)]
    dat_use_highlight <- dat_use[cells.highlight_use, , drop = FALSE]
  }

  p <- plotly::plot_ly(data = dat_use, width = width, height = height)
  p <- plotly::add_trace(
    p = p,
    data = dat_use,
    x = dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
    y = dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
    z = dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
    text = paste0(
      "Cell:", rownames(dat_use),
      "\ngroup.by:", dat_use[["group.by"]],
      "\ncolor:", dat_use[["color"]]
    ),
    type = "scatter3d",
    mode = "markers",
    color = dat_use[[group.by]],
    colors = colors,
    marker = list(size = pt.size),
    showlegend = TRUE,
    visible = TRUE
  )
  if (!is.null(cells.highlight_use)) {
    p <- plotly::add_trace(
      p = p,
      x = dat_use_highlight[[paste0(reduction_key, dims[1], "All_cells")]],
      y = dat_use_highlight[[paste0(reduction_key, dims[2], "All_cells")]],
      z = dat_use_highlight[[paste0(reduction_key, dims[3], "All_cells")]],
      text = paste0(
        "Cell:", rownames(dat_use_highlight),
        "\ngroup.by:", dat_use_highlight[["group.by"]],
        "\ncolor:", dat_use_highlight[["color"]]
      ),
      type = "scatter3d",
      mode = "markers",
      marker = list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight),
      name = "highlight",
      showlegend = TRUE,
      visible = TRUE
    )
  }
  if (!is.null(lineages)) {
    for (l in lineages) {
      dat_sub <- dat_use[!is.na(dat_use[[l]]), , drop = FALSE]
      dat_sub <- dat_sub[order(dat_sub[[l]]), , drop = FALSE]

      xlo <- loess(formula(paste(paste0(reduction_key, dims[1], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      ylo <- loess(formula(paste(paste0(reduction_key, dims[2], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      zlo <- loess(formula(paste(paste0(reduction_key, dims[3], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      dat_smooth <- data.frame(x = xlo$fitted, y = ylo$fitted, z = zlo$fitted)
      dat_smooth <- dat_smooth[dat_smooth[["x"]] <= max(dat_use[[paste0(reduction_key, dims[1], "All_cells")]], na.rm = TRUE) & dat_smooth[["x"]] >= min(dat_use[[paste0(reduction_key, dims[1], "All_cells")]], na.rm = TRUE), , drop = FALSE]
      dat_smooth <- dat_smooth[dat_smooth[["y"]] <= max(dat_use[[paste0(reduction_key, dims[2], "All_cells")]], na.rm = TRUE) & dat_smooth[["y"]] >= min(dat_use[[paste0(reduction_key, dims[2], "All_cells")]], na.rm = TRUE), , drop = FALSE]
      dat_smooth <- dat_smooth[dat_smooth[["z"]] <= max(dat_use[[paste0(reduction_key, dims[3], "All_cells")]], na.rm = TRUE) & dat_smooth[["z"]] >= min(dat_use[[paste0(reduction_key, dims[3], "All_cells")]], na.rm = TRUE), , drop = FALSE]
      dat_smooth <- unique(na.omit(dat_smooth))
      p <- plotly::add_trace(
        p = p,
        x = dat_smooth[, "x"],
        y = dat_smooth[, "y"],
        z = dat_smooth[, "z"],
        text = paste0(
          "Lineage:", l
        ),
        type = "scatter3d",
        mode = "lines",
        line = list(width = 6, color = palette_scp(x = lineages, palette = lineage_palette)[l], reverscale = FALSE),
        name = l,
        showlegend = TRUE,
        visible = TRUE
      )
    }
  }

  p <- plotly::layout(
    p = p,
    title = list(
      text = paste0("Total", " (nCells:", nrow(dat_use), ")"),
      font = list(size = 16, color = "black"),
      y = 0.95
    ),
    font = list(size = 12, color = "black"),
    showlegend = TRUE,
    legend = list(
      itemsizing = "constant",
      y = 0.5,
      x = 1,
      xanchor = "left",
      alpha = 1
    ),
    scene = list(
      xaxis = list(title = xlab, range = c(min(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE))),
      yaxis = list(title = ylab, range = c(min(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE))),
      zaxis = list(title = zlab, range = c(min(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE))),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    autosize = FALSE
  )

  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    htmlwidgets::saveWidget(
      widget = plotly::as_widget(p),
      file = save
    )
    unlink(gsub("\\.html", "_files", save), recursive = TRUE)
  }

  return(p)
}

#' 3D-Dimensional reduction plot for gene expression visualization.
#'
#' Plotting cell points on a reduced 3D space and coloring according to the gene expression in the cells.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' FeatureDimPlot3D(pancreas_sub, features = c("Ghrl", "Ins1", "Gcg", "Ins2"), reduction = "StandardpcaUMAP3D")
#' @inheritParams FeatureDimPlot
#' @inheritParams CellDimPlot3D
#'
#' @importFrom methods slot
#' @importFrom rlang "%||%"
#' @export
FeatureDimPlot3D <- function(srt, features = NULL, reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                             split.by = NULL, slot = "data", assay = NULL,
                             calculate_coexp = FALSE,
                             pt.size = 1.5, cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                             width = NULL, height = NULL, save = NULL, force = FALSE) {
  check_R("plotly")
  cols.highlight <- col2hex(cols.highlight)

  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt@meta.data[[split.by]] <- factor("All_cells")
  }
  for (i in split.by) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt, min_dim = 3)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction, min_dim = 3)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (ncol(srt@reductions[[reduction]]@cell.embeddings) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- colnames(srt@assays[[1]])
  }
  reduction_key <- srt@reductions[[reduction]]@key
  if (is.null(axis_labs) || length(axis_labs) != 3) {
    xlab <- paste0(reduction_key, dims[1])
    ylab <- paste0(reduction_key, dims[2])
    zlab <- paste0(reduction_key, dims[3])
  } else {
    xlab <- axis_labs[1]
    ylab <- axis_labs[2]
    zlab <- axis_labs[3]
  }
  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    check_R("htmlwidgets")
    if (!grepl(".html$", save)) {
      stop("'save' must be a string with .html as a suffix.")
    }
  }

  features <- unique(features)
  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression", immediate. = TRUE)
    }
    status <- check_DataType(srt, slot = slot, assay = assay)
    message("Data type detected in ", slot, " slot: ", status)
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      srt@meta.data[["CoExp"]] <- apply(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE], 2, function(x) exp(mean(log(x))))
    } else if (status == "log_normalized_counts") {
      srt@meta.data[["CoExp"]] <- apply(expm1(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE]), 2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene > 0)) {
    dat_gene <- t(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE])
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta > 0)) {
    dat_meta <- as.matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  dat_sp <- srt@meta.data[, split.by, drop = FALSE]
  dat_dim <- srt@reductions[[reduction]]@cell.embeddings
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt@assays[[1]])
  dat_use <- cbind(dat_exp, dat_dim[rownames(dat_exp), , drop = FALSE], dat_sp[rownames(dat_exp), , drop = FALSE])

  dat_use[[paste0(reduction_key, dims[1], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[1])]]
  dat_use[[paste0(reduction_key, dims[2], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[2])]]
  dat_use[[paste0(reduction_key, dims[3], "All_cells")]] <- dat_use[[paste0(reduction_key, dims[3])]]
  for (i in levels(dat_use[[split.by]])) {
    dat_use[[paste0(reduction_key, dims[1], i)]] <- ifelse(dat_use[[split.by]] == i, dat_use[[paste0(reduction_key, dims[1])]], NA)
    dat_use[[paste0(reduction_key, dims[2], i)]] <- ifelse(dat_use[[split.by]] == i, dat_use[[paste0(reduction_key, dims[2])]], NA)
    dat_use[[paste0(reduction_key, dims[3], i)]] <- ifelse(dat_use[[split.by]] == i, dat_use[[paste0(reduction_key, dims[3])]], NA)
  }
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(dat_use)]
    dat_use_highlight <- dat_use[cells.highlight, , drop = FALSE]
    for (i in levels(dat_use_highlight[[split.by]])) {
      dat_use_highlight[[paste0(reduction_key, dims[1], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[1])]], NA)
      dat_use_highlight[[paste0(reduction_key, dims[2], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[2])]], NA)
      dat_use_highlight[[paste0(reduction_key, dims[3], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[3])]], NA)
    }
  }

  p <- plotly::plot_ly(data = dat_use, width = width, height = height)
  p <- plotly::add_trace(
    p = p,
    data = dat_use,
    x = dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
    y = dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
    z = dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
    text = paste0(
      "Cell:", rownames(dat_use),
      "\nExp:", round(dat_use[[features[1]]], 3),
      "\nsplit.by:", dat_use[[split.by]]
    ),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      color = dat_use[[features[1]]],
      colorbar = list(title = list(text = features[1], font = list(color = "black", size = 14)), len = 0.5),
      size = pt.size,
      showscale = TRUE
    ),
    name = "All_cells",
    showlegend = TRUE,
    visible = TRUE
  )
  if (!is.null(cells.highlight)) {
    p <- plotly::add_trace(
      p = p,
      x = dat_use_highlight[[paste0(reduction_key, dims[1], "All_cells")]],
      y = dat_use_highlight[[paste0(reduction_key, dims[2], "All_cells")]],
      z = dat_use_highlight[[paste0(reduction_key, dims[3], "All_cells")]],
      text = paste0(
        "Cell:", rownames(dat_use_highlight),
        "\nExp:", round(dat_use_highlight[[features[1]]], 3),
        "\nsplit.by:", dat_use_highlight[[split.by]]
      ),
      type = "scatter3d",
      mode = "markers",
      marker = list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight),
      name = "highlight",
      showlegend = TRUE,
      visible = TRUE
    )
  }

  split_option <- list()
  genes_option <- list()
  for (i in 0:nlevels(dat_use[[split.by]])) {
    sp <- ifelse(i == 0, "All_cells", levels(dat_use[[split.by]])[i])
    ncells <- ifelse(i == 0, nrow(dat_use), table(dat_use[[split.by]])[sp])
    if (i != 0 && sp == "All_cells") {
      next
    }
    x <- list(dat_use[[paste0(reduction_key, dims[1], sp)]])
    y <- list(dat_use[[paste0(reduction_key, dims[2], sp)]])
    z <- list(dat_use[[paste0(reduction_key, dims[3], sp)]])
    name <- sp
    if (!is.null(cells.highlight)) {
      x <- c(x, list(dat_use_highlight[[paste0(reduction_key, dims[1], sp)]]))
      y <- c(y, list(dat_use_highlight[[paste0(reduction_key, dims[2], sp)]]))
      z <- c(z, list(dat_use_highlight[[paste0(reduction_key, dims[3], sp)]]))
      name <- c(sp, "highlight")
    }
    split_option[[i + 1]] <- list(
      method = "update",
      args = list(list(
        x = x,
        y = y,
        z = z,
        name = name,
        visible = TRUE
      ), list(title = list(
        text = paste0(sp, " (nCells:", ncells, ")"),
        font = list(size = 16, color = "black"),
        y = 0.95
      ))),
      label = sp
    )
  }
  for (j in seq_along(features)) {
    marker <- list(
      color = dat_use[[features[j]]],
      colorbar = list(title = list(text = features[j], font = list(color = "black", size = 14)), len = 0.5),
      size = pt.size,
      showscale = TRUE
    )

    if (!is.null(cells.highlight)) {
      marker <- list(marker, list(size = sizes.highlight, color = cols.highlight, symbol = shape.highlight))
    }
    genes_option[[j]] <- list(
      method = "update",
      args = list(list(
        text = list(paste0(
          "Cell:", rownames(dat_use),
          "\nExp:", round(dat_use[[features[j]]], 3),
          "\nsplit.by:", dat_use[[split.by]]
        )),
        marker = marker
      )),
      label = features[j]
    )
  }

  p <- plotly::layout(
    p = p,
    title = list(
      text = paste0("All_cells", " (nCells:", nrow(dat_use), ")"),
      font = list(size = 16, color = "black"),
      y = 0.95
    ),
    showlegend = TRUE,
    legend = list(
      itemsizing = "constant",
      y = -0.2,
      x = 0.5,
      xanchor = "center"
    ),
    scene = list(
      xaxis = list(title = xlab, range = c(min(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[1])]], na.rm = TRUE))),
      yaxis = list(title = ylab, range = c(min(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[2])]], na.rm = TRUE))),
      zaxis = list(title = zlab, range = c(min(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE), max(dat_use[[paste0(reduction_key, dims[3])]], na.rm = TRUE))),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    updatemenus = list(
      list(
        y = 0.67,
        buttons = split_option
      ),
      list(
        y = 0.33,
        buttons = genes_option
      )
    ),
    autosize = FALSE
  )

  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    htmlwidgets::saveWidget(
      widget = plotly::as_widget(p),
      file = save
    )
    unlink(gsub("\\.html", "_files", save), recursive = TRUE)
  }

  return(p)
}

#' Statistical plot of features
#'
#' @param srt
#' @param group.by
#' @param split.by
#' @param bg.by
#' @param cells
#' @param keep_empty
#' @param slot
#' @param assay
#' @param palette
#' @param palcolor
#' @param alpha
#' @param bg_palette
#' @param bg_palcolor
#' @param bg_apha
#' @param add_box
#' @param box_width
#' @param add_point
#' @param pt.color
#' @param pt.size
#' @param pt.alpha
#' @param jitter.width
#' @param cells.highlight
#' @param cols.highlight
#' @param sizes.highlight
#' @param alpha.highlight
#' @param calculate_coexp
#' @param compare_features
#' @param y.max
#' @param same.y.lims
#' @param y.trans
#' @param sort
#' @param stack
#' @param fill.by
#' @param flip
#' @param comparisons
#' @param ref_group
#' @param pairwise_method
#' @param multiplegroup_comparisons
#' @param multiple_method
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param force
#' @param individual
#' @param plot_type
#' @param y.nbreaks
#' @param y.min
#' @param stat.by
#' @param box_color
#' @param box_ptsize
#' @param add_trend
#' @param trend_color
#' @param trend_linewidth
#' @param trend_ptsize
#' @param sig_label
#' @param theme_args
#' @param seed
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType") %>% panel_fix(height = 1, width = 2)
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", plot_type = "box")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", plot_type = "bar")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", plot_type = "dot")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", plot_type = "col")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", add_box = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", add_point = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", add_trend = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", split.by = "Phase")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", split.by = "Phase", add_box = TRUE, add_trend = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", split.by = "Phase", comparisons = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("Rbp4", "Pyy"), group.by = "SubCellType", fill.by = "expression", palette = "Blues", same.y.lims = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("Rbp4", "Pyy"), group.by = "SubCellType", multiplegroup_comparisons = TRUE)
#' FeatureStatPlot(pancreas_sub, stat.by = c("Rbp4", "Pyy"), group.by = "SubCellType", comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")))
#' FeatureStatPlot(pancreas_sub, stat.by = c("Rbp4", "Pyy"), group.by = "SubCellType", comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")), sig_label = "p.format")
#' FeatureStatPlot(pancreas_sub, stat.by = c("Rbp4", "Pyy"), group.by = "SubCellType", bg.by = "CellType", add_box = TRUE, stack = TRUE)
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   legend.position = "top", legend.direction = "horizontal",
#'   group.by = "SubCellType", bg.by = "CellType", stack = TRUE
#' )
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   fill.by = "feature", plot_type = "box",
#'   group.by = "SubCellType", bg.by = "CellType", stack = TRUE, flip = TRUE
#' ) %>% panel_fix_single(width = 8, height = 5) # Because the plot is made by combining, we want to adjust the overall height and width
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom gtable gtable_add_cols gtable_add_rows gtable_add_grob gtable_add_padding
#' @importFrom ggplot2 geom_blank geom_violin geom_rect geom_boxplot geom_count geom_col geom_vline geom_hline layer_data layer_scales position_jitterdodge position_dodge stat_summary scale_x_discrete element_line element_text element_blank annotate mean_sdl after_stat scale_shape_identity
#' @importFrom grid grobHeight grobWidth
#' @importFrom rlang %||%
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix rowSums
#' @export
FeatureStatPlot <- function(srt, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, fill.by = c("group", "feature", "expression"),
                            cells = NULL, slot = c("data", "counts"), assay = NULL, keep_empty = FALSE, individual = FALSE,
                            plot_type = c("violin", "box", "bar", "dot", "col"),
                            palette = "Paired", palcolor = NULL, alpha = 1,
                            bg_palette = "Paired", bg_palcolor = NULL, bg_apha = 0.2,
                            add_box = FALSE, box_color = "black", box_width = 0.1, box_ptsize = 2,
                            add_point = FALSE, pt.color = "grey30", pt.size = NULL, pt.alpha = 1, jitter.width = 0.5,
                            add_trend = FALSE, trend_color = "black", trend_linewidth = 1, trend_ptsize = 2,
                            add_stat = c("none", "mean", "median"), stat_color = "black", stat_size = 1,
                            cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1, alpha.highlight = 1,
                            calculate_coexp = FALSE, compare_features = FALSE,
                            same.y.lims = FALSE, y.min = NULL, y.max = NULL, y.trans = "identity", y.nbreaks = 5,
                            sort = FALSE, stack = FALSE, flip = FALSE,
                            comparisons = NULL, ref_group = NULL, pairwise_method = "wilcox.test", sig_label = c("p.signif", "p.format"),
                            multiplegroup_comparisons = FALSE, multiple_method = "kruskal.test",
                            aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = "Expression level",
                            legend.position = "right", legend.direction = "vertical",
                            theme_use = "theme_scp", theme_args = list(),
                            combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  plot_type <- match.arg(plot_type)
  fill.by <- match.arg(fill.by)
  slot <- match.arg(slot)
  sig_label <- match.arg(sig_label)
  add_stat <- match.arg(add_stat)

  meta.data <- srt@meta.data
  assay <- assay %||% DefaultAssay(srt)
  exp.data <- slot(srt@assays[[assay]], slot)

  if (plot_type == "col") {
    if (isTRUE(add_box) || isTRUE(add_point) || isTRUE(add_trend) || isTRUE(add_stat != "none")) {
      warning("Cannot add other layers when plot_type is 'col'", immediate. = TRUE)
      add_box <- add_point <- add_trend <- FALSE
    }
  }
  if ((isTRUE(multiplegroup_comparisons) || length(comparisons) > 0) && plot_type %in% c("col")) {
    warning("Cannot add comparison when plot_type is 'col'", immediate. = TRUE)
    multiplegroup_comparisons <- FALSE
    comparisons <- NULL
  }
  if (isTRUE(comparisons) && is.null(split.by)) {
    stop("'split.by' must provided when comparisons=TRUE")
  }

  if (nrow(meta.data) == 0) {
    stop("meta.data is empty.")
  }
  if (is.null(group.by)) {
    group.by <- "No.group.by"
    xlab <- ""
    meta.data[[group.by]] <- factor("All")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    meta.data[[split.by]] <- factor("")
  }
  for (i in unique(c(group.by, split.by, bg.by))) {
    if (!i %in% colnames(meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }
  bg_map <- NULL
  if (!is.null(bg.by)) {
    for (g in group.by) {
      df_table <- table(meta.data[[g]], meta.data[[bg.by]])
      if (max(rowSums(df_table > 0), na.rm = TRUE) > 1) {
        stop("'group.by' must be a part of 'bg.by'")
      } else {
        bg_map[[g]] <- setNames(colnames(df_table)[apply(df_table, 1, function(x) which(x > 0))], rownames(df_table))
      }
    }
  } else {
    for (g in group.by) {
      bg_map[[g]] <- setNames(levels(meta.data[[g]]), levels(meta.data[[g]]))
    }
  }
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- colnames(srt@assays[[1]])
  }
  if (!is.null(cells.highlight) && isFALSE(add_point)) {
    warning("'cells.highlight' is valid only when add_point=TRUE.", immediate. = TRUE)
  }
  if (isTRUE(stack) & isTRUE(sort)) {
    message("Set sort to FALSE when stack is TRUE")
    sort <- FALSE
  }
  if (isTRUE(multiplegroup_comparisons) || length(comparisons) > 0) {
    check_R("ggpubr")
    require("ggpubr", quietly = TRUE)
    ncomp <- sapply(comparisons, length)
    if (any(ncomp > 2)) {
      stop("'comparisons' must be a list in which all elements must be vectors of length 2")
    }
  }

  stat.by <- unique(stat.by)
  features_drop <- stat.by[!stat.by %in% c(rownames(exp.data), colnames(meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    stat.by <- stat.by[!stat.by %in% features_drop]
  }

  features_gene <- stat.by[stat.by %in% rownames(exp.data)]
  features_meta <- stat.by[stat.by %in% colnames(meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 1) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression", immediate. = TRUE)
    }
    status <- check_DataType(data = exp.data)
    message("Data type detected in ", slot, " slot: ", status)
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      meta.data[["CoExp"]] <- apply(exp.data[features_gene, , drop = FALSE], 2, function(x) exp(mean(log(x))))
    } else if (status == "log_normalized_counts") {
      meta.data[["CoExp"]] <- apply(expm1(exp.data[features_gene, , drop = FALSE]), 2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    stat.by <- c(stat.by, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene) > 0) {
    if (all(rownames(srt@assays[[assay]]) %in% features_gene)) {
      dat_gene <- t(exp.data)
    } else {
      dat_gene <- t(exp.data[features_gene, , drop = FALSE])
    }
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as.matrix(meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  stat.by <- unique(stat.by[stat.by %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'stat.by' must be type of numeric variable.")
  }
  dat_group <- meta.data[, unique(c(group.by, bg.by, split.by)), drop = FALSE]
  dat_use <- cbind(dat_group, dat_exp[row.names(dat_group), , drop = FALSE])
  if (!is.null(cells)) {
    dat_group <- dat_group[intersect(rownames(dat_group), cells), , drop = FALSE]
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_group), 0.5)
  }

  nlev <- sapply(dat_group, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 100 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  if (isTRUE(same.y.lims) && is.null(y.max)) {
    y.max <- max(as.matrix(dat_use[, stat.by, drop = FALSE])[is.finite(as.matrix(dat_use[, stat.by, drop = FALSE]))], na.rm = TRUE)
  }
  if (isTRUE(same.y.lims) && is.null(y.min)) {
    y.min <- min(as.matrix(dat_use[, stat.by, drop = FALSE])[is.finite(as.matrix(dat_use[, stat.by, drop = FALSE]))], na.rm = TRUE)
  }

  plist <- list()
  plist_stack <- list()

  comb_list <- list()
  comb <- expand.grid(group_name = group.by, stat_name = stat.by, stringsAsFactors = FALSE)
  if (isTRUE(individual)) {
    for (g in group.by) {
      comb_list[[g]] <- merge(comb, expand.grid(
        group_name = g, group_element = levels(dat_use[[g]]),
        split_name = levels(dat_use[[split.by]]), stringsAsFactors = FALSE
      ),
      by = "group_name", all = FALSE
      )
    }
  } else {
    for (g in group.by) {
      comb_list[[g]] <- merge(comb, expand.grid(
        group_name = g, group_element = list(levels(dat_use[[g]])),
        split_name = list(levels(dat_use[[split.by]])), stringsAsFactors = FALSE
      ),
      by = "group_name", all = FALSE
      )
    }
  }
  comb <- do.call(rbind, comb_list)
  rownames(comb) <- paste0(
    comb[["stat_name"]], ":", comb[["group_name"]], ":",
    sapply(comb[["group_element"]], function(x) paste0(x, collapse = ",")), ":",
    sapply(comb[["split_name"]], function(x) paste0(x, collapse = ","))
  )

  plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
    g <- comb[i, "group_name"]
    f <- comb[i, "stat_name"]
    single_group <- comb[[i, "group_element"]]
    sp <- comb[[i, "split_name"]]
    xlab <- xlab %||% g
    ylab <- ylab %||% "Expression level"
    if (identical(theme_use, "theme_blank")) {
      theme_args[["xlab"]] <- xlab
      theme_args[["ylab"]] <- ylab
    }
    if (fill.by == "feature") {
      colors <- palette_scp(stat.by, palette = palette, palcolor = palcolor)
    }
    if (fill.by == "group") {
      if (split.by != "All_cells") {
        colors <- palette_scp(levels(dat_use[[split.by]]), palette = palette, palcolor = palcolor)
      } else {
        colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor)
      }
    }
    if (fill.by == "expression") {
      median_values <- aggregate(dat_use[, stat.by, drop = FALSE], by = list(dat_use[[g]], dat_use[[split.by]]), FUN = median)
      rownames(median_values) <- paste0(median_values[, 1], "-", median_values[, 2])
      colors <- palette_scp(unlist(median_values[, stat.by]), type = "continuous", palette = palette, palcolor = palcolor)
      colors_limits <- range(median_values[, stat.by])
    }

    dat <- dat_use[dat_use[[g]] %in% single_group & dat_use[[split.by]] %in% sp, c(colnames(dat_group), f)]
    dat[[g]] <- factor(dat[[g]], levels = levels(dat[[g]])[levels(dat[[g]]) %in% dat[[g]]])
    if (!is.null(bg.by)) {
      bg <- bg.by
      bg_color <- palette_scp(levels(dat[[bg]]), palette = bg_palette, palcolor = bg_palcolor)
    } else {
      bg <- g
      bg_color <- palette_scp(levels(dat[[bg]]), palcolor = bg_palcolor %||% rep(c("transparent", "grey85"), nlevels(dat[[bg]])))
    }
    dat[["bg.by"]] <- dat[[bg]]
    var_nm <- setNames(
      object = c("value", "group.by", "split.by"),
      nm = c(f, g, split.by)
    )
    colnames(dat)[colnames(dat) %in% names(var_nm)] <- var_nm[colnames(dat)[colnames(dat) %in% names(var_nm)]]
    stat <- table(dat[, "group.by"], dat[, "split.by"])
    stat_drop <- which(stat == 1, arr.ind = TRUE)
    if (nrow(stat_drop) > 0) {
      for (j in 1:nrow(stat_drop)) {
        dat <- dat[!(dat[, "group.by"] == rownames(stat)[stat_drop[j, 1]] & dat[, "split.by"] == colnames(stat)[stat_drop[j, 2]]), , drop = FALSE]
        rownames(stat)[stat_drop[j, 1]]
      }
    }
    dat[, "features"] <- rep(f, nrow(dat))
    if (nrow(dat) > 0 && ((is.character(x = sort) && nchar(x = sort) > 0) || sort)) {
      df_sort <- aggregate(dat[, "value", drop = FALSE], by = list(dat[["group.by"]]), median)
      if (is.character(sort) && sort == "increasing") {
        decreasing <- FALSE
      } else {
        decreasing <- TRUE
      }
      sortlevel <- as.character(df_sort[order(df_sort[["value"]], decreasing = decreasing), 1])
      dat[, "group.by"] <- factor(dat[, "group.by"], levels = sortlevel)
    }
    if (fill.by == "feature") {
      dat[, "fill.by"] <- rep(f, nrow(dat))
      keynm <- "Features"
    }
    if (fill.by == "group") {
      dat[, "fill.by"] <- if (split.by == "All_cells") dat[, "group.by"] else dat[, "split.by"]
      keynm <- ifelse(split.by == "All_cells", g, split.by)
    }
    if (fill.by == "expression") {
      dat[, "fill.by"] <- median_values[paste0(dat[["group.by"]], "-", dat[["split.by"]]), f]
      keynm <- "Median expression"
    }
    if (split.by != "All_cells") {
      levels_order <- levels(dat[["split.by"]])
    } else {
      levels_order <- levels(dat[["group.by"]])
    }
    if (fill.by == "feature") {
      levels_order <- unique(stat.by)
    }

    group_comb <- expand.grid(x = levels(dat[["split.by"]]), y = levels(dat[["group.by"]]))
    dat[["group.unique"]] <- head(
      factor(paste("sp", dat[["split.by"]], "gp", dat[["group.by"]], sep = "-"),
        levels = paste("sp", group_comb[[1]], "gp", group_comb[[2]], sep = "-")
      ),
      nrow(dat)
    )
    y_max_use <- y.max %||% suppressWarnings(max(dat[, "value"][is.finite(x = dat[, "value"])], na.rm = TRUE))
    y_min_use <- y.min %||% suppressWarnings(min(dat[, "value"][is.finite(x = dat[, "value"])], na.rm = TRUE))
    dat <- dat[dat[["value"]] >= y_min_use & dat[["value"]] <= y_max_use, , drop = FALSE]
    dat <- dat[order(dat[["group.unique"]]), , drop = FALSE]

    if (isTRUE(flip)) {
      dat[["group.by"]] <- factor(dat[["group.by"]], levels = rev(levels(dat[["group.by"]])))
      aspect.ratio <- 1 / aspect.ratio
      if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
        aspect.ratio <- NULL
      }
    }

    if (plot_type == "col") {
      if (isTRUE(flip)) {
        dat[["cell"]] <- rev(seq_len(nrow(dat)))
      } else {
        dat[["cell"]] <- seq_len(nrow(dat))
      }
      p <- ggplot(dat, aes(
        x = .data[["cell"]], y = .data[["value"]], fill = .data[["fill.by"]]
      ))
    } else {
      p <- ggplot(dat, aes(
        x = .data[["group.by"]], y = .data[["value"]], fill = .data[["fill.by"]]
      ))
    }
    if (isTRUE(flip)) {
      p <- p + do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          strip.text.x = element_text(angle = 45),
          panel.grid.major.x = element_line(color = "grey", linetype = 2),
          legend.position = legend.position,
          legend.direction = legend.direction
        ) + coord_flip()
    } else {
      p <- p + do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text.y = element_text(angle = 0),
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          legend.position = legend.position,
          legend.direction = legend.direction
        )
    }

    if (isFALSE(individual)) {
      if (plot_type == "col") {
        x_index <- split(dat[["cell"]], dat[["group.by"]])
        bg_data <- as.data.frame(t(sapply(x_index, range)))
        colnames(bg_data) <- c("xmin", "xmax")
        bg_data[["group.by"]] <- names(x_index)
        bg_data[["xmin"]] <- ifelse(bg_data[["xmin"]] == min(bg_data[["xmax"]]), -Inf, bg_data[["xmin"]] - 0.5)
        bg_data[["xmax"]] <- ifelse(bg_data[["xmax"]] == max(bg_data[["xmax"]]), Inf, bg_data[["xmax"]] + 0.5)
        bg_data[["ymin"]] <- -Inf
        bg_data[["ymax"]] <- Inf
        bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[["group.by"]])]]
      } else {
        bg_data <- unique(dat[, "group.by", drop = FALSE])
        bg_data[["x"]] <- as.numeric(bg_data[["group.by"]])
        bg_data[["xmin"]] <- ifelse(bg_data[["x"]] == min(bg_data[["x"]]), -Inf, bg_data[["x"]] - 0.5)
        bg_data[["xmax"]] <- ifelse(bg_data[["x"]] == max(bg_data[["x"]]), Inf, bg_data[["x"]] + 0.5)
        bg_data[["ymin"]] <- -Inf
        bg_data[["ymax"]] <- Inf
        bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[["group.by"]])]]
      }
      bg_layer <- geom_rect(data = bg_data, xmin = bg_data[["xmin"]], xmax = bg_data[["xmax"]], ymin = bg_data[["ymin"]], ymax = bg_data[["ymax"]], fill = bg_data[["fill"]], alpha = bg_apha, inherit.aes = FALSE)
      p <- p + bg_layer
    }

    if (plot_type %in% c("bar", "col")) {
      p <- p + geom_hline(yintercept = 0, linetype = 2)
    }
    if (plot_type == "violin") {
      p <- p + geom_violin(scale = "width", trim = TRUE, alpha = alpha, position = position_dodge())
    }
    if (plot_type == "box") {
      add_box <- FALSE
      p <- p + geom_boxplot(
        mapping = aes(group = .data[["group.unique"]]),
        position = position_dodge(width = 0.9), color = "black", width = 0.8, outlier.shape = NA
      ) +
        stat_summary(
          fun = median, geom = "point", mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9), color = "black", fill = "white", size = 1.5, shape = 21,
        )
    }
    if (plot_type == "bar") {
      p <- p + stat_summary(
        fun = mean, geom = "col", mapping = aes(group = .data[["split.by"]]),
        position = position_dodge(width = 0.9), width = 0.8, color = "black"
      ) +
        stat_summary(
          fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9), width = 0.2, color = "black"
        )
      y_min_use <- layer_scales(p)$y$range$range[1]
    }
    if (plot_type == "dot") {
      bins <- cut(dat$value, breaks = seq(min(dat$value), max(dat$value), length.out = 15), include.lowest = TRUE)
      bins_median <- sapply(strsplit(levels(bins), ","), function(x) median(as.numeric(gsub("\\(|\\)|\\[|\\]", "", x)), na.rm = TRUE))
      names(bins_median) <- levels(bins)
      dat[["bins"]] <- bins_median[bins]
      p <- p + geom_count(data = dat, aes(y = bins), shape = 21, alpha = alpha, position = position_dodge(width = 0.9)) +
        scale_size_area(name = "Count", max_size = 6, n.breaks = 4) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 2))
    }
    if (plot_type == "col") {
      p <- p + geom_col()
      if (flip) {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      } else {
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (isFALSE(individual) && isTRUE(nlevels(dat[["group.by"]]) > 1)) {
        x_index <- split(dat[["cell"]], dat[["group.by"]])
        border_data <- as.data.frame(sapply(x_index, min) - 0.5)
        colnames(border_data) <- "xintercept"
        border_data <- border_data[2:nrow(border_data), , drop = FALSE]
        border_layer <- geom_vline(xintercept = border_data[["xintercept"]], linetype = 2)
        p <- p + border_layer
      }
    }

    if (length(comparisons) > 0) {
      if (isTRUE(comparisons)) {
        group_use <- names(which(rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 2))
        if (any(rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 3)) {
          warning("Detected more than 2 groups. Use multiple_method for comparison", immediate. = TRUE)
          method <- multiple_method
        } else {
          method <- pairwise_method
        }
        if (sig_label == "p.format") {
          p <- p + stat_compare_means(
            data = dat[dat[["group.by"]] %in% group_use, , drop = FALSE],
            mapping = aes(x = .data[["group.by"]], y = .data[["value"]], label = after_stat(p.format)),
            size = 3.5,
            step.increase = 0.1,
            tip.length = 0.03,
            vjust = 1,
            method = method
          )
        } else {
          p <- p + stat_compare_means(
            data = dat[dat[["group.by"]] %in% group_use, , drop = FALSE],
            mapping = aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]], label = after_stat(p.signif)),
            size = 3.5,
            step.increase = 0.1,
            tip.length = 0.03,
            vjust = 1,
            method = method
          )
        }
        y_max_use <- layer_scales(p)$y$range$range[2]
      } else {
        p <- p + stat_compare_means(
          aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]]),
          label = sig_label,
          size = 3.5,
          step.increase = 0.1,
          tip.length = 0.03,
          vjust = 0,
          comparisons = comparisons, ref.group = ref_group, method = pairwise_method
        )
        y_max_use <- layer_scales(p)$y$range$range[1] + (layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]) * 1.15
      }
    }
    if (isTRUE(multiplegroup_comparisons)) {
      p <- p + stat_compare_means(
        aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]]),
        method = multiple_method,
        label = sig_label,
        label.y = Inf,
        size = 3.5,
        vjust = 1.2,
        hjust = 0
      )
      y_max_use <- layer_scales(p)$y$range$range[1] + (layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]) * 1.15
    }

    if (isTRUE(add_point)) {
      suppressWarnings(p <- p + geom_point(
        aes(x = .data[["group.by"]], y = .data[["value"]], linetype = rep(f, nrow(dat)), group = .data[["group.unique"]]),
        inherit.aes = FALSE,
        color = pt.color, size = pt.size, alpha = pt.alpha,
        position = position_jitterdodge(jitter.width = jitter.width, dodge.width = 0.9, seed = 11), show.legend = FALSE
      ))
      if (!is.null(cells.highlight)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight)
        if (nrow(cell_df) > 0) {
          p <- p + geom_point(
            data = cell_df, aes(x = .data[["group.by"]], y = .data[["value"]], linetype = rep(f, nrow(cell_df)), group = .data[["group.unique"]]), inherit.aes = FALSE,
            color = cols.highlight, size = sizes.highlight, alpha = alpha.highlight,
            position = position_jitterdodge(jitter.width = jitter.width, dodge.width = 0.9, seed = 11), show.legend = FALSE
          )
        }
      }
    }
    if (isTRUE(add_box)) {
      p <- p + geom_boxplot(aes(group = .data[["group.unique"]]),
        position = position_dodge(width = 0.9), color = box_color, fill = box_color, width = box_width, show.legend = FALSE, outlier.shape = NA
      ) +
        stat_summary(
          fun = median, geom = "point", mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9), color = "black", fill = "white", size = box_ptsize, shape = 21,
        )
    }
    if (isTRUE(add_trend)) {
      if (plot_type %in% c("violin", "box")) {
        if (nlevels(dat[["split.by"]]) > 1) {
          point_layer <- stat_summary(
            fun = median, geom = "point", mapping = aes(group = .data[["split.by"]], color = .data[["group.by"]]),
            position = position_dodge(width = 0.9), fill = "white", size = trend_ptsize, shape = 21
          )
          p_data <- p + point_layer
          p <- p + geom_line(
            data = layer_data(p_data, length(p_data$layers)),
            aes(x = x, y = y, group = colour),
            color = trend_color, linewidth = trend_linewidth, inherit.aes = FALSE
          ) +
            stat_summary(
              fun = median, geom = "point", mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9), color = "black", fill = "white", size = trend_ptsize, shape = 21
            )
        } else {
          p <- p + stat_summary(
            fun = median, geom = "line", mapping = aes(group = .data[["split.by"]]),
            position = position_dodge(width = 0.9), color = trend_color, linewidth = trend_linewidth
          ) +
            stat_summary(
              fun = median, geom = "point", mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9), color = "black", fill = "white", size = trend_ptsize, shape = 21
            )
        }
      }
      if (plot_type %in% c("bar")) {
        if (nlevels(dat[["split.by"]]) > 1) {
          point_layer <- stat_summary(
            fun = mean, geom = "point", mapping = aes(group = .data[["split.by"]], color = .data[["group.by"]]),
            position = position_dodge(width = 0.9), fill = "white", size = trend_ptsize, shape = 21
          )
          p_data <- p + point_layer
          p <- p + geom_line(
            data = layer_data(p_data, length(p_data$layers)),
            aes(x = x, y = y, group = colour),
            color = trend_color, linewidth = trend_linewidth, inherit.aes = FALSE
          ) +
            stat_summary(
              fun = mean, geom = "point", mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9), color = "black", fill = "white", size = trend_ptsize, shape = 21
            )
        } else {
          p <- p + stat_summary(
            fun = mean, geom = "line", mapping = aes(group = .data[["split.by"]]),
            position = position_dodge(width = 0.9), color = trend_color, linewidth = trend_linewidth,
          ) +
            stat_summary(
              fun = mean, geom = "point", mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9), color = "black", fill = "white", size = trend_ptsize, shape = 21
            )
        }
      }
    }
    if (add_stat != "none") {
      p <- p + stat_summary(
        fun = add_stat, geom = "point", mapping = aes(group = .data[["split.by"]], shape = 95),
        position = position_dodge(width = 0.9), color = stat_color, fill = "white", size = stat_size, stroke = 10,
      ) + scale_shape_identity()
    }

    if (nrow(dat) == 0) {
      p <- p + facet_null()
    } else {
      if (isTRUE(stack) && !isTRUE(flip)) {
        p <- p + facet_grid(features ~ .) + theme(strip.text.y = element_text(angle = 0))
      } else {
        p <- p + facet_grid(. ~ features)
      }
    }
    p <- p + labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
    if (nrow(dat) != 0) {
      p <- p + scale_x_discrete(drop = !keep_empty)
    }

    if (isTRUE(stack)) {
      p <- p + scale_y_continuous(
        limits = c(y_min_use, y_max_use), trans = y.trans,
        breaks = c(y_min_use, y_max_use), labels = c(round(y_min_use, 1), round(y_max_use, 1))
      )
    } else {
      p <- p + scale_y_continuous(limits = c(y_min_use, y_max_use), trans = y.trans, n.breaks = y.nbreaks)
    }

    if (fill.by != "expression") {
      if (isTRUE(stack)) {
        p <- p + scale_fill_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, limits = levels_order, drop = FALSE) +
          scale_color_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, limits = levels_order, drop = FALSE)
      } else {
        p <- p + scale_fill_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, drop = FALSE) +
          scale_color_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, drop = FALSE)
      }
      p <- p + guides(fill = guide_legend(
        title.hjust = 0,
        keywidth = 0.05,
        keyheight = 0.05,
        default.unit = "inch",
        order = 1,
        override.aes = list(size = 4.5, color = "black", alpha = 1)
      ))
    } else {
      p <- p + scale_fill_gradientn(
        name = paste0(keynm, ":"), colours = colors, limits = colors_limits
      ) + guides(
        fill = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)
      )
    }
    plist[[paste0(f, ":", g, ":", paste0(single_group, collapse = ","), ":", paste0(sp, collapse = ","))]] <- p
  })

  if (isTRUE(stack) && length(stat.by) > 1 && isFALSE(individual)) {
    for (g in group.by) {
      plist_g <- plist[sapply(strsplit(names(plist), ":"), function(x) x[2]) == g]
      legend <- get_legend(plist_g[[1]])
      if (isTRUE(flip)) {
        lab <- textGrob(label = ifelse(is.null(ylab), "Expression level", ylab), hjust = 0.5)
        plist_g <- lapply(seq_along(plist_g), FUN = function(i) {
          p <- plist_g[[i]]
          if (i != 1) {
            suppressWarnings(p <- p + theme(
              legend.position = "none",
              panel.grid = element_blank(),
              plot.title = element_blank(),
              plot.subtitle = element_blank(),
              axis.title = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(hjust = c(0, 1)),
              axis.ticks.length.y = unit(0, "pt"),
              plot.margin = unit(c(0, -0.5, 0, 0), "mm")
            ))
          } else {
            suppressWarnings(p <- p + theme(
              legend.position = "none",
              panel.grid = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(hjust = c(0, 1)),
              axis.ticks.length.y = unit(0, "pt"),
              plot.margin = unit(c(0, -0.5, 0, 0), "mm")
            ))
          }
          return(as_grob(p))
        })
        gtable <- do.call(cbind, plist_g)
        gtable <- add_grob(gtable, lab, "bottom", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      } else {
        lab <- textGrob(label = ifelse(is.null(ylab), "Expression level", ylab), rot = 90, hjust = 0.5)
        plist_g <- lapply(seq_along(plist_g), FUN = function(i) {
          p <- plist_g[[i]]
          if (i != length(plist_g)) {
            suppressWarnings(p <- p + theme(
              legend.position = "none",
              panel.grid = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(vjust = c(0, 1)),
              axis.ticks.length.x = unit(0, "pt"),
              plot.margin = unit(c(-0.5, 0, 0, 0), "mm")
            ))
            if (i == 1) {
              p <- p + theme(plot.title = element_blank(), plot.subtitle = element_blank())
            }
          } else {
            suppressWarnings(p <- p + theme(
              legend.position = "none",
              panel.grid = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(vjust = c(0, 1)),
              axis.ticks.length.x = unit(0, "pt"),
              plot.margin = unit(c(-0.5, 0, 0, 0), "mm")
            ))
          }
          return(as_grob(p))
        })
        gtable <- do.call(rbind, plist_g)
        gtable <- add_grob(gtable, lab, "left", clip = "off")
        gtable <- add_grob(gtable, legend, legend.position)
      }
      gtable <- gtable_add_padding(gtable, unit(c(1, 1, 1, 1), units = "cm"))
      plot <- wrap_plots(gtable)
      plist_stack[[g]] <- plot
    }
  }

  if (length(plist_stack) > 0) {
    plist <- plist_stack
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' Statistical plot of cells
#'
#' @param srt A Seurat object.
#' @param stat.by The name of a metadata column to be counted.
#' @param group.by
#' @param split.by
#' @param cells
#' @param keep_empty
#' @param individual
#' @param plot_type
#' @param stat_type
#' @param position
#' @param palette
#' @param palcolor
#' @param alpha
#' @param label
#' @param label.size
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param force
#' @param flip
#' @param NA_color
#' @param NA_stat
#' @param stat_level
#' @param bg.by
#' @param bg_palette
#' @param bg_palcolor
#' @param bg_apha
#' @param theme_args
#' @param seed
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "SubCellType", label = TRUE)
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "SubCellType", label = TRUE) %>% panel_fix(height = 2, width = 3)
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "SubCellType", stat_type = "count", position = "dodge", label = TRUE)
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "SubCellType", bg.by = "CellType", palette = "Set1", stat_type = "count", position = "dodge")
#'
#' CellStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "bar")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "rose")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "ring")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "pie")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "dot")
#'
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "bar")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "rose")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "ring")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "area")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "dot")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "trend")
#'
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "bar", individual = TRUE)
#'
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "bar")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "rose")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "ring")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "area")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "dot")
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "trend")
#'
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "bar", position = "dodge", label = TRUE)
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "rose", position = "dodge", label = TRUE)
#' CellStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "ring", position = "dodge", label = TRUE)
#'
#' CellStatPlot(pancreas_sub, stat.by = c("CellType", "Phase"), plot_type = "sankey")
#' CellStatPlot(pancreas_sub, stat.by = c("CellType", "Phase"), plot_type = "chord")
#'
#' CellStatPlot(pancreas_sub,
#'   stat.by = c("CellType", "Phase"), plot_type = "venn",
#'   stat_level = list(CellType = c("Ductal", "Ngn3 low EP"), Phase = "S")
#' )
#' pancreas_sub$Progenitor <- pancreas_sub$CellType %in% c("Ngn3 low EP", "Ngn3 high EP")
#' pancreas_sub$G2M <- pancreas_sub$Phase == "G2M"
#' pancreas_sub$Sox9_Expressed <- pancreas_sub[["RNA"]]@counts["Sox9", ] > 0
#' pancreas_sub$Neurog3_Expressed <- pancreas_sub[["RNA"]]@counts["Neurog3", ] > 0
#' CellStatPlot(pancreas_sub, stat.by = c("Progenitor", "G2M", "Sox9_Expressed", "Neurog3_Expressed"), plot_type = "venn", stat_level = "TRUE")
#' CellStatPlot(pancreas_sub, stat.by = c("Progenitor", "G2M", "Sox9_Expressed", "Neurog3_Expressed"), plot_type = "upset", stat_level = "TRUE")
#' sum(pancreas_sub$Progenitor == "FALSE" &
#'   pancreas_sub$G2M == "FALSE" &
#'   pancreas_sub$Sox9_Expressed == "TRUE" &
#'   pancreas_sub$Neurog3_Expressed == "FALSE")
#' @importFrom Seurat Cells
#' @export
CellStatPlot <- function(srt, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, cells = NULL, flip = FALSE,
                         NA_color = "grey", NA_stat = TRUE, keep_empty = FALSE, individual = FALSE, stat_level = NULL,
                         plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", "upset"),
                         stat_type = c("percent", "count"), position = c("stack", "dodge"),
                         palette = "Paired", palcolor = NULL, alpha = 1,
                         bg_palette = "Paired", bg_palcolor = NULL, bg_apha = 0.2,
                         label = FALSE, label.size = 3.5, label.fg = "black", label.bg = "white", label.bg.r = 0.1,
                         aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                         legend.position = "right", legend.direction = "vertical",
                         theme_use = "theme_scp", theme_args = list(),
                         combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  cells <- cells %||% colnames(srt@assays[[1]])
  meta.data <- srt@meta.data[cells, , drop = FALSE]

  plot <- StatPlot(
    meta.data = meta.data, stat.by = stat.by, group.by = group.by, split.by = split.by, bg.by = bg.by, flip = flip,
    NA_color = NA_color, NA_stat = NA_stat, keep_empty = keep_empty, individual = individual, stat_level = stat_level,
    plot_type = plot_type, stat_type = stat_type, position = position,
    palette = palette, palcolor = palcolor, alpha = alpha,
    bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_apha = bg_apha,
    label = label, label.size = label.size, label.fg = label.fg, label.bg = label.bg, label.bg.r = label.bg.r,
    aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
    legend.position = legend.position, legend.direction = legend.direction,
    theme_use = theme_use, theme_args = theme_args,
    combine = combine, nrow = nrow, ncol = ncol, byrow = byrow, force = force, seed = seed
  )
  return(plot)
}

#' StatPlot
#'
#' @param meta.data
#'
#' @param stat.by
#' @param group.by
#' @param split.by
#' @param flip
#' @param NA_color
#' @param NA_stat
#' @param keep_empty
#' @param individual
#' @param stat_level
#' @param plot_type
#' @param stat_type
#' @param position
#' @param palette
#' @param palcolor
#' @param alpha
#' @param label
#' @param label.size
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param force
#' @param bg.by
#' @param bg_palette
#' @param bg_palcolor
#' @param bg_apha
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' head(pancreas_sub@meta.data)
#' StatPlot(pancreas_sub@meta.data, stat.by = "Phase", group.by = "CellType", plot_type = "bar", label = TRUE)
#'
#' head(pancreas_sub[["RNA"]]@meta.features)
#' StatPlot(pancreas_sub[["RNA"]]@meta.features, stat.by = "highly_variable_genes", plot_type = "ring", label = TRUE)
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", IDtype = "symbol", db = "GeneType")
#' head(pancreas_sub[["RNA"]]@meta.features)
#' StatPlot(pancreas_sub[["RNA"]]@meta.features,
#'   stat.by = "highly_variable_genes", group.by = "GeneType",
#'   stat_type = "count", plot_type = "bar", position = "dodge", label = TRUE, NA_stat = FALSE
#' )
#'
#' @importFrom dplyr group_by across all_of mutate "%>%" .data summarise
#' @importFrom stats quantile xtabs
#' @importFrom ggplot2 ggplot aes labs position_identity position_stack position_dodge2 scale_x_continuous scale_y_continuous geom_col geom_area geom_vline scale_fill_manual scale_fill_identity scale_color_identity scale_fill_gradientn guides guide_legend element_line coord_polar annotate geom_sf theme_void after_stat scale_size_area
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom ggrepel geom_text_repel
#' @importFrom circlize chordDiagram circos.clear
#' @importFrom patchwork wrap_plots
#' @importFrom gtable gtable_add_rows gtable_add_cols gtable_add_grob
#' @importFrom rlang %||%
#' @export
StatPlot <- function(meta.data, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, flip = FALSE,
                     NA_color = "grey", NA_stat = TRUE, keep_empty = FALSE, individual = FALSE, stat_level = NULL,
                     plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", "upset"),
                     stat_type = c("percent", "count"), position = c("stack", "dodge"),
                     palette = "Paired", palcolor = NULL, alpha = 1,
                     bg_palette = "Paired", bg_palcolor = NULL, bg_apha = 0.2,
                     label = FALSE, label.size = 3.5, label.fg = "black", label.bg = "white", label.bg.r = 0.1,
                     aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                     legend.position = "right", legend.direction = "vertical",
                     theme_use = "theme_scp", theme_args = list(),
                     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  stat_type <- match.arg(stat_type)
  plot_type <- match.arg(plot_type)
  position <- match.arg(position)

  if (nrow(meta.data) == 0) {
    stop("meta.data is empty.")
  }
  if (is.null(group.by)) {
    group.by <- "No.group.by"
    xlab <- ""
    meta.data[[group.by]] <- factor("All")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    meta.data[[split.by]] <- factor("")
  }

  for (i in unique(c(group.by, split.by, bg.by))) {
    if (!i %in% colnames(meta.data)) {
      stop(paste0(i, " is not in the meta.data."))
    }
    if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }
  bg_map <- NULL
  if (!is.null(bg.by)) {
    for (g in group.by) {
      df_table <- table(meta.data[[g]], meta.data[[bg.by]])
      if (max(rowSums(df_table > 0), na.rm = TRUE) > 1) {
        stop("'group.by' must be a part of 'bg.by'")
      } else {
        bg_map[[g]] <- setNames(colnames(df_table)[apply(df_table, 1, function(x) which(x > 0))], rownames(df_table))
      }
    }
  } else {
    for (g in group.by) {
      bg_map[[g]] <- setNames(levels(meta.data[[g]]), levels(meta.data[[g]]))
    }
  }
  for (i in unique(stat.by)) {
    if (!i %in% colnames(meta.data)) {
      stop(paste0(i, " is not in the meta.data."))
    }
    if (plot_type %in% c("venn", "upset")) {
      if (!is.factor(meta.data[[i]]) && !is.logical(meta.data[[i]])) {
        meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
      }
    } else if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }

  if (length(stat.by) >= 2) {
    if (!plot_type %in% c("sankey", "chord", "venn", "upset")) {
      stop("plot_type must be one of 'sankey', 'chord', 'venn' and 'upset' whtn multiple 'stat.by' provided.")
    }
    if (length(stat.by) > 2 && plot_type == "chord") {
      stop("'stat.by' can only be a vector of length 2 when 'plot_type' is 'chord'.")
    }
    if (length(stat.by) > 7 && plot_type == "venn") {
      stop("'stat.by' can only be a vector of length <= 7 when 'plot_type' is 'venn'.")
    }
  }
  levels <- unique(unlist(lapply(meta.data[, stat.by, drop = FALSE], function(x) {
    if (is.factor(x)) {
      return(levels(x))
    }
    if (is.logical(x)) {
      return(as.character(unique(x)))
    }
  })))

  if (plot_type %in% c("venn", "upset")) {
    if (is.null(stat_level)) {
      stat_level <- lapply(stat.by, function(stat) {
        levels(meta.data[[stat]])[1] %||% sort(unique(meta.data[[stat]]))[1]
      })
      message("stat_level is set to ", paste0(stat_level, collapse = ","))
    } else {
      if (length(stat_level) == 1) {
        stat_level <- rep(stat_level, length(stat.by))
      }
      if (length(stat_level) != length(stat.by)) {
        stop("'stat_level' must be of length 1 or the same length as 'stat.by'")
      }
    }
    if (is.null(names(stat_level))) {
      names(stat_level) <- stat.by
    }
    for (i in stat.by) {
      meta.data[[i]] <- meta.data[[i]] %in% stat_level[[i]]
    }
  }

  if (plot_type %in% c("rose", "ring", "pie")) {
    aspect.ratio <- 1
  }

  if (any(group.by != "No.group.by") && plot_type %in% c("sankey", "chord", "venn", "upset")) {
    warning("group.by is not used when plot sankey, chord, venn or upset", immediate. = TRUE)
  }
  dat_all <- meta.data[, unique(c(stat.by, group.by, split.by, bg.by)), drop = FALSE]
  nlev <- sapply(dat_all, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 100 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  dat_split <- split.data.frame(dat_all, dat_all[[split.by]])

  plist <- list()
  if (plot_type %in% c("bar", "rose", "ring", "pie", "trend", "area", "dot")) {
    xlab <- xlab %||% group.by
    ylab <- ylab %||% ifelse(stat_type == "count", "Count", "Percentage")
    if (identical(theme_use, "theme_blank")) {
      theme_args[["xlab"]] <- xlab
      theme_args[["ylab"]] <- ylab
      if (plot_type %in% c("rose", "ring", "pie")) {
        theme_args[["add_coord"]] <- FALSE
      }
    }
    colors <- palette_scp(dat_all[[stat.by]], palette = palette, palcolor = palcolor, NA_color = NA_color, NA_keep = TRUE)

    comb_list <- list()
    comb <- expand.grid(stat_name = stat.by, group_name = group.by, stringsAsFactors = FALSE)
    if (isTRUE(individual)) {
      for (g in group.by) {
        comb_list[[g]] <- merge(comb, expand.grid(
          group_name = g, group_element = levels(dat_all[[g]]),
          split_name = levels(dat_all[[split.by]]), stringsAsFactors = FALSE
        ),
        by = "group_name"
        )
      }
    } else {
      for (g in group.by) {
        comb_list[[g]] <- merge(comb, expand.grid(
          group_name = g, group_element = list(levels(dat_all[[g]])),
          split_name = levels(dat_all[[split.by]]), stringsAsFactors = FALSE
        ),
        by = "group_name"
        )
      }
    }
    comb <- do.call(rbind, comb_list)
    rownames(comb) <- paste0(
      comb[["group_name"]], ":",
      sapply(comb[["group_element"]], function(x) paste0(x, collapse = ",")), ":",
      comb[["split_name"]]
    )

    plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
      stat.by <- comb[i, "stat_name"]
      sp <- comb[i, "split_name"]
      g <- comb[i, "group_name"]
      single_group <- comb[[i, "group_element"]]
      colors_use <- colors[names(colors) %in% dat_split[[ifelse(split.by == "All_cells", 1, sp)]][[stat.by]]]
      if (any(is.na(dat_split[[ifelse(split.by == "All_cells", 1, sp)]][[stat.by]])) && isTRUE(NA_stat)) {
        colors_use <- c(colors_use, colors["NA"])
      }
      if (stat_type == "percent") {
        dat_use <- dat_split[[ifelse(split.by == "All_cells", 1, sp)]] %>%
          xtabs(formula = paste0("~", stat.by, "+", g), addNA = NA_stat) %>%
          as.data.frame() %>%
          group_by(across(all_of(g)), .drop = FALSE) %>%
          mutate(groupn = sum(Freq)) %>%
          group_by(across(all_of(c(stat.by, g))), .drop = FALSE) %>%
          mutate(value = Freq / groupn) %>%
          as.data.frame()
      } else {
        dat_use <- dat_split[[ifelse(split.by == "All_cells", 1, sp)]] %>%
          xtabs(formula = paste0("~", stat.by, "+", g), addNA = NA_stat) %>%
          as.data.frame() %>%
          mutate(value = Freq)
      }
      dat <- dat_use[dat_use[[g]] %in% single_group, , drop = FALSE]
      dat[[g]] <- factor(dat[[g]], levels = levels(dat[[g]])[levels(dat[[g]]) %in% dat[[g]]])
      dat <- dat[!is.na(dat[["value"]]), , drop = FALSE]
      if (!is.null(bg.by)) {
        bg <- bg.by
        bg_color <- palette_scp(levels(dat_all[[bg]]), palette = bg_palette, palcolor = bg_palcolor)
      } else {
        bg <- g
        bg_color <- palette_scp(levels(dat_all[[bg]]), palcolor = bg_palcolor %||% rep(c("transparent", "grey85"), nlevels(dat_all[[bg]])))
      }

      if (isTRUE(flip)) {
        dat[[g]] <- factor(dat[[g]], levels = rev(levels(dat[[g]])))
        aspect.ratio <- 1 / aspect.ratio
        if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
          aspect.ratio <- NULL
        }
      }
      if (plot_type == "ring") {
        dat[[g]] <- factor(dat[[g]], levels = c("   ", levels(dat[[g]])))
        dat <- rbind(dat, dat[nrow(dat) + 1, , drop = FALSE])
        dat[nrow(dat), g] <- "   "
      }
      if (plot_type == "dot") {
        position_use <- position_identity()
        scalex <- scale_x_discrete(drop = !keep_empty)
      } else {
        if (position == "stack") {
          position_use <- position_stack(vjust = 0.5)
          scalex <- scale_x_discrete(drop = !keep_empty, expand = c(0, 0))
          scaley <- scale_y_continuous(
            labels = if (stat_type == "count") scales::number else scales::percent,
            expand = c(0, 0)
          )
        } else if (position == "dodge") {
          if (plot_type == "area") {
            position_use <- position_dodge2(width = 0.9, preserve = "total")
          } else {
            position_use <- position_dodge2(width = 0.9, preserve = "single")
          }
          scalex <- scale_x_discrete(drop = !keep_empty)
          scaley <- scale_y_continuous(
            limits = c(0, max(dat[["value"]], na.rm = TRUE) * 1.1),
            labels = if (stat_type == "count") scales::number else scales::percent,
            expand = c(0, 0)
          )
        }
      }
      if (position == "stack") {
        bg_layer <- NULL
      } else {
        bg_data <- na.omit(unique(dat[, g, drop = FALSE]))
        bg_data[["x"]] <- as.numeric(bg_data[[g]])
        bg_data[["xmin"]] <- ifelse(bg_data[["x"]] == min(bg_data[["x"]]), -Inf, bg_data[["x"]] - 0.5)
        bg_data[["xmax"]] <- ifelse(bg_data[["x"]] == max(bg_data[["x"]]), Inf, bg_data[["x"]] + 0.5)
        bg_data[["ymin"]] <- -Inf
        bg_data[["ymax"]] <- Inf
        bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[[g]])]]
        bg_layer <- geom_rect(data = bg_data, xmin = bg_data[["xmin"]], xmax = bg_data[["xmax"]], ymin = bg_data[["ymin"]], ymax = bg_data[["ymax"]], fill = bg_data[["fill"]], alpha = bg_apha, inherit.aes = FALSE)
      }

      if (plot_type == "bar") {
        p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
          bg_layer +
          geom_col(aes(fill = .data[[stat.by]]),
            width = 0.8,
            color = "black",
            alpha = alpha,
            position = position_use
          ) +
          scalex +
          scaley
      }
      if (plot_type == "trend") {
        dat_area <- dat[rep(seq_len(nrow(dat)), each = 2), , drop = FALSE]
        dat_area[[g]] <- as.numeric(dat_area[[g]])
        dat_area[seq(1, nrow(dat_area), 2), g] <- dat_area[seq(1, nrow(dat_area), 2), g] - 0.3
        dat_area[seq(2, nrow(dat_area), 2), g] <- dat_area[seq(2, nrow(dat_area), 2), g] + 0.3
        p <- ggplot(dat, aes(x = .data[[g]], y = value, fill = .data[[stat.by]])) +
          bg_layer +
          geom_area(
            data = dat_area, mapping = aes(x = .data[[g]], fill = .data[[stat.by]]),
            alpha = alpha / 2, color = "black", linetype = 2, position = position_use
          ) +
          geom_col(aes(fill = .data[[stat.by]]),
            width = 0.6,
            color = "black",
            alpha = alpha,
            position = position_use
          ) +
          scalex +
          scaley
      }
      if (plot_type == "rose") {
        p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
          bg_layer +
          geom_col(aes(fill = .data[[stat.by]]),
            width = 0.8,
            color = "black",
            alpha = alpha,
            position = position_use
          ) +
          scalex +
          scaley +
          coord_polar(theta = "x", start = ifelse(flip, pi / 2, 0))
      }
      if (plot_type == "ring" || plot_type == "pie") {
        p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
          bg_layer +
          geom_col(aes(fill = .data[[stat.by]]),
            width = 0.8,
            color = "black",
            alpha = alpha,
            position = position_use
          ) +
          scalex +
          scaley +
          coord_polar(theta = "y", start = ifelse(flip, pi / 2, 0))
      }
      if (plot_type == "area") {
        p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
          bg_layer +
          geom_area(aes(fill = .data[[stat.by]]),
            color = "black",
            alpha = alpha,
            position = position_use
          ) +
          scalex +
          scaley
      }
      if (plot_type == "dot") {
        p <- ggplot(dat, aes(x = .data[[g]], y = .data[[stat.by]])) +
          bg_layer +
          geom_point(aes(fill = .data[[stat.by]], size = value),
            color = "black",
            alpha = alpha,
            shape = 21,
            position = position_use
          ) +
          scalex +
          scale_size_area(name = capitalize(stat_type), max_size = 12) +
          guides(size = guide_legend(override.aes = list(fill = "grey30")))
      }
      if (isTRUE(label)) {
        if (plot_type == "dot") {
          p <- p + geom_text_repel(
            aes(
              x = .data[[g]], y = .data[[stat.by]],
              label = if (stat_type == "count") value else paste0(round(value * 100, 1), "%"),
            ),
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100, min.segment.length = 0,
            position = position_use
          )
        } else {
          p <- p + geom_text_repel(
            aes(
              label = if (stat_type == "count") value else paste0(round(value * 100, 1), "%"),
              y = value
            ),
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100, min.segment.length = 0,
            position = position_use
          )
        }
      }
      if (plot_type %in% c("rose")) {
        # angle <- 360 / (2 * pi) * rev(seq(pi / nlevels(dat[[g]]), 2 * pi - pi / nlevels(dat[[g]]), len = nlevels(dat[[g]])))
        # axis.text.x <- element_text(angle = angle)
        axis.text.x <- element_text()
      } else if (plot_type %in% c("ring", "pie")) {
        axis.text.x <- element_text()
      } else {
        axis.text.x <- element_text(angle = 45, hjust = 1, vjust = 1)
      }
      p <- p + labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
        scale_fill_manual(
          name = paste0(stat.by, ":"), values = colors_use, na.value = colors_use["NA"], drop = FALSE,
          limits = names(colors_use), na.translate = T
        ) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          axis.text.x = axis.text.x,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = if (plot_type == "trend" & stat_type == "percent") element_blank() else element_line(colour = "grey80", linetype = 2)
        ) + guides(fill = guide_legend(
          title.hjust = 0,
          keywidth = 0.05,
          keyheight = 0.05,
          default.unit = "inch",
          order = 1,
          override.aes = list(size = 4.5, color = "black", alpha = 1)
        ))
      if (isTRUE(flip) && !plot_type %in% c("pie", "rose")) {
        p <- p + coord_flip()
      }
      return(p)
    })
  } else if (plot_type %in% c("chord", "sankey", "venn", "upset")) {
    colors <- palette_scp(stat.by, palette = palette, palcolor = palcolor)
    if (plot_type == "chord" && isTRUE(combine)) {
      temp <- tempfile(fileext = "png")
      grDevices::png(temp)
      grDevices::dev.control("enable")
      nlev <- nlevels(dat_all[[split.by]])
      if (is.null(nrow) && is.null(ncol)) {
        nrow <- ceiling(sqrt(nlev))
        ncol <- ceiling(nlev / nrow)
      }
      if (is.null(nrow)) {
        nrow <- ceiling(sqrt(ncol))
      }
      if (is.null(ncol)) {
        ncol <- ceiling(sqrt(nrow))
      }
      par(mfrow = c(nrow, ncol))
    }
    for (sp in levels(dat_all[[split.by]])) {
      dat_use <- dat_split[[ifelse(split.by == "All_cells", 1, sp)]]
      if (plot_type == "venn") {
        check_R("ggVennDiagram")
        dat_list <- as.list(dat_use[, stat.by])
        dat_list <- lapply(setNames(names(dat_list), names(dat_list)), function(x) {
          lg <- dat_list[[x]]
          names(lg) <- rownames(dat_use)
          cellkeep <- names(lg)[lg]
          return(cellkeep)
        })
        venn <- ggVennDiagram::Venn(dat_list)
        data <- ggVennDiagram::process_data(venn)
        dat_venn_region <- ggVennDiagram::venn_region(data)
        idname <- dat_venn_region[["name"]][dat_venn_region[["name"]] %in% stat.by]
        names(idname) <- dat_venn_region[["id"]][dat_venn_region[["name"]] %in% stat.by]
        idcomb <- strsplit(dat_venn_region[["id"]], split = "")
        colorcomb <- lapply(idcomb, function(x) colors[idname[as.character(x)]])
        dat_venn_region[["colors"]] <- sapply(colorcomb, function(x) blendcolors(x, mode = "blend"))
        dat_venn_region[["label"]] <- paste0(
          dat_venn_region[["count"]], "\n",
          round(dat_venn_region[["count"]] / sum(dat_venn_region[["count"]]) * 100, 1), "%"
        )
        dat_venn_setedge <- ggVennDiagram::venn_setedge(data)
        dat_venn_setedge[["colors"]] <- colors[dat_venn_setedge[["name"]]]
        dat_venn_setlabel <- ggVennDiagram::venn_setlabel(data)
        dat_title <- as.data.frame(do.call(rbind, dat_venn_setlabel$geometry))
        colnames(dat_title) <- c("x", "y")
        dat_title[["label"]] <- paste0(
          dat_venn_setlabel[["name"]], "\n",
          "(", sapply(dat_list, length)[dat_venn_setlabel[["name"]]], ")"
        )
        dat_stat <- as.data.frame(do.call(rbind, lapply(dat_venn_region$geometry, function(x) sf::st_centroid(x))))
        colnames(dat_stat) <- c("x", "y")
        dat_stat[["label"]] <- dat_venn_region[["label"]]
        p <- ggplot() +
          geom_sf(data = dat_venn_region, aes(fill = colors), alpha = alpha) +
          geom_sf(data = dat_venn_setedge, aes(color = colors), size = 1) +
          geom_text_repel(
            data = dat_title, aes(label = label, x = x, y = y),
            fontface = "bold",
            colour = label.fg, size = label.size + 0.5,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100,
            min.segment.length = 0, segment.colour = "black"
          ) +
          geom_text_repel(
            data = dat_stat, aes(label = label, x = x, y = y),
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100,
            min.segment.length = 0, segment.colour = "black"
          ) +
          scale_fill_identity() +
          scale_color_identity() +
          theme(
            plot.title = element_text(hjust = 0.5),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()
          )
        p <- p + labs(x = sp, title = title, subtitle = subtitle)
      }
      if (plot_type == "upset") {
        check_R("ggupset")
        for (n in seq_len(nrow(dat_use))) {
          dat_use[["intersection"]][n] <- list(stat.by[unlist(dat_use[n, stat.by])])
        }
        dat_use <- dat_use[sapply(dat_use[["intersection"]], length) > 0, , drop = FALSE]
        p <- ggplot(dat_use, aes(x = intersection)) +
          geom_bar(aes(fill = after_stat(count)), color = "black", width = 0.5, show.legend = FALSE) +
          geom_text_repel(aes(label = after_stat(count)),
            stat = "count",
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100,
            min.segment.length = 0, segment.colour = "black"
          ) +
          labs(title = title, subtitle = subtitle, x = sp, y = "Intersection size") +
          ggupset::scale_x_upset(sets = stat.by, n_intersections = 20) +
          scale_fill_gradientn(colors = palette_scp(palette = "material-indigo")) +
          theme_scp(
            aspect.ratio = 0.6,
            panel.grid.major = element_line(colour = "grey80", linetype = 2)
          ) +
          ggupset::theme_combmatrix(
            combmatrix.label.text = element_text(size = 12, color = "black"),
            combmatrix.label.extra_spacing = 6
          )
        p <- p + labs(title = title, subtitle = subtitle)
      }
      if (plot_type == "sankey") {
        colors <- palette_scp(c(unique(unlist(lapply(dat_all[, stat.by, drop = FALSE], levels))), NA), palette = palette, palcolor = palcolor, NA_keep = TRUE, NA_color = NA_color)
        legend_list <- list()
        for (l in stat.by) {
          df <- data.frame(levels(dat_use[[l]]))
          colnames(df) <- l
          legend_list[[l]] <- get_legend(ggplot(data = df) +
            geom_col(aes(x = 1, y = 1, fill = .data[[l]]), color = "black") +
            scale_fill_manual(values = colors[levels(dat_use[[l]])]) +
            guides(fill = guide_legend(
              title.hjust = 0,
              title.vjust = 0,
              keywidth = 0.05,
              keyheight = 0.05,
              default.unit = "inch",
              order = 1,
              override.aes = list(size = 4.5, color = "black", alpha = 1)
            )) +
            theme_scp(
              legend.position = "bottom",
              legend.direction = legend.direction
            ))
          if (any(is.na(dat_use[[l]]))) {
            raw_levels <- levels(dat_use[[l]])
            dat_use[[l]] <- as.character(dat_use[[l]])
            dat_use[[l]][is.na(dat_use[[l]])] <- "NA"
            dat_use[[l]] <- factor(dat_use[[l]], levels = c(raw_levels, "NA"))
          }
        }
        if (legend.direction == "vertical") {
          legend <- do.call(cbind, legend_list)
        } else {
          legend <- do.call(rbind, legend_list)
        }
        dat <- suppressWarnings(make_long(dat_use, all_of(stat.by)))
        dat$node <- factor(dat$node, levels = rev(names(colors)))
        p0 <- ggplot(dat, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node)) +
          geom_sankey(color = "black", flow.alpha = alpha, show.legend = FALSE, na.rm = FALSE) +
          scale_fill_manual(values = colors, drop = FALSE) +
          scale_x_discrete(expand = c(0, 0.2)) +
          theme_void() +
          theme(axis.text.x = element_text())
        gtable <- as_grob(p0)
        gtable <- add_grob(gtable, legend, legend.position)
        p <- wrap_plots(gtable)
      }
      if (plot_type == "chord") {
        colors <- palette_scp(c(unique(unlist(lapply(dat_all[, stat.by, drop = FALSE], levels))), NA), palette = palette, palcolor = palcolor, NA_keep = TRUE, NA_color = NA_color)
        M <- table(dat_use[[stat.by[1]]], dat_use[[stat.by[2]]], useNA = "ifany")
        m <- matrix(M, ncol = ncol(M), dimnames = dimnames(M))
        colnames(m)[is.na(colnames(m))] <- "NA"
        chordDiagram(m,
          grid.col = colors,
          transparency = 0.2,
          link.lwd = 1,
          link.lty = 1,
          link.border = 1
        )
        circos.clear()
        p <- grDevices::recordPlot()

        # library(grid)
        # library(gridBase)
        # plot.new()
        # pushViewport(
        #   viewport(x = 0.5, y = 0.5, width = unit(1, "snpc"), height = unit(1, "snpc"), just = c("left", "center"))
        # )
        # par(omi = gridOMI(), new = TRUE)
        # chord()
      }

      plist[[sp]] <- p
    }
  }
  if (isTRUE(combine) && plot_type == "chord") {
    plot <- grDevices::recordPlot()
    grDevices::dev.off()
    unlink(temp)
    return(plot)
  }
  if (isTRUE(combine) && plot_type != "chord") {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' Features correlation plot
#'
#' @param srt
#' @param features
#' @param group.by
#' @param split.by
#' @param slot
#' @param assay
#' @param cor_method
#' @param adjust
#' @param margin
#' @param reverse
#' @param add_equation
#' @param add_r2
#' @param add_pvalue
#' @param add_smooth
#' @param palette
#' @param palcolor
#' @param bg_color
#' @param pt.size
#' @param pt.alpha
#' @param cells.highlight
#' @param cols.highlight
#' @param sizes.highlight
#' @param alpha.highlight
#' @param stroke.highlight
#' @param calculate_coexp
#' @param raster
#' @param raster.dpi
#' @param theme_use
#' @param title
#' @param subtitle
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param force
#' @param cells
#' @param aspect.ratio
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' FeatureCorPlot(pancreas_sub, features = c("Ghrl", "Gcg", "Ins1", "Ins2"), group.by = "SubCellType")
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom SeuratObject as.sparse
#' @importFrom dplyr group_by "%>%" .data
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_density_2d stat_density_2d labs scale_x_continuous scale_y_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colourbar geom_hex stat_summary_hex
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom ggrepel geom_text_repel
#' @importFrom gtable gtable_add_cols
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @importFrom methods slot
#' @export
FeatureCorPlot <- function(srt, features, group.by = NULL, split.by = NULL, cells = NULL, slot = "data", assay = NULL,
                           cor_method = "pearson", adjust = 1, margin = 1, reverse = FALSE,
                           add_equation = FALSE, add_r2 = TRUE, add_pvalue = TRUE, add_smooth = TRUE,
                           palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = NULL, pt.alpha = 1,
                           cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                           calculate_coexp = FALSE,
                           raster = NULL, raster.dpi = c(512, 512),
                           aspect.ratio = 1, title = NULL, subtitle = NULL,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  require("ggrepel", quietly = TRUE)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt@meta.data[[split.by]] <- factor("")
  }
  if (is.null(group.by)) {
    group.by <- "All_cells"
    srt@meta.data[[group.by]] <- factor("All_cells")
  }
  for (i in c(split.by, group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
  }
  if (!is.null(cells.highlight) & !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt@assays[[1]]))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt@assays[[1]]))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt@assays[[1]]))
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- colnames(srt@assays[[1]])
  }

  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  status <- check_DataType(srt, slot = slot, assay = assay)
  if (slot == "counts" && status != "raw_counts") {
    stop("Data in the 'counts' slot is not raw counts.")
  }
  if (slot == "data" && status != "log_normalized_counts") {
    if (status == "raw_counts") {
      warning("Data in the 'data' slot is raw counts. Perform NormalizeData(LogNormalize) on the data.", immediate. = TRUE)
      srt <- suppressWarnings(NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE))
    }
    if (status == "raw_normalized_counts") {
      warning("Data in the 'data' slot is raw_normalized_counts. Perform NormalizeData(LogNormalize) on the data.", immediate. = TRUE)
      srt <- suppressWarnings(NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE))
    }
    if (status == "unknown") {
      stop("Data in the 'data' slot is unknown. Please check the data type.")
    }
  }

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression", immediate. = TRUE)
    }
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      srt@meta.data[["CoExp"]] <- apply(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE], 2, function(x) exp(mean(log(x))))
    } else if (status == "log_normalized_counts") {
      srt@meta.data[["CoExp"]] <- apply(expm1(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE]), 2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene) > 0) {
    dat_gene <- t(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE])
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as.matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])
  if (length(features) < 2) {
    stop("features must be a vector of length at least 2.")
  }

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  if (!inherits(dat_exp, "dgCMatrix")) {
    dat_exp <- as.sparse(as.matrix(dat_exp))
  }
  if (length(features) > 10 && !isTRUE(force)) {
    warning("More than 10 features to be paired compared which will generate more than 50 plots.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  dat_use <- srt@meta.data[, unique(c(split.by, group.by)), drop = FALSE]
  dat_use <- cbind(dat_use, dat_exp[row.names(dat_use), , drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) * ncol(combn(features, m = 2)) > 1e5)
  if (isTRUE(raster)) {
    check_R("exaexa/scattermore")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector")
    }
  }

  plist <- list()
  colors <- palette_scp(levels(dat_use[[group.by]]), palette = palette, palcolor = palcolor)
  cor_palette <- palette_scp(x = seq(-1, 1, length.out = 200), palette = "RdBu")
  bound <- strsplit(gsub("\\(|\\)|\\[|\\]", "", names(cor_palette)), ",")
  bound <- lapply(bound, as.numeric)
  df_bound <- do.call(rbind, bound)
  rownames(df_bound) <- cor_palette
  df_bound[1, 1] <- df_bound[1, 1] - 0.01

  pair <- as.data.frame(t(combn(features, m = 2)))
  colnames(pair) <- c("feature1", "feature2")
  pair_expand <- expand.grid(features, features, stringsAsFactors = TRUE)
  colnames(pair_expand) <- c("feature1", "feature2")
  pair_expand[["feature1"]] <- factor(pair_expand[["feature1"]], levels = levels(pair_expand[["feature2"]]))

  for (s in levels(dat_use[[split.by]])) {
    dat <- dat_use[dat_use[[split.by]] == s, , drop = FALSE]
    feature_mat <- t(dat_exp[rownames(dat), features])
    if (cor_method %in% c("pearson", "spearman")) {
      if (cor_method == "spearman") {
        feature_mat <- t(apply(feature_mat, 1, rank))
      }
      cor_method <- "correlation"
    }
    pair_sim <- proxyC::simil(
      x = feature_mat,
      method = cor_method
    )
    if (isTRUE(reverse)) {
      order1 <- rev(pair_expand[, 1])
      order2 <- rev(pair_expand[, 2])
      levels(order1) <- rev(levels(order1))
      levels(order2) <- rev(levels(order2))
    } else {
      order1 <- pair_expand[, 1]
      order2 <- pair_expand[, 2]
    }
    plotlist <- mapply(FUN = function(x, y) {
      f1 <- as.character(x)
      f2 <- as.character(y)
      f1_index <- as.numeric(x)
      f2_index <- as.numeric(y)
      p <- ggplot(data = dat) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(margin, margin, margin, margin),
          legend.position = "none"
        )
      if (f1_index == f2_index) {
        p <- p + geom_violin(aes(x = .data[[group.by]], y = .data[[f1]], fill = .data[[group.by]]),
          scale = "width", adjust = adjust, trim = TRUE, na.rm = TRUE
        ) + scale_x_discrete(position = ifelse(isTRUE(reverse), "top", "bottom")) +
          scale_y_continuous(position = ifelse(isTRUE(reverse), "right", "left"))
      } else {
        p <- p + scale_x_continuous(
          n.breaks = 3, labels = scales::number_format(),
          limits = c(min(dat_exp[rownames(dat), ], na.rm = TRUE), max(dat_exp[rownames(dat), ], na.rm = TRUE)),
          position = ifelse(isTRUE(reverse), "top", "bottom")
        ) +
          scale_y_continuous(
            n.breaks = 3, labels = scales::number_format(),
            limits = c(min(dat_exp[rownames(dat), ], na.rm = TRUE), max(dat_exp[rownames(dat), ], na.rm = TRUE)),
            position = ifelse(isTRUE(reverse), "right", "left")
          )
      }
      if (f1_index < f2_index) {
        if (isTRUE(raster)) {
          p <- p + scattermore::geom_scattermore(
            mapping = aes(x = .data[[f1]], y = .data[[f2]], color = .data[[group.by]]),
            pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
          )
        } else {
          p <- p + geom_point(aes(x = .data[[f1]], y = .data[[f2]], color = .data[[group.by]]),
            alpha = pt.alpha, size = pt.size
          )
        }
        if (isTRUE(add_smooth)) {
          p <- p + geom_smooth(aes(x = .data[[f1]], y = .data[[f2]]),
            alpha = 0.5, method = "lm", color = "red", formula = y ~ x, na.rm = TRUE
          )
        }
        if (any(isTRUE(add_equation), isTRUE(add_r2), isTRUE(add_pvalue))) {
          m <- lm(dat[[f2]] ~ dat[[f1]])
          if (coef(m)[2] >= 0) {
            eq1 <- substitute(
              italic(y) == a + b %.% italic(x),
              list(
                a = format(as.numeric(coef(m)[1]), digits = 2),
                b = format(as.numeric(coef(m)[2]), digits = 2)
              )
            )
          } else {
            eq1 <- substitute(
              italic(y) == a - b %.% italic(x),
              list(
                a = format(as.numeric(coef(m)[1]), digits = 2),
                b = format(-as.numeric(coef(m)[2]), digits = 2)
              )
            )
          }
          eq1 <- as.character(as.expression(eq1))
          eq2 <- substitute(
            italic(r)^2 ~ "=" ~ r2,
            list(
              r2 = format(summary(m)$r.squared, digits = 2)
            )
          )
          eq2 <- as.character(as.expression(eq2))
          eq3 <- substitute(
            italic(p) ~ "=" ~ pvalue,
            list(
              pvalue = format(summary(m)$coefficients[2, 4], digits = 2)
            )
          )
          eq3 <- as.character(as.expression(eq3))
          eqs <- c(eq1, eq2, eq3)
          vjusts <- c(1.3, 1.3 * 2, 1.3 * 2^2)
          i <- c(isTRUE(add_equation), isTRUE(add_r2), isTRUE(add_pvalue))
          p <- p + annotate(
            geom = "text", x = -Inf, y = Inf, label = eqs[i], size = 3.5,
            hjust = -0.05, vjust = vjusts[1:sum(i)], parse = TRUE
          )
        }
        if (!is.null(cells.highlight)) {
          cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight_use)
          if (nrow(cell_df) > 0) {
            # point_size <- p$layers[[1]]$aes_params$size
            if (isTRUE(raster)) {
              p <- p + scattermore::geom_scattermore(
                data = cell_df, aes(x = .data[[f1]], y = .data[[f2]]), color = cols.highlight,
                pointsize = floor(sizes.highlight) + stroke.highlight, alpha = alpha.highlight, pixels = raster.dpi
              ) +
                scattermore::geom_scattermore(
                  data = cell_df, aes(x = .data[[f1]], y = .data[[f2]], color = .data[[group.by]]),
                  pointsize = floor(sizes.highlight), alpha = alpha.highlight, pixels = raster.dpi
                )
            } else {
              p <- p +
                suppressWarnings(geom_point(
                  data = cell_df, aes(x = .data[[f1]], y = .data[[f2]]), color = cols.highlight,
                  size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
                )) +
                suppressWarnings(geom_point(
                  data = cell_df, aes(x = .data[[f1]], y = .data[[f2]], color = .data[[group.by]]),
                  size = sizes.highlight, alpha = alpha.highlight
                ))
            }
          }
        }
      }
      if (f1_index > f2_index) {
        label <- paste0(f1, "\n", f2, "\nCor: ", round(pair_sim[f1, f2], 3)) # "\n","f1_index:",f1_index," ","f2_index:",f2_index
        label_pos <- (max(dat_exp[rownames(dat), ], na.rm = TRUE) - min(dat_exp[rownames(dat), ], na.rm = TRUE)) / 2
        fill <- rownames(df_bound)[df_bound[, 1] < pair_sim[f1, f2] & df_bound[, 2] >= pair_sim[f1, f2]]
        p <- p + annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = fill) +
          annotate(
            geom = "text_repel", x = label_pos, y = label_pos, label = label,
            fontface = "bold", color = "black", bg.color = "white", bg.r = 0.1, size = 3.5, point.size = NA
          )
      }

      if (f1_index == 1 & f2_index != 1) {
        p <- p + theme(
          axis.ticks.y = element_line(),
          axis.text.y = element_text(size = 10)
        )
      }
      if (f2_index == length(features) & f1_index != length(features)) {
        p <- p + theme(
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size = 10)
        )
      }
      if (f1_index == 1) {
        p <- p + labs(y = f2) + theme(axis.title.y = element_text(size = 12))
      }
      if (f2_index == length(features)) {
        p <- p + labs(x = f1) + theme(axis.title.x = element_text(size = 12))
      }
      p <- p + scale_color_manual(
        name = paste0(group.by, ":"),
        values = colors,
        labels = names(colors),
        na.value = bg_color
      ) + scale_fill_manual(
        name = paste0(group.by, ":"),
        values = colors,
        labels = names(colors),
        na.value = bg_color
      )
      return(p)
    }, x = order1, y = order2, SIMPLIFY = FALSE)
    if (nlevels(dat[[group.by]]) > 1) {
      legend <- suppressWarnings(get_legend(plotlist[[1]] +
        guides(fill = guide_legend(
          title.hjust = 0,
          keywidth = 0.05,
          keyheight = 0.05,
          default.unit = "inch",
          order = 1,
          override.aes = list(size = 4.5, color = "black", alpha = 1)
        )) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )))
    } else {
      legend <- NULL
    }
    grob_row <- list()
    plotlist <- suppressWarnings(lapply(plotlist, as_grob))
    for (i in seq(1, length(plotlist), length(features))) {
      grob_row[[paste0(i:(i + length(features) - 1), collapse = "-")]] <- do.call(cbind, plotlist[i:(i + length(features) - 1)])
    }
    gtable <- do.call(rbind, grob_row)
    if (!is.null(legend)) {
      gtable <- add_grob(gtable, legend, legend.position)
    }
    if (nlevels(dat_use[[split.by]]) > 1) {
      split_grob <- textGrob(s, just = "center", gp = gpar(fontface = "bold", fontsize = 13))
      gtable <- add_grob(gtable, split_grob, "top")
    }
    if (!is.null(subtitle)) {
      subtitle_grob <- textGrob(subtitle, x = 0, hjust = 0, gp = gpar(fontface = "italic", fontsize = 13))
      gtable <- add_grob(gtable, subtitle_grob, "top")
    }
    if (!is.null(title)) {
      title_grob <- textGrob(title, x = 0, hjust = 0, gp = gpar(fontsize = 14))
      gtable <- add_grob(gtable, title_grob, "top", 2 * grobHeight(title_grob))
    }
    p <- wrap_plots(gtable)
    plist[[paste0(s)]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' CellDensityPlot
#'
#' @param srt
#' @param features
#' @param group.by
#' @param split.by
#' @param flip
#' @param reverse
#' @param x_order
#' @param decreasing
#' @param palette
#' @param palcolor
#' @param cells
#' @param assay
#' @param slot
#' @param keep_empty
#' @param y.nbreaks
#' @param y.min
#' @param y.max
#' @param same.y.lims
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param force
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' CellDensityPlot(pancreas_sub, features = "Sox9", group.by = "SubCellType")
#'
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' CellDensityPlot(pancreas_sub, features = "Lineage1", group.by = "SubCellType", aspect.ratio = 1)
#' CellDensityPlot(pancreas_sub, features = "Lineage1", group.by = "SubCellType", flip = TRUE)
#'
#' @importFrom stats median
#' @importFrom dplyr %>% group_by_at summarise_at arrange_at pull desc
#' @importFrom ggplot2 ggplot scale_fill_manual labs scale_y_discrete scale_x_continuous facet_grid labs coord_flip element_text element_line
#' @importFrom patchwork wrap_plots
#' @importFrom methods slot
#' @export
CellDensityPlot <- function(srt, features, group.by, split.by = NULL, assay = NULL, slot = "data",
                            flip = FALSE, reverse = FALSE, x_order = c("value", "rank"),
                            decreasing = NULL, palette = "Paired", palcolor = NULL,
                            cells = NULL, keep_empty = FALSE,
                            y.nbreaks = 4, y.min = NULL, y.max = NULL, same.y.lims = FALSE,
                            aspect.ratio = NULL, title = NULL, subtitle = NULL,
                            legend.position = "right", legend.direction = "vertical",
                            theme_use = "theme_scp", theme_args = list(),
                            combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  check_R("ggridges")
  assay <- assay %||% DefaultAssay(srt)
  x_order <- match.arg(x_order)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in c(group.by, split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
  }

  features <- unique(features)
  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))]
  # print(colnames(srt@meta.data))
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (length(features_gene) > 0) {
    dat_gene <- t(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE])
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as.matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  dat_use <- cbind(dat_exp, srt@meta.data[row.names(dat_exp), c(group.by, split.by), drop = FALSE])
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  if (isTRUE(same.y.lims) && is.null(y.max)) {
    y.max <- max(as.matrix(dat_exp[, features])[is.finite(as.matrix(dat_exp[, features]))], na.rm = TRUE)
  }
  if (isTRUE(same.y.lims) && is.null(y.min)) {
    y.min <- min(as.matrix(dat_exp[, features])[is.finite(as.matrix(dat_exp[, features]))], na.rm = TRUE)
  }

  plist <- list()
  for (f in features) {
    for (g in group.by) {
      colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor)
      for (s in levels(dat_use[[split.by]])) {
        dat <- dat_use[dat_use[[split.by]] == s, , drop = FALSE]
        if (any(is.infinite(dat[, f]))) {
          dat[, f][dat[, f] == max(dat[, f], na.rm = TRUE)] <- max(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
          dat[, f][dat[, f] == min(dat[, f], na.rm = TRUE)] <- min(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
        }
        dat[, "cell"] <- rownames(dat)
        if (x_order == "value") {
          dat[, "value"] <- dat[, f]
        } else {
          dat[, "value"] <- rank(dat[, f])
        }
        dat[, "features"] <- f
        dat[, "split.by"] <- s
        dat <- dat[!is.na(dat[[f]]), , drop = FALSE]
        stat <- table(dat[[g]])
        stat_drop <- names(which(stat <= 2))
        if (length(stat_drop) > 0) {
          dat <- dat[!dat[[g]] %in% stat_drop, , drop = FALSE]
        }
        y_max_use <- y.max %||% suppressWarnings(max(dat[, "value"][is.finite(x = dat[, "value"])], na.rm = TRUE))
        y_min_use <- y.min %||% suppressWarnings(min(dat[, "value"][is.finite(x = dat[, "value"])], na.rm = TRUE))

        if (!is.null(decreasing)) {
          levels <- dat %>%
            group_by_at(g) %>%
            summarise_at(.funs = median, .vars = f, na.rm = TRUE) %>%
            arrange_at(.vars = f, .funs = if (decreasing) desc else list()) %>%
            pull(g) %>%
            as.character()
          dat[["order"]] <- factor(dat[[g]], levels = levels)
        } else {
          dat[["order"]] <- factor(dat[[g]], levels = rev(levels(dat[[g]])))
        }
        if (flip) {
          dat[["order"]] <- factor(dat[[g]], levels = levels(dat[[g]]))
          aspect.ratio <- 1 / aspect.ratio
          if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
            aspect.ratio <- NULL
          }
        }
        p <- ggplot(dat, aes(x = .data[["value"]], y = .data[["order"]], fill = .data[[g]])) +
          ggridges::geom_density_ridges()
        p <- p + scale_fill_manual(
          name = paste0(g, ":"),
          values = colors
        )
        y.trans <- ifelse(flip, "reverse", "identity")
        y.trans <- ifelse(reverse, setdiff(c("reverse", "identity"), y.trans), y.trans)

        limits <- if (y.trans == "reverse") c(y_max_use, y_min_use) else c(y_min_use, y_max_use)
        p <- p +
          scale_y_discrete(drop = !keep_empty, expand = c(0, 0)) +
          scale_x_continuous(
            limits = limits, trans = y.trans, n.breaks = y.nbreaks,
            expand = c(0, 0)
          )
        if (split.by != "All_cells") {
          p <- p + facet_grid(. ~ split.by)
        }
        p <- p + labs(title = title, subtitle = subtitle, x = f, y = g)
        if (isTRUE(flip)) {
          p <- p + do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              strip.text.x = element_text(angle = 0),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.ticks.x = element_line(),
              panel.grid.major.x = element_line(color = "grey", linetype = 2),
              legend.position = legend.position,
              legend.direction = legend.direction
            ) + coord_flip()
        } else {
          p <- p + do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              strip.text.y = element_text(angle = 0),
              axis.text.x = element_text(),
              axis.text.y = element_text(hjust = 1),
              axis.ticks.y = element_line(),
              panel.grid.major.y = element_line(color = "grey", linetype = 2),
              legend.position = legend.position,
              legend.direction = legend.direction
            )
        }
        plist[[paste0(f, ":", g, ":", paste0(s, collapse = ","))]] <- p
      }
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
  return(p)
}

#' LineagePlot
#'
#' @param srt
#' @param lineages
#' @param reduction
#' @param dims
#' @param trim
#' @param span
#' @param palette
#' @param palcolor
#' @param lineages_arrow
#' @param linewidth
#' @param line_bg
#' @param line_bg_stroke
#' @param whiskers
#' @param whiskers_linewidth
#' @param whiskers_alpha
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#' @param cells
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", show_plot = FALSE)
#' LineagePlot(pancreas_sub, lineages = paste0("Lineage", 1:3))
#' LineagePlot(pancreas_sub, lineages = paste0("Lineage", 1:3), whiskers = TRUE)
#' @importFrom Seurat Key Embeddings
#' @importFrom ggplot2 aes geom_path geom_segment labs
#' @importFrom stats loess quantile
#' @export
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
    fitted <- lapply(axes, function(x) {
      loess(formula(paste(x, l, sep = "~")), data = dat_sub, span = span, degree = 2)$fitted
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
#'
#' @param srt
#' @param paga
#' @param reduction
#' @param dims
#' @param show_transition
#' @param node_palette
#' @param node_size
#' @param node_alpha
#' @param node_highlight
#' @param node_highlight_color
#' @param label
#' @param label.size
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param label_insitu
#' @param label_repel
#' @param label_repulsion
#' @param label_point_size
#' @param label_point_color
#' @param label_segment_color
#' @param edge_threshold
#' @param edge_line
#' @param edge_line_curvature
#' @param edge_line_angle
#' @param edge_size
#' @param edge_color
#' @param edge_alpha
#' @param edge_shorten
#' @param edge_offset
#' @param edge_highlight
#' @param edge_highlight_color
#' @param transition_threshold
#' @param transition_line
#' @param transition_line_curvature
#' @param transition_line_angle
#' @param transition_size
#' @param transition_color
#' @param transition_alpha
#' @param transition_arrow_type
#' @param transition_arrow_angle
#' @param transition_arrow_length
#' @param transition_shorten
#' @param transition_offset
#' @param transition_highlight
#' @param transition_highlight_color
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#' @param type
#' @param cells
#' @param node_palcolor
#' @param theme_args
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP", return_seurat = TRUE)
#' PAGAPlot(pancreas_sub)
#' PAGAPlot(pancreas_sub, type = "connectivities_tree")
#' PAGAPlot(pancreas_sub, reduction = "PCA")
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
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
    node = dat, edge = as.matrix(connectivities), node_coord = paste0(reduction_key, dims),
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
#' @param node
#'
#' @param edge
#' @param transition
#' @param node_coord
#' @param node_group
#' @param node_palette
#' @param node_size
#' @param node_alpha
#' @param node_highlight
#' @param node_highlight_color
#' @param label
#' @param label.size
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param label_insitu
#' @param label_repel
#' @param label_repulsion
#' @param label_point_size
#' @param label_point_color
#' @param label_segment_color
#' @param edge_threshold
#' @param use_triangular
#' @param edge_line
#' @param edge_line_curvature
#' @param edge_line_angle
#' @param edge_color
#' @param edge_size
#' @param edge_alpha
#' @param edge_shorten
#' @param edge_offset
#' @param edge_highlight
#' @param edge_highlight_color
#' @param transition_threshold
#' @param transition_line
#' @param transition_line_curvature
#' @param transition_line_angle
#' @param transition_color
#' @param transition_size
#' @param transition_alpha
#' @param transition_arrow_type
#' @param transition_arrow_angle
#' @param transition_arrow_length
#' @param transition_shorten
#' @param transition_offset
#' @param transition_highlight
#' @param transition_highlight_color
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#' @param node_palcolor
#' @param theme_args
#'
#' @importFrom ggplot2 scale_size_identity scale_size_continuous scale_size_discrete scale_alpha_identity scale_alpha_continuous scale_alpha_discrete geom_curve geom_segment geom_point scale_color_manual guide_legend guides labs aes scale_linewidth_continuous
#' @importFrom ggnewscale new_scale
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow
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
  edge_df <- reshape2::melt(edge, na.rm = TRUE, stringsAsFactors = FALSE)
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
    trans2 <- trans1 <- as.matrix(transition)
    trans1[lower.tri(trans1)] <- 0
    trans2[upper.tri(trans2)] <- 0
    trans <- t(trans1) - trans2
    trans[abs(trans) <= transition_threshold] <- NA
    trans_df <- reshape2::melt(trans, na.rm = TRUE, stringsAsFactors = FALSE)
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
        keywidth = 0,
        keyheight = 0,
        default.unit = "inch",
        order = 1,
        override.aes = list(size = 4.5, alpha = 1)
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
        point.size = NA, max.overlaps = 100,
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
#' @param data
#' @param shorten_start
#' @param shorten_end
#' @param offset
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
#' @param srt
#' @param reduction
#' @param dims
#' @param velocity
#' @param plot_type
#' @param group_by
#' @param group_palette
#' @param n_neighbors
#' @param density
#' @param smooth
#' @param scale
#' @param min_mass
#' @param cutoff_perc
#' @param arrow_angle
#' @param arrow_flank
#' @param arrow_color
#' @param streamline_L
#' @param streamline_minL
#' @param streamline_res
#' @param streamline_n
#' @param streamlinewidth
#' @param streamline_alpha
#' @param streamline_color
#' @param streamline_palette
#' @param streamline_palcolor
#' @param streamline_bg_color
#' @param streamline_bg_stroke
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#' @param cells
#' @param group_palcolor
#' @param seed
#' @param theme_args
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
#' @export
VelocityPlot <- function(srt, reduction, dims = c(1, 2), cells = NULL, velocity = "stochastic", plot_type = c("raw", "grid", "stream"),
                         group_by = NULL, group_palette = "Paired", group_palcolor = NULL,
                         n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50), density = 1, smooth = 0.5, scale = 1, min_mass = 1, cutoff_perc = 5,
                         arrow_angle = 20, arrow_flank = 0.8, arrow_color = "black",
                         streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
                         streamlinewidth = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
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
          guide = guide_legend(
            title.hjust = 0,
            keywidth = 0,
            keyheight = 0,
            default.unit = "inch",
            order = 1,
            override.aes = list(size = 4.5, alpha = 1)
          )
        )
      )
    } else {
      velocity_layer <- list(
        geom_segment(
          data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v),
          color = arrow_color,
          arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle),
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        geom_segment(
          data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v),
          color = adjcolors(arrow_color, 0.2),
          arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle * (1 - (arrow_flank))),
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
      ),
      geom_segment(
        data = df_field, aes(x = x, y = y, xend = x + u, yend = y + v),
        color = adjcolors(arrow_color, 0.2),
        arrow = arrow(length = unit(df_field[["length_perc"]], "npc"), type = "closed", angle = arrow_angle * (1 - (arrow_flank))),
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
    u <- reshape2::melt(t(V_grid[1, , ]))
    v <- reshape2::melt(t(V_grid[2, , ]))
    df_field[, "u"] <- u$value
    df_field[, "v"] <- v$value
    df_field[is.na(df_field)] <- 0

    if (!is.null(streamline_color)) {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, size = max(streamlinewidth, na.rm = TRUE) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          n = streamline_n, size = max(streamlinewidth, na.rm = TRUE), color = streamline_color, alpha = streamline_alpha,
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
          n = streamline_n, size = max(streamlinewidth, na.rm = TRUE) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v, size = after_stat(step), color = sqrt(..dx..^2 + ..dy..^2)),
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
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)
        ),
        scale_size(range = range(streamlinewidth), guide = "none")
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
#' @param X_emb
#' @param V_emb
#' @param density
#' @param smooth
#' @param n_neighbors
#' @param min_mass
#' @param scale
#' @param adjust_for_stream
#' @param cutoff_perc
#'
#' @importFrom SeuratObject as.sparse
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
  X_grid <- as.matrix(expand.grid(grs))

  d <- proxyC::dist(
    x = as.sparse(X_emb),
    y = as.sparse(X_grid),
    method = "euclidean",
    use_nan = TRUE
  )
  neighbors <- t(as.matrix(apply(d, 2, function(x) order(x, decreasing = FALSE)[1:n_neighbors])))
  dists <- t(as.matrix(apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1:n_neighbors]])))

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
#' @param srt
#' @param group_by
#' @param test.use
#' @param DE_threshold
#' @param x_metric
#' @param palette
#' @param palcolor
#' @param pt.size
#' @param pt.alpha
#' @param cols.highlight
#' @param sizes.highlight
#' @param alpha.highlight
#' @param stroke.highlight
#' @param nlabel
#' @param features_label
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param label.size
#' @param aspect.ratio
#' @param xlab
#' @param ylab
#' @param combine
#' @param nrow
#' @param ncol
#' @param theme_use
#' @param theme_args
#' @param byrow
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
VolcanoPlot <- function(srt, group_by = NULL, test.use = "wilcox", DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                        x_metric = "diff_pct", palette = "RdBu", palcolor = NULL, pt.size = 1, pt.alpha = 1,
                        cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                        nlabel = 5, features_label = NULL, label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                        aspect.ratio = NULL, xlab = x_metric, ylab = "-log10(p-adjust)",
                        theme_use = "theme_scp", theme_args = list(),
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  if (is.null(group_by)) {
    group_by <- "custom"
  }
  slot <- paste0("DEtest_", group_by)
  if (!slot %in% names(srt@tools) || length(grep(pattern = "AllMarkers", names(srt@tools[[slot]]))) == 0) {
    stop("Cannot find the DEtest result for the group '", group_by, "'. You may perform RunDEtest first.")
  }
  index <- grep(pattern = paste0("AllMarkers_", test.use), names(srt@tools[[slot]]))[1]
  if (is.na(index)) {
    stop("Cannot find the 'AllMarkers_", test.use, "' in the DEtest result.")
  }
  de <- names(srt@tools[[slot]])[index]
  de_df <- srt@tools[[slot]][[de]]
  de_df[, "diff_pct"] <- de_df[, "pct.1"] - de_df[, "pct.2"]
  de_df[, "-log10padj"] <- -log10(de_df[, "p_val_adj"])
  de_df[, "DE"] <- FALSE
  de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), "DE"] <- TRUE

  x_upper <- quantile(de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])], c(0.99, 1))
  x_lower <- quantile(de_df[["avg_log2FC"]][is.finite(de_df[["avg_log2FC"]])], c(0.01, 0))
  x_upper <- ifelse(x_upper[1] > 0, x_upper[1], x_upper[2])
  x_lower <- ifelse(x_lower[1] < 0, x_lower[1], x_lower[2])
  if (x_upper > 0 & x_lower < 0) {
    value_range <- min(abs(c(x_upper, x_lower)), na.rm = TRUE)
    x_upper <- value_range
    x_lower <- -value_range
  }

  de_df[, "border"] <- FALSE
  de_df[de_df[["avg_log2FC"]] > x_upper, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] > x_upper, "avg_log2FC"] <- x_upper
  de_df[de_df[["avg_log2FC"]] < x_lower, "border"] <- TRUE
  de_df[de_df[["avg_log2FC"]] < x_lower, "avg_log2FC"] <- x_lower

  de_df[, "y"] <- -log10(de_df[, "p_val_adj"])
  if (x_metric == "diff_pct") {
    de_df[, "x"] <- de_df[, "diff_pct"]
    de_df[de_df[, "avg_log2FC"] < 0, "y"] <- -de_df[de_df[, "avg_log2FC"] < 0, "y"]
    de_df <- de_df[order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), , drop = FALSE]
  } else if (x_metric == "avg_log2FC") {
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df[de_df[, "diff_pct"] < 0, "y"] <- -de_df[de_df[, "diff_pct"] < 0, "y"]
    de_df <- de_df[order(abs(de_df[, "diff_pct"]), decreasing = FALSE, na.last = FALSE), , drop = FALSE]
  }
  de_df[, "distance"] <- de_df[, "x"]^2 + de_df[, "y"]^2

  plist <- list()
  for (group in levels(de_df[["group1"]])) {
    df <- de_df[de_df[["group1"]] == group, , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }
    x_nudge <- diff(range(df$x)) * 0.05
    df[, "label"] <- FALSE
    if (is.null(features_label)) {
      df[df[["y"]] >= 0, ][head(order(df[df[["y"]] >= 0, "distance"], decreasing = TRUE), nlabel), "label"] <- TRUE
      df[df[["y"]] < 0, ][head(order(df[df[["y"]] < 0, "distance"], decreasing = TRUE), nlabel), "label"] <- TRUE
    } else {
      df[df[["gene"]] %in% features_label, "label"] <- TRUE
    }
    jitter <- position_jitter(width = 0.2, height = 0.2, seed = 11)
    color_by <- ifelse(x_metric == "diff_pct", "avg_log2FC", "diff_pct")
    p <- ggplot() +
      geom_point(data = df[!df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[!df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight) +
      geom_point(data = df[df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[df[["DE"]] & df[["border"]], , drop = FALSE], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      geom_vline(xintercept = 0, color = "grey", linetype = 2) +
      geom_text_repel(
        data = df[df[["label"]], , drop = FALSE], aes(x = x, y = y, label = gene),
        min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40",
        color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, force = 20,
        nudge_x = ifelse(df[df[["label"]], "y"] >= 0, -x_nudge, x_nudge)
      ) +
      labs(x = xlab, y = ylab) +
      scale_color_gradientn(
        name = ifelse(x_metric == "diff_pct", "log2FC", "diff_pct"), colors = palette_scp(palette = palette, palcolor = palcolor),
        values = rescale(unique(c(min(c(df[, color_by], 0), na.rm = TRUE), 0, max(df[, color_by], na.rm = TRUE)))),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0, order = 1)
      ) +
      scale_y_continuous(labels = abs) +
      facet_wrap(~group1) +
      do.call(theme_use, theme_args) +
      theme(aspect.ratio = aspect.ratio)
    plist[[group]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' SankeyPlot
#'
#' @param node
#' @param edge
#' @param node_group
#' @param group_list
#'
#' @export
SankeyPlot <- function(node, edge, node_group = NULL, group_list = NULL) {
  check_R("networkD3")

  if (is.null(group_list)) {
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
    node[["name"]] <- node[[node_group]]
    edge_df <- reshape2::melt(edge, na.rm = TRUE, stringsAsFactors = FALSE)
    colnames(edge_df) <- c("from", "to", "size")
    edge_df[["from"]] <- match(edge_df[["from"]], node[["name"]]) - 1
    edge_df[["to"]] <- match(edge_df[["to"]], node[["name"]]) - 1
    edge_df <- edge_df[edge_df[["size"]] > 0, , drop = FALSE]
    edge_df <- aggregate(size ~ from + to, data = edge_df, sum)
  } else {
    if (length(unique(unlist(lapply(group_list, length)))) != 1) {
      stop("Elements in the group_list must be the same length.")
    }
    group_list <- lapply(group_list, as.factor)
    all_levels <- unique(unlist(lapply(group_list, levels)))
    group_list <- lapply(group_list, as.character)
    edge_df <- do.call(cbind.data.frame, group_list)
    if (ncol(edge_df) >= 3) {
      edge_df_raw <- edge_df
      edge_df <- edge_df_raw[, 1:2]
      for (i in 3:ncol(edge_df_raw)) {
        edge_df <- rbind(edge_df, edge_df_raw[, c(i - 1, i)])
      }
    }
    colnames(edge_df) <- c("from", "to")
    node <- data.frame(name = unique(unlist(group_list)))
    node[["name"]] <- factor(node[["name"]], all_levels)
    edge_df[["from"]] <- match(edge_df[["from"]], node[["name"]]) - 1
    edge_df[["to"]] <- match(edge_df[["to"]], node[["name"]]) - 1
    edge_df[["size"]] <- 1
    edge_df <- aggregate(size ~ from + to, data = edge_df, sum)
  }

  p <- networkD3::sankeyNetwork(
    Links = edge_df, Nodes = node,
    Source = "from", Target = "to",
    Value = "size", NodeID = "name",
    sinksRight = FALSE, nodeWidth = 40, fontSize = 15, nodePadding = 20
  )

  return(p)
}

fc_matrix <- function(matrix) {
  matrix / rowMeans(matrix)
}
zscore_matrix <- function(matrix, ...) {
  t(scale(t(matrix), ...))
}
log2fc_matrix <- function(matrix) {
  log2(matrix / rowMeans(matrix))
}
log1p_matrix <- function(matrix) {
  log1p(matrix)
}
matrix_process <- function(matrix, method = c("raw", "zscore", "fc", "log2fc", "log1p"), ...) {
  if (is.function(method)) {
    matrix_processed <- method(matrix, ...)
  } else if (method == "raw") {
    matrix_processed <- matrix
  } else if (method == "fc") {
    matrix_processed <- fc_matrix(matrix)
  } else if (method == "zscore") {
    matrix_processed <- zscore_matrix(matrix, ...)
  } else if (method == "log2fc") {
    matrix_processed <- log2fc_matrix(matrix)
  } else if (method == "log1p") {
    matrix_processed <- log1p_matrix(matrix)
  }
  if (!identical(dim(matrix_processed), dim(matrix))) {
    stop("The dimensions of the matrix are changed after processing")
  }
  return(matrix_processed)
}

extractgrobs <- function(vlnplots, x_nm, y_nm, x, y) {
  grobs <- vlnplots[paste0(x_nm[x], ":", y_nm[y])]
  if (length(grobs) == 1) {
    grobs <- grobs[[1]]
  }
  return(grobs)
}

#' @importFrom grid viewport grid.draw is.grob
grid_draw <- function(groblist, x, y, width, height) {
  if (is.grob(groblist)) {
    groblist <- list(groblist)
  }
  for (i in seq_along(groblist)) {
    groblist[[i]]$vp <- viewport(x = x[i], y = y[i], width = width[i], height = height[i])
    grid.draw(groblist[[i]])
  }
}

#' @importFrom stats hclust dist as.dendrogram order.dendrogram
#' @importFrom ComplexHeatmap merge_dendrogram
cluster_within_group2 <- function(mat, factor) {
  require("dendextend", quietly = TRUE)
  if (!is.factor(factor)) {
    factor <- factor(factor, levels = unique(factor))
  }
  dend_list <- list()
  order_list <- list()
  for (le in unique(levels(factor))) {
    m <- mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] <- which(factor == le)
      dend_list[[le]] <- structure(which(factor == le),
        class = "dendrogram", leaf = TRUE, # height = 0,
        label = 1, members = 1
      )
    } else if (ncol(m) > 1) {
      hc1 <- hclust(dist(t(m)))
      dend_list[[le]] <- as.dendrogram(hc1)
      order_list[[le]] <- which(factor == le)[order.dendrogram(dend_list[[le]])]
      order.dendrogram(dend_list[[le]]) <- order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") <- le
  }
  parent <- as.dendrogram(hclust(dist(t(sapply(
    order_list,
    function(x) rowMeans(mat[, x, drop = FALSE])
  )))))
  dend_list <- lapply(dend_list, function(dend) {
    dendrapply(
      dend,
      function(node) {
        if (is.null(attr(node, "height"))) {
          attr(node, "height") <- 0
        }
        node
      }
    )
  })
  # print(sapply(dend_list, function(x) attr(x, "height")))
  dend <- merge_dendrogram(parent, dend_list)
  order.dendrogram(dend) <- unlist(order_list[order.dendrogram(parent)])
  return(dend)
}

#' @importFrom ComplexHeatmap HeatmapAnnotation anno_empty anno_block anno_textbox
#' @importFrom grid gpar
#' @importFrom dplyr %>% filter group_by arrange desc across reframe mutate distinct n .data "%>%"
heatmap_enrichment <- function(geneID, geneID_groups, feature_split_palette = "jama", feature_split_palcolor = NULL, ha_right = NULL, flip = FALSE,
                               anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                               terms_width = unit(4, "in"), terms_fontsize = 8,
                               keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                               features_width = unit(2, "in"), features_fontsize = c(6, 10),
                               IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_combine = FALSE, db_version = "latest", convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                               db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                               GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                               pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, min_word_length = 3,
                               exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "organization", "system", "regulation", "positive", "negative", "response", "process")) {
  res <- NULL
  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    if (isTRUE(flip)) {
      stop("anno_keys, anno_features and anno_terms can only be used when flip is FALSE.")
    }
    if (all(is.na(geneID_groups))) {
      geneID_groups <- rep(1, length(geneID))
    }
    if (!is.factor(geneID_groups)) {
      geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
    }
    fill_split <- palette_scp(levels(geneID_groups), type = "discrete", palette = feature_split_palette, palcolor = feature_split_palcolor)[levels(geneID_groups) %in% geneID_groups]
    res <- RunEnrichment(
      geneID = geneID, geneID_groups = geneID_groups, IDtype = IDtype, species = species,
      db_update = db_update, db_version = db_version, db_combine = db_combine, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror,
      db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize,
      GO_simplify = GO_simplify, GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff
    )
    if (isTRUE(db_combine)) {
      db <- "Combined"
    }
    if (isTRUE(GO_simplify) && any(db %in% c("GO_BP", "GO_CC", "GO_MF"))) {
      db[db %in% c("GO_BP", "GO_CC", "GO_MF")] <- paste0(db[db %in% c("GO_BP", "GO_CC", "GO_MF")], "_sim")
    }
    if (nrow(res$enrichment) == 0) {
      warning("No enrichment result found.", immediate. = TRUE)
    } else {
      metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
      metric_value <- ifelse(is.null(padjustCutoff), pvalueCutoff, padjustCutoff)
      pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
      padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

      df <- res$enrichment
      df <- df[df[["Database"]] %in% db, , drop = FALSE]
      df <- df[df[[metric]] < metric_value, , drop = FALSE]
      df <- df[order(df[[metric]]), , drop = FALSE]
      if (nrow(df) == 0) {
        warning(
          "No term enriched using the threshold: ",
          paste0("pvalueCutoff = ", pvalueCutoff), "; ",
          paste0("padjustCutoff = ", padjustCutoff),
          immediate. = TRUE
        )
      } else {
        df_list <- split.data.frame(df, ~ Database + Groups)
        df_list <- df_list[lapply(df_list, nrow) > 0]

        for (enrich in db) {
          nm <- strsplit(names(df_list), "\\.")
          subdf_list <- df_list[unlist(lapply(nm, function(x) x[[1]])) %in% enrich]
          if (length(subdf_list) == 0) {
            warning(
              "No ", enrich, " term enriched using the threshold: ",
              paste0("pvalueCutoff = ", pvalueCutoff), "; ",
              paste0("padjustCutoff = ", padjustCutoff),
              immediate. = TRUE
            )
            next
          }
          nm <- strsplit(names(subdf_list), "\\.")

          ha_terms <- NULL
          if (isTRUE(anno_terms)) {
            terms_list <- lapply(subdf_list, function(df) {
              if (isTRUE(show_termid)) {
                terms <- paste(head(df$ID, topTerm), head(df$Description, topTerm))
              } else {
                terms <- head(df$Description, topTerm)
                terms <- capitalize(terms)
              }
              df_out <- data.frame(keyword = terms)
              df_out[["col"]] <- palette_scp(-log10(head(df[, metric], topTerm)), type = "continuous", palette = "Spectral", matched = TRUE)
              df_out[["col"]] <- sapply(df_out[["col"]], function(x) blendcolors(c(x, "black")))
              df_out[["fontsize"]] <- rep(terms_fontsize, nrow(df_out))
              return(df_out)
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(terms_list))) > 0) {
              ha_terms <- HeatmapAnnotation(
                "terms_empty" = anno_empty(width = unit(0.05, "in"), border = FALSE, which = "row"),
                "terms_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "terms" = anno_textbox(
                  align_to = geneID_groups, text = terms_list, max_width = terms_width,
                  word_wrap = TRUE, add_new_line = TRUE,
                  background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE,
                  which = "row"
                ),
                which = "row", gap = unit(0, "points")
              )
              names(ha_terms) <- paste0(names(ha_terms), "_", enrich)
            }
          }

          ha_keys <- NULL
          if (isTRUE(anno_keys)) {
            check_R("jokergoo/simplifyEnrichment")
            keys_list <- lapply(subdf_list, function(df) {
              if (df$Database[1] %in% c("GO_BP", "GO_CC", "GO_MF")) {
                df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df[["ID"]])
                if (nrow(df0) > 0) {
                  df <- df0 %>%
                    reframe(
                      keyword = .data[["keyword"]],
                      score = -(log10(.data[["padj"]])),
                      count = .data[["n_term"]],
                      Database = df[["Database"]][1],
                      Groups = df[["Groups"]][1]
                    ) %>%
                    filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
                    filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
                    filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
                    distinct() %>%
                    mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                    as.data.frame()
                  df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
                } else {
                  df <- NULL
                }
              } else {
                df <- df %>%
                  mutate(keyword = strsplit(tolower(as.character(.data[["Description"]])), " ")) %>%
                  unnest(cols = "keyword") %>%
                  group_by(.data[["keyword"]], Database, Groups) %>%
                  reframe(
                    keyword = .data[["keyword"]],
                    score = sum(-(log10(.data[[metric]]))),
                    count = n(),
                    Database = .data[["Database"]],
                    Groups = .data[["Groups"]],
                    .groups = "keep"
                  ) %>%
                  filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
                  filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
                  filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
                  distinct() %>%
                  mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                  as.data.frame()
                df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
              }
              if (isTRUE(nrow(df) > 0)) {
                df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
                df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
                df[["fontsize"]] <- rescale(df[, "count"], to = keys_fontsize)
                return(df)
              } else {
                return(NULL)
              }
            })
            names(keys_list) <- unlist(lapply(nm, function(x) x[[2]]))
            keys_list <- keys_list[lapply(keys_list, length) > 0]
            if (length(intersect(geneID_groups, names(keys_list))) > 0) {
              ha_keys <- HeatmapAnnotation(
                "keys_empty" = anno_empty(width = unit(0.05, "in"), border = FALSE, which = "row"),
                "keys_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "keys" = anno_textbox(
                  align_to = geneID_groups, text = keys_list, max_width = keys_width,
                  background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE,
                  which = "row"
                ),
                which = "row", gap = unit(0, "points")
              )
              names(ha_keys) <- paste0(names(ha_keys), "_", enrich)
            }
          }

          ha_features <- NULL
          if (isTRUE(anno_features)) {
            features_list <- lapply(subdf_list, function(df) {
              df <- df %>%
                mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
                unnest(cols = "keyword") %>%
                group_by(.data[["keyword"]], Database, Groups) %>%
                reframe(
                  keyword = .data[["keyword"]],
                  score = sum(-(log10(.data[[metric]]))),
                  count = n(),
                  Database = .data[["Database"]],
                  Groups = .data[["Groups"]],
                  .groups = "keep"
                ) %>%
                distinct() %>%
                mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                as.data.frame()
              df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
              df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
              df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
              df[["fontsize"]] <- rescale(df[, "count"], to = features_fontsize)
              return(df)
            })
            names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(features_list))) > 0) {
              ha_features <- HeatmapAnnotation(
                "features_empty" = anno_empty(width = unit(0.05, "in"), border = FALSE, which = "row"),
                "features_split" = anno_block(
                  gp = gpar(fill = fill_split),
                  width = unit(0.1, "in"),
                  which = "row"
                ),
                "features" = anno_textbox(
                  align_to = geneID_groups, text = features_list, max_width = features_width,
                  background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE,
                  which = "row"
                ),
                which = "row", gap = unit(0, "points")
              )
              names(ha_features) <- paste0(names(ha_features), "_", enrich)
            }
          }

          ha_enrichment <- list(ha_terms, ha_keys, ha_features)
          ha_enrichment <- ha_enrichment[sapply(ha_enrichment, length) > 0]
          ha_enrichment <- do.call(c, ha_enrichment)

          if (is.null(ha_right)) {
            ha_right <- ha_enrichment
          } else {
            ha_right <- c(ha_right, ha_enrichment)
          }
        }
      }
    }
  }
  return(list(ha_right = ha_right, res = res))
}

#' @importFrom grid convertWidth convertHeight unit
#' @importFrom ComplexHeatmap width.HeatmapAnnotation height.HeatmapAnnotation width.Legends
heatmap_rendersize <- function(width, height, units, ha_top_list, ha_left, ha_right, ht_list, legend_list, flip) {
  width_annotation <- height_annotation <- 0
  if (isTRUE(flip)) {
    width_sum <- width[1] %||% convertWidth(unit(1, "in"), units, valueOnly = TRUE)
    height_sum <- sum(height %||% convertHeight(unit(1, "in"), units, valueOnly = TRUE))
    if (length(ha_top_list) > 0) {
      width_annotation <- convertWidth(unit(width_annotation, units) + width.HeatmapAnnotation(ha_top_list[[1]]), units, valueOnly = TRUE)
    }
    if (!is.null(ha_left)) {
      height_annotation <- convertHeight(unit(height_annotation, units) + height.HeatmapAnnotation(ha_left), units, valueOnly = TRUE)
    }
    if (!is.null(ha_right)) {
      height_annotation <- convertHeight(unit(height_annotation, units) + height.HeatmapAnnotation(ha_right), units, valueOnly = TRUE)
    }
  } else {
    width_sum <- sum(width %||% convertWidth(unit(1, "in"), units, valueOnly = TRUE))
    height_sum <- height[1] %||% convertHeight(unit(1, "in"), units, valueOnly = TRUE)
    if (length(ha_top_list) > 0) {
      height_annotation <- convertHeight(unit(height_annotation, units) + height.HeatmapAnnotation(ha_top_list[[1]]), units, valueOnly = TRUE)
    }
    if (!is.null(ha_left)) {
      width_annotation <- convertWidth(unit(width_annotation, units) + width.HeatmapAnnotation(ha_left), units, valueOnly = TRUE)
    }
    if (!is.null(ha_right)) {
      width_annotation <- convertWidth(unit(width_annotation, units) + width.HeatmapAnnotation(ha_right), units, valueOnly = TRUE)
    }
  }
  dend_width <- name_width <- NULL
  dend_height <- name_height <- NULL
  if (inherits(ht_list, "HeatmapList")) {
    for (nm in names(ht_list@ht_list)) {
      ht <- ht_list@ht_list[[nm]]
      dend_width <- max(ht@row_dend_param$width, dend_width)
      dend_height <- max(ht@column_dend_param$height, dend_height)
      name_width <- max(ht@row_names_param$max_width, name_width)
      name_height <- max(ht@column_names_param$max_height, name_height)
    }
  } else if (inherits(ht_list, "Heatmap")) {
    ht <- ht_list
    dend_width <- max(ht@row_dend_param$width, dend_width)
    dend_height <- max(ht@column_dend_param$height, dend_height)
    name_width <- max(ht@row_names_param$max_width, name_width)
    name_height <- max(ht@column_names_param$max_height, name_height)
  } else {
    stop("ht_list is not a class of HeatmapList or Heatmap.")
  }

  lgd_width <- convertWidth(unit(unlist(lapply(legend_list, width.Legends)), unitType(width.Legends(legend_list[[1]]))), unitTo = units, valueOnly = TRUE)
  width_sum <- convertWidth(unit(width_sum, units) +
    unit(width_annotation, units) +
    dend_width +
    name_width, units, valueOnly = TRUE) + sum(lgd_width)
  height_sum <- max(
    convertHeight(unit(height_sum, units) +
      unit(height_annotation, units) +
      dend_height +
      name_height, units, valueOnly = TRUE),
    convertHeight(unit(0.95, "npc"), units, valueOnly = TRUE)
  )
  return(list(width_sum = width_sum, height_sum = height_sum))
}

#' @importFrom grid convertWidth convertHeight convertUnit unit grid.grabExpr
#' @importFrom ComplexHeatmap draw
#' @importFrom methods slotNames
heatmap_fixsize <- function(width, width_sum, height, height_sum, units, ht_list, legend_list, verbose = TRUE) {
  if (verbose) {
    message("The size of the heatmap will be fixed as some elements are not scalable.")
  }
  gTree <- grid.grabExpr(
    {
      ht <- draw(ht_list, annotation_legend_list = legend_list)
      ht_width <- ComplexHeatmap:::width(ht)
      ht_height <- ComplexHeatmap:::height(ht)
      if (inherits(ht_list, "HeatmapList")) {
        for (nm in names(ht_list@ht_list)) {
          if (is.null(names(width))) {
            width_fix <- width[1]
          } else {
            width_fix <- width[nm]
          }
          if (is.null(names(height))) {
            height_fix <- height[1]
          } else {
            height_fix <- height[nm]
          }
          ht_list@ht_list[[nm]]@matrix_param$width <- unit(width_fix %||% dim(ht_list@ht_list[[nm]]@matrix)[1], units = "null")
          ht_list@ht_list[[nm]]@matrix_param$height <- unit(height_fix %||% dim(ht_list@ht_list[[nm]]@matrix)[2], units = "null")
        }
      } else if (inherits(ht_list, "Heatmap")) {
        ht_list@matrix_param$width <- unit(width[1] %||% dim(ht_list@matrix)[1], units = "null")
        ht_list@matrix_param$height <- unit(height[1] %||% dim(ht_list@matrix)[2], units = "null")
      } else {
        stop("ht_list is not a class of HeatmapList or Heatmap.")
      }
    },
    width = unit(width_sum, units = units),
    height = unit(height_sum, units = units),
    wrap = TRUE,
    wrap.grobs = TRUE
  )
  if (unitType(ht_width) == "npc") {
    ht_width <- unit(width_sum, units = units)
  }
  if (unitType(ht_height) == "npc") {
    ht_height <- unit(height_sum, units = units)
  }
  if (is.null(width)) {
    ht_width <- max(
      convertWidth(ht@layout$max_left_component_width, units, valueOnly = TRUE) +
        convertWidth(ht@layout$max_right_component_width, units, valueOnly = TRUE) +
        convertWidth(sum(ht@layout$max_title_component_width), units, valueOnly = TRUE) +
        convertWidth(ht@annotation_legend_param$size[1], units, valueOnly = TRUE) +
        convertWidth(unit(1, "in"), units, valueOnly = TRUE),
      convertWidth(unit(0.95, "npc"), units, valueOnly = TRUE)
    )
    ht_width <- unit(ht_width, units)
  }
  if (is.null(height)) {
    ht_height <- max(
      convertHeight(ht@layout$max_top_component_height, units, valueOnly = TRUE) +
        convertHeight(ht@layout$max_bottom_component_height, units, valueOnly = TRUE) +
        convertHeight(sum(ht@layout$max_title_component_height), units, valueOnly = TRUE) +
        convertHeight(unit(1, "in"), units, valueOnly = TRUE),
      convertHeight(ht@annotation_legend_param$size[2], units, valueOnly = TRUE),
      convertHeight(unit(0.95, "npc"), units, valueOnly = TRUE)
    )
    ht_height <- unit(ht_height, units)
  }
  ht_width <- convertUnit(ht_width, unitTo = units)
  ht_height <- convertUnit(ht_height, unitTo = units)
  return(list(ht_width = ht_width, ht_height = ht_height))
}

#' GroupHeatmap
#'
#' @param srt A \code{Seurat} object.
#' @param features A vector of gene names to plot.
#' @param feature_split A vector of group names for features.
#' @param group.by Columns used to calculate cell expression. One heatmap per column name.
#' @param exp_method Method used to calculate cell expression.
#' @param assay Assay used to calculate the expression.
#' @param heatmap_palette Heatmap expression palette.
#' @param feature_annotation_palette Feature groups palette.
#' @param cell_annotation_palette Column palette.
#' @param aggregate_fun
#' @param slot
#' @param lib_normalize
#' @param libsize
#' @param add_reticle
#' @param cluster_rows
#' @param cluster_columns
#' @param split.by
#' @param cells
#' @param exp_cutoff
#' @param border
#' @param flip
#' @param feature_split_by
#' @param n_split
#' @param split_method
#' @param decreasing
#' @param cluster_features_by
#' @param cluster_row_slices
#' @param cluster_column_slices
#' @param show_row_names
#' @param show_column_names
#' @param row_names_side
#' @param column_names_side
#' @param row_names_rot
#' @param column_names_rot
#' @param row_title_side
#' @param column_title_side
#' @param row_title_rot
#' @param column_title_rot
#' @param anno_terms
#' @param anno_keys
#' @param anno_features
#' @param terms_width
#' @param terms_fontsize
#' @param keys_width
#' @param keys_fontsize
#' @param features_width
#' @param features_fontsize
#' @param IDtype
#' @param species
#' @param db_update
#' @param db_version
#' @param convert_species
#' @param Ensembl_version
#' @param mirror
#' @param db
#' @param TERM2GENE
#' @param TERM2NAME
#' @param minGSSize
#' @param maxGSSize
#' @param GO_simplify
#' @param GO_simplify_cutoff
#' @param simplify_method
#' @param simplify_similarityCutoff
#' @param pvalueCutoff
#' @param padjustCutoff
#' @param topTerm
#' @param show_termid
#' @param topWord
#' @param min_word_length
#' @param exclude_words
#' @param nlabel
#' @param features_label
#' @param label_size
#' @param label_color
#' @param add_bg
#' @param bg_alpha
#' @param add_dot
#' @param dot_size
#' @param reticle_color
#' @param add_violin
#' @param fill.by
#' @param fill_palette
#' @param fill_palcolor
#' @param heatmap_palcolor
#' @param group_palette
#' @param group_palcolor
#' @param cell_split_palette
#' @param cell_split_palcolor
#' @param feature_split_palette
#' @param feature_split_palcolor
#' @param cell_annotation
#' @param cell_annotation_palcolor
#' @param cell_annotation_params
#' @param feature_annotation
#' @param feature_annotation_palcolor
#' @param feature_annotation_params
#' @param use_raster
#' @param raster_device
#' @param height
#' @param width
#' @param units
#' @param seed
#' @param ht_params
#' @param grouping.var
#' @param numerator
#' @param limits
#' @param raster_by_magick
#' @param fuzzification
#' @param show_fuzzification
#' @param db_combine
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
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
#' ht1$plot
#' panel_fix(ht1$plot, height = 4, width = 6, raster = TRUE, dpi = 50)
#'
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' de_filter <- filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#'
#' ht2 <- GroupHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, group.by = "CellType",
#'   split.by = "Phase", cell_split_palette = "Dark2",
#'   nlabel = 10, show_row_names = FALSE
#' )
#' ht2$plot
#'
#' ht3 <- GroupHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, feature_split = de_filter$group1, group.by = "CellType",
#'   nlabel = 20, show_row_names = FALSE,
#'   species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "SP"))
#' de_top <- de_filter %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ht4 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
#'   cell_annotation_params = list(height = grid::unit(0.5, "in")),
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   add_dot = TRUE, add_bg = TRUE, show_row_names = TRUE
#' )
#' ht4$plot
#'
#' ht5 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
#'   cell_annotation_params = list(width = grid::unit(0.5, "in")),
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   add_dot = TRUE, add_bg = TRUE,
#'   flip = TRUE, column_title_rot = 45, show_column_names = TRUE
#' )
#' ht5$plot
#'
#' ht6 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   add_violin = TRUE, cluster_rows = TRUE,
#'   show_row_names = TRUE
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   add_violin = TRUE, fill.by = "expression", fill_palette = "Blues", cluster_rows = TRUE,
#'   show_row_names = TRUE
#' )
#' ht7$plot
#'
#' ht8 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, group.by = "CellType", split.by = "Phase", n_split = 4,
#'   cluster_rows = TRUE, cluster_columns = TRUE, cluster_row_slices = TRUE, cluster_column_slices = TRUE,
#'   add_dot = TRUE, add_reticle = TRUE, heatmap_palette = "viridis",
#'   show_row_names = TRUE, ht_params = list(row_gap = grid::unit(0, "mm"), row_names_gp = grid::gpar(fontsize = 10))
#' )
#' ht8$plot
#'
#' @importFrom circlize colorRamp2
#' @importFrom stats aggregate formula quantile
#' @importFrom ComplexHeatmap Legend HeatmapAnnotation anno_block anno_simple anno_customize Heatmap draw pindex restore_matrix %v%
#' @importFrom grid gpar grid.grabExpr grid.lines grid.rect grid.points grid.draw
#' @importFrom ggplot2 theme_void theme facet_null
#' @importFrom patchwork wrap_plots
#' @importFrom methods getFunction
#' @importFrom dplyr %>% filter group_by arrange desc across mutate distinct n .data "%>%"
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @export
GroupHeatmap <- function(srt, features = NULL, group.by = NULL, split.by = NULL, grouping.var = NULL, numerator = NULL, cells = NULL,
                         aggregate_fun = base::mean, exp_cutoff = 0, border = TRUE, flip = FALSE,
                         slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), limits = NULL, lib_normalize = identical(slot, "counts"), libsize = NULL,
                         feature_split = NULL, feature_split_by = NULL, n_split = NULL, split_method = c("kmeans", "hclust", "mfuzz"), decreasing = FALSE, fuzzification = NULL, show_fuzzification = FALSE,
                         cluster_features_by = NULL, cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                         show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                         row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                         anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                         terms_width = unit(4, "in"), terms_fontsize = 8,
                         keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                         features_width = unit(2, "in"), features_fontsize = c(6, 10),
                         IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                         db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                         GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                         pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, min_word_length = 3,
                         exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "organization", "system", "regulation", "positive", "negative", "response", "process"),
                         nlabel = 0, features_label = NULL, label_size = 10, label_color = "black",
                         add_bg = FALSE, bg_alpha = 0.5,
                         add_dot = FALSE, dot_size = unit(8, "mm"),
                         add_reticle = FALSE, reticle_color = "grey",
                         add_violin = FALSE, fill.by = "feature", fill_palette = "Dark2", fill_palcolor = NULL,
                         heatmap_palette = "RdBu", heatmap_palcolor = NULL, group_palette = "Paired", group_palcolor = NULL,
                         cell_split_palette = "jama", cell_split_palcolor = NULL, feature_split_palette = "jama", feature_split_palcolor = NULL,
                         cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = grid::unit(1, "cm")) else list(height = grid::unit(1, "cm")),
                         feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = list(),
                         use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, height = NULL, width = NULL, units = "inch",
                         seed = 11, ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }
  if (is.null(features)) {
    stop("No feature provided.")
  }

  split_method <- match.arg(split_method)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), "(", data_nm, ")")
  } else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }
  if (!is.null(grouping.var) && exp_method != "log2fc") {
    warning("When 'grouping.var' is specified, 'exp_method' can only be 'log2fc'", immediate. = TRUE)
    exp_method <- "log2fc"
  }
  if (!is.null(grouping.var)) {
    if (identical(split.by, grouping.var)) {
      stop("'grouping.var' must be different from 'split.by'")
    }
    if (!is.factor(srt@meta.data[[grouping.var]])) {
      srt@meta.data[[grouping.var]] <- factor(srt@meta.data[[grouping.var]], levels = unique(srt@meta.data[[grouping.var]]))
    }
    if (is.null(numerator)) {
      numerator <- levels(srt@meta.data[[grouping.var]])[1]
      warning("'numerator' is not specified. Use the first level in 'grouping.var': ", numerator, immediate. = TRUE)
    } else {
      if (!numerator %in% levels(srt@meta.data[, grouping.var])) {
        stop("'", numerator, "' is not an element of the '", grouping.var, "'")
      }
    }
    srt@meta.data[["grouping.var.use"]] <- srt@meta.data[[grouping.var]] == numerator
    add_dot <- FALSE
    exp_name <- paste0(numerator, "/", "other\n", exp_method, "(", data_nm, ")")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (is.null(group.by)) {
    srt@meta.data[["No.group.by"]] <- factor("")
    group.by <- "No.group.by"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(group.by[!group.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  group_elements <- unlist(lapply(srt@meta.data[, group.by, drop = FALSE], function(x) length(unique(x))))
  if (any(group_elements == 1) && exp_method == "zscore") {
    stop(
      "'zscore' cannot be applied to the group(s) consisting of one element: ",
      paste0(names(group_elements)[group_elements == 1], collapse = ",")
    )
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    stop("'group_palette' must be the same length as 'group.by'")
  }
  group_palette <- setNames(group_palette, nm = group.by)
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(split.by[!split.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- group.by
  }
  if (any(!feature_split_by %in% group.by)) {
    stop("feature_split_by must be a subset of group.by")
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(cell_annotation_palette, length(cell_annotation))
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(cell_annotation_palcolor, length(cell_annotation))
    }
    npal <- unique(c(length(cell_annotation_palette), length(cell_annotation_palcolor), length(cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]])))) {
      stop("cell_annotation: ", paste0(cell_annotation[!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(feature_annotation_palette, length(feature_annotation))
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(feature_annotation_palcolor, length(feature_annotation))
    }
    npal <- unique(c(length(feature_annotation_palette), length(feature_annotation_palcolor), length(feature_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features))) {
      stop("feature_annotation: ", paste0(feature_annotation[!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features)], collapse = ","), " is not in the meta data of the ", assay, " assay in the Seurat object.")
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(group.by))
  }
  if (length(height) == 1) {
    height <- rep(height, length(group.by))
  }
  if (length(width) >= 1) {
    names(width) <- group.by
  }
  if (length(height) >= 1) {
    names(height) <- group.by
  }

  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }

  if (is.null(cells)) {
    cells <- colnames(srt@assays[[1]])
  }
  if (all(!cells %in% colnames(srt@assays[[1]]))) {
    stop("No cells found.")
  }
  if (!all(cells %in% colnames(srt@assays[[1]]))) {
    warning("Some cells not found.", immediate. = TRUE)
  }
  cells <- intersect(cells, colnames(srt@assays[[1]]))

  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  index <- features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  features <- features[index]
  features_unique <- make.unique(features)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[index]
    names(feature_split) <- features_unique
  }

  cell_groups <- list()
  for (cell_group in group.by) {
    if (!is.factor(srt@meta.data[[cell_group]])) {
      srt@meta.data[[cell_group]] <- factor(srt@meta.data[[cell_group]], levels = unique(srt@meta.data[[cell_group]]))
    }
    cell_groups[[cell_group]] <- setNames(srt@meta.data[cells, cell_group], cells)
    cell_groups[[cell_group]] <- na.omit(cell_groups[[cell_group]])
    cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], levels = levels(cell_groups[[cell_group]])[levels(cell_groups[[cell_group]]) %in% cell_groups[[cell_group]]])

    if (!is.null(split.by)) {
      if (!is.factor(srt@meta.data[[split.by]])) {
        srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], levels = unique(srt@meta.data[[split.by]]))
      }
      levels <- apply(expand.grid(levels(srt@meta.data[[split.by]]), levels(cell_groups[[cell_group]])), 1, function(x) paste0(x[2:1], collapse = " : "))
      cell_groups[[cell_group]] <- setNames(paste0(cell_groups[[cell_group]][cells], " : ", srt@meta.data[cells, split.by]), cells)
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], levels = levels[levels %in% cell_groups[[cell_group]]])
    }
    if (!is.null(grouping.var)) {
      levels <- apply(expand.grid(c("TRUE", "FALSE"), levels(cell_groups[[cell_group]])), 1, function(x) paste0(x[2:1], collapse = " ; "))
      cell_groups[[cell_group]] <- setNames(paste0(cell_groups[[cell_group]][cells], " ; ", srt@meta.data[cells, "grouping.var.use"]), cells)
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], levels = levels[levels %in% cell_groups[[cell_group]]])
    }
  }

  gene <- features[features %in% rownames(srt@assays[[assay]])]
  gene_unique <- features_unique[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]

  mat_raw <- as.matrix(rbind(slot(srt@assays[[assay]], slot)[gene, cells, drop = FALSE], t(srt@meta.data[cells, meta, drop = FALSE])))[features, , drop = FALSE]
  rownames(mat_raw) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_raw, na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- colSums(slot(srt@assays[[assay]], "counts")[, colnames(mat_raw), drop = FALSE])
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", immediate. = TRUE)
        if (!is.null(grouping.var)) {
          exp_name <- paste0(numerator, "/", "other\n", exp_method, "(", slot, ")")
        } else {
          exp_name <- paste0(exp_method, "(", slot, ")")
        }
      }
    }
    mat_raw[gene_unique, ] <- t(t(mat_raw[gene_unique, , drop = FALSE]) / libsize_use * median(libsize_use))
  }

  mat_raw_list <- list()
  mat_perc_list <- list()
  for (cell_group in names(cell_groups)) {
    mat_tmp <- t(aggregate(t(mat_raw[features_unique, , drop = FALSE]), by = list(cell_groups[[cell_group]][colnames(mat_raw)]), FUN = aggregate_fun))
    colnames(mat_tmp) <- mat_tmp[1, , drop = FALSE]
    mat_tmp <- mat_tmp[-1, , drop = FALSE]
    class(mat_tmp) <- "numeric"
    mat_raw_list[[cell_group]] <- mat_tmp

    mat_perc <- t(aggregate(t(mat_raw[features_unique, , drop = FALSE]), by = list(cell_groups[[cell_group]][colnames(mat_raw)]), FUN = function(x) {
      sum(x > exp_cutoff) / length(x)
    }))
    colnames(mat_perc) <- mat_perc[1, , drop = FALSE]
    mat_perc <- mat_perc[-1, , drop = FALSE]
    class(mat_perc) <- "numeric"
    if (isTRUE(flip)) {
      mat_perc <- t(mat_perc)
    }
    mat_perc_list[[cell_group]] <- mat_perc
  }

  # data used to plot heatmap
  mat_list <- list()
  for (cell_group in group.by) {
    mat_tmp <- mat_raw_list[[cell_group]]
    if (is.null(grouping.var)) {
      mat_tmp <- matrix_process(mat_tmp, method = exp_method)
      mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)]), na.rm = TRUE) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
      mat_list[[cell_group]] <- mat_tmp
    } else {
      compare_groups <- strsplit(colnames(mat_tmp), " ; ")
      names_keep <- names(which(table(sapply(compare_groups, function(x) x[[1]])) == 2))
      group_keep <- which(sapply(compare_groups, function(x) x[[1]] %in% names_keep))
      group_TRUE <- intersect(group_keep, which(sapply(compare_groups, function(x) x[[2]]) == "TRUE"))
      group_FALSE <- intersect(group_keep, which(sapply(compare_groups, function(x) x[[2]]) == "FALSE"))
      mat_tmp <- log2(mat_tmp[, group_TRUE] / mat_tmp[, group_FALSE])
      colnames(mat_tmp) <- gsub(" ; .*", "", colnames(mat_tmp))
      mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)]), na.rm = TRUE) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
      mat_tmp[is.na(mat_tmp)] <- 0
      mat_list[[cell_group]] <- mat_tmp
      cell_groups[[cell_group]] <- factor(gsub(" ; .*", "", cell_groups[[cell_group]]), levels = unique(gsub(" ; .*", "", levels(cell_groups[[cell_group]]))))
    }
  }

  # data used to do clustering
  # if (length(feature_split_by) == 1) {
  #   mat_split <- mat_list[[feature_split_by]]
  # } else {
  #   # mat_split <- do.call(cbind, mat_raw_list[feature_split_by])
  #   # mat_split <- matrix_process(mat_split, method = exp_method)
  #   mat_split <- do.call(cbind, mat_list[feature_split_by])
  #   mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
  #   mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  # }
  mat_split <- do.call(cbind, mat_list[feature_split_by])

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(min(abs(quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 2) / 2
      colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    } else {
      b <- quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)
      colors <- colorRamp2(seq(b[1], b[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    }
  } else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }

  cell_metadata <- cbind.data.frame(
    data.frame(row.names = cells, cells = cells),
    cbind.data.frame(
      srt@meta.data[cells, c(group.by, intersect(cell_annotation, colnames(srt@meta.data))), drop = FALSE],
      t(srt@assays[[assay]]@data[intersect(cell_annotation, rownames(srt@assays[[assay]])) %||% integer(), cells, drop = FALSE])
    )
  )
  feature_metadata <- cbind.data.frame(
    data.frame(row.names = features_unique, features = features, features_uique = features_unique),
    srt@assays[[assay]]@meta.features[features, intersect(feature_annotation, colnames(srt@assays[[assay]]@meta.features)), drop = FALSE]
  )
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in% features[duplicated(features)]

  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)
  if (isTRUE(add_dot)) {
    lgd[["point"]] <- Legend(
      labels = paste0(seq(20, 100, length.out = 5), "%"),
      title = "Percent",
      type = "points",
      pch = 21,
      size = dot_size * seq(0.2, 1, length.out = 5), # unit(pi * grid_size^2 * seq(0.2, 1, length.out = 5), "cm"),
      grid_height = dot_size * seq(0.2, 1, length.out = 5) * 0.8,
      grid_width = dot_size,
      legend_gp = gpar(fill = "grey30"),
      border = FALSE,
      background = "transparent",
      direction = "vertical"
    )
  }

  ha_top_list <- NULL
  cluster_columns_list <- list()
  column_split_list <- list()
  for (i in seq_along(group.by)) {
    cell_group <- group.by[i]
    cluster_columns_list[[cell_group]] <- cluster_columns
    if (is.null(split.by)) {
      column_split_list[[cell_group]] <- NULL
    } else {
      column_split_list[[cell_group]] <- factor(gsub(" : .*", "", levels(cell_groups[[cell_group]])), levels = levels(srt@meta.data[[cell_group]]))
    }
    if (isTRUE(cluster_column_slices) && !is.null(split.by)) {
      if (!isTRUE(cluster_columns)) {
        if (nlevels(column_split_list[[cell_group]]) == 1) {
          stop("cluster_column_slices=TRUE can not be used when there is only one group.")
        }
        dend <- cluster_within_group(mat_list[[cell_group]], column_split_list[[cell_group]])
        cluster_columns_list[[cell_group]] <- dend
        column_split_list[[cell_group]] <- length(unique(column_split_list[[cell_group]]))
      }
    }
    if (cell_group != "No.group.by") {
      funbody <- paste0(
        "
        grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt@meta.data[[cell_group]]), collapse = "','"), "')"), ",palette = '", group_palette[i], "',palcolor=c(", paste0("'", paste0(group_palcolor[[i]], collapse = "','"), "'"), "))[nm]))
      "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      anno[[cell_group]] <- anno_block(
        align_to = split(seq_along(levels(cell_groups[[cell_group]])), gsub(pattern = " : .*", replacement = "", x = levels(cell_groups[[cell_group]]))),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = FALSE
      )
      ha_cell_group <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(flip, "row", "column"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "top", "left"), border = TRUE))
      ha_top_list[[cell_group]] <- ha_cell_group
      lgd[[cell_group]] <- Legend(
        title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      )
    }

    if (!is.null(split.by)) {
      funbody <- paste0(
        "
      grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt@meta.data[[split.by]]), collapse = "','"), "')"), ",palette = '", cell_split_palette, "',palcolor=c(", paste0("'", paste0(unlist(cell_split_palcolor), collapse = "','"), "'"), "))[nm]))
    "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      anno[[split.by]] <- anno_block(
        align_to = split(seq_along(levels(cell_groups[[cell_group]])), gsub(pattern = ".* : ", replacement = "", x = levels(cell_groups[[cell_group]]))),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = i == 1
      )
      ha_split_by <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(flip, "row", "column"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "top", "left"), border = TRUE))
      if (is.null(ha_top_list[[cell_group]])) {
        ha_top_list[[cell_group]] <- ha_split_by
      } else {
        ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_split_by)
      }
    }
  }
  if (!is.null(split.by)) {
    lgd[[split.by]] <- Legend(
      title = split.by, labels = levels(srt@meta.data[[split.by]]),
      legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[split.by]]), palette = cell_split_palette, palcolor = cell_split_palcolor)), border = TRUE
    )
  }

  if (!is.null(cell_annotation)) {
    subplots_list <- list()
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (cell_group in group.by) {
          subplots <- CellStatPlot(srt,
            flip = flip,
            cells = names(cell_groups[[cell_group]]), plot_type = "pie",
            group.by = cell_group, stat.by = cellan, split.by = split.by,
            palette = palette, palcolor = palcolor,
            individual = TRUE, combine = FALSE
          )
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", cell_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(cell_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(cell_group, x[1], ""), collapse = ":")
            }
          })

          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "row", "column"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = cell_group == group.by[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
          }
        }
        lgd[[cellan]] <- Legend(
          title = cellan, labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        for (cell_group in group.by) {
          subplots <- FeatureStatPlot(srt,
            assay = assay, slot = "data", flip = flip,
            stat.by = cellan, cells = names(cell_groups[[cell_group]]),
            group.by = cell_group, split.by = split.by,
            palette = palette, palcolor = palcolor,
            fill.by = "group", same.y.lims = TRUE,
            individual = TRUE, combine = FALSE
          )
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", cell_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(cellan, cell_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(cellan, cell_group, x[1], ""), collapse = ":")
            }
          })
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "row", "column"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = cell_group == group.by[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
          }
        }
        # lgd[[cellan]] <- Legend(
        #   title = cellan, labels = levels(cell_anno),
        #   legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        # )
      }
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat_split) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    } else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- setNames(rep(1, nrow(mat_split)), rownames(mat_split))
      } else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_R("Mfuzz"), error = identity)
          if (inherits(status, "error")) {
            warning("The Mfuzz package was not found. Switch split_method to 'kmeans'", immediate. = TRUE)
            split_method <- "kmeans"
          } else {
            require("Mfuzz", quietly = TRUE)
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            eset <- new("ExpressionSet", exprs = mat_split_tmp)
            eset <- Mfuzz::standardise(eset)
            min_fuzzification <- Mfuzz::mestimate(eset)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            } else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- Mfuzz::mfuzz(eset, c = n_split, m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
            }
            if (isTRUE(show_fuzzification)) {
              message("fuzzification: ", fuzzification)
            }
            # mfuzz.plot(eset, cl,new.window = FALSE)
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, iter.max = 1e4, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(stats::dist(mat_split))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      groupmean <- aggregate(t(mat_split), by = list(unlist(lapply(cell_groups[feature_split_by], levels))), mean)
      maxgroup <- groupmean[, 1][apply(groupmean[, names(row_split)], 2, which.max)]
      maxgroup <- factor(maxgroup, levels = levels(unlist(cell_groups[feature_split_by])))
      df <- data.frame(row_split = row_split, order_by = maxgroup)
      df_order <- aggregate(df[["order_by"]], by = list(df[, "row_split"]), FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
      df_order[, "row_split"] <- df_order[, "Group.1"]
      df_order[["order_by"]] <- as.numeric(factor(df_order[["x"]], levels = levels(maxgroup)))
      df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), , drop = FALSE]
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))
    }
  } else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat_split)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  } else {
    feature_metadata[["feature_split"]] <- NA
  }

  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- cluster_within_group(t(mat_split), row_split_raw)
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0(
      "
      grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(row_split_raw), collapse = "','"), "')"), ",palette = '", feature_split_palette, "',palcolor=c(", paste0("'", paste0(unlist(feature_split_palcolor), collapse = "','"), "'"), "))[nm]))
    "
    )
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())
    ha_clusters <- HeatmapAnnotation(
      features_split = anno_block(
        align_to = split(seq_along(row_split_raw), row_split_raw),
        panel_fun = getFunction("panel_fun", where = environment()),
        width = unit(0.1, "in"),
        height = unit(0.1, "in"),
        show_name = FALSE,
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"),
      border = TRUE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(
      title = "Cluster", labels = intersect(levels(row_split_raw), row_split_raw),
      legend_gp = gpar(fill = palette_scp(intersect(levels(row_split_raw), row_split_raw), type = "discrete", palette = feature_split_palette, palcolor = feature_split_palcolor, matched = TRUE)), border = TRUE
    )
  }

  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- do.call(cbind, mat_list[cluster_features_by])
    if (is.null(row_split)) {
      dend <- as.dendrogram(hclust(dist(mat_cluster)))
      dend_ordered <- reorder(dend, wts = colMeans(mat_cluster), agglo.FUN = mean)
      cluster_rows <- dend_ordered
    } else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }

  cell_group <- group.by[1]
  ht_args <- list(
    matrix = mat_list[[cell_group]],
    col = colors,
    row_split = row_split,
    column_split = column_split_list[[cell_group]],
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns_list[[cell_group]],
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_column_slices,
    use_raster = TRUE
  )
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  ht_list <- do.call(Heatmap, args = ht_args)
  features_ordered <- rownames(mat_list[[1]])[unlist(suppressWarnings(row_order(ht_list)))]
  feature_metadata[["index"]] <- setNames(object = seq_along(features_ordered), nm = features_ordered)[rownames(feature_metadata)]

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered) / nlabel) / 2)
        index_to <- length(features_ordered)
        index <- unique(round(seq(from = index_from, to = index_to, length.out = nlabel)))
      } else {
        index <- seq_along(features_ordered)
      }
    } else {
      index <- NULL
    }
  } else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the features"), immediate. = TRUE)
    }
  }
  if (length(index) > 0) {
    ha_mark <- HeatmapAnnotation(
      gene = anno_mark(
        at = which(rownames(feature_metadata) %in% features_ordered[index]),
        labels = feature_metadata[which(rownames(feature_metadata) %in% features_ordered[index]), "features"],
        side = ifelse(flip, "top", "left"),
        labels_gp = gpar(fontsize = label_size, col = label_color),
        link_gp = gpar(fontsize = label_size, col = label_color),
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"), show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    } else {
      ha_left <- c(ha_mark, ha_left)
    }
  }

  ha_right <- NULL
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, FALSE))
        } else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = as.character(featan_values),
          col = palette_scp(featan_values, palette = palette, palcolor = palcolor),
          na_col = "transparent",
          which = ifelse(flip, "column", "row"),
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, labels = levels(featan_values),
          legend_gp = gpar(fill = palette_scp(featan_values, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        col_fun <- colorRamp2(
          breaks = seq(min(featan_values, na.rm = TRUE), max(featan_values, na.rm = TRUE), length = 100),
          colors = palette_scp(palette = palette, palcolor = palcolor)
        )
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = featan_values,
          col = col_fun,
          na_col = "transparent",
          which = ifelse(flip, "column", "row"),
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, col_fun = col_fun, border = TRUE
        )
      }
    }
  }

  enrichment <- heatmap_enrichment(
    geneID = feature_metadata[["features"]], geneID_groups = feature_metadata[["feature_split"]],
    feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor,
    ha_right = ha_right, flip = flip,
    anno_terms = anno_terms, anno_keys = anno_keys, anno_features = anno_features,
    terms_width = terms_width, terms_fontsize = terms_fontsize,
    keys_width = keys_width, keys_fontsize = keys_fontsize,
    features_width = features_width, features_fontsize = features_fontsize,
    IDtype = IDtype, species = species, db_update = db_update, db_version = db_version, db_combine = db_combine, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror,
    db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize,
    GO_simplify = GO_simplify, GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff,
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, min_word_length = min_word_length,
    exclude_words = exclude_words
  )
  res <- enrichment$res
  ha_right <- enrichment$ha_right

  ht_list <- NULL
  vlnplots_list <- NULL
  x_nm_list <- NULL
  y_nm_list <- NULL
  if (fill.by == "group") {
    palette <- group_palette
    palcolor <- group_palcolor
  } else {
    palette <- feature_annotation_palette
    palcolor <- feature_annotation_palcolor
  }
  for (cell_group in group.by) {
    if (cell_group == group.by[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (cell_group == group.by[length(group.by)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }
    if (isTRUE(add_violin)) {
      vlnplots <- FeatureStatPlot(srt,
        assay = assay, slot = "data", flip = flip,
        stat.by = rownames(mat_list[[cell_group]]),
        cells = names(cell_groups[[cell_group]]),
        group.by = cell_group, split.by = split.by,
        palette = fill_palette, palcolor = fill_palcolor,
        fill.by = fill.by, same.y.lims = TRUE,
        individual = TRUE, combine = FALSE
      )
      lgd[["ht"]] <- NULL
      # legend <- get_legend(vlnplots[[1]])
      # funbody <- paste0(
      #   "
      #         g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
      #         grid.draw(g)
      #         "
      # )
      # funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      # eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
      # graphics = list(
      #   " " = function(x, y, w, h) {
      #     grid_draw(legend,x = x,y = y,width = w,height = h)
      # })
      # lgd[["ht"]] <- Legend(title = " ", at = names(graphics), graphics = graphics, border = TRUE)

      for (nm in names(vlnplots)) {
        gtable <- as_grob(vlnplots[[nm]] + facet_null() + theme_void() + theme(legend.position = "none"))
        gtable$name <- paste0(cell_group, "-", nm)
        vlnplots[[nm]] <- gtable
      }
      vlnplots_list[[paste0("heatmap_group:", cell_group)]] <- vlnplots
      x_nm <- rownames(mat_list[[cell_group]])
      x_nm_list[[paste0("heatmap_group:", cell_group)]] <- x_nm
      y_nm <- sapply(strsplit(levels(cell_groups[[cell_group]]), " : "), function(x) {
        if (length(x) == 2) {
          paste0(c(cell_group, x[1], x[2]), collapse = ":")
        } else {
          paste0(c(cell_group, x[1], ""), collapse = ":")
        }
      })
      y_nm_list[[paste0("heatmap_group:", cell_group)]] <- y_nm

      # popViewport()
      # grid.draw(roundrectGrob())
      # groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:' ,cell_group)]],
      #                          x_nm =  x_nm_list[[paste0('heatmap_group:', cell_group)]],
      #                          y_nm= y_nm_list[[paste0('heatmap_group:',cell_group)]],
      #                          x = 1:4,y = 1:4);
      # grid_draw(groblist,
      #   x = unit(c(0.33, 0.67, 0.33, 0.67), "npc"), y = unit(c(0.33, 0.33, 0.67, 0.67), "npc"),
      #   width = rep(unit(1, "in"), 4), height = rep(unit(1, "in"), 4)
      # )
    }

    funbody <- paste0(
      if (isTRUE(add_dot) || isTRUE(add_violin)) {
        "grid.rect(x, y,
          width = width, height = height,
          gp = gpar(col = 'white', lwd = 1, fill = 'white')
        );"
      },
      if (isTRUE(add_bg)) {
        paste0("
        grid.rect(x, y,
          width = width, height = height,
          gp = gpar(col = fill, lwd = 1, fill = adjcolors(fill, ", bg_alpha, "))
        );
        ")
      },
      if (isTRUE(add_reticle)) {
        paste0("
        ind_mat = restore_matrix(j, i, x, y);
        ind_top = ind_mat[1,];
        ind_left = ind_mat[,1];
        for(col in seq_len(ncol(ind_mat))){
          grid.lines(x = unit(rep(x[ind_top[col]],each=2),'npc'),y = unit(c(0,1),'npc'),gp = gpar(col = '", reticle_color, "', lwd = 1.5));
        };
        for(row in seq_len(nrow(ind_mat))){
          grid.lines(x = unit(c(0,1),'npc'),y = unit(rep(y[ind_left[row]],each=2),'npc'),gp = gpar(col = '", reticle_color, "', lwd = 1.5));
        };
        ")
      },
      if (isTRUE(add_dot)) {
        paste0("perc <- pindex(mat_perc_list[['", cell_group, "']]", ", i, j);
        grid.points(x, y,
          pch = 21,
          size = dot_size*perc,
          gp = gpar(col = 'black', lwd = 1, fill = fill)
        );
        ")
      },
      if (isTRUE(add_violin)) {
        if (isTRUE(flip)) {
          paste0("
        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '", cell_group, "')]],
               x_nm =  x_nm_list[[paste0('heatmap_group:', '", cell_group, "')]],
               y_nm= y_nm_list[[paste0('heatmap_group:', '", cell_group, "')]],
               x = j,y = i);
        grid_draw(groblist, x = x, y = y, width = width, height = height);
        ")
        } else {
          paste0("
        groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0('heatmap_group:', '", cell_group, "')]],
               x_nm =  x_nm_list[[paste0('heatmap_group:', '", cell_group, "')]],
               y_nm= y_nm_list[[paste0('heatmap_group:', '", cell_group, "')]],
               x = i,y = j);
        grid_draw(groblist, x = x, y = y, width = width, height = height);
        ")
        }
      }
    )
    # groblist <- extractgrobs(vlnplots = vlnplots_list[[paste0("heatmap_group:", cell_group)]],
    #                          x_nm =  x_nm_list[[paste0("heatmap_group:", cell_group)]],
    #                          y_nm= y_nm_list[[paste0("heatmap_group:", cell_group)]],
    #                          x = 1:3,y = 1:3)
    #
    #
    #
    # gtable <- groblist[[paste0(
    #   x_nm_list[[paste0("heatmap_group:", cell_group)]][1],
    #   ":",
    #   y_nm_list[[paste0("heatmap_group:", cell_group)]][1]
    # )]]
    # ind_mat = restore_matrix(j, i, x, y);
    # ind = unique(ind_mat);
    # grid_draw(groblist,x=x[ind],y= y[ind]);

    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("layer_fun <- function(j, i, x, y, width, height, fill) {", funbody, "}", sep = "")), envir = environment())
    ht_args <- list(
      name = cell_group,
      matrix = if (flip) t(mat_list[[cell_group]]) else mat_list[[cell_group]],
      col = colors,
      layer_fun = getFunction("layer_fun", where = environment()),
      row_title = if (flip) ifelse(cell_group != "No.group.by", cell_group, "") else character(0),
      row_title_side = row_title_side,
      column_title = if (flip) character(0) else ifelse(cell_group != "No.group.by", cell_group, ""),
      column_title_side = column_title_side,
      row_title_rot = row_title_rot,
      column_title_rot = column_title_rot,
      row_split = if (flip) column_split_list[[cell_group]] else row_split,
      column_split = if (flip) row_split else column_split_list[[cell_group]],
      cluster_rows = if (flip) cluster_columns_list[[cell_group]] else cluster_rows,
      cluster_columns = if (flip) cluster_rows else cluster_columns_list[[cell_group]],
      cluster_row_slices = if (flip) cluster_column_slices else cluster_row_slices,
      cluster_column_slices = if (flip) cluster_row_slices else cluster_column_slices,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      top_annotation = if (flip) left_annotation else ha_top_list[[cell_group]],
      left_annotation = if (flip) ha_top_list[[cell_group]] else left_annotation,
      bottom_annotation = if (flip) right_annotation else NULL,
      right_annotation = if (flip) NULL else right_annotation,
      show_heatmap_legend = FALSE,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      width = if (is.numeric(width[cell_group])) unit(width[cell_group], units = units) else NULL,
      height = if (is.numeric(height[cell_group])) unit(height[cell_group], units = units) else NULL
    )
    if (any(names(ht_params) %in% names(ht_args))) {
      warning("ht_params: ", paste0(intersect(names(ht_params), names(ht_args)), collapse = ","), " were duplicated and will not be used.", immediate. = TRUE)
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(Heatmap, args = ht_args)
      } else {
        ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
      }
    } else {
      ht_list <- ht_list + do.call(Heatmap, args = ht_args)
    }
  }

  if ((!is.null(row_split) && length(index) > 0) || any(c(anno_terms, anno_keys, anno_features)) || !is.null(width) || !is.null(height)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width, height = height, units = units,
    ha_top_list = ha_top_list, ha_left = ha_left, ha_right = ha_right,
    ht_list = ht_list, legend_list = lgd, flip = flip
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  # cat("width:", width, "\n")
  # cat("height:", height, "\n")
  # cat("width_sum:", width_sum, "\n")
  # cat("height_sum:", height_sum, "\n")

  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(
      width = width, width_sum = width_sum, height = height, height_sum = height_sum, units = units,
      ht_list = ht_list, legend_list = lgd
    )
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    # cat("ht_width:", ht_width, "\n")
    # cat("ht_height:", ht_height, "\n")

    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  } else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  }

  if (isTRUE(fix)) {
    p <- panel_fix_single(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
  } else {
    p <- wrap_plots(gTree)
  }

  return(list(
    plot = p,
    matrix_list = mat_list,
    feature_split = feature_split,
    cell_metadata = cell_metadata,
    feature_metadata = feature_metadata,
    enrichment = res
  ))

  return(p)
}

#' FeatureHeatmap
#'
#' @param srt
#' @param features
#' @param feature_split
#' @param cluster_rows
#' @param group.by
#' @param cluster_columns
#' @param max_cells
#' @param slot
#' @param assay
#' @param exp_method
#' @param n_split
#' @param feature_split_by
#' @param split_method
#' @param decreasing
#' @param lib_normalize
#' @param libsize
#' @param anno_keys
#' @param anno_features
#' @param IDtype
#' @param species
#' @param db_update
#' @param db_version
#' @param Ensembl_version
#' @param mirror
#' @param db
#' @param TERM2GENE
#' @param TERM2NAME
#' @param minGSSize
#' @param maxGSSize
#' @param GO_simplify
#' @param GO_simplify_cutoff
#' @param simplify_method
#' @param simplify_similarityCutoff
#' @param pvalueCutoff
#' @param padjustCutoff
#' @param topWord
#' @param min_word_length
#' @param exclude_words
#' @param nlabel
#' @param features_label
#' @param label_size
#' @param label_color
#' @param heatmap_palette
#' @param group_palette
#' @param feature_split_palette
#' @param cell_annotation
#' @param cell_annotation_palette
#' @param cell_annotation_palcolor
#' @param feature_annotation
#' @param feature_annotation_palette
#' @param feature_annotation_palcolor
#' @param use_raster
#' @param height
#' @param width
#' @param units
#' @param seed
#' @param cells
#' @param split.by
#' @param cell_order
#' @param border
#' @param flip
#' @param cluster_features_by
#' @param cluster_row_slices
#' @param cluster_column_slices
#' @param show_row_names
#' @param show_column_names
#' @param row_names_side
#' @param column_names_side
#' @param row_names_rot
#' @param column_names_rot
#' @param row_title_side
#' @param column_title_side
#' @param row_title_rot
#' @param column_title_rot
#' @param anno_terms
#' @param terms_width
#' @param terms_fontsize
#' @param keys_width
#' @param keys_fontsize
#' @param features_width
#' @param features_fontsize
#' @param convert_species
#' @param topTerm
#' @param show_termid
#' @param heatmap_palcolor
#' @param group_palcolor
#' @param cell_split_palette
#' @param cell_split_palcolor
#' @param feature_split_palcolor
#' @param cell_annotation_params
#' @param feature_annotation_params
#' @param raster_device
#' @param ht_params
#' @param limits
#' @param raster_by_magick
#' @param fuzzification
#' @param show_fuzzification
#' @param db_combine
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' de_filter <- filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' ht1 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, group.by = "CellType",
#'   split.by = "Phase", cell_split_palette = "Dark2",
#' )
#' ht1$plot
#' panel_fix(ht1$plot, height = 4, width = 6, raster = TRUE, dpi = 50)
#'
#' ht2 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, group.by = c("CellType", "SubCellType"), n_split = 4,
#'   cluster_rows = TRUE, cluster_row_slices = TRUE, cluster_columns = TRUE, cluster_column_slices = TRUE,
#'   ht_params = list(row_gap = grid::unit(0, "mm"))
#' )
#' ht2$plot
#'
#' ht3 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, feature_split = de_filter$group1, group.by = "CellType",
#'   species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "SP"))
#' ht4 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, n_split = 4, group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   cell_annotation = c("Phase", "G2M_score"), cell_annotation_palette = c("Dark2", "Purples")
#' )
#' ht4$plot
#'
#' ht5 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, n_split = 4, group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   cell_annotation = c("Phase", "G2M_score"), cell_annotation_palette = c("Dark2", "Purples"),
#'   flip = TRUE, column_title_rot = 45
#' )
#' ht5$plot
#'
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub, assay_X = "RNA", group_by = "SubCellType",
#'   linear_reduction = "PCA", nonlinear_reduction = "UMAP", infer_pseudotime = TRUE, root_group = "Ductal"
#' )
#' ht6 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, nlabel = 10,
#'   cell_order = names(sort(pancreas_sub$dpt_pseudotime)),
#'   cell_annotation = c("CellType", "dpt_pseudotime"),
#'   cell_annotation_palette = c("Paired", "cividis")
#' )
#' ht6$plot
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap Legend HeatmapAnnotation anno_empty anno_mark anno_simple anno_textbox draw decorate_heatmap_body width.HeatmapAnnotation height.HeatmapAnnotation width.Legends height.Legends cluster_within_group decorate_annotation row_order %v%
#' @importFrom grid gpar grid.grabExpr grid.text
#' @importFrom gtable gtable_add_padding
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat GetAssayData
#' @importFrom stats hclust order.dendrogram as.dendrogram reorder
#' @importFrom dplyr %>% filter group_by arrange desc across mutate distinct n .data "%>%"
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @export
FeatureHeatmap <- function(srt, features = NULL, cells = NULL, group.by = NULL, split.by = NULL, max_cells = 100, cell_order = NULL, border = TRUE, flip = FALSE,
                           slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), limits = NULL, lib_normalize = identical(slot, "counts"), libsize = NULL,
                           feature_split = NULL, feature_split_by = NULL, n_split = NULL, split_method = c("kmeans", "hclust", "mfuzz"), decreasing = FALSE, fuzzification = NULL, show_fuzzification = FALSE,
                           cluster_features_by = NULL, cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                           row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                           anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                           terms_width = unit(4, "in"), terms_fontsize = 8,
                           keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                           features_width = unit(2, "in"), features_fontsize = c(6, 10),
                           IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                           db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                           GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, min_word_length = 3,
                           exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "organization", "system", "regulation", "positive", "negative", "response", "process"),
                           nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                           heatmap_palette = "RdBu", heatmap_palcolor = NULL, group_palette = "Paired", group_palcolor = NULL,
                           cell_split_palette = "jama", cell_split_palcolor = NULL, feature_split_palette = "jama", feature_split_palcolor = NULL,
                           cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = list(),
                           feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = list(),
                           use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, height = NULL, width = NULL, units = "inch",
                           seed = 11, ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }

  split_method <- match.arg(split_method)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), "(", data_nm, ")")
  } else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (is.null(group.by)) {
    srt@meta.data[["No.group.by"]] <- factor("")
    group.by <- "No.group.by"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(group.by[!group.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    stop("'group_palette' must be the same length as 'group.by'")
  }
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(split.by[!split.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- group.by
  }
  if (any(!feature_split_by %in% group.by)) {
    stop("feature_split_by must be a subset of group.by")
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(cell_annotation_palette, length(cell_annotation))
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(cell_annotation_palcolor, length(cell_annotation))
    }
    npal <- unique(c(length(cell_annotation_palette), length(cell_annotation_palcolor), length(cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]])))) {
      stop("cell_annotation: ", paste0(cell_annotation[!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(feature_annotation_palette, length(feature_annotation))
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(feature_annotation_palcolor, length(feature_annotation))
    }
    npal <- unique(c(length(feature_annotation_palette), length(feature_annotation_palcolor), length(feature_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features))) {
      stop("feature_annotation: ", paste0(feature_annotation[!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features)], collapse = ","), " is not in the meta data of the ", assay, " assay in the Seurat object.")
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(group.by))
  }
  if (length(height) == 1) {
    height <- rep(height, length(group.by))
  }
  if (length(width) >= 1) {
    names(width) <- group.by
  }
  if (length(height) >= 1) {
    names(height) <- group.by
  }

  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }

  if (is.null(cells)) {
    cells <- colnames(srt@assays[[1]])
  }
  if (all(!cells %in% colnames(srt@assays[[1]]))) {
    stop("No cells found.")
  }
  if (!all(cells %in% colnames(srt@assays[[1]]))) {
    warning("Some cells not found.", immediate. = TRUE)
  }
  cells <- intersect(cells, colnames(srt@assays[[1]]))

  if (is.null(features)) {
    features <- VariableFeatures(srt, assay = assay)
  }
  index <- features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data))
  features <- features[index]
  features_unique <- make.unique(features)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[index]
    names(feature_split) <- features_unique
  }

  cell_groups <- list()
  for (cell_group in group.by) {
    if (!is.factor(srt@meta.data[[cell_group]])) {
      srt@meta.data[[cell_group]] <- factor(srt@meta.data[[cell_group]], levels = unique(srt@meta.data[[cell_group]]))
    }
    if (is.null(split.by)) {
      cell_groups[[cell_group]] <- unlist(lapply(levels(srt@meta.data[[cell_group]]), function(x) {
        cells_sub <- colnames(srt@assays[[1]])[which(srt@meta.data[[cell_group]] == x)]
        cells_sub <- intersect(cells, cells_sub)
        size <- ifelse(length(cells_sub) > max_cells, max_cells, length(cells_sub))
        cells_sample <- sample(cells_sub, size)
        out <- setNames(rep(x, size), cells_sample)
        return(out)
      }), use.names = TRUE)
      levels <- levels(srt@meta.data[[cell_group]])
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], levels = levels[levels %in% cell_groups[[cell_group]]])
    } else {
      if (!is.factor(srt@meta.data[[split.by]])) {
        srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], levels = unique(srt@meta.data[[split.by]]))
      }
      cell_groups[[cell_group]] <- unlist(lapply(levels(srt@meta.data[[cell_group]]), function(x) {
        cells_sub <- colnames(srt@assays[[1]])[srt@meta.data[[cell_group]] == x]
        cells_sub <- intersect(cells, cells_sub)
        cells_tmp <- NULL
        for (sp in levels(srt@meta.data[[split.by]])) {
          cells_sp <- cells_sub[srt@meta.data[cells_sub, split.by] == sp]
          size <- ifelse(length(cells_sp) > max_cells, max_cells, length(cells_sp))
          cells_sample <- sample(cells_sp, size)
          cells_tmp <- c(cells_tmp, setNames(rep(paste0(x, " : ", sp), size), cells_sample))
        }
        size <- ifelse(length(cells_tmp) > max_cells, max_cells, length(cells_tmp))
        out <- sample(cells_tmp, size)
        return(out)
      }), use.names = TRUE)
      levels <- apply(expand.grid(levels(srt@meta.data[[split.by]]), levels(srt@meta.data[[cell_group]])), 1, function(x) paste0(x[2:1], collapse = " : "))
      cell_groups[[cell_group]] <- factor(cell_groups[[cell_group]], levels = levels[levels %in% cell_groups[[cell_group]]])
    }
    if (!is.null(cell_order)) {
      cell_groups[[cell_group]] <- cell_groups[[cell_group]][intersect(cell_order, names(cell_groups[[cell_group]]))]
    }
  }

  gene <- features[features %in% rownames(srt@assays[[assay]])]
  gene_unique <- features_unique[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  all_cells <- unique(unlist(lapply(cell_groups, names)))
  mat_raw <- as.matrix(rbind(slot(srt@assays[[assay]], slot)[gene, all_cells, drop = FALSE], t(srt@meta.data[all_cells, meta, drop = FALSE])))[features, , drop = FALSE]
  rownames(mat_raw) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_raw, na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- colSums(slot(srt@assays[[assay]], "counts")[, colnames(mat_raw), drop = FALSE])
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", immediate. = TRUE)
      }
    }
    mat_raw[gene_unique, ] <- t(t(mat_raw[gene_unique, , drop = FALSE]) / libsize_use * median(libsize_use))
  }

  # data used to plot heatmap
  mat_list <- list()
  for (cell_group in group.by) {
    mat_tmp <- mat_raw[, names(cell_groups[[cell_group]])]
    mat_tmp <- matrix_process(mat_tmp, method = exp_method)
    mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)]), na.rm = TRUE) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
    mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
    mat_list[[cell_group]] <- mat_tmp
  }

  # data used to do clustering
  # if (length(feature_split_by) == 1) {
  #   mat_split <- mat_list[[feature_split_by]]
  # } else {
  #   # mat_split <- mat_list[, unlist(lapply(cell_groups[feature_split_by], names))]
  #   # mat_split <- matrix_process(mat_split, method = exp_method)
  #   mat_split <- do.call(cbind, mat_list[feature_split_by])
  #   mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
  #   mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  # }
  mat_split <- do.call(cbind, mat_list[feature_split_by])

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(min(abs(quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 2) / 2
      colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    } else {
      b <- quantile(do.call(cbind, mat_list), c(0.01, 0.99), na.rm = TRUE)
      colors <- colorRamp2(seq(b[1], b[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    }
  } else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }

  cell_metadata <- cbind.data.frame(
    data.frame(row.names = colnames(mat_raw), cells = colnames(mat_raw)),
    cbind.data.frame(
      srt@meta.data[colnames(mat_raw), c(group.by, intersect(cell_annotation, colnames(srt@meta.data))), drop = FALSE],
      t(srt@assays[[assay]]@data[intersect(cell_annotation, rownames(srt@assays[[assay]])) %||% integer(), colnames(mat_raw), drop = FALSE])
    )
  )
  feature_metadata <- cbind.data.frame(
    data.frame(row.names = features_unique, features = features, features_uique = features_unique),
    srt@assays[[assay]]@meta.features[features, c(feature_annotation), drop = FALSE]
  )
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in% features[duplicated(features)]

  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_top_list <- NULL
  cluster_columns_list <- list()
  column_split_list <- list()
  for (i in seq_along(group.by)) {
    cell_group <- group.by[i]
    cluster_columns_list[[cell_group]] <- cluster_columns
    column_split_list[[cell_group]] <- cell_groups[[cell_group]]
    if (isTRUE(cluster_column_slices)) {
      if (!isTRUE(cluster_columns)) {
        if (nlevels(column_split_list[[cell_group]]) == 1) {
          stop("cluster_column_slices=TRUE can not be used when there is only one group.")
        }
        dend <- cluster_within_group(mat_list[[cell_group]], column_split_list[[cell_group]])
        cluster_columns_list[[cell_group]] <- dend
        column_split_list[[cell_group]] <- length(unique(column_split_list[[cell_group]]))
      }
    }
    if (cell_group != "No.group.by") {
      funbody <- paste0(
        "
        grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt@meta.data[[cell_group]]), collapse = "','"), "')"), ",palette = '", group_palette[i], "',palcolor=c(", paste0("'", paste0(group_palcolor[[i]], collapse = "','"), "'"), "))[nm]))
      "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      anno[[cell_group]] <- anno_block(
        align_to = split(seq_along(cell_groups[[cell_group]]), gsub(pattern = " : .*", replacement = "", x = cell_groups[[cell_group]])),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = FALSE
      )
      ha_cell_group <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(flip, "row", "column"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "top", "left"), border = TRUE))
      ha_top_list[[cell_group]] <- ha_cell_group
      lgd[[cell_group]] <- Legend(
        title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      )
    }

    if (!is.null(split.by)) {
      funbody <- paste0(
        "
      grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt@meta.data[[split.by]]), collapse = "','"), "')"), ",palette = '", cell_split_palette, "',palcolor=c(", paste0("'", paste0(unlist(cell_split_palcolor), collapse = "','"), "'"), "))[nm]))
    "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      anno[[split.by]] <- anno_block(
        align_to = split(seq_along(cell_groups[[cell_group]]), gsub(pattern = ".* : ", replacement = "", x = cell_groups[[cell_group]])),
        panel_fun = getFunction("panel_fun", where = environment()),
        which = ifelse(flip, "row", "column"),
        show_name = i == 1
      )
      ha_split_by <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(flip, "row", "column"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "top", "left"), border = TRUE))
      if (is.null(ha_top_list[[cell_group]])) {
        ha_top_list[[cell_group]] <- ha_split_by
      } else {
        ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_split_by)
      }
    }
  }
  if (!is.null(split.by)) {
    lgd[[split.by]] <- Legend(
      title = split.by, labels = levels(srt@meta.data[[split.by]]),
      legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[split.by]]), palette = cell_split_palette, palcolor = cell_split_palcolor)), border = TRUE
    )
  }

  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (cell_group in group.by) {
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = as.character(cell_anno[names(cell_groups[[cell_group]])]),
            col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = cell_group == group.by[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
          }
        }
        lgd[[cellan]] <- Legend(
          title = cellan, labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        col_fun <- colorRamp2(
          breaks = seq(min(cell_anno, na.rm = TRUE), max(cell_anno, na.rm = TRUE), length = 100),
          colors = palette_scp(palette = palette, palcolor = palcolor)
        )
        for (cell_group in group.by) {
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = cell_anno[names(cell_groups[[cell_group]])],
            col = col_fun,
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = cell_group == group.by[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[cell_group]])) {
            ha_top_list[[cell_group]] <- ha_top
          } else {
            ha_top_list[[cell_group]] <- c(ha_top_list[[cell_group]], ha_top)
          }
        }
        lgd[[cellan]] <- Legend(
          title = cellan, col_fun = col_fun, border = TRUE
        )
      }
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat_split) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    } else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- setNames(rep(1, nrow(mat_split)), rownames(mat_split))
      } else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_R("Mfuzz"), error = identity)
          if (inherits(status, "error")) {
            warning("The Mfuzz package was not found. Switch split_method to 'kmeans'", immediate. = TRUE)
            split_method <- "kmeans"
          } else {
            require("Mfuzz", quietly = TRUE)
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            eset <- new("ExpressionSet", exprs = mat_split_tmp)
            eset <- Mfuzz::standardise(eset)
            min_fuzzification <- Mfuzz::mestimate(eset)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            } else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- Mfuzz::mfuzz(eset, c = n_split, m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
            }
            if (isTRUE(show_fuzzification)) {
              message("fuzzification: ", fuzzification)
            }
            # mfuzz.plot(eset, cl,new.window = FALSE)
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, iter.max = 1e4, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(stats::dist(mat_split))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      groupmean <- aggregate(t(mat_split), by = list(unlist(cell_groups[feature_split_by])), mean)
      maxgroup <- groupmean[, 1][apply(groupmean[, names(row_split)], 2, which.max)]
      df <- data.frame(row_split = row_split, order_by = maxgroup)
      df_order <- aggregate(df[["order_by"]], by = list(df[, "row_split"]), FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
      df_order[, "row_split"] <- df_order[, "Group.1"]
      df_order[["order_by"]] <- as.numeric(factor(df_order[["x"]], levels = levels(maxgroup)))
      df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), , drop = FALSE]
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))
    }
  } else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat_split)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  } else {
    feature_metadata[["feature_split"]] <- NA
  }

  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- cluster_within_group(t(mat_split), row_split_raw)
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0(
      "
      grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(row_split_raw), collapse = "','"), "')"), ",palette = '", feature_split_palette, "',palcolor=c(", paste0("'", paste0(unlist(feature_split_palcolor), collapse = "','"), "'"), "))[nm]))
    "
    )
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())
    ha_clusters <- HeatmapAnnotation(
      features_split = anno_block(
        align_to = split(seq_along(row_split_raw), row_split_raw),
        panel_fun = getFunction("panel_fun", where = environment()),
        width = unit(0.1, "in"),
        height = unit(0.1, "in"),
        show_name = FALSE,
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"),
      border = TRUE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(
      title = "Cluster", labels = intersect(levels(row_split_raw), row_split_raw),
      legend_gp = gpar(fill = palette_scp(intersect(levels(row_split_raw), row_split_raw), type = "discrete", palette = feature_split_palette, palcolor = feature_split_palcolor, matched = TRUE)), border = TRUE
    )
  }

  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- do.call(cbind, mat_list[cluster_features_by])
    if (is.null(row_split)) {
      dend <- as.dendrogram(hclust(dist(mat_cluster)))
      dend_ordered <- reorder(dend, wts = colMeans(mat_cluster), agglo.FUN = mean)
      cluster_rows <- dend_ordered
    } else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }

  cell_group <- group.by[1]
  ht_args <- list(
    name = cell_group,
    matrix = mat_list[[cell_group]],
    col = colors,
    row_split = row_split,
    column_split = column_split_list[[cell_group]],
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns_list[[cell_group]],
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_column_slices,
    use_raster = TRUE
  )
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  ht_list <- do.call(Heatmap, args = ht_args)
  features_ordered <- rownames(mat_list[[1]])[unlist(suppressWarnings(row_order(ht_list)))]
  feature_metadata[["index"]] <- setNames(object = seq_along(features_ordered), nm = features_ordered)[rownames(feature_metadata)]

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered) / nlabel) / 2)
        index_to <- length(features_ordered)
        index <- unique(round(seq(from = index_from, to = index_to, length.out = nlabel)))
      } else {
        index <- seq_along(features_ordered)
      }
    } else {
      index <- NULL
    }
  } else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the features"), immediate. = TRUE)
    }
  }
  if (length(index) > 0) {
    ha_mark <- HeatmapAnnotation(
      gene = anno_mark(
        at = which(rownames(feature_metadata) %in% features_ordered[index]),
        labels = feature_metadata[which(rownames(feature_metadata) %in% features_ordered[index]), "features"],
        side = ifelse(flip, "top", "left"),
        labels_gp = gpar(fontsize = label_size, col = label_color),
        link_gp = gpar(fontsize = label_size, col = label_color),
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"), show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    } else {
      ha_left <- c(ha_mark, ha_left)
    }
  }

  ha_right <- NULL
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, FALSE))
        } else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = as.character(featan_values),
          col = palette_scp(featan_values, palette = palette, palcolor = palcolor),
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, labels = levels(featan_values),
          legend_gp = gpar(fill = palette_scp(featan_values, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        col_fun <- colorRamp2(
          breaks = seq(min(featan_values, na.rm = TRUE), max(featan_values, na.rm = TRUE), length = 100),
          colors = palette_scp(palette = palette, palcolor = palcolor)
        )
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = featan_values,
          col = col_fun,
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, col_fun = col_fun, border = TRUE
        )
      }
    }
  }

  enrichment <- heatmap_enrichment(
    geneID = feature_metadata[["features"]], geneID_groups = feature_metadata[["feature_split"]],
    feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor,
    ha_right = ha_right, flip = flip,
    anno_terms = anno_terms, anno_keys = anno_keys, anno_features = anno_features,
    terms_width = terms_width, terms_fontsize = terms_fontsize,
    keys_width = keys_width, keys_fontsize = keys_fontsize,
    features_width = features_width, features_fontsize = features_fontsize,
    IDtype = IDtype, species = species, db_update = db_update, db_version = db_version, db_combine = db_combine, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror,
    db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize,
    GO_simplify = GO_simplify, GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff,
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, min_word_length = min_word_length,
    exclude_words = exclude_words
  )
  res <- enrichment$res
  ha_right <- enrichment$ha_right

  ht_list <- NULL
  for (cell_group in group.by) {
    if (cell_group == group.by[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (cell_group == group.by[length(group.by)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }

    ht_args <- list(
      name = cell_group,
      matrix = if (flip) t(mat_list[[cell_group]]) else mat_list[[cell_group]],
      col = colors,
      row_title = if (flip) ifelse(cell_group != "No.group.by", cell_group, "") else character(0),
      row_title_side = row_title_side,
      column_title = if (flip) character(0) else ifelse(cell_group != "No.group.by", cell_group, ""),
      column_title_side = column_title_side,
      row_title_rot = row_title_rot,
      column_title_rot = column_title_rot,
      row_split = if (flip) column_split_list[[cell_group]] else row_split,
      column_split = if (flip) row_split else column_split_list[[cell_group]],
      cluster_rows = if (flip) cluster_columns_list[[cell_group]] else cluster_rows,
      cluster_columns = if (flip) cluster_rows else cluster_columns_list[[cell_group]],
      cluster_row_slices = if (flip) cluster_column_slices else cluster_row_slices,
      cluster_column_slices = if (flip) cluster_row_slices else cluster_column_slices,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      top_annotation = if (flip) left_annotation else ha_top_list[[cell_group]],
      left_annotation = if (flip) ha_top_list[[cell_group]] else left_annotation,
      bottom_annotation = if (flip) right_annotation else NULL,
      right_annotation = if (flip) NULL else right_annotation,
      show_heatmap_legend = FALSE,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      width = if (is.numeric(width[cell_group])) unit(width[cell_group], units = units) else NULL,
      height = if (is.numeric(height[cell_group])) unit(height[cell_group], units = units) else NULL
    )
    if (!is.null(split.by) && !isTRUE(cluster_column_slices)) {
      groups_order <- sapply(strsplit(levels(column_split_list[[cell_group]]), " : "), function(x) x[[1]])
      gaps_order <- paste(groups_order[2:length(groups_order)], groups_order[1:(length(groups_order) - 1)], sep = "->")
      gaps <- rep(unit(1, "mm"), length(gaps_order))
      gaps[groups_order[2:length(groups_order)] == groups_order[1:(length(groups_order) - 1)]] <- unit(0, "mm")
      if (isTRUE(flip)) {
        ht_args[["row_gap"]] <- gaps
      } else {
        ht_args[["column_gap"]] <- gaps
      }
    }
    if (any(names(ht_params) %in% names(ht_args))) {
      warning("ht_params: ", paste0(intersect(names(ht_params), names(ht_args)), collapse = ","), " were duplicated and will not be used.", immediate. = TRUE)
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(Heatmap, args = ht_args)
      } else {
        ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
      }
    } else {
      ht_list <- ht_list + do.call(Heatmap, args = ht_args)
    }
  }

  if ((!is.null(row_split) && length(index) > 0) || any(c(anno_terms, anno_keys, anno_features)) || !is.null(width) || !is.null(height)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width, height = height, units = units,
    ha_top_list = ha_top_list, ha_left = ha_left, ha_right = ha_right,
    ht_list = ht_list, legend_list = lgd, flip = flip
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  # cat("width:", width, "\n")
  # cat("height:", height, "\n")
  # cat("width_sum:", width_sum, "\n")
  # cat("height_sum:", height_sum, "\n")

  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(
      width = width, width_sum = width_sum, height = height, height_sum = height_sum, units = units,
      ht_list = ht_list, legend_list = lgd
    )
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    # cat("ht_width:", ht_width, "\n")
    # cat("ht_height:", ht_height, "\n")

    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  } else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  }

  if (isTRUE(fix)) {
    p <- panel_fix_single(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
  } else {
    p <- wrap_plots(gTree)
  }

  return(list(
    plot = p,
    matrix_list = mat_list,
    feature_split = feature_split,
    cell_metadata = cell_metadata,
    feature_metadata = feature_metadata,
    enrichment = res
  ))
}

FeatureCorHeatmap <- function(srt, features, cells) {

}

#' CellCorHeatmap
#'
#' @param srt_query
#' @param srt_ref
#' @param bulk_ref
#' @param query_group
#' @param ref_group
#' @param query_assay
#' @param ref_assay
#' @param query_reduction
#' @param ref_reduction
#' @param query_dims
#' @param ref_dims
#' @param query_collapsing
#' @param ref_collapsing
#' @param features
#' @param features_type
#' @param feature_source
#' @param nfeatures
#' @param DEtest_param
#' @param DE_threshold
#' @param distance_metric
#' @param k
#' @param filter_lowfreq
#' @param prefix
#' @param cluster_columns
#' @param cluster_rows
#' @param nlabel
#' @param label_cutoff
#' @param label_by
#' @param border
#' @param flip
#' @param limits
#' @param show_row_names
#' @param show_column_names
#' @param row_names_side
#' @param column_names_side
#' @param row_names_rot
#' @param column_names_rot
#' @param row_title_side
#' @param column_title_side
#' @param row_title_rot
#' @param column_title_rot
#' @param heatmap_palette
#' @param heatmap_palcolor
#' @param query_group_palette
#' @param query_group_palcolor
#' @param ref_group_palette
#' @param ref_group_palcolor
#' @param query_cell_annotation
#' @param query_cell_annotation_palette
#' @param query_cell_annotation_palcolor
#' @param query_cell_annotation_params
#' @param ref_cell_annotation
#' @param ref_cell_annotation_palette
#' @param ref_cell_annotation_palcolor
#' @param ref_cell_annotation_params
#' @param use_raster
#' @param raster_device
#' @param height
#' @param width
#' @param units
#' @param seed
#' @param ht_params
#' @param label_size
#' @param raster_by_magick
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' ht1 <- CellCorHeatmap(srt_query = pancreas_sub, query_group = "SubCellType")
#' ht1$plot
#'
#' data("panc8_sub")
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(capitalize(rownames(panc8_sub), force_tolower = TRUE))
#' panc8_sub <- RenameFeatures(panc8_sub, newnames = genenames)
#' panc8_sub <- check_srtMerge(panc8_sub, batch = "tech")[["srtMerge"]]
#'
#' ht2 <- CellCorHeatmap(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub, nlabel = 3, label_cutoff = 0.6,
#'   query_group = "SubCellType", ref_group = "celltype",
#'   query_cell_annotation = "Phase", query_cell_annotation_palette = "Set2",
#'   ref_cell_annotation = "tech", ref_cell_annotation_palette = "Set3",
#'   width = 3, height = 2
#' )
#' ht2$plot
#'
#' ht3 <- CellCorHeatmap(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   query_group = "SubCellType", query_collapsing = FALSE, cluster_rows = TRUE,
#'   ref_group = "celltype", ref_collapsing = FALSE, cluster_columns = TRUE
#' )
#' ht3$plot
#'
#' ht4 <- CellCorHeatmap(
#'   srt_query = pancreas_sub, srt_ref = panc8_sub,
#'   show_row_names = TRUE, show_column_names = TRUE,
#'   query_group = "SubCellType", ref_group = "celltype",
#'   query_cell_annotation = c("Sox9", "Rbp4", "Gcg"),
#'   ref_cell_annotation = c("Sox9", "Rbp4", "Gcg")
#' )
#' ht4$plot
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Legend HeatmapAnnotation anno_block anno_simple anno_customize Heatmap draw pindex restore_matrix %v%
#' @importFrom grid gpar grid.grabExpr grid.lines grid.rect grid.points grid.draw
#' @importFrom ggplot2 theme_void theme facet_null
#' @importFrom patchwork wrap_plots
#' @importFrom methods getFunction
#' @importFrom dplyr %>% filter group_by arrange desc across mutate distinct n .data "%>%"
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @export
CellCorHeatmap <- function(srt_query, srt_ref = NULL, bulk_ref = NULL,
                           query_group = NULL, ref_group = NULL,
                           query_assay = NULL, ref_assay = NULL,
                           query_reduction = NULL, ref_reduction = NULL,
                           query_dims = 1:30, ref_dims = 1:30,
                           query_collapsing = !is.null(query_group), ref_collapsing = TRUE,
                           features = NULL, features_type = c("HVF", "DE"), feature_source = "both", nfeatures = 2000,
                           DEtest_param = list(max.cells.per.ident = 200, test.use = "wilcox"),
                           DE_threshold = "p_val_adj < 0.05",
                           distance_metric = "cosine", k = 30,
                           filter_lowfreq = 0, prefix = "KNNPredict",
                           border = TRUE, flip = FALSE, limits = NULL,
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = "left", column_names_side = "top", row_names_rot = 0, column_names_rot = 90,
                           row_title_side = "left", column_title_side = "top", row_title_rot = 90, column_title_rot = 0,
                           nlabel = 0, label_cutoff = 0, label_by = "row", label_size = 10,
                           heatmap_palette = "RdBu", heatmap_palcolor = NULL,
                           query_group_palette = "Paired", query_group_palcolor = NULL,
                           ref_group_palette = "jama", ref_group_palcolor = NULL,
                           query_cell_annotation = NULL, query_cell_annotation_palette = "Paired", query_cell_annotation_palcolor = NULL, query_cell_annotation_params = if (flip) list(height = grid::unit(1, "cm")) else list(width = grid::unit(1, "cm")),
                           ref_cell_annotation = NULL, ref_cell_annotation_palette = "Paired", ref_cell_annotation_palcolor = NULL, ref_cell_annotation_params = if (flip) list(width = grid::unit(1, "cm")) else list(height = grid::unit(1, "cm")),
                           use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, height = NULL, width = NULL, units = "inch",
                           seed = 11, ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }

  simil_method <- c(
    "cosine", "pearson", "spearman", "correlation", "jaccard", "ejaccard", "dice", "edice",
    "hamman", "simple matching", "faith"
  )
  dist_method <- c(
    "euclidean", "chisquared", "kullback", "manhattan", "maximum", "canberra",
    "minkowski", "hamming"
  )
  if (is.null(srt_ref) && is.null(bulk_ref)) {
    srt_ref <- srt_query
    ref_group <- query_group
    ref_assay <- query_assay
    ref_reduction <- query_reduction
    ref_dims <- query_dims
    ref_collapsing <- query_collapsing
    ref_group_palette <- query_group_palette
    ref_group_palcolor <- query_group_palcolor
    ref_cell_annotation <- query_cell_annotation
    ref_cell_annotation_palette <- query_cell_annotation_palette
    ref_cell_annotation_palcolor <- query_cell_annotation_palcolor
    ref_cell_annotation_params <- query_cell_annotation_params
  }
  if (!is.null(bulk_ref)) {
    srt_ref <- CreateSeuratObject(counts = bulk_ref, meta.data = data.frame(celltype = colnames(bulk_ref)), assay = "RNA")
    ref_group <- "CellType"
    ref_assay <- "RNA"
    ref_reduction <- NULL
    ref_collapsing <- FALSE
  }
  query_assay <- query_assay %||% DefaultAssay(srt_query)
  ref_assay <- ref_assay %||% DefaultAssay(srt_ref)
  other_params <- list(
    query_group = query_group, query_reduction = query_reduction, query_assay = query_assay, query_dims = query_dims, query_collapsing = query_collapsing,
    ref_group = ref_group, ref_reduction = ref_reduction, ref_assay = ref_assay, ref_dims = ref_dims, ref_collapsing = ref_collapsing
  )
  if (is.null(srt_query@tools[[paste0(prefix, "_classification")]][["distance_matrix"]]) ||
    !identical(other_params, srt_query@tools[[paste0(prefix, "_classification")]][["other_params"]])) {
    srt_query <- RunKNNPredict(
      srt_query = srt_query, srt_ref = srt_ref,
      query_group = query_group, ref_group = ref_group,
      query_assay = query_assay, ref_assay = ref_assay,
      query_reduction = query_reduction, ref_reduction = ref_reduction, query_dims = query_dims, ref_dims = ref_dims,
      query_collapsing = query_collapsing, ref_collapsing = ref_collapsing,
      features = features, features_type = features_type, feature_source = feature_source, nfeatures = nfeatures,
      DEtest_param = DEtest_param,
      DE_threshold = DE_threshold,
      distance_metric = distance_metric, k = k,
      filter_lowfreq = filter_lowfreq, prefix = prefix,
      nn_method = "raw", return_full_distance_matrix = TRUE
    )
  }
  distance_matrix <- srt_query@tools[[paste0(prefix, "_classification")]][["distance_matrix"]]
  distance_metric <- srt_query@tools[[paste0(prefix, "_classification")]][["distance_metric"]]
  if (distance_metric %in% simil_method) {
    simil_matrix <- t(as.matrix(1 - distance_matrix))
    simil_name <- paste0(capitalize(distance_metric), " similarity")
  } else if (distance_metric %in% dist_method) {
    simil_matrix <- t(as.matrix(1 - distance_matrix / max(distance_matrix, na.rm = TRUE)))
    simil_name <- paste0("1-dist[", distance_metric, "]/max(dist[", distance_metric, "])")
  }
  simil_matrix[is.infinite(simil_matrix)] <- max(abs(simil_matrix[!is.infinite(simil_matrix)]), na.rm = TRUE) * ifelse(simil_matrix[is.infinite(simil_matrix)] > 0, 1, -1)
  simil_matrix[is.na(simil_matrix)] <- 0

  cell_groups <- list()
  if (is.null(query_group)) {
    srt_query@meta.data[["No.group.by"]] <- factor("")
    query_group <- "No.group.by"
  }
  if (is.null(ref_group)) {
    srt_ref@meta.data[["No.group.by"]] <- factor("")
    ref_group <- "No.group.by"
  }
  if (!is.factor(srt_query[[query_group, drop = TRUE]])) {
    srt_query@meta.data[[query_group]] <- factor(srt_query[[query_group, drop = TRUE]], levels = unique(srt_query[[query_group, drop = TRUE]]))
  }
  cell_groups[["query_group"]] <- unlist(lapply(levels(srt_query[[query_group, drop = TRUE]]), function(x) {
    cells_sub <- colnames(srt_query)[which(srt_query[[query_group, drop = TRUE]] == x)]
    out <- setNames(object = rep(x, length(cells_sub)), nm = cells_sub)
    return(out)
  }), use.names = TRUE)
  levels <- levels(srt_query[[query_group, drop = TRUE]])
  cell_groups[["query_group"]] <- factor(cell_groups[["query_group"]], levels = levels[levels %in% cell_groups[["query_group"]]])


  if (!is.factor(srt_ref[[ref_group, drop = TRUE]])) {
    srt_ref@meta.data[[ref_group]] <- factor(srt_ref[[ref_group, drop = TRUE]], levels = unique(srt_ref[[ref_group, drop = TRUE]]))
  }
  cell_groups[["ref_group"]] <- unlist(lapply(levels(srt_ref[[ref_group, drop = TRUE]]), function(x) {
    cells_sub <- colnames(srt_ref)[which(srt_ref[[ref_group, drop = TRUE]] == x)]
    out <- setNames(object = rep(x, length(cells_sub)), nm = cells_sub)
    return(out)
  }), use.names = TRUE)
  levels <- levels(srt_ref[[ref_group, drop = TRUE]])
  cell_groups[["ref_group"]] <- factor(cell_groups[["ref_group"]], levels = levels[levels %in% cell_groups[["ref_group"]]])

  if (isTRUE(query_collapsing)) {
    simil_matrix <- simil_matrix[levels(cell_groups[["query_group"]]), , drop = FALSE]
  } else {
    simil_matrix <- simil_matrix[names(cell_groups[["query_group"]]), , drop = FALSE]
  }
  if (isTRUE(ref_collapsing)) {
    simil_matrix <- simil_matrix[, levels(cell_groups[["ref_group"]]), drop = FALSE]
  } else {
    simil_matrix <- simil_matrix[, names(cell_groups[["ref_group"]]), drop = FALSE]
  }

  if (!is.null(query_cell_annotation)) {
    if (length(query_cell_annotation_palette) == 1) {
      query_cell_annotation_palette <- rep(query_cell_annotation_palette, length(query_cell_annotation))
    }
    if (length(query_cell_annotation_palcolor) == 1) {
      query_cell_annotation_palcolor <- rep(query_cell_annotation_palcolor, length(query_cell_annotation))
    }
    npal <- unique(c(length(query_cell_annotation_palette), length(query_cell_annotation_palcolor), length(query_cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("query_cell_annotation_palette and query_cell_annotation_palcolor must be the same length as query_cell_annotation")
    }
    if (any(!query_cell_annotation %in% c(colnames(srt_query@meta.data), rownames(srt_query[[query_assay]])))) {
      stop("query_cell_annotation: ", paste0(query_cell_annotation[!query_cell_annotation %in% c(colnames(srt_query@meta.data), rownames(srt_query[[query_assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(ref_cell_annotation)) {
    if (length(ref_cell_annotation_palette) == 1) {
      ref_cell_annotation_palette <- rep(ref_cell_annotation_palette, length(ref_cell_annotation))
    }
    if (length(ref_cell_annotation_palcolor) == 1) {
      ref_cell_annotation_palcolor <- rep(ref_cell_annotation_palcolor, length(ref_cell_annotation))
    }
    npal <- unique(c(length(ref_cell_annotation_palette), length(ref_cell_annotation_palcolor), length(ref_cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("ref_cell_annotation_palette and ref_cell_annotation_palcolor must be the same length as ref_cell_annotation")
    }
    if (any(!ref_cell_annotation %in% c(colnames(srt_ref@meta.data), rownames(srt_ref[[ref_assay]])))) {
      stop("ref_cell_annotation: ", paste0(ref_cell_annotation[!ref_cell_annotation %in% c(colnames(srt_ref@meta.data), rownames(srt_ref[[ref_assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }

  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
  }
  if (is.null(limits)) {
    colors <- colorRamp2(seq(min(simil_matrix, na.rm = TRUE), max(simil_matrix, na.rm = TRUE), length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  } else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }

  cell_metadata <- data.frame(
    row.names = c(paste0("query_", colnames(srt_query)), paste0("ref_", colnames(srt_ref))),
    cells = c(colnames(srt_query), colnames(srt_ref))
  )
  query_metadata <- cbind.data.frame(
    srt_query@meta.data[cell_metadata[["cells"]], c(query_group, intersect(query_cell_annotation, colnames(srt_query@meta.data))), drop = FALSE],
    as.data.frame(t(srt_query[[query_assay]]@data[intersect(query_cell_annotation, rownames(srt_query[[query_assay]])) %||% integer(), , drop = FALSE]))[cell_metadata[["cells"]], , drop = FALSE]
  )
  colnames(query_metadata) <- paste0("query_", colnames(query_metadata))
  ref_metadata <- cbind.data.frame(
    srt_ref@meta.data[cell_metadata[["cells"]], c(ref_group, intersect(ref_cell_annotation, colnames(srt_ref@meta.data))), drop = FALSE],
    as.data.frame(t(srt_ref[[ref_assay]]@data[intersect(ref_cell_annotation, rownames(srt_ref[[ref_assay]])) %||% integer(), , drop = FALSE]))[cell_metadata[["cells"]], , drop = FALSE]
  )
  colnames(ref_metadata) <- paste0("ref_", colnames(ref_metadata))
  cell_metadata <- cbind.data.frame(cell_metadata, cbind.data.frame(query_metadata, ref_metadata))

  lgd <- list()
  lgd[["ht"]] <- Legend(title = simil_name, col_fun = colors, border = TRUE)

  ha_query_list <- NULL
  if (query_group != "No.group.by") {
    if (isFALSE(query_collapsing) && ((isFALSE(flip) & isTRUE(cluster_rows)) || (isTRUE(flip) & isTRUE(cluster_columns)))) {
      query_cell_annotation <- c(query_group, query_cell_annotation)
    } else {
      funbody <- paste0(
        "
        grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt_query[[query_group, drop = TRUE]]), collapse = "','"), "')"), ",palette = '", query_group_palette, "',palcolor=c(", paste0("'", paste0(query_group_palcolor, collapse = "','"), "'"), "))[nm]))
      "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      if (isTRUE(query_collapsing)) {
        anno[[paste0(c("Query", query_group), collapse = ":")]] <- anno_block(
          align_to = split(seq_along(levels(cell_groups[["query_group"]])), levels(cell_groups[["query_group"]])),
          panel_fun = getFunction("panel_fun", where = environment()),
          which = ifelse(flip, "column", "row"),
          show_name = FALSE
        )
      } else {
        anno[[paste0(c("Query", query_group), collapse = ":")]] <- anno_block(
          align_to = split(seq_along(cell_groups[["query_group"]]), cell_groups[["query_group"]]),
          panel_fun = getFunction("panel_fun", where = environment()),
          which = ifelse(flip, "column", "row"),
          show_name = FALSE
        )
      }
      ha_cell_group <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "bottom"), border = TRUE))
      ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- ha_cell_group
      lgd[[paste0(c("Query", query_group), collapse = ":")]] <- Legend(
        title = paste0(c("Query", query_group), collapse = ":"), labels = levels(srt_query[[query_group, drop = TRUE]]),
        legend_gp = gpar(fill = palette_scp(levels(srt_query[[query_group, drop = TRUE]]), palette = query_group_palette, palcolor = query_group_palcolor)), border = TRUE
      )
    }
  }

  ha_ref_list <- NULL
  if (ref_group != "No.group.by") {
    if (isFALSE(ref_collapsing) && ((isFALSE(flip) & isTRUE(cluster_columns)) || (isTRUE(flip) & isTRUE(cluster_rows)))) {
      ref_cell_annotation <- c(ref_group, ref_cell_annotation)
    } else {
      funbody <- paste0(
        "
        grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(srt_ref[[ref_group, drop = TRUE]]), collapse = "','"), "')"), ",palette = '", ref_group_palette, "',palcolor=c(", paste0("'", paste0(ref_group_palcolor, collapse = "','"), "'"), "))[nm]))
      "
      )
      funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
      eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())

      anno <- list()
      if (isTRUE(ref_collapsing)) {
        anno[[paste0(c("Ref", ref_group), collapse = ":")]] <- anno_block(
          align_to = split(seq_along(levels(cell_groups[["ref_group"]])), levels(cell_groups[["ref_group"]])),
          panel_fun = getFunction("panel_fun", where = environment()),
          which = ifelse(!flip, "column", "row"),
          show_name = FALSE
        )
      } else {
        anno[[paste0(c("Ref", ref_group), collapse = ":")]] <- anno_block(
          align_to = split(seq_along(cell_groups[["ref_group"]]), cell_groups[["ref_group"]]),
          panel_fun = getFunction("panel_fun", where = environment()),
          which = ifelse(!flip, "column", "row"),
          show_name = FALSE
        )
      }
      ha_cell_group <- do.call("HeatmapAnnotation", args = c(anno, which = ifelse(!flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(!flip, "left", "bottom"), border = TRUE))
      ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_cell_group
      lgd[[paste0(c("Ref", ref_group), collapse = ":")]] <- Legend(
        title = paste0(c("Ref", ref_group), collapse = ":"), labels = levels(srt_ref[[ref_group, drop = TRUE]]),
        legend_gp = gpar(fill = palette_scp(levels(srt_ref[[ref_group, drop = TRUE]]), palette = ref_group_palette, palcolor = ref_group_palcolor)), border = TRUE
      )
    }
  }

  if (!is.null(query_cell_annotation)) {
    query_subplots_list <- list()
    for (i in seq_along(query_cell_annotation)) {
      cellan <- query_cell_annotation[i]
      palette <- query_cell_annotation_palette[i]
      palcolor <- query_cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, paste0("query_", cellan)]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        if (isTRUE(query_collapsing)) {
          subplots <- CellStatPlot(srt_query,
            flip = !flip,
            cells = gsub("query_", "", names(cell_groups[["query_group"]])), plot_type = "pie",
            group.by = query_group, stat.by = cellan,
            palette = palette, palcolor = palcolor,
            individual = TRUE, combine = FALSE
          )
          query_subplots_list[[paste0(cellan, ":", query_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(query_subplots_list[['", cellan, ":", query_group, "']]", "[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", query_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[["query_group"]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(query_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(query_group, x[1], ""), collapse = ":")
            }
          })

          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(anno_args, query_cell_annotation_params[setdiff(names(query_cell_annotation_params), names(anno_args))])
          ha_query <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]])) {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- ha_query
          } else {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- c(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]], ha_query)
          }
        } else {
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = as.character(cell_anno[paste0("query_", names(cell_groups[["query_group"]]))]),
            col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
            which = ifelse(flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(anno_args, query_cell_annotation_params[setdiff(names(query_cell_annotation_params), names(anno_args))])
          ha_query <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]])) {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- ha_query
          } else {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- c(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]], ha_query)
          }
        }
        lgd[[paste0(c("Query", cellan), collapse = ":")]] <- Legend(
          title = paste0(c("Query", cellan), collapse = ":"), labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        if (isTRUE(query_collapsing)) {
          subplots <- FeatureStatPlot(srt_query,
            assay = query_assay, slot = "data", flip = !flip,
            stat.by = cellan, cells = gsub("query_", "", names(cell_groups[["query_group"]])),
            group.by = query_group,
            palette = query_group_palette,
            palcolor = query_group_palcolor,
            fill.by = "group", same.y.lims = TRUE,
            individual = TRUE, combine = FALSE
          )
          query_subplots_list[[paste0(cellan, ":", query_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(query_subplots_list[['", cellan, ":", query_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", query_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[["query_group"]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(cellan, query_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(cellan, query_group, x[1], ""), collapse = ":")
            }
          })
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(anno_args, query_cell_annotation_params[setdiff(names(query_cell_annotation_params), names(anno_args))])
          ha_query <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]])) {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- ha_query
          } else {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- c(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]], ha_query)
          }
        } else {
          col_fun <- colorRamp2(
            breaks = seq(min(cell_anno, na.rm = TRUE), max(cell_anno, na.rm = TRUE), length = 100),
            colors = palette_scp(palette = palette, palcolor = palcolor)
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = cell_anno[paste0("query_", names(cell_groups[["query_group"]]))],
            col = col_fun,
            which = ifelse(flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(flip, "left", "bottom")
          )
          anno_args <- c(anno_args, query_cell_annotation_params[setdiff(names(query_cell_annotation_params), names(anno_args))])
          ha_query <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]])) {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- ha_query
          } else {
            ha_query_list[[paste0(c("Query", query_group), collapse = ":")]] <- c(ha_query_list[[paste0(c("Query", query_group), collapse = ":")]], ha_query)
          }
          lgd[[paste0(c("Query", cellan), collapse = ":")]] <- Legend(
            title = paste0(c("Query", cellan), collapse = ":"), col_fun = col_fun, border = TRUE
          )
        }
      }
    }
  }

  if (!is.null(ref_cell_annotation)) {
    ref_subplots_list <- list()
    for (i in seq_along(ref_cell_annotation)) {
      cellan <- ref_cell_annotation[i]
      palette <- ref_cell_annotation_palette[i]
      palcolor <- ref_cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, paste0("ref_", cellan)]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        if (isTRUE(ref_collapsing)) {
          subplots <- CellStatPlot(srt_ref,
            flip = flip,
            cells = gsub("ref_", "", names(cell_groups[["ref_group"]])), plot_type = "pie",
            group.by = ref_group, stat.by = cellan,
            palette = palette, palcolor = palcolor,
            individual = TRUE, combine = FALSE
          )
          ref_subplots_list[[paste0(cellan, ":", ref_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(ref_subplots_list[['", cellan, ":", ref_group, "']]", "[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", ref_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[["ref_group"]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(ref_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(ref_group, x[1], ""), collapse = ":")
            }
          })

          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(!flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(anno_args, ref_cell_annotation_params[setdiff(names(ref_cell_annotation_params), names(anno_args))])
          ha_ref <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]], ha_ref)
          }
        } else {
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = as.character(cell_anno[paste0("ref_", names(cell_groups[["ref_group"]]))]),
            col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
            which = ifelse(!flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(anno_args, ref_cell_annotation_params[setdiff(names(ref_cell_annotation_params), names(anno_args))])
          ha_ref <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]], ha_ref)
          }
        }
        lgd[[paste0(c("Ref", cellan), collapse = ":")]] <- Legend(
          title = paste0(c("Ref", cellan), collapse = ":"), labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        if (isTRUE(ref_collapsing)) {
          subplots <- FeatureStatPlot(srt_ref,
            assay = ref_assay, slot = "data", flip = flip,
            stat.by = cellan, cells = gsub("ref_", "", names(cell_groups[["ref_group"]])),
            group.by = ref_group,
            palette = ref_group_palette,
            palcolor = ref_group_palcolor,
            fill.by = "group", same.y.lims = TRUE,
            individual = TRUE, combine = FALSE
          )
          ref_subplots_list[[paste0(cellan, ":", ref_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(ref_subplots_list[['", cellan, ":", ref_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(legend.position = 'none'));
              g$name <- '", paste0(cellan, ":", ref_group, "-", nm), "';
              grid.draw(g)
              "
            )
            funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
            eval(parse(text = paste("graphics[[nm]] <- function(x, y, w, h) {", funbody, "}", sep = "")), envir = environment())
          }
          x_nm <- sapply(strsplit(levels(cell_groups[["ref_group"]]), " : "), function(x) {
            if (length(x) == 2) {
              paste0(c(cellan, ref_group, x[1], x[2]), collapse = ":")
            } else {
              paste0(c(cellan, ref_group, x[1], ""), collapse = ":")
            }
          })
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_customize(
            x = x_nm,
            graphics = graphics,
            which = ifelse(!flip, "column", "row"),
            border = TRUE,
            verbose = FALSE
          )
          anno_args <- c(ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(anno_args, ref_cell_annotation_params[setdiff(names(ref_cell_annotation_params), names(anno_args))])
          ha_ref <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]], ha_ref)
          }
        } else {
          col_fun <- colorRamp2(
            breaks = seq(min(cell_anno, na.rm = TRUE), max(cell_anno, na.rm = TRUE), length = 100),
            colors = palette_scp(palette = palette, palcolor = palcolor)
          )
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = cell_anno[paste0("ref_", names(cell_groups[["ref_group"]]))],
            col = col_fun,
            which = ifelse(!flip, "column", "row"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell,
            which = ifelse(!flip, "column", "row"),
            show_annotation_name = TRUE,
            annotation_name_side = ifelse(!flip, "left", "bottom")
          )
          anno_args <- c(anno_args, ref_cell_annotation_params[setdiff(names(ref_cell_annotation_params), names(anno_args))])
          ha_ref <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]])) {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- ha_ref
          } else {
            ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]] <- c(ha_ref_list[[paste0(c("Ref", ref_group), collapse = ":")]], ha_ref)
          }
          lgd[[paste0(c("Ref", cellan), collapse = ":")]] <- Legend(
            title = paste0(c("Ref", cellan), collapse = ":"), col_fun = col_fun, border = TRUE
          )
        }
      }
    }
  }

  layer_fun <- function(j, i, x, y, w, h, fill) {
    if (nlabel > 0) {
      if (flip) {
        mat <- t(simil_matrix)
      } else {
        mat <- simil_matrix
      }
      value <- pindex(mat, i, j)
      ind_mat <- restore_matrix(j, i, x, y)

      inds <- NULL
      if (label_by %in% c("row", "both")) {
        for (row in 1:nrow(ind_mat)) {
          ind <- ind_mat[row, ]
          ind <- ind[which(value[ind] >= max(c(sort(value[ind], decreasing = TRUE)[nlabel]), na.rm = TRUE) & value[ind] >= label_cutoff)]
          inds <- c(inds, ind)
        }
      }
      if (label_by %in% c("column", "both")) {
        for (column in 1:ncol(ind_mat)) {
          ind <- ind_mat[, column]
          ind <- ind[which(value[ind] >= max(c(sort(value[ind], decreasing = TRUE)[nlabel]), na.rm = TRUE) & value[ind] >= label_cutoff)]
          inds <- c(inds, ind)
        }
      }
      if (label_by == "both") {
        inds <- inds[duplicated(inds)]
      }
      if (length(inds) > 0) {
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        lapply(theta, function(i) {
          x_out <- x[inds] + unit(cos(i) * label_size / 30, "mm")
          y_out <- y[inds] + unit(sin(i) * label_size / 30, "mm")
          grid.text(round(value[inds], 2), x = x_out, y = y_out, gp = gpar(fontsize = label_size, col = "white"))
        })
        grid.text(round(value[inds], 2), x[inds], y[inds], gp = gpar(fontsize = label_size, col = "black"))
      }
    }
  }

  ht_list <- NULL
  ht_args <- list(
    name = simil_name,
    matrix = if (flip) t(simil_matrix) else simil_matrix,
    col = colors,
    layer_fun = layer_fun,
    row_title = if (flip) paste0(c("Ref", ref_group), collapse = ":") else paste0(c("Query", query_group), collapse = ":"),
    row_title_side = row_title_side,
    column_title = if (flip) paste0(c("Query", query_group), collapse = ":") else paste0(c("Ref", ref_group), collapse = ":"),
    column_title_side = column_title_side,
    row_title_rot = row_title_rot,
    column_title_rot = column_title_rot,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_side = row_names_side,
    column_names_side = column_names_side,
    row_names_rot = row_names_rot,
    column_names_rot = column_names_rot,
    top_annotation = if (flip) ha_query_list[[1]] else ha_ref_list[[1]],
    left_annotation = if (flip) ha_ref_list[[1]] else ha_query_list[[1]],
    show_heatmap_legend = FALSE,
    border = border,
    use_raster = use_raster,
    raster_device = raster_device,
    raster_by_magick = raster_by_magick,
    width = if (is.numeric(width)) unit(width, units = units) else NULL,
    height = if (is.numeric(height)) unit(height, units = units) else NULL
  )
  if (any(names(ht_params) %in% names(ht_args))) {
    warning("ht_params: ", paste0(intersect(names(ht_params), names(ht_args)), collapse = ","), " were duplicated and will not be used.", immediate. = TRUE)
  }
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  if (isTRUE(flip)) {
    if (is.null(ht_list)) {
      ht_list <- do.call(Heatmap, args = ht_args)
    } else {
      ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
    }
  } else {
    ht_list <- ht_list + do.call(Heatmap, args = ht_args)
  }

  if (!is.null(width) || !is.null(height)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width, height = height, units = units,
    ha_top_list = ha_ref_list, ha_left = ha_query_list[[1]], ha_right = NULL,
    ht_list = ht_list, legend_list = lgd, flip = flip
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  # cat("width:", width, "\n")
  # cat("height:", height, "\n")
  # cat("width_sum:", width_sum, "\n")
  # cat("height_sum:", height_sum, "\n")

  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(
      width = width, width_sum = width_sum, height = height, height_sum = height_sum, units = units,
      ht_list = ht_list, legend_list = lgd
    )
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    # cat("ht_width:", ht_width, "\n")
    # cat("ht_height:", ht_height, "\n")

    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  } else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  }

  if (isTRUE(fix)) {
    p <- panel_fix_single(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
  } else {
    p <- wrap_plots(gTree)
  }

  return(list(
    plot = p,
    simil_matrix = simil_matrix,
    simil_name = simil_name,
    cell_metadata = cell_metadata
  ))
}

#' Heatmap plot for dynamic features along lineages
#'
#' @param srt
#' @param lineages
#' @param feature_from
#' @param exp_method
#' @param slot
#' @param assay
#' @param use_fitted
#' @param lib_normalize
#' @param libsize
#' @param min_expcells
#' @param r.sq
#' @param dev.expl
#' @param padjust
#' @param cell_density
#' @param order_by
#' @param decreasing
#' @param feature_split
#' @param n_split
#' @param feature_split_by
#' @param split_method
#' @param fuzzification
#' @param show_fuzzification
#' @param nlabel
#' @param features_label
#' @param label_size
#' @param label_color
#' @param pseudotime_label
#' @param pseudotime_label_color
#' @param pseudotime_label_linetype
#' @param pseudotime_label_linewidth
#' @param heatmap_palette
#' @param pseudotime_palette
#' @param cell_annotation
#' @param cell_annotation_palette
#' @param cell_annotation_palcolor
#' @param feature_annotation
#' @param feature_annotation_palette
#' @param feature_annotation_palcolor
#' @param reverse_ht
#' @param use_raster
#' @param height
#' @param width
#' @param units
#' @param seed
#' @param features
#' @param border
#' @param flip
#' @param family
#' @param cluster_features_by
#' @param cluster_rows
#' @param cluster_row_slices
#' @param cluster_columns
#' @param cluster_column_slices
#' @param show_row_names
#' @param show_column_names
#' @param row_names_side
#' @param column_names_side
#' @param row_names_rot
#' @param column_names_rot
#' @param row_title_side
#' @param column_title_side
#' @param row_title_rot
#' @param column_title_rot
#' @param anno_terms
#' @param anno_keys
#' @param anno_features
#' @param terms_width
#' @param terms_fontsize
#' @param keys_width
#' @param keys_fontsize
#' @param features_width
#' @param features_fontsize
#' @param IDtype
#' @param species
#' @param db_update
#' @param db_version
#' @param convert_species
#' @param Ensembl_version
#' @param mirror
#' @param db
#' @param TERM2GENE
#' @param TERM2NAME
#' @param minGSSize
#' @param maxGSSize
#' @param GO_simplify
#' @param GO_simplify_cutoff
#' @param simplify_method
#' @param simplify_similarityCutoff
#' @param pvalueCutoff
#' @param padjustCutoff
#' @param topTerm
#' @param show_termid
#' @param topWord
#' @param min_word_length
#' @param exclude_words
#' @param heatmap_palcolor
#' @param pseudotime_palcolor
#' @param feature_split_palette
#' @param feature_split_palcolor
#' @param cell_annotation_params
#' @param feature_annotation_params
#' @param separate_annotation
#' @param separate_annotation_palette
#' @param separate_annotation_palcolor
#' @param separate_annotation_params
#' @param raster_device
#' @param ht_params
#' @param limits
#' @param raster_by_magick
#' @param db_combine
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#' ht1 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   n_split = 5,
#'   split_method = "kmeans-peaktime",
#'   cell_annotation = "SubCellType"
#' )
#' ht1$plot
#' panel_fix(ht1$plot, raster = TRUE, dpi = 50)
#'
#' ht2 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   features = c("Sox9", "Neurod2", "Isl1", "Rbp4", "Pyy", "S_score", "G2M_score"),
#'   cell_annotation = "SubCellType"
#' )
#' ht2$plot
#' panel_fix(ht2$plot, height = 5, width = 5, raster = TRUE, dpi = 50)
#'
#' ht3 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_split = 5,
#'   split_method = "kmeans",
#'   cluster_rows = TRUE,
#'   cell_annotation = "SubCellType"
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "SP"))
#' ht4 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   use_fitted = TRUE,
#'   n_split = 6,
#'   split_method = "mfuzz",
#'   heatmap_palette = "viridis",
#'   cell_annotation = c("SubCellType", "Phase", "G2M_score"),
#'   cell_annotation_palette = c("Paired", "jama", "Purples"),
#'   separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
#'   separate_annotation_palette = c("Paired", "Set1"),
#'   separate_annotation_params = list(height = grid::unit(2, "cm")),
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   pseudotime_label = 25,
#'   pseudotime_label_color = "red"
#' )
#' ht4$plot
#'
#' ht5 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   use_fitted = TRUE,
#'   n_split = 6,
#'   split_method = "mfuzz",
#'   heatmap_palette = "viridis",
#'   cell_annotation = c("SubCellType", "Phase", "G2M_score"),
#'   cell_annotation_palette = c("Paired", "jama", "Purples"),
#'   separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
#'   separate_annotation_palette = c("Paired", "Set1"),
#'   separate_annotation_params = list(width = grid::unit(2, "cm")),
#'   feature_annotation = c("TF", "SP"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   pseudotime_label = 25,
#'   pseudotime_label_color = "red",
#'   flip = TRUE, column_title_rot = 45
#' )
#' ht5$plot
#'
#' ht6 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   cell_annotation = "SubCellType",
#'   n_split = 5, split_method = "mfuzz",
#'   species = "Mus_musculus", db = "GO_BP",
#'   anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht6$plot
#'
#' @importFrom Seurat GetAssayData NormalizeData DefaultAssay
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap Legend HeatmapAnnotation anno_empty anno_mark anno_simple anno_textbox draw decorate_heatmap_body width.HeatmapAnnotation height.HeatmapAnnotation width.Legends height.Legends decorate_annotation row_order %v%
#' @importFrom stats kmeans
#' @importFrom patchwork wrap_plots
#' @importFrom grid gpar grid.lines grid.text
#' @importFrom gtable gtable_add_padding
#' @importFrom Matrix t
#' @importFrom dplyr %>% filter group_by arrange desc across mutate reframe distinct n .data "%>%"
#' @importFrom rlang %||%
#' @export
DynamicHeatmap <- function(srt, lineages, features = NULL, feature_from = lineages, use_fitted = FALSE, border = TRUE, flip = FALSE,
                           min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05, cell_density = 1, order_by = c("peaktime", "valleytime"),
                           slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), limits = NULL,
                           lib_normalize = identical(slot, "counts"), libsize = NULL, family = NULL,
                           cluster_features_by = NULL, cluster_rows = FALSE, cluster_row_slices = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                           row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                           feature_split = NULL, feature_split_by = NULL, n_split = NULL,
                           split_method = c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime"), decreasing = FALSE,
                           fuzzification = NULL, show_fuzzification = FALSE,
                           anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                           terms_width = unit(4, "in"), terms_fontsize = 8,
                           keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                           features_width = unit(2, "in"), features_fontsize = c(6, 10),
                           IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                           db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                           GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, min_word_length = 3,
                           exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "organization", "system", "regulation", "positive", "negative", "response", "process"),
                           nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                           pseudotime_label = NULL, pseudotime_label_color = "black", pseudotime_label_linetype = 2, pseudotime_label_linewidth = 3,
                           heatmap_palette = "RdBu", heatmap_palcolor = NULL,
                           pseudotime_palette = "cividis", pseudotime_palcolor = NULL,
                           feature_split_palette = "jama", feature_split_palcolor = NULL,
                           cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = list(),
                           feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = list(),
                           separate_annotation = NULL, separate_annotation_palette = "Paired", separate_annotation_palcolor = NULL, separate_annotation_params = if (flip) list(width = grid::unit(2, "cm")) else list(height = grid::unit(2, "cm")),
                           reverse_ht = NULL, use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, height = NULL, width = NULL, units = "inch",
                           seed = 11, ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }

  split_method <- match.arg(split_method)
  order_by <- match.arg(order_by)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), "(", data_nm, ")")
  } else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (any(!lineages %in% colnames(srt@meta.data))) {
    lineages_missing <- lineages[!lineages %in% colnames(srt@meta.data)]
    for (l in lineages_missing) {
      if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
        pseudotime <- srt@tools[[paste0("DynamicFeatures_", l)]][["lineages"]]
        srt@meta.data[[l]] <- srt@meta.data[[pseudotime]]
      } else {
        stop("lineages: ", l, " is not in the meta data of the Seurat object")
      }
    }
  }
  if (any(!feature_from %in% lineages)) {
    stop("feature_from must be a subset of the lineages")
  }
  if (is.null(feature_split_by)) {
    feature_split_by <- lineages
  }
  if (any(!feature_split_by %in% lineages)) {
    stop("'feature_split_by' must be a subset of lineages.")
  }
  if (!split_method %in% c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime")) {
    stop("'split_method' must be one of 'mfuzz', 'kmeans', 'kmeans-peaktime', 'hclust', 'hclust-peaktime'.")
  }
  if (!is.null(feature_split) && is.null(names(feature_split))) {
    stop("'feature_split' must be named.")
  }
  if (!is.null(feature_split) && !is.factor(feature_split)) {
    feature_split <- factor(feature_split, levels = unique(feature_split))
  }
  if (is.numeric(pseudotime_label)) {
    if (length(pseudotime_label_color) == 1) {
      pseudotime_label_color <- rep(pseudotime_label_color, length(pseudotime_label))
    }
    if (length(pseudotime_label_linetype) == 1) {
      pseudotime_label_linetype <- rep(pseudotime_label_linetype, length(pseudotime_label))
    }
    if (length(pseudotime_label_linewidth) == 1) {
      pseudotime_label_linewidth <- rep(pseudotime_label_linewidth, length(pseudotime_label))
    }
    npal <- unique(c(length(pseudotime_label), length(pseudotime_label_color), length(pseudotime_label_linetype), length(pseudotime_label_linewidth)))
    if (length(npal[npal != 0]) > 1) {
      stop("Parameters for the pseudotime_label must be the same length!")
    }
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_annotation_palette) == 1) {
      cell_annotation_palette <- rep(cell_annotation_palette, length(cell_annotation))
    }
    if (length(cell_annotation_palcolor) == 1) {
      cell_annotation_palcolor <- rep(cell_annotation_palcolor, length(cell_annotation))
    }
    npal <- unique(c(length(cell_annotation_palette), length(cell_annotation_palcolor), length(cell_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("cell_annotation_palette and cell_annotation_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]])))) {
      stop("cell_annotation: ", paste0(cell_annotation[!cell_annotation %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_annotation_palette) == 1) {
      feature_annotation_palette <- rep(feature_annotation_palette, length(feature_annotation))
    }
    if (length(feature_annotation_palcolor) == 1) {
      feature_annotation_palcolor <- rep(feature_annotation_palcolor, length(feature_annotation))
    }
    npal <- unique(c(length(feature_annotation_palette), length(feature_annotation_palcolor), length(feature_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("feature_annotation_palette and feature_annotation_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features))) {
      stop("feature_annotation: ", paste0(feature_annotation[!feature_annotation %in% colnames(srt@assays[[assay]]@meta.features)], collapse = ","), " is not in the meta data of the ", assay, " assay in the Seurat object.")
    }
  }
  if (!is.null(separate_annotation)) {
    if (length(separate_annotation_palette) == 1) {
      separate_annotation_palette <- rep(separate_annotation_palette, length(separate_annotation))
    }
    if (length(separate_annotation_palcolor) == 1) {
      separate_annotation_palcolor <- rep(separate_annotation_palcolor, length(separate_annotation))
    }
    npal <- unique(c(length(separate_annotation_palette), length(separate_annotation_palcolor), length(separate_annotation)))
    if (length(npal[npal != 0]) > 1) {
      stop("separate_annotation_palette and separate_annotation_palcolor must be the same length as separate_annotation")
    }
    if (any(!unique(unlist(separate_annotation)) %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]])))) {
      stop("separate_annotation: ", paste0(unique(unlist(separate_annotation))[!unique(unlist(separate_annotation)) %in% c(colnames(srt@meta.data), rownames(srt@assays[[assay]]))], collapse = ","), " is not in the Seurat object.")
    }
  }
  if (length(width) == 1) {
    width <- rep(width, length(lineages))
  }
  if (length(height) == 1) {
    height <- rep(height, length(lineages))
  }
  if (length(width) >= 1) {
    names(width) <- lineages
  }
  if (length(height) >= 1) {
    names(height) <- lineages
  }

  column_split <- NULL
  if (isTRUE(flip)) {
    cluster_rows_raw <- cluster_rows
    cluster_columns_raw <- cluster_columns
    cluster_row_slices_raw <- cluster_row_slices
    cluster_column_slices_raw <- cluster_column_slices
    cluster_rows <- cluster_columns_raw
    cluster_columns <- cluster_rows_raw
    cluster_row_slices <- cluster_column_slices_raw
    cluster_column_slices <- cluster_row_slices_raw
  }

  cell_union <- unique(colnames(srt@assays[[1]])[apply(srt@meta.data[, lineages, drop = FALSE], 1, function(x) !all(is.na(x)))])
  Pseudotime_assign <- rowMeans(srt@meta.data[cell_union, lineages, drop = FALSE], na.rm = TRUE)
  cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union, cells = cell_union),
    Pseudotime_assign = Pseudotime_assign,
    srt@meta.data[cell_union, lineages, drop = FALSE]
  )
  if (cell_density != 1) {
    cell_bin <- cut(Pseudotime_assign, breaks = seq(min(Pseudotime_assign, na.rm = TRUE), max(Pseudotime_assign, na.rm = TRUE), length.out = 100), include.lowest = TRUE)
    ncell_bin <- ceiling(max(table(cell_bin), na.rm = TRUE) * cell_density)
    message("ncell/bin=", ncell_bin, "(100bins)")
    cell_keep <- unlist(sapply(levels(cell_bin), function(x) {
      cells <- names(Pseudotime_assign)[cell_bin == x]
      out <- sample(cells, size = min(length(cells), ncell_bin))
      return(out)
    }))
    cell_metadata <- cell_metadata[cell_keep, , drop = FALSE]
  }
  cell_order_list <- list()
  for (l in lineages) {
    if (!is.null(reverse_ht) && (l %in% lineages[reverse_ht] || l %in% reverse_ht)) {
      cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], decreasing = TRUE), , drop = FALSE]
    } else {
      cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
      cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], decreasing = FALSE), , drop = FALSE]
    }
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), l)
  }
  if (!is.null(cell_annotation)) {
    cell_metadata <- cbind.data.frame(
      cell_metadata,
      cbind.data.frame(
        srt@meta.data[rownames(cell_metadata), c(intersect(cell_annotation, colnames(srt@meta.data))), drop = FALSE],
        t(srt@assays[[assay]]@data[intersect(cell_annotation, rownames(srt@assays[[assay]])) %||% integer(), rownames(cell_metadata), drop = FALSE])
      )
    )
  }

  dynamic <- list()
  if (is.null(features)) {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
      if (is.null(DynamicFeatures)) {
        stop("DynamicFeatures result for ", l, " found in the srt object. Should perform RunDynamicFeatures first!")
      }
      DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > min_expcells & DynamicFeatures$r.sq > r.sq & DynamicFeatures$dev.expl > dev.expl & DynamicFeatures$padjust < padjust, , drop = FALSE]
      dynamic[[l]] <- DynamicFeatures
      if (l %in% feature_from) {
        features <- c(features, DynamicFeatures[["features"]])
      }
    }
    message(length(unique(features)), " features from ", paste0(lineages, collapse = ","), " passed the threshold (exp_ncells>", min_expcells, " & r.sq>", r.sq, " & dev.expl>", dev.expl, " & padjust<", padjust, "): \n", paste0(head(features, 10), collapse = ","), "...")
  } else {
    for (l in lineages) {
      DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
      if (is.null(DynamicFeatures)) {
        srt <- RunDynamicFeatures(srt, lineages = l, features = features, assay = assay, slot = slot, family = family, libsize = libsize)
        DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
      }
      if (any(!features %in% rownames(DynamicFeatures))) {
        srt <- RunDynamicFeatures(srt, lineages = l, features = features[!features %in% rownames(DynamicFeatures)], assay = assay, slot = slot, family = family, libsize = libsize)
        DynamicFeatures <- rbind(DynamicFeatures, srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]])
      }
      DynamicFeatures <- DynamicFeatures[features, , drop = FALSE]
      dynamic[[l]] <- DynamicFeatures
    }
  }

  features <- unique(features)
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  if (length(gene) == 0 && length(meta) == 0) {
    stop("No dynamic features found in the meta.data or in the assay: ", assay)
  }
  feature_metadata <- data.frame(row.names = features, features = features)
  for (l in lineages) {
    feature_metadata[rownames(dynamic[[l]]), paste0(l, order_by)] <- dynamic[[l]][, order_by]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "exp_ncells")] <- dynamic[[l]][, "exp_ncells"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "r.sq")] <- dynamic[[l]][, "r.sq"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "dev.expl")] <- dynamic[[l]][, "dev.expl"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "padjust")] <- dynamic[[l]][, "padjust"]
  }
  feature_metadata[, order_by] <- apply(feature_metadata[, paste0(lineages, order_by), drop = FALSE], 1, max, na.rm = TRUE)
  feature_metadata <- feature_metadata[order(feature_metadata[, order_by], decreasing = decreasing), , drop = FALSE]
  feature_metadata <- feature_metadata[rownames(feature_metadata) %in% features, , drop = FALSE]
  features <- rownames(feature_metadata)
  if (!is.null(feature_annotation)) {
    feature_metadata <- cbind.data.frame(feature_metadata, srt@assays[[assay]]@meta.features[rownames(feature_metadata), feature_annotation, drop = FALSE])
  }

  if (isTRUE(use_fitted)) {
    mat_list <- list()
    for (l in lineages) {
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, -1]
      rownames(fitted_matrix) <- paste0(rownames(fitted_matrix), l)
      mat_list[[l]] <- t(fitted_matrix[, features])
    }
    mat_raw <- do.call(cbind, mat_list)
  } else {
    mat_list <- list()
    Y_libsize <- colSums(slot(srt@assays[[assay]], "counts"))
    for (l in lineages) {
      cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
      mat_tmp <- as.matrix(rbind(slot(srt@assays[[assay]], slot)[gene, cells, drop = FALSE], t(srt@meta.data[cells, meta, drop = FALSE])))[features, , drop = FALSE]
      if (isTRUE(lib_normalize) && min(mat_tmp, na.rm = TRUE) >= 0) {
        if (!is.null(libsize)) {
          libsize_use <- libsize
        } else {
          libsize_use <- Y_libsize[colnames(mat_tmp)]
          isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
          if (isTRUE(isfloat)) {
            libsize_use <- rep(1, length(libsize_use))
            warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", immediate. = TRUE)
          }
        }
        mat_tmp[gene, ] <- t(t(mat_tmp[gene, , drop = FALSE]) / libsize_use * median(Y_libsize))
      }
      colnames(mat_tmp) <- paste0(colnames(mat_tmp), l)
      mat_list[[l]] <- mat_tmp
    }
    mat_raw <- do.call(cbind, mat_list)
  }

  # data used to plot heatmap
  mat <- matrix_process(mat_raw, method = exp_method)
  mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)]), na.rm = TRUE) * ifelse(mat[is.infinite(mat)] > 0, 1, -1)
  mat[is.na(mat)] <- mean(mat, na.rm = TRUE)

  # data used to do spliting
  # if ((!identical(sort(feature_split_by), sort(lineages)) && is.null(feature_split) && n_split < nrow(mat) && n_split > 1) || cell_density != 1) {
  #   mat_split <- mat_raw[, unlist(cell_order_list[feature_split_by]), drop = FALSE]
  #   mat_split <- matrix_process(mat_split, method = exp_method)
  #   mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
  #   mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  # } else {
  #   mat_split <- mat[, unlist(cell_order_list[feature_split_by]), drop = FALSE]
  # }
  mat_split <- mat[, unlist(cell_order_list[feature_split_by]), drop = FALSE]

  if (is.null(limits)) {
    if (!is.function(exp_method) && exp_method %in% c("zscore", "log2fc")) {
      b <- ceiling(min(abs(quantile(mat, c(0.01, 0.99), na.rm = TRUE)), na.rm = TRUE) * 2) / 2
      colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    } else {
      b <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
      colors <- colorRamp2(seq(b[1], b[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
    }
  } else {
    colors <- colorRamp2(seq(limits[1], limits[2], length = 100), palette_scp(palette = heatmap_palette, palcolor = heatmap_palcolor))
  }

  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_top_list <- list()
  pseudotime <- na.omit(unlist(cell_metadata[, lineages]))
  pseudotime_col <- colorRamp2(
    breaks = seq(min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE), length = 100),
    colors = palette_scp(palette = pseudotime_palette, palcolor = pseudotime_palcolor)
  )
  for (l in lineages) {
    ha_top_list[[l]] <- HeatmapAnnotation(
      Pseudotime = anno_simple(
        x = cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l],
        col = pseudotime_col,
        which = ifelse(flip, "row", "column"),
        na_col = "transparent",
        border = TRUE
      ), which = ifelse(flip, "row", "column"), show_annotation_name = l == lineages[1], annotation_name_side = ifelse(flip, "top", "left")
    )
  }
  lgd[["pseudotime"]] <- Legend(title = "Pseudotime", col_fun = pseudotime_col, border = TRUE)

  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_annotation_palette[i]
      palcolor <- cell_annotation_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = as.character(cell_anno[lineage_cells]),
            col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell, which = ifelse(flip, "row", "column"), show_annotation_name = l == lineages[1], annotation_name_side = ifelse(flip, "top", "left"))
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[cellan]] <- Legend(
          title = cellan, labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        col_fun <- colorRamp2(
          breaks = seq(min(cell_anno, na.rm = TRUE), max(cell_anno, na.rm = TRUE), length = 100),
          colors = palette_scp(palette = palette, palcolor = palcolor)
        )
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
          ha_cell <- list()
          ha_cell[[cellan]] <- anno_simple(
            x = cell_anno[lineage_cells],
            col = col_fun,
            which = ifelse(flip, "row", "column"),
            na_col = "transparent",
            border = TRUE
          )
          anno_args <- c(ha_cell, which = ifelse(flip, "row", "column"), show_annotation_name = l == lineages[1], annotation_name_side = ifelse(flip, "top", "left"))
          anno_args <- c(anno_args, cell_annotation_params[setdiff(names(cell_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[cellan]] <- Legend(
          title = cellan, col_fun = col_fun, border = TRUE
        )
      }
    }
  }

  if (!is.null(separate_annotation)) {
    subplots_list <- list()
    for (i in seq_along(separate_annotation)) {
      cellan <- separate_annotation[[i]]
      palette <- separate_annotation_palette[i]
      palcolor <- separate_annotation_palcolor[[i]]
      if (length(cellan) == 1 && cellan %in% colnames(srt@meta.data)) {
        cell_anno <- srt@meta.data[[cellan]]
      } else {
        cell_anno <- numeric()
      }
      if (!is.numeric(cell_anno)) {
        if (is.logical(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = c(TRUE, FALSE))
        } else if (!is.factor(cell_anno)) {
          cell_anno <- factor(cell_anno, levels = unique(cell_anno))
        }
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
          subplots <- CellDensityPlot(
            srt = srt, cells = lineage_cells, group.by = cellan, features = l,
            decreasing = TRUE, x_order = "rank",
            palette = palette, palcolor = palcolor,
            flip = flip, reverse = l %in% lineages[reverse_ht] || l %in% reverse_ht
          ) + theme_void()
          subplots_list[[paste0(cellan, ":", l)]] <- subplots
          graphics <- list()
          nm <- paste0(cellan, ":", l)
          funbody <- paste0(
            "
            g <- as_grob(subplots_list[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
            g$name <- '", nm, "';
            grid.draw(g);
            grid.rect(gp = gpar(fill = 'transparent', col = 'black'));
            "
          )
          funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
          eval(parse(text = paste("block_graphics <- function(index, levels) {", funbody, "}", sep = "")), envir = environment())

          ha_cell <- list()
          ha_cell[[paste0(cellan, "\n(separate)")]] <- anno_block(
            panel_fun = block_graphics,
            which = ifelse(flip, "row", "column"),
            show_name = l == lineages[1]
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, separate_annotation_params[setdiff(names(separate_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[paste0("separate:", cellan)]] <- Legend(
          title = paste0(cellan, "\n(separate)"), labels = levels(cell_anno),
          legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        for (l in lineages) {
          lineage_cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
          subplots <- DynamicPlot(
            srt = srt, cells = lineage_cells, lineages = l, group.by = NULL, features = cellan,
            line_palette = palette, line_palcolor = palcolor,
            add_rug = FALSE, legend.position = "none", compare_features = TRUE, x_order = "rank",
            flip = flip, reverse = l %in% lineages[reverse_ht] || l %in% reverse_ht
          ) + theme_void()
          subplots_list[[paste0(paste0(cellan, collapse = ","), ":", l)]] <- subplots
          graphics <- list()
          nm <- paste0(paste0(cellan, collapse = ","), ":", l)
          funbody <- paste0(
            "
            g <- as_grob(subplots_list[['", nm, "']] + theme_void() + theme(legend.position = 'none'));
            g$name <- '", nm, "';
            grid.draw(g);
            grid.rect(gp = gpar(fill = 'transparent', col = 'black'));
            "
          )
          funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
          eval(parse(text = paste("block_graphics <- function(index, levels) {", funbody, "}", sep = "")), envir = environment())

          ha_cell <- list()
          ha_cell[[paste0(paste0(cellan, collapse = ","), "\n(separate)")]] <- anno_block(
            panel_fun = block_graphics,
            which = ifelse(flip, "row", "column"),
            show_name = l == lineages[1]
          )
          anno_args <- c(ha_cell,
            which = ifelse(flip, "row", "column"),
            show_annotation_name = l == lineages[1],
            annotation_name_side = ifelse(flip, "top", "left")
          )
          anno_args <- c(anno_args, separate_annotation_params[setdiff(names(separate_annotation_params), names(anno_args))])
          ha_top <- do.call(HeatmapAnnotation, args = anno_args)
          if (is.null(ha_top_list[[l]])) {
            ha_top_list[[l]] <- ha_top
          } else {
            ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
          }
        }
        lgd[[paste0("separate:", paste0(cellan, collapse = ","))]] <- Legend(
          title = "Features\n(separate)", labels = cellan,
          legend_gp = gpar(fill = palette_scp(cellan, palette = palette, palcolor = palcolor)), border = TRUE
        )
      }
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat) <= n_split)) {
      row_split_raw <- row_split <- feature_split <- NULL
    } else {
      if (n_split == 1) {
        row_split_raw <- row_split <- feature_split <- setNames(rep(1, nrow(mat_split)), rownames(mat_split))
      } else {
        if (split_method == "mfuzz") {
          status <- tryCatch(check_R("Mfuzz"), error = identity)
          if (inherits(status, "error")) {
            warning("The Mfuzz package was not found. Switch split_method to 'kmeans'", immediate. = TRUE)
            split_method <- "kmeans"
          } else {
            require("Mfuzz", quietly = TRUE)
            eset <- new("ExpressionSet", exprs = mat_split)
            eset <- Mfuzz::standardise(eset)
            min_fuzzification <- Mfuzz::mestimate(eset)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            } else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- Mfuzz::mfuzz(eset, c = n_split, m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
            }
            if (isTRUE(show_fuzzification)) {
              message("fuzzification: ", fuzzification)
            }
            # mfuzz.plot(eset, cl,new.window = FALSE)
            row_split <- feature_split <- cl$cluster
          }
        }
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, iter.max = 1e4, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "kmeans-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          km <- kmeans(feature_y, centers = n_split, iter.max = 1e4, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(stats::dist(mat_split))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
        if (split_method == "hclust-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          hc <- hclust(stats::dist(feature_y))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      df <- data.frame(row_split = row_split, order_by = feature_metadata[names(row_split), order_by])
      df_order <- aggregate(df, by = list(row_split), FUN = mean)
      df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), , drop = FALSE]
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split_raw <- row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))
    }
  } else {
    row_split_raw <- row_split <- feature_split <- feature_split[row.names(mat)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
  } else {
    feature_metadata[["feature_split"]] <- NA
  }

  ha_left <- NULL
  if (!is.null(row_split)) {
    if (isTRUE(cluster_row_slices)) {
      if (!isTRUE(cluster_rows)) {
        dend <- cluster_within_group(t(mat_split), row_split_raw)
        cluster_rows <- dend
        row_split <- length(unique(row_split_raw))
      }
    }
    funbody <- paste0(
      "
      grid.rect(gp = gpar(fill = palette_scp(", paste0("c('", paste0(levels(row_split_raw), collapse = "','"), "')"), ",palette = '", feature_split_palette, "',palcolor=c(", paste0("'", paste0(unlist(feature_split_palcolor), collapse = "','"), "'"), "))[nm]))
    "
    )
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)
    eval(parse(text = paste("panel_fun <- function(index, nm) {", funbody, "}", sep = "")), envir = environment())
    ha_clusters <- HeatmapAnnotation(
      features_split = anno_block(
        align_to = split(seq_along(row_split_raw), row_split_raw),
        panel_fun = getFunction("panel_fun", where = environment()),
        width = unit(0.1, "in"),
        height = unit(0.1, "in"),
        show_name = FALSE,
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"),
      border = TRUE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(
      title = "Cluster", labels = intersect(levels(row_split_raw), row_split_raw),
      legend_gp = gpar(fill = palette_scp(intersect(levels(row_split_raw), row_split_raw), type = "discrete", palette = feature_split_palette, palcolor = feature_split_palcolor, matched = TRUE)), border = TRUE
    )
  }

  if (isTRUE(cluster_rows) && !is.null(cluster_features_by)) {
    mat_cluster <- mat[, unlist(cell_order_list[cluster_features_by]), drop = FALSE]
    if (is.null(row_split)) {
      dend <- as.dendrogram(hclust(dist(mat_cluster)))
      dend_ordered <- reorder(dend, wts = colMeans(mat_cluster), agglo.FUN = mean)
      cluster_rows <- dend_ordered
    } else {
      row_split <- length(unique(row_split_raw))
      dend <- cluster_within_group2(t(mat_cluster), row_split_raw)
      cluster_rows <- dend
    }
  }

  l <- lineages[1]
  ht_args <- list(
    matrix = mat[, cell_order_list[[l]], drop = FALSE],
    col = colors,
    row_split = row_split,
    cluster_rows = cluster_rows,
    cluster_row_slices = cluster_row_slices,
    column_split = column_split,
    cluster_columns = cluster_columns,
    cluster_column_slices = cluster_column_slices,
    use_raster = TRUE
  )
  ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
  ht_list <- do.call(Heatmap, args = ht_args)
  features_ordered <- rownames(mat)[unlist(suppressWarnings(row_order(ht_list)))]
  feature_metadata[["index"]] <- setNames(object = seq_along(features_ordered), nm = features_ordered)[rownames(feature_metadata)]

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index_from <- ceiling((length(features_ordered) / nlabel) / 2)
        index_to <- length(features_ordered)
        index <- unique(round(seq(from = index_from, to = index_to, length.out = nlabel)))
      } else {
        index <- seq_along(features_ordered)
      }
    } else {
      index <- NULL
    }
  } else {
    index <- which(features_ordered %in% features_label)
    drop <- setdiff(features_label, features_ordered)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the features"), immediate. = TRUE)
    }
  }
  if (length(index) > 0) {
    ha_mark <- HeatmapAnnotation(
      gene = anno_mark(
        at = which(rownames(feature_metadata) %in% features_ordered[index]),
        labels = feature_metadata[which(rownames(feature_metadata) %in% features_ordered[index]), "features"],
        side = ifelse(flip, "top", "left"),
        labels_gp = gpar(fontsize = label_size, col = label_color),
        link_gp = gpar(fontsize = label_size, col = label_color),
        which = ifelse(flip, "column", "row")
      ),
      which = ifelse(flip, "column", "row"), show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_mark
    } else {
      ha_left <- c(ha_mark, ha_left)
    }
  }

  ha_right <- NULL
  if (length(lineages) > 1) {
    ha_list <- list()
    for (l in lineages) {
      ha_list[[l]] <- anno_simple(
        x = is.na(feature_metadata[, paste0(l, order_by)]) + 0,
        col = c("0" = "#181830", "1" = "transparent"),
        width = unit(0.5, "cm"),
        height = unit(0.5, "cm"),
        which = ifelse(flip, "column", "row")
      )
    }
    ha_lineage <- do.call("HeatmapAnnotation", args = c(ha_list, which = ifelse(flip, "column", "row"), annotation_name_side = ifelse(flip, "left", "top"), border = TRUE))
    if (is.null(ha_right)) {
      ha_right <- ha_lineage
    } else {
      ha_right <- c(ha_right, ha_lineage)
    }
  }
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featan <- feature_annotation[i]
      palette <- feature_annotation_palette[i]
      palcolor <- feature_annotation_palcolor[[i]]
      featan_values <- feature_metadata[, featan]
      if (!is.numeric(featan_values)) {
        if (is.logical(featan_values)) {
          featan_values <- factor(featan_values, levels = c(TRUE, FALSE))
        } else if (!is.factor(featan_values)) {
          featan_values <- factor(featan_values, levels = unique(featan_values))
        }
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = as.character(featan_values),
          col = palette_scp(featan_values, palette = palette, palcolor = palcolor),
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, labels = levels(featan_values),
          legend_gp = gpar(fill = palette_scp(featan_values, palette = palette, palcolor = palcolor)), border = TRUE
        )
      } else {
        col_fun <- colorRamp2(
          breaks = seq(min(featan_values, na.rm = TRUE), max(featan_values, na.rm = TRUE), length = 100),
          colors = palette_scp(palette = palette, palcolor = palcolor)
        )
        ha_feature <- list()
        ha_feature[[featan]] <- anno_simple(
          x = featan_values,
          col = col_fun,
          which = ifelse(flip, "column", "row"),
          na_col = "transparent",
          border = TRUE
        )
        anno_args <- c(ha_feature, which = ifelse(flip, "column", "row"), show_annotation_name = TRUE, annotation_name_side = ifelse(flip, "left", "top"), border = TRUE)
        anno_args <- c(anno_args, feature_annotation_params[setdiff(names(feature_annotation_params), names(anno_args))])
        ha_feature <- do.call(HeatmapAnnotation, args = anno_args)
        if (is.null(ha_right)) {
          ha_right <- ha_feature
        } else {
          ha_right <- c(ha_right, ha_feature)
        }
        lgd[[featan]] <- Legend(
          title = featan, col_fun = col_fun, border = TRUE
        )
      }
    }
  }

  enrichment <- heatmap_enrichment(
    geneID = feature_metadata[["features"]], geneID_groups = feature_metadata[["feature_split"]],
    feature_split_palette = feature_split_palette, feature_split_palcolor = feature_split_palcolor,
    ha_right = ha_right, flip = flip,
    anno_terms = anno_terms, anno_keys = anno_keys, anno_features = anno_features,
    terms_width = terms_width, terms_fontsize = terms_fontsize,
    keys_width = keys_width, keys_fontsize = keys_fontsize,
    features_width = features_width, features_fontsize = features_fontsize,
    IDtype = IDtype, species = species, db_update = db_update, db_version = db_version, db_combine = db_combine, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror,
    db = db, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize,
    GO_simplify = GO_simplify, GO_simplify_cutoff = GO_simplify_cutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff,
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, min_word_length = min_word_length,
    exclude_words = exclude_words
  )
  res <- enrichment$res
  ha_right <- enrichment$ha_right

  ht_list <- NULL
  for (l in lineages) {
    if (l == lineages[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (l == lineages[length(lineages)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }
    ht_args <- list(
      name = l,
      matrix = if (flip) t(mat[, cell_order_list[[l]], drop = FALSE]) else mat[, cell_order_list[[l]], drop = FALSE],
      col = colors,
      row_title = if (flip) l else character(0),
      row_title_side = row_title_side,
      column_title = if (flip) character(0) else l,
      column_title_side = if (flip) "top" else column_title_side,
      row_title_rot = row_title_rot,
      column_title_rot = column_title_rot,
      row_split = if (flip) column_split else row_split,
      column_split = if (flip) row_split else column_split,
      cluster_rows = if (flip) cluster_columns else cluster_rows,
      cluster_columns = if (flip) cluster_rows else cluster_columns,
      cluster_row_slices = if (flip) cluster_column_slices else cluster_row_slices,
      cluster_column_slices = if (flip) cluster_row_slices else cluster_column_slices,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_rot = row_names_rot,
      column_names_rot = column_names_rot,
      top_annotation = if (flip) left_annotation else ha_top_list[[l]],
      left_annotation = if (flip) ha_top_list[[l]] else left_annotation,
      bottom_annotation = if (flip) right_annotation else NULL,
      right_annotation = if (flip) NULL else right_annotation,
      show_heatmap_legend = FALSE,
      border = border,
      use_raster = use_raster,
      raster_device = raster_device,
      raster_by_magick = raster_by_magick,
      width = if (is.numeric(width[l])) unit(width[l], units = units) else NULL,
      height = if (is.numeric(height[l])) unit(height[l], units = units) else NULL
    )
    if (any(names(ht_params) %in% names(ht_args))) {
      warning("ht_params: ", paste0(intersect(names(ht_params), names(ht_args)), collapse = ","), " were duplicated and will not be used.", immediate. = TRUE)
    }
    ht_args <- c(ht_args, ht_params[setdiff(names(ht_params), names(ht_args))])
    if (isTRUE(flip)) {
      if (is.null(ht_list)) {
        ht_list <- do.call(Heatmap, args = ht_args)
      } else {
        ht_list <- ht_list %v% do.call(Heatmap, args = ht_args)
      }
    } else {
      ht_list <- ht_list + do.call(Heatmap, args = ht_args)
    }
  }

  if ((!is.null(row_split) && length(index) > 0) || any(c(anno_terms, anno_keys, anno_features)) || !is.null(width) || !is.null(height)) {
    fix <- TRUE
  } else {
    fix <- FALSE
  }
  rendersize <- heatmap_rendersize(
    width = width, height = height, units = units,
    ha_top_list = ha_top_list, ha_left = ha_left, ha_right = ha_right,
    ht_list = ht_list, legend_list = lgd, flip = flip
  )
  width_sum <- rendersize[["width_sum"]]
  height_sum <- rendersize[["height_sum"]]
  # cat("width:", width, "\n")
  # cat("height:", height, "\n")
  # cat("width_sum:", width_sum, "\n")
  # cat("height_sum:", height_sum, "\n")

  if (isTRUE(fix)) {
    fixsize <- heatmap_fixsize(
      width = width, width_sum = width_sum, height = height, height_sum = height_sum, units = units,
      ht_list = ht_list, legend_list = lgd
    )
    ht_width <- fixsize[["ht_width"]]
    ht_height <- fixsize[["ht_height"]]
    # cat("ht_width:", ht_width, "\n")
    # cat("ht_height:", ht_height, "\n")

    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
        if (is.numeric(pseudotime_label)) {
          for (n in seq_along(pseudotime_label)) {
            pse <- pseudotime_label[n]
            col <- pseudotime_label_color[n]
            lty <- pseudotime_label_linetype[n]
            lwd <- pseudotime_label_linewidth[n]
            for (l in lineages) {
              for (slice in 1:max(nlevels(row_split), 1)) {
                decorate_heatmap_body(l,
                  {
                    pseudotime <- cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l]
                    i <- which.min(abs(pseudotime - pse))
                    if (flip) {
                      x <- 1 - (i / length(pseudotime))
                      grid.lines(c(0, 1), c(x, x), gp = gpar(lty = lty, lwd = lwd, col = col))
                    } else {
                      i <- which.min(abs(pseudotime - pse))
                      x <- i / length(pseudotime)
                      grid.lines(c(x, x), c(0, 1), gp = gpar(lty = lty, lwd = lwd, col = col))
                    }
                  },
                  row_slice = ifelse(flip, 1, slice),
                  column_slice = ifelse(flip, slice, 1)
                )
              }
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  } else {
    ht_width <- unit(width_sum, units = units)
    ht_height <- unit(height_sum, units = units)
    gTree <- grid.grabExpr(
      {
        draw(ht_list, annotation_legend_list = lgd)
        for (enrich in db) {
          enrich_anno <- names(ha_right)[grep(paste0("_split_", enrich), names(ha_right))]
          if (length(enrich_anno) > 0) {
            for (enrich_anno_element in enrich_anno) {
              enrich_obj <- strsplit(enrich_anno_element, "_split_")[[1]][1]
              decorate_annotation(enrich_anno_element, slice = 1, {
                grid.text(paste0(enrich, " (", enrich_obj, ")"), x = unit(1, "npc"), y = unit(1, "npc") + unit(2.5, "mm"), just = c("left", "bottom"))
              })
            }
          }
        }
        if (is.numeric(pseudotime_label)) {
          for (n in seq_along(pseudotime_label)) {
            pse <- pseudotime_label[n]
            col <- pseudotime_label_color[n]
            lty <- pseudotime_label_linetype[n]
            lwd <- pseudotime_label_linewidth[n]
            for (l in lineages) {
              for (slice in 1:max(nlevels(row_split), 1)) {
                decorate_heatmap_body(l,
                  {
                    pseudotime <- cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l]
                    i <- which.min(abs(pseudotime - pse))
                    if (flip) {
                      x <- 1 - (i / length(pseudotime))
                      grid.lines(c(0, 1), c(x, x), gp = gpar(lty = lty, lwd = lwd, col = col))
                    } else {
                      i <- which.min(abs(pseudotime - pse))
                      x <- i / length(pseudotime)
                      grid.lines(c(x, x), c(0, 1), gp = gpar(lty = lty, lwd = lwd, col = col))
                    }
                  },
                  row_slice = ifelse(flip, 1, slice),
                  column_slice = ifelse(flip, slice, 1)
                )
              }
            }
          }
        }
      },
      width = ht_width,
      height = ht_height,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
  }

  if (isTRUE(fix)) {
    p <- panel_fix_single(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
  } else {
    p <- wrap_plots(gTree)
  }

  return(list(
    plot = p,
    matrix = mat,
    cell_order = cell_order_list,
    feature_split = feature_split,
    cell_metadata = cell_metadata,
    feature_metadata = feature_metadata,
    enrichment = res
  ))
}

#' Plot for dynamic features along lineages.
#'
#' @param srt
#' @param features
#' @param lineages
#' @param slot
#' @param assay
#' @param family
#' @param libsize
#' @param exp_method
#' @param lib_normalize
#' @param group.by
#' @param compare_lineages
#' @param compare_features
#' @param add_line
#' @param add_interval
#' @param line.size
#' @param line_palette
#' @param line_palcolor
#' @param add_point
#' @param pt.size
#' @param point_palette
#' @param point_palcolor
#' @param add_rug
#' @param aspect.ratio
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param cells
#' @param flip
#' @param reverse
#' @param x_order
#' @param theme_use
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   features = c("Nnat", "Irx1", "G2M_score"),
#'   group.by = "SubCellType",
#'   compare_features = TRUE
#' )
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1", "G2M_score"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1", "G2M_score"),
#'   group.by = "SubCellType",
#'   compare_lineages = FALSE,
#'   compare_features = FALSE
#' )
#' @importFrom ggplot2 geom_line geom_point geom_ribbon geom_rug stat_density2d facet_grid expansion
#' @importFrom reshape2 melt
#' @importFrom patchwork wrap_plots
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom grDevices colorRampPalette
#' @importFrom stats runif
#' @export
DynamicPlot <- function(srt, features, lineages, group.by = NULL, cells = NULL, slot = "counts", assay = NULL, family = NULL,
                        exp_method = c("log1p", "raw", "zscore", "fc", "log2fc"), lib_normalize = identical(slot, "counts"), libsize = NULL,
                        compare_lineages = TRUE, compare_features = FALSE,
                        add_line = TRUE, add_interval = TRUE, line.size = 1, line_palette = "Dark2", line_palcolor = NULL,
                        add_point = TRUE, pt.size = 1, point_palette = "Paired", point_palcolor = NULL,
                        add_rug = TRUE, flip = FALSE, reverse = FALSE, x_order = c("value", "rank"),
                        aspect.ratio = NULL,
                        legend.position = "right", legend.direction = "vertical",
                        theme_use = "theme_scp", theme_args = list(),
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11) {
  set.seed(seed)

  check_R("MatrixGenerics")
  x_order <- match.arg(x_order)
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    stop(group.by, " is not in the meta.data of srt object.")
  }

  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  if (length(exp_method) == 1 && is.function(exp_method)) {
    exp_name <- paste0(as.character(x = formals()$exp_method), "(", data_nm, ")")
  } else {
    exp_method <- match.arg(exp_method)
    exp_name <- paste0(exp_method, "(", data_nm, ")")
  }

  assay <- assay %||% DefaultAssay(srt)
  gene <- features[features %in% rownames(srt@assays[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  features <- c(gene, meta)
  if (length(features) == 0) {
    stop("No feature found in the srt object.")
  }

  cell_union <- c()
  raw_matrix_list <- list()
  fitted_matrix_list <- list()
  upr_matrix_list <- list()
  lwr_matrix_list <- list()
  for (l in lineages) {
    features_exist <- c()
    raw_matrix <- NULL
    if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      raw_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]][, -1]
      fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, -1]
      upr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["upr_matrix"]][, -1]
      lwr_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["lwr_matrix"]][, -1]
      features_exist <- colnames(raw_matrix)
    }
    feature_calcu <- features[!features %in% features_exist]
    if (length(feature_calcu) > 0) {
      srt_tmp <- RunDynamicFeatures(srt, lineages = l, features = feature_calcu, assay = assay, slot = slot, family = family, libsize = libsize)
      if (is.null(raw_matrix)) {
        raw_matrix <- fitted_matrix <- upr_matrix <- lwr_matrix <- matrix(NA, nrow = nrow(srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]), ncol = 0)
      }
      raw_matrix <- cbind(raw_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]][, feature_calcu, drop = FALSE])
      fitted_matrix <- cbind(fitted_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, feature_calcu, drop = FALSE])
      upr_matrix <- cbind(upr_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["upr_matrix"]][, feature_calcu, drop = FALSE])
      lwr_matrix <- cbind(lwr_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["lwr_matrix"]][, feature_calcu, drop = FALSE])
    }
    raw_matrix_list[[l]] <- as.matrix(raw_matrix[, features, drop = FALSE])
    fitted_matrix_list[[l]] <- as.matrix(fitted_matrix[, features, drop = FALSE])
    upr_matrix_list[[l]] <- as.matrix(upr_matrix[, features, drop = FALSE])
    lwr_matrix_list[[l]] <- as.matrix(lwr_matrix[, features, drop = FALSE])
    cell_union <- unique(c(cell_union, rownames(raw_matrix)))
  }

  x_assign <- rowMeans(srt@meta.data[cell_union, lineages, drop = FALSE], na.rm = TRUE)
  cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union),
    x_assign = x_assign,
    srt@meta.data[cell_union, lineages, drop = FALSE]
  )

  cell_order_list <- list()
  for (l in lineages) {
    cell_metadata_sub <- na.omit(cell_metadata[, l, drop = FALSE])
    cell_metadata_sub <- cell_metadata_sub[order(cell_metadata_sub[[l]], decreasing = FALSE), , drop = FALSE]
    cell_order_list[[l]] <- paste0(rownames(cell_metadata_sub), l)
  }

  df_list <- list()
  Y_libsize <- colSums(slot(srt@assays[[assay]], "counts"))
  for (l in lineages) {
    raw_matrix <- raw_matrix_list[[l]]
    fitted_matrix <- fitted_matrix_list[[l]]
    upr_matrix <- upr_matrix_list[[l]]
    lwr_matrix <- lwr_matrix_list[[l]]
    if (isTRUE(lib_normalize) && min(raw_matrix[, gene], na.rm = TRUE) >= 0) {
      if (!is.null(libsize)) {
        libsize_use <- libsize
      } else {
        libsize_use <- Y_libsize[rownames(raw_matrix)]
        isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
        if (isTRUE(isfloat)) {
          libsize_use <- rep(1, length(libsize_use))
          warning("The values in the 'counts' slot are non-integer. Set the library size to 1.", immediate. = TRUE)
        }
      }
      raw_matrix[, gene] <- raw_matrix[, gene, drop = FALSE] / libsize_use * median(Y_libsize)
    }

    if (is.function(exp_method)) {
      raw_matrix <- t(exp_method(t(raw_matrix)))
      fitted_matrix <- t(exp_method(t(fitted_matrix)))
      upr_matrix <- t(exp_method(t(upr_matrix)))
      lwr_matrix <- t(exp_method(t(lwr_matrix)))
    } else if (exp_method == "raw") {
      raw_matrix <- raw_matrix
      fitted_matrix <- fitted_matrix
      upr_matrix <- upr_matrix
      lwr_matrix <- lwr_matrix
    } else if (exp_method == "zscore") {
      center <- colMeans(raw_matrix)
      sd <- MatrixGenerics::colSds(raw_matrix)
      raw_matrix <- scale(raw_matrix, center = center, scale = sd)
      fitted_matrix <- scale(fitted_matrix, center = center, scale = sd)
      upr_matrix <- scale(upr_matrix, center = center, scale = sd)
      lwr_matrix <- scale(lwr_matrix, center = center, scale = sd)
    } else if (exp_method == "fc") {
      colm <- colMeans(raw_matrix)
      raw_matrix <- t(t(raw_matrix) / colm)
      fitted_matrix <- t(t(fitted_matrix) / colm)
      upr_matrix <- t(t(upr_matrix) / colm)
      lwr_matrix <- t(t(lwr_matrix) / colm)
    } else if (exp_method == "log2fc") {
      colm <- colMeans(raw_matrix)
      raw_matrix <- t(log2(t(raw_matrix) / colm))
      fitted_matrix <- t(log2(t(fitted_matrix) / colm))
      upr_matrix <- t(log2(t(upr_matrix) / colm))
      lwr_matrix <- t(log2(t(lwr_matrix) / colm))
    } else if (exp_method == "log1p") {
      raw_matrix <- log1p(raw_matrix)
      fitted_matrix <- log1p(fitted_matrix)
      upr_matrix <- log1p(upr_matrix)
      lwr_matrix <- log1p(lwr_matrix)
    }
    raw_matrix[is.infinite(raw_matrix)] <- max(abs(raw_matrix[!is.infinite(raw_matrix)]), na.rm = TRUE) * ifelse(raw_matrix[is.infinite(raw_matrix)] > 0, 1, -1)
    fitted_matrix[is.infinite(fitted_matrix)] <- max(abs(fitted_matrix[!is.infinite(fitted_matrix)])) * ifelse(fitted_matrix[is.infinite(fitted_matrix)] > 0, 1, -1)
    upr_matrix[is.infinite(upr_matrix)] <- max(abs(upr_matrix[!is.infinite(upr_matrix)]), na.rm = TRUE) * ifelse(upr_matrix[is.infinite(upr_matrix)] > 0, 1, -1)
    lwr_matrix[is.infinite(lwr_matrix)] <- max(abs(lwr_matrix[!is.infinite(lwr_matrix)]), na.rm = TRUE) * ifelse(lwr_matrix[is.infinite(lwr_matrix)] > 0, 1, -1)

    raw <- as.data.frame(cbind(cell_metadata[rownames(raw_matrix), c(l, "x_assign")], raw_matrix))
    colnames(raw)[1] <- "Pseudotime"
    raw[["Cell"]] <- rownames(raw)
    raw[["Value"]] <- "raw"
    raw <- reshape2::melt(raw, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    fitted <- as.data.frame(cbind(cell_metadata[rownames(fitted_matrix), c(l, "x_assign")], fitted_matrix))
    colnames(fitted)[1] <- "Pseudotime"
    fitted[["Cell"]] <- rownames(fitted)
    fitted[["Value"]] <- "fitted"
    fitted <- reshape2::melt(fitted, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    upr <- as.data.frame(cbind(cell_metadata[rownames(upr_matrix), c(l, "x_assign")], upr_matrix))
    colnames(upr)[1] <- "Pseudotime"
    upr[["Cell"]] <- rownames(upr)
    upr[["Value"]] <- "upr"
    upr <- reshape2::melt(upr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    lwr <- as.data.frame(cbind(cell_metadata[rownames(lwr_matrix), c(l, "x_assign")], lwr_matrix))
    colnames(lwr)[1] <- "Pseudotime"
    lwr[["Cell"]] <- rownames(lwr)
    lwr[["Value"]] <- "lwr"
    lwr <- reshape2::melt(lwr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    raw[["upr"]] <- NA
    raw[["lwr"]] <- NA
    fitted[["upr"]] <- upr[["exp"]]
    fitted[["lwr"]] <- lwr[["exp"]]

    df_tmp <- rbind(raw, fitted)
    df_tmp[["Lineages"]] <- factor(l, levels = lineages)
    df_list[[l]] <- df_tmp
  }
  df_all <- do.call(rbind, df_list)
  rownames(df_all) <- NULL

  if (!is.null(group.by)) {
    cell_group <- srt@meta.data[df_all[["Cell"]], group.by, drop = FALSE]
    if (!is.factor(cell_group[, group.by])) {
      cell_group[, group.by] <- factor(cell_group[, group.by], levels = unique(cell_group[, group.by]))
    }
    df_all <- cbind(df_all, cell_group)
  }
  df_all[["LineagesFeatures"]] <- paste(df_all[["Lineages"]], df_all[["Features"]], sep = "-")

  if (!is.null(cells)) {
    df_all <- df_all[df_all[["Cell"]] %in% cells, , drop = FALSE]
  }
  df_all <- df_all[sample(seq_len(nrow(df_all))), , drop = FALSE]

  plist <- list()
  legend <- NULL
  if (isTRUE(compare_lineages)) {
    lineages_use <- list(lineages)
    lineages_formula <- "."
  } else {
    lineages_use <- lineages
    lineages_formula <- "Lineages"
  }
  if (isTRUE(compare_features)) {
    features_use <- list(features)
    features_formula <- "."
  } else {
    features_use <- features
    features_formula <- "Features"
  }
  formula <- paste(lineages_formula, "~", features_formula)
  fill_by <- "Lineages"
  if (lineages_formula == "." && length(lineages) > 1) {
    lineages_guide <- TRUE
  } else {
    lineages_guide <- FALSE
    if (isTRUE(compare_features)) {
      fill_by <- "Features"
    }
  }
  if (features_formula == "." && length(features) > 1) {
    features_guide <- TRUE
  } else {
    features_guide <- FALSE
  }

  for (l in lineages_use) {
    for (f in features_use) {
      df <- subset(df_all, df_all[["Lineages"]] %in% l & df_all[["Features"]] %in% f)
      if (x_order == "rank") {
        df[, "x_assign"] <- rank(df[, "x_assign"])
        df[, "Pseudotime"] <- rank(df[, "Pseudotime"])
      }
      df_point <- unique(df[df[["Value"]] == "raw", c("Cell", "x_assign", "exp", group.by)])
      if (isTRUE(compare_features)) {
        raw_point <- NULL
      } else {
        if (isTRUE(add_point)) {
          if (is.null(group.by)) {
            raw_point <- geom_point(data = df_point, mapping = aes(x = .data[["x_assign"]], y = .data[["exp"]]), size = pt.size, alpha = 0.8)
          } else {
            raw_point <- list(
              geom_point(data = df_point, mapping = aes(x = .data[["x_assign"]], y = .data[["exp"]], color = .data[[group.by]]), size = pt.size, alpha = 0.8),
              scale_color_manual(
                values = palette_scp(df[[group.by]], palette = point_palette, palcolor = point_palcolor)
              ),
              scale_fill_manual(
                values = palette_scp(df[[group.by]], palette = point_palette, palcolor = point_palcolor),
                guide = guide_legend(override.aes = list(alpha = 1, size = 3), order = 1)
              ),
              new_scale_color(),
              new_scale_fill()
            )
          }
        } else {
          raw_point <- NULL
        }
      }
      if (isTRUE(add_rug)) {
        if (is.null(group.by)) {
          rug <- list(geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]]), alpha = 1, length = unit(0.05, "npc"), show.legend = FALSE))
        } else {
          rug <- list(
            geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]], color = .data[[group.by]]), alpha = 1, length = unit(0.05, "npc"), show.legend = isTRUE(compare_features)),
            scale_color_manual(
              values = palette_scp(df[[group.by]], palette = point_palette, palcolor = point_palcolor)
            ),
            new_scale_color()
          )
        }
      } else {
        rug <- NULL
      }

      if (isTRUE(add_interval)) {
        interval <- list(
          geom_ribbon(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(
              x = .data[["Pseudotime"]], y = .data[["exp"]], ymin = .data[["lwr"]], ymax = .data[["upr"]],
              fill = .data[[fill_by]], group = .data[["LineagesFeatures"]]
            ),
            alpha = 0.4, color = "grey90"
          ),
          scale_fill_manual(
            values = palette_scp(df[[fill_by]], palette = line_palette, palcolor = line_palcolor),
            guide = if (fill_by == "Features" || lineages_guide || length(l) == 1) "none" else guide_legend()
          ),
          new_scale_fill()
        )
      } else {
        interval <- NULL
      }
      if (isTRUE(compare_features)) {
        line <- list(
          geom_line(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(
              x = .data[["Pseudotime"]], y = .data[["exp"]],
              color = .data[["Features"]], group = .data[["LineagesFeatures"]]
            ),
            linewidth = line.size, alpha = 0.8
          ),
          scale_color_manual(
            values = palette_scp(df[["Features"]], palette = line_palette, palcolor = line_palcolor),
            guide = if (features_guide) guide_legend(override.aes = list(alpha = 1, size = 2), order = 2) else "none"
          ),
          new_scale_color()
        )
      } else {
        if (isTRUE(add_line)) {
          line <- list(
            geom_line(
              data = subset(df, df[["Value"]] == "fitted"),
              mapping = aes(x = .data[["Pseudotime"]], y = .data[["exp"]], color = .data[["Lineages"]], group = .data[["LineagesFeatures"]]), linewidth = line.size, alpha = 0.8
            ),
            scale_color_manual(
              values = palette_scp(df[["Lineages"]], palette = line_palette, palcolor = line_palcolor),
              guide = if (lineages_guide) guide_legend(override.aes = list(alpha = 1, size = 2), order = 2) else "none"
            ),
            new_scale_color()
          )
        } else {
          line <- NULL
        }
      }

      x_trans <- ifelse(flip, "reverse", "identity")
      x_trans <- ifelse(reverse, setdiff(c("reverse", "identity"), x_trans), x_trans)
      p <- ggplot() +
        scale_x_continuous(trans = x_trans, expand = expansion(c(0, 0))) +
        scale_y_continuous(expand = expansion(c(0.1, 0.05))) +
        raw_point +
        rug +
        interval +
        line +
        labs(x = ifelse(x_order == "rank", "Pseudotime(rank)", "Pseudotime"), y = exp_name) +
        facet_grid(formula(formula), scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )

      if (isTRUE(flip)) {
        p <- p + coord_flip()
      }
      if (is.null(legend)) {
        legend <- get_legend(p + theme(legend.position = "bottom"))
      }
      plist[[paste(paste0(l, collapse = "_"), paste0(f, collapse = "_"), sep = ".")]] <- p + theme(legend.position = "none")
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    if (legend.position != "none") {
      gtable <- as_grob(plot)
      gtable <- add_grob(gtable, legend, legend.position)
      plot <- wrap_plots(gtable)
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' Projection Plot
#'
#' @param srt_query
#' @param srt_ref
#' @param query_group
#' @param ref_group
#' @param query_reduction
#' @param ref_reduction
#' @param query_param
#' @param ref_param
#' @param xlim
#' @param ylim
#' @param pt.size
#' @param stroke.highlight
#'
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous ggplot_build theme geom_point aes scale_fill_identity facet_null
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom grid grob
#' @importFrom rlang  %||%
#' @importFrom patchwork wrap_plots
#' @export
ProjectionPlot <- function(srt_query, srt_ref,
                           query_group = "orig.ident", ref_group = "orig.ident",
                           query_reduction = "ref.umap", ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||% NULL,
                           query_param = list(palette = "Dark2"), ref_param = list(palette = "Paired"),
                           xlim = NULL, ylim = NULL, pt.size = 0.8, stroke.highlight = 0.5) {
  if (is.null(ref_reduction)) {
    stop("Please specify the ref_reduction.")
  }
  if (is.null(xlim)) {
    ref_xlim <- range(srt_ref[[ref_reduction]]@cell.embeddings[, 1])
    query_xlim <- range(srt_query[[query_reduction]]@cell.embeddings[, 1])
    xlim <- range(c(ref_xlim, query_xlim))
  }
  if (is.null(ylim)) {
    ref_ylim <- range(srt_ref[[ref_reduction]]@cell.embeddings[, 2])
    query_ylim <- range(srt_query[[query_reduction]]@cell.embeddings[, 2])
    ylim <- range(c(ref_ylim, query_ylim))
  }

  p1 <- do.call(CellDimPlot, args = c(
    srt = srt_ref, reduction = ref_reduction, group.by = ref_group,
    ref_param
  )) +
    guides(color = guide_legend(title = paste0("Ref: ", ref_group), override.aes = list(size = 4.5)))
  p1legend <- get_legend(p1)
  # p1legend <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")

  suppressMessages(expr = {
    p2 <- do.call(CellDimPlot, args = c(
      srt = srt_query, reduction = query_reduction, group.by = query_group,
      query_param
    )) +
      scale_x_continuous(limits = xlim) +
      scale_y_continuous(limits = ylim)
  })
  p2data <- ggplot_build(p2)$data[[1]]
  color <- p2data$colour
  names(color) <- p2$data$group.by
  p2 <- p2 + guides(color = guide_legend(
    title = paste0("Query: ", query_group),
    override.aes = list(size = 4.5, shape = 21, color = "black", fill = na.omit(color[levels(p2$data$group.by)]))
  ))
  p2legend <- get_legend(p2)
  # p2legend <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box")

  if (!is.null(p1legend) && !is.null(p2legend)) {
    legend <- cbind(p1legend, p2legend)
  } else {
    legend <- p1legend %||% p2legend
  }

  p3 <- p1 + new_scale_fill() + new_scale_color() +
    geom_point(
      data = p2data, aes(x = x, y = y), color = "black",
      size = pt.size + stroke.highlight
    ) +
    geom_point(
      data = p2data, aes(x = x, y = y, color = colour),
      size = pt.size
    ) +
    scale_color_identity() +
    facet_null() + theme(legend.position = "none")
  if (is.null(legend)) {
    return(p3)
  } else {
    gtable <- as_grob(p3)
    gtable <- add_grob(gtable, legend, "right")
    p <- wrap_plots(gtable)
    return(p)
  }
}

#' EnrichmentPlot
#'
#' @param plot_type Type of plot to be visualized.
#' @param topTerm Number of terms to plot if \code{plot_type="bar"} or \code{plot_type="lollipop"}.
#' @param topWord Number of features to plot if \code{plot_type="wordcloud"}.
#' @param pvalueCutoff pvalueCutoff
#' @param padjustCutoff padjustCutoff
#' @param palette palette
#' @param combine combine
#' @param nrow nrow
#' @param ncol ncol
#' @param byrow byrow
#' @param db
#' @param character_width
#' @param lineheight
#' @param srt
#' @param group_by
#' @param test.use
#' @param res
#' @param group_use
#' @param word_type
#' @param word_size
#' @param min_word_length
#' @param exclude_words
#' @param aspect.ratio
#' @param legend.position
#' @param legend.direction
#' @param palcolor
#' @param only_sig
#' @param theme_use
#' @param theme_args
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, db = "GO_BP", group_by = "CellType", species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", plot_type = "bar")
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"),
#'   plot_type = "bar", character_width = 30,
#'   theme_use = ggplot2::theme_classic, theme_args = list(base_size = 10)
#' )
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", only_sig = TRUE)
#'
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "bar", ncol = 1)
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "lollipop", ncol = 1)
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "wordcloud")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "wordcloud", word_type = "feature")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "comparison")
#'
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, db = c("MP", "DO"), group_by = "CellType", convert_species = TRUE, species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, db = c("MP", "DO"), group_by = "CellType", group_use = "Endocrine", ncol = 1)
#'
#' @importFrom ggplot2 ggplot geom_bar geom_text labs scale_fill_manual scale_y_continuous facet_grid coord_flip scale_color_gradientn scale_fill_gradientn scale_size guides geom_segment expansion guide_colorbar scale_color_manual guide_none
#' @importFrom dplyr group_by filter arrange desc across mutate reframe distinct n .data
#' @importFrom stats formula
#' @importFrom patchwork wrap_plots
#' @export
#'
EnrichmentPlot <- function(srt, db = "GO_BP", group_by = NULL, group_use = NULL, test.use = "wilcox",
                           res = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                           plot_type = c("bar", "lollipop", "wordcloud", "comparison"), palette = "Spectral", palcolor = NULL,
                           topTerm = 6, topWord = 100, word_type = c("term", "feature"), word_size = c(2, 8), min_word_length = 3, only_sig = FALSE,
                           exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "organization", "system", "regulation", "positive", "negative", "response", "process"),
                           aspect.ratio = 1, character_width = 50, lineheight = 0.7,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  word_type <- match.arg(word_type)

  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    slot <- paste("Enrichment", group_by, test.use, sep = "_")
    if (!slot %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunEnrichment first.")
    }
    res <- srt@tools[[slot]][["enrichment"]]
  } else {
    res <- res[["enrichment"]]
  }

  if (is.null(pvalueCutoff) && is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(res[["Database"]])) {
    res[["Database"]] <- factor(res[["Database"]], levels = unique(res[["Database"]]))
  }
  if (!is.factor(res["Groups"])) {
    res[["Groups"]] <- factor(res[["Groups"]], levels = unique(res[["Groups"]]))
  }
  if (length(db[!db %in% res[["Database"]]]) > 0) {
    stop(paste0(db[!db %in% res[["Database"]]], " is not in the enrichment result."))
  }
  if (!is.null(group_use)) {
    res <- res[res[["Groups"]] %in% group_use, , drop = FALSE]
  }

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  metric_value <- ifelse(is.null(padjustCutoff), pvalueCutoff, padjustCutoff)
  pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

  if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
    res_sim <- res[res[["Database"]] %in% gsub("_sim", "", db), , drop = FALSE]
  }
  res <- res[res[["Database"]] %in% db, , drop = FALSE]
  res_sig <- res[res[[metric]] < metric_value, , drop = FALSE]
  res_sig <- res_sig[order(res_sig[[metric]]), , drop = FALSE]
  if (nrow(res_sig) == 0) {
    stop(
      "No term enriched using the threshold: ",
      paste0("pvalueCutoff = ", pvalueCutoff), "; ",
      paste0("padjustCutoff = ", padjustCutoff)
    )
  }
  df_list <- split.data.frame(res_sig, ~ Database + Groups)
  df_list <- df_list[lapply(df_list, nrow) > 0]

  if (plot_type == "comparison") {
    ids <- NULL
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      ids <- unique(c(ids, df[head(seq_len(nrow(df)), topTerm), "ID"]))
    }
    if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
      res_sub <- subset(res_sim, ID %in% ids)
      res_sub[["Database"]] <- paste0(res_sub[["Database"]], "_sim")
    } else {
      res_sub <- subset(res, ID %in% ids)
    }
    res_sub[["GeneRatio"]] <- sapply(res_sub[["GeneRatio"]], function(x) {
      sp <- strsplit(x, "/")[[1]]
      GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
    })
    res_sub[["BgRatio"]] <- sapply(res_sub[["BgRatio"]], function(x) {
      sp <- strsplit(x, "/")[[1]]
      BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      return(BgRatio)
    })
    res_sub[["EnrichmentScore"]] <- res_sub[["GeneRatio"]] / res_sub[["BgRatio"]]
    res_sub[["Description"]] <- capitalize(res_sub[["Description"]])
    res_sub[["Description"]] <- str_wrap(res_sub[["Description"]], width = character_width)
    terms <- setNames(res_sub[["Description"]], res_sub[["ID"]])
    res_sub[["Description"]] <- factor(res_sub[["Description"]], levels = rev(terms[ids]))
    if (isTRUE(only_sig)) {
      res_sub <- res_sub[res_sub[[metric]] < metric_value, , drop = FALSE]
    }
    p <- ggplot(res_sub, aes(x = Groups, y = Description)) +
      geom_point(aes(size = GeneRatio, fill = .data[[metric]], color = ""), shape = 21) +
      scale_size_area(name = "GeneRatio", max_size = 6, n.breaks = 4) +
      guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
      scale_fill_gradientn(
        name = paste0(metric),
        limits = c(0, min(pvalueCutoff, padjustCutoff)),
        n.breaks = 3,
        colors = palette_scp(palette = palette, palcolor = palcolor, reverse = TRUE),
        na.value = "grey80",
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, order = 2)
      ) +
      scale_color_manual(values = NA, na.value = "black") +
      guides(colour = if (isTRUE(only_sig)) guide_none() else guide_legend("Non-sig", override.aes = list(colour = "black", fill = "grey80", size = 3))) +
      facet_grid(Database ~ ., scales = "free") +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        panel.grid.major = element_line(colour = "grey80", linetype = 2),
        strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, linewidth = 1),
        legend.position = legend.position,
        legend.direction = legend.direction,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(
          lineheight = lineheight, hjust = 1,
          face = ifelse(grepl("\n", levels(res_sub[["Description"]])), "italic", "plain")
        )
      )
    plist <- list(p)
  } else if (plot_type == "bar") {
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df <- df[head(seq_len(nrow(df)), topTerm), , drop = FALSE]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = rev(df[["Description"]]))

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["metric"]],
        fill = .data[["Database"]],
        label = .data[["Count"]]
      )) +
        geom_bar(width = 0.9, stat = "identity", color = "black") +
        geom_text(hjust = -0.5, size = 3.5) +
        labs(x = "", y = paste0("-log10(", metric, ")")) +
        scale_fill_manual(
          values = palette_scp(levels(df[["Database"]]), palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = "none"
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(df[["metric"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        facet_grid(Database ~ Groups, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, linewidth = 1),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.position = legend.position,
          legend.direction = legend.direction,
          axis.text.y = element_text(
            lineheight = lineheight, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        )
      return(p)
    }))
  } else if (plot_type == "lollipop") {
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df <- df[head(seq_len(nrow(df)), topTerm), , drop = FALSE]
      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df[["BgRatio"]] <- sapply(df[["BgRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
        return(BgRatio)
      })
      df[["EnrichmentScore"]] <- df[["GeneRatio"]] / df[["BgRatio"]]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = df[order(df[["EnrichmentScore"]]), "Description"])

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["EnrichmentScore"]],
        fill = .data[["metric"]]
      )) +
        geom_point(aes(size = .data[["GeneRatio"]]), shape = 21, color = "white", stroke = 1) +
        geom_segment(aes(yend = 0, xend = .data[["Description"]], color = .data[["metric"]]), linewidth = 2, lineend = "butt") +
        scale_size(name = "GeneRatio", range = c(3, 6), scales::breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_y_continuous(limits = c(0, 1.2 * max(df[["EnrichmentScore"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        labs(x = "", y = "Enrichment Score") +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = c("gold", "red3"),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1)
        ) +
        scale_color_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = c("gold", "red3"),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1)
        ) +
        facet_grid(Database ~ Groups, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          panel.background = element_rect(fill = "#1A365C"),
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, linewidth = 1),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.position = legend.position,
          legend.direction = legend.direction,
          axis.text.y = element_text(
            lineheight = lineheight, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        )
      return(p)
    }))
  } else if (plot_type == "wordcloud") {
    check_R("ggwordcloud")
    check_R("jokergoo/simplifyEnrichment")
    plist <- lapply(df_list, function(df) {
      if (word_type == "term") {
        if (df$Database[1] %in% c("GO_BP", "GO_CC", "GO_MF")) {
          df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df[["ID"]])
          if (nrow(df0 > 0)) {
            df <- df0 %>%
              reframe(
                keyword = .data[["keyword"]],
                score = -(log10(.data[["padj"]])),
                count = .data[["n_term"]],
                Database = df[["Database"]][1],
                Groups = df[["Groups"]][1]
              ) %>%
              filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
              filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
              filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
              distinct() %>%
              mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
              as.data.frame()
            df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
          } else {
            return(NULL)
          }
        } else {
          df <- df %>%
            mutate(keyword = strsplit(tolower(as.character(.data[["Description"]])), " ")) %>%
            unnest(cols = "keyword") %>%
            group_by(.data[["keyword"]], Database, Groups) %>%
            reframe(
              keyword = .data[["keyword"]],
              score = sum(-(log10(.data[[metric]]))),
              count = n(),
              Database = .data[["Database"]],
              Groups = .data[["Groups"]],
              .groups = "keep"
            ) %>%
            filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
            filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
            filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
          df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
        }
      } else {
        df <- df %>%
          mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Database, Groups) %>%
          reframe(
            keyword = .data[["keyword"]],
            score = sum(-(log10(.data[[metric]]))),
            count = n(),
            Database = .data[["Database"]],
            Groups = .data[["Groups"]],
            .groups = "keep"
          ) %>%
          distinct() %>%
          mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
          as.data.frame()
        df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
      }
      colors <- palette_scp(df[["score"]], type = "continuous", palette = palette, palcolor = palcolor, matched = FALSE)
      colors_value <- seq(min(df[["score"]], na.rm = TRUE), quantile(df[["score"]], 0.99, na.rm = TRUE) + 0.001, length.out = 100)
      p <- ggplot(df, aes(label = .data[["keyword"]], size = .data[["count"]], color = .data[["score"]], angle = .data[["angle"]])) +
        ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, eccentricity = 1, shape = "square", show.legend = TRUE, grid_margin = 3) +
        scale_color_gradientn(
          name = "Score:", colours = colors, values = rescale(colors_value),
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, title.hjust = 0)
        ) +
        scale_size(name = "Count", range = word_size, breaks = ceiling(seq(min(df[["count"]], na.rm = TRUE), max(df[["count"]], na.rm = TRUE), length.out = 3))) +
        guides(size = guide_legend(override.aes = list(colour = "black", label = "G"), order = 1)) +
        facet_grid(Database ~ Groups, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    })
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' GSEA Plot
#'
#' @param res
#' @param geneSetID
#' @param base_size
#' @param rel_heights
#' @param subplots
#' @param n_coregene
#' @param sample_coregene
#' @param features_label
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param label.size
#' @param srt
#' @param group_by
#' @param test.use
#' @param db
#' @param pvalueCutoff
#' @param padjustCutoff
#' @param topTerm
#' @param rel_width
#' @param palette
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param group_use
#' @param plot_type
#' @param palcolor
#' @param only_pos
#' @param only_sig
#' @param linewidth
#' @param line_alpha
#' @param line_color
#' @param aspect.ratio
#' @param character_width
#' @param lineheight
#' @param seed
#'
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
#' pancreas_sub <- RunGSEA(pancreas_sub, group_by = "CellType", db = "GO_BP", species = "Mus_musculus")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", geneSetID = "GO:0006412")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", geneSetID = c("GO:0046903", "GO:0015031", "GO:0007600")) %>%
#'   panel_fix_single(width = 5) # Because the plot is made by combining, we want to adjust the overall height and width
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", only_sig = TRUE)
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", pvalueCutoff = 0.05, padjustCutoff = NULL, only_pos = TRUE, only_sig = TRUE)
#' @importFrom ggplot2 ggplot aes theme theme_classic alpha element_blank element_rect margin geom_line geom_point geom_rect geom_linerange geom_hline geom_vline geom_segment annotate ggtitle labs xlab ylab scale_x_continuous scale_y_continuous scale_color_manual scale_alpha_manual guides guide_legend guide_none
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRamp
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr case_when filter pull %>%
#' @importFrom stats quantile
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid textGrob
#' @export
GSEAPlot <- function(srt, db = "GO_BP", group_by = NULL, group_use = NULL, test.use = "wilcox",
                     res = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                     plot_type = c("line", "comparison"), palette = "Spectral", palcolor = NULL,
                     topTerm = 6, geneSetID = NULL, only_pos = FALSE, only_sig = FALSE,
                     subplots = 1:3, rel_heights = c(1.5, 0.5, 1), rel_width = 3,
                     linewidth = 1.5, line_alpha = 1, line_color = "#6BB82D",
                     n_coregene = 10, sample_coregene = FALSE, features_label = NULL,
                     label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                     aspect.ratio = NULL, base_size = 12, character_width = 50, lineheight = 0.7,
                     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    slot <- paste("GSEA", group_by, test.use, sep = "_")
    if (!slot %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunGSEA first.")
    }
    enrichment <- srt@tools[[slot]][["enrichment"]]
    res <- srt@tools[[slot]][["results"]]
  } else {
    enrichment <- res[["enrichment"]]
    res <- res[["results"]]
  }
  group_use <- group_use %||% unique(enrichment[["Groups"]])
  comb <- expand.grid(group_use, db)
  use <- names(res)[names(res) %in% paste(comb$Var1, comb$Var2, sep = "-")]
  if (length(use) == 0) {
    stop(paste0(db, " is not in the enrichment result."))
  }
  res <- res[use]
  enrichment <- enrichment[enrichment[["Groups"]] %in% group_use, , drop = FALSE]

  if (is.null(pvalueCutoff) && is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(enrichment[["Database"]])) {
    enrichment[["Database"]] <- factor(enrichment[["Database"]], levels = unique(enrichment[["Database"]]))
  }
  if (!is.factor(enrichment["Groups"])) {
    enrichment[["Groups"]] <- factor(enrichment[["Groups"]], levels = unique(enrichment[["Groups"]]))
  }
  if (length(db[!db %in% enrichment[["Database"]]]) > 0) {
    stop(paste0(db[!db %in% enrichment[["Database"]]], " is not in the enrichment result."))
  }

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  metric_value <- ifelse(is.null(padjustCutoff), pvalueCutoff, padjustCutoff)
  pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

  if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
    enrichment_sim <- enrichment[enrichment[["Database"]] %in% gsub("_sim", "", db), , drop = FALSE]
  }
  enrichment <- enrichment[enrichment[["Database"]] %in% db, , drop = FALSE]

  plist <- NULL
  if (plot_type == "comparison") {
    ids <- NULL
    for (i in group_use) {
      df <- enrichment[enrichment[["Groups"]] == i, , drop = FALSE]
      df <- df[df[["pvalue"]] <= pvalueCutoff & df[["p.adjust"]] <= padjustCutoff, , drop = FALSE]
      df <- df[order(df[[metric]]), , drop = FALSE]
      df_up <- df[df[["NES"]] > 0, , drop = FALSE]
      ID_up <- df_up[head(order(df_up[[metric]]), topTerm), "ID"]
      df_down <- df[df[["NES"]] < 0, , drop = FALSE]
      ID_down <- df_down[head(order(df_down[[metric]]), topTerm), "ID"]
      if (isTRUE(only_pos)) {
        ids <- unique(c(ids, head(c(ID_up), topTerm)))
      } else {
        ids <- unique(c(ids, head(c(ID_up, ID_down), topTerm)))
      }
    }
    if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
      res_sub <- subset(enrichment_sim, ID %in% ids)
      res_sub[["Database"]] <- paste0(res_sub[["Database"]], "_sim")
    } else {
      res_sub <- subset(enrichment, ID %in% ids)
    }
    res_sub[["Description"]] <- capitalize(res_sub[["Description"]])
    res_sub[["Description"]] <- str_wrap(res_sub[["Description"]], width = character_width)
    terms <- setNames(res_sub[["Description"]], res_sub[["ID"]])
    res_sub[["Description"]] <- factor(res_sub[["Description"]], levels = rev(terms[ids]))
    res_sub[["Significant"]] <- res_sub[[metric]] < metric_value
    res_sub[["Significant"]] <- factor(res_sub[["Significant"]], levels = c("TRUE", "FALSE"))
    if (isTRUE(only_sig)) {
      res_sub <- res_sub[res_sub[["Significant"]] == "TRUE", , drop = FALSE]
    }
    if (isTRUE(only_pos)) {
      res_sub <- res_sub[res_sub[["NES"]] > 0, , drop = FALSE]
    }
    suppressWarnings(p <- ggplot(res_sub, aes(x = Groups, y = Description)) +
      geom_point(aes(size = setSize, fill = NES, color = Significant), shape = 21, stroke = 0.8) +
      scale_size_area(name = "setSize", max_size = 6, n.breaks = 4) +
      guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 2)) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
      scale_fill_gradientn(
        name = "NES",
        n.breaks = 4,
        limits = c(-max(abs(res_sub[["NES"]])), max(abs(res_sub[["NES"]]))),
        colors = palette_scp(palette = palette, palcolor = palcolor),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4, barwidth = 1, order = 1)
      ) +
      scale_color_manual(
        name = paste0("Significant\n(", metric, "<", metric_value, ")", collapse = ""), values = c("TRUE" = "black", "FALSE" = "grey90"),
        guide = if (isTRUE(only_sig)) guide_none() else guide_legend()
      ) +
      facet_grid(Database ~ ., scales = "free") +
      theme_scp(
        aspect.ratio = aspect.ratio,
        base_size = base_size,
        panel.grid.major = element_line(colour = "grey80", linetype = 2),
        strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, linewidth = 1),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(
          lineheight = lineheight, hjust = 1,
          face = ifelse(grepl("\n", levels(res_sub[["Description"]])), "italic", "plain")
        )
      ))
    plist <- list(p)
  } else if (plot_type == "line") {
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(geneSetID)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, , drop = FALSE]
        geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, , drop = FALSE]
        geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
        if (isTRUE(only_pos)) {
          geneSetID_use <- head(c(geneSetID_up), topTerm)
        } else {
          geneSetID_use <- head(c(geneSetID_up, geneSetID_down), topTerm)
        }
      } else {
        geneSetID_use <- geneSetID
      }
      if (length(geneSetID_use) == 1) {
        gsdata <- gsInfo(object = res_enrich, geneSetID = geneSetID_use)
      } else {
        gsdata <- do.call(rbind, lapply(geneSetID_use, gsInfo, object = res_enrich))
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      stat <- res_enrich[geneSetID_use, c("Description", "NES", metric)]
      rownames(stat) <- stat[, "Description"]
      stat$p.sig <- case_when(
        stat[[metric]] > 0.05 ~ "ns  ",
        stat[[metric]] <= 0.05 & stat[[metric]] > 0.01 ~ "*   ",
        stat[[metric]] <= 0.01 & stat[[metric]] > 0.001 ~ "**  ",
        stat[[metric]] <= 0.001 & stat[[metric]] > 0.0001 ~ "*** ",
        stat[[metric]] <= 0.0001 ~ "****"
      )
      gsdata$p.sig <- stat[gsdata$Description, "p.sig"]
      gsdata[["Description"]] <- capitalize(gsdata[["Description"]])
      gsdata$DescriptionP <- paste(gsdata$p.sig, gsdata$Description)
      gsdata$DescriptionP <- factor(gsdata$DescriptionP, levels = unique(gsdata$DescriptionP))
      p <- ggplot(gsdata, aes(x = x)) +
        xlab(NULL) +
        theme_classic(base_size) +
        theme(
          panel.grid.major = element_line(colour = "grey90", linetype = 2),
          panel.grid.minor = element_line(colour = "grey90", linetype = 2)
        ) +
        scale_x_continuous(expand = c(0.01, 0))
      es_layer <- geom_line(aes(y = runningScore, color = DescriptionP),
        linewidth = linewidth, alpha = line_alpha
      )
      bg_dat <- data.frame(xmin = -Inf, xmax = Inf, ymin = c(0, -Inf), ymax = c(Inf, 0), fill = c(alpha("#C40003", 0.2), alpha("#1D008F", 0.2)))
      p1 <- p +
        geom_rect(data = bg_dat, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)), inherit.aes = FALSE) +
        geom_hline(yintercept = 0, linetype = 1, color = "grey40") +
        es_layer +
        ylab("Enrichment Score") +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(color = "black", fill = "transparent", linewidth = 1),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent")
        )

      i <- 0
      for (term in rev(levels(gsdata$DescriptionP))) {
        idx <- which(gsdata$ymin != 0 & gsdata$DescriptionP ==
          term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
      }
      p2 <- ggplot(gsdata, aes(x = x)) +
        geom_linerange(aes(ymin = ymin, ymax = ymax, color = DescriptionP), alpha = line_alpha) +
        xlab(NULL) +
        ylab(NULL) +
        theme_classic(base_size) +
        theme(
          legend.position = "none",
          plot.margin = margin(t = -0.1, b = 0, r = 0.2, l = 0.2, unit = "cm"),
          panel.border = element_rect(color = "black", fill = "transparent", linewidth = 1),
          axis.line.y = element_blank(), axis.line.x = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()
        ) +
        scale_x_continuous(expand = c(0.01, 0)) +
        scale_y_continuous(expand = c(0, 0))
      if (length(geneSetID_use) == 1) {
        subtitle_use <- paste0("(NES=", round(stat[["NES"]], 3), ", ", metric, "=", format(stat[[metric]], digits = 3, scientific = TRUE), ", ", stat[["p.sig"]], ")")
        p1 <- p1 +
          annotate(
            geom = "segment", x = 0, xend = p$data$x[which.max(abs(p$data$runningScore))],
            y = p$data$runningScore[which.max(abs(p$data$runningScore))], yend = p$data$runningScore[which.max(abs(p$data$runningScore))], linetype = 2
          ) +
          annotate(
            geom = "segment", x = p$data$x[which.max(abs(p$data$runningScore))], xend = p$data$x[which.max(abs(p$data$runningScore))],
            y = 0, yend = p$data$runningScore[which.max(abs(p$data$runningScore))], linetype = 2
          ) +
          annotate(
            geom = "point", x = p$data$x[which.max(abs(p$data$runningScore))],
            y = p$data$runningScore[which.max(abs(p$data$runningScore))],
            fill = ifelse(stat[["NES"]] < 0, "#5E34F5", "#F52323"), color = "black", size = 2.5,
            shape = ifelse(stat[["NES"]] < 0, 25, 24)
          ) +
          labs(subtitle = subtitle_use) +
          theme(plot.subtitle = element_text(face = "italic"))

        if ((is.numeric(n_coregene) && n_coregene > 1) || length(features_label) > 0) {
          if (length(features_label) == 0) {
            features_label <- unlist(strsplit(gsdata$CoreGene[1], "/"))
            n_coregene <- min(n_coregene, length(features_label))
            if (isTRUE(sample_coregene)) {
              features_label <- sample(features_label, n_coregene, replace = FALSE)
            } else {
              features_label <- gsdata$GeneName[gsdata$GeneName %in% features_label][1:n_coregene]
            }
          }
          df_gene <- gsdata[gsdata$position == 1 & gsdata$GeneName %in% features_label, , drop = FALSE]
          gene_drop <- features_label[!features_label %in% df_gene$GeneName]
          if (length(gene_drop) > 0) {
            warning("Gene ", paste(gene_drop, collapse = ","), " is not in the geneset of the ", gsdata$Description[1], immediate. = TRUE)
          }
          x_nudge <- diff(range(gsdata$x)) * 0.05
          y_nudge <- diff(range(gsdata$runningScore)) * 0.05
          p1 <- p1 + geom_point(
            data = df_gene,
            mapping = aes(y = runningScore), color = "black"
          ) +
            geom_text_repel(
              data = df_gene,
              mapping = aes(y = runningScore, label = GeneName),
              min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40",
              color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size,
              nudge_x = ifelse(df_gene$runningScore >= 0, x_nudge, -x_nudge),
              nudge_y = ifelse(df_gene$runningScore > 0, -y_nudge, y_nudge)
            )
        }

        x <- p$data$x
        y <- y_raw <- p$data$geneList
        y[y > quantile(y_raw, 0.98)] <- quantile(y_raw, 0.98)
        y[y < quantile(y_raw, 0.02)] <- quantile(y_raw, 0.02)
        col <- rep("white", length(y))
        y_pos <- which(y > 0)
        if (length(y_pos) > 0) {
          y_pos_i <- cut(y[y_pos],
            breaks = seq(min(y[y_pos], na.rm = TRUE), max(y[y_pos], na.rm = TRUE), len = 100),
            include.lowest = TRUE
          )
          col[y_pos] <- colorRampPalette(c("#F5DCDC", "#C40003"))(100)[y_pos_i]
        }

        y_neg <- which(y < 0)
        if (length(y_neg) > 0) {
          y_neg_i <- cut(y[y_neg],
            breaks = seq(min(y[y_neg], na.rm = TRUE), max(y[y_neg], na.rm = TRUE), len = 100),
            include.lowest = TRUE
          )
          col[y_neg] <- colorRampPalette(c("#1D008F", "#DDDCF5"))(100)[y_neg_i]
        }

        ymin <- min(p2$data$ymin, na.rm = TRUE)
        ymax <- max(p2$data$ymax - p2$data$ymin, na.rm = TRUE) * 0.3
        xmin <- which(!duplicated(col))
        xmax <- xmin + as.numeric(table(col)[as.character(unique(col))])
        d <- data.frame(
          ymin = ymin, ymax = ymax, xmin = xmin,
          xmax = xmax, col = unique(col)
        )
        p2 <- p2 + geom_rect(
          aes(
            xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax, fill = I(col)
          ),
          data = d,
          alpha = 0.95, inherit.aes = FALSE
        )
      }
      df2 <- p$data
      df2$y <- p$data$geneList[df2$x]
      min_y <- df2$y[which.min(abs(df2$y))]
      corss_x <- median(df2$x[df2$y == min_y])
      p3 <- p + geom_segment(data = df2, aes(
        x = x, xend = x,
        y = y, yend = 0
      ), color = "grey30")

      if (max(df2$y) > 0) {
        p3 <- p3 + annotate(geom = "text", x = 0, y = Inf, vjust = 1.3, hjust = 0, color = "#C81A1F", size = 4, label = " Positively correlated")
      }
      if (min(df2$y) < 0) {
        p3 <- p3 + annotate(geom = "text", x = Inf, y = -Inf, vjust = -0.3, hjust = 1, color = "#3C298C", size = 4, label = "Negtively correlated ")
      }
      if (max(df2$y) > 0 && min(df2$y) < 0) {
        p3 <- p3 + geom_vline(xintercept = corss_x, linetype = 2, color = "black") +
          annotate(geom = "text", y = 0, x = corss_x, vjust = ifelse(diff(abs(range(df2$y))) > 0, -0.3, 1.3), size = 4, label = paste0("Zero cross at ", corss_x))
      }
      p3 <- p3 + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") +
        theme(
          plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"),
          axis.line = element_blank(), axis.line.x = element_blank(),
          panel.border = element_rect(color = "black", fill = "transparent", linewidth = 1)
        )
      if (length(geneSetID_use) == 1) {
        p1 <- p1 + ggtitle(gsdata$Description[1], subtitle = subtitle_use)
      }
      if (length(line_color) != length(geneSetID_use)) {
        color_use <- palette_scp(levels(gsdata$DescriptionP), palette = palette, palcolor = palcolor)
      } else {
        color_use <- line_color
      }
      p1 <- p1 + scale_color_manual(values = color_use)
      if (length(color_use) == 1) {
        p1 <- p1 + theme(legend.position = "none")
        p2 <- p2 + scale_color_manual(values = "black")
      } else {
        p2 <- p2 + scale_color_manual(values = color_use)
      }
      legend <- get_legend(
        p1 +
          guides(color = guide_legend(
            ncol = 2,
            label.theme = element_text(size = base_size),
            override.aes = list(size = 3)
          )) +
          theme(
            legend.position = "bottom",
            legend.spacing = unit(0, "cm"),
            legend.background = element_blank(),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0)
          )
      )
      plotlist <- list(p1 + theme(legend.position = "none"), p2, p3)[subplots]
      if (length(subplots) == 1) {
        plist[[nm]] <- plotlist[[1]] + theme(
          aspect.ratio = rel_heights[subplots] / rel_width,
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
        )
      } else {
        plotlist <- lapply(plotlist[subplots], as_grob)
        rel_heights <- rel_heights[subplots]
        for (i in seq_along(plotlist)) {
          plotlist[[i]] <- panel_fix_single(plotlist[[i]], height = rel_heights[i], units = "null", margin = 0, respect = TRUE, return_grob = TRUE)
          plotlist[[i]] <- panel_fix_single(plotlist[[i]], width = rel_width, units = "null", margin = 0, respect = TRUE, return_grob = TRUE)
        }
        p_out <- do.call(rbind, c(plotlist, size = "first"))

        if (length(geneSetID_use) > 1) {
          p_out <- add_grob(p_out, legend, "top")
        }
        lab <- textGrob(label = nm, rot = -90, hjust = 0.5)
        p_out <- add_grob(p_out, lab, "right", clip = "off")
        p_out <- wrap_plots(p_out)
        plist[[nm]] <- p_out
      }
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  df$CoreGene <- object@result[geneSetID, "core_enrichment"]
  if (length(object@gene2Symbol) == length(object@geneList)) {
    df$GeneName <- object@gene2Symbol
  } else {
    df$GeneName <- df$gene
  }
  return(df)
}

gseaScores <- function(geneList, geneSet, exponent = 1) {
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit / NR)
  Pmiss[!hits] <- 1 / (N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES, na.rm = TRUE)
  min.ES <- min(runningES, na.rm = TRUE)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  df <- data.frame(
    x = seq_along(runningES), runningScore = runningES,
    position = as.integer(hits), gene = names(geneList)
  )
  return(df)
}

#' @importFrom ggplot2 ggplot_build ggplot_gtable panel_rows panel_cols wrap_dims
#' @importFrom gtable gtable
#' @importFrom grid unit unit.pmax is.unit
#' @importFrom utils modifyList
#' @importFrom stats na.omit
#' @importFrom BiocParallel bplapply
build_patchwork <- function(x, guides = "auto", BPPARAM = BiocParallel::SerialParam()) {
  x$layout <- modifyList(patchwork:::default_layout, x$layout[!vapply(x$layout, is.null, logical(1))])

  guides <- if (guides == "collect" && x$layout$guides != "keep") {
    "collect"
  } else {
    x$layout$guides
  }
  # bpprogressbar(BPPARAM) <- TRUE
  gt <- bplapply(x$plots, patchwork:::plot_table, guides = guides, BPPARAM = BPPARAM)
  fixed_asp <- vapply(gt, function(x) isTRUE(x$respect), logical(1))
  guide_grobs <- unlist(lapply(gt, `[[`, "collected_guides"), recursive = FALSE)
  gt <- bplapply(gt, patchwork:::simplify_gt, BPPARAM = BPPARAM)
  gt <- patchwork:::add_insets(gt)
  if (is.null(x$layout$design)) {
    if (is.null(x$layout$ncol) && !is.null(x$layout$widths) && length(x$layout$widths) > 1) {
      x$layout$ncol <- length(x$layout$widths)
    }
    if (is.null(x$layout$nrow) && !is.null(x$layout$heights) && length(x$layout$heights) > 1) {
      x$layout$nrow <- length(x$layout$heights)
    }
    dims <- wrap_dims(length(gt), nrow = x$layout$nrow, ncol = x$layout$ncol)
    x$layout$design <- patchwork:::create_design(dims[2], dims[1], x$layout$byrow)
  } else {
    dims <- c(
      max(x$layout$design$b),
      max(x$layout$design$r)
    )
  }

  TABLE_COLS <- patchwork:::TABLE_COLS
  TABLE_ROWS <- patchwork:::TABLE_ROWS
  PANEL_ROW <- patchwork:::PANEL_ROW
  PANEL_COL <- patchwork:::PANEL_COL

  gt_new <- gtable(
    unit(rep(0, TABLE_COLS * dims[2]), "null"),
    unit(rep(0, TABLE_ROWS * dims[1]), "null")
  )
  design <- as.data.frame(unclass(x$layout$design))
  if (nrow(design) < length(gt)) {
    warning("Too few patch areas to hold all plots. Dropping plots", call. = FALSE)
    gt <- gt[seq_len(nrow(design))]
    fixed_asp <- fixed_asp[seq_len(nrow(design))]
  } else {
    design <- design[seq_along(gt), ]
  }
  if (any(design$t < 1)) design$t[design$t < 1] <- 1
  if (any(design$l < 1)) design$l[design$l < 1] <- 1
  if (any(design$b > dims[1])) design$b[design$b > dims[1]] <- dims[1]
  if (any(design$r > dims[2])) design$r[design$r > dims[2]] <- dims[2]
  max_z <- lapply(gt, function(x) max(x$layout$z))
  max_z <- c(0, cumsum(max_z))
  gt_new$layout <- do.call(rbind, lapply(seq_along(gt), function(i) {
    loc <- design[i, ]
    lay <- gt[[i]]$layout
    lay$name <- paste0(lay$name, "-", i)
    lay$t <- lay$t + ifelse(lay$t <= PANEL_ROW, (loc$t - 1) * TABLE_ROWS, (loc$b - 1) * TABLE_ROWS)
    lay$l <- lay$l + ifelse(lay$l <= PANEL_COL, (loc$l - 1) * TABLE_COLS, (loc$r - 1) * TABLE_COLS)
    lay$b <- lay$b + ifelse(lay$b < PANEL_ROW, (loc$t - 1) * TABLE_ROWS, (loc$b - 1) * TABLE_ROWS)
    lay$r <- lay$r + ifelse(lay$r < PANEL_COL, (loc$l - 1) * TABLE_COLS, (loc$r - 1) * TABLE_COLS)
    lay$z <- lay$z + max_z[i]
    lay
  }))
  table_dimensions <- patchwork:::table_dims(
    lapply(gt, `[[`, "widths"),
    lapply(gt, `[[`, "heights"),
    design,
    dims[2],
    dims[1]
  )
  gt_new$grobs <- patchwork:::set_grob_sizes(gt, table_dimensions$widths, table_dimensions$heights, design)
  gt_new$widths <- table_dimensions$widths
  gt_new$heights <- table_dimensions$heights
  widths <- rep(x$layout$widths, length.out = dims[2])
  heights <- rep(x$layout$heights, length.out = dims[1])
  gt_new <- patchwork:::set_panel_dimensions(gt_new, gt, widths, heights, fixed_asp, design)
  if (x$layout$guides == "collect") {
    guide_grobs <- patchwork:::collapse_guides(guide_grobs)
    if (length(guide_grobs) != 0) {
      theme <- x$annotation$theme
      if (!attr(theme, "complete")) {
        theme <- theme_get() + theme
      }
      guide_grobs <- patchwork:::assemble_guides(guide_grobs, theme)
      gt_new <- patchwork:::attach_guides(gt_new, guide_grobs, theme)
    }
  } else {
    gt_new$collected_guides <- guide_grobs
  }

  class(gt_new) <- c("gtable_patchwork", class(gt_new))
  gt_new
}

#' @importFrom utils modifyList
patchworkGrob <- function(x, BPPARAM = BiocParallel::SerialParam(), ...) {
  annotation <- modifyList(patchwork:::default_annotation, x$patches$annotation[!vapply(x$patches$annotation, is.null, logical(1))])
  x <- patchwork:::recurse_tags(x, annotation$tag_levels, annotation$tag_prefix, annotation$tag_suffix, annotation$tag_sep)$patches
  plot <- patchwork:::get_patches(x)
  gtable <- build_patchwork(plot, BPPARAM = BPPARAM)
  gtable <- patchwork:::annotate_table(gtable, annotation)
  class(gtable) <- setdiff(class(gtable), "gtable_patchwork")
  gtable
}

#' @importFrom gridGraphics echoGrob
#' @importFrom grid grobTree
#' @importFrom ggplot2 ggplotGrob
as_grob <- function(plot, ...) {
  if (inherits(plot, "recordedplot")) {
    gridGraphics::echoGrob(plot)
  } else if (inherits(plot, "gList")) {
    grobTree(plot)
  } else if (inherits(plot, "patchwork")) {
    patchworkGrob(plot, ...)
  } else if (inherits(plot, "ggplot")) {
    ggplotGrob(plot)
  } else {
    warning("Cannot convert object of class ", paste0(class(plot), collapse = ","), " into a grob.")
  }
}

#' @importFrom grid unit
#' @importFrom gtable gtable_col
as_gtable <- function(plot, ...) {
  if (inherits(plot, "gtable")) {
    return(plot)
  }
  if (inherits(plot, "grob")) {
    u <- unit(1, "null")
    gt <- gtable_col(NULL, list(plot), u, u)
    gt$layout$clip <- "inherit"
    return(gt)
  } else {
    grob <- as_grob(plot, ...)
    if (inherits(grob, "gtable")) {
      return(grob)
    } else {
      return(as_gtable(grob, ...))
    }
  }
}

get_legend <- function(plot) {
  plot <- as_gtable(plot)
  grob_names <- plot$layout$name
  grobs <- plot$grobs
  grobIndex <- which(grepl("guide-box", grob_names))
  grobIndex <- grobIndex[1]
  matched_grobs <- grobs[[grobIndex]]
  return(matched_grobs)
}

#' @importFrom grid is.grob grobWidth grobHeight
#' @importFrom gtable is.gtable gtable_add_rows gtable_add_cols gtable_add_grob
add_grob <- function(gtable, grob, position = c("top", "bottom", "left", "right", "none"), space = NULL, clip = "on") {
  position <- match.arg(position)
  if (position == "none") {
    return(gtable)
  }

  if (is.null(space)) {
    if (is.gtable(grob)) {
      if (position %in% c("top", "bottom")) {
        space <- sum(grob$heights)
      } else {
        space <- sum(grob$widths)
      }
    } else if (is.grob(grob)) {
      if (position %in% c("top", "bottom")) {
        space <- grobHeight(grob)
      } else {
        space <- grobWidth(grob)
      }
    }
  }

  if (position == "top") {
    gtable <- gtable_add_rows(gtable, space, 0)
    gtable <- gtable_add_grob(gtable, grob, t = 1, l = mean(gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]), clip = clip)
  }
  if (position == "bottom") {
    gtable <- gtable_add_rows(gtable, space, -1)
    gtable <- gtable_add_grob(gtable, grob, t = dim(gtable)[1], l = mean(gtable$layout[grepl(pattern = "panel", x = gtable$layout$name), "l"]), clip = clip)
  }
  if (position == "left") {
    gtable <- gtable_add_cols(gtable, space, 0)
    gtable <- gtable_add_grob(gtable, grob, t = mean(gtable$layout[grep("panel", gtable$layout$name), "t"]), l = 1, clip = clip)
  }
  if (position == "right") {
    gtable <- gtable_add_cols(gtable, space, -1)
    gtable <- gtable_add_grob(gtable, grob, t = mean(gtable$layout[grep("panel", gtable$layout$name), "t"]), l = dim(gtable)[2], clip = clip)
  }
  return(gtable)
}
