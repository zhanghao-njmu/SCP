#' SCP theme
#'
#' The default theme for SCP plot function.
#'
#' @param aspect.ratio Aspect ratio of the panel.
#' @param ... Arguments passed to the \code{\link[ggplot2]{theme}}.
#' @param base_size
#'
#' @importFrom ggplot2 theme theme_classic element_blank element_text element_rect margin
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
    panel.border = element_rect(fill = "transparent", colour = "black", size = 1),
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
  return(theme_classic() + out)
}

#' Blank theme
#'
#' @param add_coord
#'
#' @param xlab
#' @param ylab
#' @param xlen_npc
#' @param ylen_npc
#' @param lab_cex
#' @param ...
#'
#' @importFrom ggplot2 theme element_blank margin annotation_custom coord_cartesian coord_sf
#' @importFrom grid grobTree gList linesGrob textGrob arrow gpar
#' @export
theme_blank <- function(add_coord = TRUE, xlab = "Dim1", ylab = "Dim2", xlen_npc = 0.15, ylen_npc = 0.15, lab_cex = 1, ...) {
  if (length(xlab) == 0) {
    xlab <- "Dim1"
  }
  if (length(ylab) == 0) {
    ylab <- "Dim2"
  }
  if (length(xlen_npc) == 0) {
    xlen_npc <- 0.15
  }
  if (length(ylen_npc) == 0) {
    ylen_npc <- 0.15
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
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.tag = element_blank(),
    plot.margin = margin(0, 0, 15, 15),
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
      textGrob(label = xlab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = 1.5, hjust = 0, gp = gpar(cex = lab_cex)),
      linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = ylab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = -0.75, hjust = 0, rot = 90, gp = gpar(cex = lab_cex))
    ))
    return(list(
      list(annotation_custom(g)),
      list(theme_scp() + out),
      list(coord_sf(clip = "off"))
    ))
  } else {
    return(list(
      list(theme_scp() + out)
    ))
  }
}

#' SCP palette
#'
#' The default palette for SCP plot function.
#'
#' @param x A vector of character/factor or numeric values for color mapping.
#' @param palette Palette name. All palette names can be queried with \code{names(SCP:::palette_list)}.
#' @param type Type of \code{x}.Can be one of "auto","discrete" or "continuous". The default is "auto", which automatically detects if \code{x} is a numeric value.
#' If it is numeric, then \code{x} is interpolated to the palette color of length 100 in ascending order.
#' Otherwise, the color is assigned to \code{x} according to its type number.
#' @param matched Whether to return a matched color vector of length \code{x}
#' @param reverse Whether to invert the colors.
#' @param NA_keep Whether to keep the color assignment to NA in \code{x}.
#' @param NA_color NA colors if NA_keep is \code{TRUE}.
#' @param n
#' @param palcolor
#'
#' @examples
#' x <- c(1:3, NA, 3:5)
#' pal1 <- palette_scp(x, palette = "Spectral")
#' pal2 <- palette_scp(x, palette = "Spectral", matched = TRUE)
#' pal3 <- palette_scp(x, palette = "Spectral", matched = TRUE, NA_keep = TRUE)
#' pal4 <- palette_scp(x, palette = "Paired", type = "discrete")
#' list_palette(list(pal1, pal2, pal3, pal4))
#'
#' names(SCP:::palette_list)
#' list_palette(SCP:::palette_list)
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
#'   ocean_names <- names(pals:::syspals)[str_detect(names(pals:::syspals), "ocean")]
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
    stop(paste0("Invalid palette Must be one of ", paste0(names(palette_list), collapse = ",")))
  }
  if (is.null(palcolor)) {
    palcolor <- palette_list[[palette]]
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
      values <- cut(x, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
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

#' Set the panel width/height of a plot object to a fixed value.
#'
#' The ggplot object, when stored, can only specify the height and width of the plot, not the panel.
#' The latter is obviously more important to control the final result of a plot.
#' This function can set the panel width/height of plot to a fixed value.
#' The plot to be fixed can be a ggplot or a plot combined using plot_grid or patchwork functions.
#'
#' @param x A ggplot or grob object.
#' @param width The width of the panel.
#' @param height The height of the panel.
#' @param margin Margins around the plot when \code{file} is specified.
#' @param units The units in which \code{height}, \code{width} and \code{margin} are given. Can be \code{mm}, \code{cm}, \code{inches}, etc. See \code{\link[grid]{unit}}.
#' @param save NULL or the file name to save.
#' @param space
#' @param bg_color
#' @param raster
#' @param dpi
#' @param return_grob
#' @param verbose
#' @param ...
#'
#' @return If \code{filename} is not specified, no return; otherwise return a list objects:
#' \itemize{
#' \item{\code{grob} An R object of class "grob", a graphical object.}
#' \item{\code{plot_size} Plot size after fixing the panel size and the margin.}
#' \item{\code{units} The units used.}
#' }
#'
#' @examples
#' library(ggplot2)
#' p <- qplot(mpg, wt, data = mtcars, colour = cyl) + facet_wrap(~gear, nrow = 2)
#' p_fix1 <- panel_fix(p, width = 5, height = 3, units = "cm")
#' p_fix2 <- panel_fix(p, width = 5, height = 3, units = "cm", raster = TRUE, dpi = 72)
#' p_fix1
#' p_fix2
#'
#' ## Save the plot with appropriate size
#' # plot_size <- attr(p_fix1, "size")
#' # ggsave(
#' #   filename = "p_fix.png", plot = p_fix1,
#' #   units = plot_size$units, width = plot_size$width, height = plot_size$height
#' # )
#'
#' ## or save the plot directly
#' # p_fix1 <- panel_fix(p, width = 5, height = 3, units = "cm", save = "p_fix.png")
#'
#' data("pancreas_sub")
#' p1 <- ClassDimPlot(pancreas_sub, "Phase") # ggplot object
#' p2 <- ExpDimPlot(pancreas_sub, "Ins1") # ggplot object
#' p <- p1 | p2 # plot is generated by patchwork
#' panel_fix(p, height = 1)
#'
#' library(cowplot)
#' p1 <- ClassDimPlot(pancreas_sub, c("Phase", "SubCellType")) # plot is generated by plot_grid
#' p2 <- ExpDimPlot(pancreas_sub, "Ins1") | ExpDimPlot(pancreas_sub, "Ins2") # plot is generated by patchwork
#' p <- plot_grid(p1, p2, nrow = 2) # plot is generated by plot_grid
#' panel_fix(p, height = 1)
#' @importFrom ggplot2 ggsave
#' @importFrom gtable gtable_add_padding
#' @importFrom grid grob unit convertWidth convertHeight
#' @importFrom cowplot plot_grid as_gtable
#' @export
#'
panel_fix <- function(x = NULL,
                      width = NULL, height = NULL, margin = 0.5, space = 0.1, units = "in", bg_color = "white",
                      raster = FALSE, dpi = 300,
                      return_grob = FALSE, save = NULL, verbose = TRUE, ...) {
  if (!"gtable" %in% class(x)) {
    tryCatch(
      {
        grob <- as_gtable(x)
      },
      error = function(error) {
        stop(error, "\nCannot convert the x to a gtable object.")
      }
    )
  } else {
    grob <- x
  }
  args <- as.list(match.call())[-1]
  depth <- args[["depth"]]
  if (is.null(depth)) {
    depth <- 1
  }

  panel_index <- grep("panel", grob[["layout"]][["name"]])
  if (length(panel_index) > 1) {
    message("More than 2 panels detected. panel_fix may not work as expected.")
  }
  for (i in panel_index) {
    geom_index <- grep("GeomDrawGrob", names(grob$grobs[[i]][["children"]]))
    if (length(geom_index) > 0) {
      if (isTRUE(verbose)) {
        message("panel ", i, " is detected as generated by plot_grid.")
      }
      for (j in geom_index) {
        subgrob <- grob$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]]
        subgrob <- panel_fix(subgrob, width = width, height = height, margin = space, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
        grob$grobs[[i]][["children"]][[j]][["children"]][[1]][["children"]][[1]] <- subgrob
        # print(paste0("plot_width:",plot_width," plot_height:",plot_height))
      }
      panel_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE) / as.numeric(grob$grobs[[i]][["children"]][[j]]$vp$width)
      panel_height <- convertWidth(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE) / as.numeric(grob$grobs[[i]][["children"]][[j]]$vp$height)
      grob <- panel_fix_single(grob, panel_index = i, width = panel_width, height = panel_height, margin = ifelse(depth == 1, margin, 0), units = units, raster = FALSE)
    } else if (inherits(grob$grobs[[i]], "gtable")) {
      if (isTRUE(verbose)) {
        message("panel ", i, " is detected as generated by patchwork.")
      }
      subgrob <- grob$grobs[[i]]
      subgrob <- panel_fix(subgrob, width = width, height = height, margin = space, units = units, raster = raster, dpi = dpi, return_grob = TRUE, verbose = verbose, depth = depth + 1)
      grob$grobs[[i]] <- subgrob
      panel_width <- convertWidth(sum(subgrob[["widths"]]), unitTo = units, valueOnly = TRUE)
      panel_height <- convertWidth(sum(subgrob[["heights"]]), unitTo = units, valueOnly = TRUE)
      grob <- panel_fix_single(grob, panel_index = i, width = panel_width, height = panel_height, margin = ifelse(depth == 1, margin, 0), units = units, raster = FALSE, respect = TRUE)
    } else {
      grob <- panel_fix_single(grob, panel_index = i, width = width, height = height, margin = margin, units = units, raster = raster, dpi = dpi)
    }
  }
  plot_width <- convertWidth(sum(grob[["widths"]]), unitTo = units, valueOnly = TRUE)
  plot_height <- convertHeight(sum(grob[["heights"]]), unitTo = units, valueOnly = TRUE)
  p <- plot_grid(grob) + theme(plot.background = element_rect(fill = bg_color, color = bg_color))
  attr(p, "size") <- list(width = plot_width, height = plot_height, units = units)
  if (!is.null(save) && is.character(save) && nchar(save) > 0) {
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
  if (isTRUE(return_grob)) {
    return(grob)
  } else {
    return(p)
  }
}

#' Set the panel width/height for the single plot.
#'
#' @param x
#'
#' @param panel_index
#' @param width
#' @param height
#' @param margin
#' @param units
#' @param raster
#' @param dpi
#'
#' @importFrom ggplot2 ggsave zeroGrob
#' @importFrom gtable gtable_add_padding
#' @importFrom grid grob unit unitType convertWidth convertHeight viewport grid.draw rasterGrob grobTree
#' @importFrom cowplot plot_grid as_gtable
#' @export
panel_fix_single <- function(x, panel_index = NULL, respect = NULL,
                             width = NULL, height = NULL, margin = 0.5, units = "in",
                             raster = FALSE, dpi = 300) {
  if (!inherits(x, "gtable")) {
    tryCatch(
      {
        grob <- as_gtable(x)
      },
      error = function(error) {
        stop(error, "\nCannot convert the x to a gtable object")
      }
    )
  } else {
    grob <- x
  }

  if (is.null(panel_index)) {
    panel_index <- grep("panel", grob[["layout"]][["name"]])
  }
  if (length(panel_index) == 0) {
    warning("No panel detected.")
    return(x)
  }

  panel_index_h <- sort(unique(grob[["layout"]][["t"]][panel_index]))
  panel_index_w <- sort(unique(grob[["layout"]][["l"]][panel_index]))
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  w <- grob[["widths"]][panel_index_w]
  h <- grob[["heights"]][panel_index_h]

  if (units != "null") {
    raw_w <- as.numeric(convertWidth(w, unitTo = units))
    raw_h <- as.numeric(convertHeight(h, unitTo = units))
    if (raw_w < 1e10) {
      raw_w <- 0
    }
    if (raw_h < 1e10) {
      raw_h <- 0
    }
    if (unitType(w) == unitType(h)) {
      raw_aspect <- as.vector(h) / as.vector(w)
    } else {
      if (raw_w != 0) {
        raw_aspect <- raw_h / raw_w
      } else {
        raw_aspect <- 1
      }
    }
    if (is.null(width) && is.null(height)) {
      width <- raw_w
      height <- raw_h
      if (width == 0 && height == 0) {
        stop("Unable to calculate size from the plot. 'width' or 'height' must be provided with at least one.")
      }
    }
    for (i in seq_along(raw_aspect)) {
      if (is.finite(raw_aspect[i]) && raw_aspect[i] != 0) {
        if (is.null(width) || is.na(width)) {
          width <- height / raw_aspect[i]
        }
        if (is.null(height) || is.na(height)) {
          height <- width * raw_aspect[i]
        }
      }
    }
  }
  # print(paste0("width:",width))
  # print(paste0("height:",height))

  if (!length(width) %in% c(0, 1, length(panel_index)) || !length(height) %in% c(0, 1, length(panel_index))) {
    stop("The length of 'width' and 'height' must be 1 or the length of panels.")
  }
  if (length(width) == 1) {
    width <- rep(width, nw)
  }
  if (length(height) == 1) {
    height <- rep(height, nh)
  }
  if (!is.null(width)) {
    grob[["widths"]][panel_index_w] <- unit(width, units = units)
  }
  if (!is.null(height)) {
    grob[["heights"]][panel_index_h] <- unit(height, units = units)
  }
  grob <- gtable_add_padding(grob, unit(margin, units = units))

  if (isTRUE(raster)) {
    check_R(c("png", "grDevices"))
    for (i in panel_index) {
      g <- grob$grobs[[i]]
      child_list <- list()
      for (j in seq_along(g[["children"]])) {
        child <- g[["children"]][[j]]
        if (!is.null(child$vp) || isTRUE(grepl("(text)|(label)", names(g[["children"]])[j]))) {
          g[["children"]][[j]] <- zeroGrob()
          child_list <- c(child_list, list(child))
        }
      }
      g$vp <- viewport()
      temp <- tempfile(fileext = "png")
      grDevices::png(temp, width = width[1], height = height[1], bg = "transparent", res = dpi, units = units)
      grid.draw(g)
      grDevices::dev.off()
      g_ras <- rasterGrob(png::readPNG(temp, native = TRUE))
      unlink(temp)
      grob$grobs[[i]] <- do.call(grobTree, c(list(g_ras), child_list))
      # grid.draw(grob)
    }
  }

  if (!is.null(respect)) {
    grob$respect <- respect
  }

  return(grob)
}

#' Convert a color with arbitrary transparency to a fixed color
#' @param colors
#'
#' @param alpha
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
#' @param colors
#'
#' @param mode
#'
#' @export
blendcolors <- function(colors, mode = "blend") {
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
  if (mode == "mix") {
    Result <- (c1 + c2) / 2
    Result[Result > 1] <- 1
  }
  if (mode == "blend") {
    Result <- (c1 * c1a + c2 * c2a * (1 - c1a)) / A
    A <- 1
  }
  if (mode == "screen") {
    Result <- 1 - (1 - c1) * (1 - c2)
  }
  if (mode == "multiply") {
    Result <- c1 * c2
  }

  return(list(Result, A))
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


#' Find the default reduction name in the srt object.
#'
#' @param srt A Seurat object.
#' @param pattern Name of the pattern to search for.
#' @param min_dim Minimum dimension threshold.
#' @param max_distance
#'
#' @return Default reduction name.
#'
#' @export
DefaultReduction <- function(srt, pattern = NULL, min_dim = 2, max_distance = 0.1) {
  if (length(Reductions(srt)) == 0) {
    stop("Unable to find any reductions.")
  }
  pattern_default <- c("umap", "tsne", "dm", "pca", "ica", "nmf", "mds", "glmpca")
  pattern_dim <- c("2D", "3D")
  reduc_all <- names(srt@reductions)
  reduc_all <- reduc_all[unlist(lapply(reduc_all, function(x) {
    dim(srt@reductions[[x]]@cell.embeddings)[2] >= min_dim
  }))]
  if (length(reduc_all) == 0) {
    stop("No dimensional reduction found in the srt object.")
  }
  if (length(reduc_all) == 1) {
    return(reduc_all)
  }
  if (is.null(pattern)) {
    if (("Default_reduction" %in% names(srt@misc))) {
      pattern <- srt@misc[["Default_reduction"]]
    } else {
      pattern <- pattern_default
    }
  }

  pattern <- c(pattern, paste0(pattern, min_dim, "D"))
  if (any(pattern %in% reduc_all)) {
    return(pattern[pattern %in% reduc_all][1])
  }
  index <- c(unlist(sapply(pattern, function(pat) {
    grep(pattern = pat, x = reduc_all, ignore.case = TRUE)
  })))
  if (length(index) > 0) {
    default_reduc <- reduc_all[index]
  } else {
    index <- c(unlist(sapply(pattern, function(pat) {
      agrep(pattern = pat, x = reduc_all, max.distance = max_distance, ignore.case = TRUE)
    })))
    if (length(index) > 0) {
      default_reduc <- reduc_all[index]
    } else {
      default_reduc <- reduc_all
    }
  }
  if (length(default_reduc) > 1) {
    default_reduc <- default_reduc[unlist(sapply(c(pattern_default, pattern_dim), function(pat) {
      grep(pattern = pat, x = default_reduc, ignore.case = TRUE)
    }))]
    default_reduc <- default_reduc[which.min(sapply(default_reduc, function(x) dim(srt[[x]])[2]))]
  }
  return(default_reduc)
}

#' 2D-Dimensional reduction plot for classification.
#'
#' Plotting cell points on a reduced 2D plane and coloring according to the groups of the cells.
#'
#' @param srt A Seurat object.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' @param reduction Which dimensionality reduction to use.
#' @param split.by Name of a metadata column to split plot by.
#' @param palette Name of palette to use. Default is "Paired".
#' @param bg_color Color value for NA points.
#' @param pt.size Point size for plotting.
#' @param pt.alpha Point transparency.
#' @param label Whether to label the groups.
#' @param label_insitu Whether place the raw label text in the original position. Default is FALSE, which using numbers instead of raw labels.
#' @param label.size Plot resolution. Only valid when \code{units} is not \code{px}
#' @param cells.highlight A vector of names of cells to highlight.
#' @param cols.highlight A vector of colors to highlight the cells.
#' @param sizes.highlight Size of highlighted cells.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical")
#' @param combine Whether to arrange multiple plots into a grid.
#' @param nrow Number of rows in the plot grid.
#' @param ncol Number of columns in the plot grid.
#' @param byrow Logical value indicating if the plots should be arrange by row (default) or by column.
#' @param dims
#' @param show_na
#' @param show_stat
#' @param palcolor
#' @param label.fg
#' @param label.bg
#' @param label.bg.r
#' @param label_repel
#' @param label_repulsion
#' @param label_point_size
#' @param label_point_color
#' @param label_segment_color
#' @param alpha.highlight
#' @param stroke.highlight
#' @param add_density
#' @param density_color
#' @param density_filled
#' @param density_filled_palette
#' @param density_filled_color
#' @param lineages
#' @param lineages_trim
#' @param lineages_span
#' @param lineages_palette
#' @param lineages_palcolor
#' @param lineages_arrow
#' @param lineages_line_size
#' @param lineages_line_bg
#' @param lineages_line_bg_r
#' @param lineages_whiskers
#' @param lineages_whiskers_size
#' @param lineages_whiskers_alpha
#' @param stat.by
#' @param stat_type
#' @param stat_plot_type
#' @param stat_plot_size
#' @param stat_plot_palette
#' @param stat_plot_alpha
#' @param stat_plot_label
#' @param stat_plot_label_size
#' @param graph
#' @param edge_size
#' @param edge_alpha
#' @param edge_color
#' @param paga
#' @param paga_node_palette
#' @param paga_node_size
#' @param paga_edge_threshold
#' @param paga_edge_size
#' @param paga_edge_color
#' @param paga_edge_alpha
#' @param paga_transition_threshold
#' @param paga_transition_size
#' @param paga_transition_color
#' @param paga_transition_alpha
#' @param paga_show_transition
#' @param velocity
#' @param velocity_plot_type
#' @param velocity_n_neighbors
#' @param velocity_density
#' @param velocity_smooth
#' @param velocity_scale
#' @param velocity_min_mass
#' @param velocity_cutoff_perc
#' @param velocity_arrow_color
#' @param velocity_arrow_angle
#' @param velocity_arrow_flank
#' @param streamline_L
#' @param streamline_minL
#' @param streamline_res
#' @param streamline_n
#' @param streamline_jitter
#' @param streamline_size
#' @param streamline_alpha
#' @param streamline_color
#' @param streamline_palette
#' @param streamline_palcolor
#' @param streamline_bg_color
#' @param streamline_bg_stroke
#' @param hex
#' @param hex.size
#' @param hex.count
#' @param hex.bins
#' @param hex.binwidth
#' @param raster
#' @param raster.dpi
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param align
#' @param axis
#' @param force
#'
#' @return A single ggplot object if combine = TRUE; otherwise, a list of ggplot objects
#'
#' @examples
#' data("pancreas_sub")
#' ClassDimPlot(pancreas_sub,
#'   group.by = "CellType", reduction = "UMAP", label = TRUE,
#'   cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Delta"]
#' )
#'
#' ClassDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", stat.by = "Phase")
#' ClassDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", stat.by = "Phase", stat_plot_type = "bar", stat_plot_label = TRUE, stat_plot_size = 0.15)
#'
#' ClassDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", hex = TRUE)
#'
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' ClassDimPlot(pancreas_sub, group.by = "CellType", reduction = "UMAP", graph = "Standardpca_SNN")
#'
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", show_plot = FALSE)
#' ExpDimPlot(pancreas_sub, features = paste0("Lineage", 1:3), reduction = "UMAP")
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3))
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_whiskers = TRUE)
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
#'
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga)
#' ClassDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2,
#'   label = TRUE, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
#'   paga = pancreas_sub@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
#'   show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
#' )
#'
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", mode = "stochastic", return_seurat = TRUE)
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga, paga_show_transition = TRUE)
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = NA, velocity = "stochastic")
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "grid")
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "grid", velocity_scale = 1.5)
#' ClassDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2, velocity = "stochastic", velocity_plot_type = "stream")
#' ClassDimPlot(pancreas_sub,
#'   group.by = "SubCellType", reduction = "UMAP", pt.size = 5, pt.alpha = 0.2,
#'   label = TRUE, label_insitu = TRUE,
#'   velocity = "stochastic", velocity_plot_type = "stream", velocity_arrow_color = "yellow",
#'   velocity_density = 2, velocity_smooth = 1, streamline_n = 20, streamline_color = "black",
#'   show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
#' )
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom dplyr group_by summarize "%>%" .data
#' @importFrom ggplot2 ggplot ggplotGrob aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_manual scale_fill_manual guides guide_legend geom_hex geom_path theme_void annotation_custom
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color new_scale_fill new_scale
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom cowplot plot_grid get_legend
#' @importFrom stats median loess aggregate
#' @importFrom utils askYesNo
#' @importFrom rlang %||%
#' @export
#'
ClassDimPlot <- function(srt, group.by = "orig.ident", reduction = NULL, dims = c(1, 2), split.by = NULL,
                         show_na = FALSE, show_stat = TRUE,
                         palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = NULL, pt.alpha = 1,
                         label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                         label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                         label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                         cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                         add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_color = NULL,
                         lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75,
                         lineages_palette = "Dark2", lineages_palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                         lineages_line_size = 1, lineages_line_bg = "white", lineages_line_bg_r = 0.5,
                         lineages_whiskers = FALSE, lineages_whiskers_size = 0.5, lineages_whiskers_alpha = 0.5,
                         stat.by = NULL, stat_type = "percent", stat_plot_type = "pie", stat_plot_size = 0.1,
                         stat_plot_palette = "Set1", stat_palcolor = NULL, stat_plot_alpha = 1, stat_plot_label = FALSE, stat_plot_label_size = 3,
                         graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
                         paga = NULL, paga_node_palette = "Paired", paga_node_size = 4,
                         paga_edge_threshold = 0.01, paga_edge_size = c(0.2, 1), paga_edge_color = "grey40", paga_edge_alpha = 0.5,
                         paga_transition_threshold = 0.01, paga_transition_size = c(0.2, 1), paga_transition_color = "black", paga_transition_alpha = 1, paga_show_transition = FALSE,
                         velocity = NULL, velocity_plot_type = "raw", velocity_n_neighbors = ceiling(ncol(srt) / 50),
                         velocity_density = 1, velocity_smooth = 0.5, velocity_scale = 1, velocity_min_mass = 1, velocity_cutoff_perc = 5,
                         velocity_arrow_color = "black", velocity_arrow_angle = 20, velocity_arrow_flank = 0.8,
                         streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15, streamline_jitter = 1,
                         streamline_size = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
                         streamline_bg_color = "white", streamline_bg_stroke = 0.5,
                         hex = FALSE, hex.size = 0.5, hex.count = TRUE, hex.bins = 50, hex.binwidth = NULL,
                         raster = NULL, raster.dpi = c(512, 512),
                         theme_use = "theme_scp", aspect.ratio = 1, title = NULL, subtitle = NULL,
                         xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                         legend.position = "right", legend.direction = "vertical",
                         combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr", force = FALSE) {
  set.seed(11)
  check_R("exaexa/scattermore")
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("")
  }
  for (i in unique(c(group.by, split.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
    if (isTRUE(show_na) && any(is.na(srt[[i, drop = TRUE]]))) {
      raw_levels <- unique(c(levels(srt[[i, drop = TRUE]]), "NA"))
      srt[[i, drop = TRUE]] <- as.character(srt[[i, drop = TRUE]])
      srt[[i, drop = TRUE]][is.na(srt[[i, drop = TRUE]])] <- "NA"
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = raw_levels)
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
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }

  dat_meta <- srt@meta.data[, unique(c(group.by, split.by)), drop = FALSE]
  nlev <- sapply(dat_meta, nlevels)
  nlev <- nlev[nlev > 50]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 50 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  reduction_key <- Key(srt[[reduction]])
  dat_dim <- Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
  dat_use <- cbind(dat_dim, dat_meta[row.names(dat_dim), , drop = FALSE])
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) > 1e5)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector")
    }
  }
  if (!is.null(stat.by)) {
    subplots <- ClassStatPlot(srt,
      stat.by = stat.by, group.by = group.by, split.by = split.by,
      stat_type = stat_type, plot_type = stat_plot_type,
      palette = stat_plot_palette, palcolor = stat_palcolor, alpha = stat_plot_alpha,
      label = stat_plot_label, label.size = stat_plot_label_size,
      legend.position = "bottom", legend.direction = legend.direction,
      stat_single = TRUE, combine = FALSE
    )
  }
  if (!is.null(lineages)) {
    lineages_layers <- LineagePlot(srt,
      lineages = lineages, reduction = reduction, dims = dims,
      trim = lineages_trim, span = lineages_span,
      palette = lineages_palette, palcolor = lineages_palcolor, lineages_arrow = lineages_arrow,
      line_size = lineages_line_size, line_bg = lineages_line_bg, line_bg_r = lineages_line_bg_r,
      whiskers = lineages_whiskers, whiskers_size = lineages_whiskers_size, whiskers_alpha = lineages_whiskers_alpha,
      theme_use = theme_use, aspect.ratio = aspect.ratio, title = title, subtitle = subtitle,
      xlab = xlab, ylab = ylab, lab_cex = lab_cex, xlen_npc = xlen_npc, ylen_npc = ylen_npc,
      legend.position = legend.position, legend.direction = legend.direction,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(paga)) {
    if (split.by != "All_cells") {
      stop("paga can only plot on the non-split data")
    }
    paga_layers <- PAGAPlot(srt,
      paga = paga, reduction = reduction, dims = dims,
      node_palette = paga_node_palette, node_size = paga_node_size,
      edge_threshold = paga_edge_threshold, edge_size = paga_edge_size, edge_color = paga_edge_color, edge_alpha = paga_edge_alpha,
      transition_threshold = paga_transition_threshold, transition_size = paga_transition_size, transition_color = paga_transition_color, transition_alpha = paga_transition_alpha, show_transition = paga_show_transition,
      theme_use = theme_use, aspect.ratio = aspect.ratio, title = title, subtitle = subtitle,
      xlab = xlab, ylab = ylab, lab_cex = lab_cex, xlen_npc = xlen_npc, ylen_npc = ylen_npc,
      legend.position = "bottom", legend.direction = legend.direction, return_layer = TRUE
    )
    paga_layers <- paga_layers[!names(paga_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(velocity)) {
    if (split.by != "All_cells") {
      stop("velocity can only plot on the non-split data")
    }
    velocity_layers <- VelocityPlot(srt,
      reduction = reduction, dims = dims, velocity = velocity, plot_type = velocity_plot_type, group_by = group.by, group_palette = palette,
      n_neighbors = velocity_n_neighbors, density = velocity_density, smooth = velocity_smooth, scale = velocity_scale, min_mass = velocity_min_mass, cutoff_perc = velocity_cutoff_perc,
      arrow_color = velocity_arrow_color, arrow_angle = velocity_arrow_angle, arrow_flank = velocity_arrow_flank,
      streamline_L = streamline_L, streamline_minL = streamline_minL, streamline_res = streamline_res, streamline_jitter = streamline_jitter, streamline_n = streamline_n,
      streamline_size = streamline_size, streamline_alpha = streamline_alpha, streamline_color = streamline_color, streamline_palette = streamline_palette, streamline_palcolor = streamline_palcolor,
      streamline_bg_color = streamline_bg_color, streamline_bg_stroke = streamline_bg_stroke,
      theme_use = theme_use, aspect.ratio = aspect.ratio, title = title, subtitle = subtitle,
      xlab = xlab, ylab = ylab, lab_cex = lab_cex, xlen_npc = xlen_npc, ylen_npc = ylen_npc,
      legend.position = "bottom", legend.direction = legend.direction, return_layer = TRUE
    )
    velocity_layers <- velocity_layers[!names(velocity_layers) %in% c("lab_layer", "theme_layer")]
  }
  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  for (g in group.by) {
    colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor)
    for (s in levels(dat_use[[split.by]])) {
      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      legend_list <- list()
      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
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
      dat[, "cell"] <- rownames(dat)
      dat[, "x"] <- dat[, paste0(reduction_key, dims[1])]
      dat[, "y"] <- dat[, paste0(reduction_key, dims[2])]
      dat[, "group.by"] <- dat[, g]
      dat[, "split.by"] <- s
      dat <- dat[order(dat[, "group.by"], decreasing = FALSE, na.last = FALSE), ]
      naindex <- which(is.na(dat[, "group.by"]))
      naindex <- ifelse(length(naindex) > 0, max(naindex), 1)
      dat <- dat[c(1:naindex, sample((min(naindex + 1, nrow(dat))):nrow(dat))), ]
      if (isTRUE(show_stat)) {
        subtitle_use <- subtitle %||% paste0(s, " nCell:", sum(!is.na(dat[["group.by"]])))
      } else {
        subtitle_use <- subtitle
      }
      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_color)
          density <- list(
            stat_density_2d(
              geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = ..density..),
              contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
            ),
            scale_fill_gradientn(name = "Density", colours = filled_color)
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
        net_mat <- as.matrix(x = srt[[graph]])[rownames(dat), rownames(dat)]
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
            data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = value),
            color = edge_color, alpha = edge_alpha, show.legend = FALSE
          ),
          scale_size_continuous(range = edge_size)
        )
      } else {
        net <- NULL
      }
      p <- ggplot(dat) +
        net +
        density +
        labs(title = title, subtitle = subtitle_use, x = xlab, y = ylab) +
        scale_x_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[1])]), max(dat_dim[, paste0(reduction_key, dims[1])]))) +
        scale_y_continuous(limits = c(min(dat_dim[, paste0(reduction_key, dims[2])]), max(dat_dim[, paste0(reduction_key, dims[2])]))) +
        facet_grid(. ~ split.by) +
        do.call(theme_use, list(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
        ))

      if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
          pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
        )
      } else if (isTRUE(hex)) {
        check_R("hexbin")
        if (isTRUE(hex.count)) {
          p <- p + geom_hex(
            mapping = aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["group.by"]], color = .data[["group.by"]], alpha = ..count..),
            size = hex.size, bins = hex.bins, binwidth = hex.binwidth
          )
        } else {
          p <- p + geom_hex(
            mapping = aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["group.by"]], color = .data[["group.by"]]),
            size = hex.size, bins = hex.bins, binwidth = hex.binwidth
          )
        }
      } else {
        p <- p + suppressWarnings(geom_point(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]], cell = .data[["cell"]]),
          size = pt.size, alpha = pt.alpha
        ))
      }
      if (!is.null(cells.highlight) && !isTRUE(hex)) {
        p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
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
            p <- p +
              suppressWarnings(geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], cell = .data[["cell"]]), color = cols.highlight,
                size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
              )) +
              suppressWarnings(geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]], cell = .data[["cell"]]),
                size = sizes.highlight, alpha = alpha.highlight
              ))
          }
        }
      }
      p <- p + scale_color_manual(
        name = paste0(g, ":"),
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color, guide = guide_legend(
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
        stat_plot <- subplots[paste0(s, ":", levels(dat[, "group.by"]))]
        names(stat_plot) <- levels(dat[, "group.by"])

        stat_plot_list <- list()
        for (i in seq_len(nrow(coor_df))) {
          stat_plot_list[[i]] <- annotation_custom(ggplotGrob(stat_plot[[coor_df[i, "group"]]] + theme_void() + theme(legend.position = "none")),
            x = coor_df[i, "x"] - x_range * stat_plot_size / 2, y = coor_df[i, "y"] - y_range * stat_plot_size / 2,
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
        label_df <- label_df[!is.na(label_df[, "label"]), ]
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
        grob <- ggplotGrob(p + theme(legend.position = "none"))
        if (legend.position == "bottom") {
          grob <- gtable_add_rows(grob, sum(legend$heights), -1)
          grob <- gtable_add_grob(grob, legend, t = dim(grob)[1], l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
        }
        if (legend.position == "top") {
          grob <- gtable_add_rows(grob, sum(legend$heights), 0)
          grob <- gtable_add_grob(grob, legend, t = 1, l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
        }
        if (legend.position == "right") {
          grob <- gtable_add_cols(grob, sum(legend$widths), -1)
          grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])
        }
        if (legend.position == "left") {
          grob <- gtable_add_cols(grob, sum(legend$widths), 0)
          grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1)
        }
        p <- plot_grid(grob)
      }
      plist[[paste0(s, ":", g)]] <- p
    }
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' Dimensional reduction plot for expression.
#'
#' Plotting cell points on a reduced 2D plane and coloring according to the expression of the features.
#'
#' @param srt A Seurat object.
#' @param features Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' @param reduction Which dimensionality reduction to use.
#' @param split.by Name of a metadata column to split plot by.
#' @param palette Name of palette to use. Default is "Paired".
#' @param bg_color Color value for NA points.
#' @param pt.size Point size for plotting.
#' @param pt.alpha Point transparency.
#' @param keep_scale How to handle the color scale across multiple plots. Options are:
#' \itemize{
#'   \item{"feature" (default; by row/feature scaling):}{ The plots for each individual feature are scaled to the maximum expression of the feature across the conditions provided to 'split.by'.}
#'   \item{"all" (universal scaling):}{ The plots for all features and conditions are scaled to the maximum expression value for the feature with the highest overall expression.}
#'   \item{NULL (no scaling):}{ Each individual plot is scaled to the maximum expression value of the feature in the condition provided to 'split.by'. Be aware setting NULL will result in color scales that are not comparable between plots.}
#' }
#' @param cells.highlight A vector of names of cells to highlight.
#' @param cols.highlight A vector of colors to highlight the cells.
#' @param sizes.highlight Size of highlighted cells.
#' @param legend.position The position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.direction Layout of items in legends ("horizontal" or "vertical")
#' @param combine Whether to arrange multiple plots into a grid.
#' @param nrow Number of rows in the plot grid.
#' @param ncol Number of columns in the plot grid.
#' @param byrow Logical value indicating if the plots should be arrange by row (default) or by column.
#' @param dims
#' @param slot
#' @param assay
#' @param show_stat
#' @param palcolor
#' @param bg_cutoff
#' @param lower_quantile
#' @param upper_quantile
#' @param add_density
#' @param density_color
#' @param density_filled
#' @param density_filled_palette
#' @param density_filled_color
#' @param alpha.highlight
#' @param stroke.highlight
#' @param calculate_coexp
#' @param compare_features
#' @param color_blend_mode
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
#' @param lineages
#' @param lineages_trim
#' @param lineages_span
#' @param lineages_palette
#' @param lineages_palcolor
#' @param lineages_arrow
#' @param lineages_line_size
#' @param lineages_line_bg
#' @param lineages_line_bg_r
#' @param lineages_whiskers
#' @param lineages_whiskers_size
#' @param lineages_whiskers_alpha
#' @param graph
#' @param edge_size
#' @param edge_alpha
#' @param edge_color
#' @param hex
#' @param hex.size
#' @param hex.color
#' @param hex.bins
#' @param hex.binwidth
#' @param raster
#' @param raster.dpi
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param align
#' @param axis
#' @param force
#'
#' @return A single ggplot object if combine = TRUE; otherwise, a list of ggplot objects
#'
#' @examples
#' data("pancreas_sub")
#' ExpDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP")
#' ExpDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", hex = TRUE)
#'
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' ExpDimPlot(pancreas_sub, features = "Lineage3", reduction = "UMAP", lineages = "Lineage3")
#' ExpDimPlot(pancreas_sub, features = "Lineage3", reduction = "UMAP", lineages = "Lineage3", lineages_span = 0.1)
#'
#' ExpDimPlot(pancreas_sub,
#'   features = c("Ghrl", "Ins1", "Gcg", "Ins2"), pt.size = 1,
#'   compare_features = TRUE, color_blend_mode = "blend",
#'   label = TRUE, label_insitu = TRUE
#' )
#'
#' ExpDimPlot(pancreas_sub,
#'   features = c("Ghrl", "Ins1", "Gcg", "Ins2"), pt.size = 1,
#'   compare_features = TRUE, color_blend_mode = "blend",
#'   label = TRUE, label_insitu = TRUE, label_repel = TRUE
#' )
#'
#' ExpDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "mix", title = "mix",
#'   label = TRUE, label_insitu = TRUE
#' )
#' ExpDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "blend", title = "blend",
#'   label = TRUE, label_insitu = TRUE
#' )
#' ExpDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "screen", title = "screen",
#'   label = TRUE, label_insitu = TRUE
#' )
#' ExpDimPlot(pancreas_sub,
#'   features = c("S_score", "G2M_score"), pt.size = 1, palcolor = c("red", "green"),
#'   compare_features = TRUE, color_blend_mode = "multiply", title = "multiply",
#'   label = TRUE, label_insitu = TRUE
#' )
#'
#' library(scales)
#' par(mfrow = c(2, 2))
#' show_col(c("red", "green", blendcolors(c("red", "green"), mode = "mix")), ncol = 4)
#' text(3.5, -0.5, "mix", cex = 1.2)
#' show_col(c("red", "green", blendcolors(c("red", "green"), mode = "blend")), ncol = 4)
#' text(3.5, -0.5, "blend", cex = 1.2)
#' show_col(c("red", "green", blendcolors(c("red", "green"), mode = "screen")), ncol = 4)
#' text(3.5, -0.5, "screen", cex = 1.2)
#' show_col(c("red", "green", blendcolors(c("red", "green"), mode = "multiply")), ncol = 4)
#' text(3.5, -0.5, "multiply", cex = 1.2)
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom dplyr group_by summarize "%>%" .data
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colourbar geom_hex stat_summary_hex geom_path
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom gtable gtable_add_cols
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid get_legend draw_grob
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @export
#'
ExpDimPlot <- function(srt, features, reduction = NULL, dims = c(1, 2), split.by = NULL, slot = "data", assay = DefaultAssay(srt),
                       show_stat = TRUE,
                       palette = ifelse(isTRUE(compare_features), "Set1", "Spectral"), palcolor = NULL,
                       bg_cutoff = 0, bg_color = "grey80", pt.size = NULL, pt.alpha = 1,
                       keep_scale = NULL, lower_quantile = 0, upper_quantile = 0.99,
                       add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_color = NULL,
                       cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                       calculate_coexp = FALSE, compare_features = FALSE, color_blend_mode = "blend",
                       label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                       label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                       label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                       lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75,
                       lineages_palette = "Dark2", lineages_palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                       lineages_line_size = 1, lineages_line_bg = "white", lineages_line_bg_r = 0.5,
                       lineages_whiskers = FALSE, lineages_whiskers_size = 0.5, lineages_whiskers_alpha = 0.5,
                       graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
                       hex = FALSE, hex.size = 0.5, hex.color = "grey90", hex.bins = 50, hex.binwidth = NULL,
                       raster = NULL, raster.dpi = c(512, 512),
                       theme_use = "theme_scp", aspect.ratio = 1, title = NULL, subtitle = NULL,
                       xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                       legend.position = "right", legend.direction = "vertical",
                       combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr", force = FALSE) {
  check_R("exaexa/scattermore")
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("")
  }
  for (i in c(split.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
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
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }

  features <- unique(features)
  features_drop <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt[[assay]]@counts)]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression")
    }
    status <- check_DataType(srt, slot = slot)
    message("Data type detected in ", slot, " slot: ", status)
    if (status == "raw_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "log_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), expm1) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "raw_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), log1p) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene > 0)) {
    dat_gene <- t(GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE])[colnames(srt), , drop = FALSE]
  } else {
    dat_gene <- matrix(nrow = ncol(srt), ncol = 0)
  }
  if (length(features_meta > 0)) {
    dat_meta <- srt@meta.data[, features_meta, drop = FALSE]
  } else {
    dat_meta <- matrix(nrow = ncol(srt), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!all(sapply(dat_exp, is.numeric))) {
    stop("'features' must be type of numeric variable.")
  }
  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  reduction_key <- Key(srt[[reduction]])
  dat_dim <- Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
  dat_sp <- srt@meta.data[, split.by, drop = FALSE]
  dat_use <- cbind(dat_dim, dat_sp[row.names(dat_dim), , drop = FALSE])

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) > 1e5)
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
      line_size = lineages_line_size, line_bg = lineages_line_bg, line_bg_r = lineages_line_bg_r,
      whiskers = lineages_whiskers, whiskers_size = lineages_whiskers_size, whiskers_alpha = lineages_whiskers_alpha,
      theme_use = theme_use, aspect.ratio = aspect.ratio, title = title, subtitle = subtitle,
      xlab = xlab, ylab = ylab, lab_cex = lab_cex, xlen_npc = xlen_npc, ylen_npc = ylen_npc,
      legend.position = legend.position, legend.direction = legend.direction,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
  }

  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  if (isTRUE(compare_features) && length(features) > 1) {
    f <- features
    for (s in levels(dat_sp[[split.by]])) {
      dat <- cbind(dat_use, dat_exp[row.names(dat_use), f, drop = FALSE])
      dat[, f][dat[, f] <= bg_cutoff] <- NA
      dat[, "cell"] <- rownames(dat)
      dat[, "x"] <- dat[, paste0(reduction_key, dims[1])]
      dat[, "y"] <- dat[, paste0(reduction_key, dims[2])]
      dat[, "split.by"] <- s
      dat[, "features"] <- paste(f, collapse = "|")
      cells_keep <- dat[[split.by]] == s
      dat <- dat[cells_keep, ]
      subtitle_use <- subtitle %||% s

      colors <- palette_scp(f, type = "discrete", palette = palette, palcolor = palcolor)
      colors_list <- list()
      value_list <- list()
      pal_list <- list()
      temp_geom <- list()
      legend_list <- list()
      for (i in seq_along(colors)) {
        lower_col <- adjcolors(colors[i], 0.1)
        colors_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = TRUE, matched = TRUE, palcolor = c(lower_col, colors[i]))
        pal_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = FALSE, matched = FALSE, palcolor = c(lower_col, colors[i]))
        value_list[[i]] <- seq(min(dat[, names(colors)[i]], na.rm = TRUE), max(dat[, names(colors)[i]], na.rm = TRUE), length.out = 100)
        temp_geom[[i]] <- list(
          suppressWarnings(geom_point(data = dat, mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[[names(colors)[i]]], cell = .data[["cell"]]))),
          scale_color_gradientn(
            colours = pal_list[[i]],
            values = rescale(value_list[[i]]), na.value = bg_color,
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
          ),
          new_scale_color()
        )
        legend_list[[i]] <- get_legend(ggplot(dat, aes(x = .data[["x"]], y = .data[["y"]])) +
          temp_geom[[i]] +
          do.call(theme_use, list(
            legend.spacing = unit(0, "cm"),
            legend.margin = margin(5, 5, 5, 5),
            xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
          )))
      }
      for (j in seq_len(nrow(dat))) {
        dat[j, "color_blend"] <- blendcolors(sapply(colors_list, function(x) x[j]), mode = color_blend_mode)
      }
      dat["color_value"] <- colSums(col2rgb(dat[, "color_blend"]))
      dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_value"] <- NA
      dat <- dat[order(dat[, "color_value"], decreasing = TRUE, na.last = FALSE), ]
      dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_blend"] <- bg_color
      if (!is.null(graph)) {
        net_mat <- as.matrix(x = srt[[graph]])[rownames(dat), rownames(dat)]
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
            data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = value),
            color = edge_color, alpha = edge_alpha, show.legend = FALSE
          ),
          scale_size_continuous(range = edge_size)
        )
      } else {
        net <- NULL
      }
      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_color)
          density <- list(
            stat_density_2d(
              geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = ..density..),
              contour = FALSE,
              inherit.aes = FALSE
            ),
            scale_fill_gradientn(name = "Density", colours = filled_color)
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
        scale_x_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[1])]), max(dat_use[, paste0(reduction_key, dims[1])]))) +
        scale_y_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[2])]), max(dat_use[, paste0(reduction_key, dims[2])]))) +
        facet_grid(split.by ~ features) +
        do.call(theme_use, list(
          aspect.ratio = aspect.ratio,
          legend.position = "none", legend.direction = legend.direction,
          xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
        ))
      if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
          pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
        ) +
          scale_color_identity() +
          new_scale_color()
      } else {
        p <- p + suppressWarnings(geom_point(
          mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]], cell = .data[["cell"]]),
          size = pt.size, alpha = pt.alpha
        )) +
          scale_color_identity() +
          new_scale_color()
      }

      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        if (nrow(cell_df) > 0) {
          if (isTRUE(raster)) {
            p <- p + scattermore::geom_scattermore(
              data = cell_df, aes(x = .data[["x"]], y = .data[["y"]]), color = cols.highlight,
              pointsize = floor(sizes.highlight) + stroke.highlight, alpha = alpha.highlight, pixels = raster.dpi
            ) +
              scattermore::geom_scattermore(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
                pointsize = floor(sizes.highlight), alpha = alpha.highlight, pixels = raster.dpi
              ) + scale_color_identity() +
              new_scale_color()
          } else {
            p <- p +
              suppressWarnings(geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], cell = .data[["cell"]]), color = cols.highlight,
                size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
              )) +
              suppressWarnings(geom_point(
                data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]], cell = .data[["cell"]]),
                size = sizes.highlight, alpha = alpha.highlight
              )) +
              scale_color_identity() +
              new_scale_color()
          }
        }
      }

      legend2 <- NULL
      if (isTRUE(label)) {
        label_df <- p$data %>%
          reshape2::melt(measure.vars = f) %>%
          group_by(variable) %>%
          filter(value >= quantile(value, 0.95, na.rm = TRUE) & value <= quantile(value, 0.99, na.rm = TRUE)) %>%
          summarize(x = median(.data[["x"]]), y = median(.data[["y"]])) %>%
          as.data.frame()
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[, "label"]), ]
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
          legend2 <- get_legend(p + guides(fill = guide_colourbar(label.theme = element_text(angle = 45, hjust = 1, vjust = 1))) +
            do.call(theme_use, list(
              legend.direction = "horizontal",
              legend.spacing = unit(0, "cm"),
              legend.margin = margin(5, 5, 5, 5),
              xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
            )))
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
        legend <- gtable_add_rows(legend, sum(legend_curve$heights), 0)
        legend <- gtable_add_grob(legend, legend_curve,
          t = 1,
          l = min(legend$layout[grepl(pattern = "guides", x = legend$layout$name), "l"])
        )
        p <- suppressMessages({
          p + lineages_layers + theme(legend.position = "none")
        })
      }

      grob <- ggplotGrob(p)
      grob <- gtable_add_cols(grob, sum(legend$widths), -1)
      grob <- gtable_add_grob(grob, legend, t = grob$layout[grepl(pattern = "panel", x = grob$layout$name), "t"], l = dim(grob)[2])
      if (!is.null(legend2)) {
        grob <- gtable_add_rows(grob, sum(legend2$heights), -1)
        grob <- gtable_add_grob(grob, legend2,
          t = dim(grob)[1],
          l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"])
        )
      }
      p <- plot_grid(grob)
      plist[[paste0(s, ":", paste0(f, collapse = "|"))]] <- p
    }
  } else {
    for (f in features) {
      for (s in levels(dat_sp[[split.by]])) {
        dat <- cbind(dat_use, dat_exp[row.names(dat_use), f, drop = FALSE])
        dat[, f][dat[, f] <= bg_cutoff] <- NA
        dat[, "cell"] <- rownames(dat)
        dat[, "x"] <- dat[, paste0(reduction_key, dims[1])]
        dat[, "y"] <- dat[, paste0(reduction_key, dims[2])]
        dat[, "value"] <- dat[, f]
        dat[, "features"] <- f
        cells_keep <- dat[[split.by]] == s
        dat <- dat[cells_keep, ]
        dat <- dat[order(dat[, "value"], decreasing = FALSE, na.last = FALSE), ]
        colors <- palette_scp(dat[, f], type = "continuous", palette = palette, palcolor = palcolor)
        legend_list <- list()
        if (isTRUE(show_stat)) {
          subtitle_use <- subtitle %||% paste0(s, " nPos:", sum(dat[["value"]] > 0, na.rm = TRUE), ", ", round(sum(dat[["value"]] > 0, na.rm = TRUE) / nrow(dat) * 100, 2), "%")
        } else {
          subtitle_use <- subtitle
        }
        if (all(is.na(dat[, f]))) {
          colors_value <- rep(0, 100)
        } else {
          if (is.null(keep_scale)) {
            colors_value <- seq(quantile(dat[, f], lower_quantile, na.rm = TRUE), quantile(dat[, f], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
          } else {
            if (keep_scale == "feature") {
              colors_value <- seq(quantile(dat_exp[, f], lower_quantile, na.rm = TRUE), quantile(dat_exp[, f], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
            }
            if (keep_scale == "all") {
              colors_value <- seq(quantile(dat_exp[, features], lower_quantile, na.rm = TRUE), quantile(dat_exp[, features], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
            }
          }
        }
        dat[which(dat[, "value"] > max(colors_value)), "value"] <- max(colors_value)
        if (!is.null(graph)) {
          net_mat <- as.matrix(x = srt[[graph]])[rownames(dat), rownames(dat)]
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
              data = net_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = value),
              color = edge_color, alpha = edge_alpha, show.legend = FALSE
            ),
            scale_size_continuous(range = edge_size)
          )
        } else {
          net <- NULL
        }
        if (isTRUE(add_density)) {
          if (isTRUE(density_filled)) {
            filled_color <- palette_scp(palette = density_filled_palette, palcolor = density_filled_color)
            density <- list(
              stat_density_2d(
                geom = "raster", aes(x = .data[["x"]], y = .data[["y"]], fill = ..density..),
                contour = FALSE,
                inherit.aes = FALSE
              ),
              scale_fill_gradientn(name = "Density", colours = filled_color)
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
          scale_x_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[1])]), max(dat_use[, paste0(reduction_key, dims[1])]))) +
          scale_y_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[2])]), max(dat_use[, paste0(reduction_key, dims[2])]))) +
          guides(color = guide_colourbar(
            barwidth = 0.9,
            barheight = 4,
            frame.colour = "black",
            ticks.colour = "black"
          )) +
          do.call(theme_use, list(
            aspect.ratio = aspect.ratio,
            legend.position = legend.position,
            legend.direction = legend.direction,
            xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
          ))
        if (isTRUE(raster)) {
          p <- p + scattermore::geom_scattermore(
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
            pointsize = floor(pt.size), alpha = pt.alpha, pixels = raster.dpi
          )
        } else if (isTRUE(hex)) {
          check_R("hexbin")
          dat_na <- dat[is.na(dat[["value"]]), ]
          dat_hex <- dat[!is.na(dat[["value"]]), ]
          if (nrow(dat_na) > 0) {
            p <- p + geom_hex(
              data = dat[is.na(dat[["value"]]), ],
              mapping = aes(x = .data[["x"]], y = .data[["y"]]),
              fill = bg_color, color = hex.color,
              size = hex.size, bins = hex.bins, binwidth = hex.binwidth
            ) + new_scale_fill()
          }
          if (nrow(dat_hex) > 0) {
            p <- p + stat_summary_hex(
              data = dat_hex,
              mapping = aes(x = .data[["x"]], y = .data[["y"]], z = .data[["value"]]),
              color = hex.color, size = hex.size, bins = hex.bins, binwidth = hex.binwidth
            )
          }
        } else {
          p <- p + suppressWarnings(geom_point(
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]], cell = .data[["cell"]]),
            size = pt.size, alpha = pt.alpha
          ))
        }
        if (!is.null(cells.highlight) && !isTRUE(hex)) {
          p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
          cell_df <- subset(p$data, cells.highlight == TRUE)
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
                suppressWarnings(geom_point(
                  data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], cell = .data[["cell"]]), color = cols.highlight,
                  size = sizes.highlight + stroke.highlight, alpha = alpha.highlight
                )) +
                suppressWarnings(geom_point(
                  data = cell_df, aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]], cell = .data[["cell"]]),
                  size = sizes.highlight, alpha = alpha.highlight
                ))
            }
          }
        }
        if (nrow(dat) > 0) {
          p <- p + facet_grid(formula(paste0(split.by, "~features")))
        }
        if (all(is.na(dat[, f]))) {
          p <- p + scale_colour_gradient(
            name = "", na.value = bg_color
          ) + scale_fill_gradient(
            name = "", na.value = bg_color
          )
        } else {
          p <- p + scale_color_gradientn(
            name = "", colours = colors, values = rescale(colors_value), na.value = bg_color
          ) + scale_fill_gradientn(
            name = "", colours = colors, values = rescale(colors_value), na.value = bg_color
          )
        }
        p <- p + guides(
          color = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1),
          fill = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
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
            summarize(x = median(.data[["x"]]), y = median(.data[["y"]])) %>%
            as.data.frame()
          label_df[, "label"] <- f
          label_df[, "rank"] <- seq_len(nrow(label_df))
          require(ggrepel)
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
          grob <- ggplotGrob(p + theme(legend.position = "none"))
          if (legend.position == "bottom") {
            grob <- gtable_add_rows(grob, sum(legend$heights), -1)
            grob <- gtable_add_grob(grob, legend, t = dim(grob)[1], l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
          }
          if (legend.position == "top") {
            grob <- gtable_add_rows(grob, sum(legend$heights), 0)
            grob <- gtable_add_grob(grob, legend, t = 1, l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
          }
          if (legend.position == "right") {
            grob <- gtable_add_cols(grob, sum(legend$widths), -1)
            grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])
          }
          if (legend.position == "left") {
            grob <- gtable_add_cols(grob, sum(legend$widths), 0)
            grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1)
          }
          p <- plot_grid(grob)
        }
        plist[[paste0(s, ":", f)]] <- p
      }
    }
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
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
#' @inheritParams ClassDimPlot
#' @param shape.highlight Shape of the cell to highlight. See \href{https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol}{scattergl-marker-symbol}
#' @param width Width in pixels, defaults to automatic sizing.
#' @param height Height in pixels, defaults to automatic sizing.
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' ClassDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#'
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#' ClassDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D", lineages = paste0("Lineage", 1:2))
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom utils askYesNo
#' @importFrom stringr str_replace
#' @export
ClassDimPlot3D <- function(srt, group.by = "orig.ident", reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                           palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = 1.5,
                           cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                           lineages = NULL, lineage_palette = "Dark2", span = 0.75, arrow_reverse = FALSE,
                           width = NULL, height = NULL, save = NULL, force = FALSE) {
  bg_color <- col2hex(bg_color)
  cols.highlight <- col2hex(cols.highlight)

  for (i in c(group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
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
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (ncol(Embeddings(srt, reduction = reduction)) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }
  reduction_key <- Key(srt[[reduction]])
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

  dat_dim <- Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
  dat_use <- cbind(dat_dim[colnames(srt), ], srt@meta.data[colnames(srt), , drop = FALSE])
  nlev <- sapply(dat_use[, group.by, drop = FALSE], nlevels)
  nlev <- nlev[nlev > 50]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 50 levels.", immediate. = TRUE)
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
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(dat_use)]
    dat_use_highlight <- dat_use[cells.highlight, ]
    dat_use_highlight[["group.by"]] <- "highlight"
  }
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
  if (!is.null(cells.highlight)) {
    cells.highlight <- cells.highlight[cells.highlight %in% rownames(dat_use)]
    dat_use_highlight <- dat_use[cells.highlight, ]
  }

  p <- plot_ly(data = dat_use, width = width, height = height)
  p <- p %>% add_trace(
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
  if (!is.null(cells.highlight)) {
    p <- p %>% add_trace(
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
      dat_sub <- dat_use[!is.na(dat_use[[l]]), ]
      dat_sub <- dat_sub[order(dat_sub[[l]]), ]

      xlo <- loess(formula(paste(paste0(reduction_key, dims[1], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      ylo <- loess(formula(paste(paste0(reduction_key, dims[2], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      zlo <- loess(formula(paste(paste0(reduction_key, dims[3], "All_cells"), l, sep = "~")), data = dat_sub, span = span, degree = 2)
      dat_smooth <- data.frame(x = xlo$fitted, y = ylo$fitted, z = zlo$fitted)
      dat_smooth <- dat_smooth[dat_smooth[["x"]] <= max(dat_use[[paste0(reduction_key, dims[1], "All_cells")]]) & dat_smooth[["x"]] >= min(dat_use[[paste0(reduction_key, dims[1], "All_cells")]]), ]
      dat_smooth <- dat_smooth[dat_smooth[["y"]] <= max(dat_use[[paste0(reduction_key, dims[2], "All_cells")]]) & dat_smooth[["y"]] >= min(dat_use[[paste0(reduction_key, dims[2], "All_cells")]]), ]
      dat_smooth <- dat_smooth[dat_smooth[["z"]] <= max(dat_use[[paste0(reduction_key, dims[3], "All_cells")]]) & dat_smooth[["z"]] >= min(dat_use[[paste0(reduction_key, dims[3], "All_cells")]]), ]
      dat_smooth <- unique(na.omit(dat_smooth))
      p <- add_trace(
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

  p <- p %>% plotly::layout(
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
      xaxis = list(title = xlab, range = c(min(dat_use[[paste0(reduction_key, dims[1])]]), max(dat_use[[paste0(reduction_key, dims[1])]]))),
      yaxis = list(title = ylab, range = c(min(dat_use[[paste0(reduction_key, dims[2])]]), max(dat_use[[paste0(reduction_key, dims[2])]]))),
      zaxis = list(title = zlab, range = c(min(dat_use[[paste0(reduction_key, dims[3])]]), max(dat_use[[paste0(reduction_key, dims[3])]]))),
      aspectratio = list(x = 1, y = 1, z = 1)
    ),
    autosize = FALSE
  )

  if ((!is.null(save) && is.character(save) && nchar(save) > 0)) {
    htmlwidgets::saveWidget(
      widget = plotly::as_widget(p),
      file = save
    )
    unlink(str_replace(save, "\\.html", "_files"), recursive = TRUE)
  }

  return(p)
}

#' 3D-Dimensional reduction plot for gene expression visualization.
#'
#' Plotting cell points on a reduced 3D space and coloring according to the gene expression in the cells.
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' ExpDimPlot3D(pancreas_sub, features = c("Ghrl", "Ins1", "Gcg", "Ins2"), reduction = "StandardpcaUMAP3D")
#' @inheritParams ExpDimPlot
#' @inheritParams ClassDimPlot3D
#'
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom rlang "%||%"
#' @importFrom stringr str_replace
#' @export
ExpDimPlot3D <- function(srt, features = NULL, reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                         split.by = NULL, slot = "data", assay = "RNA",
                         calculate_coexp = FALSE,
                         pt.size = 1.5, cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                         width = NULL, height = NULL, save = NULL, force = FALSE) {
  cols.highlight <- col2hex(cols.highlight)

  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("All_cells")
  }
  for (i in split.by) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt, min_dim = 3)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction, min_dim = 3)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  if (ncol(Embeddings(srt, reduction = reduction)) < 3) {
    stop("Reduction must be in three dimensions or higher.")
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }
  reduction_key <- Key(srt[[reduction]])
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
  features_drop <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt[[assay]]@counts)]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression")
    }
    status <- check_DataType(srt, slot = slot)
    message("Data type detected in ", slot, " slot: ", status)
    if (status == "raw_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "log_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), expm1) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "raw_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), log1p) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene > 0)) {
    dat_gene <- t(GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE])[colnames(srt), , drop = FALSE]
  } else {
    dat_gene <- matrix(nrow = ncol(srt), ncol = 0)
  }
  if (length(features_meta > 0)) {
    dat_meta <- srt@meta.data[, features_meta, drop = FALSE]
  } else {
    dat_meta <- matrix(nrow = ncol(srt), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])

  if (!all(sapply(dat_exp, is.numeric))) {
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
  dat_dim <- Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
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
    dat_use_highlight <- dat_use[cells.highlight, ]
    for (i in levels(dat_use_highlight[[split.by]])) {
      dat_use_highlight[[paste0(reduction_key, dims[1], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[1])]], NA)
      dat_use_highlight[[paste0(reduction_key, dims[2], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[2])]], NA)
      dat_use_highlight[[paste0(reduction_key, dims[3], i)]] <- ifelse(dat_use_highlight[[split.by]] == i, dat_use_highlight[[paste0(reduction_key, dims[3])]], NA)
    }
  }

  p <- plot_ly(data = dat_use, width = width, height = height)
  p <- p %>% add_trace(
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
    p <- p %>% add_trace(
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

  p <- p %>% plotly::layout(
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
      xaxis = list(title = xlab, range = c(min(dat_use[[paste0(reduction_key, dims[1])]]), max(dat_use[[paste0(reduction_key, dims[1])]]))),
      yaxis = list(title = ylab, range = c(min(dat_use[[paste0(reduction_key, dims[2])]]), max(dat_use[[paste0(reduction_key, dims[2])]]))),
      zaxis = list(title = zlab, range = c(min(dat_use[[paste0(reduction_key, dims[3])]]), max(dat_use[[paste0(reduction_key, dims[3])]]))),
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
    unlink(str_replace(save, "\\.html", "_files"), recursive = TRUE)
  }

  return(p)
}


# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", ggplot2::GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- dplyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
        1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

#' Violin plot for cell expression.
#'
#' @param srt
#' @param features
#' @param group.by
#' @param split.by
#' @param bg.by
#' @param cells_subset
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
#' @param box_fill
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
#' @param adjust
#' @param split.plot
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
#' @param align
#' @param axis
#' @param force
#'
#' @examples
#' data("pancreas_sub")
#' ExpVlnPlot(pancreas_sub,
#'   features = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   group.by = "SubCellType", bg.by = "CellType", stack = TRUE, flip = TRUE
#' )
#' ExpVlnPlot(pancreas_sub, features = c("G2M_score", "S_score"), group.by = "SubCellType")
#' ExpVlnPlot(pancreas_sub, features = c("Neurog3", "Fev"), group.by = "SubCellType", bg.by = "CellType", stack = TRUE)
#' ExpVlnPlot(pancreas_sub, features = c("Rbp4", "Pyy"), group.by = "SubCellType", comparisons = list(c("Alpha", "Beta"), c("Alpha", "Delta")), multiplegroup_comparisons = TRUE)
#' @importFrom Seurat DefaultAssay
#' @importFrom gtable gtable_add_cols gtable_add_rows gtable_add_grob gtable_add_padding
#' @importFrom ggplot2 geom_blank geom_violin geom_rect geom_boxplot layer_scales position_jitterdodge position_dodge2 stat_summary scale_x_discrete element_line annotate
#' @importFrom grid grobHeight grobWidth
#' @importFrom rlang %||%
#' @importFrom cowplot plot_grid
#' @export
ExpVlnPlot <- function(srt, features = NULL, group.by = NULL, split.by = NULL, bg.by = NULL,
                       cells_subset = NULL, keep_empty = FALSE,
                       slot = "data", assay = DefaultAssay(srt),
                       palette = "Paired", palcolor = NULL, alpha = 1,
                       bg_palette = "Paired", bg_palcolor = NULL, bg_apha = 0.5,
                       add_box = TRUE, box_fill = "black", box_width = 0.2,
                       add_point = FALSE, pt.color = "black", pt.size = NULL, pt.alpha = 1, jitter.width = 1,
                       cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1, alpha.highlight = 1,
                       calculate_coexp = FALSE, compare_features = FALSE,
                       y.nbreaks = 3, y.max = NULL, same.y.lims = FALSE, y.trans = "identity",
                       sort = FALSE, adjust = 1, split.plot = FALSE,
                       stack = FALSE, fill.by = "ident", flip = FALSE,
                       comparisons = list(), ref_group = NULL, pairwise_method = "wilcox.test",
                       multiplegroup_comparisons = FALSE, multiple_method = "kruskal.test",
                       theme_use = "theme_scp", aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = group.by, ylab = "Expression level",
                       legend.position = "right", legend.direction = "vertical",
                       combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr", force = FALSE) {
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(group.by)) {
    stop("'group.by' must be provided.")
  }
  if (!(fill.by %in% c("feature", "ident"))) {
    stop("`fill.by` must be either `feature` or `ident`")
  }
  if (isTRUE(stack)) {
    if (fill.by == "feature" && !is.null(split.by)) {
      warning("split.by is not used when stacking violin plots with fill.by='feature'")
    }
    split.by <- NULL
    if (isTRUE(sort)) {
      sort <- FALSE
    }
  }
  if (is.null(split.by)) {
    split.by <- group.by
  }
  for (i in unique(c(group.by, split.by, bg.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
  }
  if (!is.null(bg.by)) {
    df_table <- table(srt[[group.by, drop = TRUE]], srt[[bg.by, drop = TRUE]])
    if (max(rowSums(df_table > 0)) > 1) {
      stop("'group.by' must be a division of 'bg.by'")
    }
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }
  if (!is.null(cells_subset)) {
    if (!any(cells_subset %in% colnames(srt))) {
      stop("No cells in 'cells_subset' found in srt.")
    }
    if (!all(cells_subset %in% colnames(srt))) {
      warning("Some cells in 'cells_subset' not found in srt.", immediate. = TRUE)
    }
    cells_subset <- intersect(cells_subset, colnames(srt))
    srt <- srt[, cells_subset]
  }
  if (length(comparisons) > 0) {
    l <- sapply(comparisons, length)
    if (any(l > 2)) {
      stop("'comparisons' must be a list of length-2 vectors.")
    }
  }

  features <- unique(features)
  features_drop <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression")
    }
    status <- check_DataType(srt, slot = slot)
    message("Data type detected in ", slot, " slot: ", status)
    if (status == "raw_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "log_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), expm1) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "raw_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), log1p) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene > 0)) {
    dat_gene <- t(GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE])[colnames(srt), , drop = FALSE]
  } else {
    dat_gene <- matrix(nrow = ncol(srt), ncol = 0)
  }
  if (length(features_meta > 0)) {
    dat_meta <- srt@meta.data[, features_meta, drop = FALSE]
  } else {
    dat_meta <- matrix(nrow = ncol(srt), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])


  if (!all(sapply(dat_exp, is.numeric))) {
    stop("'features' must be type of numeric variable.")
  }
  dat_use <- srt@meta.data[, unique(c(group.by, bg.by, split.by)), drop = FALSE]
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }

  nlev <- sapply(dat_use, nlevels)
  nlev <- nlev[nlev > 50]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 50 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }
  if (isTRUE(same.y.lims) && is.null(y.max)) {
    y.max <- max(unlist(dat_exp[, features])[is.finite(unlist(dat_exp[, features]))])
  }
  if ((isTRUE(stack) || split.by == group.by) && (fill.by == "feature")) {
    colors <- palette_scp(features, palette = palette, palcolor = palcolor)
  } else {
    colors <- palette_scp(levels(dat_use[[split.by]]), palette = palette, palcolor = palcolor)
  }
  if (!is.null(bg.by)) {
    bg_color <- palette_scp(levels(dat_use[[bg.by]]), palette = bg_palette, palcolor = bg_palcolor)
  }

  plist <- list()
  xlab <- xlab %||% group.by
  ylab <- ylab %||% "Expression level"
  for (f in features) {
    dat <- cbind(dat_use, dat_exp[row.names(dat_use), f, drop = FALSE])
    dat[, "cell"] <- rownames(dat)
    dat[, "features"] <- f
    dat[, "value"] <- dat[, f]
    dat[, "group.by"] <- dat[, group.by]
    dat[, "split.by"] <- dat[, split.by]
    dat[, "bg.by"] <- dat[, bg.by]
    if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
      df_sort <- aggregate(dat[, "value", drop = FALSE], by = list(dat[["group.by"]]), mean)
      if (is.character(sort) && sort == "increasing") {
        decreasing <- FALSE
      } else {
        decreasing <- TRUE
      }
      sortlevel <- as.character(df_sort[order(df_sort[["value"]], decreasing = decreasing), 1])
      dat[, "group.by"] <- factor(dat[, "group.by"], levels = sortlevel)
    }
    if (fill.by == "feature") {
      dat[, "split.by"] <- f
      keynm <- "Features"
    } else {
      keynm <- split.by
    }

    comb <- expand.grid(x = levels(dat[["split.by"]]), y = levels(dat[["group.by"]]))
    dat[, "group"] <- factor(paste("a", dat[["split.by"]], "b", dat[["group.by"]], sep = "-"), levels = paste("a", comb[[1]], "b", comb[[2]], sep = "-"))

    y_max_use <- y.max %||% max(dat[, "value"][is.finite(x = dat[, "value"])])
    y_min_use <- min(dat[, "value"][is.finite(x = dat[, "value"])])

    levels_order <- levels(dat[["group.by"]])
    if (isTRUE(flip)) {
      dat[["group.by"]] <- factor(dat[["group.by"]], levels = rev(levels(dat[["group.by"]])))
      aspect.ratio <- 1 / aspect.ratio
      if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
        aspect.ratio <- NULL
      }
    }
    p <- ggplot(dat, aes(
      x = .data[["group.by"]], y = .data[["value"]], fill = .data[["split.by"]]
    )) +
      geom_blank()
    if (!is.null(bg.by)) {
      bg_df <- unique(dat[, c("group.by", "bg.by")])
      bg_list <- list()
      for (i in seq_len(nrow(bg_df))) {
        x <- as.numeric(bg_df[i, "group.by"])
        bg_list[[i]] <- annotate(
          geom = "rect",
          xmin = ifelse(x == 1, -Inf, x - 0.5),
          xmax = ifelse(x == nlevels(bg_df[, "group.by"]), Inf, x + 0.5),
          ymin = -Inf,
          ymax = Inf,
          fill = bg_color[as.character(bg_df[i, "bg.by"])], alpha = bg_apha
        )
      }
      p <- p + bg_list
    }
    if (isTRUE(split.plot)) {
      p <- p + geom_split_violin(scale = "width", adjust = adjust, trim = TRUE, alpha = alpha)
    } else {
      p <- p + geom_violin(scale = "width", adjust = adjust, trim = TRUE, alpha = alpha)
    }
    if (length(comparisons) > 0) {
      check_R("ggpubr")
      p <- p + ggpubr::stat_compare_means(
        symnum.args = list(
          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 2),
          symbols = c("****", "***", "**", "*", "ns")
        ),
        size = 3.5,
        step.increase = 0.1,
        tip.length = 0.03,
        comparisons = comparisons, ref.group = ref_group, method = pairwise_method
      )
      y_max_use <- layer_scales(p)$y$range$range[2]
    }
    if (isTRUE(multiplegroup_comparisons)) {
      check_R("ggpubr")
      p <- p + ggpubr::stat_compare_means(method = multiple_method, label.y = Inf, vjust = 1.3, hjust = 0, size = 3.5)
      y_max_use <- y_min_use + (y_max_use - y_min_use) * 1.2
    }

    if (isTRUE(add_point)) {
      p <- p + geom_point(
        color = pt.color,
        size = pt.size, alpha = pt.alpha,
        position = position_jitterdodge(jitter.width = jitter.width, dodge.width = 0.9, seed = 11), show.legend = FALSE
      )
      if (!is.null(cells.highlight)) {
        p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
        cell_df <- subset(p$data, cells.highlight == TRUE)
        if (nrow(cell_df) > 0) {
          p <- p + geom_point(
            data = cell_df,
            color = cols.highlight, size = sizes.highlight, alpha = alpha.highlight,
            position = position_jitterdodge(jitter.width = jitter.width, dodge.width = 0.9, seed = 11), show.legend = FALSE
          )
        }
      }
    }
    if (isTRUE(add_box)) {
      p <- p + geom_boxplot(aes(group = .data[["group"]]),
        position = position_dodge2(width = 0.9), color = "black", fill = box_fill, width = box_width, show.legend = FALSE, outlier.shape = NA
      ) +
        stat_summary(
          fun = median, geom = "point", mapping = aes(group = .data[["split.by"]]),
          position = position_dodge2(width = 0.9), color = "black", fill = "white", size = 1.5, shape = 21,
        )
    }
    if (isTRUE(stack) && !isTRUE(flip)) {
      p <- p + facet_grid(features ~ .) + theme(strip.text.y = element_text(angle = 0))
    } else {
      p <- p + facet_grid(. ~ features)
    }

    p <- p + labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
      scale_x_discrete(drop = !keep_empty) +
      scale_y_continuous(limits = c(y_min_use, y_max_use), trans = y.trans, n.breaks = y.nbreaks) +
      scale_fill_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, drop = FALSE) +
      scale_color_manual(name = paste0(keynm, ":"), values = colors, breaks = levels_order, drop = FALSE) +
      do.call(theme_use, list(
        aspect.ratio = aspect.ratio,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.grid.major.y = element_line(color = "grey", linetype = 2),
        legend.position = legend.position,
        legend.direction = legend.direction
      )) + guides(fill = guide_legend(
        title.hjust = 0,
        keywidth = 0.05,
        keyheight = 0.05,
        default.unit = "inch",
        order = 1,
        override.aes = list(size = 4.5, color = "black", alpha = 1)
      ))

    if (isTRUE(flip)) {
      p <- p + coord_flip()
    }
    plist[[f]] <- p
  }

  if (isTRUE(stack)) {
    legend <- cowplot::get_plot_component(plist[[1]], pattern = "guide-box")
    if (isTRUE(flip)) {
      lab <- textGrob(label = ifelse(is.null(ylab), "Expression level", ylab), hjust = 0.5)
      plist <- lapply(seq_along(plist), FUN = function(i) {
        p <- plist[[i]]
        if (i != 1) {
          p <- p + theme(
            legend.position = "none",
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(
              t = 0.2, r = -0.06, b = 0.2, l = -0.06,
              unit = "cm"
            )
          )
        } else {
          p <- p + theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            plot.margin = margin(
              t = 0.2, r = -0.06, b = 0.2, l = 0.2,
              unit = "cm"
            )
          )
        }
        return(ggplotGrob(p))
      })
      grob <- do.call(cbind, plist)
      grob <- gtable_add_rows(grob, grobHeight(lab), -1)
      grob <- gtable_add_grob(grob, lab, t = dim(grob)[1], l = mean(grob$layout[grep("panel", grob$layout$name), "l"]), clip = FALSE)

      grob <- gtable_add_cols(grob, sum(legend$widths), -1)
      grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])

      plot <- plot_grid(grob)
      return(plot)
    } else {
      lab <- textGrob(label = ifelse(is.null(ylab), "Expression level", ylab), rot = 90, hjust = 0.5)
      plist <- lapply(seq_along(plist), FUN = function(i) {
        p <- plist[[i]]
        if (i != length(plist)) {
          p <- p + theme(
            legend.position = "none",
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(
              t = -0.06, r = 0.2, b = -0.06, l = 0.2,
              unit = "cm"
            )
          )
          if (i == 1) {
            p <- p + theme(plot.title = element_blank(), plot.subtitle = element_blank())
          }
        } else {
          p <- p + theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            plot.margin = margin(
              t = -0.06, r = 0.2, b = 0.2, l = 0.2,
              unit = "cm"
            )
          )
        }
        return(ggplotGrob(p))
      })
      grob <- do.call(rbind, plist)
      grob <- gtable_add_cols(grob, grobWidth(lab), 0)
      grob <- gtable_add_grob(grob, lab, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1, clip = FALSE)

      grob <- gtable_add_cols(grob, sum(legend$widths), -1)
      grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])

      grob <- gtable_add_padding(grob, unit(c(0.2, 0, 0, 0), units = "cm"))
      plot <- plot_grid(grob)
    }
    return(plot)
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

#' ExpDotPlot
#'
#' @param srt A \code{Seurat} object.
#' @param features A vector of gene names to plot.
#' @param feature_split A vector of group names for features.
#' @param cell_split_by Columns used to calculate cell expression. One heatmap per column name.
#' @param exp_method Method used to calculate cell expression.
#' @param assay Assay used to calculate the expression.
#' @param heatmap_palette Heatmap expression palette.
#' @param feature_palette Feature groups palette.
#' @param cell_palette Column palette.
#' @param grid_size size for each dot.
#' @param aggregate_fun
#' @param slot
#' @param lib_normalize
#' @param libsize
#' @param add_reticle
#' @param cluster_features
#' @param cluster_cells
#'
#' @examples
#' library(dplyr)
#' data(pancreas_sub)
#' ExpDotPlot(pancreas_sub,
#'   features = c(
#'     "Sox9", "Anxa2", "Bicc1", # Ductal
#'     "Neurog3", "Hes6", # EPs
#'     "Fev", "Neurod1", # Pre-endocrine
#'     "Rbp4", "Pyy", # Endocrine
#'     "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#'   ),
#'   cell_split_by = c("CellType", "SubCellType")
#' )
#'
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' de_filter <- filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' de_top <- de_filter %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ExpDotPlot(pancreas_sub, features = de_top$gene, feature_split = de_top$group1, cell_split_by = "CellType")
#' @importFrom circlize colorRamp2
#' @importFrom stats aggregate formula quantile
#' @importFrom ComplexHeatmap Legend HeatmapAnnotation anno_block anno_simple Heatmap draw
#' @importFrom grid gpar grid.grabExpr grid.rect grid.circle
#' @importFrom cowplot plot_grid
#' @importFrom scales alpha
#' @importFrom methods getFunction
#' @export
ExpDotPlot <- function(srt, features = NULL, feature_split = NULL, cell_split_by = NULL, aggregate_fun = mean,
                       slot = "counts", assay = "RNA", exp_method = c("zscore", "raw", "log2fc", "log1p"),
                       lib_normalize = TRUE, libsize = NULL,
                       heatmap_palette = "YlOrRd", feature_palette = "Paired", cell_palette = "Paired",
                       add_reticle = FALSE,
                       cluster_features = FALSE, cluster_cells = FALSE, grid_size = 0.4) {
  exp_method <- match.arg(exp_method)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  exp_name <- switch(exp_method,
    "raw" = data_nm,
    "zscore" = paste0("Z-score(", data_nm, ")"),
    "log2fc" = paste0("Log2(", data_nm, "FC)"),
    "log1p" = paste0("Log(", data_nm, "+1)")
  )

  if (is.null(features)) {
    features <- VariableFeatures(srt)
  }
  if (is.null(cell_split_by)) {
    stop("'cell_split_by' must be provided.")
  }
  cell_split_by <- cell_split_by[cell_split_by %in% colnames(srt@meta.data)]
  if (length(cell_split_by) == 0) {
    stop("Stop plot! 'cell_split_by' is invalid!")
  }
  columns_length <- lapply(srt[[cell_split_by, drop = FALSE]], function(x) length(unique(x))) %>% unlist()
  if (any(columns_length < 2) && exp_method == "zscore") {
    stop(paste0("'cell_split_by' ", cell_split_by[columns_length < 2], " has only one group."))
  }
  if (!is.null(feature_split)) {
    if (length(feature_split) != length(features)) {
      stop("length(feature_split)!=length(features)")
    }
  }

  index <- features %in% rownames(srt)
  features <- features[index]
  feature_split <- feature_split[index]
  if (length(features) > 500) {
    stop("Too many features, suggest reducing the number of features.")
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }

  dat <- cbind.data.frame(srt@meta.data[, cell_split_by, drop = FALSE], t(GetAssayData(srt, slot = slot, assay = assay)[unique(features), ]))
  if (isTRUE(lib_normalize) && min(dat[, features], na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- colSums(GetAssayData(srt, slot = "counts", assay = assay)[, rownames(dat)])
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning("Values in 'counts' slot is non-integer. Set the libsize to 1.", immediate. = TRUE)
      }
    }
    dat[, features] <- dat[, features] / libsize_use * median(libsize_use)
  }
  dat <- dat[, c(cell_split_by, features)]
  dotHT_list <- list()

  for (i in cell_split_by) {
    mat <- t(aggregate(formula(paste0(".~", i)), dat[, c(i, features)], FUN = aggregate_fun))
    colnames(mat) <- mat[i, ]
    mat <- mat[-which(rownames(mat) == i), ]
    mat <- apply(mat, c(1, 2), as.numeric)

    if (exp_method == "raw") {
      mat <- mat
    } else if (exp_method == "zscore") {
      mat <- t(scale(t(mat)))
    } else if (exp_method == "log2fc") {
      mat <- log2(mat / rowMeans(mat))
    } else if (exp_method == "log1p") {
      mat <- log1p(mat)
    }
    mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)])) * ifelse(mat[is.infinite(mat)] > 0, 1, -1)
    mat[is.na(mat)] <- mean(mat, na.rm = TRUE)

    color_palette <- palette_scp(palette = heatmap_palette)
    colors <- colorRamp2(seq(quantile(mat, 0.01, na.rm = TRUE), quantile(mat, 0.99, na.rm = TRUE), length = 100), color_palette)
    legend_color <- Legend(
      col_fun = colors,
      labels_gp = gpar(fontsize = 12),
      title = exp_name, title_position = "topleft",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      type = "grid",
      border = TRUE,
      background = "transparent",
      direction = "horizontal"
    )
    legend_point <- Legend(
      labels = paste0(seq(20, 100, length.out = 5), "%"), labels_gp = gpar(fontsize = 12),
      title = "Percent", title_position = "topleft",
      title_gp = gpar(fontsize = 12, fontfamily = "sans"),
      type = "points", pch = 21,
      size = unit(pi * grid_size^2 * seq(0.2, 1, length.out = 5), "cm"),
      grid_height = unit(grid_size, "cm"),
      grid_width = unit(grid_size, "cm"),
      legend_gp = gpar(fill = "grey80"), border = FALSE,
      background = "transparent",
      direction = "vertical", nrow = 1
    )

    dat_pec <- aggregate(formula(paste0(".~", i)), dat, FUN = function(x) {
      sum(x > 0) / length(x)
    }) %>% t()
    colnames(dat_pec) <- dat_pec[i, ]
    dat_pec <- dat_pec[-which(rownames(dat_pec) == i), ]
    dat_pec <- apply(dat_pec, c(1, 2), as.numeric)
    dat_pec <- dat_pec[rownames(mat), ]
    assign(paste0(i, ".dat_pec"), dat_pec)

    if (!is.null(feature_split)) {
      if (!is.factor(feature_split)) {
        feature_split <- factor(feature_split, levels = unique(feature_split))
      }
      row_color <- palette_scp(levels(feature_split), palette = feature_palette, matched = TRUE)
      row_split <- feature_split
      names(row_split) <- row_split
      df_left_annotation <- HeatmapAnnotation(
        foo1 = anno_block(
          gp = gpar(fill = row_color)
        ),
        width = unit(0.5, "cm"), which = "row"
      )
    } else {
      df_left_annotation <- row_split <- NULL
    }
    df_top_annotation <- HeatmapAnnotation(
      Cell = anno_simple(
        x = colnames(mat),
        col = palette_scp(colnames(mat), palette = cell_palette),
        gp = gpar(col = "black"),
      ),
      height = unit(0.5, "cm"), which = "column", show_annotation_name = FALSE
    )
    use_raster <- length(features) > 2000
    funbody <- paste0(
      "
      pec <- ", i, ".dat_pec", "[i, j];
      grid.rect(x, y,
                width = w, height = h,
                gp = gpar(col = 'white', lwd = 1, fill = 'white')
      );
      grid.rect(x, y,
                width = w, height = h,
                gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
      );", if (isTRUE(add_reticle)) {
        v <- "vcolor=palette_scp(colnames(mat), palette = cell_palette)[j];
      grid.lines(x = c(x,x),y = c(y-h/2,y+h/2),
                gp = gpar(col = vcolor, lwd = 1.5)
      );"
        if (!is.null(feature_split)) {
          h <- "hcolor=palette_scp(levels(feature_split), palette = feature_palette)[feature_split[i]];
      grid.lines(x = c(x-w/2,x+w/2),y = c(y,y),
                gp = gpar(col = hcolor, lwd = 1.5)
      );"
        }
        paste0(v, h)
      },
      "grid.circle(x, y,
                  r = h * pec / 2,
                  gp = gpar(col = 'black', lwd = 1, fill = fill)
      );
    "
    )
    funbody <- gsub(pattern = "\n", replacement = "", x = funbody)

    eval(parse(text = paste("cell_fun <- function(j, i, x, y, w, h, fill) {", funbody, "}", sep = "")), envir = environment())

    dotHT_list[[i]] <- Heatmap(mat,
      col = colors,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10),
      column_names_side = "top",
      column_names_rot = 90,
      cluster_columns = cluster_cells,
      cluster_rows = cluster_features,
      row_split = row_split,
      column_split = NULL,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      row_title_rot = 0,
      cell_fun = getFunction("cell_fun", where = environment()),
      left_annotation = df_left_annotation,
      top_annotation = df_top_annotation,
      width = unit(ncol(mat) * grid_size, "cm"),
      height = unit(length(features) * grid_size, "cm"),
      show_heatmap_legend = FALSE,
      use_raster = use_raster,
      raster_device = "png"
    )
  }
  ht_list <- NULL
  for (ht in dotHT_list) {
    ht_list <- ht_list + ht
  }
  gTree <- grid.grabExpr({
    draw(ht_list,
      heatmap_legend_list = list(legend_color, legend_point), heatmap_legend_side = "bottom",
      padding = unit(c(1, max(c(nchar(features) * 0.15, 1)), max(c(nchar(as.character(srt[[cell_split_by, drop = TRUE]])) * 0.15, 1)), 1), "cm")
    ) # bottom, left, top and right
  })
  p <- plot_grid(gTree)

  return(p)
}

#' ExpHeatmap
#'
#' @param srt
#' @param features
#' @param feature_split
#' @param cluster_features
#' @param cell_split_by
#' @param cluster_cells
#' @param max_cells
#' @param slot
#' @param assay
#' @param exp_method
#' @param n_split
#' @param heatmap_split_by
#' @param split_method
#' @param row_title_size
#' @param decreasing
#' @param lib_normalize
#' @param libsize
#' @param anno_keys
#' @param anno_features
#' @param IDtype
#' @param species
#' @param db_IDtype
#' @param db_update
#' @param db_version
#' @param Ensembl_version
#' @param mirror
#' @param enrichment
#' @param TERM2GENE
#' @param TERM2NAME
#' @param minGSSize
#' @param maxGSSize
#' @param universe
#' @param GO_simplify
#' @param GO_simplify_padjustCutoff
#' @param simplify_method
#' @param simplify_similarityCutoff
#' @param pvalueCutoff
#' @param padjustCutoff
#' @param topWord
#' @param min_word_length
#' @param exclude_words
#' @param anno_width
#' @param anno_size
#' @param nlabel
#' @param features_label
#' @param label_size
#' @param label_color
#' @param heatmap_palette
#' @param cell_split_palette
#' @param feature_split_palette
#' @param cell_annotation
#' @param cell_palette
#' @param cell_palcolor
#' @param feature_annotation
#' @param feature_palette
#' @param feature_palcolor
#' @param use_raster
#' @param height
#' @param width
#' @param units
#' @param seed
#'
#' @examples
#' library(dplyr)
#' data(pancreas_sub)
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' de_filter <- filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' ht1 <- ExpHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, cell_split_by = "CellType"
#' )
#' ht1$plot
#'
#' ht2 <- ExpHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, n_split = 4, cell_split_by = "CellType"
#' )
#' ht2$plot
#'
#' ht3 <- ExpHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "CellType",
#'   species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht3$plot
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap Legend HeatmapAnnotation anno_empty anno_mark anno_simple anno_textbox draw decorate_heatmap_body width.HeatmapAnnotation height.HeatmapAnnotation
#' @importFrom grid gpar grid.grabExpr
#' @importFrom cowplot plot_grid
#' @importFrom dplyr "%>%"
#' @importFrom Seurat GetAssayData
#'
#' @export
ExpHeatmap <- function(srt, features = NULL, feature_split = NULL, cluster_features = FALSE,
                       cell_split_by = "orig.ident", cluster_cells = FALSE, max_cells = 100,
                       slot = "counts", assay = "RNA", exp_method = c("zscore", "raw", "log2fc", "log1p"),
                       n_split = NULL, heatmap_split_by = cell_split_by[1], split_method = c("kmeans", "hclust", "mfuzz"),
                       row_title_size = 12, decreasing = TRUE,
                       lib_normalize = TRUE, libsize = NULL,
                       anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                       IDtype = "symbol", species = "Homo_sapiens", db_IDtype = "symbol", db_update = FALSE, db_version = "latest", Ensembl_version = 103, mirror = NULL,
                       enrichment = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, universe = NULL,
                       GO_simplify = FALSE, GO_simplify_padjustCutoff = 0.2, simplify_method = "Rel", simplify_similarityCutoff = 0.7,
                       pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_id = TRUE, topWord = 20, min_word_length = 3, exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "system", "regulation", "positive", "negative", "response", "process"),
                       anno_width = unit(c(4, 2, 2), "in"), anno_size = c(6, 10),
                       nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                       heatmap_palette = "RdBu", cell_split_palette = "Paired", feature_split_palette = "jama",
                       cell_annotation = NULL, cell_palette = "Paired", cell_palcolor = NULL,
                       feature_annotation = NULL, feature_palette = "Dark2", feature_palcolor = NULL,
                       use_raster = NULL, height = NULL, width = NULL, units = "inch",
                       seed = 11) {
  set.seed(seed)
  exp_method <- match.arg(exp_method)
  split_method <- match.arg(split_method)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  exp_name <- switch(exp_method,
    "raw" = data_nm,
    "zscore" = paste0("Z-score(", data_nm, ")"),
    "log2fc" = paste0("Log2(", data_nm, "FC)"),
    "log1p" = paste0("Log(", data_nm, "+1)")
  )
  if (missing(srt)) {
    stop("srt must be provided.")
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features or zero")
  }

  if (length(cell_split_palette) == 1) {
    cell_split_palette <- rep(cell_split_palette, length(cell_split_by))
  }
  if (length(cell_split_palette) != length(cell_split_by)) {
    stop("'cell_split_palette' must be the same length as 'cell_split_by'")
  }
  if (any(!cell_split_by %in% colnames(srt@meta.data))) {
    stop("cell_split_by: ", cell_split_by[!cell_split_by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!heatmap_split_by %in% cell_split_by) {
    stop("'heatmap_split_by' must be a subset of 'cell_split_by'")
  }

  if (!is.null(cell_annotation)) {
    if (length(cell_palette) == 1) {
      cell_palette <- rep(cell_palette, length(cell_annotation))
    }
    if (length(cell_palcolor) == 1) {
      cell_palcolor <- rep(cell_palcolor, length(cell_annotation))
    }
    if (length(unique(length(cell_palette), length(cell_palcolor), length(cell_annotation))) != 1) {
      stop("cell_palette and cell_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% colnames(srt@meta.data))) {
      stop("cell_annotation: ", cell_annotation[!cell_annotation %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_palette) == 1) {
      feature_palette <- rep(feature_palette, length(feature_annotation))
    }
    if (length(feature_palcolor) == 1) {
      feature_palcolor <- rep(feature_palcolor, length(feature_annotation))
    }
    if (length(unique(length(feature_palette), length(feature_palcolor), length(feature_annotation))) != 1) {
      stop("feature_palette and feature_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt[[assay]]@meta.features))) {
      stop("feature_annotation: ", feature_annotation[!feature_annotation %in% colnames(srt[[assay]]@meta.features)], " is not in the meta data of the ", assay, " assay in the Seurat object.")
    }
  }

  if (is.null(features)) {
    features <- VariableFeatures(srt)
  }
  index <- features %in% rownames(srt)
  features <- features[index]
  features_unique <- make.unique(features)
  if (!is.null(feature_split)) {
    feature_split <- feature_split[index]
    names(feature_split) <- features_unique
  }

  cell_groups <- list()
  for (cell_split_var in cell_split_by) {
    if (!is.factor(srt[[cell_split_var, drop = TRUE]])) {
      srt[[cell_split_var, drop = TRUE]] <- factor(srt[[cell_split_var, drop = TRUE]], levels = unique(srt[[cell_split_var, drop = TRUE]]))
    }
    cell_groups[[cell_split_var]] <- lapply(levels(srt[[cell_split_var, drop = TRUE]]), function(x) {
      cells <- colnames(srt)[srt[[cell_split_var, drop = TRUE]] == x]
      size <- ifelse(length(cells) > max_cells, max_cells, length(cells))
      cells_sample <- sample(cells, size)
      out <- setNames(rep(x, size), cells_sample)
      return(out)
    }) %>% unlist(use.names = TRUE)
    cell_groups[[cell_split_var]] <- factor(cell_groups[[cell_split_var]], levels = levels(srt[[cell_split_var, drop = TRUE]]))
  }

  mat_all <- as.matrix(GetAssayData(srt, slot = slot, assay = assay)[features, unique(unlist(lapply(cell_groups, names)))])
  rownames(mat_all) <- features_unique
  if (isTRUE(lib_normalize) && min(mat_all, na.rm = TRUE) >= 0) {
    if (!is.null(libsize)) {
      libsize_use <- libsize
    } else {
      libsize_use <- colSums(GetAssayData(srt, slot = "counts", assay = assay)[, colnames(mat_all)])
      isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
      if (isTRUE(isfloat)) {
        libsize_use <- rep(1, length(libsize_use))
        warning("Values in 'counts' slot is non-integer. Set the libsize to 1.", immediate. = TRUE)
      }
    }
    mat_all <- t(t(mat_all) / libsize_use * median(libsize_use))
  }

  mat_list <- list()
  for (cell_split_var in cell_split_by) {
    mat_tmp <- mat_all[, names(cell_groups[[cell_split_var]])]
    if (exp_method == "raw") {
      mat_tmp <- mat_tmp
    } else if (exp_method == "zscore") {
      mat_tmp <- t(scale(t(mat_tmp)))
    } else if (exp_method == "log2fc") {
      mat_tmp <- log2(mat_tmp / rowMeans(mat_tmp))
    } else if (exp_method == "log1p") {
      mat_tmp <- log1p(mat_tmp)
    }
    mat_tmp[is.infinite(mat_tmp)] <- max(abs(mat_tmp[!is.infinite(mat_tmp)])) * ifelse(mat_tmp[is.infinite(mat_tmp)] > 0, 1, -1)
    mat_tmp[is.na(mat_tmp)] <- mean(mat_tmp, na.rm = TRUE)
    mat_list[[cell_split_var]] <- mat_tmp
  }

  if (length(heatmap_split_by) > 1 && is.null(feature_split) && n_split < nrow(mat_all) && n_split > 1) {
    mat_split <- mat_all[, unique(unlist(lapply(cell_groups[heatmap_split_by], names)))]
    colnames(mat_split) <- make.unique(mat_split)
    if (exp_method == "raw") {
      mat_split <- mat_split
    } else if (exp_method == "zscore") {
      mat_split <- t(scale(t(mat_split)))
    } else if (exp_method == "log2fc") {
      mat_split <- log2(mat_split / rowMeans(mat_split))
    } else if (exp_method == "log1p") {
      mat_split <- log1p(mat_split)
    }
    mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
    mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  } else {
    mat_split <- mat_list[[heatmap_split_by]]
  }

  if (exp_method %in% c("zscore", "log2fc")) {
    b <- ceiling(min(abs(quantile(mat_split, c(0.01, 0.99), na.rm = TRUE))) * 2) / 2
    colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette))
  } else if (exp_method %in% c("raw", "log1p")) {
    b <- quantile(mat_split, c(0.01, 0.99), na.rm = TRUE)
    colors <- colorRamp2(seq(b[1], b[2], length = 100), palette_scp(palette = heatmap_palette))
  }

  cell_metadata <- cbind.data.frame(
    data.frame(row.names = colnames(mat_all), cells = colnames(mat_all)),
    srt@meta.data[colnames(mat_all), c(cell_split_by, cell_annotation), drop = FALSE]
  )
  feature_metadata <- cbind.data.frame(
    data.frame(row.names = features_unique, features = features, features_uique = features_unique),
    srt[[assay]]@meta.features[features, c(feature_annotation), drop = FALSE]
  )
  feature_metadata[, "duplicated"] <- feature_metadata[["features"]] %in% features[duplicated(features)]

  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_top_list <- list()
  for (i in seq_along(cell_split_by)) {
    cell_split_var <- cell_split_by[i]
    palette <- cell_split_palette[i]
    ha_top_list[[cell_split_var]] <- HeatmapAnnotation(
      cell_split = anno_block(
        gp = gpar(fill = palette_scp(levels(cell_groups[[cell_split_var]]), palette = palette)),
        show_name = FALSE
      ),
      which = "column", border = TRUE
    )
    lgd[[cell_split_var]] <- Legend(
      title = cell_split_var, labels = levels(cell_groups[[cell_split_var]]),
      legend_gp = gpar(fill = palette_scp(levels(cell_groups[[cell_split_var]]), palette = palette)), border = TRUE
    )
  }

  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_palette[i]
      palcolor <- cell_palcolor[[i]]
      cell_anno <- cell_metadata[, cellan]
      names(cell_anno) <- rownames(cell_metadata)
      if (!is.factor(cell_anno)) {
        cell_anno <- factor(cell_anno, levels = unique(cell_anno))
      }
      for (cell_split_var in cell_split_by) {
        ha_top_tmp <- list()
        ha_top_tmp[[cellan]] <- anno_simple(
          x = as.character(cell_anno[names(cell_groups[[cell_split_var]])]),
          col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
          border = TRUE
        )
        ha_top <- do.call("HeatmapAnnotation", args = c(ha_top_tmp, which = "column", show_annotation_name = FALSE))
        ha_top_list[[cell_split_var]] <- c(ha_top_list[[cell_split_var]], ha_top)
      }
      lgd[[cellan]] <- Legend(
        title = cellan, labels = levels(cell_anno),
        legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
      )
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat_split) <= n_split)) {
      feature_split <- row_split <- NULL
    } else {
      if (n_split == 1) {
        feature_split <- row_split <- setNames(rep(1, nrow(mat_split)), rownames(mat_split))
      } else {
        if (split_method == "mfuzz") {
          check_R("Mfuzz")
          require("Mfuzz")
          eset <- new("ExpressionSet", exprs = mat_split)
          eset <- Mfuzz::standardise(eset)
          min_fuzzification <- Mfuzz::mestimate(eset)
          if (is.null(fuzzification)) {
            fuzzification <- min_fuzzification + 0.1
          } else {
            if (fuzzification <= min_fuzzification) {
              warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2))
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
        if (split_method == "kmeans") {
          km <- kmeans(mat_split, centers = n_split, iter.max = 1e4, nstart = 20)
          row_split <- feature_split <- km$cluster
        }
        if (split_method == "hclust") {
          hc <- hclust(dist(mat_split))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      groupmean <- aggregate(t(mat_split), by = list(unlist(cell_groups[heatmap_split_by])), mean)
      maxgroup <- groupmean[, 1][apply(groupmean[, names(row_split)], 2, which.max)]

      df <- data.frame(row_split = row_split, order_by = as.numeric(maxgroup))
      df_order <- aggregate(df, by = list(row_split), FUN = mean)
      df_order <- df_order[order(df_order[["order_by"]], decreasing = !decreasing), ]
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))
    }
  } else {
    row_split <- feature_split <- feature_split[row.names(mat_split)]
  }

  if (isTRUE(cluster_features)) {
    clusters <- row_split %||% setNames(rep(1, nrow(mat_split)), rownames(mat_split))
    dend_list <- lapply(levels(factor(clusters)), function(x) {
      feat <- names(clusters)[clusters == x]
      row_dend <- as.dendrogram(hclust(dist(mat_split[feat, ])))
      return(row_dend)
    })
    order_list <- split(clusters, clusters)
    reorder <- -rowMeans(mat_split, na.rm = TRUE)
    for (i in seq_along(dend_list)) {
      sub_ind <- names(order_list[[i]])
      dend_list[[i]] <- reorder(dend_list[[i]], reorder[sub_ind], mean)
      order_list[[i]] <- sub_ind[order.dendrogram(dend_list[[i]])]
    }
    features_ordered <- unlist(order_list)
  } else {
    if (!is.null(row_split)) {
      features_ordered <- names(row_split)[order(row_split)]
    } else {
      features_ordered <- rownames(mat_split)
    }
  }

  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split
    feature_metadata[["index"]] <- order(feature_metadata[["feature_split"]])
  } else {
    feature_metadata[["feature_split"]] <- NA
  }

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(features) > nlabel) {
        index <- seq(floor(length(features_ordered) / nlabel) - 1, length(features_ordered), ceiling(length(features_ordered) / nlabel))
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
      warning(paste0(paste0(drop, collapse = ","), "was not found in the dynamic features"))
    }
  }

  ha_left <- NULL
  if (length(index) > 0) {
    ha_left <- HeatmapAnnotation(
      gene = anno_mark(
        at = which(rownames(feature_metadata) %in% features_ordered[index]),
        labels = feature_metadata[which(rownames(feature_metadata) %in% features_ordered[index]), "features"],
        side = "left",
        labels_gp = gpar(fontsize = label_size, col = label_color),
        link_gp = gpar(fontsize = label_size, col = label_color)
      ),
      which = "row", show_annotation_name = FALSE
    )
  }
  if (!is.null(row_split)) {
    ha_clusters <- HeatmapAnnotation(
      feat_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
      feat_cluster = anno_block(
        gp = gpar(fill = palette_scp(row_split, type = "discrete", palette = feature_split_palette)),
        width = unit(0.1, "in")
      ),
      which = "row", show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(
      title = "Cluster", labels = levels(factor(row_split)),
      legend_gp = gpar(fill = palette_scp(row_split, type = "discrete", palette = feature_split_palette)), border = TRUE
    )
  }

  ha_right <- NULL
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featurean <- feature_annotation[i]
      palette <- feature_palette[i]
      palcolor <- feature_palcolor[[i]]
      feature_class <- feature_metadata[features_unique, featurean]
      if (!is.factor(feature_class)) {
        feature_class <- factor(feature_class, levels = unique(feature_class))
      }
      ha_feature <- list()
      ha_feature[[featurean]] <- anno_simple(
        x = as.character(feature_class),
        col = palette_scp(feature_class, palette = palette, palcolor = palcolor),
        width = unit(0.5, "cm"), which = "row"
      )
      ha_feature <- do.call("HeatmapAnnotation", args = c(ha_feature, which = "row", annotation_name_side = "top", border = TRUE))
      if (is.null(ha_right)) {
        ha_right <- ha_feature
      } else {
        ha_right <- c(ha_right, ha_feature)
      }
      lgd[[featurean]] <- Legend(
        title = featurean, labels = levels(feature_class),
        legend_gp = gpar(fill = palette_scp(feature_class, palette = palette, palcolor = palcolor)), border = TRUE
      )
    }
  }

  res <- NULL
  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    check_R("ggwordcloud")
    check_R("simplifyEnrichment")
    geneID <- feature_metadata[features_unique, "features"]
    geneID_groups <- feature_metadata[features_unique, "feature_split"]
    if (all(is.na(geneID_groups))) {
      geneID_groups <- rep(1, length(geneID))
    }
    res <- RunEnrichment(
      geneID = geneID, geneID_groups = geneID_groups, IDtype = IDtype, species = species,
      db_IDtype = db_IDtype, db_update = db_update, db_version = db_version, Ensembl_version = Ensembl_version, mirror = mirror,
      enrichment = enrichment, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize, universe = universe,
      GO_simplify = GO_simplify, GO_simplify_padjustCutoff = GO_simplify_padjustCutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff
    )
    if (nrow(res$enrichment) == 0) {
      stop("No enrichment result found.")
    }
    metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
    pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
    padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

    df <- res$enrichment %>%
      filter(Enrichment %in% enrichment) %>%
      group_by(Enrichment, Groups) %>%
      filter(.data[["pvalue"]] <= pvalueCutoff & .data[["p.adjust"]] <= padjustCutoff) %>%
      arrange(desc(-.data[["pvalue"]])) %>%
      as.data.frame()
    df_list <- split.data.frame(df, ~ Enrichment + Groups)
    df_list <- df_list[lapply(df_list, nrow) > 0]

    for (enrich in enrichment) {
      nm <- strsplit(names(df_list), "\\.")
      subdf_list <- df_list[unlist(lapply(nm, function(x) x[[1]])) %in% enrich]

      ha_terms <- NULL
      if (isTRUE(anno_terms)) {
        ids_list <- lapply(subdf_list, function(df) {
          if (isTRUE(show_id)) {
            ids <- paste(head(df$ID, topTerm), head(df$Description, topTerm))
          } else {
            ids <- head(df$Description, topTerm)
            ids <- paste(toupper(substr(ids, 1, 1)), substr(ids, 2, nchar(ids)), sep = "")
          }
          df_out <- data.frame(keyword = ids)
          df_out[["col"]] <- palette_scp(-log10(head(df[, "p.adjust"], topTerm)), type = "continuous", palette = "Spectral", matched = TRUE)
          df_out[["col"]] <- sapply(df_out[["col"]], function(x) blendcolors(c(x, "black")))
          df_out[["fontsize"]] <- mean(anno_size)
          return(df_out)
        })
        names(ids_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_terms <- HeatmapAnnotation(
          id_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          id_cluster = anno_block(
            gp = gpar(fill = palette_scp(geneID_groups, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          id = anno_textbox(
            align_to = geneID_groups, text = ids_list, max_width = anno_width[1],
            word_wrap = TRUE, add_new_line = TRUE,
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
      }

      ha_keys <- NULL
      if (isTRUE(anno_keys)) {
        term_list <- lapply(subdf_list, function(df) {
          if (df$Enrichment[1] %in% c("GO_BP", "GO_CC", "GO_MF")) {
            df <- simplifyEnrichment::keyword_enrichment_from_GO(df[["ID"]]) %>%
              summarise(
                keyword = .data[["keyword"]],
                score = -(log10(.data[["padj"]])),
                count = .data[["n_term"]],
                Enrichment = df[["Enrichment"]][1],
                Groups = df[["Groups"]][1]
              ) %>%
              filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
              filter(!.data[["keyword"]] %in% exclude_words) %>%
              distinct() %>%
              mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
              as.data.frame()
            df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          } else {
            df <- df %>%
              mutate(keyword = strsplit(as.character(.data[["Description"]]), " ")) %>%
              unnest(cols = "keyword") %>%
              group_by(.data[["keyword"]], Enrichment, Groups) %>%
              summarise(
                keyword = .data[["keyword"]],
                score = sum(-(log10(.data[[metric]]))),
                count = n(),
                Enrichment = .data[["Enrichment"]],
                Groups = .data[["Groups"]],
                .groups = "keep"
              ) %>%
              filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
              filter(!.data[["keyword"]] %in% exclude_words) %>%
              distinct() %>%
              mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
              as.data.frame()
            df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          }
          df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
          df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
          df[["fontsize"]] <- rescale(df[, "count"], to = anno_size)
          return(df)
        })
        names(term_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_keys <- HeatmapAnnotation(
          terms_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          terms_cluster = anno_block(
            gp = gpar(fill = palette_scp(geneID_groups, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          terms = anno_textbox(
            align_to = geneID_groups, text = term_list, max_width = anno_width[2],
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
      }

      ha_features <- NULL
      if (isTRUE(anno_features)) {
        features_list <- lapply(subdf_list, function(df) {
          df <- df %>%
            mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
            unnest(cols = "keyword") %>%
            group_by(.data[["keyword"]], Enrichment, Groups) %>%
            summarise(
              keyword = .data[["keyword"]],
              score = sum(-(log10(.data[[metric]]))),
              count = n(),
              Enrichment = .data[["Enrichment"]],
              Groups = .data[["Groups"]],
              .groups = "keep"
            ) %>%
            filter(!.data[["keyword"]] %in% exclude_words) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
          df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
          df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
          df[["fontsize"]] <- rescale(df[, "count"], to = anno_size)
          return(df)
        })
        names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_features <- HeatmapAnnotation(
          feat_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          feat_cluster = anno_block(
            gp = gpar(fill = palette_scp(geneID_groups, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          feat = anno_textbox(
            align_to = geneID_groups, text = features_list, max_width = anno_width[3],
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
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

  if (is.null(use_raster)) {
    use_raster <- ifelse(max(sapply(mat_list, ncol)) * length(features) > 1e7, TRUE, FALSE)
  }
  ht_list <- c()
  for (cell_split_var in cell_split_by) {
    if (cell_split_var == cell_split_by[1]) {
      left_annotation <- ha_left
    } else {
      left_annotation <- NULL
    }
    if (cell_split_var == cell_split_by[length(cell_split_by)]) {
      right_annotation <- ha_right
    } else {
      right_annotation <- NULL
    }
    ht_list <- ht_list + Heatmap(
      name = cell_split_var,
      matrix = mat_list[[cell_split_var]],
      col = colors,
      column_title = cell_split_var,
      column_split = cell_groups[[cell_split_var]],
      row_split = row_split,
      row_title_gp = gpar(fontsize = row_title_size),
      cluster_rows = cluster_features,
      cluster_columns = cluster_cells,
      cluster_row_slices = FALSE,
      cluster_column_slices = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      top_annotation = ha_top_list[[cell_split_var]],
      left_annotation = left_annotation,
      right_annotation = right_annotation,
      show_heatmap_legend = FALSE,
      border = TRUE,
      use_raster = use_raster,
      raster_device = "png"
    )
  }

  if (length(index) == 0 && is.null(anno_keys) && is.null(anno_features) && is.null(width) && is.null(height)) {
    fix <- FALSE
  } else {
    fix <- TRUE
  }
  if (is.null(height)) {
    height <- max(convertHeight(unit(1, "npc"), units, valueOnly = TRUE), 7)
  }
  if (length(ha_top_list) > 0) {
    height_top <- c()
    for (ha_top in ha_top_list) {
      height_top <- max(height.HeatmapAnnotation(ha_top), height_top)
    }
    height <- as.numeric(convertUnit(unit(height, units = units) + height_top, units))
  }
  if (is.null(width)) {
    width <- max(convertWidth(unit(1, "npc"), units, valueOnly = TRUE), 7)
  }
  if (!is.null(ha_left)) {
    width <- as.numeric(convertUnit(unit(width, units = units) + width.HeatmapAnnotation(ha_left), units))
  }
  if (!is.null(ha_right)) {
    width <- as.numeric(convertUnit(unit(width, units = units) + width.HeatmapAnnotation(ha_right), units))
  }

  gTree <- grid.grabExpr(
    {
      draw(ht_list,
        annotation_legend_list = lgd,
        padding = unit(c(1, 1, 1, 1), "cm") # bottom, left, top and right
      )
    },
    height = height,
    width = width
  )
  if (!isTRUE(fix)) {
    p <- plot_grid(gTree)
  } else {
    p <- plot_grid(panel_fix_single(plot_grid(gTree), width = width, height = height, units = units))
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
#' @param align
#' @param axis
#' @param force
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' ExpCorPlot(pancreas_sub, features = c("Ghrl", "Gcg", "Ins1", "Ins2"), group.by = "SubCellType")
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom SeuratObject as.sparse
#' @importFrom dplyr group_by summarize "%>%" .data
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_density_2d stat_density_2d labs scale_x_continuous scale_y_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colourbar geom_hex stat_summary_hex
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom gtable gtable_add_cols
#' @importFrom cowplot plot_grid get_legend draw_grob
#' @importFrom Matrix t
#' @importFrom rlang %||%
#' @export
ExpCorPlot <- function(srt, features, group.by = NULL, split.by = NULL, slot = "data", assay = DefaultAssay(srt),
                       cor_method = "pearson", adjust = 1, margin = 1, reverse = FALSE,
                       add_equation = FALSE, add_r2 = TRUE, add_pvalue = TRUE, add_smooth = TRUE,
                       palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = NULL, pt.alpha = 1,
                       cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                       calculate_coexp = FALSE,
                       raster = NULL, raster.dpi = c(512, 512),
                       theme_use = "theme_scp", title = NULL, subtitle = NULL,
                       legend.position = "right", legend.direction = "vertical",
                       combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr", force = FALSE) {
  check_R("exaexa/scattermore")
  require(ggrepel)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("")
  }
  if (is.null(group.by)) {
    group.by <- "All_cells"
    srt[[group.by]] <- factor("All_cells")
  }
  for (i in c(split.by, group.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
  }
  if (!is.null(cells.highlight)) {
    if (!any(cells.highlight %in% colnames(srt))) {
      stop("No cells in 'cells.highlight' found in srt.")
    }
    if (!all(cells.highlight %in% colnames(srt))) {
      warning("Some cells in 'cells.highlight' not found in srt.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, colnames(srt))
  }

  features_drop <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt[[assay]]@counts)]
  features_meta <- features[features %in% colnames(srt@meta.data)]

  if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression")
    }
    status <- check_DataType(srt, slot = slot)
    message("Data type detected in ", slot, " slot: ", status)
    if (status == "raw_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "log_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), expm1) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else if (status == "raw_normalized_counts") {
      srt[["CoExp"]] <- GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE] %>%
        apply(c(1, 2), log1p) %>%
        apply(2, function(x) log1p(exp(mean(log(x)))))
    } else {
      stop("Can not determine the data type.")
    }
    features <- c(features, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene > 0)) {
    dat_gene <- t(GetAssayData(srt, slot = slot, assay = assay)[features_gene, , drop = FALSE])[colnames(srt), , drop = FALSE]
  } else {
    dat_gene <- matrix(nrow = ncol(srt), ncol = 0)
  }
  if (length(features_meta > 0)) {
    dat_meta <- srt@meta.data[, features_meta, drop = FALSE]
  } else {
    dat_meta <- matrix(nrow = ncol(srt), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  features <- unique(features[features %in% c(features_gene, features_meta)])
  if (length(features) < 2) {
    stop("features must be a vector of length at least 2.")
  }

  if (!inherits(dat_exp, "dgCMatrix")) {
    dat_exp <- as.sparse(as.matrix(dat_exp))
  }
  if (!all(sapply(dat_exp, is.numeric))) {
    stop("'features' must be type of numeric variable.")
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

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster <- raster %||% (nrow(dat_use) * ncol(combn(features, m = 2)) > 1e5)
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
    dat <- dat_use[dat_use[[split.by]] == s, ]
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
        do.call(theme_use, list(
          aspect.ratio = 1,
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(margin, margin, margin, margin),
          legend.position = "none"
        ))
      if (f1_index == f2_index) {
        p <- p + geom_violin(aes(x = .data[[group.by]], y = .data[[f1]], fill = .data[[group.by]]),
          scale = "width", adjust = adjust, trim = TRUE, na.rm = TRUE
        ) + scale_x_discrete(position = ifelse(isTRUE(reverse), "top", "bottom")) +
          scale_y_continuous(position = ifelse(isTRUE(reverse), "right", "left"))
      } else {
        p <- p + scale_x_continuous(
          n.breaks = 4, labels = scales::number_format(),
          limits = c(min(dat_exp[rownames(dat), ], na.rm = TRUE), max(dat_exp[rownames(dat), ], na.rm = TRUE)),
          position = ifelse(isTRUE(reverse), "top", "bottom")
        ) +
          scale_y_continuous(
            n.breaks = 4, labels = scales::number_format(),
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
          p$data[, "cells.highlight"] <- rownames(p$data) %in% cells.highlight
          cell_df <- subset(p$data, cells.highlight == TRUE)
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
        do.call(theme_use, list(
          legend.position = legend.position,
          legend.direction = legend.direction
        ))))
    } else {
      legend <- NULL
    }
    grob_row <- list()
    plotlist <- suppressWarnings(lapply(plotlist, ggplotGrob))
    for (i in seq(1, length(plotlist), length(features))) {
      grob_row[[paste0(i:(i + length(features) - 1), collapse = "-")]] <- do.call(cbind, plotlist[i:(i + length(features) - 1)])
    }
    grob <- do.call(rbind, grob_row)
    if (!is.null(legend)) {
      if (legend.position == "bottom") {
        grob <- gtable_add_rows(grob, sum(legend$heights), -1)
        grob <- gtable_add_grob(grob, legend, t = dim(grob)[1], l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
      }
      if (legend.position == "top") {
        grob <- gtable_add_rows(grob, sum(legend$heights), 0)
        grob <- gtable_add_grob(grob, legend, t = 1, l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
      }
      if (legend.position == "right") {
        grob <- gtable_add_cols(grob, sum(legend$widths), -1)
        grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])
      }
      if (legend.position == "left") {
        grob <- gtable_add_cols(grob, sum(legend$widths), 0)
        grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1)
      }
    }
    if (nlevels(dat_use[[split.by]]) > 1) {
      split_grob <- textGrob(s, just = "center", gp = gpar(fontface = "bold", fontsize = 13))
      grob <- gtable_add_rows(grob, grobHeight(split_grob), 0)
      grob <- gtable_add_grob(grob, split_grob,
        t = 1,
        l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"])
      )
    }
    if (!is.null(subtitle)) {
      subtitle_grob <- textGrob(subtitle, x = 0, hjust = 0, gp = gpar(fontface = "italic", fontsize = 13))
      grob <- gtable_add_rows(grob, grobHeight(subtitle_grob), 0)
      grob <- gtable_add_grob(grob, subtitle_grob,
        t = 1,
        l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"])
      )
    }
    if (!is.null(title)) {
      title_grob <- textGrob(title, x = 0, hjust = 0, gp = gpar(fontsize = 14))
      grob <- gtable_add_rows(grob, 2 * grobHeight(title_grob), 0)
      grob <- gtable_add_grob(grob, title_grob,
        t = 1,
        l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"])
      )
    }
    p <- plot_grid(grob)
    plist[[paste0(s)]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}


#' Statistical plot of cell classification
#'
#' @param srt
#' @param stat.by
#' @param group.by
#' @param split.by
#' @param cells_subset
#' @param keep_empty
#' @param stat_single
#' @param plot_type
#' @param stat_type
#' @param position
#' @param width
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
#' @param align
#' @param axis
#' @param force
#'
#' @examples
#' data("pancreas_sub")
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType")
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", position = "dodge", label = TRUE)
#'
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "rose", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "ring", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", plot_type = "pie", label = TRUE)
#'
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType")
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "rose", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "ring", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "area", label = TRUE)
#'
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", plot_type = "bar", label = TRUE, stat_single = TRUE)
#'
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "rose", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "ring", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "area", label = TRUE)
#'
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "rose", position = "dodge", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "ring", position = "dodge", label = TRUE)
#' ClassStatPlot(pancreas_sub, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "area", position = "dodge", label = TRUE, alpha = 0.8)
#'
#' ClassStatPlot(pancreas_sub, stat.by = c("CellType", "Phase"), plot_type = "chord")
#' ClassStatPlot(pancreas_sub, stat.by = c("CellType", "Phase"), plot_type = "sankey")
#'
#' ClassStatPlot(pancreas_sub,
#'   stat.by = c("CellType", "Phase"), plot_type = "venn",
#'   stat_level = list(CellType = c("Ductal", "Ngn3 low EP", "Ngn3 high EP"), Phase = "G2M")
#' )
#'
#' pancreas_sub$Ductal <- ifelse(pancreas_sub$CellType %in% c("Ductal", "Ngn3 low EP", "Ngn3 high EP"), TRUE, FALSE)
#' pancreas_sub$Endocrine <- ifelse(pancreas_sub$CellType %in% c("Pre-endocrine", "Endocrine"), TRUE, FALSE)
#' pancreas_sub$G1S <- ifelse(pancreas_sub$Phase %in% c("G1", "S"), TRUE, FALSE)
#' pancreas_sub$G2M <- ifelse(pancreas_sub$Phase %in% c("G2M"), TRUE, FALSE)
#' head(pancreas_sub@meta.data[, c("Ductal", "Endocrine", "G1S", "G2M")])
#' ClassStatPlot(pancreas_sub, stat.by = c("Ductal", "Endocrine", "G1S", "G2M"), plot_type = "venn")
#' ClassStatPlot(pancreas_sub, stat.by = c("Ductal", "Endocrine", "G1S", "G2M"), plot_type = "upset")
#'
#' @importFrom dplyr group_by across all_of mutate "%>%" .data summarise
#' @importFrom stats quantile xtabs
#' @importFrom ggplot2 ggplot ggplotGrob aes labs position_stack position_dodge2 scale_x_continuous scale_y_continuous geom_col geom_area geom_vline scale_fill_manual scale_fill_identity scale_color_identity scale_fill_gradientn guides guide_legend element_line coord_polar annotate geom_sf theme_void
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid get_legend draw_grob
#' @importFrom gtable gtable_add_rows gtable_add_cols gtable_add_grob
#' @importFrom rlang %||%
#' @export
ClassStatPlot <- function(srt, stat.by = "orig.ident", group.by = NULL, split.by = NULL,
                          cells_subset = NULL, keep_empty = FALSE, stat_single = FALSE, stat_level = NULL,
                          plot_type = c("bar", "rose", "ring", "pie", "area", "sankey", "chord", "venn", "upset"),
                          stat_type = c("percent", "count"), position = c("stack", "dodge"), width = 0.8,
                          palette = "Paired", palcolor = NULL, alpha = 1,
                          label = FALSE, label.size = 3.5, label.fg = "black", label.bg = "white", label.bg.r = 0.1,
                          theme_use = "theme_scp", aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                          legend.position = "right", legend.direction = "vertical",
                          combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr", force = FALSE) {
  stat_type <- match.arg(stat_type)
  plot_type <- match.arg(plot_type)
  position <- match.arg(position)

  if (is.null(stat.by)) {
    stop("'stat.by' must be provided.")
  }
  if (is.null(group.by)) {
    group.by <- "ClassStat_group"
    srt[["ClassStat_group"]] <- factor(paste0(stat.by, collapse = ","))
  }
  if (is.null(split.by)) {
    split.by <- "All_cells"
    srt[[split.by]] <- factor("")
  }

  for (i in unique(c(group.by, split.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
  }
  for (i in unique(stat.by)) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (plot_type %in% c("venn", "upset")) {
      if (!is.factor(srt[[i, drop = TRUE]]) && !is.logical(srt[[i, drop = TRUE]])) {
        srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
      }
    } else if (!is.factor(srt[[i, drop = TRUE]])) {
      srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
    }
  }

  if (length(stat.by) >= 2) {
    if (!plot_type %in% c("sankey", "chord", "venn", "upset")) {
      stop("plot_type must be one of 'sankey', 'chord', 'venn' and 'upset' whtn multiple 'stat.by' provided.")
    }
    if (length(stat.by) >= 3 && plot_type == "chord") {
      stop("'plot_type' is 'chord' and 'stat.by' can only be a vector of length 2.")
    }
    if (length(stat.by) >= 8 && plot_type == "venn") {
      stop("'plot_type' is 'venn' and 'stat.by' can only be a vector of length less than 8.")
    }
  }
  levels <- unique(unlist(lapply(srt@meta.data[, stat.by, drop = FALSE], function(x) {
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
        levels(srt[[stat, drop = TRUE]])[1]
      })
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
      if (!is.logical(srt[[i, drop = TRUE]])) {
        srt[[i, drop = TRUE]] <- srt[[i, drop = TRUE]] %in% stat_level[[i]]
      }
    }
  }

  if (any(group.by != "ClassStat_group") && plot_type %in% c("sankey", "chord", "venn", "upset")) {
    warning("group.by is not used when plot sankey, chord, venn or upset")
  }

  if (!is.null(cells_subset)) {
    if (!any(cells_subset %in% colnames(srt))) {
      stop("No cells in 'cells_subset' found in srt.")
    }
    if (!all(cells_subset %in% colnames(srt))) {
      warning("Some cells in 'cells_subset' not found in srt.", immediate. = TRUE)
    }
    cells_subset <- intersect(cells_subset, colnames(srt))
    srt <- srt[, cells_subset]
  }

  dat_all <- srt@meta.data[, unique(c(stat.by, group.by, split.by)), drop = FALSE]
  nlev <- sapply(dat_all, nlevels)
  nlev <- nlev[nlev > 50]
  if (length(nlev) > 0 && !isTRUE(force)) {
    warning(paste(names(nlev), sep = ","), " have more than 50 levels.", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  plist <- list()
  if (plot_type %in% c("bar", "rose", "ring", "pie", "area")) {
    xlab <- xlab %||% group.by
    ylab <- ylab %||% ifelse(stat_type == "count", "Count", "Percentage")
    colors <- palette_scp(levels(dat_all[[stat.by]]), palette = palette, palcolor = palcolor)
    for (g in group.by) {
      for (s in levels(dat_all[[split.by]])) {
        if (stat_type == "percent") {
          dat_use <- dat_all[dat_all[[split.by]] == s, ] %>%
            xtabs(formula = paste0("~", stat.by, "+", g)) %>%
            as.data.frame() %>%
            group_by(across(all_of(g)), .drop = FALSE) %>%
            mutate(groupn = sum(Freq)) %>%
            group_by(across(all_of(c(stat.by, g))), .drop = FALSE) %>%
            mutate(value = Freq / groupn) %>%
            as.data.frame()
        } else {
          dat_use <- dat_all[dat_all[[split.by]] == s, ] %>%
            xtabs(formula = paste0("~", stat.by, "+", g)) %>%
            as.data.frame() %>%
            mutate(value = Freq)
        }
        if (isTRUE(stat_single)) {
          groups <- levels(dat_use[[g]])
        } else {
          groups <- list(levels(dat_use[[g]]))
        }
        for (single_group in groups) {
          dat <- dat_use[dat_use[[g]] %in% single_group, ]
          dat[[g]] <- factor(dat[[g]], levels = levels(dat[[g]])[levels(dat[[g]]) %in% dat[[g]]])
          if (plot_type == "ring") {
            dat[[g]] <- factor(dat[[g]], levels = c("   ", levels(dat[[g]])))
            dat <- rbind(dat, dat[nrow(dat) + 1, ])
            dat[nrow(dat), g] <- "   "
          }
          if (position == "stack") {
            position_use <- position_stack(vjust = 0.5)
            scalex <- scale_x_discrete(drop = !keep_empty, expand = c(0, 0))
            scaley <- scale_y_continuous(
              labels = if (stat_type == "count") scales::number else scales::percent,
              expand = c(0, 0)
            )
          } else {
            position_use <- position_dodge2(width = 0.9, preserve = "single")
            scalex <- scale_x_discrete(drop = !keep_empty)
            scaley <- scale_y_continuous(
              limits = c(0, max(dat[["value"]], na.rm = TRUE) * 1.1),
              labels = if (stat_type == "count") scales::number else scales::percent,
              expand = c(0, 0)
            )
          }
          if (position == "stack") {
            bg_list <- NULL
          } else {
            bg_df <- unique(dat[, g, drop = FALSE])
            bg_df[, "bg_color"] <- rep_len(c("transparent", "grey80"), nrow(bg_df))
            bg_list <- list()
            for (i in seq_len(nrow(bg_df))) {
              x <- as.numeric(bg_df[i, g])
              bg_list[[i]] <- annotate(
                geom = "rect",
                xmin = ifelse(x == 1, -Inf, x - 0.5),
                xmax = ifelse(x == nlevels(bg_df[, g]), Inf, x + 0.5),
                ymin = -Inf,
                ymax = Inf,
                fill = as.character(bg_df[i, "bg_color"]), alpha = 0.2
              )
            }
          }

          if (plot_type == "bar") {
            p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
              bg_list +
              geom_col(aes(fill = .data[[stat.by]]),
                width = width,
                color = "black",
                alpha = alpha,
                position = position_use
              ) +
              scalex +
              scaley
          }
          if (plot_type == "rose") {
            p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
              bg_list +
              geom_col(aes(fill = .data[[stat.by]]),
                width = width,
                color = "black",
                alpha = alpha,
                position = position_use
              ) +
              scalex +
              scaley +
              coord_polar(theta = "x")
          }
          if (plot_type == "ring" || plot_type == "pie") {
            p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
              bg_list +
              geom_col(aes(fill = .data[[stat.by]]),
                width = width,
                color = "black",
                alpha = alpha,
                position = position_use
              ) +
              scalex +
              scaley +
              coord_polar(theta = "y")
          }
          if (plot_type == "area") {
            p <- ggplot(dat, aes(x = .data[[g]], y = value, group = .data[[stat.by]])) +
              bg_list +
              geom_area(aes(fill = .data[[stat.by]]),
                color = "black",
                alpha = alpha,
                position = position_use
              ) +
              scalex +
              scaley
          }
          if (isTRUE(label)) {
            p <- p + geom_text_repel(aes(
              label = if (stat_type == "count") value else paste0(round(value * 100, 1), "%"),
              y = value # if (position == "stack") value else value / 2
            ),
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100, min.segment.length = 0,
            position = position_use
            )
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
            scale_fill_manual(name = paste0(stat.by, ":"), values = colors, drop = FALSE) +
            do.call(theme_use, list(
              aspect.ratio = aspect.ratio,
              axis.text.x = axis.text.x,
              legend.position = legend.position,
              legend.direction = legend.direction,
              panel.grid.major = element_line(colour = "grey80", linetype = 2)
            )) + guides(fill = guide_legend(
              title.hjust = 0,
              keywidth = 0.05,
              keyheight = 0.05,
              default.unit = "inch",
              order = 1,
              override.aes = list(size = 4.5, color = "black", alpha = 1)
            ))
          plist[[paste0(s, ":", paste0(single_group, collapse = ","))]] <- p
        }
      }
    }
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
    for (s in levels(dat_all[[split.by]])) {
      dat_use <- dat_all[dat_all[[split.by]] == s, ]
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
            colour = label.fg, size = label.size + 1,
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
          theme_void()
        p <- p + labs(title = title, subtitle = subtitle)
      }
      if (plot_type == "upset") {
        check_R("ggupset")
        for (n in seq_len(nrow(dat_use))) {
          dat_use[["intersection"]][n] <- list(stat.by[unlist(dat_use[n, stat.by])])
        }
        dat_use <- dat_use[sapply(dat_use[["intersection"]], length) > 0, ]
        p <- ggplot(dat_use, aes(x = intersection)) +
          geom_bar(aes(fill = ..count..), color = "black", width = 0.5, show.legend = FALSE) +
          geom_text_repel(aes(label = ..count..),
            stat = "count",
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100,
            min.segment.length = 0, segment.colour = "black"
          ) +
          labs(title = title, subtitle = subtitle, x = "", y = "Intersection size") +
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
        check_R("davidsjoberg/ggsankey")
        require("dplyr")
        colors <- palette_scp(unique(unlist(lapply(dat_all[, stat.by, drop = FALSE], levels))), palette = palette, palcolor = palcolor)
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
        }
        if (legend.direction == "vertical") {
          legend <- do.call(cbind, legend_list)
        } else {
          legend <- do.call(rbind, legend_list)
        }
        dat <- suppressWarnings(ggsankey::make_long(dat_use, all_of(stat.by)))
        dat$node <- factor(dat$node, levels = rev(names(colors)))
        p0 <- ggplot(dat, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node)) +
          ggsankey::geom_sankey(color = "black", flow.alpha = alpha, show.legend = FALSE) +
          scale_fill_manual(values = colors, drop = FALSE) +
          scale_x_discrete(expand = c(0, 0.2)) +
          theme_void() +
          theme(axis.text.x = element_text())
        grob <- ggplotGrob(p0)
        if (legend.position == "bottom") {
          grob <- gtable_add_rows(grob, sum(legend$heights), -1)
          grob <- gtable_add_grob(grob, legend, t = dim(grob)[1], l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
        }
        if (legend.position == "top") {
          grob <- gtable_add_rows(grob, sum(legend$heights), 0)
          grob <- gtable_add_grob(grob, legend, t = 1, l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
        }
        if (legend.position == "right") {
          grob <- gtable_add_cols(grob, sum(legend$widths), -1)
          grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])
        }
        if (legend.position == "left") {
          grob <- gtable_add_cols(grob, sum(legend$widths), 0)
          grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1)
        }
        p <- plot_grid(grob)
      }
      if (plot_type == "chord") {
        check_R("circlize")
        colors <- palette_scp(unique(unlist(lapply(dat_all[, stat.by, drop = FALSE], levels))), palette = palette, palcolor = palcolor)
        M <- table(dat_use[[stat.by[1]]], dat_use[[stat.by[2]]])
        m <- matrix(M, ncol = ncol(M), dimnames = dimnames(M))
        circlize::chordDiagram(m,
          grid.col = colors,
          transparency = 0.2,
          link.lwd = 1,
          link.lty = 1,
          link.border = 1
        )
        circlize::circos.clear()
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

      plist[[s]] <- p
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
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
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
#' @param line_size
#' @param line_bg
#' @param line_bg_r
#' @param whiskers
#' @param whiskers_size
#' @param whiskers_alpha
#' @param theme_use
#' @param aspect.ratio
#' @param title
#' @param subtitle
#' @param xlab
#' @param ylab
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param legend.position
#' @param legend.direction
#' @param return_layer
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
LineagePlot <- function(srt, lineages, reduction = NULL, dims = c(1, 2),
                        trim = c(0.01, 0.99), span = 0.75,
                        palette = "Dark2", palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                        line_size = 1, line_bg = "white", line_bg_r = 0.5,
                        whiskers = FALSE, whiskers_size = 0.5, whiskers_alpha = 0.5,
                        theme_use = "theme_scp", aspect.ratio = 1, title = NULL, subtitle = NULL,
                        xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                        legend.position = "right", legend.direction = "vertical",
                        return_layer = FALSE) {
  if (missing(lineages)) {
    stop("'lineages' must be provided")
  }

  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  reduction_key <- Key(srt[[reduction]])
  dat_dim <- Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
  dat_lineages <- srt@meta.data[, unique(lineages), drop = FALSE]
  dat <- cbind(dat_dim, dat_lineages[row.names(dat_dim), , drop = FALSE])
  dat[, "cell"] <- rownames(dat)

  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  colors <- palette_scp(lineages, palette = palette, palcolor = palcolor)
  axes <- paste0(reduction_key, dims)
  fitted_list <- lapply(lineages, function(l) {
    trim_pass <- dat[[l]] > quantile(dat[[l]], trim[1], na.rm = TRUE) & dat[[l]] < quantile(dat[[l]], trim[2], na.rm = TRUE)
    na_pass <- !is.na(dat[[l]])
    index <- which(trim_pass & na_pass)
    index <- index[order(dat[index, l])]
    dat_sub <- dat[index, ]
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
        size = whiskers_size, alpha = whiskers_alpha,
        show.legend = TRUE, inherit.aes = FALSE
      ))
    }
    curve <- c(
      curve,
      geom_path(
        data = dat_smooth, mapping = aes(x = Axis_1, y = Axis_2), color = line_bg,
        size = line_size + line_bg_r, arrow = lineages_arrow,
        show.legend = TRUE, inherit.aes = FALSE
      ),
      geom_path(
        data = dat_smooth, mapping = aes(x = Axis_1, y = Axis_2, color = Lineages),
        size = line_size, arrow = lineages_arrow,
        show.legend = TRUE, inherit.aes = FALSE
      )
    )
    return(curve)
  })
  curve_layer <- c(unlist(curve_layer), list(scale_color_manual(values = colors)))

  lab_layer <- list(labs(title = title, subtitle = subtitle, x = xlab, y = ylab))
  theme_layer <- list(do.call(theme_use, list(
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction,
    xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
  )))

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
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#' @param type
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
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
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
#' PAGAPlot(pancreas_sub, show_transition = TRUE)
#'
#' @importFrom Seurat Reductions Key Embeddings
#' @export
PAGAPlot <- function(srt, paga = srt@misc$paga, type = "connectivities",
                     reduction = NULL, dims = c(1, 2), show_transition = FALSE,
                     node_palette = "Paired", node_size = 4, node_alpha = 1,
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
                     theme_use = "theme_scp", aspect.ratio = 1, title = "PAGA", subtitle = NULL,
                     xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                     legend.position = "right", legend.direction = "vertical",
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
  if (!is.factor(srt[[groups, drop = TRUE]])) {
    srt[[groups, drop = TRUE]] <- factor(srt[[groups, drop = TRUE]])
  }
  if (nlevels(srt[[groups, drop = TRUE]]) != nrow(connectivities)) {
    stop("nlevels in ", groups, " is not identical with the group in paga")
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  reduction_key <- Key(srt[[reduction]])
  dat_dim <- as.data.frame(Embeddings(srt, reduction = reduction))
  colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
  rownames(dat_dim) <- rownames(dat_dim) %||% colnames(srt)
  dat_dim <- dat_dim[, paste0(reduction_key, dims)]
  dat_dim[[groups]] <- srt[[groups, drop = TRUE]][rownames(dat_dim)]
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

  out <- GraphPlot(
    node = dat, edge = as.matrix(connectivities), node_coord = paste0(reduction_key, dims),
    node_group = groups, node_palette = node_palette, node_size = node_size, node_alpha = node_alpha,
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
    theme_use = theme_use, aspect.ratio = aspect.ratio, title = title, subtitle = subtitle,
    xlab = xlab, ylab = ylab, lab_cex = lab_cex, xlen_npc = xlen_npc, ylen_npc = ylen_npc,
    legend.position = legend.position, legend.direction = legend.direction,
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
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#'
#' @importFrom ggplot2 scale_size_identity scale_size_continuous scale_size_discrete scale_alpha_identity scale_alpha_continuous scale_alpha_discrete geom_curve geom_segment geom_point scale_color_manual guide_legend guides labs aes
#' @importFrom ggnewscale new_scale
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow
#' @export
GraphPlot <- function(node, edge, transition = NULL,
                      node_coord = c("x", "y"), node_group = NULL, node_palette = "Paired", node_size = 4, node_alpha = 1,
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
                      theme_use = "theme_scp", aspect.ratio = 1, title = NULL, subtitle = NULL,
                      xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                      legend.position = "right", legend.direction = "vertical",
                      return_layer = FALSE) {
  use_triangular <- match.arg(use_triangular)
  edge_line <- match.arg(edge_line)
  transition_line <- match.arg(transition_line)
  if (missing(node) || missing(edge)) {
    stop("Both nodes and edges must be provided.")
  }
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
    warning("rownames of node is not identical with edge matrix. They will correspond according to the order.")
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
      warning("rownames of node is not identical with transition matrix. They will correspond according to the order.")
      colnames(transition) <- rownames(transition) <- rownames(node) <- rownames(node) %||% colnames(transition) %||% rownames(transition)
    }
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
          data = edge_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
          lineend = "round", linejoin = "mitre", linetype = linetype, color = edge_color, alpha = edge_alpha,
          inherit.aes = FALSE, show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[edge_df[["edge_name"]] %in% edge_highlight, ]
        edge_layer <- c(
          edge_layer,
          list(
            geom_segment(
              data = edge_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
              lineend = "round", linejoin = "mitre", linetype = linetype, color = edge_highlight_color, alpha = 1,
              inherit.aes = FALSE, show.legend = FALSE
            )
          )
        )
      }
    } else {
      edge_layer <- list(
        geom_curve(
          data = edge_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
          curvature = edge_line_curvature, angle = edge_line_angle,
          lineend = "round", linetype = linetype, color = edge_color, alpha = edge_alpha,
          inherit.aes = FALSE, show.legend = FALSE
        )
      )
      if (!is.null(edge_highlight)) {
        edge_df_highlight <- edge_df[edge_df[["edge_name"]] %in% edge_highlight, ]
        edge_layer <- c(
          edge_layer,
          list(
            geom_curve(
              data = edge_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
              curvature = edge_line_curvature, angle = edge_line_angle,
              lineend = "round", linetype = linetype, color = edge_highlight_color, alpha = 1,
              inherit.aes = FALSE, show.legend = FALSE
            )
          )
        )
      }
    }
    edge_layer <- c(edge_layer, list(
      scale_size_continuous(range = range(edge_size), guide = "none"),
      new_scale("size")
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
            data = trans_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
            arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
            lineend = "round", linejoin = "mitre", color = transition_color, alpha = transition_alpha,
            inherit.aes = FALSE, show.legend = FALSE
          )
        )
        if (!is.null(transition_highlight)) {
          trans_df_highlight <- trans_df[trans_df[["trans_name"]] %in% transition_highlight, ]
          trans_layer <- c(
            trans_layer,
            list(
              geom_segment(
                data = trans_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
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
            data = trans_df, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
            arrow = arrow(angle = transition_arrow_angle, type = transition_arrow_type, length = transition_arrow_length),
            curvature = transition_line_curvature, angle = transition_line_angle,
            lineend = "round", color = transition_color, alpha = transition_alpha,
            inherit.aes = FALSE, show.legend = FALSE
          )
        )
        if (!is.null(edge_highlight)) {
          trans_df_highlight <- trans_df[trans_df[["trans_name"]] %in% transition_highlight, ]
          trans_layer <- c(
            trans_layer,
            list(
              geom_curve(
                data = trans_df_highlight, mapping = aes(x = x, y = y, xend = xend, yend = yend, size = size),
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
        scale_size_continuous(range = range(transition_size), guide = "none"),
        new_scale("size")
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
    node_highlight <- node[node[["node_name"]] %in% node_highlight, ]
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
      name = node_group, values = palette_scp(node[["node_group"]], palette = node_palette), labels = label_use,
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
  theme_layer <- list(do.call(theme_use, list(
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction,
    xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
  )))

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
#' @param streamline_jitter
#' @param streamline_size
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
#' @param lab_cex
#' @param xlen_npc
#' @param ylen_npc
#' @param legend.position
#' @param legend.direction
#' @param return_layer
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
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
VelocityPlot <- function(srt, reduction, dims = c(1, 2), velocity = "stochastic", plot_type = c("raw", "grid", "stream"), group_by = NULL, group_palette = "Paired",
                         n_neighbors = ceiling(ncol(srt) / 50), density = 1, smooth = 0.5, scale = 1, min_mass = 1, cutoff_perc = 5,
                         arrow_angle = 20, arrow_flank = 0.8, arrow_color = "black",
                         streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15, streamline_jitter = 1,
                         streamline_size = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
                         streamline_bg_color = "white", streamline_bg_stroke = 0.5,
                         theme_use = "theme_scp", aspect.ratio = 1, title = "Cell velocity", subtitle = NULL,
                         xlab = NULL, ylab = NULL, lab_cex = 1, xlen_npc = 0.15, ylen_npc = 0.15,
                         legend.position = "right", legend.direction = "vertical",
                         return_layer = FALSE, seed = 11) {
  plot_type <- match.arg(plot_type)
  set.seed(seed)
  check_R("metR")

  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  V_reduction <- paste0(velocity, "_", reduction)
  if (!V_reduction %in% Reductions(srt)) {
    stop("Cannot find the velocity embedding ", V_reduction, ".")
  }
  X_emb <- Embeddings(srt, reduction)[, dims]
  V_emb <- Embeddings(srt, V_reduction)[, dims]

  reduction_key <- Key(srt[[reduction]])
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])

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
          name = group_by, values = palette_scp(df_field[["group_by"]], palette = group_palette),
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
          jitter = streamline_jitter, n = streamline_n, size = max(streamline_size) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          jitter = streamline_jitter, n = streamline_n, size = max(streamline_size), color = streamline_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          jitter = streamline_jitter, n = streamline_n, linetype = 0, color = arrow_color,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        )
      )
    } else {
      velocity_layer <- list(
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          jitter = streamline_jitter, n = streamline_n, size = max(streamline_size) + streamline_bg_stroke, color = streamline_bg_color, alpha = streamline_alpha,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v, size = ..step.., color = sqrt(..dx..^2 + ..dy..^2)),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          jitter = streamline_jitter, n = streamline_n, alpha = streamline_alpha,
          arrow = NULL, lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        metR::geom_streamline(
          data = df_field, aes(x = x, y = y, dx = u, dy = v),
          L = streamline_L, min.L = streamline_minL, res = streamline_res,
          jitter = streamline_jitter, n = streamline_n, linetype = 0, color = arrow_color,
          arrow.type = "closed", arrow.angle = arrow_angle,
          lineend = "round", linejoin = "mitre", inherit.aes = FALSE
        ),
        scale_color_gradientn(
          name = "Velocity", colors = palette_scp(palette = streamline_palette, palcolor = streamline_palcolor),
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
        ),
        scale_size(range = range(streamline_size), guide = "none")
      )
    }
  }

  lab_layer <- list(labs(title = title, subtitle = subtitle, x = xlab, y = ylab))
  theme_layer <- list(do.call(theme_use, list(
    aspect.ratio = aspect.ratio,
    legend.position = legend.position,
    legend.direction = legend.direction,
    xlab = xlab, ylab = ylab, xlen_npc = xlen_npc, ylen_npc = ylen_npc, lab_cex = lab_cex
  )))


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
#' @importFrom cowplot plot_grid
#' @export
VolcanoPlot <- function(srt, group_by = NULL, test.use = "wilcox", DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                        x_metric = "diff_pct", palette = "RdBu", palcolor = NULL, pt.size = 1, pt.alpha = 1,
                        cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                        nlabel = 5, features_label = NULL, label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                        aspect.ratio = NULL, xlab = x_metric, ylab = "-log10(p-adjust)",
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr") {
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
    value_range <- min(abs(c(x_upper, x_lower)))
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
    de_df <- de_df[order(abs(de_df[, "avg_log2FC"]), decreasing = FALSE, na.last = FALSE), ]
  } else if (x_metric == "avg_log2FC") {
    de_df[, "x"] <- de_df[, "avg_log2FC"]
    de_df[de_df[, "diff_pct"] < 0, "y"] <- -de_df[de_df[, "diff_pct"] < 0, "y"]
    de_df <- de_df[order(abs(de_df[, "diff_pct"]), decreasing = FALSE, na.last = FALSE), ]
  }
  de_df[, "distance"] <- de_df[, "x"]^2 + de_df[, "y"]^2

  plist <- list()
  for (group in levels(de_df[["group1"]])) {
    df <- de_df[de_df[["group1"]] == group, ]
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
      geom_point(data = df[!df[["DE"]] & !df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[!df[["DE"]] & df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], ], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight) +
      geom_point(data = df[df[["DE"]] & df[["border"]], ], aes(x = x, y = y), color = cols.highlight, size = sizes.highlight + stroke.highlight, alpha = alpha.highlight, position = jitter) +
      geom_point(data = df[df[["DE"]] & !df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha) +
      geom_point(data = df[df[["DE"]] & df[["border"]], ], aes(x = x, y = y, color = .data[[color_by]]), size = pt.size, alpha = pt.alpha, position = jitter) +
      geom_hline(yintercept = 0, color = "black", linetype = 1) +
      geom_vline(xintercept = 0, color = "grey", linetype = 2) +
      geom_text_repel(
        data = df[df[["label"]], ], aes(x = x, y = y, label = gene),
        min.segment.length = 0, max.overlaps = 100, segment.colour = "grey40",
        color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, force = 20,
        nudge_x = ifelse(df[df[["label"]], "y"] >= 0, -x_nudge, x_nudge)
      ) +
      labs(x = xlab, y = ylab) +
      scale_color_gradientn(
        name = ifelse(x_metric == "diff_pct", "log2FC", "diff_pct"), colors = palette_scp(palette = palette, palcolor = palcolor),
        values = rescale(unique(c(min(df[, color_by], 0), 0, max(df[, color_by])))),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
      ) +
      scale_y_continuous(labels = abs) +
      facet_wrap(~group1) +
      theme_scp(aspect.ratio = aspect.ratio)
    plist[[group]] <- p
  }
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
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
      warning("rownames of node is not identical with edge matrix. They will correspond according to the order.")
      colnames(edge) <- rownames(edge) <- rownames(node) <- rownames(node) %||% colnames(edge) %||% rownames(edge)
    }
    node[["name"]] <- node[[node_group]]
    edge_df <- reshape2::melt(edge, na.rm = TRUE, stringsAsFactors = FALSE)
    colnames(edge_df) <- c("from", "to", "size")
    edge_df[["from"]] <- match(edge_df[["from"]], node[["name"]]) - 1
    edge_df[["to"]] <- match(edge_df[["to"]], node[["name"]]) - 1
    edge_df <- edge_df[edge_df[["size"]] > 0, ]
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

#' SummaryPlot
#' @inheritParams ClassDimPlot
#' @inheritParams ExpDimPlot
#' @inheritParams panel_fix
#' @param Class_palette Name of palette to use for ClassDimPlot.
#' @param Exp_palette Name of palette to use for ExpDimPlot.
#' @param size Size of each plot.
#'
#' @return A grob obejct with size attributes.
#' @examples
#' data("pancreas_sub")
#' p <- SummaryPlot(pancreas_sub,
#'   group.by = c("CellType", "SubCellType"), group.split.by = "SubCellType",
#'   features = c("S_score", "G2M_score", "nFeature_RNA", "nCount_RNA")
#' )
#' p
#'
#' ## Save the plot with appropriate size
#' # ggplot2::ggsave(
#' #   filename = "summaryplot.png", plot = p,
#' #   units = attr(p, "size")$units, width = attr(p, "size")$width, height = attr(p, "size")$height
#' # )
#'
#' @export
SummaryPlot <- function(srt,
                        group.by = NULL,
                        group.split.by = NULL,
                        features = NULL,
                        reduction = NULL,
                        Class_palette = "Paired",
                        Exp_palette = "Spectral",
                        size = 2, margin = 0.1, units = "in") {
  if (is.null(group.by) && is.null(features)) {
    stop("features or 'group.by' must be provided.")
  }
  if (is.null(group.by) || is.null(group.split.by)) {
    for (i in c(group.by, group.split.by)) {
      if (!i %in% colnames(srt@meta.data)) {
        stop(i, " is not in the meta.data of srt object.")
      }
      if (!is.factor(srt[[i, drop = TRUE]])) {
        srt[[i, drop = TRUE]] <- factor(srt[[i, drop = TRUE]], levels = unique(srt[[i, drop = TRUE]]))
      }
    }
  }
  features_drop <- features[!features %in% c(rownames(srt), colnames(srt@meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  class_features <- c()
  exp_features <- c()
  for (i in features) {
    if (i %in% colnames(srt@meta.data)) {
      if (!is.numeric(srt@meta.data[[i]])) {
        class_features <- c(class_features, i)
      } else {
        exp_features <- c(exp_features, i)
      }
    } else {
      exp_features <- c(exp_features, i)
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% Reductions(srt)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }

  p1_list <- p2_list <- p3_list <- list()
  if (!is.null(group.by)) {
    p1_list <- ClassDimPlot(srt,
      reduction = reduction,
      group.by = group.by,
      palette = Class_palette,
      force = TRUE,
      combine = FALSE
    )
  }
  if (!is.null(group.split.by)) {
    for (sp in group.split.by) {
      p2 <- ClassDimPlot(srt,
        reduction = reduction,
        group.by = sp,
        split.by = sp,
        palette = Class_palette,
        legend.position = "none",
        label = FALSE,
        force = TRUE,
        combine = FALSE
      )
      p2_list <- c(p2_list, p2)
    }
  }
  if (length(class_features) > 0) {
    p3 <- ClassDimPlot(
      srt = srt,
      reduction = reduction,
      group.by = class_features,
      palette = Class_palette,
      force = TRUE,
      combine = FALSE
    )
    p3_list <- c(p3_list, p3)
  }
  if (length(exp_features) > 0) {
    p3 <- ExpDimPlot(
      srt = srt,
      reduction = reduction,
      features = exp_features,
      palette = Exp_palette,
      force = TRUE,
      combine = FALSE,
    )
    p3_list <- c(p3_list, p3)
  }

  p1_ncol <- ceiling(sqrt(length(p1_list)))
  p2_ncol <- ceiling(sqrt(length(p2_list)))
  p3_ncol <- ceiling(sqrt(length(p3_list)))
  ncol <- max(p1_ncol, p2_ncol, p3_ncol)
  p1_ncol <- min(p1_ncol, ncol)
  p2_ncol <- min(p2_ncol, ncol)
  p3_ncol <- min(p3_ncol, ncol)
  p1_nrow <- ceiling(length(p1_list) / p1_ncol)
  p2_nrow <- ceiling(length(p2_list) / p2_ncol)
  p3_nrow <- ceiling(length(p3_list) / p3_ncol)

  p1_height <- p1_nrow
  p2_height <- p2_nrow
  p3_height <- p3_nrow

  if (!is.null(size)) {
    p1_width <- p1_height <- p2_width <- p2_height <- p3_width <- p3_height <- 0
    if (length(p1_list) > 0) {
      p1_list <- lapply(p1_list, function(p) {
        panel_fix(p, width = size, height = size, units = units, margin = margin)
      })
      p1_size <- do.call(rbind.data.frame, lapply(p1_list, function(p) unlist(attr(p, "size"))))
      colnames(p1_size) <- c("width", "height", "units")
      p1_size[, "nrow"] <- ceiling(seq_len(nrow(p1_size)) / p1_ncol)
      p1_size[, "ncol"] <- ceiling(seq_len(nrow(p1_size)) %% p1_ncol)
      p1_size[, "ncol"][p1_size[, "ncol"] == 0] <- p1_ncol
      p1_width <- max(aggregate(p1_size[, "width", drop = FALSE], by = list(p1_size$nrow), FUN = function(x) {
        sum(as.numeric(x))
      })[["width"]])
      p1_height <- max(aggregate(p1_size[, "height", drop = FALSE], by = list(p1_size$ncol), FUN = function(x) {
        sum(as.numeric(x))
      })[["height"]])
    }

    if (length(p2_list) > 0) {
      p2_list <- lapply(p2_list, function(p) {
        panel_fix(p, width = size, height = size, units = units, margin = margin)
      })
      p2_size <- do.call(rbind.data.frame, lapply(p2_list, function(p) unlist(attr(p, "size"))))
      colnames(p2_size) <- c("width", "height", "units")
      p2_size[, "nrow"] <- ceiling(seq_len(nrow(p2_size)) / p2_ncol)
      p2_size[, "ncol"] <- ceiling(seq_len(nrow(p2_size)) %% p2_ncol)
      p2_size[, "ncol"][p2_size[, "ncol"] == 0] <- p2_ncol
      p2_width <- max(aggregate(p2_size[, "width", drop = FALSE], by = list(p2_size$nrow), FUN = function(x) {
        sum(as.numeric(x))
      })[["width"]])
      p2_height <- max(aggregate(p2_size[, "height", drop = FALSE], by = list(p2_size$ncol), FUN = function(x) {
        sum(as.numeric(x))
      })[["height"]])
    }

    if (length(p3_list) > 0) {
      p3_list <- lapply(p3_list, function(p) {
        panel_fix(p, width = size, height = size, units = units, margin = margin)
      })
      p3_size <- do.call(rbind.data.frame, lapply(p3_list, function(p) unlist(attr(p, "size"))))
      colnames(p3_size) <- c("width", "height", "units")
      p3_size[, "nrow"] <- ceiling(seq_len(nrow(p3_size)) / p3_ncol)
      p3_size[, "ncol"] <- ceiling(seq_len(nrow(p3_size)) %% p3_ncol)
      p3_size[, "ncol"][p3_size[, "ncol"] == 0] <- p3_ncol
      p3_width <- max(aggregate(p3_size[, "width", drop = FALSE], by = list(p3_size$nrow), FUN = function(x) {
        sum(as.numeric(x))
      })[["width"]])
      p3_height <- max(aggregate(p3_size[, "height", drop = FALSE], by = list(p3_size$ncol), FUN = function(x) {
        sum(as.numeric(x))
      })[["height"]])
    }

    height <- sum(p1_height, p2_height, p3_height)
    width <- max(p1_width, p2_width, p3_width)
  }

  p1 <- p2 <- p3 <- NULL
  if (length(p1_list) > 0) {
    p1 <- plot_grid(plotlist = p1_list, align = "hv", axis = "tblr", nrow = p1_nrow)
  }
  if (length(p2_list) > 0) {
    p2 <- plot_grid(plotlist = p2_list, align = "hv", axis = "tblr", nrow = p2_nrow)
  }
  if (length(p3_list) > 0) {
    p3 <- plot_grid(plotlist = p3_list, align = "hv", axis = "tblr", nrow = p3_nrow)
  }
  p <- plot_grid(p1, p2, p3, align = "hv", axis = "tblr", ncol = 1, rel_heights = c(p1_height, p2_height, p3_height)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

  if (!is.null(size)) {
    p <- plot_grid(panel_fix_single(p, width = width, height = height, units = units))
    attr(p, "size") <- list(width = width, height = height, units = units)
  }
  return(p)
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
#' @param libsize
#' @param order.by
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
#' @param add_density
#' @param aspect.ratio
#' @param legend.position
#' @param legend.direction
#' @param combine
#' @param nrow
#' @param ncol
#' @param byrow
#' @param align
#' @param axis
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#' @importFrom ggplot2 geom_line geom_point geom_ribbon geom_rug stat_density2d facet_grid expansion
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid get_legend
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom grDevices colorRampPalette
#' @importFrom stats runif
#' @export
DynamicPlot <- function(srt, features, lineages, slot = "counts", assay = "RNA", family = NULL,
                        exp_method = c("log1p", "zscore", "log2fc", "raw"), lib_normalize = TRUE, libsize = NULL,
                        order.by = "pseudotime", group.by = NULL, compare_lineages = TRUE, compare_features = FALSE,
                        add_line = TRUE, add_interval = TRUE, line.size = 1, line_palette = "Dark2", line_palcolor = NULL,
                        add_point = TRUE, pt.size = 1, point_palette = "Paired", point_palcolor = NULL,
                        add_rug = TRUE, add_density = TRUE,
                        aspect.ratio = 0.5, legend.position = "right", legend.direction = "vertical",
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr") {
  check_R("MatrixGenerics")
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    stop(group.by, " is not in the meta.data of srt object.")
  }

  gene <- features[features %in% rownames(srt[[assay]]@counts)]
  meta <- features[features %in% colnames(srt@meta.data)]
  features <- c(gene, meta)
  if (length(features) == 0) {
    stop("No feature found in the srt object.")
  }
  exp_method <- match.arg(exp_method)

  # if (isTRUE(compare_lineages) && exp_method != "raw") {
  #   stop("exp_method must be 'raw' when compare_lineages is TRUE")
  # }
  # if (length(meta) > 0 && exp_method != "raw") {
  #   stop("exp_method must be 'raw' when features are not all features.")
  # }

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
  if (order.by == "rank") {
    x_assign <- rank(x_assign)
  }
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
  Y_libsize <- colSums(GetAssayData(srt, slot = "counts", assay = assay))
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
          warning("Values in 'counts' slot is non-integer. Set the libsize to 1.", immediate. = TRUE)
        }
      }
      raw_matrix[, gene] <- raw_matrix[, gene] / libsize_use * median(Y_libsize)
    }
    if (exp_method == "raw") {
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
    raw_matrix[is.infinite(raw_matrix)] <- max(abs(raw_matrix[!is.infinite(raw_matrix)])) * ifelse(raw_matrix[is.infinite(raw_matrix)] > 0, 1, -1)
    fitted_matrix[is.infinite(fitted_matrix)] <- max(abs(fitted_matrix[!is.infinite(fitted_matrix)])) * ifelse(fitted_matrix[is.infinite(fitted_matrix)] > 0, 1, -1)
    upr_matrix[is.infinite(upr_matrix)] <- max(abs(upr_matrix[!is.infinite(upr_matrix)])) * ifelse(upr_matrix[is.infinite(upr_matrix)] > 0, 1, -1)
    lwr_matrix[is.infinite(lwr_matrix)] <- max(abs(lwr_matrix[!is.infinite(lwr_matrix)])) * ifelse(lwr_matrix[is.infinite(lwr_matrix)] > 0, 1, -1)

    raw <- as.data.frame(cbind(cell_metadata[rownames(raw_matrix), c(l, "x_assign")], raw_matrix))
    colnames(raw)[1] <- "Pseudotime"
    raw[["Cell"]] <- rownames(raw)
    raw[["Value"]] <- "raw"
    raw <- reshape2::melt(raw, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = exp_method, variable.name = "Features")

    fitted <- as.data.frame(cbind(cell_metadata[rownames(fitted_matrix), c(l, "x_assign")], fitted_matrix))
    colnames(fitted)[1] <- "Pseudotime"
    fitted[["Cell"]] <- rownames(fitted)
    fitted[["Value"]] <- "fitted"
    fitted <- reshape2::melt(fitted, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = exp_method, variable.name = "Features")

    upr <- as.data.frame(cbind(cell_metadata[rownames(upr_matrix), c(l, "x_assign")], upr_matrix))
    colnames(upr)[1] <- "Pseudotime"
    upr[["Cell"]] <- rownames(upr)
    upr[["Value"]] <- "upr"
    upr <- reshape2::melt(upr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = exp_method, variable.name = "Features")

    lwr <- as.data.frame(cbind(cell_metadata[rownames(lwr_matrix), c(l, "x_assign")], lwr_matrix))
    colnames(lwr)[1] <- "Pseudotime"
    lwr[["Cell"]] <- rownames(lwr)
    lwr[["Value"]] <- "lwr"
    lwr <- reshape2::melt(lwr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = exp_method, variable.name = "Features")

    raw[["upr"]] <- NA
    raw[["lwr"]] <- NA
    fitted[["upr"]] <- upr[[exp_method]]
    fitted[["lwr"]] <- lwr[[exp_method]]

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

  df_all <- df_all[sample(seq_len(nrow(df_all))), ]

  plist <- list()
  legend <- NULL
  if (isTRUE(compare_lineages)) {
    lineages_use <- list(lineages)
    formula <- ". ~ Features"
  } else {
    lineages_use <- lineages
    formula <- "Lineages ~ Features"
  }
  if (isTRUE(compare_features)) {
    features_use <- list(features)
    formula <- ".~."
  } else {
    features_use <- features
  }
  for (l in lineages_use) {
    for (f in features_use) {
      if (all(f %in% meta)) {
        data_nm <- f
        exp_name <- NULL
      } else if (all(!f %in% meta)) {
        data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
        data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
        exp_name <- NULL
      } else {
        exp_name <- exp_method
      }
      if (is.null(exp_name)) {
        exp_name <- switch(exp_method,
          "raw" = data_nm,
          "zscore" = paste0("Z-score(", data_nm, ")"),
          "log2fc" = paste0("Log2(", data_nm, "FC)"),
          "log1p" = paste0("Log(", data_nm, "+1)")
        )
      }

      df <- subset(df_all, df_all[["Lineages"]] %in% l & df_all[["Features"]] %in% f)
      random_noise <- runif(nrow(df), -0.01 * diff(range(df[, exp_method], na.rm = TRUE)), 0.01 * diff(range(df[, exp_method], na.rm = TRUE)))
      df[, "random_noise"] <- random_noise
      df_point <- unique(df[df[["Value"]] == "raw", c("Cell", "x_assign", exp_method, group.by)])
      if (isTRUE(compare_features)) {
        raw_point <- NULL
      } else {
        if (isTRUE(add_point)) {
          if (is.null(group.by)) {
            raw_point <- geom_point(data = df_point, mapping = aes(x = .data[["x_assign"]], y = .data[[exp_method]]), size = pt.size, alpha = 0.8)
          } else {
            raw_point <- list(
              geom_point(data = df_point, mapping = aes(x = .data[["x_assign"]], y = .data[[exp_method]], color = .data[[group.by]]), size = pt.size, alpha = 0.8),
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
          rug <- list(geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]]), alpha = 0.1, length = unit(0.05, "npc"), show.legend = FALSE))
        } else {
          rug <- list(
            geom_rug(data = df_point, mapping = aes(x = .data[["x_assign"]], color = .data[[group.by]]), alpha = 0.5, length = unit(0.05, "npc"), show.legend = isTRUE(compare_features)),
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
            mapping = aes(x = .data[["Pseudotime"]], y = .data[[exp_method]], ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["Lineages"]], group = .data[["LineagesFeatures"]]),
            alpha = 0.4, color = "grey90", show.legend = isTRUE(compare_features)
          ),
          scale_fill_manual(values = palette_scp(df[["Lineages"]], palette = line_palette, palcolor = line_palcolor)),
          new_scale_fill()
        )
      } else {
        interval <- NULL
      }
      if (isTRUE(compare_features)) {
        density <- NULL
      } else {
        if (isTRUE(add_density)) {
          density <- list(
            stat_density2d(data = subset(df, df[["Value"]] == "raw"), geom = "raster", aes(x = .data[["x_assign"]], y = .data[[exp_method]] + .data[["random_noise"]], fill = ..density..^0.25, alpha = 1), contour = FALSE, show.legend = FALSE),
            stat_density2d(data = subset(df, df[["Value"]] == "raw"), geom = "raster", aes(x = .data[["x_assign"]], y = .data[[exp_method]] + .data[["random_noise"]], fill = ..density..^0.25, alpha = ifelse(..density..^0.25 < 0.4, 0, 1)), contour = FALSE, show.legend = FALSE),
            scale_fill_gradientn(colours = colorRampPalette(c("white", "grey90"))(256)),
            new_scale_fill()
          )
        } else {
          density <- NULL
        }
      }
      if (isTRUE(compare_features)) {
        line <- list(
          geom_line(
            data = subset(df, df[["Value"]] == "fitted"),
            mapping = aes(x = .data[["Pseudotime"]], y = .data[[exp_method]], color = .data[["Features"]], group = .data[["LineagesFeatures"]]), size = line.size, alpha = 0.8
          ),
          scale_color_manual(
            values = palette_scp(df[["Features"]], palette = line_palette, palcolor = line_palcolor),
            guide = guide_legend(override.aes = list(alpha = 1, size = 2), order = 2)
          ),
          new_scale_color()
        )
      } else {
        if (isTRUE(add_line)) {
          line <- list(
            geom_line(
              data = subset(df, df[["Value"]] == "fitted"),
              mapping = aes(x = .data[["Pseudotime"]], y = .data[[exp_method]], color = .data[["Lineages"]], group = .data[["LineagesFeatures"]]), size = line.size, alpha = 0.8
            ),
            scale_color_manual(
              values = palette_scp(df[["Lineages"]], palette = line_palette, palcolor = line_palcolor),
              guide = guide_legend(override.aes = list(alpha = 1, size = 2), order = 2)
            ),
            new_scale_color()
          )
        } else {
          line <- NULL
        }
      }

      p <- ggplot() +
        scale_y_continuous(expand = expansion(c(0.1, 0.05))) +
        density +
        raw_point +
        rug +
        interval +
        line +
        labs(x = "Pseudotime", y = exp_name) +
        facet_grid(formula(formula), scales = "free") +
        theme_scp(
          aspect.ratio = aspect.ratio,
          legend.position = "bottom",
          legend.direction = legend.direction
        )
      if (is.null(legend)) {
        legend <- get_legend(p)
      }
      plist[[paste(paste0(l, collapse = "_"), paste0(f, collapse = "_"), sep = ".")]] <- p + theme(legend.position = "none")
    }
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
    } else {
      plot <- plist[[1]]
    }
    if (!is.null(legend)) {
      grob <- ggplotGrob(plot)
      if (legend.position == "bottom") {
        grob <- gtable_add_rows(grob, sum(legend$heights), -1)
        grob <- gtable_add_grob(grob, legend, t = dim(grob)[1], l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
      }
      if (legend.position == "top") {
        grob <- gtable_add_rows(grob, sum(legend$heights), 0)
        grob <- gtable_add_grob(grob, legend, t = 1, l = min(grob$layout[grepl(pattern = "panel", x = grob$layout$name), "l"]))
      }
      if (legend.position == "right") {
        grob <- gtable_add_cols(grob, sum(legend$widths), -1)
        grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = dim(grob)[2])
      }
      if (legend.position == "left") {
        grob <- gtable_add_cols(grob, sum(legend$widths), 0)
        grob <- gtable_add_grob(grob, legend, t = mean(grob$layout[grep("panel", grob$layout$name), "t"]), l = 1)
      }
      plot <- plot_grid(grob)
    }
    return(plot)
  } else {
    return(plist)
  }
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
#' @param row_title_size
#' @param n_split
#' @param heatmap_split_by
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
#' @param cell_palette
#' @param cell_palcolor
#' @param feature_annotation
#' @param feature_palette
#' @param feature_palcolor
#' @param reverse_ht
#' @param use_raster
#' @param height
#' @param width
#' @param units
#' @param seed
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#' ht_result1 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1"), order_by = "peaktime",
#'   cell_annotation = "SubCellType",
#'   n_split = 5, split_method = "kmeans-peaktime",
#'   height = 6, width = 7, use_raster = FALSE
#' )
#' ht_result1$plot
#'
#' ht_result2 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 5, reverse_ht = "Lineage1",
#'   height = 6, width = 7, use_raster = FALSE
#' )
#' ht_result2$plot
#'
#' ht_result3 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 5, reverse_ht = "Lineage1",
#'   species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
#'   height = 6, width = 7, use_raster = FALSE
#' )
#' ht_result3$plot
#' @importFrom Seurat GetAssayData NormalizeData DefaultAssay
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap Legend HeatmapAnnotation anno_empty anno_mark anno_simple anno_textbox draw decorate_heatmap_body width.HeatmapAnnotation height.HeatmapAnnotation
#' @importFrom stats kmeans
#' @importFrom cowplot ggdraw draw_grob
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid grid.lines convertUnit
#' @importFrom dplyr group_by filter arrange desc across mutate summarise distinct n .data
#' @export
DynamicHeatmap <- function(srt, lineages, feature_from = lineages,
                           exp_method = c("zscore", "raw", "log2fc", "log1p"),
                           slot = "counts", assay = "RNA",
                           use_fitted = FALSE, lib_normalize = TRUE, libsize = NULL,
                           min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05,
                           cell_density = 1, order_by = c("peaktime", "valleytime"), decreasing = FALSE,
                           feature_split = NULL, row_title_size = 12,
                           n_split = NULL, heatmap_split_by = lineages,
                           split_method = c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime"), fuzzification = NULL, show_fuzzification = FALSE,
                           anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                           IDtype = "symbol", species = "Homo_sapiens", db_IDtype = "symbol", db_update = FALSE, db_version = "latest", Ensembl_version = 103, mirror = NULL,
                           enrichment = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, universe = NULL,
                           GO_simplify = FALSE, GO_simplify_padjustCutoff = 0.2, simplify_method = "Rel", simplify_similarityCutoff = 0.7,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_id = TRUE, topWord = 20, min_word_length = 3, exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "system", "regulation", "positive", "negative", "response", "process"),
                           anno_width = unit(c(4, 2, 2), "in"), anno_size = c(6, 10),
                           nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                           pseudotime_label = NULL, pseudotime_label_color = "black",
                           pseudotime_label_linetype = 2, pseudotime_label_linewidth = 2,
                           heatmap_palette = "RdBu", pseudotime_palette = "cividis", feature_split_palette = "jama",
                           cell_annotation = NULL, cell_palette = "Paired", cell_palcolor = NULL,
                           feature_annotation = NULL, feature_palette = "Dark2", feature_palcolor = NULL,
                           reverse_ht = NULL, use_raster = NULL, height = NULL, width = NULL, units = "inch",
                           seed = 11) {
  set.seed(seed)
  exp_method <- match.arg(exp_method)
  split_method <- match.arg(split_method)
  order_by <- match.arg(order_by)
  data_nm <- c(ifelse(isTRUE(lib_normalize), "normalized", ""), slot)
  data_nm <- paste(data_nm[data_nm != ""], collapse = " ")
  exp_name <- switch(exp_method,
    "raw" = data_nm,
    "zscore" = paste0("Z-score(", data_nm, ")"),
    "log2fc" = paste0("Log2(", data_nm, "FC)"),
    "log1p" = paste0("Log(", data_nm, "+1)")
  )
  if (missing(srt)) {
    stop("srt must be provided.")
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(srt)
  }
  if (missing(lineages)) {
    stop("lineages must be provided.")
  }
  if (any(!feature_from %in% lineages)) {
    stop("feature_from must be a subset of the lineages")
  }
  if (any(!heatmap_split_by %in% lineages)) {
    stop("'heatmap_split_by' must be a subset of lineages.")
  }
  if (!split_method %in% c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime")) {
    stop("'split_method' must be one of 'mfuzz', 'kmeans', 'kmeans-peaktime', 'hclust', 'hclust-peaktime'.")
  }
  if (!is.null(feature_split) && is.null(names(feature_split))) {
    stop("'feature_split' must be named with names.")
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
    if (length(unique(length(pseudotime_label), length(pseudotime_label_color), length(pseudotime_label_linetype), length(pseudotime_label_linewidth))) != 1) {
      stop("Parameters for the pseudotime_label must be the same length!")
    }
  }
  if (!is.null(cell_annotation)) {
    if (length(cell_palette) == 1) {
      cell_palette <- rep(cell_palette, length(cell_annotation))
    }
    if (length(cell_palcolor) == 1) {
      cell_palcolor <- rep(cell_palcolor, length(cell_annotation))
    }
    if (length(unique(length(cell_palette), length(cell_palcolor), length(cell_annotation))) != 1) {
      stop("cell_palette and cell_palcolor must be the same length as cell_annotation")
    }
    if (any(!cell_annotation %in% colnames(srt@meta.data))) {
      stop("cell_annotation: ", cell_annotation[!cell_annotation %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
    }
  }
  if (!is.null(feature_annotation)) {
    if (length(feature_palette) == 1) {
      feature_palette <- rep(feature_palette, length(feature_annotation))
    }
    if (length(feature_palcolor) == 1) {
      feature_palcolor <- rep(feature_palcolor, length(feature_annotation))
    }
    if (length(unique(length(feature_palette), length(feature_palcolor), length(feature_annotation))) != 1) {
      stop("feature_palette and feature_palcolor must be the same length as feature_annotation")
    }
    if (any(!feature_annotation %in% colnames(srt[[assay]]@meta.data))) {
      stop("feature_annotation: ", feature_annotation[!feature_annotation %in% colnames(srt[[assay]]@meta.data)], " is not in the meta data of the ", assay, " assay in the Seurat object.")
    }
  }

  feature_union <- c()
  cell_union <- c()
  dynamic <- list()
  for (l in lineages) {
    if (!paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      stop(l, " info not found in the srt object. Should perform RunDynamicFeatures first!")
    }
    DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
    DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > min_expcells & DynamicFeatures$r.sq > r.sq & DynamicFeatures$dev.expl > dev.expl & DynamicFeatures$padjust < padjust, ]
    dynamic[[l]] <- DynamicFeatures
    if (l %in% feature_from) {
      feature_union <- c(feature_union, DynamicFeatures[, "features"])
    }
    cell_union <- c(cell_union, rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]))
    if ("lineages" %in% names(srt@tools[[paste0("DynamicFeatures_", l)]])) {
      pseudotime <- srt@tools[[paste0("DynamicFeatures_", l)]][["lineages"]]
      if (!l %in% colnames(srt@meta.data)) {
        srt[[l]] <- srt[[pseudotime]]
      }
    }
  }
  feature_union <- unique(feature_union)
  gene <- feature_union[feature_union %in% rownames(srt[[assay]])]
  meta <- feature_union[feature_union %in% colnames(srt@meta.data)]
  cell_union <- unique(cell_union)
  Pseudotime_assign <- rowMeans(srt@meta.data[cell_union, lineages, drop = FALSE], na.rm = TRUE)

  cell_metadata <- cbind.data.frame(data.frame(row.names = cell_union, cells = cell_union),
    Pseudotime_assign = Pseudotime_assign,
    srt@meta.data[cell_union, lineages, drop = FALSE]
  )
  if (cell_density != 1) {
    cell_bin <- cut(Pseudotime_assign, breaks = seq(min(Pseudotime_assign), max(Pseudotime_assign), length.out = 100), include.lowest = TRUE)
    ncell_bin <- ceiling(max(table(cell_bin)) * cell_density)
    message("ncell/bin=", ncell_bin, "(100bins)")
    cell_keep <- unlist(sapply(levels(cell_bin), function(x) {
      cells <- names(Pseudotime_assign)[cell_bin == x]
      out <- sample(cells, size = min(length(cells), ncell_bin))
      return(out)
    }))
    cell_metadata <- cell_metadata[cell_keep, ]
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
    cell_metadata <- cbind.data.frame(cell_metadata, srt@meta.data[rownames(cell_metadata), cell_annotation, drop = FALSE])
  }

  feature_metadata <- data.frame(row.names = feature_union, features = feature_union)
  for (l in lineages) {
    feature_metadata[rownames(dynamic[[l]]), paste0(l, order_by)] <- dynamic[[l]][, order_by]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "exp_ncells")] <- dynamic[[l]][, "exp_ncells"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "r.sq")] <- dynamic[[l]][, "r.sq"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "dev.expl")] <- dynamic[[l]][, "dev.expl"]
    feature_metadata[rownames(dynamic[[l]]), paste0(l, "padjust")] <- dynamic[[l]][, "padjust"]
  }
  feature_metadata[, order_by] <- apply(feature_metadata[, paste0(lineages, order_by), drop = FALSE], 1, max, na.rm = TRUE)
  feature_metadata <- feature_metadata[order(feature_metadata[, order_by], decreasing = decreasing), ]
  feature_metadata <- feature_metadata[rownames(feature_metadata) %in% feature_union, ]
  feature_union <- rownames(feature_metadata)
  if (!is.null(feature_annotation)) {
    feature_metadata <- cbind.data.frame(feature_metadata, srt[[assay]]@meta.features[rownames(feature_metadata), feature_annotation, drop = FALSE])
  }

  if (isTRUE(use_fitted)) {
    mat_list <- list()
    for (l in lineages) {
      fitted_matrix <- NULL
      if (paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
        fitted_matrix <- srt@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, -1]
      }
      feature_calcu <- feature_union[!feature_union %in% colnames(fitted_matrix)]
      if (length(feature_calcu) > 0) {
        srt_tmp <- RunDynamicFeatures(srt, lineages = l, features = feature_calcu, assay = assay, slot = slot, family = family, libsize = libsize)
        if (is.null(fitted_matrix)) {
          fitted_matrix <- matrix(NA, nrow = nrow(srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]), ncol = 0)
        }
        fitted_matrix <- cbind(fitted_matrix, srt_tmp@tools[[paste0("DynamicFeatures_", l)]][["fitted_matrix"]][, feature_calcu, drop = FALSE])
      }
      rownames(fitted_matrix) <- paste0(rownames(fitted_matrix), l)
      mat_list[[l]] <- t(fitted_matrix[, feature_union])
    }
    mat <- do.call(cbind, mat_list)
  } else {
    mat_list <- list()
    Y_libsize <- colSums(GetAssayData(srt, slot = "counts", assay = assay))
    for (l in lineages) {
      cells <- gsub(pattern = l, replacement = "", x = cell_order_list[[l]])
      mat_tmp <- as.matrix(rbind(GetAssayData(srt, assay = assay, slot = slot)[gene, cells], t(srt[[meta]])[, cells]))[feature_union, ]
      if (isTRUE(lib_normalize) && min(mat_tmp, na.rm = TRUE) >= 0) {
        if (!is.null(libsize)) {
          libsize_use <- libsize
        } else {
          libsize_use <- Y_libsize[colnames(mat_tmp)]
          isfloat <- any(libsize_use %% 1 != 0, na.rm = TRUE)
          if (isTRUE(isfloat)) {
            libsize_use <- rep(1, length(libsize_use))
            warning("Values in 'counts' slot is non-integer. Set the libsize to 1.", immediate. = TRUE)
          }
        }
        mat_tmp <- t(t(mat_tmp) / libsize_use * median(Y_libsize))
      }
      colnames(mat_tmp) <- paste0(colnames(mat_tmp), l)
      mat_list[[l]] <- mat_tmp
    }
    mat <- do.call(cbind, mat_list)
  }
  mat_split <- mat[, unlist(cell_order_list[heatmap_split_by]), drop = FALSE]

  if (exp_method == "raw") {
    mat <- mat
  } else if (exp_method == "zscore") {
    mat <- t(scale(t(mat)))
  } else if (exp_method == "log2fc") {
    mat <- log2(mat / rowMeans(mat))
  } else if (exp_method == "log1p") {
    mat <- log1p(mat)
  }
  mat[is.infinite(mat)] <- max(abs(mat[!is.infinite(mat)])) * ifelse(mat[is.infinite(mat)] > 0, 1, -1)
  mat[is.na(mat)] <- mean(mat, na.rm = TRUE)

  if ((!identical(sort(heatmap_split_by), sort(lineages)) && is.null(feature_split) && n_split < nrow(mat) && n_split > 1) || cell_density != 1) {
    if (exp_method == "raw") {
      mat_split <- mat_split
    } else if (exp_method == "zscore") {
      mat_split <- t(scale(t(mat_split)))
    } else if (exp_method == "log2fc") {
      mat_split <- log2(mat_split / rowMeans(mat_split))
    } else if (exp_method == "log1p") {
      mat_split <- log1p(mat_split)
    }
    mat_split[is.infinite(mat_split)] <- max(abs(mat_split[!is.infinite(mat_split)])) * ifelse(mat_split[is.infinite(mat_split)] > 0, 1, -1)
    mat_split[is.na(mat_split)] <- mean(mat_split, na.rm = TRUE)
  } else {
    mat_split <- mat[, unlist(cell_order_list[heatmap_split_by]), drop = FALSE]
  }

  if (exp_method %in% c("zscore", "log2fc")) {
    b <- ceiling(min(abs(quantile(mat, c(0.01, 0.99), na.rm = TRUE))) * 2) / 2
    colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = heatmap_palette))
  } else if (exp_method %in% c("raw", "log1p")) {
    b <- quantile(mat, c(0.01, 0.99), na.rm = TRUE)
    colors <- colorRamp2(seq(b[1], b[2], length = 100), palette_scp(palette = heatmap_palette))
  }

  lgd <- list()
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_top_list <- list()
  pseudotime <- na.omit(unlist(cell_metadata[, lineages]))
  pseudotime_col <- colorRamp2(
    breaks = seq(min(pseudotime), max(pseudotime), length = 100),
    colors = palette_scp(palette = pseudotime_palette)
  )
  for (l in lineages) {
    ha_top_list[[l]] <- HeatmapAnnotation(
      Pseudotime = anno_simple(
        x = cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l],
        col = pseudotime_col,
        border = TRUE
      ), which = "column", show_annotation_name = FALSE
    )
  }
  lgd[["pseudotime"]] <- Legend(title = "Pseudotime", col_fun = pseudotime_col, border = TRUE)

  if (!is.null(cell_annotation)) {
    for (i in seq_along(cell_annotation)) {
      cellan <- cell_annotation[i]
      palette <- cell_palette[i]
      palcolor <- cell_palcolor[[i]]
      cell_anno <- cell_metadata[cell_union, cellan]
      names(cell_anno) <- cell_union
      if (!is.factor(cell_anno)) {
        cell_anno <- factor(cell_anno, levels = unique(cell_anno))
      }
      for (l in lineages) {
        ha_top_tmp <- list()
        ha_top_tmp[[cellan]] <- anno_simple(
          x = as.character(cell_anno[gsub(pattern = l, replacement = "", x = cell_order_list[[l]])]),
          col = palette_scp(cell_anno, palette = palette, palcolor = palcolor),
          border = TRUE
        )
        ha_top <- do.call("HeatmapAnnotation", args = c(ha_top_tmp, which = "column", show_annotation_name = FALSE))
        ha_top_list[[l]] <- c(ha_top_list[[l]], ha_top)
      }
      lgd[[cellan]] <- Legend(
        title = cellan, labels = levels(cell_anno),
        legend_gp = gpar(fill = palette_scp(cell_anno, palette = palette, palcolor = palcolor)), border = TRUE
      )
    }
  }

  if (is.null(feature_split)) {
    if (is.null(n_split) || isTRUE(nrow(mat) <= n_split)) {
      feature_split <- row_split <- NULL
    } else {
      if (n_split == 1) {
        feature_split <- row_split <- setNames(rep(1, nrow(mat_split)), rownames(mat_split))
      } else {
        if (split_method == "mfuzz") {
          check_R("Mfuzz")
          require("Mfuzz")
          eset <- new("ExpressionSet", exprs = mat_split)
          eset <- Mfuzz::standardise(eset)
          min_fuzzification <- Mfuzz::mestimate(eset)
          if (is.null(fuzzification)) {
            fuzzification <- min_fuzzification + 0.1
          } else {
            if (fuzzification <= min_fuzzification) {
              warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2))
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
          hc <- hclust(dist(mat_split))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
        if (split_method == "hclust-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          hc <- hclust(dist(feature_y))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }

      df <- data.frame(row_split = row_split, order_by = feature_metadata[names(row_split), order_by])
      df_order <- aggregate(df, by = list(row_split), FUN = mean)
      df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), ]
      split_levels <- c()
      for (i in seq_len(nrow(df_order))) {
        raw_nm <- df_order[i, "row_split"]
        feature_split[feature_split == raw_nm] <- paste0("C", i)
        level <- paste0("C", i, "(", sum(row_split == raw_nm), ")")
        row_split[row_split == raw_nm] <- level
        split_levels <- c(split_levels, level)
      }
      row_split <- factor(row_split, levels = split_levels)
      feature_split <- factor(feature_split, levels = paste0("C", seq_len(nrow(df_order))))
    }
  } else {
    feature_split <- feature_split[row.names(mat)]
    row_split <- feature_split[row.names(mat)]
  }
  if (!is.null(feature_split)) {
    feature_metadata[["feature_split"]] <- feature_split[rownames(feature_metadata)]
    feature_metadata <- feature_metadata[order(feature_metadata[["feature_split"]]), ]
    feature_metadata[["index"]] <- seq_len(nrow(feature_metadata))
  }

  if (is.null(features_label)) {
    if (nlabel > 0) {
      if (length(feature_union) > nlabel) {
        index <- seq(floor(length(feature_union) / nlabel) - 1, length(feature_union), ceiling(length(feature_union) / nlabel))
        if (!is.null(row_split)) {
          index <- order(row_split)[index]
        }
      } else {
        index <- seq_along(feature_union)
      }
    } else {
      index <- NULL
    }
  } else {
    index <- which(feature_union %in% features_label)
    drop <- setdiff(features_label, feature_union)
    if (length(drop) > 0) {
      warning(paste0(paste0(drop, collapse = ","), "was not found in the dynamic features"))
    }
  }

  ha_left <- NULL
  if (length(index) > 0) {
    ha_left <- HeatmapAnnotation(
      gene = anno_mark(
        at = index, labels = feature_union[index], side = "left",
        labels_gp = gpar(fontsize = label_size, col = label_color),
        link_gp = gpar(fontsize = label_size, col = label_color)
      ),
      which = "row", show_annotation_name = FALSE
    )
  }
  if (!is.null(row_split)) {
    ha_clusters <- HeatmapAnnotation(
      feat_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
      feat_cluster = anno_block(
        gp = gpar(fill = palette_scp(row_split, type = "discrete", palette = feature_split_palette)),
        width = unit(0.1, "in")
      ),
      which = "row", show_annotation_name = FALSE
    )
    if (is.null(ha_left)) {
      ha_left <- ha_clusters
    } else {
      ha_left <- c(ha_left, ha_clusters)
    }
    lgd[["Cluster"]] <- Legend(
      title = "Cluster", labels = levels(factor(row_split)),
      legend_gp = gpar(fill = palette_scp(row_split, type = "discrete", palette = feature_split_palette)), border = TRUE
    )
  }

  ha_right <- NULL
  if (length(lineages) > 1) {
    ha_list <- list()
    for (l in lineages) {
      ha_list[[l]] <- anno_simple(x = is.na(feature_metadata[rownames(mat), paste0(l, order_by)]) + 0, col = c("0" = "#181830", "1" = "white"), width = unit(0.5, "cm"), which = "row")
    }
    ha_lineage <- do.call("HeatmapAnnotation", args = c(ha_list, which = "row", annotation_name_side = "top", border = TRUE))
    if (is.null(ha_right)) {
      ha_right <- ha_lineage
    } else {
      ha_right <- c(ha_right, ha_lineage)
    }
  }
  if (!is.null(feature_annotation)) {
    for (i in seq_along(feature_annotation)) {
      featurean <- feature_annotation[i]
      palette <- feature_palette[i]
      palcolor <- feature_palcolor[[i]]
      feature_class <- feature_metadata[rownames(mat), featurean]
      if (!is.factor(feature_class)) {
        feature_class <- factor(feature_class, levels = unique(feature_class))
      }
      ha_feature <- list()
      ha_feature[[featurean]] <- anno_simple(
        x = as.character(feature_class),
        col = palette_scp(feature_class, palette = palette, palcolor = palcolor),
        width = unit(0.5, "cm"), which = "row"
      )
      ha_feature <- do.call("HeatmapAnnotation", args = c(ha_feature, which = "row", annotation_name_side = "top", border = TRUE))
      if (is.null(ha_right)) {
        ha_right <- ha_feature
      } else {
        ha_right <- c(ha_right, ha_feature)
      }
      lgd[[featurean]] <- Legend(
        title = featurean, labels = levels(feature_class),
        legend_gp = gpar(fill = palette_scp(feature_class, palette = palette, palcolor = palcolor)), border = TRUE
      )
    }
  }

  res <- NULL
  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    check_R("ggwordcloud")
    check_R("simplifyEnrichment")
    clusters <- row_split %||% setNames(rep(1, nrow(mat_split)), rownames(mat_split))
    res <- RunEnrichment(
      geneID = names(clusters), geneID_groups = clusters, IDtype = IDtype, species = species,
      db_IDtype = db_IDtype, db_update = db_update, db_version = db_version, Ensembl_version = Ensembl_version, mirror = mirror,
      enrichment = enrichment, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize, universe = universe,
      GO_simplify = GO_simplify, GO_simplify_padjustCutoff = GO_simplify_padjustCutoff, simplify_method = simplify_method, simplify_similarityCutoff = simplify_similarityCutoff
    )
    if (nrow(res$enrichment) == 0) {
      stop("No enrichment result found.")
    }
    metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
    pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
    padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

    df <- res$enrichment %>%
      filter(Enrichment %in% enrichment) %>%
      group_by(Enrichment, Groups) %>%
      filter(.data[["pvalue"]] <= pvalueCutoff & .data[["p.adjust"]] <= padjustCutoff) %>%
      arrange(desc(-.data[["pvalue"]])) %>%
      as.data.frame()
    df_list <- split.data.frame(df, ~ Enrichment + Groups)
    df_list <- df_list[lapply(df_list, nrow) > 0]

    for (enrich in enrichment) {
      nm <- strsplit(names(df_list), "\\.")
      subdf_list <- df_list[unlist(lapply(nm, function(x) x[[1]])) %in% enrich]

      ha_terms <- NULL
      if (isTRUE(anno_terms)) {
        ids_list <- lapply(subdf_list, function(df) {
          if (isTRUE(show_id)) {
            ids <- paste(head(df$ID, topTerm), head(df$Description, topTerm))
          } else {
            ids <- head(df$Description, topTerm)
            ids <- paste(toupper(substr(ids, 1, 1)), substr(ids, 2, nchar(ids)), sep = "")
          }
          df_out <- data.frame(keyword = ids)
          df_out[["col"]] <- palette_scp(-log10(head(df[, "p.adjust"], topTerm)), type = "continuous", palette = "Spectral", matched = TRUE)
          df_out[["col"]] <- sapply(df_out[["col"]], function(x) blendcolors(c(x, "black")))
          df_out[["fontsize"]] <- mean(anno_size)
          return(df_out)
        })
        names(ids_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_terms <- HeatmapAnnotation(
          id_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          id_cluster = anno_block(
            gp = gpar(fill = palette_scp(clusters, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          id = anno_textbox(
            align_to = clusters, text = ids_list, max_width = anno_width[1],
            word_wrap = TRUE, add_new_line = TRUE,
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
      }

      ha_keys <- NULL
      if (isTRUE(anno_keys)) {
        term_list <- lapply(subdf_list, function(df) {
          if (df$Enrichment[1] %in% c("GO_BP", "GO_CC", "GO_MF")) {
            df <- simplifyEnrichment::keyword_enrichment_from_GO(df[["ID"]]) %>%
              summarise(
                keyword = .data[["keyword"]],
                score = -(log10(.data[["padj"]])),
                count = .data[["n_term"]],
                Enrichment = df[["Enrichment"]][1],
                Groups = df[["Groups"]][1]
              ) %>%
              filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
              filter(!.data[["keyword"]] %in% exclude_words) %>%
              distinct() %>%
              mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
              as.data.frame()
            df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          } else {
            df <- df %>%
              mutate(keyword = strsplit(as.character(.data[["Description"]]), " ")) %>%
              unnest(cols = "keyword") %>%
              group_by(.data[["keyword"]], Enrichment, Groups) %>%
              summarise(
                keyword = .data[["keyword"]],
                score = sum(-(log10(.data[[metric]]))),
                count = n(),
                Enrichment = .data[["Enrichment"]],
                Groups = .data[["Groups"]],
                .groups = "keep"
              ) %>%
              filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
              filter(!.data[["keyword"]] %in% exclude_words) %>%
              distinct() %>%
              mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
              as.data.frame()
            df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          }
          df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
          df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
          df[["fontsize"]] <- rescale(df[, "count"], to = anno_size)
          return(df)
        })
        names(term_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_keys <- HeatmapAnnotation(
          terms_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          terms_cluster = anno_block(
            gp = gpar(fill = palette_scp(clusters, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          terms = anno_textbox(
            align_to = clusters, text = term_list, max_width = anno_width[2],
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
      }

      ha_features <- NULL
      if (isTRUE(anno_features)) {
        features_list <- lapply(subdf_list, function(df) {
          df <- df %>%
            mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
            unnest(cols = "keyword") %>%
            group_by(.data[["keyword"]], Enrichment, Groups) %>%
            summarise(
              keyword = .data[["keyword"]],
              score = sum(-(log10(.data[[metric]]))),
              count = n(),
              Enrichment = .data[["Enrichment"]],
              Groups = .data[["Groups"]],
              .groups = "keep"
            ) %>%
            filter(!.data[["keyword"]] %in% exclude_words) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
          df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
          df[["col"]] <- palette_scp(df[, "score"], type = "continuous", palette = "Spectral", matched = TRUE)
          df[["col"]] <- sapply(df[["col"]], function(x) blendcolors(c(x, "black")))
          df[["fontsize"]] <- rescale(df[, "count"], to = anno_size)
          return(df)
        })
        names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
        ha_features <- HeatmapAnnotation(
          feat_empty = anno_empty(width = unit(0.05, "in"), border = FALSE),
          feat_cluster = anno_block(
            gp = gpar(fill = palette_scp(clusters, type = "discrete", palette = feature_split_palette)),
            width = unit(0.1, "in")
          ),
          feat = anno_textbox(
            align_to = clusters, text = features_list, max_width = anno_width[3],
            background_gp = gpar(fill = "grey98", col = "black"), round_corners = TRUE
          ),
          which = "row", gap = unit(0, "points")
        )
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

  if (is.null(use_raster)) {
    use_raster <- ifelse(max(sapply(cell_order_list, length)) * length(feature_union) > 1e7, TRUE, FALSE)
  }
  ht_list <- c()
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
    ht_list <- ht_list + Heatmap(
      name = l,
      matrix = mat[, cell_order_list[[l]], drop = FALSE],
      col = colors,
      column_title = l,
      row_split = row_split,
      row_title_gp = gpar(fontsize = row_title_size),
      cluster_rows = FALSE,
      cluster_row_slices = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      top_annotation = ha_top_list[[l]],
      left_annotation = left_annotation,
      right_annotation = right_annotation,
      show_heatmap_legend = FALSE,
      border = TRUE,
      use_raster = use_raster,
      raster_device = "png"
    )
  }

  if (length(index) == 0 && is.null(anno_keys) && is.null(anno_features) && is.null(width) && is.null(height)) {
    fix <- FALSE
  } else {
    fix <- TRUE
  }
  if (is.null(height)) {
    height <- max(convertHeight(unit(1, "npc"), units, valueOnly = TRUE), 7)
  }
  if (length(ha_top_list) > 0) {
    height_top <- c()
    for (ha_top in ha_top_list) {
      height_top <- max(height.HeatmapAnnotation(ha_top), height_top)
    }
    height <- as.numeric(convertUnit(unit(height, units = units) + height_top, units))
  }
  if (is.null(width)) {
    width <- max(convertWidth(unit(1, "npc"), units, valueOnly = TRUE), 7)
  }
  if (!is.null(ha_left)) {
    width <- as.numeric(convertUnit(unit(width, units = units) + width.HeatmapAnnotation(ha_left), units))
  }
  if (!is.null(ha_right)) {
    width <- as.numeric(convertUnit(unit(width, units = units) + width.HeatmapAnnotation(ha_right), units))
  }

  gTree <- grid.grabExpr(
    {
      draw(ht_list,
        annotation_legend_list = lgd,
        padding = unit(c(1, 1, 1, 1), "cm") # bottom, left, top and right
      )
      # list_components()
      if (is.numeric(pseudotime_label)) {
        for (n in seq_along(pseudotime_label)) {
          pse <- pseudotime_label[n]
          col <- pseudotime_label_color[n]
          lty <- pseudotime_label_linetype[n]
          lwd <- pseudotime_label_linewidth[n]
          for (l in lineages) {
            for (slice in 1:nlevels(row_split)) {
              decorate_heatmap_body(l,
                {
                  pseudotime <- cell_metadata[gsub(pattern = l, replacement = "", x = cell_order_list[[l]]), l]
                  i <- which.min(abs(pseudotime - pse))
                  x <- i / length(pseudotime)
                  grid.lines(c(x, x), c(0, 1), gp = gpar(lty = lty, lwd = lwd, col = col))
                },
                slice = slice
              )
            }
          }
        }
      }
    },
    height = height,
    width = width
  )
  if (!isTRUE(fix)) {
    p <- plot_grid(gTree)
  } else {
    p <- plot_grid(panel_fix_single(plot_grid(gTree), width = width, height = height, units = units))
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
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous ggplot_build theme geom_point aes scale_fill_identity facet_null ggplotGrob
#' @importFrom cowplot get_legend ggdraw draw_grob
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom grid grob
#' @importFrom rlang  %||%
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

  p1 <- do.call(ClassDimPlot, args = c(
    srt = srt_ref, reduction = ref_reduction, group.by = ref_group,
    ref_param
  )) +
    guides(color = guide_legend(title = paste0("Ref: ", ref_group), override.aes = list(size = 4.5)))
  p1legend <- get_legend(p1)
  # p1legend <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")

  suppressMessages(expr = {
    p2 <- do.call(ClassDimPlot, args = c(
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
    override.aes = list(size = 4.5, shape = 21, color = "black", fill = color[levels(p2$data$group.by)])
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
    grob <- ggplotGrob(p3)
    grob <- gtable_add_cols(grob, sum(legend$widths), -1)
    grob <- gtable_add_grob(grob, legend, t = grob$layout[grob$layout$name == "panel", "t"], l = dim(grob)[2])
    p <- plot_grid(grob)
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
#' @param enrichment
#' @param base_size
#' @param character_width
#' @param line_height
#' @param panel_fix
#' @param panel_height_scale
#' @param panel_width
#' @param align
#' @param axis
#' @param srt
#' @param group_by
#' @param test.use
#' @param res
#' @param group_use
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, group_by = "CellType", enrichment = "GO_BP", species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, group_by = "CellType", plot_type = "bar")
#' EnrichmentPlot(pancreas_sub, group_by = "CellType", plot_type = "bar", character_width = 30, base_size = 10)
#'
#' EnrichmentPlot(pancreas_sub, group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "lollipop", ncol = 1)
#' EnrichmentPlot(pancreas_sub, group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "wordcloud")
#' EnrichmentPlot(pancreas_sub, group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "wordcloud", word_type = "feature")
#' @importFrom ggplot2 ggplot geom_bar geom_text labs scale_fill_manual scale_y_continuous facet_grid coord_flip scale_color_gradientn scale_fill_gradientn scale_size guides geom_segment expansion guide_colorbar
#' @importFrom scales breaks_extended
#' @importFrom dplyr group_by filter arrange desc across mutate summarise distinct n .data
#' @importFrom stringr str_extract str_wrap str_split
#' @importFrom stats formula
#' @export
#'
EnrichmentPlot <- function(srt, enrichment = "GO_BP", group_by = NULL, group_use = NULL, test.use = "wilcox",
                           res = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                           plot_type = c("bar", "lollipop", "wordcloud"), palette = "Spectral",
                           topTerm = 6, topWord = 100, word_type = c("term", "feature"), word_size = c(2, 8),
                           min_word_length = 3, exclude_words = c("cell", "cellular", "dna", "rna", "protein", "development", "system", "regulation", "positive", "negative", "response", "process"),
                           aspect.ratio = 1, base_size = 12, character_width = 50, line_height = 0.7,
                           legend.position = "right", legend.direction = "vertical",
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, align = "hv", axis = "lr") {
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
  if (!is.factor(res[["Enrichment"]])) {
    res[["Enrichment"]] <- factor(res[["Enrichment"]], levels = unique(res[["Enrichment"]]))
  }
  if (!is.factor(res["Groups"])) {
    res[["Groups"]] <- factor(res[["Groups"]], levels = unique(res[["Groups"]]))
  }
  if (length(enrichment[!enrichment %in% res[["Enrichment"]]]) > 0) {
    stop(paste0(enrichment[!enrichment %in% res[["Enrichment"]]], " is not in the enrichment result."))
  }
  if (!is.null(group_use)) {
    res <- res[res[["Groups"]] %in% group_use, ]
  }

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

  df <- res %>%
    filter(Enrichment %in% enrichment) %>%
    group_by(Enrichment, Groups) %>%
    filter(.data[["pvalue"]] <= pvalueCutoff & .data[["p.adjust"]] <= padjustCutoff) %>%
    arrange(desc(-.data[["pvalue"]])) %>%
    as.data.frame()
  df_list <- split.data.frame(df, ~ Enrichment + Groups)
  df_list <- df_list[lapply(df_list, nrow) > 0]

  if (plot_type == "bar") {
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df <- df[head(order(df[[metric]]), topTerm), ]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = df[order(df[, "metric"]), "Description"])

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["metric"]],
        fill = .data[["Enrichment"]],
        label = str_extract(.data[["GeneRatio"]], pattern = "\\d+(?=\\/)")
      )) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(hjust = -0.5, size = 3.5 * base_size / 12) +
        labs(x = "", y = paste0("-log10(", metric, ")")) +
        scale_fill_manual(
          name = "Enrichment:",
          values = palette_scp(levels(df[["Enrichment"]]), palette = palette),
          na.value = "grey80",
          guide = "none"
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(df[["metric"]])), expand = expansion(0, 0)) +
        facet_grid(Enrichment ~ Groups, scales = "free") +
        theme_scp(
          aspect.ratio = aspect.ratio,
          base_size = base_size,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, size = 1),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.position = legend.position,
          legend.direction = legend.direction,
          plot.margin = margin(0, 0, 0, 0),
          axis.text.y = element_text(
            lineheight = line_height, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        ) +
        coord_flip()
      return(p)
    }))
  } else if (plot_type == "lollipop") {
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df <- df[head(order(df[[metric]]), topTerm), ]
      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- str_split(string = x, pattern = "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df[["BgRatio"]] <- sapply(df[["BgRatio"]], function(x) {
        sp <- str_split(string = x, pattern = "/")[[1]]
        BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
        return(BgRatio)
      })
      df[["EnrichmentScore"]] <- df[["GeneRatio"]] / df[["BgRatio"]]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = df[order(df[, "EnrichmentScore"]), "Description"])

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["EnrichmentScore"]],
        fill = .data[["metric"]]
      )) +
        geom_point(aes(size = .data[["Count"]]), shape = 21, color = "white", stroke = 1) +
        geom_segment(aes(yend = 0, xend = .data[["Description"]], color = .data[["metric"]]), size = 2, lineend = "butt") +
        scale_size(name = "Count", range = c(3, 6), breaks = ceiling(seq(min(df[["Count"]]), max(df[["Count"]]), length.out = 3))) +
        guides(size = guide_legend(override.aes = list(colour = "black", shape = 16), order = 1)) +
        scale_y_continuous(limits = c(0, 1.2 * max(df[["EnrichmentScore"]])), expand = expansion(0, 0)) +
        labs(x = "", y = "Enrichment Score") +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          breaks = breaks_extended(n = 3),
          colors = c("gold", "red3"),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4)
        ) +
        scale_color_gradientn(
          name = paste0("-log10(", metric, ")"),
          breaks = breaks_extended(n = 3),
          colors = c("gold", "red3"),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barheight = 4)
        ) +
        facet_grid(Enrichment ~ Groups, scales = "free") +
        theme_scp(
          aspect.ratio = aspect.ratio,
          panel.background = element_rect(fill = "#1A365C"),
          base_size = base_size,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          strip.background.y = element_rect(fill = "white", color = "black", linetype = 1, size = 1),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.position = legend.position,
          legend.direction = legend.direction,
          plot.margin = margin(0, 0, 0, 0),
          axis.text.y = element_text(
            lineheight = line_height, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        ) +
        coord_flip()
      return(p)
    }))
  } else if (plot_type == "wordcloud") {
    check_R("ggwordcloud")
    check_R("simplifyEnrichment")
    plist <- lapply(df_list, function(df) {
      if (word_type == "term") {
        if (df$Enrichment[1] %in% c("GO_BP", "GO_CC", "GO_MF")) {
          df <- simplifyEnrichment::keyword_enrichment_from_GO(df[["ID"]]) %>%
            summarise(
              keyword = .data[["keyword"]],
              score = -(log10(.data[["padj"]])),
              count = .data[["n_term"]],
              Enrichment = df[["Enrichment"]][1],
              Groups = df[["Groups"]][1]
            ) %>%
            filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
            filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
          df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
        } else {
          df <- df %>%
            mutate(keyword = strsplit(as.character(.data[["Description"]]), " ")) %>%
            unnest(cols = "keyword") %>%
            group_by(.data[["keyword"]], Enrichment, Groups) %>%
            summarise(
              keyword = .data[["keyword"]],
              score = sum(-(log10(.data[[metric]]))),
              count = n(),
              Enrichment = .data[["Enrichment"]],
              Groups = .data[["Groups"]],
              .groups = "keep"
            ) %>%
            filter(nchar(.data[["keyword"]]) >= min_word_length) %>%
            filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
            distinct() %>%
            mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
            as.data.frame()
          df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
        }
      } else {
        df <- df %>%
          mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Enrichment, Groups) %>%
          summarise(
            keyword = .data[["keyword"]],
            score = sum(-(log10(.data[[metric]]))),
            count = n(),
            Enrichment = .data[["Enrichment"]],
            Groups = .data[["Groups"]],
            .groups = "keep"
          ) %>%
          filter(!tolower(.data[["keyword"]]) %in% tolower(exclude_words)) %>%
          distinct() %>%
          mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
          as.data.frame()
        df <- df[head(order(df[["score"]], decreasing = TRUE), topWord), ]
      }
      colors <- palette_scp(df[, "score"], type = "continuous", palette = palette, matched = FALSE)
      colors_value <- seq(min(df[, "score"], na.rm = TRUE), quantile(df[, "score"], 0.99, na.rm = TRUE) + 0.001, length.out = 100)
      p <- ggplot(df, aes(label = .data[["keyword"]], size = .data[["count"]], color = .data[["score"]], angle = .data[["angle"]])) +
        ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, eccentricity = 1, shape = "square", show.legend = TRUE, grid_margin = 3) +
        scale_color_gradientn(
          name = "Score:", colours = colors, values = rescale(colors_value),
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        scale_size(name = "Count", range = word_size, breaks = ceiling(seq(min(df[["count"]]), max(df[["count"]]), length.out = 3))) +
        guides(size = guide_legend(override.aes = list(colour = "black", label = "G"), order = 1)) +
        facet_grid(Enrichment ~ Groups, scales = "free") +
        coord_flip() +
        theme_scp(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    })
  }

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow, align = align, axis = axis)
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
#' @param title
#' @param subtitle
#' @param base_size
#' @param rel_heights
#' @param subplots
#' @param ES_geom
#' @param color
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
#' @param group_use
#' @param term_use
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
#' pancreas_sub <- RunGSEA(pancreas_sub, group_by = "CellType", enrichment = "GO_BP", species = "Mus_musculus")
#' GSEAPlot(pancreas_sub, group_by = "CellType", group_use = c("Ngn3 low EP", "Endocrine"))
#' GSEAPlot(pancreas_sub, group_by = "CellType", group_use = "Ductal", geneSetID = "GO:0006412")
#' GSEAPlot(pancreas_sub, group_by = "CellType", group_use = "Endocrine", geneSetID = c("GO:0046903", "GO:0015031", "GO:0007600"))
#' @importFrom ggplot2 ggplot aes_ theme theme_classic alpha element_blank element_rect margin geom_line geom_point geom_rect geom_linerange geom_hline geom_vline geom_segment annotate ggtitle labs xlab ylab scale_x_continuous scale_y_continuous scale_color_manual guides guide_legend
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRamp
#' @importFrom cowplot plot_grid get_legend
#' @importFrom dplyr case_when filter pull
#' @importFrom stringr str_split
#' @importFrom stats quantile
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid textGrob
#' @export
GSEAPlot <- function(srt, enrichment = "GO_BP", group_by = NULL, group_use = NULL, test.use = "wilcox",
                     res = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                     topTerm = 6, geneSetID = NULL, ES_geom = c("line", "dot"),
                     subplots = 1:3, rel_heights = c(1.5, 0.5, 1), rel_width = 3,
                     size = 1.5, alpha = 1, color = "#6BB82D", palette = "Set1",
                     n_coregene = 10, sample_coregene = FALSE, features_label = NULL,
                     label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                     base_size = 12, title = NULL, subtitle = NULL,
                     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  ES_geom <- match.arg(ES_geom)
  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    slot <- paste("GSEA", group_by, test.use, sep = "_")
    if (!slot %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunGSEA first.")
    }
    res <- srt@tools[[slot]][["results"]]
  } else {
    res <- res[["results"]]
  }
  if (is.null(group_use)) {
    use <- grep(pattern = paste0("-", paste(enrichment, collapse = "|"), "$"), x = names(res))
  } else {
    comb <- expand.grid(group_use, enrichment)
    use <- names(res)[names(res) %in% paste(comb$Var1, comb$Var2, sep = "-")]
  }
  if (length(use) == 0) {
    stop(paste0(enrichment, " is not in the enrichment result."))
  }
  res <- res[use]

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

  plist <- suppressWarnings(lapply(names(res), function(nm) {
    res_enrich <- res[[nm]]
    if (is.null(geneSetID)) {
      geneSetID_filter <- res_enrich@result %>%
        filter(.data[["pvalue"]] <= pvalueCutoff & .data[["p.adjust"]] <= padjustCutoff) %>%
        arrange(.data[[metric]])
      geneSetID_up <- geneSetID_filter %>%
        filter(.data[["enrichmentScore"]] >= 0)
      geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
      geneSetID_down <- geneSetID_filter %>%
        filter(.data[["enrichmentScore"]] < 0)
      geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
      geneSetID_use <- head(c(geneSetID_up, geneSetID_down), topTerm)
    } else {
      geneSetID_use <- geneSetID
    }
    if (length(geneSetID_use) == 1) {
      gsdata <- gsInfo(object = res_enrich, geneSetID = geneSetID_use)
    } else {
      gsdata <- do.call(rbind, lapply(geneSetID_use, gsInfo, object = res_enrich))
    }
    padj <- res_enrich[geneSetID_use, c("Description", "p.adjust", "enrichmentScore")]
    rownames(padj) <- padj[, "Description"]
    padj$p.sig <- case_when(
      padj$p.adjust > 0.05 ~ "ns  ",
      padj$p.adjust <= 0.05 & padj$p.adjust > 0.01 ~ "*   ",
      padj$p.adjust <= 0.01 & padj$p.adjust > 0.001 ~ "**  ",
      padj$p.adjust <= 0.001 & padj$p.adjust > 0.0001 ~ "*** ",
      padj$p.adjust <= 0.0001 ~ "****"
    )
    gsdata$p.sig <- padj[gsdata$Description, "p.sig"]
    gsdata$DescriptionP <- paste(gsdata$p.sig, gsdata$Description)
    gsdata$DescriptionP <- factor(gsdata$DescriptionP, levels = unique(gsdata$DescriptionP))
    p <- ggplot(gsdata, aes_(x = ~x)) +
      xlab(NULL) +
      theme_classic(base_size) +
      theme(
        panel.grid.major = element_line(colour = "grey90", linetype = 2),
        panel.grid.minor = element_line(colour = "grey90", linetype = 2)
      ) +
      scale_x_continuous(expand = c(0.01, 0))
    if (ES_geom == "line") {
      es_layer <- geom_line(aes_(y = ~runningScore, color = ~DescriptionP),
        size = size, alpha = alpha
      )
    } else {
      es_layer <- geom_point(aes_(y = ~runningScore, color = ~DescriptionP),
        size = size, alpha = alpha, data = subset(gsdata, position == 1)
      )
    }
    bg_dat <- data.frame(xmin = -Inf, xmax = Inf, ymin = c(0, -Inf), ymax = c(Inf, 0), fill = c(alpha("#C40003", 0.2), alpha("#1D008F", 0.2)))
    p1 <- p +
      geom_rect(data = bg_dat, mapping = aes_(xmin = ~xmin, xmax = ~xmax, ymin = ~ymin, ymax = ~ymax, fill = ~ I(fill)), inherit.aes = FALSE) +
      geom_hline(yintercept = 0, linetype = 1, color = "grey40") +
      es_layer +
      ylab("Enrichment Score") +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1),
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
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
      geom_linerange(aes_(
        ymin = ~ymin,
        ymax = ~ymax, color = ~DescriptionP
      ), alpha = alpha) +
      xlab(NULL) +
      ylab(NULL) +
      theme_classic(base_size) +
      theme(
        legend.position = "none",
        plot.margin = margin(t = -0.1, b = 0, r = 0.2, l = 0.2, unit = "cm"),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1),
        axis.line.y = element_blank(), axis.line.x = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()
      ) +
      scale_x_continuous(expand = c(0.01, 0)) +
      scale_y_continuous(expand = c(0, 0))

    subtitle_use <- subtitle
    if (length(geneSetID_use) == 1) {
      subtitle_use <- paste0(ifelse(is.null(subtitle), "(ES=", paste0(subtitle, "   (ES=")), round(padj$enrichmentScore, 3), ", Sig=", padj$p.sig, ")")
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
          fill = ifelse(padj$enrichmentScore < 0, "#5E34F5", "#F52323"), color = "black", size = 2.5,
          shape = ifelse(padj$enrichmentScore < 0, 25, 24)
        ) +
        labs(subtitle = subtitle_use) +
        theme(plot.subtitle = element_text(face = "italic"))

      if ((is.numeric(n_coregene) && n_coregene > 1) || length(features_label) > 0) {
        if (length(features_label) == 0) {
          features_label <- unlist(str_split(gsdata$CoreGene[1], pattern = "/"))
          n_coregene <- min(n_coregene, length(features_label))
          if (isTRUE(sample_coregene)) {
            features_label <- sample(features_label, n_coregene, replace = FALSE)
          } else {
            features_label <- gsdata$GeneName[gsdata$GeneName %in% features_label][1:n_coregene]
          }
        }
        df_gene <- gsdata[gsdata$position == 1 & gsdata$GeneName %in% features_label, ]
        gene_drop <- features_label[!features_label %in% df_gene$GeneName]
        if (length(gene_drop) > 0) {
          warning("Gene ", paste(gene_drop, collapse = ","), " is not in the geneset of the ", gsdata$Description[1], immediate. = TRUE)
        }
        x_nudge <- diff(range(gsdata$x)) * 0.05
        y_nudge <- diff(range(gsdata$runningScore)) * 0.05
        p1 <- p1 + geom_point(
          data = df_gene,
          mapping = aes_(y = ~runningScore), color = "black"
        ) +
          geom_text_repel(
            data = df_gene,
            mapping = aes_(y = ~runningScore, label = ~GeneName),
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
      y_pos_i <- cut(y[y_pos],
        breaks = seq(min(y[y_pos]), max(y[y_pos]), len = 100),
        include.lowest = TRUE
      )
      y_neg <- which(y < 0)
      y_neg_i <- cut(y[y_neg],
        breaks = seq(min(y[y_neg]), max(y[y_neg]), len = 100),
        include.lowest = TRUE
      )
      col[y_pos] <- colorRampPalette(c("#F5DCDC", "#C40003"))(100)[y_pos_i]
      col[y_neg] <- colorRampPalette(c("#1D008F", "#DDDCF5"))(100)[y_neg_i]
      ymin <- min(p2$data$ymin)
      ymax <- max(p2$data$ymax - p2$data$ymin) * 0.3
      xmin <- which(!duplicated(col))
      xmax <- xmin + as.numeric(table(col)[as.character(unique(col))])
      d <- data.frame(
        ymin = ymin, ymax = ymax, xmin = xmin,
        xmax = xmax, col = unique(col)
      )
      p2 <- p2 + geom_rect(aes_(
        xmin = ~xmin, xmax = ~xmax,
        ymin = ~ymin, ymax = ~ymax, fill = ~ I(col)
      ),
      data = d,
      alpha = 0.95, inherit.aes = FALSE
      )
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    min_y <- df2$y[which.min(abs(df2$y))]
    corss_x <- median(df2$x[df2$y == min_y])
    p3 <- p + geom_segment(data = df2, aes_(
      x = ~x, xend = ~x,
      y = ~y, yend = 0
    ), color = "grey30") +
      geom_vline(xintercept = corss_x, linetype = 2, color = "black") +
      annotate(geom = "text", y = 0, x = corss_x, vjust = ifelse(diff(abs(range(df2$y))) > 0, -0.3, 1.3), size = 4, label = paste0("Zero cross at ", corss_x)) +
      annotate(geom = "text", x = 0, y = Inf, vjust = 1.3, hjust = 0, color = "#C81A1F", size = 4, label = " Positively correlated") +
      annotate(geom = "text", x = Inf, y = -Inf, vjust = -0.3, hjust = 1, color = "#3C298C", size = 4, label = "Negtively correlated ")
    p3 <- p3 + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") +
      theme(
        plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"),
        axis.line = element_blank(), axis.line.x = element_blank(),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1)
      )
    if (!is.null(title) && !is.na(title) && title != "") {
      p1 <- p1 + ggtitle(title, subtitle = subtitle_use)
    }
    if (is.null(title) && length(geneSetID_use) == 1) {
      p1 <- p1 + ggtitle(gsdata$Description[1], subtitle = subtitle_use)
    }
    if (length(color) != length(geneSetID_use)) {
      color_use <- palette_scp(levels(gsdata$DescriptionP), palette = palette)
    } else {
      color_use <- color
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
      return(plotlist[[1]] + theme(
        aspect.ratio = rel_heights[subplots] / rel_width,
        plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
      ))
    } else {
      plotlist <- lapply(plotlist[subplots], ggplotGrob)
      rel_heights <- rel_heights[subplots]
      for (i in seq_along(plotlist)) {
        plotlist[[i]] <- panel_fix_single(plotlist[[i]], height = rel_heights[i], units = "null", margin = 0, respect = TRUE)
        plotlist[[i]] <- panel_fix_single(plotlist[[i]], width = rel_width, units = "null", margin = 0, respect = TRUE)
      }
      p_out <- do.call(rbind, c(plotlist, size = "first"))

      if (length(geneSetID_use) > 1) {
        p_out <- gtable_add_rows(p_out, sum(legend$heights), 0)
        p_out <- gtable_add_grob(p_out, legend, t = 1, l = min(p_out$layout[grepl(pattern = "panel", x = p_out$layout$name), "l"]))
      }
      lab <- textGrob(label = nm, rot = -90, hjust = 0.5)
      p_out <- gtable_add_cols(p_out, grobWidth(lab), -1)
      p_out <- gtable_add_grob(p_out, lab, t = mean(p_out$layout[grep("panel", p_out$layout$name), "t"]), l = dim(p_out)[2], clip = FALSE)
      p_out <- plot_grid(p_out)
      return(p_out)
    }
  }))
  names(plist) <- names(res)

  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- cowplot::plot_grid(plotlist = plist, nrow = nrow, ncol = ncol)
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
  max.ES <- max(runningES)
  min.ES <- min(runningES)
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
