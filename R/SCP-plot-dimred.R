CellDimPlot <- function(srt, group.by, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL,
                        show_na = FALSE, show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
                        pt.size = NULL, pt.alpha = 1, palette = "Paired", palcolor = NULL, bg_color = "grey80",
                        label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                        label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                        label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                        cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                        add_density = FALSE, density_color = "grey80", density_filled = FALSE, density_filled_palette = "Greys", density_filled_palcolor = NULL,
                        add_mark = FALSE, mark_type = c("hull", "ellipse", "rect", "circle"), mark_expand = unit(3, "mm"), mark_alpha = 0.1, mark_linetype = 1,
                        lineages = NULL, lineages_trim = c(0.01, 0.99), lineages_span = 0.75,
                        lineages_palette = "Dark2", lineages_palcolor = NULL, lineages_arrow = arrow(length = unit(0.1, "inches")),
                        lineages_linewidth = 1, lineages_line_bg = "white", lineages_line_bg_stroke = 0.5,
                        lineages_whiskers = FALSE, lineages_whiskers_linewidth = 0.5, lineages_whiskers_alpha = 0.5,
                        stat.by = NULL, stat_type = "percent", stat_plot_type = "pie", stat_plot_position = c("stack", "dodge"), stat_plot_size = 0.15,
                        stat_plot_palette = "Set1", stat_palcolor = NULL, stat_plot_alpha = 1, stat_plot_label = FALSE, stat_plot_label_size = 3,
                        graph = NULL, edge_size = c(0.05, 0.5), edge_alpha = 0.1, edge_color = "grey40",
                        paga = NULL, paga_type = "connectivities", paga_node_size = 4,
                        paga_edge_threshold = 0.01, paga_edge_size = c(0.2, 1), paga_edge_color = "grey40", paga_edge_alpha = 0.5,
                        paga_transition_threshold = 0.01, paga_transition_size = c(0.2, 1), paga_transition_color = "black", paga_transition_alpha = 1, paga_show_transition = FALSE,
                        velocity = NULL, velocity_plot_type = "raw", velocity_n_neighbors = ceiling(ncol(srt@assays[[1]]) / 50),
                        velocity_density = 1, velocity_smooth = 0.5, velocity_scale = 1, velocity_min_mass = 1, velocity_cutoff_perc = 5,
                        velocity_arrow_color = "black", velocity_arrow_angle = 20,
                        streamline_L = 5, streamline_minL = 1, streamline_res = 1, streamline_n = 15,
                        streamline_width = c(0, 0.8), streamline_alpha = 1, streamline_color = NULL, streamline_palette = "RdYlBu", streamline_palcolor = NULL,
                        streamline_bg_color = "white", streamline_bg_stroke = 0.5,
                        hex = FALSE, hex.linewidth = 0.5, hex.count = TRUE, hex.bins = 50, hex.binwidth = NULL,
                        raster = NULL, raster.dpi = c(512, 512),
                        aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                        legend.position = "right", legend.direction = "vertical",
                        theme_use = "theme_scp", theme_args = list(),
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)
  mark_type <- match.arg(mark_type)

  if (is.null(split.by)) {
    split.by <- "All.groups"
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
      theme_use = theme_use, theme_args = theme_args,
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
      theme_use = theme_use, theme_args = theme_args,
      return_layer = TRUE
    )
    lineages_layers <- lineages_layers[!names(lineages_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(paga)) {
    if (split.by != "All.groups") {
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
      theme_use = theme_use, theme_args = theme_args,
      return_layer = TRUE
    )
    paga_layers <- paga_layers[!names(paga_layers) %in% c("lab_layer", "theme_layer")]
  }
  if (!is.null(velocity)) {
    if (split.by != "All.groups") {
      stop("velocity can only plot on the non-split data")
    }
    velocity_layers <- VelocityPlot(srt,
      cells = cells,
      reduction = reduction, dims = dims, velocity = velocity, plot_type = velocity_plot_type, group_by = group.by, group_palette = palette, group_palcolor = palcolor,
      n_neighbors = velocity_n_neighbors, density = velocity_density, smooth = velocity_smooth, scale = velocity_scale, min_mass = velocity_min_mass, cutoff_perc = velocity_cutoff_perc,
      arrow_color = velocity_arrow_color, arrow_angle = velocity_arrow_angle,
      streamline_L = streamline_L, streamline_minL = streamline_minL, streamline_res = streamline_res, streamline_n = streamline_n,
      streamline_width = streamline_width, streamline_alpha = streamline_alpha, streamline_color = streamline_color, streamline_palette = streamline_palette, streamline_palcolor = streamline_palcolor,
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
    dat[["x"]] <- dat[[paste0(reduction_key, dims[1])]]
    dat[["y"]] <- dat[[paste0(reduction_key, dims[2])]]
    dat[["group.by"]] <- dat[[g]]
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

    if (isTRUE(add_mark)) {
      mark_fun <- switch(mark_type,
        "ellipse" = "geom_mark_ellipse",
        "hull" = "geom_mark_hull",
        "rect" = "geom_mark_rect",
        "circle" = "geom_mark_circle"
      )
      mark <- list(
        do.call(
          mark_fun,
          list(
            data = dat[!is.na(dat[["group.by"]]), , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]], fill = .data[["group.by"]]),
            expand = mark_expand, alpha = mark_alpha, linetype = mark_linetype, show.legend = FALSE, inherit.aes = FALSE
          ),
        ),
        scale_fill_manual(values = colors[names(labels_tb)]),
        scale_color_manual(values = colors[names(labels_tb)]),
        new_scale_fill(),
        new_scale_color()
      )
    } else {
      mark <- NULL
    }

    if (!is.null(graph)) {
      net_mat <- as_matrix(graph)[rownames(dat), rownames(dat)]
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
          color = density_color, inherit.aes = FALSE, show.legend = FALSE
        )
      }
    } else {
      density <- NULL
    }

    p <- ggplot(dat) +
      mark +
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
    if (split.by != "All.groups") {
      p <- p + facet_grid(. ~ split.by)
    }

    if (isTRUE(raster)) {
      p <- p + scattermore::geom_scattermore(
        data = dat[is.na(dat[, "group.by"]), , drop = FALSE],
        mapping = aes(x = .data[["x"]], y = .data[["y"]]), color = bg_color,
        pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
      ) + scattermore::geom_scattermore(
        data = dat[!is.na(dat[, "group.by"]), , drop = FALSE],
        mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group.by"]]),
        pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
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
        order = 1,
        override.aes = list(size = 4, alpha = 1)
      )
    ) + scale_fill_manual(
      name = paste0(g, ":"),
      values = colors[names(labels_tb)],
      labels = label_use,
      na.value = bg_color,
      guide = guide_legend(
        title.hjust = 0,
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
          point.size = NA, max.overlaps = 100, force = 0,
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
#' @param features A character vector or a named list of features to plot. Features can be gene names in Assay or names of numeric columns in meta.data.
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
#' @seealso \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP")
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", bg_cutoff = -Inf)
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", theme_use = "theme_blank")
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP", theme_use = ggplot2::theme_classic, theme_args = list(base_size = 16))
#' FeatureDimPlot(pancreas_sub, features = "G2M_score", reduction = "UMAP") %>% panel_fix(height = 2, raster = TRUE, dpi = 30)
#'
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' FeatureDimPlot(pancreas_sub, features = c("StandardPC_1", "StandardPC_2"), reduction = "UMAP", bg_cutoff = -Inf)
#'
#' # Label and highlight cell points
#' FeatureDimPlot(pancreas_sub,
#'   features = "Rbp4", reduction = "UMAP", label = TRUE,
#'   cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Delta"]
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = "Rbp4", split.by = "Phase", reduction = "UMAP",
#'   cells.highlight = TRUE, theme_use = "theme_blank"
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
#' # Input a named feature list
#' markers <- list(
#'   "Ductal" = c("Sox9", "Anxa2", "Bicc1"),
#'   "EPs" = c("Neurog3", "Hes6"),
#'   "Pre-endocrine" = c("Fev", "Neurod1"),
#'   "Endocrine" = c("Rbp4", "Pyy"),
#'   "Beta" = "Ins1", "Alpha" = "Gcg", "Delta" = "Sst", "Epsilon" = "Ghrl"
#' )
#' FeatureDimPlot(pancreas_sub,
#'   features = markers, reduction = "UMAP",
#'   theme_use = "theme_blank",
#'   theme_args = list(plot.subtitle = ggplot2::element_text(size = 10), strip.text = ggplot2::element_text(size = 8))
#' )
#'
#' # Plot multiple features with different scales
#' endocrine_markers <- c("Beta" = "Ins1", "Alpha" = "Gcg", "Delta" = "Sst", "Epsilon" = "Ghrl")
#' FeatureDimPlot(pancreas_sub, endocrine_markers, reduction = "UMAP")
#' FeatureDimPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", lower_quantile = 0, upper_quantile = 0.8)
#' FeatureDimPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", lower_cutoff = 1, upper_cutoff = 4)
#' FeatureDimPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", bg_cutoff = 2, lower_cutoff = 2, upper_cutoff = 4)
#' FeatureDimPlot(pancreas_sub, endocrine_markers, reduction = "UMAP", keep_scale = "all")
#' FeatureDimPlot(pancreas_sub, c("Delta" = "Sst", "Epsilon" = "Ghrl"), split.by = "Phase", reduction = "UMAP", keep_scale = "feature")
#'
#' # Plot multiple features on one picture
#' FeatureDimPlot(pancreas_sub,
#'   features = endocrine_markers, pt.size = 1,
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
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d geom_segment labs scale_x_continuous scale_y_continuous scale_size_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colorbar geom_hex stat_summary_hex geom_path scale_linewidth_continuous after_stat
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom gtable gtable_add_cols
#' @importFrom ggrepel geom_text_repel GeomTextRepel
#' @importFrom grid arrow unit
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix t
#' @importFrom methods slot
#' @importFrom reshape2 melt
#' @export
FeatureDimPlot <- function(srt, features, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL, slot = "data", assay = NULL,
                           show_stat = ifelse(identical(theme_use, "theme_blank"), FALSE, TRUE),
                           palette = ifelse(isTRUE(compare_features), "Set1", "Spectral"), palcolor = NULL,
                           pt.size = NULL, pt.alpha = 1, bg_cutoff = 0, bg_color = "grey80",
                           keep_scale = "feature", lower_quantile = 0, upper_quantile = 0.99, lower_cutoff = NULL, upper_cutoff = NULL,
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
  if (!is.null(keep_scale)) {
    keep_scale <- match.arg(keep_scale, choices = c("feature", "all"))
  }

  if (is.list(features)) {
    if (is.null(names(features))) {
      features <- unlist(features)
    } else {
      features <- setNames(unlist(features), nm = rep(names(features), sapply(features, length)))
    }
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All.groups"
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

  feature_input <- features
  features <- unique(features)
  embeddings <- lapply(srt@reductions, function(x) colnames(x@cell.embeddings))
  embeddings <- setNames(rep(names(embeddings), sapply(embeddings, length)), nm = unlist(embeddings))
  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data), names(embeddings))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  features_embedding <- features[features %in% names(embeddings)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names: ", paste0(intersect(features_gene, features_meta), collapse = ","))
  }
  if (length(c(features_gene, features_meta, features_embedding)) == 0) {
    stop("There are no valid features present.")
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
      dat_gene <- t(as_matrix(slot(srt@assays[[assay]], slot)))
    } else {
      dat_gene <- t(as_matrix(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE]))
    }
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_embedding) > 0) {
    dat_embedding <- as_matrix(FetchData(srt, vars = features_embedding))
  } else {
    dat_embedding <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- as_matrix(do.call(cbind, list(dat_gene, dat_meta, dat_embedding)))
  features <- unique(features[features %in% c(features_gene, features_meta, features_embedding)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'features' must be type of numeric variable.")
  }
  dat_exp[, features][dat_exp[, features] <= bg_cutoff] <- NA

  if (length(features) > 50 && !isTRUE(force)) {
    warning("More than 50 features to be plotted", immediate. = TRUE)
    answer <- askYesNo("Are you sure to continue?", default = FALSE)
    if (!isTRUE(answer)) {
      return(invisible(NULL))
    }
  }

  if (is.null(subtitle)) {
    if (!is.null(names(feature_input))) {
      subtitle <- setNames(names(feature_input), nm = feature_input)
    }
  } else {
    if (length(subtitle) == 1) {
      subtitle <- setNames(rep(subtitle, length(features)), nm = features)
    } else if (length(subtitle) == length(features)) {
      subtitle <- setNames(subtitle, nm = features)
    } else {
      stop(paste0("Subtitle length must be 1 or length of features(", length(features), ")"))
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
      aspect.ratio = aspect.ratio, xlab = xlab, ylab = ylab,
      legend.position = legend.position, legend.direction = legend.direction,
      theme_use = theme_use, theme_args = theme_args,
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
      dat <- dat_split[[ifelse(split.by == "All.groups", 1, s)]][, , drop = FALSE]
      for (f in features) {
        if (any(is.infinite(dat[, f]))) {
          dat[, f][which(dat[, f] == max(dat[, f], na.rm = TRUE))] <- max(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
          dat[, f][which(dat[, f] == min(dat[, f], na.rm = TRUE))] <- min(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
        }
      }
      dat[["x"]] <- dat[[paste0(reduction_key, dims[1])]]
      dat[["y"]] <- dat[[paste0(reduction_key, dims[2])]]
      dat[, "split.by"] <- s
      dat[, "features"] <- paste(features, collapse = "|")
      subtitle_use <- paste0(subtitle, collapse = "|") %||% s
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
          geom_point(data = dat, mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[[names(colors)[i]]]))
        )
        if (all(is.na(colors_list[[i]]))) {
          temp_geom[[i]] <- append(temp_geom[[i]], scale_colour_gradient(
            na.value = bg_color,
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
          ))
        } else if (length(colors_list[[i]]) == 1) {
          temp_geom[[i]] <- append(temp_geom[[i]], scale_colour_gradient(
            low = colors_list[[i]], na.value = bg_color,
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
          ))
        } else {
          temp_geom[[i]] <- append(temp_geom[[i]], scale_color_gradientn(
            colours = pal_list[[i]],
            values = rescale(value_list[[i]]), na.value = bg_color,
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
          ))
        }
        legend_list[[i]] <- get_legend(
          ggplot(dat, aes(x = .data[["x"]], y = .data[["y"]])) +
            temp_geom[[i]] +
            do.call(theme_use, theme_args) +
            theme(
              aspect.ratio = aspect.ratio,
              legend.position = legend.position,
              legend.direction = legend.direction
            )
        )
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
        net_mat <- as_matrix(graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
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
        scale_y_continuous(limits = c(min(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE), max(dat_use[, paste0(reduction_key, dims[2])], na.rm = TRUE)))
      if (split.by == "All.groups") {
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
          pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
        ) +
          scattermore::geom_scattermore(
            data = dat[dat[, "color_blend"] != bg_color, , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["color_blend"]]),
            pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
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
        label_df <- melt(p$data, measure.vars = features)
        label_df <- label_df %>%
          group_by(variable) %>%
          filter(value >= quantile(value[is.finite(value)], 0.95, na.rm = TRUE) & value <= quantile(value[is.finite(value)], 0.99, na.rm = TRUE)) %>%
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
              point.size = NA, max.overlaps = 100, force = 0,
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
              point.size = NA, max.overlaps = 100, force = 0,
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
      dat <- dat_split[[ifelse(split.by == "All.groups", 1, s)]][, c(colnames(dat_use), f), drop = FALSE]
      if (any(is.infinite(dat[, f]))) {
        dat[, f][dat[, f] == max(dat[, f], na.rm = TRUE)] <- max(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
        dat[, f][dat[, f] == min(dat[, f], na.rm = TRUE)] <- min(dat[, f][is.finite(dat[, f])], na.rm = TRUE)
      }
      dat[["x"]] <- dat[[paste0(reduction_key, dims[1])]]
      dat[["y"]] <- dat[[paste0(reduction_key, dims[2])]]
      dat[["value"]] <- dat[[f]]
      dat <- dat[order(dat[, "value"], method = "radix", decreasing = FALSE, na.last = FALSE), , drop = FALSE]
      dat[, "features"] <- f
      cells.highlight_use <- cells.highlight
      if (isTRUE(cells.highlight_use)) {
        cells.highlight_use <- rownames(dat)[!is.na(dat[["value"]])]
      }
      legend_list <- list()
      if (isTRUE(show_stat)) {
        subtitle_use <- subtitle[f] %||% paste0(s, " nPos:", sum(dat[["value"]] > 0, na.rm = TRUE), ", ", round(sum(dat[["value"]] > 0, na.rm = TRUE) / nrow(dat) * 100, 2), "%")
      } else {
        subtitle_use <- subtitle[f]
      }
      if (all(is.na(dat[["value"]]))) {
        colors_value <- rep(0, 100)
      } else {
        if (is.null(keep_scale)) {
          colors_value <- seq(lower_cutoff %||% quantile(dat[is.finite(dat[, "value"]), "value"], lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(dat[is.finite(dat[, "value"]), "value"], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
        } else {
          if (keep_scale == "feature") {
            colors_value <- seq(lower_cutoff %||% quantile(dat_exp[is.finite(dat_exp[, f]), f], lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(dat_exp[is.finite(dat_exp[, f]), f], upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
          }
          if (keep_scale == "all") {
            all_values <- as_matrix(dat_exp[, features])
            colors_value <- seq(lower_cutoff %||% quantile(all_values[is.finite(all_values)], lower_quantile, na.rm = TRUE), upper_cutoff %||% quantile(all_values, upper_quantile, na.rm = TRUE) + 0.001, length.out = 100)
          }
        }
      }
      dat[which(dat[, "value"] > max(colors_value, na.rm = TRUE)), "value"] <- max(colors_value, na.rm = TRUE)
      dat[which(dat[, "value"] < min(colors_value, na.rm = TRUE)), "value"] <- min(colors_value, na.rm = TRUE)
      if (!is.null(graph)) {
        net_mat <- as_matrix(graph)[rownames(dat), rownames(dat)]
        net_mat[net_mat == 0] <- NA
        net_mat[upper.tri(net_mat)] <- NA
        net_df <- melt(net_mat, na.rm = TRUE, stringsAsFactors = FALSE)
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
          pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
        ) +
          scattermore::geom_scattermore(
            data = dat[!is.na(dat[, "value"]), , drop = FALSE],
            mapping = aes(x = .data[["x"]], y = .data[["y"]], color = .data[["value"]]),
            pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
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
        if (split.by == "All.groups") {
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
        color = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
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
          filter(value >= quantile(value[is.finite(value)], 0.95, na.rm = TRUE) & value <= quantile(value[is.finite(value)], 0.99, na.rm = TRUE)) %>%
          reframe(x = median(.data[["x"]]), y = median(.data[["y"]])) %>%
          as.data.frame()
        label_df[, "label"] <- f
        label_df[, "rank"] <- seq_len(nrow(label_df))
        if (isTRUE(label_repel)) {
          p <- p + annotate(
            geom = "point", x = label_df[["x"]], y = label_df[["y"]],
            color = "black", size = pt.size + 1
          ) + annotate(
            geom = GeomTextRepel, x = label_df[["x"]], y = label_df[["y"]], label = label_df[["label"]],
            fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
            point.size = pt.size + 1, max.overlaps = 100, force = label_repulsion,
            color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size
          )
        } else {
          p <- p + annotate(
            geom = GeomTextRepel, x = label_df[["x"]], y = label_df[["y"]], label = label_df[["label"]],
            fontface = "bold",
            point.size = NA, max.overlaps = 100, force = 0,
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
#' @param dims Dimensions to plot, must be a three-length numeric vector specifying x-, y- and z-dimensions
#' @param axis_labs A character vector of length 3 indicating the labels for the axes.
#' @param span A numeric value specifying the span of the loess smoother for lineages line.
#' @param shape.highlight Shape of the cell to highlight. See \href{https://plotly.com/r/reference/scattergl/#scattergl-marker-symbol}{scattergl-marker-symbol}
#' @param width Width in pixels, defaults to automatic sizing.
#' @param height Height in pixels, defaults to automatic sizing.
#' @param save The name of the file to save the plot to. Must end in ".html".
#' @seealso \code{\link{CellDimPlot}} \code{\link{FeatureDimPlot3D}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' CellDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#'
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#' CellDimPlot3D(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D", lineages = "Lineage1")
#'
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom utils askYesNo
#' @importFrom plotly plot_ly add_trace layout as_widget
#' @export
CellDimPlot3D <- function(srt, group.by, reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                          palette = "Paired", palcolor = NULL, bg_color = "grey80", pt.size = 1.5,
                          cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                          lineages = NULL, lineages_palette = "Dark2", span = 0.75,
                          width = NULL, height = NULL, save = NULL, force = FALSE) {
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

  p <- plot_ly(data = dat_use, width = width, height = height)
  p <- add_trace(
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
    p <- add_trace(
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
        line = list(width = 6, color = palette_scp(x = lineages, palette = lineages_palette)[l], reverscale = FALSE),
        name = l,
        showlegend = TRUE,
        visible = TRUE
      )
    }
  }

  p <- layout(
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
      widget = as_widget(p),
      file = save
    )
    unlink(gsub("\\.html", "_files", save), recursive = TRUE)
  }

  return(p)
}

#' 3D-Dimensional reduction plot for gene expression visualization.
#'
#' Plotting cell points on a reduced 3D space and coloring according to the gene expression in the cells.
#' @inheritParams FeatureDimPlot
#' @inheritParams CellDimPlot3D
#'
#' @seealso \code{\link{FeatureDimPlot}} \code{\link{CellDimPlot3D}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' FeatureDimPlot3D(pancreas_sub, features = c("Ghrl", "Ins1", "Gcg", "Ins2"), reduction = "StandardpcaUMAP3D")
#' FeatureDimPlot3D(pancreas_sub, features = c("StandardPC_1", "StandardPC_2"), reduction = "StandardpcaUMAP3D")
#'
#' @importFrom Seurat Reductions Embeddings Key FetchData
#' @importFrom Matrix t
#' @importFrom methods slot
#' @importFrom utils askYesNo
#' @importFrom plotly plot_ly add_trace layout as_widget
#' @export
FeatureDimPlot3D <- function(srt, features, reduction = NULL, dims = c(1, 2, 3), axis_labs = NULL,
                             split.by = NULL, slot = "data", assay = NULL,
                             calculate_coexp = FALSE,
                             pt.size = 1.5, cells.highlight = NULL, cols.highlight = "black", shape.highlight = "circle-open", sizes.highlight = 2,
                             width = NULL, height = NULL, save = NULL, force = FALSE) {
  cols.highlight <- col2hex(cols.highlight)

  if (is.list(features)) {
    features <- unlist(features)
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }

  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All.cells"
    srt@meta.data[[split.by]] <- factor("All.cells")
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
  embeddings <- lapply(srt@reductions, function(x) colnames(x@cell.embeddings))
  embeddings <- setNames(rep(names(embeddings), sapply(embeddings, length)), nm = unlist(embeddings))
  features_drop <- features[!features %in% c(rownames(srt@assays[[assay]]), colnames(srt@meta.data), names(embeddings))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not in the features of srt.", immediate. = TRUE)
    features <- features[!features %in% features_drop]
  }

  features_gene <- features[features %in% rownames(srt@assays[[assay]])]
  features_meta <- features[features %in% colnames(srt@meta.data)]
  features_embedding <- features[features %in% names(embeddings)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names: ", paste0(intersect(features_gene, features_meta), collapse = ","))
  }
  if (length(c(features_gene, features_meta, features_embedding)) == 0) {
    stop("There are no valid features present.")
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
      dat_gene <- t(as_matrix(slot(srt@assays[[assay]], slot)))
    } else {
      dat_gene <- t(as_matrix(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE]))
    }
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_embedding) > 0) {
    dat_embedding <- as_matrix(FetchData(srt, vars = features_embedding))
  } else {
    dat_embedding <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  dat_exp <- as_matrix(do.call(cbind, list(dat_gene, dat_meta, dat_embedding)))
  features <- unique(features[features %in% c(features_gene, features_meta, features_embedding)])

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

  p <- plot_ly(data = dat_use, width = width, height = height)
  p <- add_trace(
    p = p,
    data = dat_use,
    x = dat_use[[paste0(reduction_key, dims[1], "All_cells")]],
    y = dat_use[[paste0(reduction_key, dims[2], "All_cells")]],
    z = dat_use[[paste0(reduction_key, dims[3], "All_cells")]],
    text = paste0(
      "Cell:", rownames(dat_use),
      "\nExp:", round(dat_use[[features[1]]], 3)
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
    p <- add_trace(
      p = p,
      x = dat_use_highlight[[paste0(reduction_key, dims[1], "All_cells")]],
      y = dat_use_highlight[[paste0(reduction_key, dims[2], "All_cells")]],
      z = dat_use_highlight[[paste0(reduction_key, dims[3], "All_cells")]],
      text = paste0(
        "Cell:", rownames(dat_use_highlight),
        "\nExp:", round(dat_use_highlight[[features[1]]], 3)
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
    sp <- ifelse(i == 0, "All.cells", levels(dat_use[[split.by]])[i])
    ncells <- ifelse(i == 0, nrow(dat_use), table(dat_use[[split.by]])[sp])
    if (i != 0 && sp == "All.cells") {
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
          "\nExp:", round(dat_use[[features[j]]], 3)
        )),
        marker = marker
      )),
      label = features[j]
    )
  }

  p <- layout(
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
      widget = as_widget(p),
      file = save
    )
    unlink(gsub("\\.html", "_files", save), recursive = TRUE)
  }

  return(p)
}

#' Statistical plot of features
#'
#' This function generates a statistical plot for features.
#'
#' @param srt A Seurat object.
#' @param stat.by A character vector specifying the features to plot.
#' @param group.by A character vector specifying the groups to group by. Default is NULL.
#' @param split.by A character vector specifying the variable to split the plot by. Default is NULL.
#' @param plot.by A character vector specifying how to plot the data, by group or feature. Possible values are "group", "feature". Default is "group".
#' @param bg.by A character vector specifying the variable to use as the background color. Default is NULL.
#' @param fill.by A string specifying what to fill the plot by. Possible values are "group", "feature", or "expression". Default is "group".
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#' @param slot A string specifying which slot of the Seurat object to use. Default is "data".
#' @param assay A string specifying which assay to use. Default is NULL.
#' @param keep_empty A logical indicating whether to keep empty levels in the plot. Default is FALSE.
#' @param individual A logical indicating whether to create individual plots for each group. Default is FALSE.
#' @param plot_type A string specifying the type of plot to create. Possible values are "violin", "box", "bar", "dot", or "col". Default is "violin".
#' @param palette A string specifying the color palette to use for filling. Default is "Paired".
#' @param palcolor A character vector specifying specific colors to use for filling. Default is NULL.
#' @param alpha A numeric value specifying the transparency of the plot. Default is 1.
#' @param bg_palette A string specifying the color palette to use for the background. Default is "Paired".
#' @param bg_palcolor A character vector specifying specific colors to use for the background. Default is NULL.
#' @param bg_alpha A numeric value specifying the transparency of the background. Default is 0.2.
#' @param add_box A logical indicating whether to add a box plot to the plot. Default is FALSE.
#' @param box_color A string specifying the color of the box plot. Default is "black".
#' @param box_width A numeric value specifying the width of the box plot. Default is 0.1.
#' @param box_ptsize A numeric value specifying the size of the points of the box plot. Default is 2.
#' @param add_point A logical indicating whether to add individual data points to the plot. Default is FALSE.
#' @param pt.color A string specifying the color of the data points. Default is "grey30".
#' @param pt.size A numeric value specifying the size of the data points. If NULL, the size is automatically determined. Default is NULL.
#' @param pt.alpha A numeric value specifying the transparency of the data points. Default is 1.
#' @param jitter.width A numeric value specifying the width of the jitter. Default is 0.5.
#' @param jitter.height A numeric value specifying the height of the jitter. Default is 0.1.
#' @param add_trend A logical indicating whether to add a trend line to the plot. Default is FALSE.
#' @param trend_color A string specifying the color of the trend line. Default is "black".
#' @param trend_linewidth A numeric value specifying the width of the trend line. Default is 1.
#' @param trend_ptsize A numeric value specifying the size of the points of the trend line. Default is 2.
#' @param add_stat A string specifying which statistical summary to add to the plot. Possible values are "none", "mean", or "median". Default is "none".
#' @param stat_color A string specifying the color of the statistical summary. Default is "black".
#' @param stat_size A numeric value specifying the size of the statistical summary. Default is 1.
#' @param cells.highlight A logical or character vector specifying the cells to highlight in the plot. If TRUE, all cells are highlighted. If FALSE, no cells are highlighted. Default is NULL.
#' @param cols.highlight A string specifying the color of the highlighted cells. Default is "red".
#' @param sizes.highlight A numeric value specifying the size of the highlighted cells. Default is 1.
#' @param alpha.highlight A numeric value specifying the transparency of the highlighted cells. Default is 1.
#' @param calculate_coexp A logical indicating whether to calculate co-expression values. Default is FALSE.
#' @param same.y.lims A logical indicating whether to use the same y-axis limits for all plots. Default is FALSE.
#' @param y.min A numeric or character value specifying the minimum y-axis limit. If a character value is provided, it must be of the form "qN" where N is a number between 0 and 100 (inclusive) representing the quantile to use for the limit. Default is NULL.
#' @param y.max A numeric or character value specifying the maximum y-axis limit. If a character value is provided, it must be of the form "qN" where N is a number between 0 and 100 (inclusive) representing the quantile to use for the limit. Default is NULL.
#' @param y.trans A string specifying the transformation to apply to the y-axis. Possible values are "identity" or "log2". Default is "identity".
#' @param y.nbreaks An integer specifying the number of breaks to use for the y-axis. Default is 5.
#' @param sort A logical or character value specifying whether to sort the groups on the x-axis. If TRUE, groups are sorted in increasing order. If FALSE, groups are not sorted. If "increasing", groups are sorted in increasing order. If "decreasing", groups are sorted in decreasing order. Default is FALSE.
#' @param stack A logical specifying whether to stack the plots on top of each other. Default is FALSE.
#' @param flip A logical specifying whether to flip the plot vertically. Default is FALSE.
#' @param comparisons A list of length-2 vectors. The entries in the vector are either the names of 2 values on the x-axis or the 2 integers that correspond to the index of the groups of interest, to be compared.
#' @param ref_group A string specifying the reference group for pairwise comparisons. Default is NULL.
#' @param pairwise_method A string specifying the method to use for pairwise comparisons. Default is "wilcox.test".
#' @param multiplegroup_comparisons A logical indicating whether to add multiple group comparisons to the plot. Default is FALSE.
#' @param multiple_method A string specifying the method to use for multiple group comparisons. Default is "kruskal.test".
#' @param sig_label A string specifying the label to use for significant comparisons. Possible values are "p.signif" or "p.format". Default is "p.format".
#' @param sig_labelsize A numeric value specifying the size of the significant comparison labels. Default is 3.5.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot. Default is NULL.
#' @param title A string specifying the title of the plot. Default is NULL.
#' @param subtitle A string specifying the subtitle of the plot. Default is NULL.
#' @param xlab A string specifying the label of the x-axis. Default is NULL.
#' @param ylab A string specifying the label of the y-axis. Default is "Expression level".
#' @param legend.position A string specifying the position of the legend. Possible values are "right", "left", "top", "bottom", or "none". Default is "right".
#' @param legend.direction A string specifying the direction of the legend. Possible values are "vertical" or "horizontal". Default is "vertical".
#' @param theme_use A string specifying the theme to use for the plot. Default is "theme_scp".
#' @param theme_args A list of arguments to pass to the theme function. Default is an empty list.
#' @param combine A logical indicating whether to combine the individual plots into a single plot. Default is TRUE.
#' @param nrow An integer specifying the number of rows for the combined plot. Default is NULL.
#' @param ncol An integer specifying the number of columns for the combined plot. Default is NULL.
#' @param byrow A logical specifying whether to fill the combined plot by row or by column. Default is TRUE.
#' @param force A logical indicating whether to force the plot creation even if there are more than 100 levels in a variable. Default is FALSE.
#' @param seed An integer specifying the random seed to use for generating jitter. Default is 11.
#'
#' @examples
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
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", add_stat = "mean")
#' FeatureStatPlot(pancreas_sub, stat.by = c("G2M_score", "Fev"), group.by = "SubCellType", add_line = 0.2, line_type = 2)
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
#' ) %>% panel_fix_overall(width = 8, height = 5) # As the plot is created by combining, we can adjust the overall height and width directly.
#'
#' FeatureStatPlot(pancreas_sub, stat.by = c("Neurog3", "Rbp4", "Ins1"), group.by = "CellType", plot.by = "group")
#' FeatureStatPlot(pancreas_sub, stat.by = c("Neurog3", "Rbp4", "Ins1"), group.by = "CellType", plot.by = "feature")
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"), group.by = "CellType", plot.by = "feature",
#'   multiplegroup_comparisons = TRUE, sig_label = "p.format", sig_labelsize = 4
#' )
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"), group.by = "CellType", plot.by = "feature",
#'   comparisons = list(c("Neurog3", "Rbp4"), c("Rbp4", "Ins1")),
#'   stack = TRUE
#' )
#' FeatureStatPlot(pancreas_sub, stat.by = c(
#'   "Sox9", "Anxa2", "Bicc1", # Ductal
#'   "Neurog3", "Hes6", # EPs
#'   "Fev", "Neurod1", # Pre-endocrine
#'   "Rbp4", "Pyy", # Endocrine
#'   "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
#' ), group.by = "SubCellType", plot.by = "feature", stack = TRUE)
#'
#' library(Matrix)
#' data <- pancreas_sub@assays$RNA@data
#' pancreas_sub@assays$RNA@scale.data <- as_matrix(data / rowMeans(data))
#' FeatureStatPlot(pancreas_sub,
#'   stat.by = c("Neurog3", "Rbp4", "Ins1"), group.by = "CellType",
#'   slot = "scale.data", ylab = "FoldChange", same.y.lims = TRUE, y.max = 4
#' )
#'
#' @importFrom Seurat FetchData
#' @importFrom reshape2 melt
#' @importFrom gtable gtable_add_cols gtable_add_rows gtable_add_grob gtable_add_padding
#' @importFrom grid grobHeight grobWidth
#' @importFrom patchwork wrap_plots
#' @export
