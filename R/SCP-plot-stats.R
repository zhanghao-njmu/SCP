FeatureStatPlot <- function(srt, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, plot.by = c("group", "feature"), fill.by = c("group", "feature", "expression"),
                            cells = NULL, slot = "data", assay = NULL, keep_empty = FALSE, individual = FALSE,
                            plot_type = c("violin", "box", "bar", "dot", "col"),
                            palette = "Paired", palcolor = NULL, alpha = 1,
                            bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
                            add_box = FALSE, box_color = "black", box_width = 0.1, box_ptsize = 2,
                            add_point = FALSE, pt.color = "grey30", pt.size = NULL, pt.alpha = 1, jitter.width = 0.4, jitter.height = 0.1,
                            add_trend = FALSE, trend_color = "black", trend_linewidth = 1, trend_ptsize = 2,
                            add_stat = c("none", "mean", "median"), stat_color = "black", stat_size = 1, stat_stroke = 1, stat_shape = 25,
                            add_line = NULL, line_color = "red", line_size = 1, line_type = 1,
                            cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1, alpha.highlight = 1,
                            calculate_coexp = FALSE,
                            same.y.lims = FALSE, y.min = NULL, y.max = NULL, y.trans = "identity", y.nbreaks = 5,
                            sort = FALSE, stack = FALSE, flip = FALSE,
                            comparisons = NULL, ref_group = NULL, pairwise_method = "wilcox.test",
                            multiplegroup_comparisons = FALSE, multiple_method = "kruskal.test",
                            sig_label = c("p.signif", "p.format"), sig_labelsize = 3.5,
                            aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = "Expression level",
                            legend.position = "right", legend.direction = "vertical",
                            theme_use = "theme_scp", theme_args = list(),
                            combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  if (is.null(group.by)) {
    group.by <- "All.groups" # avoid having the same name with split.by. split.by will be All.groups by default
    xlab <- "All groups"
    srt[[group.by]] <- factor("All groups")
  }

  meta.data <- srt@meta.data
  meta.data[["cells"]] <- rownames(meta.data)
  assay <- assay %||% DefaultAssay(srt)
  exp.data <- slot(srt@assays[[assay]], slot)
  plot.by <- match.arg(plot.by)

  if (plot.by == "feature") {
    if (length(group.by) > 1) {
      stop("The 'group.by' must have a length of 1 when 'plot.by' is set to 'feature'")
    }
    if (!is.null(bg.by)) {
      message("'bg.by' is invalid when plot.by is set to 'feature'")
    }
    message("Setting 'group.by' to 'Features' as 'plot.by' is set to 'feature'")
    srt@assays[setdiff(names(srt@assays), assay)] <- NULL
    meta.reshape <- FetchData(srt, vars = c(stat.by, group.by, split.by), cells = cells %||% rownames(meta.data), slot = slot)
    meta.reshape[["cells"]] <- rownames(meta.reshape)
    meta.reshape <- melt(meta.reshape, measure.vars = stat.by, variable.name = "Features", value.name = "Stat.by")
    rownames(meta.reshape) <- paste0(meta.reshape[["cells"]], "-", meta.reshape[["Features"]])
    exp.data <- matrix(0, nrow = 1, ncol = nrow(meta.reshape), dimnames = list("Stat.by", rownames(meta.reshape)))
    plist <- list()
    for (g in unique(meta.reshape[[group.by]])) {
      if (length(rownames(meta.reshape)[meta.reshape[[group.by]] == g]) > 0) {
        meta.use <- meta.reshape
        meta.use[[group.by]] <- NULL
        colnames(meta.use)[colnames(meta.use) == "Stat.by"] <- g
        p <- ExpressionStatPlot(
          exp.data = exp.data, meta.data = meta.use, stat.by = g, group.by = "Features", split.by = split.by, bg.by = NULL, plot.by = "group", fill.by = fill.by,
          cells = rownames(meta.reshape)[meta.reshape[[group.by]] == g], keep_empty = keep_empty, individual = individual,
          plot_type = plot_type,
          palette = palette, palcolor = palcolor, alpha = alpha,
          bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_alpha = bg_alpha,
          add_box = add_box, box_color = box_color, box_width = box_width, box_ptsize = box_ptsize,
          add_point = add_point, pt.color = pt.color, pt.size = pt.size, pt.alpha = pt.alpha, jitter.width = jitter.width, jitter.height = jitter.height,
          add_trend = add_trend, trend_color = trend_color, trend_linewidth = trend_linewidth, trend_ptsize = trend_ptsize,
          add_stat = add_stat, stat_color = stat_color, stat_size = stat_size, stat_stroke = stat_stroke, stat_shape = stat_shape,
          add_line = add_line, line_color = line_color, line_size = line_size, line_type = line_type,
          cells.highlight = cells.highlight, cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, alpha.highlight = alpha.highlight,
          calculate_coexp = calculate_coexp,
          same.y.lims = same.y.lims, y.min = y.min, y.max = y.max, y.trans = y.trans, y.nbreaks = y.nbreaks,
          sort = sort, stack = stack, flip = flip,
          comparisons = comparisons, ref_group = ref_group, pairwise_method = pairwise_method,
          multiplegroup_comparisons = multiplegroup_comparisons, multiple_method = multiple_method,
          sig_label = sig_label, sig_labelsize = sig_labelsize,
          aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
          legend.position = legend.position, legend.direction = legend.direction,
          theme_use = theme_use, theme_args = theme_args,
          force = force, seed = seed
        )
        plist <- append(plist, p)
      }
    }
    group.by <- "Features"
  } else {
    plist <- ExpressionStatPlot(
      exp.data = exp.data, meta.data = meta.data, stat.by = stat.by, group.by = group.by, split.by = split.by, bg.by = bg.by, plot.by = "group", fill.by = fill.by,
      cells = cells, keep_empty = keep_empty, individual = individual,
      plot_type = plot_type,
      palette = palette, palcolor = palcolor, alpha = alpha,
      bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_alpha = bg_alpha,
      add_box = add_box, box_color = box_color, box_width = box_width, box_ptsize = box_ptsize,
      add_point = add_point, pt.color = pt.color, pt.size = pt.size, pt.alpha = pt.alpha, jitter.width = jitter.width, jitter.height = jitter.height,
      add_trend = add_trend, trend_color = trend_color, trend_linewidth = trend_linewidth, trend_ptsize = trend_ptsize,
      add_stat = add_stat, stat_color = stat_color, stat_size = stat_size, stat_stroke = stat_stroke, stat_shape = stat_shape,
      add_line = add_line, line_color = line_color, line_size = line_size, line_type = line_type,
      cells.highlight = cells.highlight, cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, alpha.highlight = alpha.highlight,
      calculate_coexp = calculate_coexp,
      same.y.lims = same.y.lims, y.min = y.min, y.max = y.max, y.trans = y.trans, y.nbreaks = y.nbreaks,
      sort = sort, stack = stack, flip = flip,
      comparisons = comparisons, ref_group = ref_group, pairwise_method = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons, multiple_method = multiple_method,
      sig_label = sig_label, sig_labelsize = sig_labelsize,
      aspect.ratio = aspect.ratio, title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = legend.position, legend.direction = legend.direction,
      theme_use = theme_use, theme_args = theme_args,
      force = force, seed = seed
    )
  }

  plist_stack <- list()
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
              axis.text.x = element_text(vjust = c(1, 0)),
              axis.ticks.length.y = unit(0, "pt"),
              plot.margin = unit(c(0, -0.5, 0, 0), "mm")
            ))
          } else {
            suppressWarnings(p <- p + theme(
              legend.position = "none",
              panel.grid = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(vjust = c(1, 0)),
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

#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom ggplot2 geom_blank geom_violin geom_rect geom_boxplot geom_count geom_col geom_vline geom_hline layer_data layer_scales position_jitterdodge position_dodge stat_summary scale_x_discrete element_line element_text element_blank annotate mean_sdl after_stat scale_shape_identity
#' @importFrom Matrix rowSums
ExpressionStatPlot <- function(exp.data, meta.data, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, plot.by = c("group", "feature"), fill.by = c("group", "feature", "expression"),
                               cells = NULL, keep_empty = FALSE, individual = FALSE,
                               plot_type = c("violin", "box", "bar", "dot", "col"),
                               palette = "Paired", palcolor = NULL, alpha = 1,
                               bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
                               add_box = FALSE, box_color = "black", box_width = 0.1, box_ptsize = 2,
                               add_point = FALSE, pt.color = "grey30", pt.size = NULL, pt.alpha = 1, jitter.width = 0.4, jitter.height = 0.1,
                               add_trend = FALSE, trend_color = "black", trend_linewidth = 1, trend_ptsize = 2,
                               add_stat = c("none", "mean", "median"), stat_color = "black", stat_size = 1, stat_stroke = 1, stat_shape = 25,
                               add_line = NULL, line_color = "red", line_size = 1, line_type = 1,
                               cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1, alpha.highlight = 1,
                               calculate_coexp = FALSE,
                               same.y.lims = FALSE, y.min = NULL, y.max = NULL, y.trans = "identity", y.nbreaks = 5,
                               sort = FALSE, stack = FALSE, flip = FALSE,
                               comparisons = NULL, ref_group = NULL, pairwise_method = "wilcox.test",
                               multiplegroup_comparisons = FALSE, multiple_method = "kruskal.test",
                               sig_label = c("p.signif", "p.format"), sig_labelsize = 3.5,
                               aspect.ratio = NULL, title = NULL, subtitle = NULL, xlab = NULL, ylab = "Expression level",
                               legend.position = "right", legend.direction = "vertical",
                               theme_use = "theme_scp", theme_args = list(),
                               force = FALSE, seed = 11) {
  set.seed(seed)

  plot.by <- match.arg(plot.by)
  plot_type <- match.arg(plot_type)
  fill.by <- match.arg(fill.by)
  sig_label <- match.arg(sig_label)
  add_stat <- match.arg(add_stat)
  if (!is.null(add_line)) {
    stopifnot(is.numeric(add_line))
  }

  if (missing(exp.data)) {
    exp.data <- matrix(0, nrow = 1, ncol = nrow(meta.data), dimnames = list("", rownames(meta.data)))
  }

  allfeatures <- rownames(exp.data)
  allcells <- rownames(meta.data)

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
    group.by <- "All.groups"
    xlab <- ""
    meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
    meta.data[[split.by]] <- factor("")
  }
  if (group.by == split.by && group.by == "All.groups") {
    legend.position <- "none"
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
  if (!is.null(cells.highlight) && !isTRUE(cells.highlight)) {
    if (!any(cells.highlight %in% allcells)) {
      stop("No cells in 'cells.highlight' found.")
    }
    if (!all(cells.highlight %in% allcells)) {
      warning("Some cells in 'cells.highlight' not found.", immediate. = TRUE)
    }
    cells.highlight <- intersect(cells.highlight, allcells)
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- allcells
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
    ncomp <- sapply(comparisons, length)
    if (any(ncomp > 2)) {
      stop("'comparisons' must be a list in which all elements must be vectors of length 2")
    }
  }

  stat.by <- unique(stat.by)
  features_drop <- stat.by[!stat.by %in% c(rownames(exp.data), colnames(meta.data))]
  if (length(features_drop) > 0) {
    warning(paste0(features_drop, collapse = ","), " are not found.", immediate. = TRUE)
    stat.by <- stat.by[!stat.by %in% features_drop]
  }

  features_gene <- stat.by[stat.by %in% rownames(exp.data)]
  features_meta <- stat.by[stat.by %in% colnames(meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names: ", paste0(intersect(features_gene, features_meta), collapse = ","))
  }

  if (isTRUE(calculate_coexp) && length(features_gene) > 1) {
    if (length(features_meta) > 0) {
      warning(paste(features_meta, collapse = ","), "is not used when calculating co-expression", immediate. = TRUE)
    }
    status <- check_DataType(data = exp.data)
    message("Data type: ", status)
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
    if (all(allfeatures %in% features_gene)) {
      dat_gene <- t(exp.data)
    } else {
      dat_gene <- t(exp.data[features_gene, , drop = FALSE])
    }
  } else {
    dat_gene <- matrix(nrow = length(allcells), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = length(allcells), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  stat.by <- unique(stat.by[stat.by %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    stop("'stat.by' must be type of numeric variable.")
  }
  dat_group <- meta.data[, unique(c("cells", group.by, bg.by, split.by)), drop = FALSE]
  dat_use <- cbind(dat_group, dat_exp[row.names(dat_group), , drop = FALSE])
  if (!is.null(cells)) {
    dat_group <- dat_group[intersect(rownames(dat_group), cells), , drop = FALSE]
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (nrow(dat_group) == 0) {
    stop("No specified cells found.")
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

  if (isTRUE(same.y.lims)) {
    valus <- as_matrix(dat_use[, stat.by, drop = FALSE])[is.finite(as_matrix(dat_use[, stat.by, drop = FALSE]))]
    if (is.null(y.max)) {
      y.max <- max(valus, na.rm = TRUE)
    } else if (is.character(y.max)) {
      q.max <- as.numeric(sub("(^q)(\\d+)", "\\2", y.max)) / 100
      y.max <- quantile(values, q.max, na.rm = TRUE)
    }
    if (is.null(y.min)) {
      y.min <- min(valus, na.rm = TRUE)
    } else if (is.character(y.min)) {
      q.min <- as.numeric(sub("(^q)(\\d+)", "\\2", y.min)) / 100
      y.min <- quantile(values, q.min, na.rm = TRUE)
    }
  }

  plist <- list()

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
      if (split.by != "All.groups") {
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
    dat[["value"]] <- dat[[f]]
    dat[["group.by"]] <- dat[[g]]
    dat[["split.by"]] <- dat[[split.by]]
    if (split.by == g) {
      dat[["split.by"]] <- dat[["group.by"]]
    }
    # stat <- table(dat[, "group.by"], dat[, "split.by"])
    # stat_drop <- which(stat == 1, arr.ind = TRUE)
    # if (nrow(stat_drop) > 0) {
    #   for (j in 1:nrow(stat_drop)) {
    #     dat <- dat[!(dat[, "group.by"] == rownames(stat)[stat_drop[j, 1]] & dat[, "split.by"] == colnames(stat)[stat_drop[j, 2]]), , drop = FALSE]
    #     rownames(stat)[stat_drop[j, 1]]
    #   }
    # }

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
      dat[, "fill.by"] <- if (split.by == "All.groups") dat[, "group.by"] else dat[, "split.by"]
      keynm <- ifelse(split.by == "All.groups", g, split.by)
    }
    if (fill.by == "expression") {
      dat[, "fill.by"] <- median_values[paste0(dat[["group.by"]], "-", dat[["split.by"]]), f]
      keynm <- "Median expression"
    }
    if (split.by != "All.groups") {
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
    dat <- dat[order(dat[["group.unique"]]), , drop = FALSE]

    values <- dat[, "value"][is.finite(x = dat[, "value"])]
    if (is.null(y.max)) {
      y_max_use <- max(values, na.rm = TRUE)
    } else if (is.character(y.max)) {
      q.max <- as.numeric(sub("(^q)(\\d+)", "\\2", y.max)) / 100
      y_max_use <- quantile(values, q.max, na.rm = TRUE)
    } else {
      y_max_use <- y.max
    }
    if (is.null(y.min)) {
      y_min_use <- min(values, na.rm = TRUE)
    } else if (is.character(y.min)) {
      q.min <- as.numeric(sub("(^q)(\\d+)", "\\2", y.min)) / 100
      y_min_use <- quantile(values, q.min, na.rm = TRUE)
    } else {
      y_min_use <- y.min
    }

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
      bg_layer <- geom_rect(data = bg_data, xmin = bg_data[["xmin"]], xmax = bg_data[["xmax"]], ymin = bg_data[["ymin"]], ymax = bg_data[["ymax"]], fill = bg_data[["fill"]], alpha = bg_alpha, inherit.aes = FALSE)
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
        border_layer <- geom_vline(xintercept = border_data[["xintercept"]], linetype = 2, alpha = 0.5)
        p <- p + border_layer
      }
    }

    if (length(comparisons) > 0) {
      if (isTRUE(comparisons)) {
        group_use <- names(which(rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 2))
        if (any(rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 3)) {
          message("Detected more than 2 groups. Use multiple_method for comparison")
          method <- multiple_method
        } else {
          method <- pairwise_method
        }
        p <- p + ggpubr::stat_compare_means(
          data = dat[dat[["group.by"]] %in% group_use, , drop = FALSE],
          mapping = aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]]),
          label = sig_label,
          label.y = y_max_use,
          size = sig_labelsize,
          step.increase = 0.1,
          tip.length = 0.03,
          vjust = 1,
          method = method
        )

        y_max_use <- layer_scales(p)$y$range$range[2]
      } else {
        p <- p + ggpubr::stat_compare_means(
          mapping = aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]]),
          label = sig_label,
          label.y = y_max_use,
          size = sig_labelsize,
          step.increase = 0.1,
          tip.length = 0.03,
          vjust = 0,
          comparisons = comparisons,
          ref.group = ref_group,
          method = pairwise_method
        )
        y_max_use <- layer_scales(p)$y$range$range[1] + (layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]) * 1.15
      }
    }
    if (isTRUE(multiplegroup_comparisons)) {
      p <- p + ggpubr::stat_compare_means(
        aes(x = .data[["group.by"]], y = .data[["value"]], group = .data[["group.unique"]]),
        method = multiple_method,
        label = sig_label,
        label.y = y_max_use,
        size = sig_labelsize,
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
        position = position_jitterdodge(jitter.width = jitter.width, jitter.height = jitter.height, dodge.width = 0.9, seed = 11), show.legend = FALSE
      ))
      if (!is.null(cells.highlight)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight)
        if (nrow(cell_df) > 0) {
          p <- p + geom_point(
            data = cell_df, aes(x = .data[["group.by"]], y = .data[["value"]], linetype = rep(f, nrow(cell_df)), group = .data[["group.unique"]]), inherit.aes = FALSE,
            color = cols.highlight, size = sizes.highlight, alpha = alpha.highlight,
            position = position_jitterdodge(jitter.width = jitter.width, jitter.height = jitter.height, dodge.width = 0.9, seed = 11), show.legend = FALSE
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
        fun = add_stat, geom = "point", mapping = aes(group = .data[["split.by"]], shape = stat_shape),
        position = position_dodge(width = 0.9), color = stat_color, fill = stat_color, size = stat_size, stroke = stat_stroke,
      ) + scale_shape_identity()
    }
    if (!is.null(add_line)) {
      p <- p + geom_hline(
        yintercept = add_line,
        color = line_color, linetype = line_type, linewidth = line_size
      )
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

    if (isTRUE(flip)) {
      if (isTRUE(stack)) {
        p <- p + do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(angle = 90),
            panel.grid.major.x = element_line(color = "grey", linetype = 2),
            legend.position = legend.position,
            legend.direction = legend.direction
          ) + coord_flip(ylim = c(y_min_use, y_max_use))
      } else {
        p <- p + do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            strip.text.y = element_text(angle = 0),
            panel.grid.major.x = element_line(color = "grey", linetype = 2),
            legend.position = legend.position,
            legend.direction = legend.direction
          ) + coord_flip(ylim = c(y_min_use, y_max_use))
      }
    } else {
      p <- p + do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text.y = element_text(angle = 0),
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          legend.position = legend.position,
          legend.direction = legend.direction
        ) + coord_cartesian(ylim = c(y_min_use, y_max_use))
    }

    if (isTRUE(stack)) {
      p <- p + scale_y_continuous(
        trans = y.trans, breaks = c(y_min_use, y_max_use), labels = c(round(y_min_use, 1), round(y_max_use, 1))
      )
    } else {
      p <- p + scale_y_continuous(trans = y.trans, n.breaks = y.nbreaks)
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
        order = 1,
        override.aes = list(size = 4, color = "black", alpha = 1)
      ))
    } else {
      p <- p + scale_fill_gradientn(
        name = paste0(keynm, ":"), colours = colors, limits = colors_limits
      ) + guides(
        fill = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
      )
    }
    # plist[[paste0(f, ":", g, ":", paste0(single_group, collapse = ","), ":", paste0(sp, collapse = ","))]] <- p
  })

  return(plist)
}

#' Statistical plot of cells
#'
#' @inheritParams StatPlot
#' @param srt A Seurat object.
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#'
#' @seealso \code{\link{StatPlot}}
#'
#' @examples
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
                         bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
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
    bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_alpha = bg_alpha,
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
#' Visualizes data using various plot types such as bar plots, rose plots, ring plots, pie charts, trend plots, area plots, dot plots, sankey plots, chord plots, venn diagrams, and upset plots.
#'
#' @param meta.data The data frame containing the data to be plotted.
#' @param stat.by The column name(s) in \code{meta.data} specifying the variable(s) to be plotted.
#' @param group.by The column name in \code{meta.data} specifying the grouping variable.
#' @param split.by The column name in \code{meta.data} specifying the splitting variable.
#' @param bg.by The column name in \code{meta.data} specifying the background variable for bar plots.
#' @param flip Logical indicating whether to flip the plot.
#' @param NA_color The color to use for missing values.
#' @param NA_stat Logical indicating whether to include missing values in the plot.
#' @param keep_empty Logical indicating whether to keep empty groups in the plot.
#' @param individual Logical indicating whether to plot individual groups separately.
#' @param stat_level The level(s) of the variable(s) specified in \code{stat.by} to include in the plot.
#' @param plot_type The type of plot to create. Can be one of "bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", or "upset".
#' @param stat_type The type of statistic to compute for the plot. Can be one of "percent" or "count".
#' @param position The position adjustment for the plot. Can be one of "stack" or "dodge".
#' @param palette The name of the color palette to use for the plot.
#' @param palcolor The color to use in the color palette.
#' @param alpha The transparency level for the plot.
#' @param bg_palette The name of the background color palette to use for bar plots.
#' @param bg_palcolor The color to use in the background color palette.
#' @param bg_alpha The transparency level for the background color in bar plots.
#' @param label Logical indicating whether to add labels on the plot.
#' @param label.size The size of the labels.
#' @param label.fg The foreground color of the labels.
#' @param label.bg The background color of the labels.
#' @param label.bg.r The radius of the rounded corners of the label background.
#' @param aspect.ratio The aspect ratio of the plot.
#' @param title The main title of the plot.
#' @param subtitle The subtitle of the plot.
#' @param xlab The x-axis label of the plot.
#' @param ylab The y-axis label of the plot.
#' @param legend.position The position of the legend in the plot. Can be one of "right", "left", "bottom", "top", or "none".
#' @param legend.direction The direction of the legend in the plot. Can be one of "vertical" or "horizontal".
#' @param theme_use The name of the theme to use for the plot. Can be one of the predefined themes or a custom theme.
#' @param theme_args A list of arguments to be passed to the theme function.
#' @param combine Logical indicating whether to combine multiple plots into a single plot.
#' @param nrow The number of rows in the combined plot.
#' @param ncol The number of columns in the combined plot.
#' @param byrow Logical indicating whether to fill the plot by row or by column.
#' @param force Logical indicating whether to force the plot even if some variables have more than 100 levels.
#' @param seed The random seed to use for reproducible results.
#'
#' @seealso \code{\link{CellStatPlot}}
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
#' @importFrom grDevices png dev.control recordPlot dev.off
#' @export
StatPlot <- function(meta.data, stat.by, group.by = NULL, split.by = NULL, bg.by = NULL, flip = FALSE,
                     NA_color = "grey", NA_stat = TRUE, keep_empty = FALSE, individual = FALSE, stat_level = NULL,
                     plot_type = c("bar", "rose", "ring", "pie", "trend", "area", "dot", "sankey", "chord", "venn", "upset"),
                     stat_type = c("percent", "count"), position = c("stack", "dodge"),
                     palette = "Paired", palcolor = NULL, alpha = 1,
                     bg_palette = "Paired", bg_palcolor = NULL, bg_alpha = 0.2,
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
    group.by <- "All.groups"
    xlab <- ""
    meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
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

  if (any(group.by != "All.groups") && plot_type %in% c("sankey", "chord", "venn", "upset")) {
    warning("group.by is not used when plot sankey, chord, venn or upset", immediate. = TRUE)
  }
  if (stat_type == "percent" && plot_type %in% c("sankey", "chord", "venn", "upset")) {
    warning("stat_type is forcibly set to 'count' when plot sankey, chord, venn or upset", immediate. = TRUE)
    stat_type <- "count"
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
      colors_use <- colors[names(colors) %in% dat_split[[ifelse(split.by == "All.groups", 1, sp)]][[stat.by]]]
      if (any(is.na(dat_split[[ifelse(split.by == "All.groups", 1, sp)]][[stat.by]])) && isTRUE(NA_stat)) {
        colors_use <- c(colors_use, colors["NA"])
      }
      if (stat_type == "percent") {
        dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]] %>%
          xtabs(formula = paste0("~", stat.by, "+", g), addNA = NA_stat) %>%
          as.data.frame() %>%
          group_by(across(all_of(g)), .drop = FALSE) %>%
          mutate(groupn = sum(Freq)) %>%
          group_by(across(all_of(c(stat.by, g))), .drop = FALSE) %>%
          mutate(value = Freq / groupn) %>%
          as.data.frame()
      } else {
        dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]] %>%
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
        bg_layer <- geom_rect(data = bg_data, xmin = bg_data[["xmin"]], xmax = bg_data[["xmax"]], ymin = bg_data[["ymin"]], ymax = bg_data[["ymax"]], fill = bg_data[["fill"]], alpha = bg_alpha, inherit.aes = FALSE)
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
            alpha = alpha / 2, color = "grey50", position = position_use
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
            point.size = NA, max.overlaps = 100, min.segment.length = 0, force = 0,
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
            point.size = NA, max.overlaps = 100, min.segment.length = 0, force = 0,
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
      title <- title %||% sp
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
          order = 1,
          override.aes = list(size = 4, color = "black", alpha = 1)
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
      png(temp)
      dev.control("enable")
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
      dat_use <- dat_split[[ifelse(split.by == "All.groups", 1, sp)]]
      if (plot_type == "venn") {
        check_R(c("ggVennDiagram", "sf"))
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
            point.size = NA, max.overlaps = 100, force = 0,
            min.segment.length = 0, segment.colour = "black"
          ) +
          geom_text_repel(
            data = dat_stat, aes(label = label, x = x, y = y),
            colour = label.fg, size = label.size,
            bg.color = label.bg, bg.r = label.bg.r,
            point.size = NA, max.overlaps = 100, force = 0,
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
            point.size = NA, max.overlaps = 100, force = 0,
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
          df <- data.frame(factor(levels(dat_use[[l]]), levels = levels(dat_use[[l]])))
          colnames(df) <- l
          legend_list[[l]] <- get_legend(ggplot(data = df) +
            geom_col(aes(x = 1, y = 1, fill = .data[[l]]), color = "black") +
            scale_fill_manual(values = colors[levels(dat_use[[l]])]) +
            guides(fill = guide_legend(
              title.hjust = 0,
              title.vjust = 0,
              order = 1,
              override.aes = list(size = 4, color = "black", alpha = 1)
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
        p <- recordPlot()

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
    plot <- recordPlot()
    dev.off()
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
#' This function creates a correlation plot to visualize the pairwise correlations between selected features in a Seurat object.
#'
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to compare. Should be present in both the assay data and the metadata of the Seurat object.
#' @param group.by A character string specifying the column in the metadata to group cells by.
#' @param split.by A character string specifying the column in the metadata to split the plot by.
#' @param cells A character vector specifying the cells to include in the plot. If NULL (default), all cells will be included.
#' @param slot A character string specifying the slot in the Seurat object to use. Defaults to "data".
#' @param assay A character string specifying the assay to use. Defaults to the default assay in the Seurat object.
#' @param cor_method A character string specifying the correlation method to use. Can be "pearson" (default) or "spearman".
#' @param adjust A numeric value specifying the adjustment factor for the width of the violin plots. Defaults to 1.
#' @param margin A numeric value specifying the margin size for the plot. Defaults to 1.
#' @param reverse A logical value indicating whether to reverse the order of the features in the plot. Defaults to FALSE.
#' @param add_equation A logical value indicating whether to add the equation of the linear regression line to each scatter plot. Defaults to FALSE.
#' @param add_r2 A logical value indicating whether to add the R-squared value of the linear regression line to each scatter plot. Defaults to TRUE.
#' @param add_pvalue A logical value indicating whether to add the p-value of the linear regression line to each scatter plot. Defaults to TRUE.
#' @param add_smooth A logical value indicating whether to add a smoothed line to each scatter plot. Defaults to TRUE.
#' @param palette A character string specifying the name of the color palette to use for the groups. Defaults to "Paired".
#' @param palcolor A character string specifying the color for the groups. Defaults to NULL.
#' @param cor_palette A character string specifying the name of the color palette to use for the correlation. Defaults to "RuBu".
#' @param cor_palcolor A character string specifying the color for the correlation. Defaults to "RuBu".
#' @param cor_range A two-length numeric vector specifying the range for the correlation.
#' @param pt.size A numeric value specifying the size of the points in the scatter plots. If NULL (default), the size will be automatically determined based on the number of cells.
#' @param pt.alpha A numeric value between 0 and 1 specifying the transparency of the points in the scatter plots. Defaults to 1.
#' @param cells.highlight A logical value or a character vector specifying the cells to highlight in the scatter plots. If TRUE, all cells will be highlighted. Defaults to NULL.
#' @param cols.highlight A character string specifying the color for the highlighted cells. Defaults to "black".
#' @param sizes.highlight A numeric value specifying the size of the highlighted cells in the scatter plots. Defaults to 1.
#' @param alpha.highlight A numeric value between 0 and 1 specifying the transparency of the highlighted cells in the scatter plots. Defaults to 1.
#' @param stroke.highlight A numeric value specifying the stroke size of the highlighted cells in the scatter plots. Defaults to 0.5.
#' @param calculate_coexp A logical value indicating whether to calculate the co-expression of selected features. Defaults to FALSE.
#' @param raster A logical value indicating whether to use raster graphics for scatter plots. Defaults to NULL.
#' @param raster.dpi A numeric vector specifying the dpi (dots per inch) resolution for raster graphics in the scatter plots. Defaults to c(512, 512).
#' @param aspect.ratio A numeric value specifying the aspect ratio of the scatter plots. Defaults to 1.
#' @param title A character string specifying the title for the correlation plot. Defaults to NULL.
#' @param subtitle A character string specifying the subtitle for the correlation plot. Defaults to NULL.
#' @param legend.position A character string specifying the position of the legend. Can be "right" (default), "left", "top", or "bottom".
#' @param legend.direction A character string specifying the direction of the legend. Can be "vertical" (default) or "horizontal".
#' @param theme_use A character string specifying the name of the theme to use for the plot. Defaults to "theme_scp".
#' @param theme_args A list of arguments to pass to the theme function. Defaults to an empty list.
#' @param combine A logical value indicating whether to combine the plots into a single plot. Defaults to TRUE.
#' @param nrow A numeric value specifying the number of rows in the combined plot. If NULL (default), the number of rows will be automatically determined.
#' @param ncol A numeric value specifying the number of columns in the combined plot. If NULL (default), the number of columns will be automatically determined.
#' @param byrow A logical value indicating whether to fill the combined plot byrow (top to bottom, left to right). Defaults to TRUE.
#' @param force A logical value indicating whether to force the creation of the plot, even if it contains more than 50 subplots. Defaults to FALSE.
#' @param seed A numeric value specifying the random seed for reproducibility. Defaults to 11.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#' FeatureCorPlot(pancreas_sub, features = c("Neurog3", "Hes6", "Fev", "Neurod1", "Rbp4", "Pyy"), group.by = "SubCellType")
#' FeatureCorPlot(pancreas_sub,
#'   features = c("nFeature_RNA", "nCount_RNA", "nFeature_spliced", "nCount_spliced", "nFeature_unspliced", "nCount_unspliced"),
#'   group.by = "SubCellType", cor_palette = "Greys", cor_range = c(0, 1)
#' )
#' FeatureCorPlot(pancreas_sub,
#'   features = c("nFeature_RNA", "nCount_RNA"),
#'   group.by = "SubCellType", add_equation = TRUE
#' )
#' @importFrom Seurat Reductions Embeddings Key
#' @importFrom SeuratObject as.sparse
#' @importFrom dplyr group_by "%>%" .data
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_density_2d stat_density_2d labs scale_x_continuous scale_y_continuous facet_grid scale_color_gradientn scale_fill_gradientn scale_colour_gradient scale_fill_gradient guide_colorbar scale_color_identity scale_fill_identity guide_colorbar geom_hex stat_summary_hex
#' @importFrom ggnewscale new_scale_color new_scale_fill
#' @importFrom ggrepel geom_text_repel GeomTextRepel
#' @importFrom gtable gtable_add_cols
#' @importFrom patchwork wrap_plots
#' @importFrom Matrix t
#' @importFrom methods slot
#' @export
FeatureCorPlot <- function(srt, features, group.by = NULL, split.by = NULL, cells = NULL, slot = "data", assay = NULL,
                           cor_method = "pearson", adjust = 1, margin = 1, reverse = FALSE,
                           add_equation = FALSE, add_r2 = TRUE, add_pvalue = TRUE, add_smooth = TRUE,
                           palette = "Paired", palcolor = NULL, cor_palette = "RdBu", cor_palcolor = NULL, cor_range = c(-1, 1),
                           pt.size = NULL, pt.alpha = 1,
                           cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                           calculate_coexp = FALSE,
                           raster = NULL, raster.dpi = c(512, 512),
                           aspect.ratio = 1, title = NULL, subtitle = NULL,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)

  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    srt@meta.data[[group.by]] <- factor("")
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
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names: ", paste0(intersect(features_gene, features_meta), collapse = ","))
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
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
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
    dat_exp <- as.sparse(as_matrix(dat_exp))
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
  cor_colors <- palette_scp(x = seq(cor_range[1], cor_range[2], length.out = 200), palette = cor_palette, palcolor = cor_palcolor)
  bound <- strsplit(gsub("\\(|\\)|\\[|\\]", "", names(cor_colors)), ",")
  bound <- lapply(bound, as.numeric)
  df_bound <- do.call(rbind, bound)
  rownames(df_bound) <- cor_colors
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
          position = ifelse(isTRUE(reverse), "top", "bottom")
        ) +
          scale_y_continuous(
            n.breaks = 3, labels = scales::number_format(),
            position = ifelse(isTRUE(reverse), "right", "left")
          )
      }
      if (f1_index < f2_index) {
        if (isTRUE(raster)) {
          p <- p + scattermore::geom_scattermore(
            mapping = aes(x = .data[[f1]], y = .data[[f2]], color = .data[[group.by]]),
            pointsize = ceiling(pt.size), alpha = pt.alpha, pixels = raster.dpi
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
            geom = GeomTextRepel, x = -Inf, y = Inf, label = eqs[i],
            color = "black", bg.color = "white", bg.r = 0.1, size = 3.5, point.size = NA,
            max.overlaps = 100, force = 0, min.segment.length = Inf,
            hjust = -0.05, vjust = vjusts[1:sum(i)], parse = TRUE
          )
        }
        if (!is.null(cells.highlight)) {
          cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight)
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
        label_pos <- (max(dat_exp[rownames(dat), ], na.rm = TRUE) + min(dat_exp[rownames(dat), ], na.rm = TRUE)) / 2
        fill <- rownames(df_bound)[df_bound[, 1] < pair_sim[f1, f2] & df_bound[, 2] >= pair_sim[f1, f2]]
        p <- p + annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = fill) +
          annotate(
            geom = GeomTextRepel, x = label_pos, y = label_pos, label = label,
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
        labels = names(colors)
      ) + scale_fill_manual(
        name = paste0(group.by, ":"),
        values = colors,
        labels = names(colors)
      )
      return(p)
    }, x = order1, y = order2, SIMPLIFY = FALSE)

    legend_list <- NULL
    if (length(features) > 1) {
      legend_list[["correlation"]] <- get_legend(ggplot(data.frame(range = cor_range, x = 1, y = 1), aes(x = x, y = y, fill = range)) +
        geom_point() +
        scale_fill_gradientn(
          name = paste0("Correlation"),
          limits = cor_range,
          n.breaks = 3,
          colors = cor_colors,
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        ))
    }
    if (nlevels(dat[[group.by]]) > 1) {
      legend_list[["group.by"]] <- suppressWarnings(get_legend(plotlist[[1]] +
        guides(fill = guide_legend(
          title.hjust = 0,
          order = 1,
          override.aes = list(size = 4, color = "black", alpha = 1)
        )) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )))
    }

    grob_row <- list()
    plotlist <- suppressWarnings(lapply(plotlist, as_grob))
    for (i in seq(1, length(plotlist), length(features))) {
      grob_row[[paste0(i:(i + length(features) - 1), collapse = "-")]] <- do.call(cbind, plotlist[i:(i + length(features) - 1)])
    }
    gtable <- do.call(rbind, grob_row)
    if (length(legend_list) > 0) {
      legend_list <- legend_list[!sapply(legend_list, is.null)]
      if (legend.direction == "vertical") {
        legend <- do.call(cbind, legend_list)
      } else {
        legend <- do.call(rbind, legend_list)
      }
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
#' Plots the density of specified features in a single or multiple groups,
#' grouped by specified variables.
#'
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to plot.
#' @param group.by A character vector specifying the variables to group the data by.
#' @param split.by A character vector specifying the variables to split the data by.
#' Default is NULL, which means no splitting is performed.
#' @param assay A character specifying the assay to use from the Seurat object.
#'   Default is NULL, which means the default assay will be used.
#' @param slot A character specifying the slot to use from the assay. Default is "data".
#' @param flip A logical indicating whether to flip the x-axis. Default is FALSE.
#' @param reverse A logical indicating whether to reverse the y-axis. Default is FALSE.
#' @param x_order A character specifying how to order the x-axis. Can be "value" or "rank". Default is "value".
#' @param decreasing A logical indicating whether to order the groups in decreasing order. Default is NULL.
#' @param palette A character specifying the color palette to use for grouping variables. Default is "Paired".
#' @param palcolor A character specifying the color to use for each group. Default is NULL.
#' @param cells A character vector specifying the cells to plot. Default is NULL, which means all cells are included.
#' @param keep_empty A logical indicating whether to keep empty groups. Default is FALSE.
#' @param y.nbreaks An integer specifying the number of breaks on the y-axis. Default is 4.
#' @param y.min A numeric specifying the minimum value on the y-axis. Default is NULL, which means the minimum value will be automatically determined.
#' @param y.max A numeric specifying the maximum value on the y-axis. Default is NULL, which means the maximum value will be automatically determined.
#' @param same.y.lims A logical indicating whether to use the same y-axis limits for all plots. Default is FALSE.
#' @param aspect.ratio A numeric specifying the aspect ratio of the plot. Default is NULL, which means the aspect ratio will be automatically determined.
#' @param title A character specifying the title of the plot. Default is NULL.
#' @param subtitle A character specifying the subtitle of the plot. Default is NULL.
#' @param legend.position A character specifying the position of the legend. Default is "right".
#' @param legend.direction A character specifying the direction of the legend. Default is "vertical".
#' @param theme_use A character specifying the theme to use. Default is "theme_scp".
#' @param theme_args A list of arguments to pass to the theme function.
#' @param combine A logical indicating whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow An integer specifying the number of rows in the combined plot.
#'   Default is NULL, which means determined automatically based on the number of plots.
#' @param ncol An integer specifying the number of columns in the combined plot.
#'   Default is NULL, which means determined automatically based on the number of plots.
#' @param byrow A logical indicating whether to add plots by row or by column in the combined plot. Default is TRUE.
#' @param force A logical indicating whether to continue plotting if there are more than 50 features. Default is FALSE.
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
CellDensityPlot <- function(srt, features, group.by = NULL, split.by = NULL, assay = NULL, slot = "data",
                            flip = FALSE, reverse = FALSE, x_order = c("value", "rank"),
                            decreasing = NULL, palette = "Paired", palcolor = NULL,
                            cells = NULL, keep_empty = FALSE,
                            y.nbreaks = 4, y.min = NULL, y.max = NULL, same.y.lims = FALSE,
                            aspect.ratio = NULL, title = NULL, subtitle = NULL,
                            legend.position = "right", legend.direction = "vertical",
                            theme_use = "theme_scp", theme_args = list(),
                            combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE) {
  check_R("ggridges")
  assay <- assay %||% DefaultAssay(srt)
  x_order <- match.arg(x_order)
  if (is.null(features)) {
    stop("'features' must be provided.")
  }
  if (!inherits(features, "character")) {
    stop("'features' is not a character vectors")
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    srt@meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  if (group.by == split.by & group.by == "All.groups") {
    legend.position <- "none"
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
  if (length(intersect(features_gene, features_meta)) > 0) {
    warning("Features appear in both gene names and metadata names: ", paste0(intersect(features_gene, features_meta), collapse = ","))
  }

  if (length(features_gene) > 0) {
    dat_gene <- t(slot(srt@assays[[assay]], slot)[features_gene, , drop = FALSE])
  } else {
    dat_gene <- matrix(nrow = ncol(srt@assays[[1]]), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(srt@meta.data[, features_meta, drop = FALSE])
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
    y.max <- max(as_matrix(dat_exp[, features])[is.finite(as_matrix(dat_exp[, features]))], na.rm = TRUE)
  }
  if (isTRUE(same.y.lims) && is.null(y.min)) {
    y.min <- min(as_matrix(dat_exp[, features])[is.finite(as_matrix(dat_exp[, features]))], na.rm = TRUE)
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
        # stat <- table(dat[[g]])
        # stat_drop <- names(which(stat <= 2))
        # if (length(stat_drop) > 0) {
        #   dat <- dat[!dat[[g]] %in% stat_drop, , drop = FALSE]
        # }
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
        if (split.by != "All.groups") {
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
#' Generate a lineage plot based on the pseudotime.
#'
#' @param srt An object of class Seurat.
#' @param lineages A character vector that specifies the lineages to be included. Typically, use the pseudotime of cells.
#' @param reduction An optional string specifying the dimensionality reduction method to use.
#' @param dims A numeric vector of length 2 specifying the dimensions to plot.
#' @param cells An optional character vector specifying the cells to include in the plot.
#' @param trim A numeric vector of length 2 specifying the quantile range of lineages to include in the plot.
#' @param span A numeric value specifying the span of the loess smoother.
#' @param palette A character string specifying the color palette to use for the lineages.
#' @param palcolor An optional string specifying the color for the palette.
#' @param lineages_arrow An arrow object specifying the arrow for lineages.
#' @param linewidth A numeric value specifying the linewidth for the lineages.
#' @param line_bg A character string specifying the color for the background lines.
#' @param line_bg_stroke A numeric value specifying the stroke width for the background lines.
#' @param whiskers A logical value indicating whether to include whiskers in the plot.
#' @param whiskers_linewidth A numeric value specifying the linewidth for the whiskers.
#' @param whiskers_alpha A numeric value specifying the transparency for the whiskers.
#' @param aspect.ratio A numeric value specifying the aspect ratio of the plot.
#' @param title An optional character string specifying the plot title.
#' @param subtitle An optional character string specifying the plot subtitle.
#' @param xlab An optional character string specifying the x-axis label.
#' @param ylab An optional character string specifying the y-axis label.
#' @param legend.position A character string specifying the position of the legend.
#' @param legend.direction A character string specifying the direction of the legend.
#' @param theme_use A character string specifying the theme to use for the plot.
#' @param theme_args A list of additional arguments to pass to the theme function.
#' @param return_layer A logical value indicating whether to return the plot as a layer.
#' @param seed An optional integer specifying the random seed for reproducibility.
#'
#' @seealso \code{\link{RunSlingshot}} \code{\link{CellDimPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", show_plot = FALSE)
#' LineagePlot(pancreas_sub, lineages = paste0("Lineage", 1:3))
#' LineagePlot(pancreas_sub, lineages = paste0("Lineage", 1:3), whiskers = TRUE)
#' @importFrom Seurat Key Embeddings
#' @importFrom ggplot2 aes geom_path geom_segment labs
#' @importFrom grid arrow unit
#' @importFrom stats loess quantile
#' @export
