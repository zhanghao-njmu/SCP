DynamicHeatmap <- function(srt, lineages, features = NULL, use_fitted = FALSE, border = TRUE, flip = FALSE,
                           min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05, num_intersections = NULL,
                           cell_density = 1, cell_bins = 100, order_by = c("peaktime", "valleytime"),
                           slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), exp_legend_title = NULL, limits = NULL,
                           lib_normalize = identical(slot, "counts"), libsize = NULL, family = NULL,
                           cluster_features_by = NULL, cluster_rows = FALSE, cluster_row_slices = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                           row_title = NULL, column_title = NULL, row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                           feature_split = NULL, feature_split_by = NULL, n_split = NULL, split_order = NULL,
                           split_method = c("mfuzz", "kmeans", "kmeans-peaktime", "hclust", "hclust-peaktime"), decreasing = FALSE,
                           fuzzification = NULL,
                           anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                           terms_width = unit(4, "in"), terms_fontsize = 8,
                           keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                           features_width = unit(2, "in"), features_fontsize = c(6, 10),
                           IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                           db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                           GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, words_excluded = NULL,
                           nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                           pseudotime_label = NULL, pseudotime_label_color = "black", pseudotime_label_linetype = 2, pseudotime_label_linewidth = 3,
                           heatmap_palette = "viridis", heatmap_palcolor = NULL,
                           pseudotime_palette = "cividis", pseudotime_palcolor = NULL,
                           feature_split_palette = "simspec", feature_split_palcolor = NULL,
                           cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = unit(5, "mm")) else list(height = unit(5, "mm")),
                           feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = if (flip) list(height = unit(5, "mm")) else list(width = unit(5, "mm")),
                           separate_annotation = NULL, separate_annotation_palette = "Paired", separate_annotation_palcolor = NULL, separate_annotation_params = if (flip) list(width = unit(10, "mm")) else list(height = unit(10, "mm")),
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
  exp_name <- exp_legend_title %||% exp_name

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
    bins <- cut(Pseudotime_assign, breaks = seq(min(Pseudotime_assign, na.rm = TRUE), max(Pseudotime_assign, na.rm = TRUE), length.out = cell_bins), include.lowest = TRUE)
    ncells <- ceiling(max(table(bins), na.rm = TRUE) * cell_density)
    message("ncell/bin=", ncells, "(", cell_bins, "bins)")
    cell_keep <- unlist(sapply(levels(bins), function(x) {
      cells <- names(Pseudotime_assign)[bins == x]
      out <- sample(cells, size = min(length(cells), ncells))
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
      features <- c(features, DynamicFeatures[["features"]])
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

  all_calculated <- Reduce(intersect, lapply(lineages, function(l) rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]])))
  features_tab <- table(features)
  features <- names(features_tab)[which(features_tab %in% (num_intersections %||% seq_along(lineages)))]
  if (!all(features %in% all_calculated)) {
    message("Some features were missing in at least one lineage: \n", paste0(head(features[!features %in% all_calculated], 10), collapse = ","), "...")
    features <- intersect(features, all_calculated)
  }
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
      mat_tmp <- as_matrix(rbind(slot(srt@assays[[assay]], slot)[gene, cells, drop = FALSE], t(srt@meta.data[cells, meta, drop = FALSE])))[features, , drop = FALSE]
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
            g <- as_grob(subplots_list[['", nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
            g <- as_grob(subplots_list[['", nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
          status <- tryCatch(check_R("e1071"), error = identity)
          if (inherits(status, "error")) {
            warning("The e1071 package was not found. Switch split_method to 'kmeans'", immediate. = TRUE)
            split_method <- "kmeans"
          } else {
            mat_split_tmp <- mat_split
            colnames(mat_split_tmp) <- make.unique(colnames(mat_split_tmp))
            mat_split_tmp <- standardise(mat_split_tmp)
            min_fuzzification <- mestimate(mat_split_tmp)
            if (is.null(fuzzification)) {
              fuzzification <- min_fuzzification + 0.1
            } else {
              if (fuzzification <= min_fuzzification) {
                warning("fuzzification value is samller than estimated:", round(min_fuzzification, 2), immediate. = TRUE)
              }
            }
            cl <- e1071::cmeans(mat_split_tmp, centers = n_split, method = "cmeans", m = fuzzification)
            if (length(cl$cluster) == 0) {
              stop("Clustering with mfuzz failed (fuzzification=", round(fuzzification, 2), "). Please set a larger fuzzification parameter manually.")
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
          hc <- hclust(as.dist(dist(mat_split)))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
        if (split_method == "hclust-peaktime") {
          feature_y <- feature_metadata[rownames(mat_split), order_by]
          names(feature_y) <- rownames(mat_split)
          hc <- hclust(as.dist(dist(feature_y)))
          row_split <- feature_split <- cutree(hc, k = n_split)
        }
      }
      df <- data.frame(row_split = row_split, order_by = feature_metadata[names(row_split), order_by])
      df_order <- aggregate(df, by = list(row_split), FUN = mean)
      df_order <- df_order[order(df_order[["order_by"]], decreasing = decreasing), , drop = FALSE]
      if (!is.null(split_order)) {
        df_order <- df_order[split_order, , drop = FALSE]
      }
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
      dend <- as.dendrogram(hclust(as.dist(dist(mat_cluster))))
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
        width = unit(5, "mm"),
        height = unit(5, "mm"),
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
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, words_excluded = words_excluded
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
      row_title = row_title %||% if (flip) l else character(0),
      row_title_side = row_title_side,
      column_title = column_title %||% if (flip) character(0) else l,
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
    if (is.null(width) || is.null(height)) {
      message("\nThe size of the heatmap is fixed because certain elements are not scalable.\nThe width and height of the heatmap are determined by the size of the current viewport.\nIf you want to have more control over the size, you can manually set the parameters 'width' and 'height'.\n")
    }
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
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
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

#' DynamicPlot
#'
#' Plot dynamic features across pseudotime.
#'
#' @param srt A Seurat object.
#' @param features A character vector specifying the features to plot.
#' @param lineages A character vector specifying the lineages to plot.
#' @param group.by A character specifying a metadata column to group the cells by. Default is NULL.
#' @param cells A character vector specifying the cells to include in the plot. Default is NULL.
#' @param slot A character string specifying the slot to use for the analysis. Default is "counts".
#' @param assay A character string specifying the assay to use for the analysis. Default is NULL.
#' @param family A character specifying the model used to calculate the dynamic features if needed. By default, this parameter is set to NULL, and the appropriate family will be automatically determined.
#' @param exp_method A character specifying the method to transform the expression values. Default is "log1p" with options "log1p", "raw", "zscore", "fc", "log2fc".
#' @param lib_normalize A boolean specifying whether to normalize the expression values using library size. By default, if the \code{slot} is counts, this parameter is set to TRUE. Otherwise, it is set to FALSE.
#' @param libsize A numeric vector specifying the library size for each cell. Default is NULL.
#' @param compare_lineages A boolean specifying whether to compare the lineages in the plot. Default is TRUE.
#' @param compare_features A boolean specifying whether to compare the features in the plot. Default is FALSE.
#' @param add_line A boolean specifying whether to add lines to the plot. Default is TRUE.
#' @param add_interval A boolean specifying whether to add confidence intervals to the plot. Default is TRUE.
#' @param line.size A numeric specifying the size of the lines. Default is 1.
#' @param line_palette A character string specifying the name of the palette to use for the line colors. Default is "Dark2".
#' @param line_palcolor A vector specifying the colors to use for the line palette. Default is NULL.
#' @param add_point A boolean specifying whether to add points to the plot. Default is TRUE.
#' @param pt.size A numeric specifying the size of the points. Default is 1.
#' @param point_palette A character string specifying the name of the palette to use for the point colors. Default is "Paired".
#' @param point_palcolor A vector specifying the colors to use for the point palette. Default is NULL.
#' @param add_rug A boolean specifying whether to add rugs to the plot. Default is TRUE.
#' @param flip A boolean specifying whether to flip the x-axis. Default is FALSE.
#' @param reverse A boolean specifying whether to reverse the x-axis. Default is FALSE.
#' @param x_order A character specifying the order of the x-axis values. Default is c("value", "rank").
#' @param aspect.ratio A numeric specifying the aspect ratio of the plot. Default is NULL.
#' @param legend.position A character string specifying the position of the legend in the plot. Default is "right".
#' @param legend.direction A character string specifying the direction of the legend in the plot. Default is "vertical".
#' @param theme_use A character string specifying the name of the theme to use for the plot. Default is "theme_scp".
#' @param theme_args A list specifying the arguments to pass to the theme function. Default is list().
#' @param combine A boolean specifying whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow A numeric specifying the number of rows in the combined plot. Default is NULL.
#' @param ncol A numeric specifying the number of columns in the combined plot. Default is NULL.
#' @param byrow A boolean specifying whether to fill plots by row in the combined plot. Default is TRUE.
#' @param seed A numeric specifying the random seed. Default is 11.
#'
#' @seealso \code{\link{RunDynamicFeatures}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
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
#' @importFrom reshape2 melt
#' @export
DynamicPlot <- function(srt, lineages, features, group.by = NULL, cells = NULL, slot = "counts", assay = NULL, family = NULL,
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
    raw_matrix_list[[l]] <- as_matrix(raw_matrix[, features, drop = FALSE])
    fitted_matrix_list[[l]] <- as_matrix(fitted_matrix[, features, drop = FALSE])
    upr_matrix_list[[l]] <- as_matrix(upr_matrix[, features, drop = FALSE])
    lwr_matrix_list[[l]] <- as_matrix(lwr_matrix[, features, drop = FALSE])
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
    raw <- melt(raw, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    fitted <- as.data.frame(cbind(cell_metadata[rownames(fitted_matrix), c(l, "x_assign")], fitted_matrix))
    colnames(fitted)[1] <- "Pseudotime"
    fitted[["Cell"]] <- rownames(fitted)
    fitted[["Value"]] <- "fitted"
    fitted <- melt(fitted, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    upr <- as.data.frame(cbind(cell_metadata[rownames(upr_matrix), c(l, "x_assign")], upr_matrix))
    colnames(upr)[1] <- "Pseudotime"
    upr[["Cell"]] <- rownames(upr)
    upr[["Value"]] <- "upr"
    upr <- melt(upr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

    lwr <- as.data.frame(cbind(cell_metadata[rownames(lwr_matrix), c(l, "x_assign")], lwr_matrix))
    colnames(lwr)[1] <- "Pseudotime"
    lwr[["Cell"]] <- rownames(lwr)
    lwr[["Value"]] <- "lwr"
    lwr <- melt(lwr, id.vars = c("Cell", "Pseudotime", "x_assign", "Value"), value.name = "exp", variable.name = "Features")

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

GroupTreePlot <- function() {

}

#' Projection Plot
#'
#' This function generates a projection plot, which can be used to compare two groups of cells in a dimensionality reduction space.
#'
#' @param srt_query An object of class Seurat storing the query group cells.
#' @param srt_ref An object of class Seurat storing the reference group cells.
#' @param query_group The grouping variable for the query group cells.
#' @param ref_group The grouping variable for the reference group cells.
#' @param query_reduction The name of the reduction in the query group cells.
#' @param ref_reduction The name of the reduction in the reference group cells.
#' @param query_param A list of parameters for customizing the query group plot. Available parameters: palette (color palette for groups) and cells.highlight (whether to highlight cells).
#' @param ref_param A list of parameters for customizing the reference group plot. Available parameters: palette (color palette for groups) and cells.highlight (whether to highlight cells).
#' @param xlim The x-axis limits for the plot. If not provided, the limits will be calculated based on the data.
#' @param ylim The y-axis limits for the plot. If not provided, the limits will be calculated based on the data.
#' @param pt.size The size of the points in the plot.
#' @param stroke.highlight The size of the stroke highlight for cells.
#'
#' @examples
#' data("panc8_sub")
#' srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
#' srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
#' srt_ref <- Integration_SCP(srt_ref, batch = "tech", integration_method = "Seurat")
#' CellDimPlot(srt_ref, group.by = c("celltype", "tech"))
#'
#' # Projection
#' srt_query <- RunKNNMap(srt_query = srt_query, srt_ref = srt_ref, ref_umap = "SeuratUMAP2D")
#' ProjectionPlot(srt_query = srt_query, srt_ref = srt_ref, query_group = "celltype", ref_group = "celltype")
#'
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous ggplot_build theme geom_point aes scale_fill_identity facet_null
#' @importFrom ggnewscale new_scale_fill new_scale_color
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom grid grob
#' @importFrom rlang  %||%
#' @importFrom patchwork wrap_plots
#' @export
ProjectionPlot <- function(srt_query, srt_ref,
                           query_group = NULL, ref_group = NULL,
                           query_reduction = "ref.embeddings", ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||% NULL,
                           query_param = list(palette = "Set1", cells.highlight = TRUE), ref_param = list(palette = "Paired"),
                           xlim = NULL, ylim = NULL, pt.size = 0.8, stroke.highlight = 0.5) {
  if (is.null(ref_reduction)) {
    stop("Please specify the ref_reduction.")
  }
  query_param[["show_stat"]] <- FALSE
  ref_param[["show_stat"]] <- FALSE

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
    guides(color = guide_legend(title = paste0("Ref: ", ref_group), override.aes = list(size = 4)))
  p1legend <- get_legend(p1)
  # p1legend <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")


  p2 <- do.call(CellDimPlot, args = c(
    srt = srt_query, reduction = query_reduction, group.by = query_group,
    query_param
  )) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim)
  p2data <- ggplot_build(p2)$data[[1]]
  color <- p2data$colour
  names(color) <- p2$data$group.by
  p2 <- p2 + guides(color = guide_legend(
    title = paste0("Query: ", query_group),
    override.aes = list(size = 4, shape = 21, color = "black", fill = na.omit(color[levels(p2$data$group.by)]))
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
#' This function generates various types of plots for enrichment (over-representation) analysis.
#'
#' @param srt A Seurat object containing the results of RunDEtest and RunEnrichment.
#' If specified, enrichment results will be extracted from the Seurat object automatically.
#' If not specified, the \code{res} arguments must be provided.
#' @param group_by A character vector specifying the grouping variable in the Seurat object. This argument is only used if \code{srt} is specified.
#' @param test.use A character vector specifying the test to be used in differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param res Enrichment results generated by RunEnrichment function. If provided, 'srt', 'test.use' and 'group_by' are ignored.
#' @param db The database to use for enrichment plot. Default is "GO_BP".
#' @param plot_type The type of plot to generate. Options are: "bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison". Default is "bar".
#' @param split_by The splitting variable(s) for the plot. Can be "Database", "Groups", or both. Default is c("Database", "Groups") for plots.
#' @param color_by The variable used for coloring. Default is "Database".
#' @param group_use The group(s) to be used for enrichment plot. Default is NULL.
#' @param id_use List of IDs to be used to display specific terms in the enrichment plot. Default value is NULL.
#' @param pvalueCutoff The p-value cutoff. Default is NULL. Only work when \code{padjustCutoff} is NULL.
#' @param padjustCutoff The p-adjusted cutoff. Default is 0.05.
#' @param topTerm The number of top terms to display. Default is 6, or 100 if 'plot_type' is "enrichmap".
#' @param compare_only_sig Whether to compare only significant terms. Default is FALSE.
#' @param topWord The number of top words to display for wordcloud. Default is 100.
#' @param word_type The type of words to display in wordcloud. Options are "term" and "feature". Default is "term".
#' @param word_size The size range for words in wordcloud. Default is c(2, 8).
#' @param words_excluded Words to be excluded from the wordcloud. The default value is NULL, which means that the built-in words (SCP::words_excluded) will be used.
#' @param network_layout The layout algorithm to use for network plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param network_labelsize The label size for network plot. Default is 5.
#' @param network_blendmode The blend mode for network plot. Default is "blend".
#' @param network_layoutadjust Whether to adjust the layout of the network plot to avoid overlapping words. Default is TRUE.
#' @param network_adjscale The scale for adjusting network plot layout. Default is 60.
#' @param network_adjiter The number of iterations for adjusting network plot layout. Default is 100.
#' @param enrichmap_layout The layout algorithm to use for enrichmap plot. Options are "fr", "kk","random", "circle", "tree", "grid", or other algorithm from 'igraph' package. Default is "fr".
#' @param enrichmap_cluster The clustering algorithm to use for enrichmap plot. Options are "walktrap", "fast_greedy", or other algorithm from 'igraph' package. Default is "fast_greedy".
#' @param enrichmap_label  The label type for enrichmap plot. Options are "term" and "feature". Default is "term".
#' @param enrichmap_labelsize The label size for enrichmap plot. Default is 5.
#' @param enrlichmap_nlabel The number of labels to display for each cluster in enrichmap plot. Default is 4.
#' @param enrichmap_show_keyword Whether to show the keyword of terms or features in enrichmap plot. Default is FALSE.
#' @param enrichmap_mark The mark shape for enrichmap plot. Options are "ellipse" and "hull". Default is "ellipse".
#' @param enrichmap_expand The expansion factor for enrichmap plot. Default is c(0.5, 0.5).
#' @param character_width  The maximum width of character of descriptions. Default is 50.
#' @param lineheight The line height for y-axis labels. Default is 0.5.
#' @param palette The color palette to use. Default is "Spectral".
#' @param palcolor Custom colors for palette. Default is NULL.
#' @param aspect.ratio The aspect ratio of the plot. Default is 1.
#' @param legend.position The position of the legend. Default is "right".
#' @param legend.direction The direction of the legend. Default is "vertical".
#' @param theme_use The theme to use for the plot. Default is "theme_scp".
#' @param theme_args The arguments to pass to the theme. Default is an empty list.
#' @param combine Whether to combine multiple plots into a single plot. Default is TRUE.
#' @param nrow The number of rows in the combined plot. Default is NULL, calculated based on the number of plots.
#' @param ncol The number of columns in the combined plot. Default is NULL, calculated based on the number of plots.
#' @param byrow  Whether to fill the combined plot by row. Default is TRUE.
#' @param seed The random seed to use. Default is 11.
#'
#' @seealso \code{\link{RunEnrichment}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, db = c("GO_BP", "GO_CC"), group_by = "CellType", species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "bar")
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = "Endocrine",
#'   plot_type = "bar", character_width = 30,
#'   theme_use = ggplot2::theme_classic, theme_args = list(base_size = 10)
#' )
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", plot_type = "bar", color_by = "Groups", ncol = 2)
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", plot_type = "bar",
#'   split_by = "Database", color_by = "Groups", palette = "Set1",
#'   id_use = list(
#'     "Ductal" = c("GO:0002181", "GO:0045787", "GO:0006260", "GO:0050679"),
#'     "Ngn3 low EP" = c("GO:0050678", "GO:0051101", "GO:0072091", "GO:0006631"),
#'     "Ngn3 high EP" = c("GO:0035270", "GO:0030325", "GO:0008637", "GO:0030856"),
#'     "Pre-endocrine" = c("GO:0090276", "GO:0031018", "GO:0030073", "GO:1903532"),
#'     "Endocrine" = c("GO:0009914", "GO:0030073", "GO:0009743", "GO:0042593")
#'   )
#' )
#'
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", compare_only_sig = TRUE)
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "comparison")
#'
#' EnrichmentPlot(pancreas_sub,
#'   db = c("GO_BP", "GO_CC"), group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "bar",
#'   split_by = "Groups"
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = c("GO_BP", "GO_CC"), group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "bar",
#'   split_by = "Database", color_by = "Groups",
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = c("GO_BP", "GO_CC"), group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "bar",
#'   split_by = c("Database", "Groups")
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = c("GO_BP", "GO_CC"), group_by = "CellType", group_use = c("Ductal", "Endocrine"), plot_type = "bar",
#'   split_by = c("Groups", "Database")
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = c("GO_BP", "GO_CC"), group_by = "CellType", plot_type = "bar",
#'   split_by = "Database", color_by = "Groups", palette = "Set1"
#' )
#'
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "dot", palette = "GdRd")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "lollipop", palette = "GdRd")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "wordcloud")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "wordcloud", word_type = "feature")
#'
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "network")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "network", id_use = c("GO:0050678", "GO:0035270", "GO:0090276", "GO:0030073"))
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "network", network_layoutadjust = FALSE)
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "network",
#'   topTerm = 4, network_blendmode = "average",
#'   theme_use = "theme_blank", theme_args = list(add_coord = FALSE)
#' ) %>% panel_fix(height = 5)
#'
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "enrichmap")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "enrichmap", enrichmap_expand = c(2, 1))
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "enrichmap",
#'   enrichmap_show_keyword = TRUE, character_width = 10
#' )
#' EnrichmentPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = "Ductal", plot_type = "enrichmap",
#'   topTerm = 200, enrichmap_mark = "hull", enrichmap_label = "feature", enrlichmap_nlabel = 3, character_width = 10,
#'   theme_use = "theme_blank", theme_args = list(add_coord = FALSE)
#' ) %>% panel_fix(height = 4)
#'
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, db = c("MP", "DO"), group_by = "CellType", convert_species = TRUE, species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, db = c("MP", "DO"), group_by = "CellType", group_use = "Endocrine", ncol = 1)
#'
#' @importFrom ggplot2 ggplot geom_bar geom_text geom_label labs scale_fill_manual scale_y_continuous scale_linewidth facet_grid coord_flip scale_color_gradientn scale_fill_gradientn scale_size guides geom_segment expansion guide_colorbar scale_color_manual guide_none draw_key_point  scale_color_identity scale_fill_identity .pt
#' @importFrom dplyr %>% group_by filter arrange desc across mutate slice_head reframe distinct n .data
#' @importFrom stats formula
#' @importFrom patchwork wrap_plots
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid grobTree convertUnit unit
#' @importFrom igraph as_data_frame graph_from_data_frame V layout_with_dh layout_with_drl layout_with_fr layout_with_gem layout_with_graphopt layout_with_kk layout_with_lgl layout_with_mds layout_in_circle layout_as_tree layout_on_grid cluster_fast_greedy cluster_infomap cluster_leiden cluster_louvain cluster_spinglass cluster_walktrap cluster_fluid_communities
#' @export
#'
