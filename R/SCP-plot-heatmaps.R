GroupHeatmap <- function(srt, features = NULL, group.by = NULL, split.by = NULL, within_groups = FALSE, grouping.var = NULL, numerator = NULL, cells = NULL,
                         aggregate_fun = base::mean, exp_cutoff = 0, border = TRUE, flip = FALSE,
                         slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), exp_legend_title = NULL, limits = NULL, lib_normalize = identical(slot, "counts"), libsize = NULL,
                         feature_split = NULL, feature_split_by = NULL, n_split = NULL, split_order = NULL,
                         split_method = c("kmeans", "hclust", "mfuzz"), decreasing = FALSE, fuzzification = NULL,
                         cluster_features_by = NULL, cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                         show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                         row_title = NULL, column_title = NULL, row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                         anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                         terms_width = unit(4, "in"), terms_fontsize = 8,
                         keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                         features_width = unit(2, "in"), features_fontsize = c(6, 10),
                         IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                         db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                         GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                         pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, words_excluded = NULL,
                         nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                         add_bg = FALSE, bg_alpha = 0.5,
                         add_dot = FALSE, dot_size = unit(8, "mm"),
                         add_reticle = FALSE, reticle_color = "grey",
                         add_violin = FALSE, fill.by = "feature", fill_palette = "Dark2", fill_palcolor = NULL,
                         heatmap_palette = "RdBu", heatmap_palcolor = NULL, group_palette = "Paired", group_palcolor = NULL,
                         cell_split_palette = "simspec", cell_split_palcolor = NULL, feature_split_palette = "simspec", feature_split_palcolor = NULL,
                         cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = unit(10, "mm")) else list(height = unit(10, "mm")),
                         feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = if (flip) list(height = unit(5, "mm")) else list(width = unit(5, "mm")),
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
  exp_name <- exp_legend_title %||% exp_name

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
    srt@meta.data[["All.groups"]] <- factor("")
    group.by <- "All.groups"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(group.by[!group.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!is.null(group.by)) {
    for (g in group.by) {
      if (!is.factor(srt@meta.data[[g]])) {
        srt@meta.data[[g]] <- factor(srt@meta.data[[g]], levels = unique(srt@meta.data[[g]]))
      }
    }
  }
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(split.by[!split.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!is.null(split.by)) {
    if (!is.factor(srt@meta.data[[split.by]])) {
      srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], levels = unique(srt@meta.data[[split.by]]))
    }
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
  raw.group.by <- group.by
  raw.group_palette <- group_palette
  if (isTRUE(within_groups)) {
    new.group.by <- c()
    new.group_palette <- group_palette
    for (g in group.by) {
      groups <- split(colnames(srt), srt[[g, drop = TRUE]])
      new.group_palette[g] <- list(rep(new.group_palette[g], length(groups)))
      for (nm in names(groups)) {
        srt[[make.names(nm)]] <- factor(NA, levels = levels(srt[[g, drop = TRUE]]))
        srt[[make.names(nm)]][colnames(srt) %in% groups[[nm]], ] <- nm
        new.group.by <- c(new.group.by, make.names(nm))
      }
    }
    group.by <- new.group.by
    group_palette <- unlist(new.group_palette)
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

  mat_raw <- as_matrix(rbind(slot(srt@assays[[assay]], slot)[gene, cells, drop = FALSE], t(srt@meta.data[cells, meta, drop = FALSE])))[features, , drop = FALSE]
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
    if (cell_group != "All.groups") {
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
      #   lgd[[cell_group]] <- Legend(
      #     title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
      #     legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      #   )
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

  for (i in seq_along(raw.group.by)) {
    cell_group <- raw.group.by[i]
    if (cell_group != "All.groups") {
      lgd[[cell_group]] <- Legend(
        title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = raw.group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      )
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
            stat.by = cellan, group.by = cell_group, split.by = split.by,
            palette = palette, palcolor = palcolor, title = NULL,
            individual = TRUE, combine = FALSE
          )
          subplots_list[[paste0(cellan, ":", cell_group)]] <- subplots
          graphics <- list()
          for (nm in names(subplots)) {
            funbody <- paste0(
              "
              g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
              g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
        if (split_method == "hclust") {
          hc <- hclust(as.dist(dist(mat_split)))
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
      dend <- as.dendrogram(hclust(as.dist(dist(mat_cluster))))
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
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, words_excluded = words_excluded
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
      #         g <- as_grob(subplots_list[['", cellan, ":", cell_group, "']]", "[['", nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
      row_title = row_title %||% if (flip) ifelse(cell_group != "All.groups", cell_group, "") else character(0),
      row_title_side = row_title_side,
      column_title = column_title %||% if (flip) character(0) else ifelse(cell_group != "All.groups", cell_group, ""),
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
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
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

#' FeatureHeatmap
#'
#' @inheritParams GroupHeatmap
#' @param max_cells An integer, maximum number of cells to sample per group, default is 100.
#' @param cell_order A vector of cell names defining the order of cells, default is NULL.
#'
#' @seealso \code{\link{RunDEtest}}
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
#'   ht_params = list(row_gap = unit(0, "mm"))
#' )
#' ht2$plot
#'
#' ht3 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, feature_split = de_filter$group1, group.by = "CellType",
#'   species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
#' ht4 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, n_split = 4, group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   cell_annotation = c("Phase", "G2M_score"), cell_annotation_palette = c("Dark2", "Purples")
#' )
#' ht4$plot
#'
#' ht5 <- FeatureHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, n_split = 4, group.by = "CellType",
#'   heatmap_palette = "viridis",
#'   feature_annotation = c("TF", "CSPA"),
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
#' @importFrom grid unit gpar grid.grabExpr grid.text
#' @importFrom gtable gtable_add_padding
#' @importFrom patchwork wrap_plots
#' @importFrom Seurat GetAssayData
#' @importFrom stats hclust order.dendrogram as.dendrogram reorder sd
#' @importFrom dplyr %>% filter group_by arrange desc across mutate distinct n .data "%>%"
#' @importFrom Matrix t
#' @importFrom proxyC dist
#' @export
FeatureHeatmap <- function(srt, features = NULL, cells = NULL, group.by = NULL, split.by = NULL, within_groups = FALSE, max_cells = 100, cell_order = NULL, border = TRUE, flip = FALSE,
                           slot = "counts", assay = NULL, exp_method = c("zscore", "raw", "fc", "log2fc", "log1p"), exp_legend_title = NULL, limits = NULL,
                           lib_normalize = identical(slot, "counts"), libsize = NULL,
                           feature_split = NULL, feature_split_by = NULL, n_split = NULL, split_order = NULL,
                           split_method = c("kmeans", "hclust", "mfuzz"), decreasing = FALSE, fuzzification = NULL,
                           cluster_features_by = NULL, cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = ifelse(flip, "left", "right"), column_names_side = ifelse(flip, "bottom", "top"), row_names_rot = 0, column_names_rot = 90,
                           row_title = NULL, column_title = NULL, row_title_side = "left", column_title_side = "top", row_title_rot = 0, column_title_rot = ifelse(flip, 90, 0),
                           anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                           terms_width = unit(4, "in"), terms_fontsize = 8,
                           keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                           features_width = unit(2, "in"), features_fontsize = c(6, 10),
                           IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                           db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                           GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                           pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, words_excluded = NULL,
                           nlabel = 20, features_label = NULL, label_size = 10, label_color = "black",
                           heatmap_palette = "RdBu", heatmap_palcolor = NULL, group_palette = "Paired", group_palcolor = NULL,
                           cell_split_palette = "simspec", cell_split_palcolor = NULL, feature_split_palette = "simspec", feature_split_palcolor = NULL,
                           cell_annotation = NULL, cell_annotation_palette = "Paired", cell_annotation_palcolor = NULL, cell_annotation_params = if (flip) list(width = unit(5, "mm")) else list(height = unit(5, "mm")),
                           feature_annotation = NULL, feature_annotation_palette = "Dark2", feature_annotation_palcolor = NULL, feature_annotation_params = if (flip) list(height = unit(5, "mm")) else list(width = unit(5, "mm")),
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
  exp_name <- exp_legend_title %||% exp_name

  assay <- assay %||% DefaultAssay(srt)
  if (length(feature_split) != 0 && length(feature_split) != length(features)) {
    stop("feature_split must be the same length as features")
  }
  if (is.null(group.by)) {
    srt@meta.data[["All.groups"]] <- factor("")
    group.by <- "All.groups"
  }
  if (any(!group.by %in% colnames(srt@meta.data))) {
    stop(group.by[!group.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!is.null(group.by)) {
    for (g in group.by) {
      if (!is.factor(srt@meta.data[[g]])) {
        srt@meta.data[[g]] <- factor(srt@meta.data[[g]], levels = unique(srt@meta.data[[g]]))
      }
    }
  }
  if (length(split.by) > 1) {
    stop("'split.by' only support one variable.")
  }
  if (any(!split.by %in% colnames(srt@meta.data))) {
    stop(split.by[!split.by %in% colnames(srt@meta.data)], " is not in the meta data of the Seurat object.")
  }
  if (!is.null(split.by)) {
    if (!is.factor(srt@meta.data[[split.by]])) {
      srt@meta.data[[split.by]] <- factor(srt@meta.data[[split.by]], levels = unique(srt@meta.data[[split.by]]))
    }
  }
  if (length(group_palette) == 1) {
    group_palette <- rep(group_palette, length(group.by))
  }
  if (length(group_palette) != length(group.by)) {
    stop("'group_palette' must be the same length as 'group.by'")
  }
  group_palette <- setNames(group_palette, nm = group.by)
  raw.group.by <- group.by
  raw.group_palette <- group_palette
  if (isTRUE(within_groups)) {
    new.group.by <- c()
    new.group_palette <- group_palette
    for (g in group.by) {
      groups <- split(colnames(srt), srt[[g, drop = TRUE]])
      new.group_palette[g] <- list(rep(new.group_palette[g], length(groups)))
      for (nm in names(groups)) {
        srt[[make.names(nm)]] <- factor(NA, levels = levels(srt[[g, drop = TRUE]]))
        srt[[make.names(nm)]][colnames(srt) %in% groups[[nm]], ] <- nm
        new.group.by <- c(new.group.by, make.names(nm))
      }
    }
    group.by <- new.group.by
    group_palette <- unlist(new.group_palette)
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
  mat_raw <- as_matrix(rbind(slot(srt@assays[[assay]], slot)[gene, all_cells, drop = FALSE], t(srt@meta.data[all_cells, meta, drop = FALSE])))[features, , drop = FALSE]
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
    if (cell_group != "All.groups") {
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
      # lgd[[cell_group]] <- Legend(
      #   title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
      #   legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      # )
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
  for (i in seq_along(raw.group.by)) {
    cell_group <- raw.group.by[i]
    if (cell_group != "All.groups") {
      lgd[[cell_group]] <- Legend(
        title = cell_group, labels = levels(srt@meta.data[[cell_group]]),
        legend_gp = gpar(fill = palette_scp(levels(srt@meta.data[[cell_group]]), palette = raw.group_palette[i], palcolor = group_palcolor[[i]])), border = TRUE
      )
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
        if (split_method == "hclust") {
          hc <- hclust(as.dist(dist(mat_split)))
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
      dend <- as.dendrogram(hclust(as.dist(dist(mat_cluster))))
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
    pvalueCutoff = pvalueCutoff, padjustCutoff = padjustCutoff, topTerm = topTerm, show_termid = show_termid, topWord = topWord, words_excluded = words_excluded
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
      row_title = row_title %||% if (flip) ifelse(cell_group != "All.groups", cell_group, "") else character(0),
      row_title_side = row_title_side,
      column_title = column_title %||% if (flip) character(0) else ifelse(cell_group != "All.groups", cell_group, ""),
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
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
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
#' This function generates a heatmap to visualize the similarity between different cell types or conditions. It takes in Seurat objects or expression matrices as input and calculates pairwise similarities or distance.
#'
#' @param srt_query A Seurat object or count matrix representing the query dataset. This dataset will be used to calculate the similarities between cells.
#' @param srt_ref A Seurat object or count matrix representing the reference dataset. If provided, the similarities will be calculated between cells from the query and reference datasets. If not provided, the similarities will be calculated within the query dataset.
#' @param bulk_ref A count matrix representing bulk data. If provided, the similarities will be calculated between cells from the query dataset and bulk data.
#' @param query_group The grouping variable in the query dataset. This variable will be used to group cells in the heatmap rows. If not provided, all cells will be treated as one group.
#' @param ref_group The grouping variable in the reference dataset. This variable will be used to group cells in the heatmap columns. If not provided, all cells will be treated as one group.
#' @param query_assay The assay to use for the query dataset. If not provided, the default assay of the query dataset will be used.
#' @param ref_assay The assay to use for the reference dataset. If not provided, the default assay of the reference dataset will be used.
#' @param query_reduction The dimensionality reduction method to use for the query dataset. If not provided, no dimensionality reduction will be applied to the query dataset.
#' @param ref_reduction The dimensionality reduction method to use for the reference dataset. If not provided, no dimensionality reduction will be applied to the reference dataset.
#' @param query_dims The dimensions to use for the query dataset. If not provided, the first 30 dimensions will be used.
#' @param ref_dims The dimensions to use for the reference dataset. If not provided, the first 30 dimensions will be used.
#' @param query_collapsing Whether to collapse cells within each query group before calculating similarities. If set to TRUE, the similarities will be calculated between query groups rather than individual cells.
#' @param ref_collapsing Whether to collapse cells within each reference group before calculating similarities. If set to TRUE, the similarities will be calculated between reference groups rather than individual cells.
#' @param features A vector of feature names to include in the heatmap. If not provided, a default set of highly variable features (HVF) will be used.
#' @param features_type The type of features to use. Options are "HVF" for highly variable features, "DE" for differentially expressed features between query and reference groups.
#' @param feature_source The source of features to use. Options are "query" to use only features from the query dataset, "ref" to use only features from the reference dataset, or "both" to use features from both datasets. If not provided or set to "both", features will be selected from both datasets.
#' @param nfeatures The maximum number of features to include in the heatmap. If not provided, the default is 2000.
#' @param DEtest_param The parameters to use for differential expression testing. This should be a list with two elements: "max.cells.per.ident" specifying the maximum number of cells per group for differential expression testing, and "test.use" specifying the statistical test to use for differential expression testing. If not provided, the default parameters will be used.
#' @param DE_threshold The threshold for differential expression. Only features with adjusted p-values below this threshold will be considered differentially expressed.
#' @param distance_metric The distance metric to use for calculating similarities between cells. This can be any of the following: "cosine", "pearson", "spearman", "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman", "simple matching", or "faith". If not provided, the default is "cosine".
#' @param k The number of nearest neighbors to use for calculating similarities. If not provided, the default is 30.
#' @param filter_lowfreq The minimum frequency threshold for selecting query dataset features. Features with a frequency below this threshold will be excluded from the heatmap. If not provided, the default is 0.
#' @param prefix The prefix to use for the KNNPredict tool slot in the query dataset. This can be used to avoid conflicts with other tools in the Seurat object. If not provided, the default is "KNNPredict".
#' @param exp_legend_title The title for the color legend in the heatmap. If not provided, a default title based on the similarity metric will be used.
#' @param border Whether to add a border around each heatmap cell. If not provided, the default is TRUE.
#' @param flip Whether to flip the orientation of the heatmap. If set to TRUE, the rows and columns of the heatmap will be swapped. This can be useful for visualizing large datasets in a more compact form. If not provided, the default is FALSE.
#' @param limits The limits for the color scale in the heatmap. If not provided, the default is to use the range of similarity values.
#' @param cluster_rows Whether to cluster the rows of the heatmap. If set to TRUE, the rows will be rearranged based on hierarchical clustering. If not provided, the default is FALSE.
#' @param cluster_columns Whether to cluster the columns of the heatmap. If set to TRUE, the columns will be rearranged based on hierarchical clustering. If not provided, the default is FALSE.
#' @param show_row_names Whether to show the row names in the heatmap. If not provided, the default is FALSE.
#' @param show_column_names Whether to show the column names in the heatmap. If not provided, the default is FALSE.
#' @param row_names_side The side of the heatmap to show the row names. Options are "left" or "right". If not provided, the default is "left".
#' @param column_names_side The side of the heatmap to show the column names. Options are "top" or "bottom". If not provided, the default is "top".
#' @param row_names_rot The rotation angle of the row names. If not provided, the default is 0 degrees.
#' @param column_names_rot The rotation angle of the column names. If not provided, the default is 90 degrees.
#' @param row_title The title for the row names in the heatmap. If not provided, the default is to use the query grouping variable.
#' @param column_title The title for the column names in the heatmap. If not provided, the default is to use the reference grouping variable.
#' @param row_title_side The side of the heatmap to show the row title. Options are "top" or "bottom". If not provided, the default is "left".
#' @param column_title_side The side of the heatmap to show the column title. Options are "left" or "right". If not provided, the default is "top".
#' @param row_title_rot The rotation angle of the row title. If not provided, the default is 90 degrees.
#' @param column_title_rot The rotation angle of the column title. If not provided, the default is 0 degrees.
#' @param nlabel The maximum number of labels to show on each side of the heatmap. If set to 0, no labels will be shown. This can be useful for reducing clutter in large heatmaps. If not provided, the default is 0.
#' @param label_cutoff The similarity cutoff for showing labels. Only cells with similarity values above this cutoff will have labels. If not provided, the default is 0.
#' @param label_by The dimension to use for labeling cells. Options are "row" to label cells by row, "column" to label cells by column, or "both" to label cells by both row and column. If not provided, the default is "row".
#' @param label_size The size of the labels in points. If not provided, the default is 10.
#' @param heatmap_palette The color palette to use for the heatmap. This can be any of the palettes available in the circlize package. If not provided, the default is "RdBu".
#' @param heatmap_palcolor The specific colors to use for the heatmap palette. This should be a vector of color names or RGB values. If not provided, the default is NULL.
#' @param query_group_palette The color palette to use for the query group legend. This can be any of the palettes available in the circlize package. If not provided, the default is "Paired".
#' @param query_group_palcolor The specific colors to use for the query group palette. This should be a vector of color names or RGB values. If not provided, the default is NULL.
#' @param ref_group_palette The color palette to use for the reference group legend. This can be any of the palettes available in the circlize package. If not provided, the default is "simspec".
#' @param ref_group_palcolor The specific colors to use for the reference group palette. This should be a vector of color names or RGB values. If not provided, the default is NULL.
#' @param query_cell_annotation A vector of cell metadata column names or assay feature names to use for highlighting specific cells in the heatmap. Each element of the vector will create a separate cell annotation track in the heatmap. If not provided, no cell annotations will be shown.
#' @param query_cell_annotation_palette The color palette to use for the query cell annotation tracks. This can be any of the palettes available in the circlize package. If a single color palette is provided, it will be used for all cell annotation tracks. If multiple color palettes are provided, each track will be assigned a separate palette. If not provided, the default is "Paired".
#' @param query_cell_annotation_palcolor The specific colors to use for the query cell annotation palettes. This should be a list of vectors, where each vector contains the colors for a specific cell annotation track. If a single color vector is provided, it will be used for all cell annotation tracks. If multiple color vectors are provided, each track will be assigned a separate color vector. If not provided, the default is NULL.
#' @param query_cell_annotation_params Additional parameters for customizing the appearance of the query cell annotation tracks. This should be a list with named elements, where the names correspond to parameter names in the heatmaps_annotation() function from the ComplexHeatmap package. If not provided, the default parameters will be used.
#' @param ref_cell_annotation A vector of cell metadata column names or assay feature names to use for highlighting specific cells in the heatmap. Each element of the vector will create a separate cell annotation track in the heatmap. If not provided, no cell annotations will be shown.
#' @param ref_cell_annotation_palette The color palette to use for the reference cell annotation tracks. This can be any of the palettes available in the circlize package. If a single color palette is provided, it will be used for all cell annotation tracks. If multiple color palettes are provided, each track will be assigned a separate palette. If not provided, the default is "Paired".
#' @param ref_cell_annotation_palcolor The specific colors to use for the reference cell annotation palettes. This should be a list of vectors, where each vector contains the colors for a specific cell annotation track. If a single color vector is provided, it will be used for all cell annotation tracks. If multiple color vectors are provided, each track will be assigned a separate color vector. If not provided, the default is NULL.
#' @param ref_cell_annotation_params Additional parameters for customizing the appearance of the reference cell annotation tracks. This should be a list with named elements, where the names correspond to parameter names in the heatmaps_annotation() function from the ComplexHeatmap package. If not provided, the default parameters will be used.
#' @param use_raster Whether to use raster images for rendering the heatmap. If set to TRUE, the heatmap will be rendered as a raster image using the raster_device argument. If not provided, the default is determined based on the number of rows and columns in the heatmap.
#' @param raster_device The raster device to use for rendering the heatmap. This should be a character string specifying the device name, such as "png", "jpeg", or "pdf". If not provided, the default is "png".
#' @param raster_by_magick Whether to use the magick package for rendering rasters. If set to TRUE, the magick package will be used instead of the raster package. This can be useful for rendering large heatmaps more efficiently. If the magick package is not installed, this argument will be ignored.
#' @param width The width of the heatmap in the specified units. If not provided, the width will be automatically determined based on the number of columns in the heatmap and the default unit.
#' @param height The height of the heatmap in the specified units. If not provided, the height will be automatically determined based on the number of rows in the heatmap and the default unit.
#' @param units The units to use for the width and height of the heatmap. Options are "mm", "cm", or "inch". If not provided, the default is "inch".
#' @param seed The random seed to use for reproducible results. If not provided, the default is 11.
#' @param ht_params Additional parameters to customize the appearance of the heatmap. This should be a list with named elements, where the names correspond to parameter names in the Heatmap() function from the ComplexHeatmap package. Any conflicting parameters will override the defaults set by this function.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item{\code{plot}}{The heatmap plot as a ggplot object.}
#'     \item{\code{features}}{The features used in the heatmap.}
#'     \item{\code{simil_matrix}}{The similarity matrix used to generate the heatmap.}
#'     \item{\code{simil_name}}{The name of the similarity metric used to generate the heatmap.}
#'     \item{\code{cell_metadata}}{The cell metadata used to generate the heatmap.}
#'   }
#'
#' @seealso \code{\link{RunKNNMap}} \code{\link{RunKNNPredict}}
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
#'   width = 4, height = 4
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
                           filter_lowfreq = 0, prefix = "KNNPredict", exp_legend_title = NULL,
                           border = TRUE, flip = FALSE, limits = NULL,
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           show_row_names = FALSE, show_column_names = FALSE, row_names_side = "left", column_names_side = "top", row_names_rot = 0, column_names_rot = 90,
                           row_title = NULL, column_title = NULL, row_title_side = "left", column_title_side = "top", row_title_rot = 90, column_title_rot = 0,
                           nlabel = 0, label_cutoff = 0, label_by = "row", label_size = 10,
                           heatmap_palette = "RdBu", heatmap_palcolor = NULL,
                           query_group_palette = "Paired", query_group_palcolor = NULL,
                           ref_group_palette = "simspec", ref_group_palcolor = NULL,
                           query_cell_annotation = NULL, query_cell_annotation_palette = "Paired", query_cell_annotation_palcolor = NULL, query_cell_annotation_params = if (flip) list(height = unit(10, "mm")) else list(width = unit(10, "mm")),
                           ref_cell_annotation = NULL, ref_cell_annotation_palette = "Paired", ref_cell_annotation_palcolor = NULL, ref_cell_annotation_params = if (flip) list(width = unit(10, "mm")) else list(height = unit(10, "mm")),
                           use_raster = NULL, raster_device = "png", raster_by_magick = FALSE, height = NULL, width = NULL, units = "inch",
                           seed = 11, ht_params = list()) {
  set.seed(seed)
  if (isTRUE(raster_by_magick)) {
    check_R("magick")
  }

  ref_legend <- TRUE
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
    ref_legend <- FALSE
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
    simil_matrix <- t(as_matrix(1 - distance_matrix))
    simil_name <- paste0(capitalize(distance_metric), " similarity")
  } else if (distance_metric %in% dist_method) {
    simil_matrix <- t(as_matrix(1 - distance_matrix / max(distance_matrix, na.rm = TRUE)))
    simil_name <- paste0("1-dist[", distance_metric, "]/max(dist[", distance_metric, "])")
  }
  simil_matrix[is.infinite(simil_matrix)] <- max(abs(simil_matrix[!is.infinite(simil_matrix)]), na.rm = TRUE) * ifelse(simil_matrix[is.infinite(simil_matrix)] > 0, 1, -1)
  simil_matrix[is.na(simil_matrix)] <- 0
  exp_name <- exp_legend_title %||% simil_name

  cell_groups <- list()
  if (is.null(query_group)) {
    srt_query@meta.data[["All.groups"]] <- factor("")
    query_group <- "All.groups"
  }
  if (is.null(ref_group)) {
    srt_ref@meta.data[["All.groups"]] <- factor("")
    ref_group <- "All.groups"
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
  lgd[["ht"]] <- Legend(title = exp_name, col_fun = colors, border = TRUE)

  ha_query_list <- NULL
  if (query_group != "All.groups") {
    if (isFALSE(query_collapsing) && ((isFALSE(flip) & isTRUE(cluster_rows)) || (isTRUE(flip) & isTRUE(cluster_columns)))) {
      query_cell_annotation <- c(query_group, query_cell_annotation)
      query_cell_annotation_palette <- c(query_group_palette, query_cell_annotation_palette)
      query_cell_annotation_palcolor <- c(list(query_group_palcolor), query_cell_annotation_palcolor %||% list(NULL))
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
  if (ref_group != "All.groups") {
    if (isFALSE(ref_collapsing) && ((isFALSE(flip) & isTRUE(cluster_columns)) || (isTRUE(flip) & isTRUE(cluster_rows)))) {
      ref_cell_annotation <- c(ref_group, ref_cell_annotation)
      ref_cell_annotation_palette <- c(ref_group_palette, ref_cell_annotation_palette)
      ref_cell_annotation_palcolor <- c(list(ref_group_palcolor), ref_cell_annotation_palcolor %||% list(NULL))
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
              g <- as_grob(query_subplots_list[['", cellan, ":", query_group, "']]", "[['", nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
              g <- as_grob(query_subplots_list[['", cellan, ":", query_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
              g <- as_grob(ref_subplots_list[['", cellan, ":", ref_group, "']]", "[['", nm, "']] + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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
              g <- as_grob(ref_subplots_list[['", cellan, ":", ref_group, "']]", "[['", nm, "']]  + facet_null() + theme_void() + theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = 'none'));
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

  if (!isTRUE(ref_legend)) {
    lgd <- lgd[grep("Ref", names(lgd), invert = TRUE)]
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
    name = exp_name,
    matrix = if (flip) t(simil_matrix) else simil_matrix,
    col = colors,
    layer_fun = layer_fun,
    row_title = row_title %||% if (flip) paste0(c("Ref", ref_group), collapse = ":") else paste0(c("Query", query_group), collapse = ":"),
    row_title_side = row_title_side,
    column_title = column_title %||% if (flip) paste0(c("Query", query_group), collapse = ":") else paste0(c("Ref", ref_group), collapse = ":"),
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
    p <- panel_fix_overall(gTree, width = as.numeric(ht_width), height = as.numeric(ht_height), units = units)
  } else {
    p <- wrap_plots(gTree)
  }

  return(list(
    plot = p,
    features = features,
    simil_matrix = simil_matrix,
    simil_name = simil_name,
    cell_metadata = cell_metadata
  ))
}

#' Heatmap plot for dynamic features along lineages
#'
#' @inheritParams GroupHeatmap
#' @param srt A Seurat object.
#' @param lineages A character vector specifying the lineages to plot.
#' @param features A character vector specifying the features to plot. By default, this parameter is set to NULL, and the dynamic features will be determined by the parameters  \code{min_expcells}, \code{r.sq}, \code{dev.expl}, \code{padjust} and \code{num_intersections}.
#' @param use_fitted A logical indicating whether to use fitted values. Default is FALSE.
#' @param border A logical indicating whether to add a border to the heatmap. Default is TRUE.
#' @param flip A logical indicating whether to flip the heatmap. Default is FALSE.
#' @param min_expcells A numeric value specifying the minimum number of expected cells. Default is 20.
#' @param r.sq A numeric value specifying the R-squared threshold. Default is 0.2.
#' @param dev.expl A numeric value specifying the deviance explained threshold. Default is 0.2.
#' @param padjust A numeric value specifying the p-value adjustment threshold. Default is 0.05.
#' @param num_intersections This parameter is a numeric vector used to determine the number of intersections among lineages. It helps in selecting which dynamic features will be used. By default, when this parameter is set to NULL, all dynamic features that pass the specified threshold will be used for each lineage.
#' @param cell_density A numeric value is used to define the cell density within each cell bin. By default, this parameter is set to 1, which means that all cells will be included within each cell bin.
#' @param cell_bins A numeric value specifying the number of cell bins. Default is 100.
#' @param order_by A character vector specifying the order of the heatmap. Default is "peaktime".
#' @param family A character specifying the model used to calculate the dynamic features if needed. By default, this parameter is set to NULL, and the appropriate family will be automatically determined.
#' @param cluster_features_by A character vector specifying which lineage to use when clustering features. By default, this parameter is set to NULL, which means that all lineages will be used.
#' @param pseudotime_label A numeric vector specifying the pseudotime label. Default is NULL.
#' @param pseudotime_label_color A character string specifying the pseudotime label color. Default is "black".
#' @param pseudotime_label_linetype A numeric value specifying the pseudotime label line type. Default is 2.
#' @param pseudotime_label_linewidth A numeric value specifying the pseudotime label line width. Default is 3.
#' @param pseudotime_palette A character vector specifying the color palette to use for pseudotime.
#' @param pseudotime_palcolor A list specifying the colors to use for the pseudotime in the heatmap.
#' @param separate_annotation A character vector of names of annotations to be displayed in separate annotation blocks. Each name should match a column name in the metadata of the Seurat object.
#' @param separate_annotation_palette A character vector specifying the color palette to use for separate annotations.
#' @param separate_annotation_palcolor A list specifying the colors to use for each level of the separate annotations.
#' @param separate_annotation_params A list of other parameters to be passed to the HeatmapAnnotation function when creating the separate annotation blocks.
#' @param reverse_ht A logical indicating whether to reverse the heatmap. Default is NULL.
#'
#' @seealso \code{\link{RunDynamicFeatures}} \code{\link{RunDynamicEnrichment}}
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
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
#' ht4 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   reverse_ht = "Lineage1",
#'   use_fitted = TRUE,
#'   n_split = 6,
#'   split_method = "mfuzz",
#'   heatmap_palette = "viridis",
#'   cell_annotation = c("SubCellType", "Phase", "G2M_score"),
#'   cell_annotation_palette = c("Paired", "simspec", "Purples"),
#'   separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
#'   separate_annotation_palette = c("Paired", "Set1"),
#'   separate_annotation_params = list(height = unit(10, "mm")),
#'   feature_annotation = c("TF", "CSPA"),
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
#'   cell_annotation_palette = c("Paired", "simspec", "Purples"),
#'   separate_annotation = list("SubCellType", c("Nnat", "Irx1")),
#'   separate_annotation_palette = c("Paired", "Set1"),
#'   separate_annotation_params = list(width = unit(10, "mm")),
#'   feature_annotation = c("TF", "CSPA"),
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
#' @importFrom stats kmeans sd
#' @importFrom patchwork wrap_plots
#' @importFrom grid gpar grid.lines grid.text unit
#' @importFrom gtable gtable_add_padding
#' @importFrom Matrix t
#' @importFrom dplyr %>% filter group_by arrange desc across mutate reframe distinct n .data
#' @importFrom proxyC dist
#' @export
