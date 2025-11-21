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
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
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

#' @importFrom stats hclust as.dendrogram order.dendrogram
#' @importFrom proxyC dist
#' @importFrom ComplexHeatmap merge_dendrogram
cluster_within_group2 <- function(mat, factor) {
  check_R("dendextend")
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
      hc1 <- hclust(as.dist(dist(t(m))))
      dend_list[[le]] <- as.dendrogram(hc1)
      order_list[[le]] <- which(factor == le)[order.dendrogram(dend_list[[le]])]
      dendextend::order.dendrogram(dend_list[[le]]) <- order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") <- le
  }
  parent <- as.dendrogram(hclust(as.dist(dist(t(sapply(
    order_list,
    function(x) rowMeans(mat[, x, drop = FALSE])
  ))))))
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
#' @importFrom grid gpar unit
#' @importFrom dplyr %>% filter group_by arrange desc across reframe mutate distinct n .data "%>%"
heatmap_enrichment <- function(geneID, geneID_groups, feature_split_palette = "simspec", feature_split_palcolor = NULL, ha_right = NULL, flip = FALSE,
                               anno_terms = FALSE, anno_keys = FALSE, anno_features = FALSE,
                               terms_width = unit(4, "in"), terms_fontsize = 8,
                               keys_width = unit(2, "in"), keys_fontsize = c(6, 10),
                               features_width = unit(2, "in"), features_fontsize = c(6, 10),
                               IDtype = "symbol", species = "Homo_sapiens", db_update = FALSE, db_combine = FALSE, db_version = "latest", convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                               db = "GO_BP", TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                               GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                               pvalueCutoff = NULL, padjustCutoff = 0.05, topTerm = 5, show_termid = FALSE, topWord = 20, words_excluded = NULL) {
  res <- NULL
  words_excluded <- words_excluded %||% SCP::words_excluded

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
              if (all(df$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
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
                    filter(nchar(.data[["keyword"]]) >= 1) %>%
                    filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) %>%
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
                  filter(nchar(.data[["keyword"]]) >= 1) %>%
                  filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) %>%
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
heatmap_fixsize <- function(width, width_sum, height, height_sum, units, ht_list, legend_list) {
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

#' @importFrom stats sd
standardise <- function(data) {
  data[] <- t(apply(data, 1, scale))
  return(data)
}

mestimate <- function(data) {
  N <- nrow(data)
  D <- ncol(data)
  m.sj <- 1 + (1418 / N + 22.05) * D^(-2) + (12.33 / N + 0.243) *
    D^(-0.0406 * log(N) - 0.1134)
  return(m.sj)
}

#' GroupHeatmap
#'
#' @param srt A Seurat object.
#' @param features The features to include in the heatmap.
#' @param group.by A character vector specifying the groups to group by. Default is NULL.
#' @param split.by A character vector specifying the variable to split the heatmap by. Default is NULL.
#' @param within_groups A logical value indicating whether to create separate heatmap scales for each group or within each group. Default is FALSE.
#' @param grouping.var A character vector that specifies another variable for grouping, such as certain conditions. The default value is NULL.
#' @param numerator A character vector specifying the value to use as the numerator in the grouping.var grouping. Default is NULL.
#' @param cells A character vector specifying the cells to include in the heatmap. Default is NULL.
#' @param aggregate_fun A function to use for aggregating data within groups. Default is base::mean.
#' @param exp_cutoff A numeric value specifying the threshold for cell counting if \code{add_dot} is TRUE. Default is 0.
#' @param border A logical value indicating whether to add a border to the heatmap. Default is TRUE.
#' @param flip A logical value indicating whether to flip the heatmap. Default is FALSE.
#' @param slot A character vector specifying the slot in the Seurat object to use. Default is "counts".
#' @param assay A character vector specifying the assay in the Seurat object to use. Default is NULL.
#' @param exp_method A character vector specifying the method for calculating expression values. Default is "zscore" with options "zscore", "raw", "fc", "log2fc", "log1p".
#' @param exp_legend_title A character vector specifying the title for the legend of expression value. Default is NULL.
#' @param limits A two-length numeric vector specifying the limits for the color scale. Default is NULL.
#' @param lib_normalize A logical value indicating whether to normalize the data by library size.
#' @param libsize A numeric vector specifying the library size for each cell. Default is NULL.
#' @param feature_split A factor specifying how to split the features. Default is NULL.
#' @param feature_split_by A character vector specifying which group.by to use when splitting features (into n_split feature clusters). Default is NULL.
#' @param n_split An integer specifying the number of feature splits (feature clusters) to create. Default is NULL.
#' @param split_order A numeric vector specifying the order of splits. Default is NULL.
#' @param split_method A character vector specifying the method for splitting features. Default is "kmeans" with options "kmeans", "hclust", "mfuzz").
#' @param decreasing A logical value indicating whether to sort feature splits in decreasing order. Default is FALSE.
#' @param fuzzification A numeric value specifying the fuzzification coefficient. Default is NULL.
#' @param cluster_features_by A character vector specifying which group.by to use when clustering features. Default is NULL. By default, this parameter is set to NULL, which means that all groups will be used.
#' @param cluster_rows A logical value indicating whether to cluster rows in the heatmap. Default is FALSE.
#' @param cluster_columns A logical value indicating whether to cluster columns in the heatmap. Default is FALSE.
#' @param cluster_row_slices A logical value indicating whether to cluster row slices in the heatmap. Default is FALSE.
#' @param cluster_column_slices A logical value indicating whether to cluster column slices in the heatmap. Default is FALSE.
#' @param show_row_names A logical value indicating whether to show row names in the heatmap. Default is FALSE.
#' @param show_column_names A logical value indicating whether to show column names in the heatmap. Default is FALSE.
#' @param row_names_side A character vector specifying the side to place row names.
#' @param column_names_side A character vector specifying the side to place column names.
#' @param row_names_rot A numeric value specifying the rotation angle for row names. Default is 0.
#' @param column_names_rot A numeric value specifying the rotation angle for column names. Default is 90.
#' @param row_title A character vector specifying the title for rows. Default is NULL.
#' @param column_title A character vector specifying the title for columns. Default is NULL.
#' @param row_title_side A character vector specifying the side to place row title. Default is "left".
#' @param column_title_side A character vector specifying the side to place column title. Default is "top".
#' @param row_title_rot A numeric value specifying the rotation angle for row title. Default is 0.
#' @param column_title_rot A numeric value specifying the rotation angle for column title.
#' @param anno_terms A logical value indicating whether to include term annotations. Default is FALSE.
#' @param anno_keys A logical value indicating whether to include key annotations. Default is FALSE.
#' @param anno_features A logical value indicating whether to include feature annotations. Default is FALSE.
#' @param terms_width A unit specifying the width of term annotations. Default is unit(4, "in").
#' @param terms_fontsize A numeric vector specifying the font size(s) for term annotations. Default is 8.
#' @param keys_width A unit specifying the width of key annotations. Default is unit(2, "in").
#' @param keys_fontsize A two-length numeric vector specifying the minimum and maximum font size(s) for key annotations. Default is c(6, 10).
#' @param features_width A unit specifying the width of feature annotations. Default is unit(2, "in").
#' @param features_fontsize A two-length numeric vector specifying the minimum and maximum font size(s) for feature annotations. Default is c(6, 10).
#' @param IDtype A character vector specifying the type of IDs for features. Default is "symbol".
#' @param species A character vector specifying the species for features. Default is "Homo_sapiens".
#' @param db_update A logical value indicating whether to update the database. Default is FALSE.
#' @param db_version A character vector specifying the version of the database. Default is "latest".
#' @param db_combine A logical value indicating whether to use a combined database. Default is FALSE.
#' @param convert_species A logical value indicating whether to use a species-converted database if annotation is missing for \code{species}. Default is FALSE.
#' @param Ensembl_version An integer specifying the Ensembl version. Default is 103.
#' @param mirror A character vector specifying the mirror for the Ensembl database. Default is NULL.
#' @param db A character vector specifying the database to use. Default is "GO_BP".
#' @param TERM2GENE A data.frame specifying the TERM2GENE mapping for the database. Default is NULL.
#' @param TERM2NAME A data.frame specifying the TERM2NAME mapping for the database. Default is NULL.
#' @param minGSSize An integer specifying the minimum gene set size for the database. Default is 10.
#' @param maxGSSize An integer specifying the maximum gene set size for the database. Default is 500.
#' @param GO_simplify A logical value indicating whether to simplify gene ontology terms. Default is FALSE.
#' @param GO_simplify_cutoff A character vector specifying the cutoff for GO simplification. Default is "p.adjust < 0.05".
#' @param simplify_method A character vector specifying the method for GO simplification. Default is "Wang".
#' @param simplify_similarityCutoff A numeric value specifying the similarity cutoff for GO simplification. Default is 0.7.
#' @param pvalueCutoff A numeric vector specifying the p-value cutoff(s) for significance. Default is NULL.
#' @param padjustCutoff A numeric value specifying the adjusted p-value cutoff for significance. Default is 0.05.
#' @param topTerm An integer specifying the number of top terms to include. Default is 5.
#' @param show_termid A logical value indicating whether to show term IDs. Default is FALSE.
#' @param topWord An integer specifying the number of top words to include. Default is 20.
#' @param words_excluded A character vector specifying the words to exclude. Default is NULL.
#' @param nlabel An integer specifying the number of labels to include. Default is 0.
#' @param features_label A character vector specifying the features to label. Default is NULL.
#' @param label_size A numeric value specifying the size of labels. Default is 10.
#' @param label_color A character vector specifying the color of labels. Default is "black".
#' @param add_bg A logical value indicating whether to add a background to the heatmap. Default is FALSE.
#' @param bg_alpha A numeric value specifying the alpha value for the background color. Default is 0.5.
#' @param add_dot A logical value indicating whether to add dots to the heatmap. The size of dot represents percentage of expressed cells based on the specified \code{exp_cutoff}. Default is FALSE.
#' @param dot_size A unit specifying the base size of the dots. Default is unit(8, "mm").
#' @param add_reticle A logical value indicating whether to add reticles to the heatmap. Default is FALSE.
#' @param reticle_color A character vector specifying the color of the reticles. Default is "grey".
#' @param add_violin A logical value indicating whether to add violins to the heatmap. Default is FALSE.
#' @param fill.by A character vector specifying what to fill the violin. Possible values are "group", "feature", or "expression". Default is "feature".
#' @param fill_palette A character vector specifying the palette to use for fill. Default is "Dark2".
#' @param fill_palcolor A character vector specifying the fill color to use. Default is NULL.
#' @param heatmap_palette A character vector specifying the palette to use for the heatmap. Default is "RdBu".
#' @param heatmap_palcolor A character vector specifying the heatmap color to use. Default is NULL.
#' @param group_palette A character vector specifying the palette to use for groups. Default is "Paired".
#' @param group_palcolor A character vector specifying the group color to use. Default is NULL.
#' @param cell_split_palette A character vector specifying the palette to use for cell splits. Default is "simspec".
#' @param cell_split_palcolor A character vector specifying the cell split color to use. Default is NULL.
#' @param feature_split_palette A character vector specifying the palette to use for feature splits. Default is "simspec".
#' @param feature_split_palcolor A character vector specifying the feature split color to use. Default is NULL.
#' @param cell_annotation A character vector specifying the cell annotation(s) to include. Default is NULL.
#' @param cell_annotation_palette A character vector specifying the palette to use for cell annotations. The length of the vector should match the number of cell_annotation. Default is "Paired".
#' @param cell_annotation_palcolor A list of character vector specifying the cell annotation color(s) to use. The length of the list should match the number of cell_annotation. Default is NULL.
#' @param cell_annotation_params A list specifying additional parameters for cell annotations. Default is a list with width = unit(1, "cm") if flip is TRUE, else a list with height = unit(1, "cm").
#' @param feature_annotation A character vector specifying the feature annotation(s) to include. Default is NULL.
#' @param feature_annotation_palette A character vector specifying the palette to use for feature annotations. The length of the vector should match the number of feature_annotation. Default is "Dark2".
#' @param feature_annotation_palcolor A list of character vector specifying the feature annotation color to use. The length of the list should match the number of feature_annotation. Default is NULL.
#' @param feature_annotation_params A list specifying additional parameters for feature annotations. Default is an empty list.
#' @param use_raster A logical value indicating whether to use a raster device for plotting. Default is NULL.
#' @param raster_device A character vector specifying the raster device to use. Default is "png".
#' @param raster_by_magick A logical value indicating whether to use the 'magick' package for raster. Default is FALSE.
#' @param height A numeric vector specifying the height(s) of the heatmap body. Default is NULL.
#' @param width A numeric vector specifying the width(s) of the heatmap body. Default is NULL.
#' @param units A character vector specifying the units for the height and width. Default is "inch".
#' @param seed An integer specifying the random seed. Default is 11.
#' @param ht_params A list specifying additional parameters passed to the ComplexHeatmap::Heatmap function. Default is an empty list.
#'
#' @seealso \code{\link{RunDEtest}}
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item{\code{plot}}{The heatmap plot.}
#'     \item{\code{matrix_list}}{A list of matrix for each \code{group.by} used in the heatmap.}
#'     \item{\code{feature_split}}{NULL or a factor if splitting is performed in the heatmap.}
#'     \item{\code{cell_metadata}}{Meta data of cells used to generate the heatmap.}
#'     \item{\code{cell_metadata}}{Meta data of features used to generate the heatmap.}
#'     \item{\code{enrichment}}{NULL or a enrichment result generated by RunEnrichment when any of the parameters \code{anno_terms}, \code{anno_keys}, or \code{anno_features} is set to TRUE.}
#'   }
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
#'   group.by = c("CellType", "SubCellType")
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
#'   cluster_rows = TRUE, cluster_columns = TRUE
#' )
#' ht2$plot
#'
#' ht3 <- GroupHeatmap(
#'   srt = pancreas_sub, features = de_filter$gene, feature_split = de_filter$group1, group.by = "CellType",
#'   species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE
#' )
#' ht3$plot
#'
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", db = c("TF", "CSPA"))
#' de_top <- de_filter %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ht4 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
#'   cell_annotation_params = list(height = unit(10, "mm")),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   add_dot = TRUE, add_bg = TRUE, nlabel = 0, show_row_names = TRUE
#' )
#' ht4$plot
#'
#' ht5 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   heatmap_palette = "YlOrRd",
#'   cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
#'   cell_annotation_params = list(width = unit(10, "mm")),
#'   feature_annotation = c("TF", "CSPA"),
#'   feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#'   add_dot = TRUE, add_bg = TRUE,
#'   flip = TRUE, column_title_rot = 45, nlabel = 0, show_row_names = TRUE
#' )
#' ht5$plot
#'
#' ht6 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   add_violin = TRUE, cluster_rows = TRUE,
#'   nlabel = 0, show_row_names = TRUE
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, feature_split = de_top$group1, group.by = "CellType",
#'   add_violin = TRUE, fill.by = "expression", fill_palette = "Blues", cluster_rows = TRUE,
#'   nlabel = 0, show_row_names = TRUE
#' )
#' ht7$plot
#'
#' ht8 <- GroupHeatmap(pancreas_sub,
#'   features = de_top$gene, group.by = "CellType", split.by = "Phase", n_split = 4,
#'   cluster_rows = TRUE, cluster_columns = TRUE, cluster_row_slices = TRUE, cluster_column_slices = TRUE,
#'   add_dot = TRUE, add_reticle = TRUE, heatmap_palette = "viridis",
#'   nlabel = 0, show_row_names = TRUE,
#'   ht_params = list(row_gap = unit(0, "mm"), row_names_gp = gpar(fontsize = 10))
#' )
#' ht8$plot
#'
#' @importFrom circlize colorRamp2
#' @importFrom stats aggregate formula quantile sd
#' @importFrom ComplexHeatmap Legend HeatmapAnnotation anno_block anno_simple anno_customize Heatmap draw pindex restore_matrix %v%
#' @importFrom grid gpar grid.grabExpr grid.lines grid.rect grid.points grid.draw
#' @importFrom ggplot2 theme_void theme facet_null
#' @importFrom patchwork wrap_plots
#' @importFrom methods getFunction
#' @importFrom dplyr %>% filter group_by arrange desc across mutate distinct n .data "%>%"
#' @importFrom Matrix t
#' @importFrom proxyC dist
#' @export
