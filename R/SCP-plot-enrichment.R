EnrichmentPlot <- function(srt, db = "GO_BP", group_by = NULL, test.use = "wilcox", res = NULL,
                           plot_type = c("bar", "dot", "lollipop", "network", "enrichmap", "wordcloud", "comparison"),
                           split_by = c("Database", "Groups"), color_by = "Database",
                           group_use = NULL, id_use = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                           topTerm = ifelse(plot_type == "enrichmap", 100, 6), compare_only_sig = FALSE,
                           topWord = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = NULL,
                           network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
                           network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
                           enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = c("term", "feature"), enrichmap_labelsize = 5,
                           enrlichmap_nlabel = 4, enrichmap_show_keyword = FALSE, enrichmap_mark = c("ellipse", "hull"), enrichmap_expand = c(0.5, 0.5),
                           character_width = 50, lineheight = 0.5,
                           palette = "Spectral", palcolor = NULL,
                           aspect.ratio = 1, legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  word_type <- match.arg(word_type)
  enrichmap_label <- match.arg(enrichmap_label)
  enrichmap_mark <- match.arg(enrichmap_mark)
  words_excluded <- words_excluded %||% SCP::words_excluded

  if (any(!split_by %in% c("Database", "Groups"))) {
    stop("'split_by' must be either 'Database', 'Groups', or both of them")
  }
  if (plot_type %in% c("network", "enrichmap") & length(split_by) == 1) {
    warning("When 'plot_type' is 'network' or 'enrichmap', the 'split_by' parameter does not take effect.", immediate. = TRUE)
    split_by <- c("Database", "Groups")
  }

  if (is.null(res)) {
    if (is.null(group_by)) {
      stop("'group_by' must be provided.")
    }
    slot <- paste("Enrichment", group_by, test.use, sep = "_")
    if (!slot %in% names(srt@tools)) {
      stop("No enrichment result found. You may perform RunEnrichment first.")
    }
    enrichment <- srt@tools[[slot]][["enrichment"]]
  } else {
    enrichment <- res[["enrichment"]]
  }

  if (is.null(pvalueCutoff) && is.null(padjustCutoff)) {
    stop("One of 'pvalueCutoff' or 'padjustCutoff' must be specified")
  }
  if (!is.factor(enrichment["Groups"])) {
    enrichment[["Groups"]] <- factor(enrichment[["Groups"]], levels = unique(enrichment[["Groups"]]))
  }
  if (length(db[!db %in% enrichment[["Database"]]]) > 0) {
    stop(paste0(db[!db %in% enrichment[["Database"]]], " is not in the enrichment result."))
  }
  if (!is.factor(enrichment[["Database"]])) {
    enrichment[["Database"]] <- factor(enrichment[["Database"]], levels = unique(enrichment[["Database"]]))
  }
  if (!is.null(group_use)) {
    enrichment <- enrichment[enrichment[["Groups"]] %in% group_use, , drop = FALSE]
  }
  if (length(id_use) > 0) {
    topTerm <- Inf
    if (is.list(id_use)) {
      if (is.null(names(id_use))) {
        stop("'id_use' must be named when it is a list.")
      }
      if (!all(names(id_use) %in% enrichment[["Groups"]])) {
        stop(paste0("Names in 'id_use' is invalid: ", paste0(names(id_use)[!names(id_use) %in% enrichment[["Groups"]]], collapse = ",")))
      }
      enrichment_list <- list()
      for (i in seq_along(id_use)) {
        enrichment_list[[i]] <- enrichment[enrichment[["ID"]] %in% id_use[[i]] & enrichment[["Groups"]] %in% names(id_use)[i], , drop = FALSE]
      }
      enrichment <- do.call(rbind, enrichment_list)
    } else {
      enrichment <- enrichment[enrichment[["ID"]] %in% unlist(id_use), , drop = FALSE]
    }
  }

  metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
  metric_value <- ifelse(is.null(padjustCutoff), pvalueCutoff, padjustCutoff)

  pvalueCutoff <- ifelse(is.null(pvalueCutoff), Inf, pvalueCutoff)
  padjustCutoff <- ifelse(is.null(padjustCutoff), Inf, padjustCutoff)

  if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
    enrichment_sim <- enrichment[enrichment[["Database"]] %in% gsub("_sim", "", db), , drop = FALSE]
  }
  enrichment <- enrichment[enrichment[["Database"]] %in% db, , drop = FALSE]

  enrichment_sig <- enrichment[enrichment[[metric]] < metric_value | enrichment[["ID"]] %in% unlist(id_use), , drop = FALSE]
  enrichment_sig <- enrichment_sig[order(enrichment_sig[[metric]]), , drop = FALSE]
  if (nrow(enrichment_sig) == 0) {
    stop(
      "No term enriched using the threshold: ",
      paste0("pvalueCutoff = ", pvalueCutoff), "; ",
      paste0("padjustCutoff = ", padjustCutoff)
    )
  }
  df_list <- split(enrichment_sig, formula(paste0("~", split_by, collapse = "+")))
  df_list <- df_list[lapply(df_list, nrow) > 0]

  facet <- switch(paste0(split_by, collapse = "~"),
    "Groups" = formula(paste0("Database ~ Groups")),
    "Database" = formula(paste0("Groups ~ Database")),
    formula(paste0(split_by, collapse = "~"))
  )

  if (plot_type == "comparison") {
    # comparison -------------------------------------------------------------------------------------------------
    ids <- NULL
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)
      ids <- unique(c(ids, df[, "ID"]))
    }
    if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
      enrichment_sub <- subset(enrichment_sim, ID %in% ids)
      enrichment_sub[["Database"]][enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")] <- paste0(enrichment_sub[["Database"]][enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")], "_sim")
    } else {
      enrichment_sub <- subset(enrichment, ID %in% ids)
    }
    enrichment_sub[["Database"]] <- factor(enrichment_sub[["Database"]], levels = db)
    enrichment_sub[["GeneRatio"]] <- sapply(enrichment_sub[["GeneRatio"]], function(x) {
      sp <- strsplit(x, "/")[[1]]
      GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
    })
    enrichment_sub[["BgRatio"]] <- sapply(enrichment_sub[["BgRatio"]], function(x) {
      sp <- strsplit(x, "/")[[1]]
      BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      return(BgRatio)
    })
    enrichment_sub[["EnrichmentScore"]] <- enrichment_sub[["GeneRatio"]] / enrichment_sub[["BgRatio"]]
    enrichment_sub[["Description"]] <- capitalize(enrichment_sub[["Description"]])
    enrichment_sub[["Description"]] <- str_wrap(enrichment_sub[["Description"]], width = character_width)
    terms <- setNames(enrichment_sub[["Description"]], enrichment_sub[["ID"]])
    enrichment_sub[["Description"]] <- factor(enrichment_sub[["Description"]], levels = unique(rev(terms[ids])))
    if (isTRUE(compare_only_sig)) {
      enrichment_sub <- enrichment_sub[enrichment_sub[[metric]] < metric_value, , drop = FALSE]
    }
    p <- ggplot(enrichment_sub, aes(x = Groups, y = Description)) +
      geom_point(aes(size = GeneRatio, fill = .data[[metric]], color = ""), shape = 21) +
      scale_size_area(name = "GeneRatio", max_size = 6, n.breaks = 4) +
      guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
      scale_fill_gradientn(
        name = paste0(metric),
        limits = c(0, min(metric_value, 1)),
        n.breaks = 3,
        colors = palette_scp(palette = palette, palcolor = palcolor, reverse = TRUE),
        na.value = "grey80",
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 2)
      ) +
      scale_color_manual(values = NA, na.value = "black") +
      guides(colour = if (isTRUE(compare_only_sig)) guide_none() else guide_legend("Non-sig", override.aes = list(colour = "black", fill = "grey80", size = 3))) +
      facet_grid(Database ~ ., scales = "free") +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction,
        panel.grid.major = element_line(colour = "grey80", linetype = 2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(
          lineheight = lineheight, hjust = 1,
          face = ifelse(grepl("\n", levels(enrichment_sub[["Description"]])), "italic", "plain")
        )
      )
    plist <- list(p)
  } else if (plot_type == "bar") {
    # bar -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = unique(rev(df[["Description"]])))

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["metric"]],
        fill = .data[[color_by]], label = .data[["Count"]]
      )) +
        geom_bar(width = 0.9, stat = "identity", color = "black") +
        geom_text(hjust = -0.5, size = 3.5, color = "white", fontface = "bold") +
        geom_text(hjust = -0.5, size = 3.5) +
        labs(x = "", y = paste0("-log10(", metric, ")")) +
        scale_fill_manual(
          values = palette_scp(levels(df[[color_by]]), palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = "none"
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(df[["metric"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.y = element_text(
            lineheight = lineheight, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        )
      return(p)
    }))
  } else if (plot_type == "dot") {
    # dot -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df <- df[order(df[["GeneRatio"]], decreasing = TRUE), ]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = unique(rev(df[["Description"]])))

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["GeneRatio"]]
      )) +
        geom_point(aes(fill = .data[["metric"]], size = .data[["Count"]]), color = "black", shape = 21) +
        labs(x = "", y = "GeneRatio") +
        scale_size(name = "Count", range = c(3, 6), scales::breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = palette_scp(palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        scale_y_continuous(limits = c(0, 1.3 * max(df[["GeneRatio"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(
            lineheight = lineheight, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        )
      return(p)
    }))
  } else if (plot_type == "lollipop") {
    # lollipop -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["GeneRatio"]] <- sapply(df[["GeneRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        GeneRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
      })
      df[["BgRatio"]] <- sapply(df[["BgRatio"]], function(x) {
        sp <- strsplit(x, "/")[[1]]
        BgRatio <- as.numeric(sp[1]) / as.numeric(sp[2])
        return(BgRatio)
      })
      df[["FoldEnrichment"]] <- df[["GeneRatio"]] / df[["BgRatio"]]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = unique(df[order(df[["FoldEnrichment"]]), "Description"]))

      p <- ggplot(df, aes(
        x = .data[["Description"]], y = .data[["FoldEnrichment"]],
        fill = .data[["metric"]]
      )) +
        geom_blank() +
        geom_segment(
          aes(y = 0, xend = .data[["Description"]], yend = .data[["FoldEnrichment"]]),
          color = "black", linewidth = 2
        ) +
        geom_segment(
          aes(y = 0, xend = .data[["Description"]], yend = .data[["FoldEnrichment"]], color = .data[["metric"]]),
          linewidth = 1
        ) +
        geom_point(aes(size = .data[["GeneRatio"]]), shape = 21, color = "black") +
        scale_size(name = "GeneRatio", range = c(3, 6), scales::breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_y_continuous(limits = c(0, 1.2 * max(df[["FoldEnrichment"]], na.rm = TRUE)), expand = expansion(0, 0)) +
        labs(x = "", y = "Fold Enrichment") +
        scale_fill_gradientn(
          name = paste0("-log10(", metric, ")"),
          n.breaks = 3,
          colors = palette_scp(palette = palette, palcolor = palcolor),
          na.value = "grey80",
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0),
          aesthetics = c("color", "fill")
        ) +
        facet_grid(facet, scales = "free") +
        coord_flip() +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          panel.grid.major = element_line(colour = "grey80", linetype = 2),
          axis.text.y = element_text(
            lineheight = lineheight, hjust = 1,
            face = ifelse(grepl("\n", levels(df[["Description"]])), "italic", "plain")
          )
        )
      return(p)
    }))
  } else if (plot_type == "network") {
    # network -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = unique(df[["Description"]]))
      df$geneID <- strsplit(df$geneID, "/")
      df_unnest <- unnest(df, cols = "geneID")

      nodes <- rbind(
        data.frame("ID" = df[["Description"]], class = "term", metric = df[["metric"]]),
        data.frame("ID" = unique(df_unnest$geneID), class = "gene", metric = 0)
      )
      nodes$Database <- df$Database[1]
      nodes$Groups <- df$Groups[1]
      edges <- as.data.frame(df_unnest[, c("Description", "geneID")])
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- 1
      graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
      if (network_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(network_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", network_layout), list(graph))
      }
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      if (isTRUE(network_layoutadjust)) {
        width <- nchar(df_nodes$name)
        width[df_nodes$class == "term"] <- 8
        layout <- adjustlayout(
          graph = graph, layout = layout, width = width, height = 2,
          scale = network_adjscale, iter = network_adjiter
        )
      }
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      colors <- palette_scp(levels(df[["Description"]]), palette = palette, palcolor = palcolor)
      df_edges[["color"]] <- colors[df_edges$from]
      node_colors <- aggregate(df_unnest$Description, by = list(df_unnest$geneID), FUN = function(x) blendcolors(colors = colors[x], mode = network_blendmode))
      colors <- c(colors, setNames(node_colors[, 2], node_colors[, 1]))
      label_colors <- ifelse(colSums(col2rgb(colors)) > 255 * 2, "black", "white")
      df_nodes[["color"]] <- colors[df_nodes$name]
      df_nodes[["label_color"]] <- label_colors[df_nodes$name]
      df_nodes[["label"]] <- NA
      df_nodes[levels(df[["Description"]]), "label"] <- seq_len(nlevels(df[["Description"]]))

      draw_key_cust <- function(data, params, size) {
        data_text <- data
        data_text$label <- which(levels(df[["Description"]]) %in% names(colors)[colors == data_text$fill])
        data_text$colour <- "black"
        data_text$alpha <- 1
        data_text$size <- 11 / .pt
        grobTree(
          draw_key_point(data, list(color = "white", shape = 21)),
          ggrepel:::shadowtextGrob(label = data_text$label, bg.colour = "black", bg.r = 0.1, gp = gpar(col = "white", fontface = "bold"))
        )
      }

      p <- ggplot() +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, color = color), alpha = 1, lineend = "round", show.legend = FALSE) +
        geom_label(data = df_nodes[df_nodes$class == "gene", ], aes(x = dim1, y = dim2, label = name, fill = color, color = label_color), size = 3, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2), size = 8, color = "black", fill = "black", stroke = 1, shape = 21, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, fill = color), size = 7, color = "white", stroke = 1, shape = 21, key_glyph = draw_key_cust) +
        geom_text_repel(
          data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, label = label),
          fontface = "bold", min.segment.length = 0, segment.color = "black",
          point.size = NA, max.overlaps = 100, force = 0, color = "white", bg.color = "black", bg.r = 0.1, size = network_labelsize
        ) +
        scale_color_identity(guide = "none") +
        scale_fill_identity(
          name = "Term:", guide = "legend",
          labels = levels(df[["Description"]]),
          breaks = colors[levels(df[["Description"]])]
        ) +
        guides(color = guide_legend(override.aes = list(color = "transparent"))) +
        labs(x = "", y = "") +
        facet_grid(facet, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    }))
  } else if (plot_type == "enrichmap") {
    # enrichmap -------------------------------------------------------------------------------------------------
    plist <- suppressWarnings(lapply(df_list, function(df) {
      df_groups <- split(df, list(df$Database, df$Groups))
      df_groups <- lapply(df_groups, function(group) {
        filtered_group <- group[head(seq_len(nrow(group)), topTerm), , drop = FALSE]
        return(filtered_group)
      })
      df <- do.call(rbind, df_groups)

      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- factor(df[["Description"]], levels = unique(df[["Description"]]))
      df$geneID <- strsplit(df$geneID, "/")
      rownames(df) <- df[["ID"]]

      nodes <- df
      edges <- as.data.frame(t(combn(nodes$ID, 2)))
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- mapply(function(x, y) length(intersect(df[[x, "geneID"]], df[[y, "geneID"]])), edges$from, edges$to)
      edges <- edges[edges[["weight"]] > 0, , drop = FALSE]
      graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
      if (enrichmap_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(enrichmap_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", enrichmap_layout), list(graph))
      }
      clusters <- do.call(paste0("cluster_", enrichmap_cluster), list(graph))
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]
      df_nodes[["clusters"]] <- factor(paste0("C", clusters$membership), paste0("C", unique(sort(clusters$membership))))

      if (isTRUE(enrichmap_show_keyword)) {
        df_keyword1 <- df_nodes %>%
          mutate(keyword = strsplit(tolower(as.character(.data[["Description"]])), "\\s|\\n", perl = TRUE)) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Database, Groups, clusters) %>%
          reframe(
            keyword = capitalize(.data[["keyword"]]),
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
          group_by(Database, Groups, clusters) %>%
          arrange(desc(score)) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["keyword"]] <- str_wrap(df_keyword1[["keyword"]], width = character_width)
        df_keyword1[["label"]] <- paste0(df_keyword1[["clusters"]], ":\n", df_keyword1[["keyword"]])
      } else {
        if (enrichmap_label == "term") {
          df_nodes[["Description"]] <- str_wrap(df_nodes[["Description"]], width = character_width)
        }
        df_keyword1 <- df_nodes %>%
          group_by(Database, Groups, clusters) %>%
          arrange(desc(metric)) %>%
          reframe(keyword = Description) %>%
          distinct() %>%
          group_by(Database, Groups, clusters) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = "\n")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["label"]] <- paste0(df_keyword1[["clusters"]], ":\n", df_keyword1[["keyword"]])
      }

      df_keyword2 <- df_nodes %>%
        mutate(keyword = .data[["geneID"]]) %>%
        unnest(cols = "keyword") %>%
        group_by(.data[["keyword"]], Database, Groups, clusters) %>%
        reframe(
          keyword = .data[["keyword"]],
          score = sum(-(log10(.data[[metric]]))),
          count = n(),
          Database = .data[["Database"]],
          Groups = .data[["Groups"]],
          .groups = "keep"
        ) %>%
        distinct() %>%
        group_by(Database, Groups, clusters) %>%
        arrange(desc(score)) %>%
        slice_head(n = enrlichmap_nlabel) %>%
        reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
        as.data.frame()
      rownames(df_keyword2) <- as.character(df_keyword2[["clusters"]])
      df_keyword2[["keyword"]] <- str_wrap(df_keyword2[["keyword"]], width = character_width)
      df_keyword2[["label"]] <- paste0(df_keyword2[["clusters"]], ":\n", df_keyword2[["keyword"]])

      df_nodes[["keyword1"]] <- df_keyword1[as.character(df_nodes$clusters), "keyword"]
      df_nodes[["keyword2"]] <- df_keyword2[as.character(df_nodes$clusters), "keyword"]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      if (enrichmap_mark == "hull") {
        check_R("concaveman")
      }
      mark_layer <- do.call(
        switch(enrichmap_mark,
          "ellipse" = "geom_mark_ellipse",
          "hull" = "geom_mark_hull"
        ),
        list(
          data = df_nodes, aes(
            x = dim1, y = dim2, color = clusters, fill = clusters,
            label = clusters, description = if (enrichmap_label == "term") keyword1 else keyword2
          ),
          expand = unit(3, "mm"),
          alpha = 0.1,
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = enrichmap_labelsize * 2,
          label.fill = "grey95",
          label.minwidth = unit(character_width, "in"),
          label.buffer = unit(0, "mm"),
          con.size = 1,
          con.cap = 0
        )
      )

      p <- ggplot() +
        mark_layer +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, linewidth = weight), alpha = 0.1, lineend = "round") +
        geom_point(data = df_nodes, aes(x = dim1, y = dim2, size = Count, fill = clusters), color = "black", shape = 21) +
        labs(x = "", y = "") +
        scale_size(name = "Count", range = c(2, 6), scales::breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_linewidth(name = "Intersection", range = c(0.3, 3), scales::breaks_extended(n = 4)) +
        guides(linewidth = guide_legend(override.aes = list(alpha = 1, color = "grey"), order = 2)) +
        scale_fill_manual(
          name = switch(enrichmap_label,
            "term" = "Feature:",
            "feature" = "Term:"
          ),
          values = palette_scp(levels(df_nodes[["clusters"]]), palette = palette, palcolor = palcolor),
          labels = if (enrichmap_label == "term") df_keyword2[levels(df_nodes[["clusters"]]), "label"] else df_keyword1[levels(df_nodes[["clusters"]]), "label"],
          na.value = "grey80",
          aesthetics = c("colour", "fill")
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black", shape = NA), byrow = TRUE, order = 3)) +
        guides(color = guide_none()) +
        scale_x_continuous(expand = expansion(c(enrichmap_expand[1], enrichmap_expand[1]), 0)) +
        scale_y_continuous(expand = expansion(c(enrichmap_expand[2], enrichmap_expand[2]), 0)) +
        facet_grid(facet, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      return(p)
    }))
  } else if (plot_type == "wordcloud") {
    # wordcloud -------------------------------------------------------------------------------------------------
    check_R("ggwordcloud")
    check_R("jokergoo/simplifyEnrichment")
    plist <- lapply(df_list, function(df) {
      if (word_type == "term") {
        df_groups <- split(df, list(df$Database, df$Groups))
        df_groups <- df_groups[sapply(df_groups, nrow) > 0]
        for (i in seq_along(df_groups)) {
          df_sub <- df_groups[[i]]
          if (all(df_sub$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
            df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df_sub[["ID"]])
            if (nrow(df0 > 0)) {
              df_sub <- df0 %>%
                reframe(
                  keyword = .data[["keyword"]],
                  score = -(log10(.data[["padj"]])),
                  count = .data[["n_term"]],
                  Database = df_sub[["Database"]][1],
                  Groups = df_sub[["Groups"]][1]
                ) %>%
                filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
                filter(nchar(.data[["keyword"]]) >= 1) %>%
                filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) %>%
                distinct() %>%
                mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                as.data.frame()
              df_sub <- df_sub[head(order(df_sub[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
            } else {
              df_sub <- NULL
            }
          } else {
            df_sub <- df_sub %>%
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
            df_sub <- df_sub[head(order(df_sub[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
          }
          df_groups[[i]] <- df_sub
        }
        df <- do.call(rbind, df_groups)
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
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
        ) +
        scale_size(name = "Count", range = word_size, breaks = ceiling(seq(min(df[["count"]], na.rm = TRUE), max(df[["count"]], na.rm = TRUE), length.out = 3))) +
        guides(size = guide_legend(override.aes = list(colour = "black", label = "G"), order = 1)) +
        facet_grid(facet, scales = "free") +
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

#' @importFrom igraph degree neighbors
adjustlayout <- function(graph, layout, width, height = 2, scale = 100, iter = 100) {
  w <- width / 2
  layout[, 1] <- layout[, 1] / diff(range(layout[, 1])) * scale
  layout[, 2] <- layout[, 2] / diff(range(layout[, 2])) * scale

  adjusted <- c()
  # for (i in seq_len(iter)) {
  for (v in order(degree(graph), decreasing = TRUE)) {
    adjusted <- c(adjusted, v)
    neighbors <- as.numeric(neighbors(graph, V(graph)[v]))
    neighbors <- setdiff(neighbors, adjusted)
    x <- layout[v, 1]
    y <- layout[v, 2]
    r <- w[v]
    for (neighbor in neighbors) {
      nx <- layout[neighbor, 1]
      ny <- layout[neighbor, 2]
      ndist <- sqrt((nx - x)^2 + (ny - y)^2)
      nr <- w[neighbor]
      expect <- r + nr
      if (ndist < expect) {
        dx <- (x - nx) * (expect - ndist) / ndist
        dy <- (y - ny) * (expect - ndist) / ndist
        layout[neighbor, 1] <- nx - dx
        layout[neighbor, 2] <- ny - dy
        adjusted <- c(adjusted, neighbor)
      }
    }
  }
  # }

  for (i in seq_len(iter)) {
    dist_matrix <- as_matrix(dist(layout))
    nearest_neighbors <- apply(dist_matrix, 2, function(x) which(x == min(x[x > 0])), simplify = FALSE)
    # nearest_neighbors <- apply(dist_matrix, 2, function(x) {
    #   head(order(x), 3)[-1]
    # }, simplify = FALSE)
    for (v in sample(seq_len(nrow(layout)))) {
      neighbors <- unique(nearest_neighbors[[v]])
      x <- layout[v, 1]
      y <- layout[v, 2]
      r <- w[v]
      for (neighbor in neighbors) {
        nx <- layout[neighbor, 1]
        ny <- layout[neighbor, 2]
        nr <- w[neighbor]
        if (abs(nx - x) < (r + nr) && abs(ny - y) < height) {
          dx <- r + nr - (nx - x)
          dy <- height - (ny - y)
          if (sample(c(1, 0), 1) == 1) {
            dx <- 0
          } else {
            dy <- 0
          }
          layout[neighbor, 1] <- nx - dx
          layout[neighbor, 2] <- ny - dy
        }
      }
    }
  }
  return(layout)
}

#' GSEA Plot
#'
#' This function generates various types of plots for Gene Set Enrichment Analysis (GSEA) results.
#'
#' @inheritParams EnrichmentPlot
#' @param srt A Seurat object containing the results of RunDEtest and RunGSEA.
#' If specified, GSEA results will be extracted from the Seurat object automatically.
#' If not specified, the \code{res} arguments must be provided.
#' @param res Enrichment results generated by RunGSEA function. If provided, 'srt', 'test.use' and 'group_by' are ignored.
#' @param plot_type The type of plot to generate. Options are: "line", "comparison", "bar", "network", "enrichmap", "wordcloud". Default is "line".
#' @param direction The direction of enrichment to include in the plot. Must be one of "pos", "neg", or "both". The default value is "both".
#' @param line_width The linewidth for the line plot.
#' @param line_alpha The alpha value for the line plot.
#' @param line_color The color for the line plot.
#' @param n_coregene The number of core genes to label in the line plot.
#' @param sample_coregene Whether to randomly sample core genes for labeling in the line plot.
#' @param features_label A character vector of feature names to include as labels in the line plot.
#' @param label.fg The color of the labels.
#' @param label.bg The background color of the labels.
#' @param label.bg.r The radius of the rounding of the label's background.
#' @param label.size The size of the labels.
#'
#' @seealso \code{\link{RunGSEA}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
#' pancreas_sub <- RunGSEA(pancreas_sub, group_by = "CellType", db = "GO_BP", species = "Mus_musculus")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", id_use = "GO:0006412")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", id_use = c("GO:0046903", "GO:0015031", "GO:0007600")) %>%
#'   panel_fix_overall(height = 6) # As the plot is created by combining, we can adjust the overall height and width directly.
#'
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", direction = "neg")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", direction = "both")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", topTerm = 3, plot_type = "comparison", compare_only_sig = TRUE)
#'
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", plot_type = "bar")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", plot_type = "bar", direction = "both")
#' GSEAPlot(pancreas_sub,
#'   db = "GO_BP", group_by = "CellType", group_use = "Ductal",
#'   plot_type = "bar", topTerm = 20, direction = "both", palcolor = c("red3", "steelblue")
#' )
#'
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", plot_type = "network")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", plot_type = "enrichmap")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", plot_type = "wordcloud")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", plot_type = "wordcloud", word_type = "feature")
#'
#' @importFrom ggplot2 ggplot aes theme theme_classic alpha element_blank element_rect margin geom_line geom_point geom_rect geom_linerange geom_hline geom_vline geom_segment annotate ggtitle labs xlab ylab scale_x_continuous scale_y_continuous scale_color_manual scale_alpha_manual guides guide_legend guide_none
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRamp
#' @importFrom patchwork wrap_plots
#' @importFrom dplyr case_when filter pull %>%
#' @importFrom stats quantile
#' @importFrom gtable gtable_add_rows gtable_add_grob
#' @importFrom grid textGrob
#' @export
GSEAPlot <- function(srt, db = "GO_BP", group_by = NULL, test.use = "wilcox", res = NULL,
                     plot_type = c("line", "bar", "network", "enrichmap", "wordcloud", "comparison"),
                     group_use = NULL, id_use = NULL, pvalueCutoff = NULL, padjustCutoff = 0.05,
                     topTerm = ifelse(plot_type == "enrichmap", 100, 6), direction = c("pos", "neg", "both"), compare_only_sig = FALSE,
                     topWord = 100, word_type = c("term", "feature"), word_size = c(2, 8), words_excluded = NULL,
                     line_width = 1.5, line_alpha = 1, line_color = "#6BB82D",
                     n_coregene = 10, sample_coregene = FALSE, features_label = NULL,
                     label.fg = "black", label.bg = "white", label.bg.r = 0.1, label.size = 4,
                     network_layout = "fr", network_labelsize = 5, network_blendmode = "blend",
                     network_layoutadjust = TRUE, network_adjscale = 60, network_adjiter = 100,
                     enrichmap_layout = "fr", enrichmap_cluster = "fast_greedy", enrichmap_label = c("term", "feature"), enrichmap_labelsize = 5,
                     enrlichmap_nlabel = 4, enrichmap_show_keyword = FALSE, enrichmap_mark = c("ellipse", "hull"), enrichmap_expand = c(0.5, 0.5),
                     character_width = 50, lineheight = 0.5,
                     palette = "Spectral", palcolor = NULL,
                     aspect.ratio = NULL, legend.position = "right", legend.direction = "vertical",
                     theme_use = "theme_scp", theme_args = list(),
                     combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, seed = 11) {
  set.seed(seed)
  plot_type <- match.arg(plot_type)
  word_type <- match.arg(word_type)
  direction <- match.arg(direction)
  enrichmap_label <- match.arg(enrichmap_label)
  enrichmap_mark <- match.arg(enrichmap_mark)
  words_excluded <- words_excluded %||% SCP::words_excluded

  subplots <- 1:3
  rel_heights <- c(1.5, 0.5, 1)
  rel_width <- 3

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
  if (length(id_use) > 0) {
    topTerm <- Inf
    if (is.list(id_use)) {
      if (is.null(names(id_use))) {
        stop("'id_use' must be named when it is a list.")
      }
      if (!all(names(id_use) %in% enrichment[["Groups"]])) {
        stop(paste0("Names in 'id_use' is invalid: ", paste0(names(id_use)[!names(id_use) %in% enrichment[["Groups"]]], collapse = ",")))
      }
      enrichment_list <- list()
      for (i in seq_along(id_use)) {
        enrichment_list[[i]] <- enrichment[enrichment[["ID"]] %in% id_use[[i]] & enrichment[["Groups"]] %in% names(id_use)[i], , drop = FALSE]
      }
      enrichment <- do.call(rbind, enrichment_list)
    } else {
      enrichment <- enrichment[enrichment[["ID"]] %in% unlist(id_use), , drop = FALSE]
    }
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
    # comparison -------------------------------------------------------------------------------------------------
    if (length(id_use) > 0) {
      ids <- unlist(id_use)
    } else {
      ids <- NULL
      for (i in group_use) {
        df <- enrichment[enrichment[["Groups"]] == i, , drop = FALSE]
        df <- df[df[[metric]] < metric_value, , drop = FALSE]
        df <- df[order(df[[metric]]), , drop = FALSE]
        df_up <- df[df[["NES"]] > 0, , drop = FALSE]
        ID_up <- df_up[head(order(df_up[[metric]]), topTerm), "ID"]
        df_down <- df[df[["NES"]] < 0, , drop = FALSE]
        ID_down <- df_down[head(order(df_down[[metric]]), topTerm), "ID"]
        ids <- switch(direction,
          "pos" = unique(c(ids, head(ID_up, topTerm))),
          "neg" = unique(c(ids, head(ID_down, topTerm))),
          "both" = unique(c(ids, head(
            c(
              head(ID_up, ceiling(topTerm / 2)),
              head(ID_down, ceiling(topTerm / 2))
            ),
            topTerm
          )))
        )
      }
    }

    if (any(db %in% c("GO_sim", "GO_BP_sim", "GO_CC_sim", "GO_MF_sim"))) {
      enrichment_sub <- subset(enrichment_sim, ID %in% ids)
      enrichment_sub[["Database"]][enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")] <- paste0(enrichment_sub[["Database"]][enrichment_sub[["Database"]] %in% c("GO", "GO_BP", "GO_CC", "GO_MF")], "_sim")
    } else {
      enrichment_sub <- subset(enrichment, ID %in% ids)
    }
    enrichment_sub[["Database"]] <- factor(enrichment_sub[["Database"]], levels = db)
    enrichment_sub[["Description"]] <- capitalize(enrichment_sub[["Description"]])
    enrichment_sub[["Description"]] <- str_wrap(enrichment_sub[["Description"]], width = character_width)
    terms <- setNames(enrichment_sub[["Description"]], enrichment_sub[["ID"]])
    enrichment_sub[["Description"]] <- factor(enrichment_sub[["Description"]], levels = unique(rev(terms[ids])))
    enrichment_sub[["Significant"]] <- enrichment_sub[[metric]] < metric_value
    enrichment_sub[["Significant"]] <- factor(enrichment_sub[["Significant"]], levels = c("TRUE", "FALSE"))
    if (isTRUE(compare_only_sig)) {
      enrichment_sub <- enrichment_sub[enrichment_sub[["Significant"]] == "TRUE", , drop = FALSE]
    }
    enrichment_sub <- switch(direction,
      "pos" = enrichment_sub[enrichment_sub[["NES"]] > 0, , drop = FALSE],
      "neg" = enrichment_sub[enrichment_sub[["NES"]] < 0, , drop = FALSE],
      "both" = enrichment_sub
    )

    p <- ggplot(enrichment_sub, aes(x = Groups, y = Description)) +
      geom_point(aes(size = setSize, fill = NES, color = Significant), shape = 21, stroke = 0.8) +
      scale_size_area(name = "setSize", max_size = 6, n.breaks = 4) +
      guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 2)) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
      scale_fill_gradientn(
        name = "NES",
        n.breaks = 4,
        limits = c(-max(abs(enrichment_sub[["NES"]])), max(abs(enrichment_sub[["NES"]]))),
        colors = palette_scp(palette = palette, palcolor = palcolor),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0, order = 1)
      ) +
      scale_color_manual(
        name = paste0("Significant\n(", metric, "<", metric_value, ")", collapse = ""), values = c("TRUE" = "black", "FALSE" = "grey90"),
        guide = if (isTRUE(compare_only_sig)) guide_none() else guide_legend()
      ) +
      facet_grid(Database ~ ., scales = "free") +
      do.call(theme_use, theme_args) +
      theme(
        aspect.ratio = aspect.ratio,
        legend.position = legend.position,
        legend.direction = legend.direction,
        panel.grid.major = element_line(colour = "grey80", linetype = 2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(
          lineheight = lineheight, hjust = 1,
          face = ifelse(grepl("\n", levels(enrichment_sub[["Description"]])), "italic", "plain")
        )
      )
    plist <- list(p)
  } else if (plot_type == "line") {
    # line -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, , drop = FALSE]
        geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, , drop = FALSE]
        geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(res_enrich@result[["ID"]], id_use[[unique(res_enrich@result$Groups)]])
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 1) {
        gsdata <- gsInfo(object = res_enrich, id_use = geneSetID_use)
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
      gsdata[["NES"]] <- stat[gsdata$Description, "NES"]
      gsdata[[metric]] <- stat[gsdata$Description, metric]
      gsdata[["p.sig"]] <- stat[gsdata$Description, "p.sig"]
      gsdata[["DescriptionP"]] <- capitalize(gsdata[["Description"]])
      gsdata[["DescriptionP"]] <- str_wrap(gsdata[["DescriptionP"]], width = character_width)
      gsdata[["DescriptionP"]] <- paste0(gsdata[["DescriptionP"]], "\n(NES=", round(gsdata[["NES"]], 3), ", ", metric, "=", format(gsdata[[metric]], digits = 3, scientific = TRUE), ", ", gsdata[["p.sig"]], ")")
      gsdata[["DescriptionP"]] <- factor(gsdata[["DescriptionP"]], levels = unique(gsdata[["DescriptionP"]]))
      p <- ggplot(gsdata, aes(x = x)) +
        xlab(NULL) +
        theme_classic(base_size = 12) +
        theme(
          panel.grid.major = element_line(colour = "grey90", linetype = 2),
          panel.grid.minor = element_line(colour = "grey90", linetype = 2)
        ) +
        scale_x_continuous(expand = c(0.01, 0))
      es_layer <- geom_line(aes(y = runningScore, color = DescriptionP),
        linewidth = line_width, alpha = line_alpha
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
        theme_classic(base_size = 12) +
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
            features_label_tmp <- unlist(strsplit(gsdata$CoreGene[1], "/"))
            n_coregene <- min(n_coregene, length(features_label_tmp))
            if (isTRUE(sample_coregene)) {
              features_label_tmp <- sample(features_label_tmp, n_coregene, replace = FALSE)
            } else {
              features_label_tmp <- gsdata$GeneName[gsdata$GeneName %in% features_label_tmp][1:n_coregene]
            }
          } else {
            features_label_tmp <- features_label
          }
          df_gene <- gsdata[gsdata$position == 1 & gsdata$GeneName %in% features_label_tmp, , drop = FALSE]
          gene_drop <- features_label_tmp[!features_label_tmp %in% df_gene$GeneName]
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
          guides(color = guide_legend(title = "Term:", byrow = TRUE)) +
          do.call(theme_use, theme_args) +
          theme(
            legend.position = legend.position,
            legend.direction = legend.direction
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
          plotlist[[i]] <- panel_fix_overall(plotlist[[i]], height = rel_heights[i], units = "null", margin = 0, respect = TRUE, return_grob = TRUE)
          plotlist[[i]] <- panel_fix_overall(plotlist[[i]], width = rel_width, units = "null", margin = 0, respect = TRUE, return_grob = TRUE)
        }
        p_out <- do.call(rbind, c(plotlist, size = "first"))

        if (length(geneSetID_use) > 1) {
          p_out <- add_grob(p_out, legend, legend.position)
        }
        lab <- textGrob(label = nm, rot = -90, hjust = 0.5)
        p_out <- add_grob(p_out, lab, "right", clip = "off")
        p_out <- wrap_plots(p_out)
        plist[[nm]] <- p_out
      }
    }
  } else if (plot_type == "bar") {
    # bar -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, , drop = FALSE]
        geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, , drop = FALSE]
        geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(res_enrich@result[["ID"]], id_use[[unique(res_enrich@result$Groups)]])
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      stat <- res_enrich[geneSetID_use, , drop = FALSE]
      stat <- stat[order(stat[["NES"]]), , drop = FALSE]
      rownames(stat) <- stat[, "Description"]
      stat[["Description"]] <- capitalize(stat[["Description"]])
      stat[["Description"]] <- str_wrap(stat[["Description"]], width = character_width)
      stat[["Description"]] <- factor(stat[["Description"]], levels = unique(stat[["Description"]]))
      stat[["Direction"]] <- ifelse(stat[["NES"]] > 0, "Pos", "Neg")
      stat[["Direction"]] <- factor(stat[["Direction"]], levels = c("Pos", "Neg"))

      p <- ggplot(stat, aes(
        x = .data[["NES"]], y = .data[["Description"]]
      )) +
        geom_vline(xintercept = 0) +
        geom_col(aes(fill = .data[["Direction"]], alpha = -log10(.data[[metric]])), color = "black") +
        geom_text(
          aes(
            x = 0, y = .data[["Description"]], label = .data[["Description"]],
            hjust = ifelse(.data[["NES"]] > 0, 1, 0),
          ),
          nudge_x = ifelse(stat[["NES"]] > 0, -0.05, 0.05),
          lineheight = lineheight,
          fontface = ifelse(grepl("\n", levels(stat[["Description"]])), "italic", "plain")
        ) +
        scale_fill_manual(
          values = palette_scp(x = rev(levels(stat[["Direction"]])), palette = palette, palcolor = rev(palcolor)),
          guide = if (direction == "both") guide_legend(order = 1) else guide_none()
        ) +
        facet_grid(Database ~ Groups, scales = "free") +
        coord_cartesian(xlim = c(-max(abs(stat[["NES"]])), max(abs(stat[["NES"]])))) +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "network") {
    # network -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, , drop = FALSE]
        geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, , drop = FALSE]
        geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(res_enrich@result[["ID"]], id_use[[unique(res_enrich@result$Groups)]])
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]
      df$p.sig <- case_when(
        df[[metric]] > 0.05 ~ "ns  ",
        df[[metric]] <= 0.05 & df[[metric]] > 0.01 ~ "*   ",
        df[[metric]] <= 0.01 & df[[metric]] > 0.001 ~ "**  ",
        df[[metric]] <= 0.001 & df[[metric]] > 0.0001 ~ "*** ",
        df[[metric]] <= 0.0001 ~ "****"
      )
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- paste0(df[["Description"]], "\n(NES=", round(df[["NES"]], 3), ", ", metric, "=", format(df[[metric]], digits = 3, scientific = TRUE), ", ", df[["p.sig"]], ")")
      df[["Description"]] <- factor(df[["Description"]], levels = unique(df[["Description"]]))
      df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
      df_unnest <- unnest(df, cols = "geneID")

      nodes <- rbind(
        data.frame("ID" = df[["Description"]], class = "term", metric = df[["metric"]]),
        data.frame("ID" = unique(df_unnest$geneID), class = "gene", metric = 0)
      )
      nodes$Database <- df$Database[1]
      nodes$Groups <- df$Groups[1]
      edges <- as.data.frame(df_unnest[, c("Description", "geneID")])
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- 1
      graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
      if (network_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(network_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", network_layout), list(graph))
      }
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      if (isTRUE(network_layoutadjust)) {
        width <- nchar(df_nodes$name)
        width[df_nodes$class == "term"] <- 8
        layout <- adjustlayout(
          graph = graph, layout = layout, width = width, height = 2,
          scale = network_adjscale, iter = network_adjiter
        )
      }
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      colors <- palette_scp(levels(df[["Description"]]), palette = palette, palcolor = palcolor)
      df_edges[["color"]] <- colors[df_edges$from]
      node_colors <- aggregate(df_unnest$Description, by = list(df_unnest$geneID), FUN = function(x) blendcolors(colors = colors[x], mode = network_blendmode))
      colors <- c(colors, setNames(node_colors[, 2], node_colors[, 1]))
      label_colors <- ifelse(colSums(col2rgb(colors)) > 255 * 2, "black", "white")
      df_nodes[["color"]] <- colors[df_nodes$name]
      df_nodes[["label_color"]] <- label_colors[df_nodes$name]
      df_nodes[["label"]] <- NA
      df_nodes[levels(df[["Description"]]), "label"] <- seq_len(nlevels(df[["Description"]]))

      draw_key_cust <- function(data, params, size) {
        data_text <- data
        data_text$label <- which(levels(df[["Description"]]) %in% names(colors)[colors == data_text$fill])
        data_text$colour <- "black"
        data_text$alpha <- 1
        data_text$size <- 11 / .pt
        grobTree(
          draw_key_point(data, list(color = "white", shape = 21)),
          ggrepel:::shadowtextGrob(label = data_text$label, bg.colour = "black", bg.r = 0.1, gp = gpar(col = "white", fontface = "bold"))
        )
      }

      p <- ggplot() +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, color = color), alpha = 1, lineend = "round", show.legend = FALSE) +
        geom_label(data = df_nodes[df_nodes$class == "gene", ], aes(x = dim1, y = dim2, label = name, fill = color, color = label_color), size = 3, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2), size = 8, color = "black", fill = "black", stroke = 1, shape = 21, show.legend = FALSE) +
        geom_point(data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, fill = color), size = 7, color = "white", stroke = 1, shape = 21, key_glyph = draw_key_cust) +
        geom_text_repel(
          data = df_nodes[df_nodes$class == "term", ], aes(x = dim1, y = dim2, label = label),
          fontface = "bold", min.segment.length = 0, segment.color = "black",
          point.size = NA, max.overlaps = 100, force = 0, color = "white", bg.color = "black", bg.r = 0.1, size = network_labelsize
        ) +
        scale_color_identity(guide = "none") +
        scale_fill_identity(
          name = "Term:", guide = "legend",
          labels = levels(df[["Description"]]),
          breaks = colors[levels(df[["Description"]])]
        ) +
        guides(fill = guide_legend(title = "Term:", byrow = TRUE)) +
        labs(x = "", y = "") +
        facet_grid(Database ~ Groups, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "enrichmap") {
    # enrichmap -------------------------------------------------------------------------------------------------
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, , drop = FALSE]
        geneSetID_up <- geneSetID_up[head(order(geneSetID_up[[metric]]), topTerm), "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, , drop = FALSE]
        geneSetID_down <- geneSetID_down[head(order(geneSetID_down[[metric]]), topTerm), "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(head(geneSetID_up, topTerm)),
          "neg" = unique(head(geneSetID_down, topTerm)),
          "both" = unique(head(
            c(
              head(geneSetID_up, ceiling(topTerm / 2)),
              head(geneSetID_down, ceiling(topTerm / 2))
            ),
            topTerm
          ))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(res_enrich@result[["ID"]], id_use[[unique(res_enrich@result$Groups)]])
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]
      df[["metric"]] <- -log10(df[[metric]])
      df[["Description"]] <- capitalize(df[["Description"]])
      df[["Description"]] <- str_wrap(df[["Description"]], width = character_width)
      df[["Description"]] <- factor(df[["Description"]], levels = unique(df[["Description"]]))
      df[["Direction"]] <- ifelse(df[["NES"]] > 0, "Pos", "Neg")
      df[["Direction"]] <- factor(df[["Direction"]], levels = c("Pos", "Neg"))
      df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
      df[["Count"]] <- sapply(df[["geneID"]], length)
      rownames(df) <- df[["ID"]]

      nodes <- df
      edges <- as.data.frame(t(combn(nodes$ID, 2)))
      colnames(edges) <- c("from", "to")
      edges[["weight"]] <- mapply(function(x, y) length(intersect(df[[x, "geneID"]], df[[y, "geneID"]])), edges$from, edges$to)
      edges <- edges[edges[["weight"]] > 0, , drop = FALSE]
      graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
      if (enrichmap_layout %in% c("circle", "tree", "grid")) {
        layout <- switch(enrichmap_layout,
          "circle" = layout_in_circle(graph),
          "tree" = layout_as_tree(graph),
          "grid" = layout_on_grid(graph)
        )
      } else {
        layout <- do.call(paste0("layout_with_", enrichmap_layout), list(graph))
      }
      clusters <- do.call(paste0("cluster_", enrichmap_cluster), list(graph))
      df_graph <- as_data_frame(graph, what = "both")

      df_nodes <- df_graph$vertices
      df_nodes[["dim1"]] <- layout[, 1]
      df_nodes[["dim2"]] <- layout[, 2]
      df_nodes[["clusters"]] <- factor(paste0("C", clusters$membership), paste0("C", unique(sort(clusters$membership))))

      if (isTRUE(enrichmap_show_keyword)) {
        df_keyword1 <- df_nodes %>%
          mutate(keyword = strsplit(tolower(as.character(.data[["Description"]])), "\\s|\\n", perl = TRUE)) %>%
          unnest(cols = "keyword") %>%
          group_by(.data[["keyword"]], Database, Groups, clusters) %>%
          reframe(
            keyword = capitalize(.data[["keyword"]]),
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
          group_by(Database, Groups, clusters) %>%
          arrange(desc(score)) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["keyword"]] <- str_wrap(df_keyword1[["keyword"]], width = character_width)
        df_keyword1[["label"]] <- paste0(df_keyword1[["clusters"]], ":\n", df_keyword1[["keyword"]])
      } else {
        if (enrichmap_label == "term") {
          df_nodes[["Description"]] <- str_wrap(df_nodes[["Description"]], width = character_width)
        }
        df_keyword1 <- df_nodes %>%
          group_by(Database, Groups, clusters) %>%
          arrange(desc(metric)) %>%
          reframe(keyword = Description) %>%
          distinct() %>%
          group_by(Database, Groups, clusters) %>%
          slice_head(n = enrlichmap_nlabel) %>%
          reframe(keyword = paste0(.data[["keyword"]], collapse = "\n")) %>%
          as.data.frame()
        rownames(df_keyword1) <- as.character(df_keyword1[["clusters"]])
        df_keyword1[["label"]] <- paste0(df_keyword1[["clusters"]], ":\n", df_keyword1[["keyword"]])
      }

      df_keyword2 <- df_nodes %>%
        mutate(keyword = .data[["geneID"]]) %>%
        unnest(cols = "keyword") %>%
        group_by(.data[["keyword"]], Database, Groups, clusters) %>%
        reframe(
          keyword = .data[["keyword"]],
          score = sum(-(log10(.data[[metric]]))),
          count = n(),
          Database = .data[["Database"]],
          Groups = .data[["Groups"]],
          .groups = "keep"
        ) %>%
        distinct() %>%
        group_by(Database, Groups, clusters) %>%
        arrange(desc(score)) %>%
        slice_head(n = enrlichmap_nlabel) %>%
        reframe(keyword = paste0(.data[["keyword"]], collapse = " ")) %>%
        as.data.frame()
      rownames(df_keyword2) <- as.character(df_keyword2[["clusters"]])
      df_keyword2[["keyword"]] <- str_wrap(df_keyword2[["keyword"]], width = character_width)
      df_keyword2[["label"]] <- paste0(df_keyword2[["clusters"]], ":\n", df_keyword2[["keyword"]])

      df_nodes[["keyword1"]] <- df_keyword1[as.character(df_nodes$clusters), "keyword"]
      df_nodes[["keyword2"]] <- df_keyword2[as.character(df_nodes$clusters), "keyword"]

      df_edges <- df_graph$edges
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]] <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]] <- df_nodes[df_edges[["to"]], "dim2"]

      if (enrichmap_mark == "hull") {
        check_R("concaveman")
      }
      mark_layer <- do.call(
        switch(enrichmap_mark,
          "ellipse" = "geom_mark_ellipse",
          "hull" = "geom_mark_hull"
        ),
        list(
          data = df_nodes, aes(
            x = dim1, y = dim2, color = clusters, fill = clusters,
            label = clusters, description = if (enrichmap_label == "term") keyword1 else keyword2
          ),
          expand = unit(3, "mm"),
          alpha = 0.1,
          label.margin = margin(1, 1, 1, 1, "mm"),
          label.fontsize = enrichmap_labelsize * 2,
          label.fill = "grey95",
          label.minwidth = unit(character_width, "in"),
          label.buffer = unit(0, "mm"),
          con.size = 1,
          con.cap = 0
        )
      )

      p <- ggplot() +
        mark_layer +
        geom_segment(data = df_edges, aes(x = from_dim1, y = from_dim2, xend = to_dim1, yend = to_dim2, linewidth = weight), alpha = 0.1, lineend = "round") +
        geom_point(data = df_nodes, aes(x = dim1, y = dim2, size = Count, fill = clusters), color = "black", shape = 21) +
        labs(x = "", y = "") +
        scale_size(name = "Count", range = c(2, 6), scales::breaks_extended(n = 4)) +
        guides(size = guide_legend(override.aes = list(fill = "grey30", shape = 21), order = 1)) +
        scale_linewidth(name = "Intersection", range = c(0.3, 3), scales::breaks_extended(n = 4)) +
        guides(linewidth = guide_legend(override.aes = list(alpha = 1, color = "grey"), order = 2)) +
        scale_fill_manual(
          name = switch(enrichmap_label,
            "term" = "Feature:",
            "feature" = "Term:"
          ),
          values = palette_scp(levels(df_nodes[["clusters"]]), palette = palette, palcolor = palcolor),
          labels = if (enrichmap_label == "term") df_keyword2[levels(df_nodes[["clusters"]]), "label"] else df_keyword1[levels(df_nodes[["clusters"]]), "label"],
          na.value = "grey80",
          aesthetics = c("colour", "fill")
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black", shape = NA), byrow = TRUE, order = 3)) +
        guides(color = guide_none()) +
        scale_x_continuous(expand = expansion(c(enrichmap_expand[1], enrichmap_expand[1]), 0)) +
        scale_y_continuous(expand = expansion(c(enrichmap_expand[2], enrichmap_expand[2]), 0)) +
        facet_grid(Database ~ Groups, scales = "free") +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )
      plist[[nm]] <- p
    }
  } else if (plot_type == "wordcloud") {
    # wordcloud -------------------------------------------------------------------------------------------------
    check_R("ggwordcloud")
    check_R("jokergoo/simplifyEnrichment")
    for (nm in names(res)) {
      res_enrich <- res[[nm]]
      if (is.null(id_use)) {
        geneSetID_filter <- res_enrich@result[res_enrich@result[[metric]] < metric_value, , drop = FALSE]
        geneSetID_filter <- geneSetID_filter[order(geneSetID_filter[[metric]]), , drop = FALSE]
        geneSetID_up <- geneSetID_filter[geneSetID_filter[["NES"]] > 0, "ID"]
        geneSetID_down <- geneSetID_filter[geneSetID_filter[["NES"]] < 0, "ID"]
        geneSetID_use <- switch(direction,
          "pos" = unique(geneSetID_up),
          "neg" = unique(geneSetID_down),
          "both" = unique(c(geneSetID_up, geneSetID_down))
        )
      } else {
        if (is.list(id_use)) {
          geneSetID_use <- intersect(res_enrich@result[["ID"]], id_use[[unique(res_enrich@result$Groups)]])
        } else {
          geneSetID_use <- id_use
        }
      }
      if (length(geneSetID_use) == 0) {
        plist[[nm]] <- NULL
        next
      }
      df <- res_enrich[geneSetID_use, , drop = FALSE]

      if (word_type == "term") {
        df_groups <- split(df, list(df$Database, df$Groups))
        df_groups <- df_groups[sapply(df_groups, nrow) > 0]
        for (i in seq_along(df_groups)) {
          df_sub <- df_groups[[i]]
          if (all(df_sub$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
            df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df_sub[["ID"]])
            if (nrow(df0 > 0)) {
              df_sub <- df0 %>%
                reframe(
                  keyword = .data[["keyword"]],
                  score = -(log10(.data[["padj"]])),
                  count = .data[["n_term"]],
                  Database = df_sub[["Database"]][1],
                  Groups = df_sub[["Groups"]][1]
                ) %>%
                filter(!grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])) %>%
                filter(nchar(.data[["keyword"]]) >= 1) %>%
                filter(!tolower(.data[["keyword"]]) %in% tolower(words_excluded)) %>%
                distinct() %>%
                mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
                as.data.frame()
              df_sub <- df_sub[head(order(df_sub[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
            } else {
              df_sub <- NULL
            }
          } else {
            df_sub <- df_sub %>%
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
            df_sub <- df_sub[head(order(df_sub[["score"]], decreasing = TRUE), topWord), , drop = FALSE]
          }
          df_groups[[i]] <- df_sub
        }
        df <- do.call(rbind, df_groups)
      } else {
        df <- df %>%
          mutate(keyword = strsplit(as.character(.data[["core_enrichment"]]), "/")) %>%
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
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
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
      plist[[nm]] <- p
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

gsInfo <- function(object, id_use) {
  geneList <- object@geneList
  if (is.numeric(id_use)) {
    id_use <- object@result[id_use, "ID"]
  }
  geneSet <- object@geneSets[[id_use]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[id_use, "Description"]
  df$CoreGene <- object@result[id_use, "core_enrichment"]
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

#' @importFrom grid grobTree
#' @importFrom ggplot2 ggplotGrob
as_grob <- function(plot, ...) {
  if (inherits(plot, "gList")) {
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
  if (position == "none" || is.null(grob)) {
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
