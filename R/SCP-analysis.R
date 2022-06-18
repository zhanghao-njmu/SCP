#' Gene ID conversion function using biomart
#'
#' This function can convert different gene ID types within one species or bewteen two species using biomart service.
#'
#' @param geneID A vector of the geneID character.
#' @param geneID_from_IDtype ID type of the geneID. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param geneID_to_IDtype ID type to convert to. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param species_from Latin names for animals of the input geneID.
#' @param species_to Latin names for animals of the output geneID. e.g. "Homo_sapiens","Mus_musculus"
#' @param Ensembl_version Ensembl database version. If NULL, use the current release version.
#' @param try_times Times when trying to connect with the biomart service.
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'.
#'
#' @return A list.
#'
#' @examples
#' if (interactive()) {
#'   res <- GeneConvert(
#'     geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"), geneID_from_IDtype = "symbol", geneID_to_IDtype = "entrez_id",
#'     species_from = "Homo_sapiens", species_to = "Mus_musculus", Ensembl_version = 103
#'   )
#'   str(res)
#' }
#' @importFrom httr set_config config
#' @importFrom dplyr "%>%" group_by mutate .data
#' @importFrom tidyr unnest
#' @importFrom tidyselect all_of
#' @importFrom reshape2 dcast melt
#' @importFrom R.cache loadCache saveCache
#' @importFrom biomaRt listEnsemblArchives useMart listDatasets useDataset getBM listAttributes
#' @export
#'
GeneConvert <- function(geneID, geneID_from_IDtype = "symbol", geneID_to_IDtype = "entrez_id",
                        species_from = "Homo_sapiens", species_to = NULL,
                        Ensembl_version = 103, try_times = 5, mirror = NULL) {
  httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))

  if (missing(geneID)) {
    stop("geneID must be provided.")
  }
  if (is.null(species_to)) {
    species_to <- species_from
  }
  if (is.null(Ensembl_version)) {
    Ensembl_version <- "current_release"
  }
  species_from_split <- unlist(strsplit(species_from, split = "_"))
  species_to_split <- unlist(strsplit(species_to, split = "_"))
  species_from_simp <- paste0(tolower(substring(species_from_split[1], 1, 1)), species_from_split[2])
  species_to_simp <- paste0(tolower(substring(species_to_split[1], 1, 1)), species_to_split[2])
  geneID_from_IDtype <- sapply(geneID_from_IDtype, tolower)
  geneID_to_IDtype <- sapply(unique(geneID_to_IDtype), tolower)

  if ("symbol" %in% geneID_from_IDtype) {
    geneID_from_IDtype <- geneID_from_IDtype[geneID_from_IDtype != "symbol"]
    geneID_from_IDtype <- c("ensembl_symbol", "entrez_symbol", "uniprot_symbol", "wiki_symbol", geneID_from_IDtype)
  }
  from_IDtype <- sapply(geneID_from_IDtype, function(x) {
    switch(x,
      "ensembl_symbol" = "external_gene_name",
      "ensembl_id" = "ensembl_gene_id",
      "entrez_symbol" = "entrezgene_accession",
      "entrez_id" = "entrezgene_id",
      "uniprot_symbol" = "uniprot_gn_symbol",
      "wiki_symbol" = "wikigene_name"
    )
  })
  names(from_IDtype) <- geneID_from_IDtype

  to_IDtype <- sapply(geneID_to_IDtype, function(x) {
    switch(x,
      "symbol" = "external_gene_name",
      "ensembl_symbol" = "external_gene_name",
      "entrez_symbol" = "external_gene_name",
      "ensembl_id" = "ensembl_gene_id",
      "entrez_id" = "entrezgene_id"
    )
  })

  if (species_from != species_to & all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    to_IDtype <- sapply(to_IDtype, function(x) {
      switch(x,
        "external_gene_name" = "associated_gene_name",
        "ensembl_gene_id" = "ensembl_gene"
      )
    })
    to_attr <- paste(species_to_simp, to_IDtype, sep = "_homolog_")
    names(to_attr) <- geneID_to_IDtype
  } else {
    to_attr <- to_IDtype
    names(to_attr) <- geneID_to_IDtype
  }

  message("Connecting to the Ensembl archives...")
  archives <- NULL
  ntry <- 0
  while (is.null(archives)) {
    ntry <- ntry + 1
    archives <- tryCatch(expr = {
      listEnsemblArchives()
    }, error = function(e) {
      message(e)
      message("Get errors when connecting with EnsemblArchives...\nRetrying...")
      Sys.sleep(1)
      return(NULL)
    })
    if (is.null(archives) & ntry >= try_times) {
      stop("Stop connecting...")
    }
  }

  Ensembl_version <- as.character(Ensembl_version)
  if (Ensembl_version == "current_release") {
    url <- archives[which(archives$current_release == "*"), "url"]
    version <- as.character(archives[which(archives$current_release == "*"), "version"])
    message("Using the ", Ensembl_version, "(", version, ")", " version of biomart...")
  } else if (Ensembl_version %in% archives$version) {
    url <- archives[which(archives$version == Ensembl_version), "url"]
    version <- as.character(archives[which(archives$version == Ensembl_version), "version"])
    message("Using the ", version, " version of biomart...")
  } else {
    stop("Ensembl_version is invalid. Must be one of current_release,", paste0(archives$version, collapse = ","))
  }

  message("Connecting to the biomart...")
  mart <- NULL
  ntry <- 0

  while (is.null(mart)) {
    ntry <- ntry + 1
    if (!is.null(mirror)) {
      mart <- tryCatch(
        expr = {
          useEnsembl("ensembl", mirror = mirror)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with mirror...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
    }
    if (is.null(mart)) {
      mart <- tryCatch(
        expr = {
          useMart("ensembl", host = url)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with ensembl mart...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
    }
    if (is.null(mart) & ntry >= try_times) {
      stop("Stop connecting...")
    }
  }
  Datasets <- listDatasets(mart)

  dataset <- paste0(species_from_simp, "_gene_ensembl")
  if (!dataset %in% Datasets$dataset) {
    warning(paste0("Can not find the dataset for the species: ", species_from, " (", dataset, ")"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = NULL))
  }
  message("Connecting to the dataset ", dataset, " ...")
  mart1 <- NULL
  ntry <- 0
  while (is.null(mart1)) {
    ntry <- ntry + 1
    mart1 <- tryCatch(
      expr = {
        useDataset(dataset = dataset, mart = mart)
      },
      error = function(e) {
        message(e)
        message("Get errors when connecting with ensembl mart...\nRetrying...")
        Sys.sleep(1)
        return(NULL)
      }
    )
    if (is.null(mart1) & ntry >= try_times) {
      stop("Stop connecting...")
    }
  }

  if (species_from != species_to & any(!geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    dataset2 <- paste0(species_to_simp, "_gene_ensembl")
    if (!dataset2 %in% Datasets$dataset) {
      warning(paste0("Can not find the dataset for the species: ", species_from, " (", dataset2, ")"), immediate. = TRUE)
      return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = NULL))
    }
    message("Connecting to the dataset ", dataset2, " ...")
    mart2 <- NULL
    ntry <- 0
    while (is.null(mart2)) {
      ntry <- ntry + 1
      mart2 <- tryCatch(
        expr = {
          useDataset(dataset = dataset2, mart = mart)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with ensembl mart...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
      if (is.null(mart2) & ntry >= try_times) {
        stop("Stop connecting...")
      }
    }
  }

  Attributes <- listAttributes(mart1)
  from_IDtype <- from_IDtype[from_IDtype %in% Attributes[, "name"]]
  geneID_res_list <- list()
  total <- length(geneID)
  geneID_res <- NULL
  if (any(!to_attr %in% Attributes$name)) {
    to_attr_drop <- to_attr[!to_attr %in% Attributes$name]
    to_attr <- to_attr[to_attr %in% Attributes$name]
    warning(paste0("Can not find the attributes for the species ", species_from, ": ", paste(to_attr_drop, collapse = ", ")), immediate. = TRUE)
    if (length(to_attr) == 0) {
      warning("No attribute found for the species ", species_from, ". Please check the 'Attributes' in the result.", immediate. = TRUE)
      return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = Attributes))
    }
  }

  message("Match and convert the geneID...\n")
  if (species_from != species_to) {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- getBM(
          mart = mart1,
          attributes = c(from_attr, "ensembl_gene_id"),
          filters = from_attr,
          values = list(geneID)
        )
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id_tmp")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1[, "from_IDtype"] <- from_name

        if (all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
          geneID_res2 <- getBM(
            mart = mart1,
            attributes = unique(c("ensembl_gene_id", to_attr)),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id_tmp"])
          )
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          geneID_res2[, "ensembl_gene_id_tmp"] <- geneID_res2[, "ensembl_gene_id"]
          geneID_res2 <- geneID_res2[, c("ensembl_gene_id_tmp", to_attr)]
          geneID_res2 <- unique(reshape2::melt(geneID_res2, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
          geneID_res2$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res2$to_IDtype]
          geneID_res_merge <- merge(x = geneID_res1, y = geneID_res2, by = "ensembl_gene_id_tmp")
        } else {
          homolog_ensembl_gene <- paste(species_to_simp, "ensembl_gene", sep = "_homolog_")
          geneID_res2 <- getBM(
            mart = mart1,
            attributes = c("ensembl_gene_id", homolog_ensembl_gene),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id_tmp"])
          )
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          colnames(geneID_res2) <- c("ensembl_gene_id_tmp", homolog_ensembl_gene)

          geneID_res3 <- getBM(
            mart = mart2,
            attributes = unique(c("ensembl_gene_id", to_attr)),
            filters = "ensembl_gene_id",
            values = list(geneID_res2[, homolog_ensembl_gene])
          )
          geneID_res3 <- geneID_res3[geneID_res3[, "ensembl_gene_id"] %in% geneID_res2[, homolog_ensembl_gene], ]
          if (nrow(geneID_res3) == 0) {
            next
          }
          geneID_res3[, "ensembl_gene_id_tmp2"] <- geneID_res3[, "ensembl_gene_id"]
          geneID_res3 <- geneID_res3[, c("ensembl_gene_id_tmp2", to_attr)]
          geneID_res3 <- unique(reshape2::melt(geneID_res3, id.vars = "ensembl_gene_id_tmp2", variable.name = "to_IDtype", value.name = "to_geneID"))
          colnames(geneID_res3)[1] <- homolog_ensembl_gene
          geneID_res3$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res3$to_IDtype]
          geneID_res_merge <- Reduce(merge, list(geneID_res1, geneID_res2, geneID_res3))
        }
        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c("from_IDtype", "from_geneID", "to_IDtype", "to_geneID")]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  } else {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- getBM(
          mart = mart1,
          attributes = unique(c("ensembl_gene_id", from_attr, to_attr)),
          filters = from_attr,
          values = list(geneID)
        )
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        geneID_res1[, "ensembl_gene_id_tmp"] <- geneID_res1[, "ensembl_gene_id"]
        geneID_res1_from <- unique(geneID_res1[, c("ensembl_gene_id_tmp", from_attr)])
        colnames(geneID_res1_from) <- c("ensembl_gene_id_tmp", "from_geneID")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1_from[, "from_IDtype"] <- from_name

        geneID_res1_to <- unique(geneID_res1[, c("ensembl_gene_id_tmp", to_attr)])
        geneID_res1_to <- unique(reshape2::melt(geneID_res1_to, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
        geneID_res1_to$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res1_to$to_IDtype]

        geneID_res_merge <- merge(x = geneID_res1_from, y = geneID_res1_to, by = "ensembl_gene_id_tmp")

        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c("from_IDtype", "from_geneID", "to_IDtype", "to_geneID")]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  }
  message(
    paste0(
      paste0(rep("=", 30), collapse = ""), "\n",
      total - length(geneID), " genes mapped\n",
      length(geneID), " genes unmapped"
    ), "\n",
    paste0(rep("=", 30), collapse = ""), "\n"
  )
  geneID_res <- bind_rows(geneID_res_list) %>% unique()
  if (is.null(geneID_res) || nrow(geneID_res) == 0) {
    warning(paste0("No gene mapped"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Datasets = Datasets, Attributes = Attributes))
  }
  geneID_collapse <- geneID_res %>%
    group_by(.data[["from_geneID"]], .data[["to_IDtype"]]) %>%
    mutate(
      from_geneID = unique(.data[["from_geneID"]]),
      to_IDtype = unique(.data[["to_IDtype"]]),
      to_geneID = list(unique(.data[["to_geneID"]][!.data[["to_geneID"]] %in% c("", NA)]))
    )
  geneID_collapse <- unique(as.data.frame(geneID_collapse[, c("from_geneID", "to_IDtype", "to_geneID")]))
  geneID_collapse <- geneID_collapse[sapply(geneID_collapse$to_geneID, length) > 0, ]
  geneID_collapse <- reshape2::dcast(geneID_collapse, formula = from_geneID ~ to_IDtype, value.var = "to_geneID")
  rownames(geneID_collapse) <- geneID_collapse[, "from_geneID"]

  geneID_expand <- geneID_collapse
  for (i in colnames(geneID_expand)[sapply(geneID_expand, class) == "list"]) {
    geneID_expand <- unnest(data = geneID_expand, cols = all_of(i))
  }
  geneID_expand <- as.data.frame(geneID_expand)

  return(list(geneID_res = geneID_res, geneID_collapse = geneID_collapse, geneID_expand = geneID_expand, Ensembl_version = version, Datasets = Datasets, Attributes = Attributes, geneID_unmapped = geneID))
}

#' Prefetch cycle gene
#'
#' Based on the human cell cycle genes, the cell cycle genes of the corresponding species were captured by homologous gene conversion.
#'
#' @inheritParams GeneConvert
#' @param species Latin names for animals,i.e., "Homo_sapiens", "Mus_musculus"
#' @param use_cached_gene Whether to use previously cached cell cycle gene conversion results for the species.
#'
#' @return A list of S-phase and G2M-phase genes.
#'
#' @examples
#' if (interactive()) {
#'   ccgenes <- CC_GenePrefetch("Mus_musculus")
#'   names(ccgenes)
#' }
#' @importFrom R.cache loadCache saveCache
#' @export
CC_GenePrefetch <- function(species, Ensembl_version = 103, mirror = NULL, try_times = 5, use_cached_gene = TRUE) {
  cc_S_genes <- Seurat::cc.genes.updated.2019$s.genes
  cc_G2M_genes <- Seurat::cc.genes.updated.2019$g2m.genes
  res <- NULL
  if (species != "Homo_sapiens") {
    if (isTRUE(use_cached_gene)) {
      res <- loadCache(key = list(species))
    }
    if (is.null(res)) {
      res <- GeneConvert(
        geneID = unique(c(cc_S_genes, cc_G2M_genes)),
        geneID_from_IDtype = "symbol",
        geneID_to_IDtype = "symbol",
        species_from = "Homo_sapiens",
        species_to = species,
        Ensembl_version = Ensembl_version,
        try_times = try_times,
        mirror = mirror
      )
      saveCache(res, key = list(species))
    } else {
      message("Using cached conversion results for ", species)
    }
    genes <- res[["geneID_collapse"]]
    cc_S_genes <- unlist(genes[cc_S_genes[cc_S_genes %in% rownames(genes)], "symbol"])
    cc_G2M_genes <- unlist(genes[cc_G2M_genes[cc_G2M_genes %in% rownames(genes)], "symbol"])
  }
  return(list(
    res = res,
    cc_S_genes = cc_S_genes,
    cc_G2M_genes = cc_G2M_genes
  ))
}

#' Cell classification
#'
#' @examples
#' data("pancreas1k")
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' pancreas1k <- CellScoring(srt = pancreas1k, features = list(S = ccgenes$cc_S_genes, G2M = ccgenes$cc_G2M_genes), nbin = 20)
#' ClassDimPlot(pancreas1k, "CC_classification")
#' ExpDimPlot(pancreas1k, "CC_G2M")
#' @importFrom Seurat AddModuleScore AddMetaData
#' @export
CellScoring <- function(srt, features, ncores = 1, method = "Seurat", classification = TRUE,
                        name = "CC", slot = "data", assay = "RNA", new_assay = FALSE, ...) {
  if (!is.list(features) || length(names(features)) == 0) {
    stop("'features' must be named list")
  }
  if (!method %in% c("Seurat", "AUCell", "UCell")) {
    stop("method must be 'Seurat', 'AUCell'or 'UCell'.")
  }
  status <- check_DataType(srt, slot = "data")
  if (status != "log_normalized_counts") {
    warning("object is not log normalized")
  }
  if (method == "Seurat") {
    ## need to add a 'slot' parameter
    srt_tmp <- AddModuleScore(
      srt,
      features = features,
      name = name,
      assay = assay,
      ...
    )
    scores <- srt_tmp[[paste0(name, 1:length(features))]]
  } else if (method == "UCell") {
    check_R(pkgs = "UCell")
    srt_tmp <- UCell::AddModuleScore_UCell(
      srt,
      features = features,
      name = name,
      ncores = ncores,
      slot = slot,
      assay = assay,
      ...
    )
    scores <- srt_tmp[[paste0(names(features), name)]]
  } else if (method == "AUCell") {
    check_R(pkgs = "AUCell")
    CellRank <- AUCell::AUCell_buildRankings(as.matrix(GetAssayData(srt, slot = slot, assay = assay)), plotStats = FALSE)
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets = features, rankings = CellRank, nCores = ncores)
    scores <- as.data.frame(t(AUCell::getAUC(cells_AUC)))[, names(features)]
  }
  if (isTRUE(new_assay)) {
    srt[[name]] <- CreateAssayObject(counts = t(scores))
  } else {
    colnames(scores) <- make.names(paste(name, names(features), sep = "_"))
    srt <- AddMetaData(object = srt, metadata = scores)
  }

  if (isTRUE(classification)) {
    assignments <- apply(scores, MARGIN = 1, FUN = function(x) {
      if (all(x < 0)) {
        return(NA)
      } else {
        if (length(which(x == max(x))) > 1) {
          return(NA)
        } else {
          return(names(features)[which(x == max(x))])
        }
      }
    })
    srt[[paste0(name, "_classification")]] <- assignments[rownames(scores)]
  }

  return(srt)
}

#' Differential gene test
#'
#' @param srt A \code{Seurat} object
#' @param FindAllMarkers Whether to find markers across all cells.
#' @param FindPairedMarkers Whether to find markers between pairs of the cell groups.
#' @param group_by group_by
#' @param cell_group1 cell_group1
#' @param cell_group2 cell_group2
#' @param fc.threshold fc.threshold
#' @param min.pct min.pct
#' @param cl cl
#' @param force force
#'
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr bind_rows
#' @importFrom Seurat FindMarkers Assays
#' @importFrom stats p.adjust
#'
#' @examples
#' library(dplyr)
#' data(pancreas1k)
#' pancreas1k <- RunDEtest(pancreas1k, group_by = "CellType")
#'
#' # Heatmap
#' de_filter <- filter(pancreas1k@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' ExpHeatmapPlot(pancreas1k, genes = de_filter$gene, gene_groups = de_filter$group1, cell_group_by = "CellType")
#'
#' # Dot plot
#' de_top <- de_filter %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ExpDotPlot(pancreas1k, genes = de_top$gene, gene_groups = de_top$group1, cell_group_by = "CellType")
#' @export
#'
RunDEtest <- function(srt, group_by = NULL, cell_group1 = NULL, cell_group2 = NULL,
                      FindAllMarkers = TRUE, FindPairedMarkers = FALSE,
                      test.use = "wilcox", only.pos = TRUE,
                      fc.threshold = 1.5, min.pct = 0.1, max.cells.per.ident = Inf, latent.vars = NULL,
                      slot = "data", assay = "RNA", BPPARAM = BiocParallel::bpparam(), progressbar = TRUE, force = FALSE, ...) {
  if ("progressbar" %in% names(BPPARAM)) {
    BPPARAM[["progressbar"]] <- progressbar
  }

  time_start <- Sys.time()
  cat(paste0("[", time_start, "] ", "Start DEtest\n"))
  message("Threads used for DE test: ", BPPARAM$workers)

  if (!is.null(cell_group1)) {
    if (!all(cell_group1 %in% colnames(srt))) {
      stop("cell_group1 has some cells not in the Seurat object.")
    }
    if (is.null(cell_group2)) {
      message("cell_group2 is not provided. Use the remaining cells.")
      cell_group2 <- colnames(srt)[!colnames(srt) %in% cell_group1]
    }
    if (!all(cell_group2 %in% colnames(srt))) {
      stop("cell_group2 has some cells not in the Seurat object.")
    }
    if (length(cell_group1) < 3 | length(cell_group2) < 3) {
      stop("Cell groups must have more than 3 cells")
    }
    cat("Perform FindMarkers(", test.use, ") for custom cell groups...\n", sep = "")
    markers <- FindMarkers(
      object = Assays(srt, assay), slot = slot,
      cells.1 = cell_group1,
      cells.2 = cell_group2,
      test.use = test.use,
      logfc.threshold = log2(fc.threshold),
      min.pct = min.pct,
      max.cells.per.ident = max.cells.per.ident,
      latent.vars = latent.vars,
      only.pos = only.pos,
      verbose = FALSE,
      ...
    )
    if (nrow(markers) > 0) {
      markers[, "gene"] <- rownames(markers)
      markers[, "group1"] <- "cell_group1"
      markers[, "group2"] <- "cell_group2"
      rownames(markers) <- NULL
      markers[, "group1"] <- factor(markers[, "group1"], levels = unique(markers[, "group1"]))
      if ("p_val" %in% colnames(markers)) {
        markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = "BH")
      }
      markers[, "DE_group_number"] <- as.integer(table(markers[["gene"]])[markers[, "gene"]])
      MarkersMatrix <- as.data.frame.matrix(table(markers[, c("gene", "group1")]))
      markers[, "DE_group"] <- apply(MarkersMatrix, 1, function(x) {
        paste0(colnames(MarkersMatrix)[x > 0], collapse = ";")
      })[markers[, "gene"]]
      srt@tools[["DEtest_custom"]][[paste0("custom_", test.use)]] <- markers
      srt@tools[["DEtest_custom"]][["cell_group1"]] <- cell_group1
      srt@tools[["DEtest_custom"]][["cell_group2"]] <- cell_group2
    } else {
      warning("No markers found.", immediate. = TRUE)
    }
  } else {
    if (is.null(group_by)) {
      cell_group <- Idents(srt)
      group_by <- "active.ident"
    } else {
      cell_group <- srt[[group_by, drop = TRUE]]
    }
    if (!is.factor(cell_group)) {
      cell_group <- factor(cell_group, levels = unique(cell_group))
    }
    cell_group <- lapply(levels(cell_group), function(x) {
      cell <- cell_group[cell_group == x]
      out <- sample(cell, size = min(max.cells.per.ident, length(cell)), replace = FALSE)
      return(out)
    })
    cell_group <- setNames(unlist(lapply(cell_group, function(x) x), use.names = F), unlist(lapply(cell_group, names)))

    if (isTRUE(FindAllMarkers)) {
      cat("Perform FindAllMarkers(", test.use, ")...\n", sep = "")
      AllMarkers <- bplapply(levels(cell_group), FUN = function(group) {
        cells.1 <- names(cell_group)[which(cell_group == group)]
        cells.2 <- names(cell_group)[which(cell_group != group)]
        # print(paste0(group," cells.1:",length(cells.1)," cells.2:",length(cells.2)))
        if (length(cells.1) < 3 | length(cells.2) < 3) {
          return(NULL)
        } else {
          markers <- FindMarkers(
            object = Assays(srt, assay), slot = slot,
            cells.1 = cells.1,
            cells.2 = cells.2,
            test.use = test.use,
            logfc.threshold = log2(fc.threshold),
            min.pct = min.pct,
            max.cells.per.ident = Inf,
            latent.vars = latent.vars,
            only.pos = only.pos,
            verbose = FALSE,
            ...
          )
          if (nrow(markers) > 0) {
            markers[, "gene"] <- rownames(markers)
            markers[, "group1"] <- as.character(group)
            markers[, "group2"] <- "others"
            return(markers)
          } else {
            return(NULL)
          }
        }
      }, BPPARAM = BPPARAM)
      AllMarkers <- do.call(rbind.data.frame, AllMarkers)
      rownames(AllMarkers) <- NULL
      AllMarkers[, "group1"] <- factor(AllMarkers[, "group1"], levels = levels(cell_group))
      if ("p_val" %in% colnames(AllMarkers)) {
        AllMarkers[, "p_val_adj"] <- p.adjust(AllMarkers[, "p_val"], method = "BH")
      }
      AllMarkers[, "DE_group_number"] <- as.integer(table(AllMarkers[["gene"]])[AllMarkers[, "gene"]])
      AllMarkersMatrix <- as.data.frame.matrix(table(AllMarkers[, c("gene", "group1")]))
      AllMarkers[, "DE_group"] <- apply(AllMarkersMatrix, 1, function(x) {
        paste0(colnames(AllMarkersMatrix)[x > 0], collapse = ";")
      })[AllMarkers[, "gene"]]
      srt@tools[[paste0("DEtest_", group_by)]][[paste0("AllMarkers_", test.use)]] <- AllMarkers
    }

    if (isTRUE(FindPairedMarkers)) {
      if (nlevels(cell_group) > 30 & (!isTRUE(force))) {
        warning("Too many groups for FindPairedMarkers function. If you want to force to run, set force=TRUE.", immediate. = TRUE)
      } else {
        cat("Perform FindPairedMarkers(", test.use, ")...\n", sep = "")
        pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
        pair <- pair[pair[, 1] != pair[, 2], ]
        PairedMarkers <- bplapply(1:nrow(pair), function(i) {
          cells.1 <- names(cell_group)[which(cell_group == pair[i, 1])]
          cells.2 <- names(cell_group)[which(cell_group == pair[i, 2])]
          if (length(cells.1) < 3 | length(cells.2) < 3) {
            return(NULL)
          } else {
            markers <- FindMarkers(
              object = Assays(srt, assay), slot = slot,
              cells.1 = cells.1,
              cells.2 = cells.2,
              test.use = test.use,
              logfc.threshold = log2(fc.threshold),
              min.pct = min.pct,
              max.cells.per.ident = Inf,
              latent.vars = latent.vars,
              only.pos = only.pos,
              verbose = FALSE,
              ...
            )
            if (nrow(markers) > 0) {
              markers[, "gene"] <- rownames(markers)
              markers[, "group1"] <- as.character(pair[i, 1])
              markers[, "group2"] <- as.character(pair[i, 2])
              return(markers)
            } else {
              return(NULL)
            }
          }
        }, BPPARAM = BPPARAM)
        PairedMarkers <- do.call(rbind.data.frame, PairedMarkers)
        rownames(PairedMarkers) <- NULL
        PairedMarkers[, "group1"] <- factor(PairedMarkers[, "group1"], levels = levels(cell_group))
        if ("p_val" %in% colnames(PairedMarkers)) {
          PairedMarkers[, "p_val_adj"] <- p.adjust(PairedMarkers[, "p_val"], method = "BH")
        }
        PairedMarkers[, "DE_group_number"] <- as.integer(table(PairedMarkers[["gene"]])[PairedMarkers[, "gene"]])
        PairedMarkersMatrix <- as.data.frame.matrix(table(PairedMarkers[, c("gene", "group1")]))
        PairedMarkers[, "DE_group"] <- apply(PairedMarkersMatrix, 1, function(x) {
          paste0(colnames(PairedMarkersMatrix)[x > 0], collapse = ";")
        })[PairedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkers_", test.use)]] <- PairedMarkers
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkersMatrix_", test.use)]] <- PairedMarkersMatrix
      }
    }
  }
  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "DEtest done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")
  return(srt)
}


#' Prepare the database for enrichment analysis
#' @param species species
#' @param enrichment Enrichment database name.
#' @param db_update Whether update the database.
#'
#' @return A list containing the database.
#'
#' @examples
#' if (interactive()) {
#'   db_list <- PrepareEnrichmentDB(species = "Homo_sapiens")
#' }
#' @importFrom R.cache loadCache saveCache readCacheHeader findCache
#' @importFrom utils packageVersion read.table
#' @importFrom stats na.omit
#' @export
#'
PrepareEnrichmentDB <- function(species = "Homo_sapiens",
                                enrichment = c(
                                  "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome",
                                  "ProteinComplex", "PFAM", "Chromosome"
                                ),
                                db_IDtypes = c("symbol", "entrez_id", "ensembl_id"), db_update = FALSE,
                                Ensembl_version = 103, mirror = NULL) {
  if (length(species) > 1) {
    stop("Only one species name can be provided at a time")
  }
  message("Species: ", species)

  db_list <- list()
  if (!isTRUE(db_update)) {
    for (term in enrichment) {
      # Try to load cached database, if already generated.
      key <- list(species, term)
      db <- loadCache(key)
      if (!is.null(db)) {
        header <- readCacheHeader(findCache(key = key))
        message("Loaded cached db: ", term, " version:", header$comment, " created:", header$timestamp)
        db_list[[species]][[term]] <- db
      }
    }
  }

  sp <- unlist(strsplit(species, split = "_"))
  org_sp <- paste0("org.", paste0(substring(sp, 1, 1), collapse = ""), ".eg.db")
  # mesh_sp <-  paste0("MeSH.", paste0(substring(sp, 1, c(1,2)), collapse = ""), ".eg.db")
  kegg_sp <- tolower(paste0(substring(sp[1], 1, 1), substring(sp[2], 1, 2), collapse = ""))
  complex_sp <- reactome_sp <- wiki_sp <- gsub(pattern = "_", replacement = " ", x = species)

  orgdb_dependent <- c("GO_BP", "GO_CC", "GO_MF", "PFAM", "Chromosome")
  if (any(orgdb_dependent %in% enrichment)) {
    status <- tryCatch(
      {
        check_R(c("AnnotationDbi", org_sp, "GO.db", "GOSemSim"))
      },
      error = function(e) {
        return("error")
      }
    )
    if (identical(status, "error")) {
      check_R("AnnotationDbi", force = TRUE)
      check_R(c(org_sp, "GO.db", "GOSemSim"))
    }
    suppressPackageStartupMessages(library(org_sp, character.only = TRUE))
    orgdb <- get(org_sp)
    orgdbCHR <- get(paste0(gsub(pattern = ".db", "", org_sp), "CHR"))
  }
  if ("PFAM" %in% enrichment) {
    check_R(c("PFAM.db", "AnnotationDbi"))
  }
  if ("Reactome" %in% enrichment) {
    check_R(c("reactome.db", "AnnotationDbi"))
  }
  if ("WikiPathway" %in% enrichment) {
    check_R(c("rWikiPathways", "clusterProfiler"))
  }
  if ("MeSH" %in% enrichment) {
    check_R(c("AHMeSHDbs", "MeSHDbi", "MeSH.db", "AnnotationHub"))
  }

  ## Prepare ----------------------------------------------------------------------
  if (any(!species %in% names(db_list)) || any(!enrichment %in% names(db_list[[species]]))) {
    ## GO ---------------------------------------------------------------------------------
    if (any(enrichment %in% c("GO_BP", "GO_CC", "GO_MF")) && any(!c("GO_BP", "GO_CC", "GO_MF") %in% names(db_list[[species]]))) {
      terms <- enrichment[enrichment %in% c("GO_BP", "GO_CC", "GO_MF")]
      for (subterm in terms) {
        message("Preparing database: ", subterm)
        sub_enrichment <- unlist(strsplit(subterm, split = "_"))[2]
        bg <- suppressMessages(AnnotationDbi::select(orgdb, keys = AnnotationDbi::keys(orgdb), columns = c("GOALL")))
        bg$EVIDENCEALL <- NULL
        bg <- unique(bg[!is.na(bg$GOALL), ])
        bg2 <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = AnnotationDbi::keys(GO.db::GO.db), columns = c("GOID", "TERM", "DEFINITION")))
        bg <- merge(x = bg, by.x = "GOALL", y = bg2, by.y = "GOID", all.x = TRUE)
        TERM2GENE <- bg[which(bg$ONTOLOGYALL == sub_enrichment), c(1, 2)]
        TERM2NAME <- bg[which(bg$ONTOLOGYALL == sub_enrichment), c(1, 4)]
        colnames(TERM2GENE) <- c("Term", "entrez_id")
        colnames(TERM2NAME) <- c("Term", "Name")
        semData <- suppressMessages(GOSemSim::godata(orgdb, ont = sub_enrichment))
        version <- packageVersion(org_sp)
        db_list[[species]][[subterm]][["TERM2GENE"]] <- unique(TERM2GENE)
        db_list[[species]][[subterm]][["TERM2NAME"]] <- unique(TERM2NAME)
        db_list[[species]][[subterm]][["semData"]] <- semData
        db_list[[species]][[subterm]][["version"]] <- version
        saveCache(db_list[[species]][[subterm]],
          key = list(species, subterm),
          comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
        )
      }
    }

    ## KEGG ---------------------------------------------------------------------------
    if (any(enrichment == "KEGG") && (!"KEGG" %in% names(db_list[[species]]))) {
      message("Preparing database: KEGG")
      kegg_db <- "pathway"
      kegg_pathwaygene_url <- paste0("http://rest.kegg.jp/link/", kegg_sp, "/", kegg_db, collapse = "")
      TERM2GENE <- kegg_get(kegg_pathwaygene_url)
      TERM2GENE[, 1] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, 1])
      TERM2GENE[, 2] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, 2])
      kegg_pathwayname_url <- paste0("http://rest.kegg.jp/list/", kegg_db, collapse = "")
      TERM2NAME <- kegg_get(kegg_pathwayname_url)
      TERM2NAME[, 1] <- gsub(pattern = "path:map", replacement = kegg_sp, TERM2NAME[, 1])
      TERM2NAME <- TERM2NAME[TERM2NAME[, 1] %in% TERM2GENE[, 1], ]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      kegg_info <- readLines("http://rest.kegg.jp/info/hsa")
      version <- kegg_info[grepl("Release", x = kegg_info)] %>% gsub(".*(?=Release)", replacement = "", x = ., perl = TRUE)
      db_list[[species]][["KEGG"]][["TERM2GENE"]] <- unique(TERM2GENE)
      db_list[[species]][["KEGG"]][["TERM2NAME"]] <- unique(TERM2NAME)
      db_list[[species]][["KEGG"]][["version"]] <- version
      saveCache(db_list[[species]][["KEGG"]],
        key = list(species, "KEGG"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    ## WikiPathway ---------------------------------------------------------------------------
    if (any(enrichment == "WikiPathway") && (!"WikiPathway" %in% names(db_list[[species]]))) {
      message("Preparing database: WikiPathway")
      tempdir <- tempdir()
      gmt_files <- list.files(tempdir)[grep(".gmt", x = list.files(tempdir))]
      if (length(gmt_files) > 0) {
        file.remove(paste0(tempdir, "/", gmt_files))
      }
      rWikiPathways::downloadPathwayArchive(organism = wiki_sp, format = "gmt", date = "current", destpath = tempdir)
      "https://wikipathways-data.wmcloud.org/current/gmt/"
      version <- list.files(tempdir)[grep(".gmt", x = list.files(tempdir))] %>%
        strsplit(split = "-") %>%
        unlist() %>%
        .[2]
      wiki_gmt <- clusterProfiler::read.gmt(paste0(tempdir, "/", list.files(tempdir)[grep(".gmt", x = list.files(tempdir))]))
      wiki_gmt <- apply(wiki_gmt, 1, function(x) {
        wikiid <- strsplit(x[["term"]], split = "%")[[1]][3]
        wikiterm <- strsplit(x[["term"]], split = "%")[[1]][1]
        gmt <- x[["gene"]]
        data.frame(v0 = wikiid, v1 = gmt, v2 = wikiterm, stringsAsFactors = FALSE)
      })
      bg <- bind_rows(wiki_gmt)
      TERM2GENE <- bg[, c(1, 2)]
      TERM2NAME <- bg[, c(1, 3)]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      db_list[[species]][["WikiPathway"]][["TERM2GENE"]] <- unique(TERM2GENE)
      db_list[[species]][["WikiPathway"]][["TERM2NAME"]] <- unique(TERM2NAME)
      db_list[[species]][["WikiPathway"]][["version"]] <- version
      saveCache(db_list[[species]][["WikiPathway"]],
        key = list(species, "WikiPathway"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    ## pathwaycommons ---------------------------------------------------------------------------
    # check_R("paxtoolsr")

    ## Reactome ---------------------------------------------------------------------------
    if (any(enrichment == "Reactome") && (!"Reactome" %in% names(db_list[[species]]))) {
      message("Preparing database: Reactome")
      bg <- suppressMessages(AnnotationDbi::select(reactome.db::reactome.db, keys = AnnotationDbi::keys(reactome.db::reactome.db), columns = c("PATHID", "PATHNAME")))
      bg <- bg[grepl(pattern = reactome_sp, x = bg$PATHNAME), ]
      bg <- na.omit(bg)
      bg$PATHNAME <- gsub(x = bg$PATHNAME, pattern = paste0("^", reactome_sp, ": "), replacement = "", perl = TRUE)
      TERM2GENE <- bg[, c(2, 1)]
      TERM2NAME <- bg[, c(2, 3)]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      version <- packageVersion("reactome.db")
      db_list[[species]][["Reactome"]][["TERM2GENE"]] <- unique(TERM2GENE)
      db_list[[species]][["Reactome"]][["TERM2NAME"]] <- unique(TERM2NAME)
      db_list[[species]][["Reactome"]][["version"]] <- version
      saveCache(db_list[[species]][["Reactome"]],
        key = list(species, "Reactome"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    ## Protein complex ---------------------------------------------------------------------------
    if (any(enrichment == "ProteinComplex") && (!"ProteinComplex" %in% names(db_list[[species]]))) {
      message("Preparing database: ProteinComplex")
      check_R("taxize")
      uid <- suppressMessages(taxize::get_uid(complex_sp, messages = FALSE))
      common_name <- suppressMessages(taxize::sci2comm(uid, db = "ncbi")) %>%
        unlist() %>%
        strsplit(" ") %>%
        unlist()
      common_name <- common_name[length(common_name)]
      common_name <- paste(toupper(substr(common_name, 1, 1)), substr(common_name, 2, nchar(common_name)), sep = "")

      temp <- tempfile()
      download.file("http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt.zip", temp)
      df <- read.table(unz(temp, "coreComplexes.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
      unlink(temp)
      df <- df[which(df$Organism == common_name), ]
      s <- strsplit(df$subunits.Entrez.IDs., split = ";")
      complex <- data.frame(V1 = rep(df$ComplexName, sapply(s, length)), V2 = unlist(s), V3 = rep(paste0("ComplexID:", df$ComplexID), sapply(s, length)))
      complex$V1 <- trimws(complex$V1) %>% gsub(pattern = "\\([^\\)]*\\)$", replacement = "", perl = FALSE)
      complex <- complex[!duplicated(complex), ]
      complex <- complex[which(complex$V2 != "NULL"), ]
      TERM2GENE <- complex[, c(3, 2)]
      TERM2NAME <- complex[, c(3, 1)]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      temp <- tempfile()
      download.file("http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.xml.zip", temp)
      con <- unz(temp, "coreComplexes.xml")
      complex_info <- readLines(con)
      close(con)
      version <- complex_info[grepl("releaseDate", x = complex_info)] %>%
        strsplit(x = ., split = "\"") %>%
        unlist() %>%
        .[2]
      unlink(temp)
      db_list[[species]][["ProteinComplex"]][["TERM2GENE"]] <- unique(TERM2GENE)
      db_list[[species]][["ProteinComplex"]][["TERM2NAME"]] <- unique(TERM2NAME)
      db_list[[species]][["ProteinComplex"]][["version"]] <- version
      saveCache(db_list[[species]][["ProteinComplex"]],
        key = list(species, "ProteinComplex"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    ## PFAM ---------------------------------------------------------------------------
    if (any(enrichment == "PFAM") && (!"PFAM" %in% names(db_list[[species]]))) {
      message("Preparing database: PFAM")
      bg <- suppressMessages(AnnotationDbi::select(orgdb, keys = AnnotationDbi::keys(orgdb), columns = "PFAM"))
      bg <- unique(bg[!is.na(bg$PFAM), ])
      bg2 <- as.data.frame(PFAM.db::PFAMDE2AC[AnnotationDbi::mappedkeys(PFAM.db::PFAMDE2AC)])
      rownames(bg2) <- bg2[["ac"]]
      bg[["PFAM_name"]] <- bg2[bg$PFAM, "de"]
      bg[is.na(bg[["PFAM_name"]]), "PFAM_name"] <- bg[is.na(bg[["PFAM_name"]]), "PFAM"]
      TERM2GENE <- bg[, c(2, 1)]
      TERM2NAME <- bg[, c(2, 3)]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      version <- packageVersion(org_sp)
      db_list[[species]][["PFAM"]][["TERM2GENE"]] <- TERM2GENE
      db_list[[species]][["PFAM"]][["TERM2NAME"]] <- TERM2NAME
      db_list[[species]][["PFAM"]][["version"]] <- version
      saveCache(db_list[[species]][["PFAM"]],
        key = list(species, "PFAM"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    ## Chromosome ---------------------------------------------------------------------------
    if (any(enrichment == "Chromosome") && (!"Chromosome" %in% names(db_list[[species]]))) {
      message("Preparing database: Chromosome")
      chr <- as.data.frame(orgdbCHR[AnnotationDbi::mappedkeys(orgdbCHR)])
      chr[, 2] <- paste0("chr", chr[, 2])
      TERM2GENE <- chr[, c(2, 1)]
      TERM2NAME <- chr[, c(2, 2)]
      colnames(TERM2GENE) <- c("Term", "entrez_id")
      colnames(TERM2NAME) <- c("Term", "Name")
      version <- packageVersion(org_sp)
      db_list[[species]][["Chromosome"]][["TERM2GENE"]] <- unique(TERM2GENE)
      db_list[[species]][["Chromosome"]][["TERM2NAME"]] <- unique(TERM2NAME)
      db_list[[species]][["Chromosome"]][["version"]] <- version
      saveCache(db_list[[species]][["Chromosome"]],
        key = list(species, "Chromosome"),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }

    # ## MeSH ---------------------------------------------------------------------------
    # if (any(enrichment == "MeSH") && (!"MeSH" %in% names(db_list[[species]]))) {
    #   message("Preparing database: MeSH")
    #   # dir.create("~/.cache/R/AnnotationHub",recursive = TRUE,showWarnings = FALSE)
    #   ### A (Anatomy);B (Organisms);C (Diseases);D (Chemicals and Drugs);
    #   ### E (Analytical Diagnostic and Therapeutic Techniques and Equipment);F (Psychiatry and Psychology);
    #   ### G (Phenomena and Processes);H (Disciplines and Occupations);
    #   ### I (Anthropology, Education, Sociology and Social Phenomena);J (Technology and Food and Beverages);
    #   ### K (Humanities);L (Information Science);M (Persons);N (Health Care);
    #   ### V (Publication Type);Z (Geographical Locations)
    #
    #   ah <- AnnotationHub()
    #   qr <- query(ah, c(paste(c("MeSHDb for", sp), collapse = " "), "v001"))
    #   if (length(qr) == 0) {
    #     stop("no MeSH records found for ", sp)
    #   }
    #   if (length(qr) == 1) {
    #     db <- qr[[1]]
    #     <- RSQLite::dbConnect(db)
    #     library("MeSHDbi")
    #     meshdb <- MeSHDbi::MeSHDb(db)
    #     MeSHDbi::dbconn(meshdb)
    #     MeSHDbi::dbschema(meshdb)
    #     MeSHDbi::dbconn(meshdb)
    #   }
    #
    #   select(file, keys = keys(file), columns = columns(file))
    #
    #
    #   bg <- AnnotationDbi::select(meshdb, keys = keys(meshdb, keytype = "GENEID"), keytype = "GENEID", columns = c("GENEID", "MESHID", "MESHCATEGORY"))
    #   # bg <-  bg[which(bg$GENEID %in% keys(orgdb)),]
    #   # bg_all <- AnnotationDbi::select(MeSH.db, keys=keys(MeSH.db,keytype = "MESHID"),columns = c("MESHID","MESHTERM"),keytype = "MESHID")
    #   # saveRDS(bg_all,"./MeSHID2Term.rds")
    #   bg_all <- readRDS("./MeSHID2Term.rds")
    #   bg <- merge(x = bg, by.x = "MESHID", y = bg_all, by.y = "MESHID", all.x = T)
    #   bg[which(is.na(bg$MESHTERM)), "MESHTERM"] <- bg[which(is.na(bg$MESHTERM)), "MESHID"]
    #   assign(paste0(species, "_meshall"), bg)
    #
    #   bg <- suppressMessages(AnnotationDbi::select(meshdb, keys = AnnotationDbi::keys(meshdb), columns = "PFAM"))
    #   bg <- unique(bg[!is.na(bg$PFAM), ])
    #   bg2 <- as.data.frame(PFAM.db::PFAMDE2AC[AnnotationDbi::mappedkeys(PFAM.db::PFAMDE2AC)])
    #   rownames(bg2) <- bg2[["ac"]]
    #   bg[["PFAM_name"]] <- bg2[bg$PFAM, "de"]
    #   bg[is.na(bg[["PFAM_name"]]), "PFAM_name"] <- bg[is.na(bg[["PFAM_name"]]), "PFAM"]
    #   TERM2GENE <- bg[, c(2, 1)]
    #   TERM2NAME <- bg[, c(2, 3)]
    #   version <- packageVersion(org_sp)
    #   db_list[[species]][["PFAM"]][["TERM2GENE"]] <- TERM2GENE
    #   db_list[[species]][["PFAM"]][["TERM2NAME"]] <- TERM2NAME
    #   db_list[[species]][["PFAM"]][["version"]] <- version
    #   saveCache(db_list[[species]][["PFAM"]],
    #     key = list(species, "PFAM"),
    #     comment = paste0(version, " nterm:", length(unique(TERM2NAME[[1]])))
    #   )
    # }
  }

  ### Convert ID
  for (term in names(db_list[[species]])) {
    IDtypes <- db_IDtypes[!db_IDtypes %in% colnames(db_list[[species]][[term]][["TERM2GENE"]])]
    if (length(IDtypes) > 0) {
      message("Convert ID types for the database: ", term)
      TERM2GENE <- db_list[[species]][[term]][["TERM2GENE"]]
      TERM2NAME <- db_list[[species]][[term]][["TERM2NAME"]]
      res <- GeneConvert(
        geneID = as.character(unique(TERM2GENE[, "entrez_id"])),
        geneID_from_IDtype = "entrez_id",
        geneID_to_IDtype = IDtypes,
        species_from = species,
        species_to = species,
        Ensembl_version = Ensembl_version,
        mirror = mirror
      )
      map <- res$geneID_collapse
      for (type in IDtypes) {
        TERM2GENE[[type]] <- map[as.character(TERM2GENE[, "entrez_id"]), type]
        TERM2GENE <- as.data.frame(tidyr::unnest(TERM2GENE, all_of(type), keep_empty = TRUE))
      }
      db_list[[species]][[term]][["TERM2GENE"]] <- TERM2GENE
      ### save cache
      version <- db_list[[species]][[term]][["version"]]
      saveCache(db_list[[species]][[term]],
        key = list(species, term),
        comment = paste0(version, " nterm:", length(TERM2NAME[[1]]))
      )
    }
  }

  return(db_list)
}

#' Perform the enrichment analysis(over-representation) on the genes
#'
#' @param geneID geneID
#' @param geneID_groups geneID_groups
#' @param IDtype IDtype
#' @param result_IDtype result_IDtype
#' @param species species
#' @param enrichment enrichment
#' @param db_update db_update
#' @param GO_simplify GO_simplify
#' @param GO_simplify_padjustCutoff GO_simplify_padjustCutoff
#' @param simplify_method simplify_method
#' @param simplify_similarityCutoff simplify_similarityCutoff
#' @param TERM2GENE TERM2GENE
#' @param TERM2NAME TERM2NAME
#' @param minGSSize minGSSize
#' @param maxGSSize maxGSSize
#' @param universe universe
#' @param cl cores
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas1k")
#'   library(dplyr)
#'   pancreas1k <- RunDEtest(pancreas1k, group_by = "CellType")
#'   de_filter <- filter(pancreas1k@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05)
#'   res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, enrichment = "GO_BP", species = "Mus_musculus")
#'   EnrichmentPlot(res$enrichment, plot_type = "bar")
#'   EnrichmentPlot(res$enrichment, plot_type = "bar", character_width = 30, text_size_scale = 0.8)
#'
#'   EnrichmentPlot(filter(res$enrichment, Groups %in% c("Ductal", "Endocrine")), plot_type = "lollipop", ncol = 1)
#'   EnrichmentPlot(filter(res$enrichment, Groups %in% c("Ductal", "Endocrine")), plot_type = "wordcloud")
#' }
#' @importFrom BiocParallel bplapply
#' @export
#'
RunEnrichment <- function(geneID = NULL, geneID_groups = NULL, IDtype = "symbol", result_IDtype = "symbol", species = "Homo_sapiens",
                          enrichment = "GO_BP", db_IDtype = "symbol", db_update = FALSE, Ensembl_version = 103, mirror = NULL,
                          TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, universe = NULL,
                          GO_simplify = TRUE, GO_simplify_padjustCutoff = 0.05, simplify_method = "Rel", simplify_similarityCutoff = 0.7,
                          BPPARAM = BiocParallel::bpparam(), progressbar = TRUE) {
  if ("progressbar" %in% names(BPPARAM)) {
    BPPARAM[["progressbar"]] <- progressbar
  }

  time_start <- Sys.time()
  cat(paste0("[", time_start, "] ", "Start Enrichment\n"))
  if (is.null(geneID_groups)) {
    geneID_groups <- rep(" ", length(geneID))
  }
  if (!is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  if (length(geneID_groups) != length(geneID)) {
    stop("length(geneID_groups)!=length(geneID)")
  }
  names(geneID_groups) <- geneID
  input <- data.frame(geneID = geneID, geneID_groups = geneID_groups)

  if (is.null(TERM2GENE)) {
    db_list <- PrepareEnrichmentDB(
      species = species, enrichment = enrichment, db_update = db_update,
      db_IDtypes = db_IDtype, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    enrichment <- "custom"
    db_list <- list()
    db_list[[species]][[enrichment]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[enrichment]][["TERM2NAME"]] <- unique(TERM2NAME)
  }

  if (length(unique(c(IDtype, db_IDtype, result_IDtype))) != 1) {
    res <- GeneConvert(
      geneID = unique(geneID),
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = unique(c(db_IDtype, result_IDtype)),
      species_from = species,
      species_to = species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
    geneMap <- res$geneID_collapse
    geneMap <- geneMap[!is.na(geneMap[, db_IDtype]), ]
  } else {
    geneMap <- data.frame(from_geneID = unique(geneID), db_IDtype = unique(geneID), row.names = unique(geneID))
    colnames(geneMap)[2] <- db_IDtype
  }

  input[[db_IDtype]] <- geneMap[as.character(input$geneID), db_IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- as.data.frame(tidyr::unnest(input, all_of(c(db_IDtype, result_IDtype))))
  input <- input[!is.na(input[[db_IDtype]]), ]

  message("Permform enrichment...")
  comb <- expand.grid(group = levels(geneID_groups), term = enrichment, stringsAsFactors = FALSE)
  res_list <- bplapply(1:nrow(comb),
    FUN = function(i, id) {
      group <- comb[i, "group"]
      term <- comb[i, "term"]
      gene <- input[input$geneID_groups == group, db_IDtype]
      gene_mapid <- input[input$geneID_groups == group, result_IDtype]
      TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", db_IDtype)]
      TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), ]
      TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], ]
      enrich_res <- suppressPackageStartupMessages({
        clusterProfiler::enricher(
          gene = gene,
          minGSSize = ifelse(term %in% c("Chromosome"), 1, minGSSize),
          maxGSSize = ifelse(term %in% c("Chromosome"), Inf, maxGSSize),
          pAdjustMethod = "BH",
          pvalueCutoff = Inf,
          qvalueCutoff = Inf,
          universe = universe,
          TERM2GENE = TERM2GENE_tmp,
          TERM2NAME = TERM2NAME_tmp
        )
      })
      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        result <- enrich_res@result
        result[, "Enrichment"] <- term
        result[, "Groups"] <- group
        result[, "BgVersion"] <- as.character(db_list[[species]][[term]][["version"]])
        IDlist <- result$geneID %>%
          strsplit("/")
        result$geneID <- lapply(IDlist, function(x) {
          result_ID <- geneMap[geneMap[, db_IDtype] %in% x, result_IDtype]
          remain_ID <- x[!x %in% geneMap[, db_IDtype]]
          paste0(c(result_ID, remain_ID), collapse = "/")
        }) %>% unlist()
        enrich_res@result <- result
        enrich_res@gene2Symbol <- gene_mapid

        if (isTRUE(GO_simplify) && term %in% c("GO_BP", "GO_CC", "GO_MF")) {
          sim_res <- enrich_res
          sim_res@ontology <- gsub(pattern = "GO_", replacement = "", x = term)
          nterm_simplify <- sum(sim_res@result$p.adjust < GO_simplify_padjustCutoff)
          if (nterm_simplify < 1) {
            warning(group, "|", term, " has no term to simplify.", immediate. = TRUE)
          } else {
            sim_res@result <- sim_res@result[sim_res@result$p.adjust < GO_simplify_padjustCutoff, ]
            semData <- db_list[[species]][[term]][["semData"]]
            BiocParallel::ipclock(id)
            sim_res <- suppressPackageStartupMessages({
              clusterProfiler::simplify(sim_res,
                measure = simplify_method,
                cutoff = simplify_similarityCutoff,
                semData = semData
              )
            })
            BiocParallel::ipcunlock(id)
            result_sim <- sim_res@result
            result_sim[, "Enrichment"] <- paste0(term, "_sim")
            result_sim[, "Groups"] <- group
            result_sim[, "BgVersion"] <- as.character(db_list[[species]][[term]][["version"]])
            sim_res@result <- result_sim
            enrich_res <- list(enrich_res, sim_res)
            names(enrich_res) <- paste(group, c(term, paste0(term, "_sim")), sep = "-")
          }
        }
        return(enrich_res)
      } else {
        return(NULL)
      }
    }, BPPARAM = BPPARAM, id = BiocParallel::ipcid()
  )
  nm <- paste(comb$group, comb$term, sep = "-")
  sim_index <- sapply(res_list, function(x) length(x) == 2)
  sim_list <- unlist(res_list[sim_index], recursive = FALSE)
  raw_list <- res_list[!sim_index]
  names(raw_list) <- nm[!sim_index]
  results <- c(raw_list, sim_list)
  results <- results[!sapply(results, is.null)]
  results <- results[intersect(c(nm, paste0(nm, "_sim")), names(results))]
  res_enrichment <- bind_rows(lapply(results, function(x) x@result))
  rownames(res_enrichment) <- NULL

  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "Enrichment done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

  return(list(enrichment = res_enrichment, results = results, geneMap = geneMap, input = input))
}

#' Perform the enrichment analysis(GSEA) on the genes
#'
#' @param geneID geneID
#' @param geneID_groups geneID_groups
#' @param IDtype IDtype
#' @param result_IDtype result_IDtype
#' @param species species
#' @param enrichment enrichment
#' @param db_update db_update
#' @param GO_simplify GO_simplify
#' @param GO_simplify_padjustCutoff GO_simplify_padjustCutoff
#' @param simplify_method simplify_method
#' @param simplify_similarityCutoff simplify_similarityCutoff
#' @param TERM2GENE TERM2GENE
#' @param TERM2NAME TERM2NAME
#' @param minGSSize minGSSize
#' @param maxGSSize maxGSSize
#' @param universe universe
#' @param cl cores
#'
#' @importFrom BiocParallel bplapply
#' @examples
#' if (interactive()) {
#'   data("pancreas1k")
#'   library(dplyr)
#'   pancreas1k <- RunDEtest(pancreas1k, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
#'   de_filter <- filter(pancreas1k@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05)
#'   res <- RunGSEA(geneID = de_filter$gene, geneScore = de_filter$avg_log2FC, geneID_groups = de_filter$group1, species = "Mus_musculus")
#'   GSEAPlot(x = res$results[[1]], geneSetID = res$results[[1]]@result$ID[1])
#' }
#' @export
#'
RunGSEA <- function(geneID = NULL, geneScore = NULL, geneID_groups = NULL, IDtype = "symbol", result_IDtype = "symbol", species = "Homo_sapiens",
                    enrichment = "GO_BP", db_IDtype = "symbol", db_update = FALSE, Ensembl_version = 103, mirror = NULL,
                    TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, scoreType = "std",
                    GO_simplify = FALSE, GO_simplify_padjustCutoff = 0.05, simplify_method = "Rel", simplify_similarityCutoff = 0.7,
                    BPPARAM = BiocParallel::bpparam(), progressbar = TRUE) {
  if ("progressbar" %in% names(BPPARAM)) {
    BPPARAM[["progressbar"]] <- progressbar
  }

  time_start <- Sys.time()
  cat(paste0("[", time_start, "] ", "Start GSEA\n"))

  if (is.null(geneID_groups)) {
    geneID_groups <- rep(" ", length(geneID))
  }
  if (!is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  if (length(geneID_groups) != length(geneID)) {
    stop("length(geneID_groups)!=length(geneID)")
  }
  if (length(geneScore) != length(geneID)) {
    stop("geneScore must be the same length with geneID")
  }
  if (all(geneScore > 0) && scoreType != "pos") {
    warning("All values in the geneScore are greater than zero and scoreType is '", scoreType, "', maybe you should switch to scoreType = 'pos'.", immediate. = TRUE)
  }
  if (all(geneScore < 0) && scoreType != "neg") {
    warning("All values in the geneScore are less than zero and scoreType is '", scoreType, "', maybe you should switch to scoreType = 'neg'.", immediate. = TRUE)
  }
  input <- data.frame(geneID = geneID, geneScore = geneScore, geneID_groups = geneID_groups)

  na_index <- which(is.na(geneScore))
  if (length(na_index) > 0) {
    message("Ignore ", length(na_index), " NA geneScore")
    input <- input[-na_index, ]
  }
  input[is.infinite(input$geneScore) & input$geneScore < 0, "geneScore"] <- min(input[!is.infinite(input$geneScore), "geneScore"])
  input[is.infinite(input$geneScore) & input$geneScore > 0, "geneScore"] <- max(input[!is.infinite(input$geneScore), "geneScore"])

  geneID <- input$geneID
  geneScore <- input$geneScore
  geneID_groups <- input$geneID_groups
  names(geneID_groups) <- geneID
  names(geneScore) <- paste(geneID, geneID_groups, sep = ".")

  if (is.null(TERM2GENE)) {
    db_list <- PrepareEnrichmentDB(
      species = species, enrichment = enrichment, db_update = db_update,
      db_IDtypes = db_IDtype, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    enrichment <- "custom"
    db_list <- list()
    db_list[[species]][[enrichment]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[enrichment]][["TERM2NAME"]] <- unique(TERM2NAME)
  }

  if (length(unique(c(IDtype, db_IDtype, result_IDtype))) != 1) {
    res <- GeneConvert(
      geneID = unique(geneID),
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = unique(c(db_IDtype, result_IDtype)),
      species_from = species,
      species_to = species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
    geneMap <- res$geneID_collapse
    geneMap <- geneMap[!is.na(geneMap[, db_IDtype]), ]
  } else {
    geneMap <- data.frame(from_geneID = unique(geneID), db_IDtype = unique(geneID), row.names = unique(geneID))
    colnames(geneMap)[2] <- db_IDtype
  }

  input[[db_IDtype]] <- geneMap[as.character(input$geneID), db_IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- as.data.frame(tidyr::unnest(input, all_of(c(db_IDtype, result_IDtype))))
  input <- input[!is.na(input[[db_IDtype]]), ]

  message("Permform enrichment...")
  comb <- expand.grid(group = levels(geneID_groups), term = enrichment, stringsAsFactors = FALSE)
  res_list <- BiocParallel::bplapply(1:nrow(comb),
    FUN = function(i, id) {
      # print(paste0("current", i))
      group <- comb[i, "group"]
      term <- comb[i, "term"]
      geneList <- input[input$geneID_groups == group, "geneScore"]
      names(geneList) <- input[input$geneID_groups == group, db_IDtype]
      gene_mapid <- input[input$geneID_groups == group, result_IDtype]
      ord <- order(geneList, decreasing = TRUE)
      geneList <- geneList[ord]
      gene_mapid <- gene_mapid[ord]
      TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", db_IDtype)]
      TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), ]
      TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], ]
      enrich_res <- tryCatch(
        {
          suppressPackageStartupMessages({
            clusterProfiler::GSEA(
              geneList = geneList,
              minGSSize = ifelse(term %in% c("Chromosome"), 1, minGSSize),
              maxGSSize = ifelse(term %in% c("Chromosome"), Inf, maxGSSize),
              nPermSimple = 10000,
              eps = 0,
              scoreType = scoreType,
              pAdjustMethod = "BH",
              pvalueCutoff = Inf,
              TERM2GENE = TERM2GENE_tmp,
              TERM2NAME = TERM2NAME_tmp,
              by = "fgsea",
              verbose = FALSE
            )
          })
        },
        error = function(error) {
          warning("fgseaMultilevel failed to run. Try using fgseaSimple.", immediate. = TRUE)
          suppressPackageStartupMessages({
            clusterProfiler::GSEA(
              geneList = geneList,
              minGSSize = ifelse(term %in% c("Chromosome"), 1, minGSSize),
              maxGSSize = ifelse(term %in% c("Chromosome"), Inf, maxGSSize),
              nPerm = 10000,
              eps = 0,
              scoreType = scoreType,
              pAdjustMethod = "BH",
              pvalueCutoff = Inf,
              TERM2GENE = TERM2GENE_tmp,
              TERM2NAME = TERM2NAME_tmp,
              by = "fgsea",
              verbose = FALSE
            )
          })
        }
      )
      if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
        result <- enrich_res@result
        result[, "Enrichment"] <- term
        result[, "Groups"] <- group
        result[, "BgVersion"] <- as.character(db_list[[species]][[term]][["version"]])
        IDlist <- result$core_enrichment %>%
          strsplit("/")
        result$core_enrichment <- lapply(IDlist, function(x) {
          result_ID <- input[input[, db_IDtype] %in% x, result_IDtype]
          remain_ID <- x[!x %in% input[, db_IDtype]]
          paste0(c(result_ID, remain_ID), collapse = "/")
        }) %>% unlist()
        enrich_res@result <- result
        enrich_res@gene2Symbol <- gene_mapid

        if (isTRUE(GO_simplify) && term %in% c("GO_BP", "GO_CC", "GO_MF")) {
          sim_res <- enrich_res
          sim_res@setType <- gsub(pattern = "GO_", replacement = "", x = term)
          nterm_simplify <- sum(sim_res@result$p.adjust < GO_simplify_padjustCutoff)
          if (nterm_simplify < 1) {
            warning(group, "|", term, " has no term to simplify.", immediate. = TRUE)
          } else {
            sim_res@result <- sim_res@result[sim_res@result$p.adjust < GO_simplify_padjustCutoff, ]
            semData <- db_list[[species]][[term]][["semData"]]
            BiocParallel::ipclock(id)
            sim_res <- clusterProfiler::simplify(sim_res,
              measure = simplify_method,
              cutoff = simplify_similarityCutoff,
              semData = semData
            )
            BiocParallel::ipcunlock(id)
            result_sim <- sim_res@result
            result_sim[, "Enrichment"] <- paste0(term, "_sim")
            result_sim[, "Groups"] <- group
            result_sim[, "BgVersion"] <- as.character(db_list[[species]][[term]][["version"]])
            sim_res@result <- result_sim
            enrich_res <- list(enrich_res, sim_res)
            names(enrich_res) <- paste(group, c(term, paste0(term, "_sim")), sep = "-")
          }
        }

        return(enrich_res)
      } else {
        return(NULL)
      }
    }, BPPARAM = BPPARAM, id = BiocParallel::ipcid()
  )
  nm <- paste(comb$group, comb$term, sep = "-")
  sim_index <- sapply(res_list, function(x) length(x) == 2)
  sim_list <- unlist(res_list[sim_index], recursive = FALSE)
  raw_list <- res_list[!sim_index]
  names(raw_list) <- nm[!sim_index]
  results <- c(raw_list, sim_list)
  results <- results[!sapply(results, is.null)]
  results <- results[intersect(c(nm, paste0(nm, "_sim")), names(results))]
  res_enrichment <- bind_rows(lapply(results, function(x) x@result))
  rownames(res_enrichment) <- NULL

  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "Enrichment done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

  return(list(enrichment = res_enrichment, results = results, geneMap = geneMap, input = input))
}

#' RunSlingshot
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- RunSlingshot(pancreas1k, group.by = "SubCellType", reduction = "UMAP")
#'
#' # 3D lineage
#' pancreas1k <- Standard_SCP(pancreas1k)
#' pancreas1k <- RunSlingshot(pancreas1k, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#' @importFrom Seurat AddMetaData as.SingleCellExperiment
#' @importFrom slingshot slingshot slingPseudotime slingBranchID
#' @export
RunSlingshot <- function(srt, group.by, reduction = NULL, start = NULL, end = NULL, prefix = NULL,
                         reverse = FALSE, align_start = FALSE, show_plot = TRUE, lineage_palette = "Dark2", seed = 11, ...) {
  if (missing(group.by)) {
    stop("group.by is missing")
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (is.null(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }
  srt_sub <- srt[, !is.na(srt[[group.by, drop = TRUE]])]

  set.seed(seed)
  sl <- tryCatch(
    {
      slingshot::slingshot(
        data = as.data.frame(srt_sub[[reduction]]@cell.embeddings),
        clusterLabels = as.character(srt_sub[[group.by, drop = TRUE]]),
        start.clus = start, end.clus = end,
        ...
      )
    },
    error = function(error) {
      sce <- as.SingleCellExperiment(srt_sub)
      sce <- slingshot::slingshot(
        data = sce, reducedDim = reduction,
        clusterLabels = as.character(srt_sub[[group.by, drop = TRUE]]),
        start.clus = start, end.clus = end,
        ...
      )
      sl <- sce@colData$slingshot
      return(sl)
    }
  )

  srt@tools[[paste("Slingshot", group.by, reduction, sep = "_")]] <- sl
  df <- as.data.frame(slingshot::slingPseudotime(sl))
  colnames(df) <- paste0(prefix, colnames(df))
  if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
      df <- max(df, na.rm = TRUE) - df
    }
  }
  srt <- AddMetaData(srt, metadata = df)
  srt <- AddMetaData(srt, metadata = slingshot::slingBranchID(sl), col.name = paste0(prefix, "BranchID"))

  if (isTRUE(show_plot)) {
    if (ncol(srt[[reduction]]@cell.embeddings) == 2) {
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 2, type = "lineages", col = "black")
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 3, col = 1:length(slingshot::SlingshotDataSet(sl)@lineages))
      p <- ClassDimPlot(srt, group.by = group.by, reduction = reduction, lineages = colnames(df))
      print(p)
    } else if (ncol(srt[[reduction]]@cell.embeddings) == 3) {
      p <- ClassDimPlot3D(srt, group.by = group.by, reduction = reduction, lineages = colnames(df))
      print(p)
    } else {
      warning("Unable to plot lineages in the specified reduction.", immediate. = TRUE)
    }
  }
  return(srt)
}

RunMoncle2 <- function(srt, features = NULL, assay = DefaultAssay(srt)) {
  library(monocle)
  expr_matrix <- as(as.matrix(srt_epi@assays$RNA@counts), "sparseMatrix")
  p_data <- srt_epi@meta.data
  f_data <- data.frame(gene_short_name = row.names(srt_epi), row.names = row.names(srt_epi))
  pd <- new("AnnotatedDataFrame", data = p_data)
  fd <- new("AnnotatedDataFrame", data = f_data)
  cds <- newCellDataSet(expr_matrix,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.1,
    expressionFamily = negbinomial.size()
  )
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds, min_expr = 0.1)
  express_genes <- VariableFeatures(srt_epi)
  cds <- setOrderingFilter(cds, express_genes)
  plot_ordering_genes(cds)
  cds <- clusterCells(cds, num_clusters = 2)
  cds <- reduceDimension(cds,
    max_components = 3,
    method = "DDRTree", norm_method = "log"
  )
  cds <- orderCells(cds)
  plot_cell_trajectory(cds, color_by = c("State"), cell_size = 0.8)
  plot_cell_trajectory(cds, color_by = c("Cluster"), cell_size = 0.8)
  plot_cell_trajectory(cds, color_by = c("Pseudotime"), cell_size = 0.8)
  cds <- orderCells(cds, root_state = 2)
  srt <- as.Seurat(cds)
  srt_epi[["DDRTree"]] <- srt[["DDRTree"]]
  srt_epi[["Cluster"]] <- srt[["Cluster"]]
  srt_epi[["Pseudotime"]] <- srt[["Pseudotime"]]
  srt_epi[["State"]] <- srt[["State"]]

  srt_epi@tools$Monocle2 <- list(
    reducedDimS = cds@reducedDimS,
    reducedDimW = cds@reducedDimW,
    reducedDimK = cds@reducedDimK
  )
  srt_epi[["State"]] <- sapply(srt_epi[["State", drop = TRUE]], function(x) {
    switch(as.character(x),
      "1" = 2,
      "2" = 1,
      "3" = 3,
      "4" = 4,
      "5" = 5
    )
  })
  ClassDimPlot(srt_epi, "State", "DDRTree")
  ClassDimPlot(srt_epi, "State", "UncorrectedUMAP2D")
  ExpDimPlot(srt_epi, "Pseudotime", "UncorrectedUMAP2D")
}

RunMoncle3 <- function(srt, group.by, reduction = NULL, start = NULL, end = NULL,
                       seed = 11, plot = TRUE, ...) {
  check_R("slingshot")
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  set.seed(seed)
  sl <- slingshot::slingshot(srt[[reduction]]@cell.embeddings,
    clusterLabels = srt[[group.by, drop = TRUE]],
    start.clus = start, end.clus = end,
    ...
  )

  srt@tools[[paste("Slingshot", group.by, reduction, sep = "_")]] <- sl
  sl@assays@data$pseudotime

  srt <- AddMetaData(srt, metadata = as.data.frame(slingshot::slingPseudotime(sl)))
  srt <- AddMetaData(srt, metadata = slingshot::slingBranchID(sl), col.name = "BranchID")
  if (isTRUE(plot)) {
    # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
    # lines(slingshot::SlingshotDataSet(sl), lwd = 3, type = "l")
    plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
    lines(slingshot::SlingshotDataSet(sl), lwd = 3, col = 1:length(slingshot::SlingshotDataSet(sl)@lineages))
  }
  return(srt)
}

#' RunDynamicFeatures
#'
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- RunSlingshot(pancreas1k, group.by = "SubCellType", reduction = "UMAP")
#' pancreas1k <- RunDynamicFeatures(pancreas1k, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#' ht_result <- DynamicHeatmap(
#'   srt = pancreas1k,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_cluster = 6, reverse_ht = "Lineage1",
#'   height = 5, width = 7, use_raster = FALSE
#' )
#' ht_result$plot
#'
#' p <- DynamicPlot(
#'   srt = pancreas1k,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#' p
#' @export
RunDynamicFeatures <- function(srt, lineages, features = NULL, suffix = lineages,
                               n_candidates = 1000, minfreq = 5,
                               family = NULL, # c("nb", "gaussian", "poisson", "binomial"),
                               slot = "counts", assay = "RNA", libsize = NULL,
                               BPPARAM = BiocParallel::bpparam(), progressbar = TRUE, seed = 11) {
  set.seed(seed)
  if ("progressbar" %in% names(BPPARAM)) {
    BPPARAM[["progressbar"]] <- progressbar
  }

  time_start <- Sys.time()
  cat(paste0("[", time_start, "] ", "Start RunDynamicFeatures\n"))

  check_R(c("mgcv"))
  meta <- c()
  gene <- c()
  if (!is.null(features)) {
    gene <- features[features %in% rownames(srt[[assay]])]
    meta <- features[features %in% colnames(srt@meta.data)]
    isnum <- sapply(srt@meta.data[, meta], is.numeric)
    if (!all(isnum)) {
      warning(paste0(meta[!isnum], collapse = ","), " is not numeric and will be dropped.", immediate. = TRUE)
      meta <- meta[isnum]
    }
    features <- c(gene, meta)
    if (length(features) == 0) {
      stop("No feature found in the srt object.")
    }
  }

  Y <- GetAssayData(srt, slot = slot, assay = assay)
  if (is.null(libsize)) {
    status <- check_DataType(srt, assay = assay, slot = "counts")
    if (status != "raw_counts") {
      Y_libsize <- setNames(rep(1, ncol(srt)), colnames(srt))
    } else {
      Y_libsize <- colSums(GetAssayData(srt, slot = "counts", assay = assay))
    }
  } else {
    if (length(libsize) == 1) {
      Y_libsize <- setNames(rep(libsize, ncol(srt)), colnames(srt))
    } else if (length(libsize) == ncol(srt)) {
      Y_libsize <- setNames(libsize, colnames(srt))
    } else {
      stop("libsize must be length of 1 or the number of cells.")
    }
  }

  if (length(meta) > 0) {
    Y <- rbind(Y, t(srt@meta.data[, meta]))
  }

  features_list <- c()
  srt_sub_list <- list()
  for (l in lineages) {
    srt_sub <- subset(srt, cell = rownames(na.omit(srt[[l]])))
    if (is.null(features)) {
      if (is.null(n_candidates)) {
        stop("'features' or 'n_candidates' must provided at least one.")
      }
      HVF <- VariableFeatures(FindVariableFeatures(srt_sub, nfeatures = n_candidates))
      HVF_counts <- srt_sub[[assay]]@counts[HVF, , drop = FALSE]
      # freq <- aggregate(HVF_counts@x,
      #   by = list(HVF_counts@i),
      #   FUN = function(x) length(unique(x))
      # )
      # sum(freq$x>=(minfreq-1))
      # sum(apply(HVF_counts, 1, function(x) {
      #   length(unique(x))
      # }) >= minfreq)
      HVF <- HVF[apply(HVF_counts, 1, function(x) {
        length(unique(x))
      }) >= minfreq]
      features_list[[l]] <- HVF
    } else {
      features_list[[l]] <- features
    }
    srt_sub_list[[l]] <- srt_sub
  }
  features <- unique(unlist(features_list))
  gene <- features[features %in% rownames(srt[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  message("Number of candidate features(union): ", length(features))

  if (slot == "counts") {
    gene_status <- status
  }
  gene_status <- status <- check_DataType(srt, assay = assay, slot = slot)
  meta_status <- sapply(meta, function(x) {
    check_DataType(data = srt[[x]])
  })
  if (is.null(family)) {
    family <- rep("gaussian", length(features))
    names(family) <- features
    family[names(meta_status)[meta_status == "raw_counts"]] <- "nb"
    if (gene_status == "raw_counts") {
      family[gene] <- "nb"
    }
  } else {
    if (length(family) == 1) {
      family <- rep(family, length(features))
      names(family) <- features
    }
    if (length(family) != length(features)) {
      stop("'family' must be 1 or a vector of the same length as the feature.")
    }
  }

  for (i in 1:length(lineages)) {
    l <- lineages[i]
    srt_sub <- srt_sub_list[[l]]
    t <- na.omit(srt_sub[[l, drop = TRUE]])
    t_ordered <- t[order(t)]
    Y_ordered <- as.matrix(Y[features, names(t_ordered), drop = FALSE])
    l_libsize <- Y_libsize[names(t_ordered)]
    raw_matrix <- as.matrix(cbind(data.frame(pseudotime = t_ordered), t(Y_ordered)))

    # df <- data.frame(x = rowMeans(Y_ordered), y = MatrixGenerics::rowVars(Y_ordered))
    # p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
    #   geom_point() +
    #   geom_abline(slope = 1, intercept = 0, color = "red") +
    #   labs(title = l, subtitle = "Var VS Mean") +
    #   theme_scp()
    # print(p)

    message("Calculate dynamic features for ", l, "...")
    system.time({
      gam_out <- BiocParallel::bplapply(1:nrow(Y_ordered), function(n, Y_ordered, t_ordered, l_libsize, family) {
        feature_nm <- rownames(Y_ordered)[n]
        family_current <- family[feature_nm]
        if (min(Y_ordered[feature_nm, ]) < 0 && family_current %in% c("nb", "poisson", "binomial")) {
          warning("Negative values detected. Replace family with 'gaussian' for the feature: ", feature_nm, immediate. = TRUE)
          family_use <- "gaussian"
        } else {
          family_use <- family_current
        }
        if (slot == "counts" && family_use != "gaussian" && !feature_nm %in% meta) {
          l_libsize <- l_libsize
        } else {
          l_libsize <- rep(median(Y_libsize), ncol(Y_ordered))
        }
        sizefactror <- median(Y_libsize) / l_libsize
        mod <- mgcv::gam(y ~ s(x, bs = "cs") + offset(log(l_libsize)),
          family = family_use,
          data = data.frame(y = Y_ordered[feature_nm, ], x = t_ordered, l_libsize = l_libsize)
        )
        pre <- predict(mod, type = "link", se.fit = TRUE)
        upr <- pre$fit + (2 * pre$se.fit)
        lwr <- pre$fit - (2 * pre$se.fit)
        upr <- mod$family$linkinv(upr)
        lwr <- mod$family$linkinv(lwr)
        res <- summary(mod)
        fitted <- fitted(mod)
        pvalue <- res$s.table[[4]]
        dev.expl <- res$dev.expl
        r.sq <- res$r.sq
        fitted.values <- fitted * sizefactror
        upr.values <- upr * sizefactror
        lwr.values <- lwr * sizefactror
        exp_ncells <- sum(Y_ordered[feature_nm, ] > min(Y_ordered[feature_nm, ]), na.rm = TRUE)
        peak_time <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
        valley_time <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])

        # ggplot(data = data.frame(
        #   x = t_ordered,
        #   raw = FetchData(srt_sub, vars = feature_nm, slot = "counts")[names(t_ordered), feature_nm, drop = TRUE],
        #   fitted = fitted.values,
        #   upr.values = upr.values,
        #   lwr.values = lwr.values,
        #   l_libsize = l_libsize
        # )) +
        #   geom_point(aes(x = x, y = raw), color = "black", size = 0.5) +
        #   geom_point(aes(x = x, y = fitted), color = "red", size = 0.5) +
        #   geom_path(aes(x = x, y = upr.values), color = "blue") +
        #   geom_path(aes(x = x, y = lwr.values), color = "green")

        # a <- data.frame(x = t_ordered, y = Y_ordered[feature_nm, ])
        # qplot(a$x, a$y)
        # length(unique(a$y) > 5)
        # a <- a[a$y > 0, ]
        # qplot(a$x, a$y)

        return(list(
          features = feature_nm, exp_ncells = exp_ncells,
          r.sq = r.sq, dev.expl = dev.expl,
          peak_time = peak_time, valley_time = valley_time,
          pvalue = pvalue, fitted.values = fitted.values,
          upr.values = upr.values, lwr.values = lwr.values
        ))
      }, BPPARAM = BPPARAM, Y_ordered = Y_ordered, t_ordered = t_ordered, l_libsize = l_libsize, family = family)
    })
    fitted_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["fitted.values"]]))
    colnames(fitted_matrix) <- rownames(Y_ordered)
    fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)

    upr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["upr.values"]]))
    colnames(upr_matrix) <- rownames(Y_ordered)
    upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)

    lwr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["lwr.values"]]))
    colnames(lwr_matrix) <- rownames(Y_ordered)
    lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)

    DynamicFeatures <- as.data.frame(do.call(rbind.data.frame, lapply(gam_out, function(x) x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")])))
    char_var <- c("features")
    numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
    DynamicFeatures[, char_var] <- lapply(DynamicFeatures[, char_var, drop = FALSE], as.character)
    DynamicFeatures[, numb_var] <- lapply(DynamicFeatures[, numb_var, drop = FALSE], as.numeric)
    rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
    DynamicFeatures[, "padjust"] <- p.adjust(DynamicFeatures[, "pvalue", drop = TRUE])

    res <- list(
      DynamicFeatures = DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = fitted_matrix,
      upr_matrix = upr_matrix,
      lwr_matrix = lwr_matrix,
      libsize = l_libsize,
      lineages = l,
      family = family
    )
    srt@tools[[paste0("DynamicFeatures_", suffix[i])]] <- res
  }

  time_end <- Sys.time()
  cat(paste0("[", time_end, "] ", "RunDynamicFeatures done\n"))
  cat("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"), "\n")

  return(srt)
}

#' RunDynamicEnrichment
#'
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas1k")
#'   pancreas1k <- RunSlingshot(pancreas1k, group.by = "SubCellType", reduction = "UMAP")
#'   pancreas1k <- RunDynamicFeatures(pancreas1k, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#'   ht_result <- DynamicHeatmap(
#'     srt = pancreas1k,
#'     lineages = c("Lineage1", "Lineage2"),
#'     cell_annotation = "SubCellType",
#'     n_cluster = 6, reverse_ht = 1, use_raster = FALSE
#'   )
#'   ht_result$plot
#'
#'   pancreas1k <- RunDynamicEnrichment(
#'     srt = pancreas1k,
#'     lineages = c("Lineage1", "Lineage2"),
#'     score_method = "AUCell",
#'     enrichment = "GO_BP",
#'     species = "Mus_musculus"
#'   )
#'   ht_result <- DynamicHeatmap(
#'     srt = pancreas1k, assay = "GO_BP", use_fitted = TRUE,
#'     lineages = c("Lineage1_GO_BP", "Lineage2_GO_BP"),
#'     cell_annotation = "SubCellType",
#'     n_cluster = 6, reverse_ht = 1,
#'     height = 5, width = 7, use_raster = FALSE
#'   )
#'   ht_result$plot
#' }
#' @export
RunDynamicEnrichment <- function(srt, lineages,
                                 score_method = "AUCell", ncore = 1,
                                 slot = "data", assay = "RNA",
                                 min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05,
                                 geneID = NULL, IDtype = "symbol", species = "Homo_sapiens",
                                 enrichment = "GO_BP", db_update = FALSE, Ensembl_version = 103, mirror = NULL,
                                 TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                                 BPPARAM = BiocParallel::bpparam(), progressbar = TRUE, seed = 11) {
  set.seed(seed)
  if ("progressbar" %in% names(BPPARAM)) {
    BPPARAM[["progressbar"]] <- progressbar
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
    feature_union <- c(feature_union, DynamicFeatures[, "features"])
    cell_union <- c(cell_union, rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]))
  }
  feature_union <- unique(feature_union)

  if (is.null(TERM2GENE)) {
    db_list <- PrepareEnrichmentDB(
      species = species, enrichment = enrichment, db_update = db_update,
      db_IDtypes = IDtype, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    enrichment <- "custom"
    db_list <- list()
    db_list[[species]][[enrichment]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[enrichment]][["TERM2NAME"]] <- unique(TERM2NAME)
  }

  for (i in 1:length(enrichment)) {
    term <- enrichment[i]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", IDtype)]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), ]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], ]

    term_use <- unique(TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% feature_union, "Term"])
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, "Term"] %in% term_use, ]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% term_use, ]
    rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% rownames(srt), ]
    feature_list <- split(TERM2GENE_tmp[, IDtype], TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"])
    GSSize <- sapply(feature_list, length)
    feature_list <- feature_list[GSSize >= minGSSize & GSSize <= maxGSSize]

    srt <- CellScoring(
      srt = srt,
      features = feature_list,
      method = score_method,
      classification = FALSE,
      ncore = ncore,
      name = term,
      slot = slot,
      assay = assay,
      new_assay = TRUE
    )
    srt <- RunDynamicFeatures(
      srt = srt,
      lineages = lineages,
      features = rownames(srt[[term]]@counts),
      suffix = paste(lineages, term, sep = "_"),
      assay = term
    )
    # ht_result <- DynamicHeatmap(
    #   srt = srt, assay = term,use_fitted = F,
    #   lineages = c("Lineage1", "Lineage2"),
    #   cell_annotation = "SubCellType",
    #     reverse_ht = 1, use_raster = F
    # )
    # ht_result$plot
  }
  return(srt)
}


#' Convert a seurat object to an anndata object using reticulate
#'
#' @param srt A \code{Seurat} object.
#' @param assay_X Assays to convert as X(main data matrix) in anndata object.
#' @param assay_layers Assays to convert as layers in anndata object.
#' @param slot slot use to convert for the slot of X and layer in anndata object.
#'
#' @return A \code{anndata} object.
#' @examples
#' data("pancreas1k")
#' adata <- srt_to_adata(pancreas1k)
#' adata
#' @importFrom reticulate import
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @export
#'
srt_to_adata <- function(srt,
                         assay_X = "RNA", slot_X = "counts",
                         assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                         convert_tools = FALSE, convert_misc = FALSE) {
  if (class(srt) != "Seurat") {
    stop("'srt' is not a Seurat object.")
  }
  if (length(slot_layers) == 1) {
    slot_layers <- rep(slot_layers, length(assay_layers))
    names(slot_layers) <- assay_layers
  } else if (length(slot_layers) != length(assay_layers)) {
    stop("slot_layers must be one character or the same length with the assay_layers")
  }

  sc <- import("scanpy", convert = FALSE)
  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in 1:ncol(obs)) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }

  var <- srt[[assay_X]]@meta.features
  if (ncol(var) > 0) {
    for (i in 1:ncol(var)) {
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(as.character(var[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }

  adata <- sc$AnnData(
    X = t(GetAssayData(srt, assay = assay_X, slot = slot_X)),
    obs = obs,
    var = cbind(data.frame(features = rownames(srt)), var)
  )
  adata$var_names <- rownames(srt)
  if (length(VariableFeatures(srt, assay = assay_X) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      adata$var <- var[, colnames(var) != "highly_variable"]
    }
    adata$var <- adata$var$join(data.frame(row.names = rownames(srt), highly_variable = rownames(srt) %in% VariableFeatures(srt, assay = assay_X)))
  }

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_X]) {
    if (assay %in% assay_layers) {
      layer_list[[assay]] <- t(GetAssayData(srt, assay = assay, slot = slot_layers[assay]))
    } else {
      message("Assay '", assay, "' is in the srt object but not converted.")
    }
  }
  if (length(layer_list) > 0) {
    adata$layers <- layer_list
  }

  reduction_list <- list()
  for (reduction in names(srt@reductions)) {
    reduction_list[[paste0(reduction)]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  obsp_list <- list()
  for (graph in names(srt@graphs)) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in names(srt@neighbors)) {
    obsp_list[[neighbor]] <- srt[[neighbor]]
  }
  if (length(obsp_list) > 0) {
    adata$obsp <- obsp_list
  }

  uns_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   uns_list[[paste0(reduction,"_stdev")]] <- srt[[reduction]]@stdev
  #   uns_list[[paste0(reduction,"_misc")]] <- srt[[reduction]]@misc
  # }
  # for (image in names(srt@images)) {
  #   uns_list[[image]] <- srt[[image]]
  # }
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@misc[[nm]]
      }
    }
  } else {
    message("'misc' slot in the srt object is not converted.")
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    message("'tools' slot in the srt object is not converted.")
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}

#' Convert an anndata object to a seurat object using reticulate
#' @param adata a connected python anndata object.
#'
#' @examples
#' data("pancreas1k")
#' adata <- srt_to_adata(pancreas1k)
#' adata <- RunPAGA(adata = adata, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP")
#' srt <- adata_to_srt(adata)
#' srt
#' @importFrom Seurat CreateSeuratObject CreateAssayObject CreateDimReducObject AddMetaData
#' @importFrom SeuratObject as.Graph
#' @importFrom reticulate iterate
#' @importFrom Matrix t
#' @export
adata_to_srt <- function(adata) {
  if (!"python.builtin.object" %in% class(adata)) {
    stop("'adata' is not a Python object.")
  }
  x <- t(adata$X)
  rownames(x) <- adata$var_names$values
  colnames(x) <- adata$obs_names$values

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(adata$obs)
  }

  srt <- CreateSeuratObject(counts = x, meta.data = metadata)

  if (length(adata$layers$keys()) > 0) {
    for (k in iterate(adata$layers$keys())) {
      layer <- t(adata$layers[[k]])
      rownames(layer) <- adata$var_names$values
      colnames(layer) <- adata$obs_names$values
      srt[[k]] <- CreateAssayObject(counts = layer)
    }
  }
  if (length(iterate(adata$obsm$keys())) > 0) {
    for (k in iterate(adata$obsm$keys())) {
      obsm <- adata$obsm[[k]]
      colnames(obsm) <- paste0(k, "_", 1:ncol(obsm))
      rownames(obsm) <- adata$obs_names$values
      srt[[k]] <- CreateDimReducObject(embeddings = obsm, assay = "RNA", key = paste0(k, "_"))
    }
  }
  if (length(iterate(adata$obsp$keys())) > 0) {
    for (k in iterate(adata$obsp$keys())) {
      obsp <- adata$obsp[[k]]
      colnames(obsp) <- adata$obs_names$values
      rownames(obsp) <- adata$obs_names$values
      obsp <- as.Graph(obsp[1:nrow(obsp), ])
      DefaultAssay(object = obsp) <- "RNA"
      srt[[k]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- AddMetaData(srt[["RNA"]], metadata = as.data.frame(adata$var))
  }
  if (length(iterate(adata$varm$keys())) > 0) {
    for (k in iterate(adata$varm$keys())) {
      varm <- adata$varm[[k]]
      colnames(varm) <- adata$var_names$values
      rownames(varm) <- adata$var_names$values
      srt[["RNA"]]@misc[["reductions"]][[k]] <- varm
    }
  }
  if (length(iterate(adata$varp$keys())) > 0) {
    for (k in iterate(adata$varp$keys())) {
      varp <- adata$varp[[k]]
      colnames(varp) <- adata$var_names$values
      rownames(varp) <- adata$var_names$values
      srt[["RNA"]]@misc[["graphs"]][[k]] <- varp
    }
  }

  if (length(iterate(adata$uns$keys())) > 0) {
    for (k in iterate(adata$uns$keys())) {
      uns <- adata$uns[[k]]
      uns <- check_python_element(uns)
      if (!"python.builtin.object" %in% class(uns)) {
        srt@misc[[k]] <- uns
      } else {
        warning("'uns: ", k, "' will not be converted. It is a Python object connected with the Python session. You may need to convert manually.", immediate. = TRUE)
      }
    }
  }
  return(srt)
}

maxDepth <- function(x, depth = 0) {
  if (is.list(x)) {
    return(max(sapply(x, maxDepth, depth + 1)))
  } else {
    return(depth)
  }
}

#' @importFrom reticulate py_to_r
check_python_element <- function(x, depth = maxDepth(x)) {
  if (depth == 0 || !is.list(x)) {
    if ("python.builtin.object" %in% class(x)) {
      x_r <- tryCatch(py_to_r(x), error = identity)
      if (inherits(x_r, "error")) {
        return(x)
      } else {
        return(x_r)
      }
    } else {
      return(x)
    }
  } else {
    raw_depth <- maxDepth(x)
    x <- lapply(x, function(element) {
      if ("python.builtin.object" %in% class(element)) {
        element_r <- tryCatch(py_to_r(element), error = identity)
        if (inherits(element_r, "error")) {
          return(element)
        } else {
          return(element_r)
        }
      } else {
        return(element)
      }
    })
    cur_depth <- maxDepth(x)
    if (cur_depth > raw_depth) {
      depth <- depth + 1
    }
    x_checked <- lapply(x, check_python_element, depth - 1)
    return(x_checked)
  }
}

#' Run PAGA analysis
#'
#' @inheritParams RunSCVELO
#' @param use_rna_velocity use_rna_velocity.
#' @param threshold edge threshold.
#' @param point_size point_size.
#' @param vkey Key for velocity. Default is "stochastic".
#' @param paga_layout Plotting layout that computes positions. See See \href{https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.paga.html}{layout} param in scanpy.pl.paga function.
#'
#' @return A \code{anndata} object.
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- RunPAGA(srt = pancreas1k, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
#' @export
#'
RunPAGA <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                    adata = NULL, h5ad = NULL, group_by = NULL,
                    liner_reduction = NULL, nonliner_reduction = NULL, basis = NULL,
                    n_pcs = 30, n_neighbors = 30, use_rna_velocity = FALSE, vkey = "stochastic",
                    embedded_with_PAGA = FALSE, paga_layout = "fr", threshold = 0.1, point_size = 20,
                    show_plot = TRUE, dpi = 300, save = FALSE, dirpath = "./", fileprefix = "",
                    return_seurat = FALSE) {
  if (all(is.null(srt), is.null(adata), is.null(h5ad))) {
    stop("One of 'srt', 'adata' or 'h5ad' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(liner_reduction) & is.null(nonliner_reduction)) {
    stop("'liner_reduction' or 'nonliner_reduction' must be provided at least one.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
    args[["adata"]] <- adata
  }
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat")]

  adata <- do.call(
    what = PAGA,
    args = args
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(paga)|(distances)|(connectivities)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}


#' Run scVelo analysis
#'
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells.
#' This function performs scVelo workflow in R by reticulate.
#'
#' @param srt A \code{Seurat} object.
#' @param adata An \code{anndata} object. Generally created through \code{\link{srt_to_adata}}
#' @param h5ad h5ad file path.
#' @param group_by group_by.
#' @param liner_reduction liner_reduction.
#' @param nonliner_reduction nonliner_reduction.
#' @param basis basis.
#' @param mode mode.
#' @param fitting_by fitting_by.
#' @param n_jobs n_jobs.
#' @param min_shared_counts min_shared_counts.
#' @param n_pcs n_pcs.
#' @param n_neighbors n_neighbors.
#' @param approx approx.
#' @param stream_smooth stream_smooth.
#' @param stream_density stream_density.
#' @param arrow_density arrow_density.
#' @param arrow_length arrow_length.
#' @param arrow_size arrow_size.
#' @param velocity_with_noise velocity_with_noise.
#' @param diff_kinetics diff_kinetics.
#' @param calculate_velocity_genes calculate_velocity_genes.
#' @param s_genes s_genes.
#' @param g2m_genes g2m_genes.
#' @param save save.
#' @param dirpath dirpath.
#' @param fileprefix fileprefix.
#' @param dpi dpi.
#'
#' @return A \code{anndata} object.
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- RunSCVELO(srt = pancreas1k, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE)
#' @export
#'
RunSCVELO <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                      adata = NULL, h5ad = NULL, group_by = NULL, n_jobs = 1,
                      liner_reduction = NULL, nonliner_reduction = NULL, basis = NULL,
                      mode = "stochastic", fitting_by = "stochastic",
                      magic_impute = FALSE, knn = 5, t = 2,
                      min_shared_counts = 30, n_pcs = 30, n_neighbors = 30, approx = TRUE,
                      stream_smooth = NULL, stream_density = 2,
                      arrow_length = 5, arrow_size = 5, arrow_density = 0.5,
                      denoise = FALSE, denoise_topn = 3, kinetics = FALSE, kinetics_topn = 100,
                      calculate_velocity_genes = FALSE, s_genes = NULL, g2m_genes = NULL,
                      show_plot = TRUE, dpi = 300, save = FALSE, dirpath = "./", fileprefix = "",
                      return_seurat = FALSE) {
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata), is.null(h5ad))) {
    stop("One of 'srt', 'adata' or 'h5ad' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(liner_reduction) & is.null(nonliner_reduction)) {
    stop("'liner_reduction' or 'nonliner_reduction' must be provided at least one.")
  }
  mode <- as.list(mode)
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
    args[["adata"]] <- adata
  }
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat")]

  adata <- do.call(
    what = SCVELO,
    args = args
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = paste0("(Ms)|(Mu)|(velocity)|(", mode, ")"), overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run Palantir analysis
#' @inheritParams RunSCVELO
#' @param point_size point_size.
#'
#' @return A \code{anndata} object.
#'
#' @examples
#' data("pancreas1k")
#' pancreas1k <- RunPalantir(
#'   srt = pancreas1k, group_by = "SubCellType", liner_reduction = "PCA", nonliner_reduction = "UMAP",
#'   early_group = "Ductal", terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon"), return_seurat = TRUE
#' )
#' @export
#'
RunPalantir <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, h5ad = NULL, group_by = NULL,
                        liner_reduction = NULL, nonliner_reduction = NULL, basis = NULL,
                        n_pcs = 30, n_neighbors = 30, dm_n_components = 10, dm_alpha = 0, dm_n_eigs = NULL,
                        early_group = NULL, terminal_groups = NULL, early_cell = NULL, terminal_cells = NULL,
                        num_waypoints = 1200, scale_components = TRUE, use_early_cell_as_start = FALSE,
                        max_iterations = 25, n_jobs = 8, point_size = 20,
                        show_plot = TRUE, dpi = 300, save = FALSE, dirpath = "./", fileprefix = "",
                        return_seurat = FALSE) {
  if (all(is.null(srt), is.null(adata), is.null(h5ad))) {
    stop("One of 'srt', 'adata' or 'h5ad' must be provided.")
  }
  if (is.null(group_by) & any(!is.null(early_group), !is.null(terminal_groups))) {
    stop("'group_by' must be provided when early_group or terminal_groups provided.")
  }
  if (is.null(liner_reduction) & is.null(nonliner_reduction)) {
    stop("'liner_reduction' or 'nonliner_reduction' must be provided at least one.")
  }
  if (is.null(early_cell) & is.null(early_group)) {
    stop("'early_cell' or 'early_group' must be provided.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
    args[["adata"]] <- adata
  }
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat")]

  adata <- do.call(
    what = Palantir,
    args = args
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(dm_kernel)|(palantir)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

RunCellRank <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, h5ad = NULL, group_by = NULL, n_jobs = 1,
                        liner_reduction = NULL, nonliner_reduction = NULL, basis = NULL,
                        mode = "stochastic", fitting_by = "stochastic",
                        magic_impute = FALSE, knn = 5, t = 2,
                        min_shared_counts = 30, n_pcs = 30, n_neighbors = 30, approx = TRUE,
                        stream_smooth = NULL, stream_density = 2,
                        arrow_size = 5, arrow_length = 5, arrow_density = 0.5,
                        s_genes = NULL, g2m_genes = NULL, calculate_velocity_genes = FALSE,
                        denoise = FALSE, kinetics = FALSE, axis = "equal",
                        show_plot = TRUE, dpi = 300, save = FALSE, dirpath = "./", fileprefix = "",
                        return_seurat = FALSE) {
  check_Python(c("pandas", "numpy", "scvelo", "cellrank[external]", "cython", "wot", "statot", "POT"))
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata), is.null(h5ad))) {
    stop("One of 'srt', 'adata' or 'h5ad' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(liner_reduction) & is.null(nonliner_reduction)) {
    stop("'liner_reduction' or 'nonliner_reduction' must be provided at least one.")
  }
  mode <- as.list(mode)
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
    args[["adata"]] <- adata
  }
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat")]

  adata <- do.call(
    what = CellRank,
    args = args
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      return(SrtAppend(srt_raw = srt, srt_append = srt_out))
    }
  } else {
    return(adata)
  }
}

RunDynamo <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                      adata = NULL, h5ad = NULL, group_by = NULL,
                      liner_reduction = NULL, nonliner_reduction = NULL, basis = NULL,
                      n_pcs = 30, n_neighbors = 30, dm_n_components = 10, dm_alpha = 0, dm_n_eigs = NULL,
                      early_group = NULL, terminal_groups = NULL, early_cell = NULL, terminal_cells = NULL,
                      num_waypoints = 1200, scale_components = TRUE, use_early_cell_as_start = FALSE,
                      max_iterations = 25, n_jobs = 1, point_size = 20,
                      show_plot = TRUE, dpi = 300, save = FALSE, dirpath = "./", fileprefix = "",
                      return_seurat = FALSE) {
  check_Python(c("matplotlib", "python-igraph", "pandas", "numpy", "dynamo-release"))
  if (all(is.null(srt), is.null(adata), is.null(h5ad))) {
    stop("One of 'srt', 'adata' or 'h5ad' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(liner_reduction) & is.null(nonliner_reduction)) {
    stop("'liner_reduction' or 'nonliner_reduction' must be provided at least one.")
  }
  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      y <- ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      y <- x
    }
    return(y)
  })
  call.envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  if (!is.null(srt)) {
    adata <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
    args[["adata"]] <- adata
  }
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat")]

  adata <- do.call(
    what = Dynamo,
    args = args
  )

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      return(SrtAppend(srt_raw = srt, srt_append = srt_out))
    }
  } else {
    return(adata)
  }
}
