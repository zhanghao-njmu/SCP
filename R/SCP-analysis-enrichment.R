RunEnrichment <- function(srt = NULL, group_by = NULL, test.use = "wilcox", DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
                          geneID = NULL, geneID_groups = NULL, geneID_exclude = NULL, IDtype = "symbol", result_IDtype = "symbol", species = "Homo_sapiens",
                          db = "GO_BP", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = TRUE, Ensembl_version = 103, mirror = NULL,
                          TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
                          GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                          BPPARAM = BiocParallel::bpparam(), seed = 11) {
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start Enrichment"))
  message("Workers: ", bpworkers(BPPARAM))

  use_srt <- FALSE
  if (is.null(geneID)) {
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
    de_df <- de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), , drop = FALSE]
    rownames(de_df) <- seq_len(nrow(de_df))

    geneID <- de_df[["gene"]]
    geneID_groups <- de_df[["group1"]]
    use_srt <- TRUE
  }

  if (is.null(geneID_groups)) {
    geneID_groups <- rep(" ", length(geneID))
  }
  if (!is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  geneID_groups <- factor(geneID_groups, levels = levels(geneID_groups)[levels(geneID_groups) %in% geneID_groups])
  if (length(geneID_groups) != length(geneID)) {
    stop("length(geneID_groups)!=length(geneID)")
  }
  names(geneID_groups) <- geneID
  input <- data.frame(geneID = geneID, geneID_groups = geneID_groups)
  input <- input[!geneID %in% geneID_exclude, , drop = FALSE]

  if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species, db = db, db_update = db_update, db_version = db_version,
      db_IDtypes = IDtype, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    db <- "custom"
    db_list <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- "custom"
  }
  if (isTRUE(db_combine)) {
    message("Create 'Combined' database ...")
    TERM2GENE <- do.call(rbind, lapply(db_list[[species]], function(x) x[["TERM2GENE"]][, c("Term", IDtype)]))
    TERM2NAME <- do.call(rbind, lapply(names(db_list[[species]]), function(x) {
      db_list[[species]][[x]][["TERM2NAME"]][["Name"]] <- paste0(db_list[[species]][[x]][["TERM2NAME"]][["Name"]], " [", x, "]")
      db_list[[species]][[x]][["TERM2NAME"]][, c("Term", "Name")]
    }))
    version <- unlist(lapply(db_list[[species]], function(x) as.character(x[["version"]])))
    version <- paste0(names(version), ":", version, collapse = ";")
    db <- "Combined"
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- unique(version)
  }

  if (length(unique(c(IDtype, result_IDtype))) != 1) {
    res <- GeneConvert(
      geneID = unique(geneID),
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = result_IDtype,
      species_from = species,
      species_to = species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
    geneMap <- res$geneID_collapse
    colnames(geneMap)[colnames(geneMap) == "from_geneID"] <- IDtype
  } else {
    geneMap <- data.frame(IDtype = unique(geneID), row.names = unique(geneID))
    colnames(geneMap)[1] <- IDtype
  }

  input[[IDtype]] <- geneMap[as.character(input$geneID), IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- unnest(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  message("Permform enrichment...")
  comb <- expand.grid(group = levels(geneID_groups), term = db, stringsAsFactors = FALSE)

  res_list <- bplapply(seq_len(nrow(comb)), function(i, id) {
    group <- comb[i, "group"]
    term <- comb[i, "term"]
    gene <- input[input$geneID_groups == group, IDtype]
    gene_mapid <- input[input$geneID_groups == group, result_IDtype]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", IDtype)]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], , drop = FALSE]
    enrich_res <- enricher(
      gene = gene,
      minGSSize = ifelse(term %in% unlimited_db, 1, minGSSize),
      maxGSSize = ifelse(term %in% unlimited_db, Inf, maxGSSize),
      pAdjustMethod = "BH",
      pvalueCutoff = Inf,
      qvalueCutoff = Inf,
      universe = NULL,
      TERM2GENE = TERM2GENE_tmp,
      TERM2NAME = TERM2NAME_tmp
    )

    if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
      result <- enrich_res@result
      result[["Groups"]] <- group
      result[["Database"]] <- term
      result[["Version"]] <- as.character(db_list[[species]][[term]][["version"]])
      IDlist <- strsplit(result$geneID, split = "/")
      result$geneID <- unlist(lapply(IDlist, function(x) {
        x_result <- NULL
        for (i in x) {
          if (i %in% geneMap[[IDtype]]) {
            x_result <- c(x_result, unique(geneMap[geneMap[[IDtype]] == i, result_IDtype]))
          } else {
            x_result <- c(x_result, i)
          }
        }
        return(paste0(x_result, collapse = "/"))
      }))
      enrich_res@result <- result
      enrich_res@gene2Symbol <- as.character(gene_mapid)

      if (isTRUE(GO_simplify) && term %in% c("GO", "GO_BP", "GO_CC", "GO_MF")) {
        sim_res <- enrich_res
        if (term == "GO") {
          sim_res@result[["ONTOLOGY"]] <- setNames(TERM2NAME_tmp[["ONTOLOGY"]], TERM2NAME_tmp[["Term"]])[sim_res@result[["ID"]]]
          sim_res@ontology <- "GOALL"
        } else {
          sim_res@ontology <- gsub(pattern = "GO_", replacement = "", x = term)
        }
        nterm_simplify <- sum(with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))))
        if (nterm_simplify <= 1) {
          warning(group, "|", term, " has no term to simplify.", immediate. = TRUE)
        } else {
          sim_res@result <- sim_res@result[with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))), , drop = FALSE]
          semData <- db_list[[species]][[term]][["semData"]]
          ipclock(id)
          sim_res <- simplify(sim_res,
            measure = simplify_method,
            cutoff = simplify_similarityCutoff,
            semData = semData
          )
          ipcunlock(id)
          result_sim <- sim_res@result
          result_sim[["Groups"]] <- group
          result_sim[["Database"]] <- paste0(term, "_sim")
          result_sim[["Version"]] <- as.character(db_list[[species]][[term]][["version"]])
          result_sim[["ONTOLOGY"]] <- NULL
          sim_res@result <- result_sim
          enrich_res <- list(enrich_res, sim_res)
          names(enrich_res) <- paste(group, c(term, paste0(term, "_sim")), sep = "-")
        }
      }
    } else {
      enrich_res <- NULL
    }
    return(enrich_res)
  }, BPPARAM = BPPARAM, id = ipcid())

  nm <- paste(comb$group, comb$term, sep = "-")
  sim_index <- sapply(res_list, function(x) length(x) == 2)
  sim_list <- unlist(res_list[sim_index], recursive = FALSE)
  raw_list <- res_list[!sim_index]
  names(raw_list) <- nm[!sim_index]
  results <- c(raw_list, sim_list)
  results <- results[!sapply(results, is.null)]
  results <- results[intersect(c(nm, paste0(nm, "_sim")), names(results))]
  enrichment <- do.call(rbind, lapply(results, function(x) x@result))
  rownames(enrichment) <- NULL

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "Enrichment done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  res <- list(enrichment = enrichment, results = results, geneMap = geneMap, input = input)
  if (isTRUE(use_srt)) {
    res[["DE_threshold"]] <- DE_threshold
    srt@tools[[paste("Enrichment", group_by, test.use, sep = "_")]] <- res
    return(srt)
  } else {
    return(res)
  }
}

#' Perform the enrichment analysis (GSEA) on the genes
#'
#' @inheritParams RunEnrichment
#' @param geneScore A numeric vector that specifies the gene scores, for example, the log2(fold change) values of gene expression.
#' @param scoreType This parameter defines the GSEA score type. Possible options are "std", "pos", "neg". By default ("std") the enrichment score is computed as in the original GSEA. The "pos" and "neg" score types are intended to be used for one-tailed tests (i.e. when one is interested only in positive ("pos") or negateive ("neg") enrichment).
#' @returns
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#'
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#'
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item{\code{enrichment}:}{ A data.frame containing all enrichment results.}
#'  \item{\code{results}:}{ A list of \code{gseaResult} objects from the DOSE package.}
#'  \item{\code{geneMap}:}{ A data.frame containing the ID mapping table for input gene IDs.}
#'  \item{\code{input}:}{ A data.frame containing the input gene IDs and gene ID groups.}
#'  \item{\code{DE_threshold}:}{ A specific threshold for differential expression analysis (only returned if input is a Seurat object).}
#' }
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}} \code{\link{GSEAPlot}} \code{\link{RunEnrichment}} \code{\link{EnrichmentPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
#' pancreas_sub <- RunGSEA(pancreas_sub,
#'   group_by = "CellType", DE_threshold = "p_val_adj < 0.05",
#'   scoreType = "std", db = "GO_BP", species = "Mus_musculus"
#' )
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", plot_type = "comparison")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Ductal", id_use = "GO:0006412")
#' GSEAPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", group_use = "Endocrine", id_use = c("GO:0046903", "GO:0015031", "GO:0007600"))
#'
#' # Remove redundant GO terms
#' pancreas_sub <- RunGSEA(srt = pancreas_sub, group_by = "CellType", db = "GO_BP", GO_simplify = TRUE, species = "Mus_musculus")
#' GSEAPlot(pancreas_sub, db = "GO_BP_sim", group_by = "CellType", plot_type = "comparison")
#'
#' # Use a combined database
#' pancreas_sub <- RunGSEA(
#'   srt = pancreas_sub, group_by = "CellType",
#'   db = c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' GSEAPlot(pancreas_sub, db = "Combined", group_by = "CellType", plot_type = "comparison")
#'
#' # Or use "geneID", "geneScore" and "geneID_groups" as input to run GSEA
#' de_df <- dplyr::filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05)
#' gsea_out <- RunGSEA(geneID = de_df[["gene"]], geneScore = de_df[["avg_log2FC"]], geneID_groups = de_df[["group1"]], db = "GO_BP", species = "Mus_musculus")
#' GSEAPlot(res = gsea_out, db = "GO_BP", plot_type = "comparison")
#'
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers ipcid ipclock ipcunlock
#' @importFrom clusterProfiler GSEA simplify
#' @export
RunGSEA <- function(srt = NULL, group_by = NULL, test.use = "wilcox", DE_threshold = "p_val_adj < 0.05", scoreType = "std",
                    geneID = NULL, geneScore = NULL, geneID_groups = NULL, geneID_exclude = NULL, IDtype = "symbol", result_IDtype = "symbol", species = "Homo_sapiens",
                    db = "GO_BP", db_update = FALSE, db_version = "latest", db_combine = FALSE, convert_species = TRUE, Ensembl_version = 103, mirror = NULL,
                    TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500, unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
                    GO_simplify = FALSE, GO_simplify_cutoff = "p.adjust < 0.05", simplify_method = "Wang", simplify_similarityCutoff = 0.7,
                    BPPARAM = BiocParallel::bpparam(), seed = 11) {
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start GSEA"))
  message("Workers: ", bpworkers(BPPARAM))

  use_srt <- FALSE
  if (is.null(geneID)) {
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
    de_df <- de_df[with(de_df, eval(rlang::parse_expr(DE_threshold))), , drop = FALSE]
    rownames(de_df) <- seq_len(nrow(de_df))

    geneID <- de_df[["gene"]]
    geneScore <- de_df[["avg_log2FC"]]
    geneID_groups <- de_df[["group1"]]
    use_srt <- TRUE
  }

  if (is.null(geneID_groups)) {
    geneID_groups <- rep(" ", length(geneID))
  }
  if (!is.factor(geneID_groups)) {
    geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
  }
  geneID_groups <- factor(geneID_groups, levels = levels(geneID_groups)[levels(geneID_groups) %in% geneID_groups])
  if (length(geneID_groups) != length(geneID)) {
    stop("length(geneID_groups)!=length(geneID)")
  }
  if (length(geneScore) != length(geneID)) {
    stop("geneScore must be the same length with geneID")
  }
  if (all(geneScore > 0) && scoreType != "pos") {
    scoreType <- "pos"
    warning("All values in the geneScore are greater than zero. Set scoreType = 'pos'.", immediate. = TRUE)
  }
  if (all(geneScore < 0) && scoreType != "neg") {
    scoreType <- "neg"
    warning("All values in the geneScore are less than zero. Set scoreType = 'neg'.", immediate. = TRUE)
  }

  input <- data.frame(geneID = geneID, geneScore = geneScore, geneID_groups = geneID_groups)
  input <- input[!geneID %in% geneID_exclude, , drop = FALSE]

  na_index <- which(is.na(geneScore))
  if (length(na_index) > 0) {
    message("Ignore ", length(na_index), " NA geneScore")
    input <- input[-na_index, , drop = FALSE]
  }
  input[is.infinite(input$geneScore) & input$geneScore < 0, "geneScore"] <- min(input[!is.infinite(input$geneScore), "geneScore"])
  input[is.infinite(input$geneScore) & input$geneScore > 0, "geneScore"] <- max(input[!is.infinite(input$geneScore), "geneScore"])

  geneID <- input$geneID
  geneScore <- input$geneScore
  geneID_groups <- input$geneID_groups
  names(geneID_groups) <- geneID
  names(geneScore) <- paste(geneID, geneID_groups, sep = ".")

  if (is.null(TERM2GENE)) {
    db_list <- PrepareDB(
      species = species, db = db, db_update = db_update, db_version = db_version,
      db_IDtypes = IDtype, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror
    )
  } else {
    colnames(TERM2GENE) <- c("Term", IDtype)
    db <- "custom"
    db_list <- list()
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    if (is.null(TERM2NAME)) {
      TERM2NAME <- unique(TERM2GENE)[, c(1, 1)]
      colnames(TERM2NAME) <- c("Term", "Name")
    }
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- "custom"
  }
  if (isTRUE(db_combine)) {
    message("Create 'Combined' database ...")
    TERM2GENE <- do.call(rbind, lapply(db_list[[species]], function(x) x[["TERM2GENE"]][, c("Term", IDtype)]))
    TERM2NAME <- do.call(rbind, lapply(names(db_list[[species]]), function(x) {
      db_list[[species]][[x]][["TERM2NAME"]][["Name"]] <- paste0(db_list[[species]][[x]][["TERM2NAME"]][["Name"]], " [", x, "]")
      db_list[[species]][[x]][["TERM2NAME"]][, c("Term", "Name")]
    }))
    version <- unlist(lapply(db_list[[species]], function(x) as.character(x[["version"]])))
    version <- paste0(names(version), ":", version, collapse = ";")
    db <- "Combined"
    db_list[[species]][[db]][["TERM2GENE"]] <- unique(TERM2GENE)
    db_list[[species]][[db]][["TERM2NAME"]] <- unique(TERM2NAME)
    db_list[[species]][[db]][["version"]] <- unique(version)
  }

  if (length(unique(c(IDtype, result_IDtype))) != 1) {
    res <- GeneConvert(
      geneID = unique(geneID),
      geneID_from_IDtype = IDtype,
      geneID_to_IDtype = result_IDtype,
      species_from = species,
      species_to = species,
      Ensembl_version = Ensembl_version,
      mirror = mirror
    )
    geneMap <- res$geneID_collapse
    colnames(geneMap)[colnames(geneMap) == "from_geneID"] <- IDtype
  } else {
    geneMap <- data.frame(from_geneID = unique(geneID), row.names = unique(geneID))
    colnames(geneMap)[1] <- IDtype
  }

  input[[IDtype]] <- geneMap[as.character(input$geneID), IDtype]
  input[[result_IDtype]] <- geneMap[as.character(input$geneID), result_IDtype]
  input <- unnest(input, cols = c(IDtype, result_IDtype))
  input <- input[!is.na(input[[IDtype]]), , drop = FALSE]

  message("Permform GSEA...")
  comb <- expand.grid(group = levels(geneID_groups), term = db, stringsAsFactors = FALSE)

  res_list <- bplapply(seq_len(nrow(comb)), function(i, id) {
    group <- comb[i, "group"]
    term <- comb[i, "term"]
    geneList <- input[input$geneID_groups == group, "geneScore"]
    names(geneList) <- input[input$geneID_groups == group, IDtype]
    gene_mapid <- input[input$geneID_groups == group, result_IDtype]
    ord <- order(geneList, decreasing = TRUE)
    geneList <- geneList[ord]
    gene_mapid <- gene_mapid[ord]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", IDtype)]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[["Term"]] %in% TERM2GENE_tmp[["Term"]], , drop = FALSE]
    enrich_res <- GSEA(
      geneList = geneList,
      minGSSize = ifelse(term %in% unlimited_db, 1, minGSSize),
      maxGSSize = ifelse(term %in% unlimited_db, Inf, maxGSSize),
      nPermSimple = 1e5, # nPermSimple:fgseaMultilevel; nperm:fgseaSimple
      eps = 0,
      scoreType = scoreType,
      pAdjustMethod = "BH",
      pvalueCutoff = Inf,
      TERM2GENE = TERM2GENE_tmp,
      TERM2NAME = TERM2NAME_tmp,
      by = "fgsea",
      verbose = FALSE
    )

    if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
      result <- enrich_res@result
      result[["Groups"]] <- group
      result[["Database"]] <- term
      result[["Version"]] <- as.character(db_list[[species]][[term]][["version"]])
      IDlist <- strsplit(result$core_enrichment, "/")
      result$core_enrichment <- unlist(lapply(IDlist, function(x) {
        x_result <- NULL
        for (i in x) {
          if (i %in% input[[IDtype]]) {
            x_result <- c(x_result, unique(geneMap[geneMap[[IDtype]] == i, result_IDtype]))
          } else {
            x_result <- c(x_result, i)
          }
        }
        return(paste0(x_result, collapse = "/"))
      }))
      enrich_res@result <- result
      enrich_res@gene2Symbol <- as.character(gene_mapid)

      if (isTRUE(GO_simplify) && term %in% c("GO", "GO_BP", "GO_CC", "GO_MF")) {
        sim_res <- enrich_res
        if (term == "GO") {
          sim_res@result[["ONTOLOGY"]] <- setNames(TERM2NAME_tmp[["ONTOLOGY"]], TERM2NAME_tmp[["Term"]])[enrich_res@result[["ID"]]]
          sim_res@setType <- "GOALL"
        } else {
          sim_res@setType <- gsub(pattern = "GO_", replacement = "", x = term)
        }
        nterm_simplify <- sum(with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))))
        if (nterm_simplify <= 1) {
          warning(group, "|", term, " has no term to simplify.", immediate. = TRUE)
        } else {
          sim_res@result <- sim_res@result[with(sim_res@result, eval(rlang::parse_expr(GO_simplify_cutoff))), , drop = FALSE]
          semData <- db_list[[species]][[term]][["semData"]]
          ipclock(id)
          sim_res <- simplify(sim_res,
            measure = simplify_method,
            cutoff = simplify_similarityCutoff,
            semData = semData
          )
          ipcunlock(id)
          result_sim <- sim_res@result
          result_sim[["Groups"]] <- group
          result_sim[["Database"]] <- paste0(term, "_sim")
          result_sim[["Version"]] <- as.character(db_list[[species]][[term]][["version"]])
          result_sim[["ONTOLOGY"]] <- NULL
          sim_res@result <- result_sim
          enrich_res <- list(enrich_res, sim_res)
          names(enrich_res) <- paste(group, c(term, paste0(term, "_sim")), sep = "-")
        }
      }
    } else {
      enrich_res <- NULL
    }
    return(enrich_res)
  }, BPPARAM = BPPARAM, id = ipcid())

  nm <- paste(comb$group, comb$term, sep = "-")
  sim_index <- sapply(res_list, function(x) length(x) == 2)
  sim_list <- unlist(res_list[sim_index], recursive = FALSE)
  raw_list <- res_list[!sim_index]
  names(raw_list) <- nm[!sim_index]
  results <- c(raw_list, sim_list)
  results <- results[!sapply(results, is.null)]
  results <- results[intersect(c(nm, paste0(nm, "_sim")), names(results))]
  enrichment <- do.call(rbind, lapply(results, function(x) x@result))
  rownames(enrichment) <- NULL

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "GSEA done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  res <- list(enrichment = enrichment, results = results, geneMap = geneMap, input = input)
  if (isTRUE(use_srt)) {
    res[["DE_threshold"]] <- DE_threshold
    srt@tools[[paste("GSEA", group_by, test.use, sep = "_")]] <- res
    return(srt)
  } else {
    return(res)
  }
}

#' RunSlingshot
#'
#' Runs the Slingshot algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param group.by The variable to group the cells by.
#' @param reduction The reduction technique to use for dimensionality reduction. Default is NULL, which uses the default reduction for the Seurat object.
#' @param dims The dimensions to use for the Slingshot algorithm. Default is NULL, which uses first two dimensions.
#' @param start The starting group for the Slingshot algorithm. Default is NULL.
#' @param end The ending group for the Slingshot algorithm. Default is NULL.
#' @param prefix The prefix to add to the column names of the resulting pseudotime variable. Default is NULL.
#' @param reverse Logical value indicating whether to reverse the pseudotime variable. Default is FALSE.
#' @param align_start Logical value indicating whether to align the starting pseudotime values at the maximum pseudotime. Default is FALSE.
#' @param show_plot Logical value indicating whether to show the dimensionality plot. Default is TRUE.
#' @param lineage_palette The color palette to use for the lineages in the plot. Default is "Dark2".
#' @param seed The random seed to use for reproducibility. Default is 11.
#' @param ... Additional arguments to be passed to the \code{\link[slingshot]{slingshot}} function.
#'
#' @seealso \code{\link{CellDimPlot}} \code{\link{RunDynamicFeatures}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "PCA", dims = 1:10)
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:2), lineages_span = 0.1)
#'
#' # 3D lineage
#' pancreas_sub <- Standard_SCP(pancreas_sub)
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "StandardpcaUMAP3D")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1, lineages_trim = c(0.05, 0.95))
#' @importFrom Seurat AddMetaData as.SingleCellExperiment
#' @importFrom slingshot slingshot slingPseudotime slingBranchID
#' @export
