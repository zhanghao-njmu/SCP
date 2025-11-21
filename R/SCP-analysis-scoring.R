CellScoring <- function(srt, features = NULL, slot = "data", assay = NULL, split.by = NULL,
                        IDtype = "symbol", species = "Homo_sapiens",
                        db = "GO_BP", termnames = NULL, db_update = FALSE, db_version = "latest", convert_species = TRUE,
                        Ensembl_version = 103, mirror = NULL, minGSSize = 10, maxGSSize = 500,
                        method = "Seurat", classification = TRUE, name = "", new_assay = FALSE,
                        BPPARAM = BiocParallel::bpparam(), seed = 11, ...) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed

  if (!method %in% c("Seurat", "AUCell", "UCell")) {
    stop("method must be 'Seurat', 'AUCell'or 'UCell'.")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (slot == "counts") {
    status <- check_DataType(srt, slot = "counts", assay = assay)
    if (status != "raw_counts") {
      warning("Data is not raw counts", immediate. = TRUE)
    }
  }
  if (slot == "data") {
    status <- check_DataType(srt, slot = "data", assay = assay)
    if (status == "raw_counts") {
      message("Data is raw counts. Perform NormalizeData(LogNormalize) on the data ...")
      srt <- NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (status == "raw_normalized_counts") {
      message("Data is normalized without log transformation. Perform NormalizeData(LogNormalize) on the data...")
      srt <- NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (status == "unknown") {
      warning("Can not determine whether data ", i, " is log-normalized...\n", immediate. = TRUE)
    }
  }
  if (name == "" && isTRUE(new_assay)) {
    stop("name must be specified when new_assay=TRUE")
  }
  if (is.null(features)) {
    for (single_db in db) {
      db_list <- PrepareDB(
        species = species, db = single_db, db_update = db_update, db_version = db_version,
        db_IDtypes = IDtype, convert_species = convert_species, Ensembl_version = Ensembl_version, mirror = mirror
      )
      TERM2GENE_tmp <- db_list[[species]][[single_db]][["TERM2GENE"]][, c("Term", IDtype)]
      TERM2NAME_tmp <- db_list[[species]][[single_db]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], , drop = FALSE]

      TERM2GENE_tmp <- unique(TERM2GENE_tmp)
      TERM2NAME_tmp <- unique(TERM2NAME_tmp)
      rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
      features_tmp <- split(TERM2GENE_tmp[, IDtype], TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"])
      if (is.null(termnames)) {
        GSSize <- sapply(features_tmp, length)
        features_tmp <- features_tmp[GSSize >= minGSSize & GSSize <= maxGSSize]
      } else {
        if (length(intersect(termnames, names(features_tmp)) > 0)) {
          features_tmp <- features_tmp[intersect(termnames, names(features_tmp))]
        } else {
          stop("None of termnames found in the db: ", single_db)
        }
      }
      features <- c(features, features_tmp)
    }
  }

  if (!is.list(features) || length(names(features)) == 0) {
    stop("'features' must be a named list")
  }
  expressed <- names(which(rowSums(GetAssayData(srt, slot = slot, assay = assay) > 0) > 0))
  features <- lapply(setNames(names(features), names(features)), function(x) features[[x]][features[[x]] %in% expressed])
  filtered_none <- names(which(sapply(features, length) == 0))
  if (length(filtered_none) > 0) {
    warning("The following list of features were filtered because none of features were found in the srt assay:\n", paste0(filtered_none, collapse = ", "), immediate. = TRUE)
  }
  features <- features[!names(features) %in% filtered_none]
  features_raw <- features
  names(features) <- make.names(names(features))
  message("Number of feature lists to be scored: ", length(features))

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start CellScoring"))
  message("Workers: ", bpworkers(BPPARAM))

  if (!is.null(split.by)) {
    split_list <- SplitObject(srt, split.by = split.by)
  } else {
    split_list <- list(srt)
  }
  scores_list <- list()
  features_nm_list <- list()
  for (i in seq_along(split_list)) {
    srt_sp <- split_list[[i]]
    if (method == "Seurat") {
      ## need to add a 'slot' parameter
      srt_tmp <- AddModuleScore2(
        srt_sp,
        features = features,
        name = name,
        slot = slot,
        assay = assay,
        BPPARAM = BPPARAM,
        ...
      )
      if (name != "") {
        scores <- srt_tmp[[paste0(name, seq_along(features))]]
      } else {
        scores <- srt_tmp[[paste0("X", seq_along(features))]]
      }
      features_nm <- features_raw
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    } else if (method == "UCell") {
      check_R("UCell")
      srt_tmp <- UCell::AddModuleScore_UCell(
        srt_sp,
        features = features,
        name = name,
        slot = slot,
        assay = assay,
        BPPARAM = BPPARAM,
        ...
      )
      filtered <- names(features)[!paste0(names(features), name) %in% colnames(srt_tmp@meta.data)]
      if (length(filtered) > 0) {
        warning("The following list of features were filtered when scoring:\n", paste0(filtered, collapse = ", "), immediate. = TRUE)
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- srt_tmp[[paste0(names(features_keep), name)]]
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    } else if (method == "AUCell") {
      check_R("AUCell")
      CellRank <- AUCell::AUCell_buildRankings(as_matrix(GetAssayData(srt_sp, slot = slot, assay = assay)), BPPARAM = BPPARAM, plotStats = FALSE)
      cells_AUC <- AUCell::AUCell_calcAUC(
        geneSets = features,
        rankings = CellRank,
        ...
      )
      filtered <- names(features)[!names(features) %in% rownames(AUCell::getAUC(cells_AUC))]
      if (length(filtered) > 0) {
        warning("The following list of features were filtered when scoring:", paste0(filtered, collapse = ", "), immediate. = TRUE)
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- as.data.frame(t(AUCell::getAUC(cells_AUC)))[, names(features_keep)]
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    }
    features_nm_list[[i]] <- setNames(object = names(features_nm), nm = colnames(scores))
    scores_list[[i]] <- scores
  }
  features_used <- Reduce(intersect, lapply(scores_list, colnames))
  features_nm_used <- Reduce(intersect, features_nm_list)
  scores_mat <- do.call(rbind, lapply(scores_list, function(x) x[intersect(rownames(x), colnames(srt)), features_used, drop = FALSE]))

  if (isTRUE(new_assay)) {
    srt[[name]] <- CreateAssayObject(counts = t(as_matrix(scores_mat[colnames(srt), , drop = FALSE])))
    srt[[name]] <- AddMetaData(object = srt[[name]], metadata = data.frame(termnames = features_nm_used[colnames(scores_mat)]))
  } else {
    srt <- AddMetaData(object = srt, metadata = scores_mat)
  }

  if (isTRUE(classification)) {
    assignments <- apply(scores_mat, MARGIN = 1, FUN = function(x) {
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
    srt[[paste0(name, "_classification")]] <- assignments[rownames(scores_mat)]
  }

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "CellScoring done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  return(srt)
}

metap <- function(p, method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"), ...) {
  method <- match.arg(method)
  res <- do.call(method, args = list(p = p, ...))
  return(res)
}

#' @importFrom stats pbeta
wilkinsonp <- function(p, r = 1, alpha = 0.05, log.p = FALSE) {
  alpha <- ifelse(alpha > 1, alpha / 100, alpha)
  stopifnot(alpha > 0, alpha < 1)
  alpha <- ifelse(alpha > 0.5, 1 - alpha, alpha)
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(
      p = NA_real_, pr = NA_real_, r = r, critp = NA_real_,
      alpha = alpha, validp = p[keep]
    )
  } else {
    pi <- p[keep]
    k <- length(pi)
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if ((r < 1) | (r > k)) {
      r <- 1
      warning("Illegal r set to 1")
    }
    pi <- sort(pi)
    pr <- pi[r]
    res <- list(
      p = pbeta(pr, r, k + 1 - r, log.p = log.p),
      pr = pr, r = r, critp = qbeta(alpha, r, k + 1 - r),
      alpha = alpha, validp = pi
    )
  }
  res
}

maximump <- function(p, alpha = 0.05, log.p = FALSE) {
  keep <- (p >= 0) & (p <= 1)
  validp <- p[keep]
  k <- length(validp)
  res <- wilkinsonp(p, r = k, alpha, log.p)
  res
}

minimump <- function(p, alpha = 0.05, log.p = FALSE) {
  res <- wilkinsonp(p, r = 1, alpha, log.p)
  res
}

#' @importFrom stats pnorm
meanp <- function(p) {
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 4
  if (invalid) {
    warning("Must have at least four valid p values")
    res <- list(z = NA_real_, p = NA_real_, validp = p[keep])
  } else {
    pi <- mean(p[keep])
    k <- length(p[keep])
    z <- (0.5 - pi) * sqrt(12 * k)
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    res <- list(
      z = z, p = pnorm(z, lower.tail = FALSE),
      validp = p[keep]
    )
  }
  res
}

#' @importFrom stats pnorm
sump <- function(p) {
  keep <- (p >= 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(p = NA_real_, conservativep = NA_real_, validp = p[keep])
  } else {
    sigmap <- sum(p[keep])
    k <- length(p[keep])
    conservativep <- exp(k * log(sigmap) - lgamma(k + 1))
    nterm <- floor(sigmap) + 1
    denom <- lfactorial(k)
    psum <- 0
    terms <- vector("numeric", nterm)
    for (i in 1:nterm) {
      terms[i] <- lchoose(k, i - 1) + k * log(sigmap -
        i + 1) - denom
      pm <- 2 * (i %% 2) - 1
      psum <- psum + pm * exp(terms[i])
    }
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if (sigmap > 20) {
      warning("Likely to be unreliable, check with another method")
    }
    res <- list(
      p = psum, conservativep = conservativep,
      validp = p[keep]
    )
  }
  res
}

#' @importFrom stats binom.test
votep <- function(p, alpha = 0.5) {
  alpha <- ifelse(alpha > 1, alpha / 100, alpha)
  stopifnot(alpha > 0, alpha < 1)
  keep <- (p >= 0) & (p <= 1)
  alp <- vector("numeric", 2)
  if (alpha <= 0.5) {
    alp[1] <- alpha
    alp[2] <- 1 - alpha
  } else {
    alp[2] <- alpha
    alp[1] <- 1 - alpha
  }
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(
      p = NA_real_, pos = NA_integer_, neg = NA_integer_,
      alpha = alpha, validp = p[keep]
    )
  } else {
    pi <- p[keep]
    k <- length(pi)
    pos <- sum(1L * (pi < alp[1]))
    neg <- sum(1L * (pi > alp[2]))
    if (k != length(p)) {
      warning("Some studies omitted")
    }
    if ((pos + neg) <= 0) {
      warning("All p values are within specified limits of alpha")
      p <- 1
    } else {
      p <- binom.test(pos, pos + neg, 0.5, alternative = "greater")$p.value
    }
    res <- list(
      p = p, pos = pos, neg = neg, alpha = alpha,
      validp = pi
    )
  }
  res
}

#' @importFrom SeuratObject PackageCheck FetchData WhichCells SetIdent Idents
#' @importFrom Seurat FindMarkers FoldChange
FindConservedMarkers2 <- function(object, grouping.var, ident.1, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL, features = NULL,
                                  test.use = "wilcox", logfc.threshold = 0.25, base = 2, pseudocount.use = 1, mean.fxn = NULL,
                                  min.pct = 0.1, min.diff.pct = -Inf, max.cells.per.ident = Inf, latent.vars = NULL, only.pos = FALSE,
                                  assay = NULL, slot = "data", min.cells.group = 3, min.cells.feature = 3,
                                  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
                                  norm.method = "LogNormalize", verbose = TRUE, ...) {
  meta.method <- match.arg(meta.method)
  object.var <- FetchData(object = object, vars = grouping.var)
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)
  assay <- assay %||% DefaultAssay(object)

  cells <- list()
  for (i in 1:num.groups) {
    cells[[i]] <- rownames(x = object.var[object.var[, 1] == levels.split[i], , drop = FALSE])
  }
  marker.test <- list()
  if (is.null(cells.1)) {
    cells.1 <- WhichCells(object = object, idents = ident.1)
    if (!is.null(ident.2)) {
      cells.2 <- cells.2 %||% WhichCells(object = object, idents = ident.2)
    } else {
      cells.2 <- setdiff(colnames(object), cells.1)
    }
    object <- SetIdent(object = object, cells = colnames(x = object), value = paste(Idents(object = object), object.var[, 1], sep = "_"))
    ident.2.save <- ident.2
    for (i in 1:num.groups) {
      level.use <- levels.split[i]
      ident.use.1 <- paste(ident.1, level.use, sep = "_")
      ident.use.1.exists <- ident.use.1 %in% Idents(object = object)
      if (!all(ident.use.1.exists)) {
        bad.ids <- ident.1[!ident.use.1.exists]
        warning("Identity: ", paste(bad.ids, collapse = ", "), " not present in group ", level.use, ". Skipping ", level.use, call. = FALSE, immediate. = TRUE)
        next
      }
      ident.2 <- ident.2.save
      cells.1.use <- WhichCells(object = object, idents = ident.use.1)
      if (length(cells.1.use) < min.cells.group) {
        warning(level.use, " has fewer than ", min.cells.group, " cells in Identity: ", paste(ident.1, collapse = ", "), ". Skipping ", level.use, call. = FALSE, immediate. = TRUE)
        next
      }
      if (is.null(x = ident.2)) {
        cells.2.use <- setdiff(x = cells[[i]], y = cells.1.use)
        ident.use.2 <- names(x = which(x = table(Idents(object = object)[cells.2.use]) > 0))
        ident.2 <- gsub(pattern = paste0("_", level.use), replacement = "", x = ident.use.2)
        if (length(x = ident.use.2) == 0) {
          stop(paste("Only one identity class present:", ident.1))
        }
      } else {
        ident.use.2 <- paste(ident.2, level.use, sep = "_")
        cells.2.use <- WhichCells(object = object, idents = ident.use.2)
      }
      if (length(cells.2.use) < min.cells.group) {
        warning(level.use, " has fewer than ", min.cells.group, " cells. Skipping ", level.use, call. = FALSE, immediate. = TRUE)
        next
      }
      if (verbose) {
        message("Testing group ", level.use, ": (", paste(ident.1, collapse = ", "), ") vs (", paste(ident.2, collapse = ", "), ")")
      }
      ident.use.2.exists <- ident.use.2 %in% Idents(object = object)
      if (!all(ident.use.2.exists)) {
        bad.ids <- ident.2[!ident.use.2.exists]
        warning("Identity: ", paste(bad.ids, collapse = ", "),
          " not present in group ", level.use, ". Skipping ",
          level.use,
          call. = FALSE, immediate. = TRUE
        )
        next
      }
      marker.test[[level.use]] <- FindMarkers(
        object = Assays(object, assay), slot = slot, cells.1 = cells.1.use, cells.2 = cells.2.use, features = features,
        test.use = test.use, logfc.threshold = logfc.threshold,
        min.pct = min.pct, min.diff.pct = min.diff.pct, max.cells.per.ident = max.cells.per.ident,
        min.cells.group = min.cells.group, min.cells.feature = min.cells.feature,
        norm.method = norm.method, base = base, pseudocount.use = pseudocount.use, mean.fxn = mean.fxn,
        latent.vars = latent.vars, only.pos = only.pos,
        verbose = verbose, ...
      )
    }
  } else {
    for (i in 1:num.groups) {
      level.use <- levels.split[i]
      cells.1.use <- intersect(cells[[i]], cells.1)
      if (length(cells.1.use) < min.cells.group) {
        warning(level.use, " has fewer than ", min.cells.group, " cells. Skipping ", level.use, call. = FALSE, immediate. = TRUE)
        next
      }
      if (is.null(x = cells.2)) {
        cells.2.use <- setdiff(x = cells[[i]], y = cells.1.use)
      } else {
        cells.2.use <- intersect(cells[[i]], cells.2)
      }
      if (length(cells.2.use) < min.cells.group) {
        warning(level.use, " has fewer than ", min.cells.group, " cells. Skipping ", level.use, call. = FALSE, immediate. = TRUE)
        next
      }
      if (verbose) {
        message("Testing group ", level.use, ": (", paste("cells.1", collapse = ", "), ") vs (", paste("cells.2", collapse = ", "), ")")
      }
      marker.test[[level.use]] <- FindMarkers(
        object = Assays(object, assay), slot = slot, cells.1 = cells.1.use, cells.2 = cells.2.use, features = features,
        test.use = test.use, logfc.threshold = logfc.threshold,
        min.pct = min.pct, min.diff.pct = min.diff.pct, max.cells.per.ident = max.cells.per.ident,
        min.cells.group = min.cells.group, min.cells.feature = min.cells.feature,
        norm.method = norm.method, base = base, pseudocount.use = pseudocount.use, mean.fxn = mean.fxn,
        latent.vars = latent.vars, only.pos = only.pos,
        verbose = verbose, ...
      )
    }
  }
  marker.test <- marker.test[!sapply(marker.test, is.null)]
  if (length(marker.test) == 0) {
    warning("No group was tested", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  genes.conserved <- Reduce(f = intersect, x = lapply(X = marker.test, FUN = function(x) {
    return(rownames(x = x))
  }))
  markers.conserved <- list()
  for (i in 1:length(x = marker.test)) {
    markers.conserved[[i]] <- marker.test[[i]][genes.conserved, , drop = FALSE]
    colnames(x = markers.conserved[[i]]) <- paste(names(x = marker.test)[i], colnames(x = markers.conserved[[i]]), sep = "_")
  }
  markers.combined <- Reduce(cbind, markers.conserved)
  fc <- FoldChange(Assays(object, assay), slot = slot, cells.1 = cells.1, cells.2 = cells.2, features = genes.conserved, norm.method = norm.method, base = base, pseudocount.use = pseudocount.use, mean.fxn = mean.fxn)
  markers.combined <- cbind(markers.combined, fc[genes.conserved, , drop = FALSE])
  logFC.codes <- colnames(x = markers.combined)[grepl(pattern = "*avg_log.*FC$", x = colnames(x = markers.combined))]
  if (isTRUE(only.pos)) {
    markers.combined <- markers.combined[apply(markers.combined[, logFC.codes] > 0, 1, all), , drop = FALSE]
  } else {
    markers.combined <- markers.combined[apply(markers.combined[, logFC.codes] < 0, 1, all) | apply(markers.combined[, logFC.codes] > 0, 1, all), , drop = FALSE]
  }
  pval.codes <- colnames(x = markers.combined)[grepl(pattern = "*_p_val$", x = colnames(x = markers.combined))]
  if (length(x = pval.codes) > 1) {
    markers.combined[["max_pval"]] <- apply(X = markers.combined[, pval.codes, drop = FALSE], MARGIN = 1, FUN = max)
    combined.pval <- data.frame(cp = apply(X = markers.combined[, pval.codes, drop = FALSE], MARGIN = 1, FUN = function(x) {
      return(metap(x, method = meta.method)$p)
    }))
    meta.method.name <- meta.method
    colnames(x = combined.pval) <- paste0(meta.method.name, "_p_val")
    markers.combined <- cbind(markers.combined, combined.pval)
    markers.combined[, "p_val"] <- markers.combined[, paste0(meta.method.name, "_p_val")]
    markers.combined <- markers.combined[order(markers.combined[, paste0(meta.method.name, "_p_val")]), , drop = FALSE]
  } else {
    markers.combined[, "max_pval"] <- markers.combined[, "p_val"] <- markers.combined[, pval.codes]
    warning("Only a single group was tested", call. = FALSE, immediate. = TRUE)
  }
  return(markers.combined)
}

#' @examples
#' markers <- FindExpressedMarkers(pancreas_sub, cells.1 = WhichCells(pancreas_sub, expression = Phase == "G2M"))
#' head(markers)
#' FeatureStatPlot(pancreas_sub, rownames(markers)[1], group.by = "Phase", add_point = TRUE)
#' @importFrom Matrix rowSums
#' @importFrom SeuratObject PackageCheck FetchData WhichCells SetIdent Idents
#' @importFrom Seurat FindMarkers FoldChange Command
#' @importFrom pbapply pbsapply
#' @export
