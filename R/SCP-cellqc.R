#' Run doublet-calling with scDblFinder
#'
#' @param srt
#' @param assay
#' @param db_rate
#' @param ...
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- db_scDblFinder(pancreas_sub)
#' CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scDblFinder_class")
#' FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scDblFinder_score")
#' @importFrom Seurat as.SingleCellExperiment
#' @export
db_scDblFinder <- function(srt, assay = "RNA", db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  check_R("scDblFinder")
  sce <- as.SingleCellExperiment(srt, assay = assay)
  sce <- scDblFinder::scDblFinder(sce, dbr = db_rate, verbose = FALSE, ...)
  srt[["db.scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["db.scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  return(srt)
}

#' Run doublet-calling with scds
#'
#' @param srt
#' @param assay
#' @param db_rate
#' @param method
#' @param ...
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#' CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scds_hybrid_class")
#' FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scds_hybrid_score")
#' @importFrom Seurat as.SingleCellExperiment
#' @export
db_scds <- function(srt, assay = "RNA", db_rate = ncol(srt) / 1000 * 0.01, method = "hybrid", ...) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  check_R("scds")
  sce <- as.SingleCellExperiment(srt, assay = assay)
  sce <- scds::cxds_bcds_hybrid(sce, ...)
  srt[["db.scds_cxds_score"]] <- sce[["cxds_score"]]
  srt[["db.scds_bcds_score"]] <- sce[["bcds_score"]]
  srt[["db.scds_hybrid_score"]] <- sce[["hybrid_score"]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_qc <- names(sort(srt[[paste0("db.scds_", method, "_score"), drop = TRUE]], decreasing = TRUE)[1:ntop])
  srt[[paste0("db.scds_", method, "_class")]] <- "singlet"
  srt[[paste0("db.scds_", method, "_class")]][db_qc, ] <- "doublet"
  return(srt)
}

#' Run doublet-calling with Scrublet
#'
#' @param srt
#' @param assay
#' @param db_rate
#' @param ...
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- db_Scrublet(pancreas_sub)
#' CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.Scrublet_class")
#' FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.Scrublet_score")
#' @importFrom reticulate import
#' @importFrom Seurat GetAssayData
#' @export
db_Scrublet <- function(srt, assay = "RNA", db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  check_Python("scrublet")
  scr <- import("scrublet")
  raw_counts <- t(as.matrix(GetAssayData(object = srt, assay = assay, slot = "counts")))
  scrub <- scr$Scrublet(raw_counts, expected_doublet_rate = db_rate, ...)
  res <- scrub$scrub_doublets()
  doublet_scores <- res[[1]]
  predicted_doublets <- res[[2]]

  srt[["db.Scrublet_score"]] <- doublet_scores
  srt[["db.Scrublet_class"]] <- sapply(predicted_doublets, function(i) {
    switch(as.character(i),
      "FALSE" = "singlet",
      "TRUE" = "doublet"
    )
  })
  return(srt)
}

#' Run doublet-calling with DoubletDetection
#'
#' @param srt
#' @param assay
#' @param db_rate
#' @param ...
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- db_DoubletDetection(pancreas_sub)
#' CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.DoubletDetection_class")
#' FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.DoubletDetection_score")
#' @importFrom reticulate import
#' @importFrom Seurat GetAssayData
#' @export
db_DoubletDetection <- function(srt, assay = "RNA", db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  check_Python("doubletdetection")
  doubletdetection <- import("doubletdetection")
  counts <- GetAssayData(object = srt, assay = assay, slot = "counts")
  clf <- doubletdetection$BoostClassifier(
    n_iters = as.integer(5),
    standard_scaling = TRUE,
    ...
  )
  labels <- clf$fit(Matrix::t(counts))$predict()
  scores <- clf$doublet_score()

  srt[["db.DoubletDetection_score"]] <- scores
  srt[["db.DoubletDetection_class"]] <- sapply(labels, function(i) {
    switch(as.character(i),
      "0" = "singlet",
      "1" = "doublet"
    )
  })
  return(srt)
}

#' Run doublet-calling for single cell RNA-seq data.
#'
#' Identification of heterotypic (or neotypic) doublets in single-cell RNAseq data.
#'
#' @param srt A \code{Seurat} object.
#' @param assay Assay to use.
#' @param db_method Doublet-calling methods used. Can be one of \code{scDblFinder}, \code{Scrublet}, \code{DoubletDetection}, \code{scds_cxds}, \code{scds_bcds}, \code{scds_hybrid}
#' @param db_rate The expected doublet rate. By default this is assumed to be 1\% per thousand cells captured (so 4\% among 4000 thousand cells), which is appropriate for 10x datasets.
#' @param ... Arguments passed to the corresponding doublet-calling method.
#'
#' @return Returns Seurat object with the doublet prediction results and prediction scores stored in the meta.data slot.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDoubletCalling(pancreas_sub, db_method = "scDblFinder")
#' CellDimPlot(pancreas_sub, reduction = "umap", group.by = "db.scDblFinder_class")
#' FeatureDimPlot(pancreas_sub, reduction = "umap", features = "db.scDblFinder_score")
#' @export
#'
RunDoubletCalling <- function(srt, assay = "RNA", db_method = "scDblFinder", db_rate = ncol(srt) / 1000 * 0.01, ...) {
  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  if (db_method %in% c("scDblFinder", "Scrublet", "DoubletDetection", "scds_cxds", "scds_bcds", "scds_hybrid")) {
    message("Run doublet-calling with ", db_method)
    methods <- unlist(strsplit(db_method, "_"))
    method1 <- methods[1]
    method2 <- methods[2]
    if (is.na(method2)) {
      args1 <- mget(names(formals()), sys.frame(sys.nframe()))
      args2 <- as.list(match.call())
    } else {
      args1 <- c(mget(names(formals()), sys.frame(sys.nframe())), method = method2)
      args2 <- c(as.list(match.call()), method = method2)
    }
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("db_method", "...")]
    tryCatch(expr = {
      srt <- do.call(
        what = paste0("db_", method1),
        args = args1
      )
    }, error = function(e) {
      message(e)
    })
    return(srt)
  } else {
    stop(paste(db_method, "is not a suppoted doublet-calling method!"))
  }
}

#' Detect outliers using MAD(Median and Median Absolute Deviation) method
#'
#' @param x
#' @param nmads
#' @param constant
#' @param type
#'
#' @importFrom stats mad
#' @export
isOutlier <- function(x, nmads = 2.5, constant = 1.4826, type = c("both", "lower", "higher")) {
  type <- match.arg(type, c("both", "lower", "higher"))
  mad <- mad(x, constant = constant, na.rm = TRUE)
  upper <- median(x, na.rm = TRUE) + nmads * mad
  lower <- median(x, na.rm = TRUE) - nmads * mad
  if (type == "both") {
    out <- which(x > upper | x < lower)
  }
  if (type == "lower") {
    out <- which(x < lower)
  }
  if (type == "higher") {
    out <- which(x > upper)
  }
  out <- c(which(is.na(x)), out)
  return(out)
}

#' Run cell-level quality control for single cell RNA-seq data.
#'
#' In CellQC, doublet-calling will be run first and then reject abnormal cell data based on median absolute deviation (MAD) outliers.
#' After doublet-calling and outlier filtering, CellQC will perform general QC (gene count, UMI count, etc.)
#' and reject cell data for non-species of interest based on the proportion and number of UMIs for the species.
#'
#' @inheritParams RunDoubletCalling
#' @param qc_metrics QC metrics applied.
#' @param return_filtered Whether to return a cell-filtered \code{Seurat} object.
#' @param outlier_cutoff MAD outlier metrics. See \link[scuttle]{isOutlier}.
#' @param outlier_n Minimum number of outlier metrics that meet the conditions for determining outlier cells. Default is 1.
#' @param UMI_threshold UMI number threshold. Cells that exceed this threshold will be considered as kept. Default is 3000.
#' @param gene_threshold Gene number threshold. Cells that exceed this threshold will be considered as kept. Default is 1000.
#' @param mito_threshold Percentage of UMI counts of mitochondrial genes. Cells that exceed this threshold will be considered as discarded. Default is 20.
#' @param mito_pattern Regex patterns to match the mitochondrial genes.
#' @param mito_gene A defined mitochondrial genes. If features provided, will ignore the \code{mito_pattern} matching.
#' @param ribo_threshold Percentage of UMI counts of ribosomal genes. Cells that exceed this threshold will be considered as discarded. Default is 50.
#' @param ribo_pattern Regex patterns to match the ribosomal genes.
#' @param ribo_gene A defined ribosomal genes. If features provided, will ignore the \code{ribo_pattern} matching.
#' @param species Species used as the suffix of the QC metrics. The first is the species of interest. Default is \code{NULL}.
#' @param species_gene_prefix Species gene prefix used to calculate QC metrics for each species. Default is \code{NULL}.
#' @param species_percent Percentage of UMI counts of the first species. Cells that exceed this threshold will be considered as kept. Default is 95.
#' @param seed Set a random seed. Default is 11.
#'
#' @note
#' General quality control usually uses metrics such as gene count, UMI count, etc.
#' In addition, ribo content and mito content can be used as QC indicators:
#' mito content can be used as an indication of apoptosis,
#' with a general threshold of 20% and less than 10% mito content in high quality cells;
#' ribo content is cell type dependent, with some cell types having even more than 50% ribo content;
#' however, ribo content in single cell data can also indicate whether the empty droplets,
#' and the ribo content in empty droplets is usually high and the mito content is low instead.
#'
#' @return Returns Seurat object with the QC results stored in the meta.data slot.
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunCellQC(pancreas_sub)
#' @importFrom Seurat Assays as.SingleCellExperiment PercentageFeatureSet WhichCells
#' @importFrom stats loess predict aggregate
#' @importFrom Matrix colSums t
#' @export
#'
RunCellQC <- function(srt, assay = "RNA",
                      qc_metrics = c("doublets", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio", "species"),
                      return_filtered = FALSE,
                      db_method = "scDblFinder", db_rate = ncol(srt) / 1000 * 0.01,
                      outlier_cutoff = c(
                        "log10_nCount:lower:2.5",
                        "log10_nCount:higher:5",
                        "log10_nFeature:lower:2.5",
                        "log10_nFeature:higher:5",
                        "featurecount_dist:lower:2.5"
                      ), outlier_n = 1,
                      UMI_threshold = 3000, gene_threshold = 1000,
                      mito_threshold = 20, mito_pattern = c("MT-", "Mt-", "mt-"), mito_gene = NULL,
                      ribo_threshold = 50, ribo_pattern = c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$", "rp[sl]\\d+\\w{0,1}\\d*$"), ribo_gene = NULL,
                      ribo_mito_ratio_range = c(1, Inf),
                      species = NULL, species_gene_prefix = NULL, species_percent = 95,
                      seed = 11) {
  set.seed(seed)

  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }
  if (!isTRUE(assay %in% Seurat::Assays(srt))) {
    stop("srt does not contain '", assay, "' assay.")
  }
  if (length(species) != length(species_gene_prefix)) {
    stop("'species_gene_prefix' must be the same length as 'species'.")
  }
  if (length(species) == 0) {
    species <- species_gene_prefix <- NULL
  }
  status <- check_DataType(srt, slot = "counts", assay = assay)
  if (status != "raw_counts") {
    stop("Data type is not raw counts!")
  }
  if (!paste0("nCount_", assay) %in% colnames(srt@meta.data)) {
    srt@meta.data[[paste0("nCount_", assay)]] <- colSums(srt[[assay]]@counts)
  }
  if (!paste0("nFeature_", assay) %in% colnames(srt@meta.data)) {
    srt@meta.data[[paste0("nFeature_", assay)]] <- colSums(srt[[assay]]@counts > 0)
  }

  ntotal <- ncol(srt)

  db_qc <- c()
  if ("doublets" %in% qc_metrics) {
    if (!is.null(db_method)) {
      for (dbm in db_method) {
        srt <- RunDoubletCalling(srt = srt, db_method = dbm, db_rate = db_rate)
        db_qc <- unique(c(db_qc, colnames(srt)[srt[[paste0("db.", dbm, "_class"), drop = TRUE]] == "doublet"]))
      }
    }
  }

  outlier_qc <- c()
  for (n in 1:length(species)) {
    if (n == 0) {
      break
    }
    sp <- species[n]
    prefix <- species_gene_prefix[n]
    sp_genes <- rownames(srt[[assay]])[grep(pattern = paste0("^", prefix), x = rownames(srt[[assay]]))]
    nCount <- srt[[paste0(c(paste0("nCount_", assay), sp), collapse = ".")]] <- colSums(srt[[assay]]@counts[sp_genes, ])
    nFeature <- srt[[paste0(c(paste0("nFeature_", assay), sp), collapse = ".")]] <- colSums(srt[[assay]]@counts[sp_genes, ] > 0)
    percent.mito <- srt[[paste0(c("percent.mito", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, assay = assay, pattern = paste0("(", paste0("^", prefix, "-*", mito_pattern), ")", collapse = "|"), features = mito_gene)[[1]]
    percent.ribo <- srt[[paste0(c("percent.ribo", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, assay = assay, pattern = paste0("(", paste0("^", prefix, "-*", ribo_pattern), ")", collapse = "|"), features = ribo_gene)[[1]]
    percent.ribo <- srt[[paste0(c("ribo.mito.ratio", sp), collapse = ".")]] <- srt[[paste0(c("percent.ribo", sp), collapse = "."), drop = TRUE]] / srt[[paste0(c("percent.mito", sp), collapse = "."), drop = TRUE]]
    percent.genome <- srt[[paste0(c("percent.genome", sp), collapse = ".")]] <- PercentageFeatureSet(object = srt, assay = assay, pattern = paste0("^", prefix))[[1]]

    if (n == 1) {
      if ("outlier" %in% qc_metrics) {
        message("Calculate outlier cells")

        # "percent.top_20:higher:5"
        # countx <- as(srt[[assay]]@counts[sp_genes, ], "sparseMatrix")
        # agg <- aggregate(x = countx@x, by = list(rep(colnames(countx), diff(countx@p))), FUN = function(x) {
        #   sum(head(x, 20))
        # })
        # rownames(agg) <- agg[[1]]
        # percent.top_20 <- srt[[paste0(c("percent.top_20", sp), collapse = ".")]] <- agg[colnames(srt), "x"]

        log10_nFeature <- srt[[paste0(c(paste0("log10_nFeature_", assay), sp), collapse = ".")]] <- log10(nFeature)
        log10_nCount <- srt[[paste0(c(paste0("log10_nCount_", assay), sp), collapse = ".")]] <- log10(nCount)
        log10_nCount[is.infinite(log10_nCount)] <- NA
        log10_nFeature[is.infinite(log10_nFeature)] <- NA
        mod <- loess(log10_nFeature ~ log10_nCount, span = 1)
        pred <- predict(mod, newdata = data.frame(log10_nCount = log10_nCount))
        featurecount_dist <- srt[[paste0(c("featurecount_dist", sp), collapse = ".")]] <- log10_nFeature - pred

        # df <- data.frame(cell = colnames(srt), ribo = srt$percent.ribo.Homo_sapiens, y = log10_nFeature, x = log10_nCount, pred = pred, featurecount_dist = featurecount_dist)
        # lower_df <- subset(df, featurecount_dist < median(df$featurecount_dist) - 2.5 * mad(df$featurecount_dist))
        # higher_df <- subset(df, featurecount_dist > median(df$featurecount_dist) + 2.5 * mad(df$featurecount_dist))
        # ggplot(df) +
        #   geom_point(aes(x = x, y = y, color = featurecount_dist)) +
        #   scale_color_gradientn(colors = c("green", "white", "orange"), values = scales::rescale(c(min(df$featurecount_dist), 0, max(df$featurecount_dist)))) +
        #   geom_point(data = lower_df, aes(x = x, y = y), shape = 21, fill = "transparent", color = "blue") +
        #   geom_point(data = higher_df, aes(x = x, y = y), shape = 21, fill = "transparent", color = "red") +
        #   geom_line(aes(x = x, y = pred), color = "black")+
        #   theme(panel.background = element_rect(fill = "grey"))
        # nrow(lower_df)

        var <- sapply(strsplit(outlier_cutoff, ":"), function(x) x[[1]])
        var_valid <- var %in% colnames(srt@meta.data) | sapply(var, FUN = function(x) exists(x, where = environment()))
        if (any(!var_valid)) {
          stop("Variable ", paste0(names(var_valid)[!var_valid], collapse = ","), " is not found in the srt object.")
        }
        outlier <- lapply(strsplit(outlier_cutoff, ":"), function(m) {
          colnames(srt)[isOutlier(get(m[1]), nmads = as.numeric(m[3]), type = m[2])]
        })
        names(outlier) <- outlier_cutoff
        print(unlist(lapply(outlier, length)))
        outlier_tb <- table(unlist(outlier))
        outlier_qc <- c(outlier_qc, names(outlier_tb)[outlier_tb >= outlier_n])
        for (nm in names(outlier)) {
          srt[[make.names(nm)]] <- colnames(srt) %in% outlier[[nm]]
        }
      }
    }
  }

  umi_qc <- gene_qc <- mito_qc <- ribo_qc <- ribo_mito_ratio_qc <- species_qc <- c()
  if ("umi" %in% qc_metrics) {
    umi_qc <- colnames(srt)[srt[[paste0(c(paste0("nCount_", assay), species[1]), collapse = "."), drop = TRUE]] < UMI_threshold]
  }
  if ("gene" %in% qc_metrics) {
    gene_qc <- colnames(srt)[srt[[paste0(c(paste0("nFeature_", assay), species[1]), collapse = "."), drop = TRUE]] < gene_threshold]
  }
  if ("mito" %in% qc_metrics) {
    mito_qc <- colnames(srt)[srt[[paste0(c("percent.mito", species[1]), collapse = "."), drop = TRUE]] > mito_threshold]
  }
  if ("ribo" %in% qc_metrics) {
    ribo_qc <- colnames(srt)[srt[[paste0(c("percent.ribo", species[1]), collapse = "."), drop = TRUE]] > ribo_threshold]
  }
  if ("ribo_mito_ratio" %in% qc_metrics) {
    ribo_mito_ratio_qc <- colnames(srt)[srt[[paste0(c("ribo.mito.ratio", species[1]), collapse = "."), drop = TRUE]] < ribo_mito_ratio_range[1] | srt[[paste0(c("ribo.mito.ratio", species[1]), collapse = "."), drop = TRUE]] > ribo_mito_ratio_range[2]]
  }
  if ("species" %in% qc_metrics) {
    species_qc <- colnames(srt)[srt[[paste0(c("percent.genome", species[1]), collapse = "."), drop = TRUE]] < species_percent]
  }

  CellQC <- unique(c(db_qc, outlier_qc, umi_qc, gene_qc, mito_qc, ribo_qc, ribo_mito_ratio_qc, species_qc))
  cat(">>>", "Total cells:", ntotal, "\n")
  cat(">>>", "Cells which are filtered out:", length(CellQC), "\n")
  cat("...", length(db_qc), "potential doublets", "\n")
  cat("...", length(outlier_qc), "outlier cells", "\n")
  cat("...", length(umi_qc), "low-UMI cells", "\n")
  cat("...", length(gene_qc), "low-gene cells", "\n")
  cat("...", length(mito_qc), "high-mito cells", "\n")
  cat("...", length(ribo_qc), "high-ribo cells", "\n")
  cat("...", length(ribo_mito_ratio_qc), " ribo_mito_ratio outlier cells", "\n")
  cat("...", length(species_qc), "species-contaminated cells", "\n")
  cat(">>>", "Remained cells after filtering:", ntotal - length(CellQC), "\n")

  qc_nm <- c("db_qc", "outlier_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc", "CellQC")
  for (qc in qc_nm) {
    srt[[qc]] <- ifelse(colnames(srt) %in% get(qc), "Fail", "Pass")
    srt[[qc]] <- factor(srt[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
  }

  if (return_filtered) {
    srt <- srt[, srt$CellQC == "Pass"]
    srt@meta.data[, intersect(qc_nm, colnames(srt@meta.data))] <- NULL
  }

  return(srt)
}
