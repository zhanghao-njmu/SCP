#' AnnotateFeatures
#'
#' @param srt
#' @param gtf
#' @param attribute
#' @param assays
#' @param overwrite
#' @param species
#' @param anno_TF Whether to add the transcription factor/cofactor annotation.
#' @param anno_LR Whether to add the ligand/receptor annotation.
#' @param merge_gtf_by
#' @param IDtype
#' @param anno_db
#' @param db
#' @param db_update
#' @param db_version
#' @param Ensembl_version
#' @param mirror
#'
#' @return
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- AnnotateFeatures(pancreas_sub, species = "Mus_musculus", anno_TF = TRUE, anno_LR = TRUE, anno_db = TRUE, db = "PFAM")
#' head(pancreas_sub[["RNA"]]@meta.features)
#'
#' ## Annotate features using a GTF file
#' # pancreas_sub <- AnnotateFeatures(pancreas_sub, gtf = "/data/reference/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf")
#'
AnnotateFeatures <- function(srt, species = "Homo_sapiens", IDtype = c("symbol", "ensembl_id", "entrez_id"),
                             anno_TF = FALSE, anno_LR = FALSE, anno_db = FALSE,
                             db = "GO_BP", db_update = FALSE, db_version = "latest", Ensembl_version = 103, mirror = NULL,
                             gtf = NULL, merge_gtf_by = "gene_name", attribute = c("gene_id", "gene_name", "gene_type"),
                             assays = "RNA", overwrite = FALSE) {
  IDtype <- match.arg(IDtype)
  merge_tf_by <- switch(IDtype,
    "symbol" = "Symbol",
    "ensembl_id" = "Ensembl",
    "entrez_id" = "Entrez.ID"
  )
  merge_lr_by <- switch(IDtype,
    "symbol" = "gene_symbol",
    "ensembl_id" = "ensembl_gene_id",
    "entrez_id" = "gene_id"
  )

  if (!is.null(gtf)) {
    gtf_all <- suppressWarnings(data.table::fread(gtf, sep = "\t"))
    colnames(gtf_all) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    for (type in c("gene", "transcript", "exon", "CDS")) {
      if (type %in% gtf_all[["feature"]]) {
        gtf_all <- gtf_all[gtf_all[["feature"]] == type, ]
        break
      }
    }
    gtf_attribute <- unique(gtf_all[["attribute"]])
    gtf_attribute <- gsub(pattern = "\"", replacement = "", x = gtf_attribute)
    gtf_attribute <- strsplit(gtf_attribute, split = "; *")
    gene_attr <- lapply(gtf_attribute, function(x) {
      detail <- strsplit(x, " ")
      out <- lapply(detail, function(x) x[2:length(x)])
      names(out) <- sapply(detail, function(x) x[1])
      out <- out[intersect(attribute, names(out))]
      return(out)
    })
    gene_attr_df <- data.table::rbindlist(unique(gene_attr), fill = TRUE)
    colnames(gene_attr_df) <- make.unique(colnames(gene_attr_df))
    gene_attr_df_collapse <- aggregate(gene_attr_df, by = list(rowid = gene_attr_df[[merge_gtf_by]]), FUN = function(x) {
      paste0(unique(x), collapse = ";")
    })
    rownames(gene_attr_df_collapse) <- gene_attr_df_collapse[["rowid"]]
    gene_attr_df_collapse[["rowid"]] <- NULL
    for (assay in assays) {
      meta.features <- srt[[assay]]@meta.features
      if (length(intersect(colnames(meta.features), colnames(gene_attr_df_collapse))) > 0 && isTRUE(overwrite)) {
        meta.features <- meta.features[, setdiff(colnames(meta.features), colnames(gene_attr_df_collapse))]
      }
      meta.features <- cbind(meta.features, gene_attr_df_collapse[rownames(meta.features), setdiff(colnames(gene_attr_df_collapse), colnames(meta.features)), drop = FALSE])
      srt[[assay]]@meta.features <- meta.features
    }
  }
  if (isTRUE(anno_TF)) {
    for (tf_type in c("TF", "TF_cofactors")) {
      url <- paste0("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/", species, "_", tf_type)
      temp <- tempfile()
      methods <- c("auto", "wget", "libcurl", "curl", "internal", "wininet")
      for (method in methods) {
        status <- tryCatch(expr = {
          download.file(url, destfile = temp, method = method, quiet = TRUE)
        }, error = function(e) {
          message("Get errors when connecting with AnimalTFDB...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        })
        if (!is.null(status)) {
          break
        }
      }
      if (is.null(status)) {
        stop("Error when try to download TF data.")
      }
      tf <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
      unlink(temp)
      for (assay in assays) {
        meta.features <- srt[[assay]]@meta.features
        if (any(c(tf_type, paste0(tf_type, "_Family")) %in% colnames(meta.features)) && isTRUE(overwrite)) {
          meta.features <- meta.features[, setdiff(colnames(meta.features), c(tf_type, paste0(tf_type, "_Family")))]
        }
        tf_sub <- unique(tf[tf[[merge_tf_by]] %in% rownames(meta.features), c(merge_tf_by, "Family")])
        if (nrow(tf_sub) == 0) {
          stop(paste0("No TF found in the seurat object. Please check if the species name is correct. The expected feature names are ", paste(head(tf_sub[[merge_tf_by]], 10), collapse = ","), "."))
        }
        rownames(tf_sub) <- tf_sub[[merge_tf_by]]
        meta.features[, tf_type] <- rownames(meta.features) %in% tf_sub[[merge_tf_by]]
        meta.features[rownames(tf_sub), paste0(tf_type, "_Family")] <- tf_sub[, "Family"]
        srt[[assay]]@meta.features <- meta.features
      }
    }
  }
  if (isTRUE(anno_LR)) {
    if (!species %in% c("Homo_sapiens", "Mus_musculus")) {
      stop("Ligand-receptor annotation is currently supported only for Homo_sapiens and Mus_musculus.")
    }
    url <- paste0("http://tcm.zju.edu.cn/celltalkdb/download/processed_data/", switch(species,
      "Homo_sapiens" = "human_lr_pair.txt",
      "Mus_musculus" = "mouse_lr_pair.txt"
    ))
    temp <- tempfile()
    methods <- c("auto", "wget", "libcurl", "curl", "internal", "wininet")
    for (method in methods) {
      status <- tryCatch(expr = {
        download.file(url, destfile = temp, method = method, quiet = TRUE)
      }, error = function(e) {
        message("Get errors when connecting with AnimalTFDB...\nRetrying...")
        Sys.sleep(1)
        return(NULL)
      })
      if (!is.null(status)) {
        break
      }
    }
    if (is.null(status)) {
      stop("Error when try to download TF data.")
    }
    lr <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
    unlink(temp)
    ligand <- aggregate(lr, by = list(rowid = lr[[paste0("ligand_", merge_lr_by)]]), FUN = function(x) {
      paste0(unique(x), collapse = ";")
    })
    rownames(ligand) <- ligand[["rowid"]]
    ligand[["LR_type"]] <- "ligand"
    ligand <- ligand[, c("LR_type", "receptor_gene_symbol", "receptor_gene_id", "receptor_ensembl_gene_id", "receptor_ensembl_protein_id", "evidence")]
    colnames(ligand) <- c("LR_type", "paired_gene_symbol", "paired_gene_id", "paired_ensembl_gene_id", "paired_ensembl_protein_id", "evidence")
    receptor <- aggregate(lr, by = list(rowid = lr[[paste0("receptor_", merge_lr_by)]]), FUN = function(x) {
      paste0(unique(x), collapse = ";")
    })
    rownames(receptor) <- receptor[["rowid"]]
    receptor[["LR_type"]] <- "receptor"
    receptor <- receptor[, c("LR_type", "ligand_gene_symbol", "ligand_gene_id", "ligand_ensembl_gene_id", "ligand_ensembl_protein_id", "evidence")]
    colnames(receptor) <- c("LR_type", "paired_gene_symbol", "paired_gene_id", "paired_ensembl_gene_id", "paired_ensembl_protein_id", "evidence")
    lr_df <- rbind(ligand, receptor)
    for (assay in assays) {
      meta.features <- srt[[assay]]@meta.features
      if (any(colnames(lr_df) %in% colnames(meta.features)) && isTRUE(overwrite)) {
        meta.features <- meta.features[, setdiff(colnames(meta.features), colnames(lr_df))]
      }
      lr_sub <- unique(lr_df[rownames(lr_df) %in% rownames(meta.features), ])
      if (nrow(lr_sub) == 0) {
        stop(paste0("No LR found in the seurat object. Please check if the species name is correct. The expected feature names are ", paste(head(lr_sub[[merge_lr_by]], 10), collapse = ","), "."))
      }
      meta.features <- cbind(meta.features, lr_sub[rownames(meta.features), setdiff(colnames(lr_sub), colnames(meta.features)), drop = FALSE])
      srt[[assay]]@meta.features <- meta.features
    }
  }
  if (isTRUE(anno_db)) {
    db_list <- PrepareEnrichmentDB(
      species = species, db = db, db_update = db_update, db_version = db_version,
      db_IDtypes = IDtype, Ensembl_version = Ensembl_version, mirror = mirror
    )
    TERM2GENE <- unique(db_list[[species]][[db]][["TERM2GENE"]])
    TERM2NAME <- unique(db_list[[species]][[db]][["TERM2NAME"]])
    rownames(TERM2NAME) <- TERM2NAME[, 1]
    TERM2GENE[, db] <- TERM2NAME[TERM2GENE[, 1], 2]
    db_df <- aggregate(
      x = TERM2GENE[, !colnames(TERM2GENE) %in% c("Term", "entrez_id", "symbol", "ensembl_id"), drop = FALSE],
      by = list(rowid = TERM2GENE[[IDtype]]), FUN = function(x) {
        paste0(unique(x), collapse = ";")
      }
    )
    rownames(db_df) <- db_df[["rowid"]]
    db_df[["rowid"]] <- NULL
    for (assay in assays) {
      meta.features <- srt[[assay]]@meta.features
      if (any(colnames(db_df) %in% colnames(meta.features)) && isTRUE(overwrite)) {
        meta.features <- meta.features[, setdiff(colnames(meta.features), colnames(db_df))]
      }
      db_sub <- unique(db_df[rownames(db_df) %in% rownames(meta.features), , drop = FALSE])
      if (nrow(db_sub) == 0) {
        stop(paste0("No db data found in the seurat object. Please check if the species name is correct. The expected feature names are ", paste(head(db_sub[[IDtype]], 10), collapse = ","), "."))
      }
      meta.features <- cbind(meta.features, db_sub[rownames(meta.features), setdiff(colnames(db_sub), colnames(meta.features)), drop = FALSE])
      srt[[assay]]@meta.features <- meta.features
    }
  }
  return(srt)
}
