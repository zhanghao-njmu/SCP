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
#' pancreas_sub <- AnnotateFeatures(pancreas_sub,
#'   species = "Mus_musculus", IDtype = "symbol",
#'   db = c("GO_BP", "Chromosome", "GeneType", "Enzyme", "TF", "SP", "CellTalk")
#' )
#' head(pancreas_sub[["RNA"]]@meta.features)
#'
#' ## Annotate features using a GTF file
#' # pancreas_sub <- AnnotateFeatures(pancreas_sub, gtf = "/data/reference/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf")
#'
AnnotateFeatures <- function(srt, species = "Homo_sapiens", IDtype = c("symbol", "ensembl_id", "entrez_id"),
                             db = "TF", db_update = FALSE, db_version = "latest", convert_species = FALSE, Ensembl_version = 103, mirror = NULL,
                             gtf = NULL, merge_gtf_by = "gene_name", attribute = c("gene_id", "gene_name", "gene_type"),
                             assays = "RNA", overwrite = FALSE) {
  IDtype <- match.arg(IDtype)
  db_list <- PrepareDB(
    species = species, db = db, db_update = db_update, db_version = db_version, convert_species = convert_species,
    db_IDtypes = IDtype, Ensembl_version = Ensembl_version, mirror = mirror
  )
  for (single_db in db) {
    TERM2GENE <- unique(db_list[[species]][[single_db]][["TERM2GENE"]])
    TERM2NAME <- unique(db_list[[species]][[single_db]][["TERM2NAME"]])
    rownames(TERM2NAME) <- TERM2NAME[, 1]
    TERM2GENE[, single_db] <- TERM2NAME[TERM2GENE[, 1], 2]
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
      db_sub <- db_df[rownames(db_df) %in% rownames(meta.features), , drop = FALSE]
      if (nrow(db_sub) == 0) {
        stop(paste0("No db data found in the seurat object. Please check if the species name is correct. The expected feature names are ", paste(head(db_sub[[IDtype]], 10), collapse = ","), "."))
      }
      meta.features <- cbind(meta.features, db_sub[rownames(meta.features), setdiff(colnames(db_sub), colnames(meta.features)), drop = FALSE])
      srt[[assay]]@meta.features <- meta.features
    }
  }
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
  return(srt)
}
