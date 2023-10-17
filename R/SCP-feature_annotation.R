#' AnnotateFeatures
#' Annotate features in a Seurat object with additional metadata from databases or a GTF file.
#' @param srt Seurat object to be annotated.
#' @param species Name of the species to be used for annotation. Default is "Homo_sapiens".
#' @param IDtype Type of identifier to use for annotation. Default is "symbol" with options "symbol", "ensembl_id", and "entrez_id".
#' @param db Vector of database names to be used for annotation. Default is NULL.
#' @param db_update Logical value indicating whether to update the database. Default is FALSE.
#' @param db_version Version of the database to use. Default is "latest".
#' @param convert_species A logical value indicating whether to use a species-converted database when the annotation is missing for the specified species. The default value is TRUE.
#' @param Ensembl_version Version of the Ensembl database to use. Default is 103.
#' @param mirror URL of the mirror to use for Ensembl database. Default is NULL.
#' @param gtf Path to the GTF file to be used for annotation. Default is NULL.
#' @param merge_gtf_by Column name to merge the GTF file by. Default is "gene_name".
#' @param columns Vector of column names to be used from the GTF file. Default is
#'     "seqname", "feature", "start", "end", "strand", "gene_id", "gene_name", "gene_type".
#' @param assays Character vector of assay names to be annotated. Default is "RNA".
#' @param overwrite Logical value indicating whether to overwrite existing metadata. Default is FALSE.
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}}
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- AnnotateFeatures(pancreas_sub,
#'   species = "Mus_musculus", IDtype = "symbol",
#'   db = c("Chromosome", "GeneType", "Enzyme", "TF", "CSPA", "VerSeDa")
#' )
#' head(pancreas_sub[["RNA"]]@meta.features)
#'
#' ## Annotate features using a GTF file
#' # pancreas_sub <- AnnotateFeatures(pancreas_sub, gtf = "/data/reference/CellRanger/refdata-gex-mm10-2020-A/genes/genes.gtf")
#' # head(pancreas_sub[["RNA"]]@meta.features)
#' @importFrom data.table fread rbindlist
#' @export
AnnotateFeatures <- function(srt, species = "Homo_sapiens", IDtype = c("symbol", "ensembl_id", "entrez_id"),
                             db = NULL, db_update = FALSE, db_version = "latest", convert_species = TRUE, Ensembl_version = 103, mirror = NULL,
                             gtf = NULL, merge_gtf_by = "gene_name", columns = c(
                               "seqname", "feature", "start", "end", "strand",
                               "gene_id", "gene_name", "gene_type"
                             ),
                             assays = "RNA", overwrite = FALSE) {
  IDtype <- match.arg(IDtype)
  if (is.null(db) && is.null(gtf)) {
    stop("Neither 'db' nor 'gtf' is specified")
  }
  if (!is.null(db)) {
    db_list <- PrepareDB(
      species = species, db = db, db_update = db_update, db_version = db_version, convert_species = convert_species,
      db_IDtypes = IDtype, Ensembl_version = Ensembl_version, mirror = mirror
    )
    db_notfound <- setdiff(db, names(db_list[[species]]))
    if (length(db_notfound) > 0) {
      warning(paste0("The following databases are not found:", paste0(db_notfound, collapse = ",")))
    }
    for (single_db in names(db_list[[species]])) {
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
          stop(paste0("No data to append was found in the Seurat object. Please check if the species name is correct. The expected feature names are ", paste(head(rownames(db_df), 10), collapse = ","), "."))
        }
        meta.features <- cbind(meta.features, db_sub[rownames(meta.features), setdiff(colnames(db_sub), colnames(meta.features)), drop = FALSE])
        srt[[assay]]@meta.features <- meta.features
      }
    }
  }

  if (!is.null(gtf)) {
    gtf_all <- suppressWarnings(fread(gtf, sep = "\t"))
    gtf_all <- gtf_all[, 1:9]
    colnames(gtf_all) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    for (type in c("gene", "transcript", "exon", "CDS")) {
      if (type %in% gtf_all[["feature"]]) {
        gtf_all <- gtf_all[gtf_all[["feature"]] == type, ]
        break
      }
    }
    columns1 <- intersect(colnames(gtf_all), columns)

    gtf_attribute <- gtf_all[["attribute"]]
    gtf_attribute <- gsub(pattern = "\"", replacement = "", x = gtf_attribute)
    gtf_attribute <- strsplit(gtf_attribute, split = "; *")
    gene_attr <- lapply(gtf_attribute, function(x) {
      detail <- strsplit(x, " ")
      out <- lapply(detail, function(x) x[2:length(x)])
      names(out) <- sapply(detail, function(x) x[1])
      out <- out[intersect(columns, names(out))]
      return(out)
    })
    gene_attr_df <- rbindlist(gene_attr, fill = TRUE)
    gtf_columns <- cbind(gtf_all[, intersect(colnames(gtf_all), columns), with = FALSE], gene_attr_df)
    colnames(gtf_columns) <- make.unique(colnames(gtf_columns))
    gtf_columns_collapse <- aggregate(gtf_columns, by = list(rowid = gtf_columns[[merge_gtf_by]]), FUN = function(x) {
      paste0(unique(x), collapse = ";")
    })
    rownames(gtf_columns_collapse) <- gtf_columns_collapse[["rowid"]]
    gtf_columns_collapse[["rowid"]] <- NULL
    for (assay in assays) {
      meta.features <- srt[[assay]]@meta.features
      if (length(intersect(colnames(meta.features), colnames(gtf_columns_collapse))) > 0 && isTRUE(overwrite)) {
        meta.features <- meta.features[, setdiff(colnames(meta.features), colnames(gtf_columns_collapse))]
      }
      meta.features <- cbind(meta.features, gtf_columns_collapse[rownames(meta.features), setdiff(colnames(gtf_columns_collapse), colnames(meta.features)), drop = FALSE])
      srt[[assay]]@meta.features <- meta.features
    }
  }
  return(srt)
}
