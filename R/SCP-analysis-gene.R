#' Gene ID conversion function using biomart
#'
#' This function can convert different gene ID types within one species or between two species using the biomart service.
#'
#' @param geneID A vector of the geneID character.
#' @param geneID_from_IDtype Gene ID type of the input \code{geneID}. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param geneID_to_IDtype Gene ID type(s) to convert to. e.g. "symbol", "ensembl_id", "entrez_id"
#' @param species_from Latin names for animals of the input geneID. e.g. "Homo_sapiens","Mus_musculus"
#' @param species_to Latin names for animals of the output geneID. e.g. "Homo_sapiens","Mus_musculus"
#' @param Ensembl_version Ensembl database version. If NULL, use the current release version.
#' @param biomart The name of the BioMart database that you want to connect to. Possible options include "ensembl", "protists_mart", "fungi_mart", and "plants_mart".
#' @param max_tries The maximum number of attempts to connect with the BioMart service.
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'.
#'
#' @return A list with the following elements:
#'   \itemize{
#'     \item{\code{geneID_res:}}{ A data.frame contains the all gene IDs mapped in the database with columns: 'from_IDtype','from_geneID','to_IDtype','to_geneID'}
#'     \item{\code{geneID_collapse:}}{ The data.frame contains all the successfully converted gene IDs, and the output gene IDs are collapsed into a list. As a result, the 'from_geneID' column (which is set as the row names) of the data.frame is unique.}
#'     \item{\code{geneID_expand:}}{ The data.frame contains all the successfully converted gene IDs, and the output gene IDs are expanded.}
#'     \item{\code{Ensembl_version:}}{ Ensembl database version.}
#'     \item{\code{Datasets:}}{ Datasets available in the selected BioMart database.}
#'     \item{\code{Attributes:}}{ Attributes available in the selected BioMart database.}
#'     \item{\code{geneID_unmapped:}}{ A character vector of gene IDs that are unmapped in the database.}
#'   }
#'
#' @examples
#' res <- GeneConvert(
#'   geneID = c("CDK1", "MKI67", "TOP2A", "AURKA", "CTCF"),
#'   geneID_from_IDtype = "symbol",
#'   geneID_to_IDtype = "entrez_id",
#'   species_from = "Homo_sapiens",
#'   species_to = "Mus_musculus",
#'   Ensembl_version = 103
#' )
#' str(res)
#'
#' # Convert the human genes to mouse homologs and replace the raw counts in a Seurat object.
#' data("pancreas_sub")
#' counts <- pancreas_sub@assays$RNA@counts
#' res <- GeneConvert(
#'   geneID = rownames(counts),
#'   geneID_from_IDtype = "symbol",
#'   geneID_to_IDtype = "symbol",
#'   species_from = "Mus_musculus",
#'   species_to = "Homo_sapiens",
#'   Ensembl_version = 103
#' )
#' # Check the number of input and converted gene IDs
#' input_genes <- length(rownames(counts))
#' db_genes <- length(unique(res$geneID_res$from_geneID))
#' converted_genes_input <- length(unique(res$geneID_collapse$from_geneID))
#' converted_genes_output <- length(unique(res$geneID_expand$symbol))
#' cat("Number of input gene IDs:", input_genes, "\n")
#' cat("Number of gene IDs mapped in the database:", db_genes, "\n")
#' cat("Number of input gene IDs that were successfully converted:", converted_genes_input, "\n")
#' cat("Number of converted gene IDs:", converted_genes_output, "\n")
#'
#' homologs_counts <- aggregate(
#'   x = counts[res$geneID_expand[, "from_geneID"], ],
#'   by = list(res$geneID_expand[, "symbol"]), FUN = sum
#' )
#' rownames(homologs_counts) <- homologs_counts[, 1]
#' homologs_counts <- as(as_matrix(homologs_counts[, -1]), "dgCMatrix")
#' homologs_counts
#'
#' @importFrom R.cache loadCache saveCache
#' @importFrom reshape2 dcast melt
#' @importFrom biomaRt listEnsemblArchives useMart listDatasets useDataset getBM listAttributes useEnsembl
#' @export
GeneConvert <- function(geneID, geneID_from_IDtype = "symbol", geneID_to_IDtype = "entrez_id",
                        species_from = "Homo_sapiens", species_to = NULL,
                        Ensembl_version = 103, biomart = NULL, mirror = NULL, max_tries = 5) {
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
  }

  if (missing(geneID)) {
    stop("'geneID' must be provided.")
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
      "wiki_symbol" = "wikigene_name",
      x
    )
  })
  names(from_IDtype) <- geneID_from_IDtype

  to_IDtype <- sapply(geneID_to_IDtype, function(x) {
    switch(x,
      "symbol" = "external_gene_name",
      "ensembl_symbol" = "external_gene_name",
      "entrez_symbol" = "external_gene_name",
      "ensembl_id" = "ensembl_gene_id",
      "entrez_id" = "entrezgene_id",
      x
    )
  })

  if (species_from != species_to && all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
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

  if (is.null(biomart)) {
    message("Connect to the Ensembl archives...")
    archives <- try_get(
      expr = {
        listEnsemblArchives()
      },
      max_tries = max_tries,
      error_message = "Get errors when connecting with EnsemblArchives..."
    )
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
    mart <- try_get(
      expr = {
        if (!is.null(mirror)) {
          useEnsembl(biomart = "ensembl", mirror = mirror)
        } else {
          useMart(biomart = "ensembl", host = url)
        }
      },
      max_tries = max_tries,
      error_message = "Get errors when connecting with ensembl mart..."
    )
    mart_from <- mart_to <- mart
  } else {
    biomart <- match.arg(biomart, choices = c("ensembl", "protists_mart", "fungi_mart", "plants_mart"), several.ok = TRUE)
    url <- setNames(
      object = c("https://ensembl.org", "https://protists.ensembl.org", "https://fungi.ensembl.org", "https://plants.ensembl.org"),
      nm = c("ensembl", "protists_mart", "fungi_mart", "plants_mart")
    )
    message("Connecting to the biomart(", biomart, ")...")
    if (length(biomart) == 1) {
      mart <- try_get(
        expr = {
          useMart(biomart = biomart, host = url[biomart])
        },
        max_tries = max_tries,
        error_message = paste0("Get errors when connecting with ensembl mart(", biomart, ")")
      )
      mart_from <- mart_to <- mart
    } else {
      stop("Supports conversion within one mart only.")
      mart_from <- try_get(
        expr = {
          useMart(biomart = biomart[1], host = url[biomart[1]])
        },
        max_tries = max_tries,
        error_message = paste0("Get errors when connecting with ensembl mart(", biomart[1], ")")
      )
      mart_to <- try_get(
        expr = {
          useMart(biomart = biomart[2], host = url[biomart[2]])
        },
        max_tries = max_tries,
        error_message = paste0("Get errors when connecting with ensembl mart(", biomart[2], ")")
      )
    }
  }

  message("Searching the dataset ", species_from_simp, " ...")
  Datasets <- try_get(
    expr = {
      listDatasets(mart_from)
    },
    max_tries = max_tries,
    error_message = paste0("Get errors when connecting with ensembl mart(", mart_from@biomart, ")")
  )
  dataset <- searchDatasets(Datasets, pattern = species_from_simp)[["dataset"]][1]
  if (is.null(dataset)) {
    warning(paste0("Can not find the dataset for the species: ", species_from, " (", species_from_simp, ")"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = NULL))
  }

  message("Connecting to the dataset ", dataset, " ...")
  mart1 <- try_get(
    expr = {
      useDataset(dataset = dataset, mart = mart_from)
    },
    max_tries = max_tries,
    error_message = paste0("Get errors when connecting with Dataset(", dataset, ")")
  )

  if (species_from != species_to && any(!geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    message("Searching the dataset ", species_to_simp, " ...")
    Datasets2 <- try_get(
      expr = {
        listDatasets(mart_to)
      },
      max_tries = max_tries,
      error_message = paste0("Get errors when connecting with ensembl mart(", mart_to@biomart, ")")
    )
    dataset2 <- searchDatasets(Datasets2, pattern = species_to_simp)[["dataset"]][1]
    if (is.null(dataset2)) {
      warning(paste0("Can not find the dataset for the species: ", species_to, " (", species_to_simp, ")"), immediate. = TRUE)
      return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = list(Datasets, Datasets2), Attributes = NULL))
    }

    message("Connecting to the dataset ", dataset2, " ...")
    mart2 <- try_get(
      expr = {
        useDataset(dataset = dataset2, mart = mart_to)
      },
      max_tries = max_tries,
      error_message = paste0("Get errors when connecting with Dataset(", dataset2, ")")
    )
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

  message("Converting the geneIDs...")
  if (species_from != species_to) {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- try_get(
          expr = {
            getBM(
              mart = mart1,
              attributes = c(from_attr, "ensembl_gene_id"),
              filters = from_attr,
              values = list(geneID)
            )
          },
          max_tries = max_tries,
          error_message = "Get errors when retrieving information from the BioMart database"
        )

        geneID_res1 <- geneID_res1[, c(from_attr, "ensembl_gene_id"), drop = FALSE]
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, , drop = FALSE]
        if (nrow(geneID_res1) == 0) {
          next
        }
        colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id_tmp")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1[, "from_IDtype"] <- from_name

        if (all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
          geneID_res2 <- try_get(
            expr = {
              getBM(
                mart = mart1,
                attributes = unique(c("ensembl_gene_id", to_attr)),
                filters = "ensembl_gene_id",
                values = list(geneID_res1[, "ensembl_gene_id_tmp"])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res2 <- geneID_res2[, unique(c("ensembl_gene_id", to_attr))]
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], , drop = FALSE]
          if (nrow(geneID_res2) == 0) {
            next
          }
          geneID_res2[, "ensembl_gene_id_tmp"] <- geneID_res2[, "ensembl_gene_id"]
          geneID_res2 <- geneID_res2[, c("ensembl_gene_id_tmp", to_attr)]
          geneID_res2 <- unique(melt(geneID_res2, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
          geneID_res2$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res2$to_IDtype]
          geneID_res_merge <- merge(x = geneID_res1, y = geneID_res2, by = "ensembl_gene_id_tmp")
        } else {
          homolog_ensembl_gene <- paste(species_to_simp, "ensembl_gene", sep = "_homolog_")
          geneID_res2 <- try_get(
            expr = {
              getBM(
                mart = mart1,
                attributes = c("ensembl_gene_id", homolog_ensembl_gene),
                filters = "ensembl_gene_id",
                values = list(geneID_res1[, "ensembl_gene_id_tmp"])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res2 <- geneID_res2[, c("ensembl_gene_id", homolog_ensembl_gene), drop = FALSE]
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], , drop = FALSE]
          if (nrow(geneID_res2) == 0) {
            next
          }
          colnames(geneID_res2) <- c("ensembl_gene_id_tmp", homolog_ensembl_gene)

          geneID_res3 <- try_get(
            expr = {
              getBM(
                mart = mart2,
                attributes = unique(c("ensembl_gene_id", to_attr)),
                filters = "ensembl_gene_id",
                values = list(geneID_res2[, homolog_ensembl_gene])
              )
            },
            max_tries = max_tries,
            error_message = "Get errors when retrieving information from the BioMart database"
          )
          geneID_res3 <- geneID_res3[, unique(c("ensembl_gene_id", to_attr)), drop = FALSE]
          geneID_res3 <- geneID_res3[geneID_res3[, "ensembl_gene_id"] %in% geneID_res2[, homolog_ensembl_gene], , drop = FALSE]
          if (nrow(geneID_res3) == 0) {
            next
          }
          geneID_res3[, "ensembl_gene_id_tmp2"] <- geneID_res3[, "ensembl_gene_id"]
          geneID_res3 <- geneID_res3[, c("ensembl_gene_id_tmp2", to_attr)]
          geneID_res3 <- unique(melt(geneID_res3, id.vars = "ensembl_gene_id_tmp2", variable.name = "to_IDtype", value.name = "to_geneID"))
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
        geneID_res1 <- try_get(
          expr = {
            getBM(
              mart = mart1,
              attributes = unique(c("ensembl_gene_id", from_attr, to_attr)),
              filters = from_attr,
              values = list(geneID)
            )
          },
          max_tries = max_tries,
          error_message = "Get errors when retrieving information from the BioMart database"
        )
        geneID_res1 <- geneID_res1[, unique(c("ensembl_gene_id", from_attr, to_attr)), drop = FALSE]
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, , drop = FALSE]
        if (nrow(geneID_res1) == 0) {
          next
        }
        geneID_res1[, "ensembl_gene_id_tmp"] <- geneID_res1[, "ensembl_gene_id"]
        geneID_res1_from <- unique(geneID_res1[, c("ensembl_gene_id_tmp", from_attr)])
        colnames(geneID_res1_from) <- c("ensembl_gene_id_tmp", "from_geneID")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1_from[, "from_IDtype"] <- from_name

        geneID_res1_to <- unique(geneID_res1[, c("ensembl_gene_id_tmp", to_attr)])
        geneID_res1_to <- unique(melt(geneID_res1_to, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
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
  geneID_res <- unique(do.call(rbind, geneID_res_list))
  rownames(geneID_res) <- NULL

  if (is.null(geneID_res) || nrow(geneID_res) == 0 || all(is.na(geneID_res[["to_geneID"]]))) {
    warning(paste0("None of the gene IDs were converted"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Datasets = Datasets, Attributes = Attributes))
  }
  geneID_res_stat <- by(geneID_res, list(geneID_res[["from_geneID"]], geneID_res[["to_IDtype"]]), function(x) {
    data.frame(
      from_IDtype = unique(x[["from_IDtype"]]),
      from_geneID = unique(x[["from_geneID"]]),
      to_IDtype = unique(x[["to_IDtype"]]),
      to_geneID = I(list(unique(x[["to_geneID"]][!x[["to_geneID"]] %in% c("", NA)])))
    )
  })
  geneID_collapse <- do.call(rbind, geneID_res_stat)
  geneID_collapse <- unique(as.data.frame(geneID_collapse[, c("from_geneID", "to_IDtype", "to_geneID")]))
  geneID_collapse <- geneID_collapse[sapply(geneID_collapse$to_geneID, length) > 0, , drop = FALSE]
  geneID_collapse <- dcast(geneID_collapse, formula = from_geneID ~ to_IDtype, value.var = "to_geneID")
  rownames(geneID_collapse) <- geneID_collapse[, "from_geneID"]
  geneID_expand <- unnest(data = geneID_collapse, cols = colnames(geneID_collapse)[sapply(geneID_collapse, class) %in% c("list", "AsIs")], keep_empty = FALSE)

  return(list(geneID_res = geneID_res, geneID_collapse = geneID_collapse, geneID_expand = geneID_expand, Ensembl_version = version, Datasets = Datasets, Attributes = Attributes, geneID_unmapped = geneID))
}

searchDatasets <- function(datasets, pattern) {
  colIdx <- vapply(datasets, FUN = function(x) grepl(pattern = pattern, x = x, ignore.case = TRUE), FUN.VALUE = logical(length = nrow(datasets)))
  rowIdx <- apply(colIdx, 1, any)
  if (any(rowIdx)) {
    return(datasets[rowIdx, , drop = FALSE])
  } else {
    message("No matching datasets found")
    return(NULL)
  }
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
#' \code{\link{GeneConvert}}
#'
#' @examples
#' ccgenes <- CC_GenePrefetch("Homo_sapiens")
#' str(ccgenes)
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' str(ccgenes)
#' @importFrom R.cache loadCache saveCache
#' @export
CC_GenePrefetch <- function(species = "Homo_sapiens", Ensembl_version = 103, mirror = NULL, max_tries = 5, use_cached_gene = TRUE) {
  S <- Seurat::cc.genes.updated.2019$s.genes
  G2M <- Seurat::cc.genes.updated.2019$g2m.genes
  res <- NULL
  if (species != "Homo_sapiens") {
    if (isTRUE(use_cached_gene)) {
      res <- loadCache(key = list(species))
    }
    if (is.null(res)) {
      res <- GeneConvert(
        geneID = unique(c(S, G2M)),
        geneID_from_IDtype = "symbol",
        geneID_to_IDtype = "symbol",
        species_from = "Homo_sapiens",
        species_to = species,
        Ensembl_version = Ensembl_version,
        max_tries = max_tries,
        mirror = mirror
      )
      saveCache(res, key = list(species))
    } else {
      message("Using cached conversion results for ", species)
    }
    genes <- res[["geneID_collapse"]]
    S <- unlist(genes[S[S %in% rownames(genes)], "symbol"])
    G2M <- unlist(genes[G2M[G2M %in% rownames(genes)], "symbol"])
  }
  return(list(
    res = res,
    S = S,
    G2M = G2M
  ))
}

LengthCheck <- function(values, cutoff = 0) {
  return(vapply(X = values, FUN = function(x) {
    return(length(x = x) > cutoff)
  }, FUN.VALUE = logical(1)))
}

#' @importFrom Seurat UpdateSymbolList CaseMatch
#' @importFrom SeuratObject DefaultAssay GetAssayData CheckGC
#' @importFrom BiocParallel bplapply bpaggregate
#' @importFrom stats rnorm
#' @importFrom Matrix rowMeans colMeans
#' @importFrom ggplot2 cut_number
AddModuleScore2 <- function(object, slot = "data", features, pool = NULL, nbin = 24, ctrl = 100,
                            k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE,
                            BPPARAM = BiocParallel::bpparam(), ...) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object, slot = slot)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster ==
        i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning("The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          ifelse(test = search, yes = ", attempting to find updated synonyms",
            no = ", not searching for symbol synonyms"
          ),
          call. = FALSE, immediate. = TRUE
        )
        if (search) {
          tryCatch(expr = {
            updated.features <- UpdateSymbolList(symbols = missing.features, ...)
            names(x = updated.features) <- missing.features
            for (miss in names(x = updated.features)) {
              index <- which(x == miss)
              x[index] <- updated.features[miss]
            }
          }, error = function(...) {
            warning("Could not reach HGNC's gene names database",
              call. = FALSE, immediate. = TRUE
            )
          })
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning("The following features are still not present in the object: ",
              paste(missing.features, collapse = ", "),
              call. = FALSE, immediate. = TRUE
            )
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      "Could not find enough features in the object from the following feature lists:",
      paste(names(x = which(x = !LengthCheck(values = features)))),
      "Attempting to match case..."
    ))
    features <- lapply(X = features.old, FUN = CaseMatch, match = rownames(x = object))
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      "The following feature lists do not have enough features present in the object:",
      paste(names(x = which(x = !LengthCheck(values = features)))),
      "exiting..."
    ))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- rowMeans(x = assay.data[pool, , drop = FALSE])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg)) / 1e+30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)

  scores <- bplapply(1:cluster.length, function(i) {
    features.use <- features[[i]]
    # ctrl.use <- data.cut[which(data.cut %in% data.cut[features.use])]
    # ctrl.use <- names(sample(ctrl.use, size = min(ctrl * length(features.use), length(ctrl.use)), replace = FALSE))
    ctrl.use <- unlist(lapply(1:length(features.use), function(j) {
      data.cut[which(data.cut == data.cut[features.use[j]])]
    }))
    ctrl.use <- names(sample(ctrl.use, size = min(ctrl * length(features.use), length(ctrl.use)), replace = FALSE))
    ctrl.scores_i <- colMeans(x = assay.data[ctrl.use, , drop = FALSE])
    features.scores_i <- colMeans(x = assay.data[features.use, , drop = FALSE])
    return(list(ctrl.scores_i, features.scores_i))
  }, BPPARAM = BPPARAM)
  ctrl.scores <- do.call(rbind, lapply(scores, function(x) x[[1]]))
  features.scores <- do.call(rbind, lapply(scores, function(x) x[[2]]))

  # features_collapse <- bplapply(1:cluster.length, function(i) {
  #   features.use <- features[[i]]
  #   # ctrl.use <- data.cut[which(data.cut %in% data.cut[features.use])]
  #   # ctrl.use <- names(sample(ctrl.use, size = min(ctrl * length(features.use), length(ctrl.use)), replace = FALSE))
  #   ctrl.use <- unlist(lapply(1:length(features.use), function(j) {
  #     data.cut[which(data.cut == data.cut[features.use[j]])]
  #   }))
  #   ctrl.use <- names(sample(ctrl.use, size = min(ctrl * length(features.use), length(ctrl.use)), replace = FALSE))
  #   return(list(features.use, ctrl.use))
  # }, BPPARAM = BPPARAM)
  # features.scores <- bpaggregate(
  #   x = as_matrix(assay.data[unlist(lapply(features_collapse, function(x) x[[1]])), ]),
  #   by = list(unlist(lapply(1:length(features_collapse), function(x) rep(x, length(features_collapse[[x]][[1]]))))),
  #   FUN = mean,
  #   BPPARAM = BPPARAM
  # )
  # features.scores <- features.scores[,-1]
  # ctrl.scores <- bpaggregate(
  #   x = as_matrix(assay.data[unlist(lapply(features_collapse, function(x) x[[2]])), ]),
  #   by = list(unlist(lapply(1:length(features_collapse), function(x) rep(x, length(features_collapse[[x]][[2]]))))),
  #   FUN = mean,
  #   BPPARAM = BPPARAM
  # )
  # ctrl.scores <- ctrl.scores[,-1]

  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

#' CellScoring
#'
#' This function performs cell scoring on a Seurat object.
#' It calculates scores for a given set of features and adds the scores as metadata to the Seurat object.
#'
#' @inheritParams RunEnrichment
#' @param srt A Seurat object
#' @param features A named list of feature lists for scoring. If NULLL, \code{db} will be used to create features sets.
#' @param slot The slot of the Seurat object to use for scoring. Defaults to "data".
#' @param assay The assay of the Seurat object to use for scoring. Defaults to NULL, in which case the default assay of the object is used.
#' @param split.by A cell metadata variable used for splitting the Seurat object into subsets and performing scoring on each subset. Defaults to NULL.
#' @param termnames A vector of term names to be used from the database. Defaults to NULL, in which case all features from the database are used.
#' @param method The method to use for scoring. Can be "Seurat", "AUCell", or "UCell". Defaults to "Seurat".
#' @param classification Whether to perform classification based on the scores. Defaults to TRUE.
#' @param name The name of the assay to store the scores in. Only used if new_assay is TRUE. Defaults to an empty string.
#' @param new_assay Whether to create a new assay for storing the scores. Defaults to FALSE.
#' @param BPPARAM The BiocParallel parameter object. Defaults to BiocParallel::bpparam().
#' @param seed The random seed for reproducibility. Defaults to 11.
#' @param ... Additional arguments to be passed to the scoring methods.
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}}
#'
#' @examples
#' data("pancreas_sub")
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub,
#'   features = list(S = ccgenes$S, G2M = ccgenes$G2M),
#'   method = "Seurat", name = "CC"
#' )
#' CellDimPlot(pancreas_sub, "CC_classification")
#' FeatureDimPlot(pancreas_sub, "CC_G2M")
#'
#' \dontrun{
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- CellScoring(
#'   srt = panc8_sub, slot = "data", assay = "RNA",
#'   db = "GO_BP", species = "Homo_sapiens",
#'   minGSSize = 10, maxGSSize = 100,
#'   method = "Seurat", name = "GO", new_assay = TRUE
#' )
#' panc8_sub <- Integration_SCP(panc8_sub,
#'   assay = "GO",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub, slot = "data", assay = "RNA",
#'   db = "GO_BP", species = "Mus_musculus",
#'   termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
#'   method = "Seurat", name = "GO", new_assay = TRUE
#' )
#' pancreas_sub <- Standard_SCP(pancreas_sub, assay = "GO")
#' CellDimPlot(pancreas_sub, "SubCellType")
#'
#' pancreas_sub[["tech"]] <- "Mouse"
#' panc_merge <- Integration_SCP(
#'   srtList = list(panc8_sub, pancreas_sub),
#'   assay = "GO",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc_merge, group.by = c("tech", "celltype", "SubCellType", "Phase"))
#'
#' genenames <- make.unique(capitalize(rownames(panc8_sub[["RNA"]]), force_tolower = TRUE))
#' panc8_sub <- RenameFeatures(panc8_sub, newnames = genenames, assay = "RNA")
#' head(rownames(panc8_sub))
#' panc_merge <- Integration_SCP(
#'   srtList = list(panc8_sub, pancreas_sub),
#'   assay = "RNA",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc_merge, group.by = c("tech", "celltype", "SubCellType", "Phase"))
#' }
#'
#' @importFrom BiocParallel bpprogressbar<- bpRNGseed<- bpworkers
#' @importFrom Seurat AddModuleScore AddMetaData
#' @export
