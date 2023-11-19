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
FindExpressedMarkers <- function(object, ident.1 = NULL, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL,
                                 features = NULL, assay = NULL, slot = "data", min.expression = 0,
                                 test.use = "wilcox", logfc.threshold = 0.25, base = 2, pseudocount.use = 1, mean.fxn = NULL, fc.name = NULL,
                                 min.pct = 0.1, min.diff.pct = -Inf, max.cells.per.ident = Inf, latent.vars = NULL, only.pos = FALSE,
                                 min.cells.group = 3, min.cells.feature = 3,
                                 norm.method = "LogNormalize", verbose = TRUE, ...) {
  ########## FindMarkers.Seurat ##########

  assay <- assay %||% DefaultAssay(object)
  data.use <- object[[assay]]

  if (!is.null(cells.1)) {
    if (is.null(cells.2)) {
      cells.2 <- setdiff(colnames(object), cells.1)
    }
  } else {
    cells.1 <- WhichCells(object = object, idents = ident.1)
    cells.2 <- WhichCells(object = object, idents = ident.2)
  }

  # fetch latent.vars
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells.1, cells.2)
    )
  }

  ########## FindMarkers.Assay ##########

  object <- data.use

  pseudocount.use <- pseudocount.use %||% 1
  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = "counts",
    no = slot
  )
  data.use <- GetAssayData(object = object, slot = data.slot)
  data.use <- data.use[rowSums(data.use) > 0, ]
  data.use <- as_matrix(data.use)
  data.use[data.use <= min.expression] <- NA
  counts <- switch(
    EXPR = data.slot,
    "scale.data" = GetAssayData(object = object, slot = "counts"),
    numeric()
  )

  ########## FoldChange.Assay ##########
  features <- features %||% rownames(x = data.use)
  slot <- data.slot

  # By default run as if LogNormalize is done
  log1pdata.mean.fxn <- function(x) {
    return(log(x = rowMeans(x = expm1(x), na.rm = TRUE) + pseudocount.use, base = base))
  }
  scaledata.mean.fxn <- function(x) {
    return(rowMeans(x, na.rm = TRUE))
  }
  counts.mean.fxn <- function(x) {
    return(log(x = rowMeans(x = x, na.rm = TRUE) + pseudocount.use, base = base))
  }
  if (!is.null(x = norm.method)) {
    # For anything apart from log normalization set to rowMeans
    if (norm.method != "LogNormalize") {
      new.mean.fxn <- counts.mean.fxn
    } else {
      new.mean.fxn <- switch(
        EXPR = slot,
        "data" = log1pdata.mean.fxn,
        "scale.data" = scaledata.mean.fxn,
        "counts" = counts.mean.fxn,
        counts.mean.fxn
      )
    }
  } else {
    # If no normalization method is passed use slots to decide mean function
    new.mean.fxn <- switch(
      EXPR = slot,
      "data" = log1pdata.mean.fxn,
      "scale.data" = scaledata.mean.fxn,
      "counts" = counts.mean.fxn,
      log1pdata.mean.fxn
    )
  }
  mean.fxn <- mean.fxn %||% new.mean.fxn
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- fc.name %||% ifelse(
    test = slot == "scale.data",
    yes = "avg_diff",
    no = paste0("avg_log", base.text, "FC")
  )
  fc.results <- FoldChange.default(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )

  ########## FindMarkers.default ##########

  object <- data.use
  slot <- data.slot

  Seurat:::ValidateCellGroups(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )

  # reset parameters so no feature filtering is performed
  if (test.use %in% "DESeq2") {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  data <- switch(
    EXPR = slot,
    "scale.data" = counts,
    object
  )
  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  # feature selection (based on logFC)
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1] # first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }
  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }

  de.results <- PerformDE(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    verbose = verbose,
    min.cells.feature = min.cells.feature,
    latent.vars = latent.vars,
    ...
  )
  de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, 1]), ]
    de.results$p_val_adj <- p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }

  return(de.results)
}

FoldChange.default <- function(object, cells.1, cells.2, mean.fxn, fc.name, features = NULL, ...) {
  features <- features %||% rownames(x = object)
  # Calculate percent expressed
  thresh.min <- 0
  pct.1 <- round(
    x = rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min, na.rm = TRUE) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min, na.rm = TRUE) /
      length(x = cells.2),
    digits = 3
  )
  # Calculate fold change
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}


PerformDE <- function(object, cells.1, cells.2, features, test.use, verbose, min.cells.feature, latent.vars, ...) {
  if (!(test.use %in% c("negbinom", "poisson", "MAST", "LR")) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for the following tests: ",
      paste(c("negbinom", "poisson", "MAST", "LR"), collapse = ", "),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
  data.use <- as_matrix(data.use)

  de.results <- switch(
    EXPR = test.use,
    "wilcox" = WilcoxDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "bimod" = Seurat:::DiffExpTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "roc" = Seurat:::MarkerTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "t" = Seurat:::DiffTTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose
    ),
    "negbinom" = Seurat:::GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "poisson" = Seurat:::GLMDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      min.cells = min.cells.feature,
      latent.vars = latent.vars,
      test.use = test.use,
      verbose = verbose
    ),
    "MAST" = Seurat:::MASTDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = Seurat:::DESeq2DETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      verbose = verbose,
      ...
    ),
    "LR" = Seurat:::LRDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    ),
    stop("Unknown test: ", test.use)
  )
  return(de.results)
}

#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
WilcoxDETest <- function(data.use, cells.1, cells.2, verbose = TRUE, ...) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  check_R("limma")
  p_val <- my.sapply(
    X = 1:nrow(x = data.use),
    FUN = function(x) {
      keep <- colnames(data.use)[!is.na(data.use[x, ])]
      j <- seq_len(length.out = length(x = intersect(cells.1, keep)))
      statistics <- data.use[x, keep]
      return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = statistics)), 1))
    }
  )

  return(data.frame(p_val, row.names = rownames(x = data.use)))
}


#' Differential gene test
#'
#' This function utilizes the Seurat package to perform a differential expression (DE) test on gene expression data.
#' Users have the flexibility to specify custom cell groups, marker types, and various options for DE analysis.
#'
#' @inheritParams  Seurat::FindMarkers
#' @param srt A Seurat object.
#' @param group_by A grouping variable in the dataset to define the groups or conditions for the differential test. If not provided, the function uses the "active.ident" variable in the Seurat object.
#' @param group1 A vector of cell IDs or a character vector specifying the cells that belong to the first group. If both group_by and group1 are provided, group1 takes precedence.
#' @param group2 A vector of cell IDs or a character vector specifying the cells that belong to the second group. This parameter is only used when group_by or group1 is provided.
#' @param cells1 A vector of cell IDs specifying the cells that belong to group1. If provided, group1 is ignored.
#' @param cells2 A vector of cell IDs specifying the cells that belong to group2. This parameter is only used when cells1 is provided.
#' @param features A vector of gene names specifying the features to consider for the differential test. If not provided, all features in the dataset are considered.
#' @param markers_type A character value specifying the type of markers to find. Possible values are "all", "paired", "conserved", and "disturbed".
#' @param grouping.var A character value specifying the grouping variable for finding conserved or disturbed markers. This parameter is only used when markers_type is "conserved" or "disturbed".
#' @param fc.threshold A numeric value used to filter genes for testing based on their average fold change between/among the two groups. By default, it is set to 1.5
#' @param meta.method A character value specifying the method to use for combining p-values in the conserved markers test. Possible values are "maximump", "minimump", "wilkinsonp", "meanp", "sump", and "votep".
#' @param norm.method Normalization method for fold change calculation when slot is 'data'. Default is "LogNormalize".
#' @param p.adjust.method A character value specifying the method to use for adjusting p-values. Default is "bonferroni".
#' @param BPPARAM A BiocParallelParam object specifying the parallelization parameters for the differential test. Default is BiocParallel::bpparam().
#' @param seed An integer value specifying the seed. Default is 11.
#' @param verbose A logical value indicating whether to display progress messages during the differential test. Default is TRUE.
#' @param ... Additional arguments to pass to the \code{\link[Seurat]{FindMarkers}} function.
#'
#' @seealso \code{\link{RunEnrichment}} \code{\link{RunGSEA}} \code{\link{GroupHeatmap}}
#' @examples
#' library(dplyr)
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "SubCellType")
#' AllMarkers <- filter(pancreas_sub@tools$DEtest_SubCellType$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' table(AllMarkers$group1)
#' ht1 <- GroupHeatmap(pancreas_sub, features = AllMarkers$gene, feature_split = AllMarkers$group1, group.by = "SubCellType")
#' ht1$plot
#'
#' TopMarkers <- AllMarkers %>%
#'   group_by(gene) %>%
#'   top_n(1, avg_log2FC) %>%
#'   group_by(group1) %>%
#'   top_n(3, avg_log2FC)
#' ht2 <- GroupHeatmap(pancreas_sub, features = TopMarkers$gene, feature_split = TopMarkers$group1, group.by = "SubCellType", show_row_names = TRUE)
#' ht2$plot
#'
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "SubCellType", markers_type = "paired")
#' PairedMarkers <- filter(pancreas_sub@tools$DEtest_SubCellType$PairedMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' table(PairedMarkers$group1)
#' ht3 <- GroupHeatmap(pancreas_sub, features = PairedMarkers$gene, feature_split = PairedMarkers$group1, group.by = "SubCellType")
#' ht3$plot
#'
#' data("panc8_sub")
#' panc8_sub <- Integration_SCP(panc8_sub, batch = "tech", integration_method = "Seurat")
#' CellDimPlot(panc8_sub, group.by = c("celltype", "tech"))
#'
#' panc8_sub <- RunDEtest(panc8_sub, group_by = "celltype", grouping.var = "tech", markers_type = "conserved")
#' ConservedMarkers1 <- filter(panc8_sub@tools$DEtest_celltype$ConservedMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' table(ConservedMarkers1$group1)
#' ht4 <- GroupHeatmap(panc8_sub,
#'   slot = "data",
#'   features = ConservedMarkers1$gene, feature_split = ConservedMarkers1$group1,
#'   group.by = "tech", split.by = "celltype", within_groups = TRUE
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(panc8_sub, group_by = "tech", grouping.var = "celltype", markers_type = "conserved")
#' ConservedMarkers2 <- filter(panc8_sub@tools$DEtest_tech$ConservedMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1)
#' table(ConservedMarkers2$group1)
#' ht4 <- GroupHeatmap(panc8_sub,
#'   slot = "data",
#'   features = ConservedMarkers2$gene, feature_split = ConservedMarkers2$group1,
#'   group.by = "tech", split.by = "celltype"
#' )
#' ht4$plot
#'
#' panc8_sub <- RunDEtest(panc8_sub, group_by = "celltype", grouping.var = "tech", markers_type = "disturbed")
#' DisturbedMarkers <- filter(panc8_sub@tools$DEtest_celltype$DisturbedMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 1 & var1 == "smartseq2")
#' table(DisturbedMarkers$group1)
#' ht5 <- GroupHeatmap(panc8_sub,
#'   slot = "data",
#'   features = DisturbedMarkers$gene, feature_split = DisturbedMarkers$group1,
#'   group.by = "celltype", split.by = "tech"
#' )
#' ht5$plot
#'
#' gene_specific <- names(which(table(DisturbedMarkers$gene) == 1))
#' DisturbedMarkers_specific <- DisturbedMarkers[DisturbedMarkers$gene %in% gene_specific, ]
#' table(DisturbedMarkers_specific$group1)
#' ht6 <- GroupHeatmap(panc8_sub,
#'   slot = "data",
#'   features = DisturbedMarkers_specific$gene, feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype", split.by = "tech"
#' )
#' ht6$plot
#'
#' ht7 <- GroupHeatmap(panc8_sub,
#'   slot = "data", aggregate_fun = function(x) mean(expm1(x)) + 1,
#'   features = DisturbedMarkers_specific$gene, feature_split = DisturbedMarkers_specific$group1,
#'   group.by = "celltype", grouping.var = "tech", numerator = "smartseq2"
#' )
#' ht7$plot
#'
#' @importFrom BiocParallel bplapply SerialParam bpprogressbar<- bpRNGseed<- bpworkers
#' @importFrom Seurat FindMarkers Assays Idents
#' @importFrom stats p.adjust
#' @export
#'
RunDEtest <- function(srt, group_by = NULL, group1 = NULL, group2 = NULL, cells1 = NULL, cells2 = NULL, features = NULL,
                      markers_type = c("all", "paired", "conserved", "disturbed"),
                      grouping.var = NULL, meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
                      test.use = "wilcox", only.pos = TRUE, fc.threshold = 1.5, base = 2, pseudocount.use = 1, mean.fxn = NULL,
                      min.pct = 0.1, min.diff.pct = -Inf, max.cells.per.ident = Inf, latent.vars = NULL,
                      min.cells.feature = 3, min.cells.group = 3,
                      norm.method = "LogNormalize", p.adjust.method = "bonferroni", slot = "data", assay = NULL,
                      BPPARAM = BiocParallel::bpparam(), seed = 11, verbose = TRUE, ...) {
  set.seed(seed)
  markers_type <- match.arg(markers_type)
  meta.method <- match.arg(meta.method)
  if (markers_type %in% c("conserved", "disturbed")) {
    if (is.null(grouping.var)) {
      stop("'grouping.var' must be provided when finding conserved or disturbed markers")
    }
  }
  assay <- assay %||% DefaultAssay(srt)

  status <- check_DataType(srt, slot = slot, assay = assay)
  if (slot == "counts" && status != "raw_counts") {
    stop("Data in the 'counts' slot is not raw counts.")
  }
  if (slot == "data" && status != "log_normalized_counts") {
    if (status == "raw_counts") {
      warning("Data in the 'data' slot is raw counts. Perform NormalizeData(LogNormalize) on the data.", immediate. = TRUE)
      srt <- NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (status == "raw_normalized_counts") {
      warning("Data in the 'data' slot is raw_normalized_counts. Perform NormalizeData(LogNormalize) on the data.", immediate. = TRUE)
      srt <- NormalizeData(object = srt, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
    }
    if (status == "unknown") {
      warning("Data in the 'data' slot is unknown. Please check the data type.")
    }
  }
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed

  time_start <- Sys.time()
  if (verbose) {
    message(paste0("[", time_start, "] ", "Start DEtest"))
    message("Workers: ", bpworkers(BPPARAM))
  }

  if (fc.threshold < 1) {
    stop("fc.threshold must be greater than or equal to 1")
  }

  if (!is.null(cells1) || !is.null(group1)) {
    if (is.null(cells1)) {
      if (is.null(group_by)) {
        stop("'group_by' must be provided when 'group1' specified")
      }
      cells1 <- colnames(srt)[srt[[group_by, drop = TRUE]] %in% group1]
    }
    if (is.null(cells2) && !is.null(group2)) {
      cells2 <- colnames(srt)[srt[[group_by, drop = TRUE]] %in% group2]
    }
    if (!all(cells1 %in% colnames(srt))) {
      stop("cells1 has some cells not in the Seurat object.")
    }
    if (is.null(cells2)) {
      cells2 <- colnames(srt)[!colnames(srt) %in% cells1]
      group2 <- "others"
    }
    if (!all(cells2 %in% colnames(srt))) {
      stop("cells2 has some cells not in the Seurat object.")
    }
    if (length(cells1) < 3 || length(cells2) < 3) {
      stop("Cell groups must have more than 3 cells")
    }
    if (verbose) {
      message("Find ", markers_type, " markers(", test.use, ") for custom cell groups...")
    }

    if (markers_type == "all") {
      markers <- FindMarkers(
        object = Assays(srt, assay), slot = slot,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        markers[, "group1"] <- group1 %||% "group1"
        markers[, "group2"] <- group2 %||% "group2"
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(markers[, "group1"], levels = unique(markers[, "group1"]))
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = p.adjust.method)
        }
        markers[, "test_group_number"] <- as.integer(table(markers[["gene"]])[markers[, "gene"]])
        MarkersMatrix <- as.data.frame.matrix(table(markers[, c("gene", "group1")]))
        markers[, "test_group"] <- apply(MarkersMatrix, 1, function(x) {
          paste0(colnames(MarkersMatrix)[x > 0], collapse = ";")
        })[markers[, "gene"]]
        srt@tools[["DEtest_custom"]][[paste0("AllMarkers_", test.use)]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        warning("No markers found.", immediate. = TRUE)
      }
    }
    if (markers_type == "conserved") {
      markers <- FindConservedMarkers2(
        object = srt, assay = assay, slot = slot,
        cells.1 = cells1,
        cells.2 = cells2,
        features = features,
        grouping.var = grouping.var,
        test.use = test.use,
        logfc.threshold = log(fc.threshold, base = base),
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        meta.method = meta.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        verbose = FALSE,
        ...
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        markers[, "gene"] <- rownames(markers)
        markers[, "group1"] <- group1 %||% "group1"
        markers[, "group2"] <- group2 %||% "group2"
        rownames(markers) <- NULL
        markers[, "group1"] <- factor(markers[, "group1"], levels = unique(markers[, "group1"]))
        if ("p_val" %in% colnames(markers)) {
          markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = p.adjust.method)
        }
        markers[, "test_group_number"] <- as.integer(table(markers[["gene"]])[markers[, "gene"]])
        MarkersMatrix <- as.data.frame.matrix(table(markers[, c("gene", "group1")]))
        markers[, "test_group"] <- apply(MarkersMatrix, 1, function(x) {
          paste0(colnames(MarkersMatrix)[x > 0], collapse = ";")
        })[markers[, "gene"]]
        srt@tools[["DEtest_custom"]][[paste0("ConservedMarkers_", test.use)]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
        srt@tools[["DEtest_custom"]][["cells2"]] <- cells2
      } else {
        warning("No markers found.", immediate. = TRUE)
      }
    }
    if (markers_type == "disturbed") {
      srt_tmp <- srt
      srt_tmp[[grouping.var, drop = TRUE]][setdiff(colnames(srt_tmp), cells1)] <- NA
      bpprogressbar(BPPARAM) <- FALSE
      srt_tmp <- RunDEtest(
        srt = srt_tmp, assay = assay, slot = slot,
        group_by = grouping.var,
        markers_type = "all",
        features = features,
        test.use = test.use,
        fc.threshold = fc.threshold,
        base = base,
        min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars,
        only.pos = only.pos,
        norm.method = norm.method,
        p.adjust.method = p.adjust.method,
        pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn,
        BPPARAM = BPPARAM,
        seed = seed,
        verbose = FALSE,
        ...
      )
      markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0("AllMarkers_", test.use)]]
      if (!is.null(markers) && nrow(markers) > 0) {
        colnames(markers) <- gsub("group", "var", colnames(markers))
        markers[["group1"]] <- group1 %||% "group1"
        srt@tools[["DEtest_custom"]][[paste0("DisturbedMarkers_", test.use)]] <- markers
        srt@tools[["DEtest_custom"]][["cells1"]] <- cells1
      } else {
        warning("No markers found.", immediate. = TRUE)
      }
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
    cell_group <- setNames(unlist(lapply(cell_group, function(x) x), use.names = FALSE), unlist(lapply(cell_group, names)))

    args1 <- list(
      object = Assays(srt, assay),
      slot = slot,
      features = features,
      test.use = test.use,
      logfc.threshold = log(fc.threshold, base = base),
      base = base,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      max.cells.per.ident = max.cells.per.ident,
      min.cells.feature = min.cells.feature,
      min.cells.group = min.cells.group,
      latent.vars = latent.vars,
      only.pos = only.pos,
      norm.method = norm.method,
      pseudocount.use = pseudocount.use,
      mean.fxn = mean.fxn,
      verbose = FALSE,
      ...
    )

    if (verbose) {
      message("Find ", markers_type, " markers(", test.use, ") among ", nlevels(cell_group), " groups...")
    }
    if (markers_type == "all") {
      AllMarkers <- bplapply(levels(cell_group), FUN = function(group) {
        cells.1 <- names(cell_group)[which(cell_group == group)]
        cells.2 <- names(cell_group)[which(cell_group != group)]
        # print(paste0(group," cells.1:",length(cells.1)," cells.2:",length(cells.2)))
        if (length(cells.1) < 3 || length(cells.2) < 3) {
          return(NULL)
        } else {
          args1[["cells.1"]] <- cells.1
          args1[["cells.2"]] <- cells.2
          markers <- do.call(FindMarkers, args1)
          if (!is.null(markers) && nrow(markers) > 0) {
            markers[, "gene"] <- rownames(markers)
            markers[, "group1"] <- as.character(group)
            markers[, "group2"] <- "others"
            if ("p_val" %in% colnames(markers)) {
              markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = p.adjust.method)
            }
            return(markers)
          } else {
            return(NULL)
          }
        }
      }, BPPARAM = BPPARAM)
      AllMarkers <- do.call(rbind.data.frame, AllMarkers)
      if (!is.null(AllMarkers) && nrow(AllMarkers) > 0) {
        rownames(AllMarkers) <- NULL
        AllMarkers[, "group1"] <- factor(AllMarkers[, "group1"], levels = levels(cell_group))
        AllMarkers[, "test_group_number"] <- as.integer(table(AllMarkers[["gene"]])[AllMarkers[, "gene"]])
        AllMarkersMatrix <- as.data.frame.matrix(table(AllMarkers[, c("gene", "group1")]))
        AllMarkers[, "test_group"] <- apply(AllMarkersMatrix, 1, function(x) {
          paste0(colnames(AllMarkersMatrix)[x > 0], collapse = ";")
        })[AllMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("AllMarkers_", test.use)]] <- AllMarkers
      } else {
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("AllMarkers_", test.use)]] <- data.frame()
      }
    }

    if (markers_type == "paired") {
      pair <- expand.grid(x = levels(cell_group), y = levels(cell_group))
      pair <- pair[pair[, 1] != pair[, 2], , drop = FALSE]
      PairedMarkers <- bplapply(seq_len(nrow(pair)), function(i) {
        cells.1 <- names(cell_group)[which(cell_group == pair[i, 1])]
        cells.2 <- names(cell_group)[which(cell_group == pair[i, 2])]
        if (length(cells.1) < 3 || length(cells.2) < 3) {
          return(NULL)
        } else {
          args1[["cells.1"]] <- cells.1
          args1[["cells.2"]] <- cells.2
          markers <- do.call(FindMarkers, args1)
          if (!is.null(markers) && nrow(markers) > 0) {
            markers[, "gene"] <- rownames(markers)
            markers[, "group1"] <- as.character(pair[i, 1])
            markers[, "group2"] <- as.character(pair[i, 2])
            if ("p_val" %in% colnames(markers)) {
              markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = p.adjust.method)
            }
            return(markers)
          } else {
            return(NULL)
          }
        }
      }, BPPARAM = BPPARAM)
      PairedMarkers <- do.call(rbind.data.frame, PairedMarkers)
      if (!is.null(PairedMarkers) && nrow(PairedMarkers) > 0) {
        rownames(PairedMarkers) <- NULL
        PairedMarkers[, "group1"] <- factor(PairedMarkers[, "group1"], levels = levels(cell_group))
        PairedMarkers[, "test_group_number"] <- as.integer(table(PairedMarkers[["gene"]])[PairedMarkers[, "gene"]])
        PairedMarkersMatrix <- as.data.frame.matrix(table(PairedMarkers[, c("gene", "group1")]))
        PairedMarkers[, "test_group"] <- apply(PairedMarkersMatrix, 1, function(x) {
          paste0(colnames(PairedMarkersMatrix)[x > 0], collapse = ";")
        })[PairedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkers_", test.use)]] <- PairedMarkers
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkersMatrix_", test.use)]] <- PairedMarkersMatrix
      } else {
        warning("No markers found.", immediate. = TRUE)
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkers_", test.use)]] <- data.frame()
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("PairedMarkersMatrix_", test.use)]] <- NULL
      }
    }

    if (markers_type == "conserved") {
      ConservedMarkers <- bplapply(levels(cell_group), FUN = function(group) {
        cells.1 <- names(cell_group)[which(cell_group == group)]
        cells.2 <- names(cell_group)[which(cell_group != group)]
        # print(paste0(group," cells.1:",length(cells.1)," cells.2:",length(cells.2)))
        if (length(cells.1) < 3 || length(cells.2) < 3) {
          return(NULL)
        } else {
          args1[["cells.1"]] <- cells.1
          args1[["cells.2"]] <- cells.2
          args1[["object"]] <- srt
          args1[["assay"]] <- assay
          args1[["grouping.var"]] <- grouping.var
          args1[["meta.method"]] <- meta.method
          markers <- do.call(FindConservedMarkers2, args1)
          if (!is.null(markers) && nrow(markers) > 0) {
            markers[, "gene"] <- rownames(markers)
            markers[, "group1"] <- as.character(group)
            markers[, "group2"] <- "others"
            if ("p_val" %in% colnames(markers)) {
              markers[, "p_val_adj"] <- p.adjust(markers[, "p_val"], method = p.adjust.method)
            }
            return(markers)
          } else {
            return(NULL)
          }
        }
      }, BPPARAM = BPPARAM)
      ConservedMarkers <- do.call(rbind.data.frame, lapply(ConservedMarkers, function(x) x[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj", "gene", "group1", "group2")]))
      if (!is.null(ConservedMarkers) && nrow(ConservedMarkers) > 0) {
        rownames(ConservedMarkers) <- NULL
        ConservedMarkers[, "group1"] <- factor(ConservedMarkers[, "group1"], levels = levels(cell_group))
        ConservedMarkers[, "test_group_number"] <- as.integer(table(ConservedMarkers[["gene"]])[ConservedMarkers[, "gene"]])
        ConservedMarkersMatrix <- as.data.frame.matrix(table(ConservedMarkers[, c("gene", "group1")]))
        ConservedMarkers[, "test_group"] <- apply(ConservedMarkersMatrix, 1, function(x) {
          paste0(colnames(ConservedMarkersMatrix)[x > 0], collapse = ";")
        })[ConservedMarkers[, "gene"]]
        ConservedMarkers <- ConservedMarkers[, c("avg_log2FC", "pct.1", "pct.2", "max_pval", "p_val", "p_val_adj", "gene", "group1", "group2", "test_group_number", "test_group")]
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("ConservedMarkers_", test.use)]] <- ConservedMarkers
      } else {
        warning("No markers found.", immediate. = TRUE)
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("ConservedMarkers_", test.use)]] <- data.frame()
      }
    }

    if (markers_type == "disturbed") {
      sub_BPPARAM <- SerialParam()
      bpprogressbar(sub_BPPARAM) <- FALSE
      DisturbedMarkers <- bplapply(levels(cell_group), FUN = function(group) {
        cells.1 <- names(cell_group)[which(cell_group == group)]
        srt_tmp <- srt
        srt_tmp[[grouping.var, drop = TRUE]][setdiff(colnames(srt_tmp), cells.1)] <- NA
        if (length(na.omit(unique(srt_tmp[[grouping.var, drop = TRUE]]))) < 2) {
          return(NULL)
        } else {
          srt_tmp <- RunDEtest(
            srt = srt_tmp, assay = assay, slot = slot,
            group_by = grouping.var,
            markers_type = "all",
            features = features,
            test.use = test.use,
            fc.threshold = fc.threshold,
            base = base,
            min.pct = min.pct,
            min.diff.pct = min.diff.pct,
            max.cells.per.ident = max.cells.per.ident,
            min.cells.feature = min.cells.feature,
            min.cells.group = min.cells.group,
            latent.vars = latent.vars,
            only.pos = only.pos,
            norm.method = norm.method,
            p.adjust.method = p.adjust.method,
            pseudocount.use = pseudocount.use,
            mean.fxn = mean.fxn,
            BPPARAM = sub_BPPARAM,
            seed = seed,
            verbose = FALSE,
            ...
          )
          markers <- srt_tmp@tools[[paste0("DEtest_", grouping.var)]][[paste0("AllMarkers_", test.use)]]
          if (!is.null(markers) && nrow(markers) > 0) {
            colnames(markers) <- gsub("group", "var", colnames(markers))
            markers[["group1"]] <- as.character(group)
            return(markers)
          } else {
            return(NULL)
          }
        }
      }, BPPARAM = BPPARAM)
      DisturbedMarkers <- do.call(rbind.data.frame, DisturbedMarkers)
      if (!is.null(DisturbedMarkers) && nrow(DisturbedMarkers) > 0) {
        rownames(DisturbedMarkers) <- NULL
        DisturbedMarkers[, "group1"] <- factor(DisturbedMarkers[, "group1"], levels = levels(cell_group))
        DisturbedMarkers[, "test_group_number"] <- as.integer(table(unique(DisturbedMarkers[, c("gene", "group1")])[["gene"]])[DisturbedMarkers[, "gene"]])
        DisturbedMarkersMatrix <- as.data.frame.matrix(table(DisturbedMarkers[, c("gene", "group1")]))
        DisturbedMarkers[, "test_group"] <- apply(DisturbedMarkersMatrix, 1, function(x) {
          paste0(colnames(DisturbedMarkersMatrix)[x > 0], collapse = ";")
        })[DisturbedMarkers[, "gene"]]
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("DisturbedMarkers_", test.use)]] <- DisturbedMarkers
      } else {
        warning("No markers found.", immediate. = TRUE)
        srt@tools[[paste0("DEtest_", group_by)]][[paste0("DisturbedMarkers_", test.use)]] <- data.frame()
      }
    }
  }
  time_end <- Sys.time()
  if (verbose) {
    message(paste0("[", time_end, "] ", "DEtest done"))
    message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))
  }
  return(srt)
}

#' ListDB
#'
#' Retrieves information about databases based on a given species and database name.
#'
#' @param species The species for which to retrieve database information. Default is "Homo_sapiens".
#' @param db The pattern to match against the database names. Default is NULL, which matches all databases.
#'
#' @seealso \code{\link{PrepareDB}}
#'
#' @importFrom R.cache getCacheRootPath readCacheHeader
#' @export
#' @examples
#' ListDB(species = "Homo_sapiens")
#' ListDB(species = "Mus_musculus", db = "GO_BP")
#'
ListDB <- function(species = "Homo_sapiens",
                   db = NULL) {
  stopifnot(length(species) == 1)
  pathnames <- dir(path = getCacheRootPath(), pattern = "[.]Rcache$", full.names = TRUE)
  if (length(pathnames) == 0) {
    return(NULL)
  }
  dbinfo <- lapply(pathnames, function(x) {
    info <- readCacheHeader(x)
    info[["date"]] <- as.character(info[["timestamp"]])
    info[["db_version"]] <- strsplit(info[["comment"]], "\\|")[[1]][1]
    info[["db_name"]] <- strsplit(info[["comment"]], "\\|")[[1]][2]
    return(info)
  })
  dbinfo <- do.call(rbind.data.frame, dbinfo)
  dbinfo[["file"]] <- pathnames

  if (is.null(db)) {
    db <- ".*"
  }
  patterns <- paste0("^", species, "-", db, "$")
  dbinfo <- dbinfo[unlist(lapply(patterns, function(pat) grep(pat, dbinfo[["db_name"]]))), , drop = FALSE]
  dbinfo <- dbinfo[order(dbinfo[["timestamp"]], decreasing = TRUE), , drop = FALSE]
  rownames(dbinfo) <- NULL
  return(dbinfo)
}

#' Prepare the gene annotation databases
#'
#' This function prepares the gene annotation databases for a given species and set of annotation sources.
#' It retrieves the necessary information from various annotation packages or external resources and organizes it into a list.
#' The list contains the annotation data for each specified annotation source.
#'
#' @inheritParams GeneConvert
#' @param species A character vector specifying the species for which the gene annotation databases should be prepared. Default is c("Homo_sapiens", "Mus_musculus").
#' @param db A character vector specifying the annotation sources to be included in the gene annotation databases.
#' Default is c("GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome",
#' "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa",
#' "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB",
#' "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF").
#' @param db_IDtypes A character vector specifying the desired ID types to be used for gene identifiers in the gene annotation databases. Default is c("symbol", "entrez_id", "ensembl_id").
#' @param db_version A character vector specifying the version of the gene annotation databases to be retrieved. Default is "latest".
#' @param db_update A logical value indicating whether the gene annotation databases should be forcefully updated. If set to FALSE, the function will attempt to load the cached databases instead. Default is FALSE.
#' @param convert_species A logical value indicating whether to use a species-converted database when the annotation is missing for the specified species. The default value is TRUE.
#' @param custom_TERM2GENE A data frame containing a custom TERM2GENE mapping for the specified species and annotation source. Default is NULL.
#' @param custom_TERM2NAME A data frame containing a custom TERM2NAME mapping for the specified species and annotation source. Default is NULL.
#' @param custom_species A character vector specifying the species name to be used in a custom database. Default is NULL.
#' @param custom_IDtype A character vector specifying the ID type to be used in a custom database. Default is NULL.
#' @param custom_version A character vector specifying the version to be used in a custom database. Default is NULL.
#'
#' @details
#' The `PrepareDB` function prepares gene annotation databases for a given species and set of annotation sources.
#' It retrieves the necessary information from various annotation packages or external resources and organizes it into a list.
#' The function also supports creating custom databases based on user-provided gene sets.
#'
#' @return A list containing the prepared gene annotation databases:
#'   \itemize{
#'     \item{\code{TERM2GENE:}}{ mapping of gene identifiers to terms}
#'     \item{\code{TERM2NAME:}}{ mapping of terms to their names}
#'     \item{\code{semData:}}{ semantic similarity data for gene sets (only for Gene Ontology terms)}
#'     }
#'
#' @seealso \code{\link{ListDB}}
#'
#' @examples
#' if (interactive()) {
#'   db_list <- PrepareDB(species = "Homo_sapiens", db = "GO_BP")
#'   ListDB(species = "Homo_sapiens", db = "GO_BP")
#'   head(db_list[["Homo_sapiens"]][["GO_BP"]][["TERM2GENE"]])
#'
#'   # Based on homologous gene conversion, prepare a gene annotation database that originally does not exist in the species.
#'   db_list <- PrepareDB(species = "Homo_sapiens", db = "MP")
#'   ListDB(species = "Homo_sapiens", db = "MP")
#'   head(db_list[["Homo_sapiens"]][["MP"]][["TERM2GENE"]])
#'
#'   # Prepare databases for other species
#'   db_list <- PrepareDB(species = "Macaca_fascicularis", db = "GO_BP")
#'   ListDB(species = "Macaca_fascicularis", db = "GO_BP")
#'   head(db_list[["Macaca_fascicularis"]][["GO_BP"]][["TERM2GENE"]])
#'
#'   db_list <- PrepareDB(species = "Saccharomyces_cerevisiae", db = "GO_BP")
#'   ListDB(species = "Saccharomyces_cerevisiae", db = "GO_BP")
#'   head(db_list[["Saccharomyces_cerevisiae"]][["GO_BP"]][["TERM2GENE"]])
#'
#'   # Prepare databases for Arabidopsis (plant)
#'   db_list <- PrepareDB(
#'     species = "Arabidopsis_thaliana",
#'     db = c(
#'       "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway",
#'       "ENZYME", "Chromosome"
#'     ),
#'     biomart = "plants_mart"
#'   )
#'   head(db_list[["Arabidopsis_thaliana"]][["KEGG"]][["TERM2GENE"]])
#'
#'   # You can also build a custom database based on the gene sets you have
#'   ccgenes <- CC_GenePrefetch("Homo_sapiens")
#'   custom_TERM2GENE <- rbind(
#'     data.frame(term = "S_genes", gene = ccgenes[["cc_S_genes"]]),
#'     data.frame(term = "G2M_genes", gene = ccgenes[["cc_G2M_genes"]])
#'   )
#'   str(custom_TERM2GENE)
#'
#'   # Set convert_species = TRUE to build a custom database for both species, with the name "CellCycle"
#'   db_list <- PrepareDB(
#'     species = c("Homo_sapiens", "Mus_musculus"), db = "CellCycle", convert_species = TRUE,
#'     custom_TERM2GENE = custom_TERM2GENE, custom_species = "Homo_sapiens", custom_IDtype = "symbol", custom_version = "Seurat_v4"
#'   )
#'   ListDB(db = "CellCycle")
#'
#'   db_list <- PrepareDB(species = "Mus_musculus", db = "CellCycle")
#'   head(db_list[["Mus_musculus"]][["CellCycle"]][["TERM2GENE"]])
#' }
#' @importFrom R.cache loadCache saveCache readCacheHeader
#' @importFrom utils packageVersion read.table
#' @importFrom stats na.omit
#' @importFrom AnnotationDbi select keys columns mappedkeys
#' @importFrom clusterProfiler read.gmt
#' @export
#'
PrepareDB <- function(species = c("Homo_sapiens", "Mus_musculus"),
                      db = c(
                        "GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome",
                        "CORUM", "MP", "DO", "HPO", "PFAM",
                        "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa",
                        "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE",
                        "MSigDB", "CellTalk", "CellChat",
                        "Chromosome", "GeneType", "Enzyme", "TF"
                      ),
                      db_IDtypes = c("symbol", "entrez_id", "ensembl_id"),
                      db_version = "latest", db_update = FALSE,
                      convert_species = TRUE,
                      Ensembl_version = 103, mirror = NULL, biomart = NULL, max_tries = 5,
                      custom_TERM2GENE = NULL, custom_TERM2NAME = NULL,
                      custom_species = NULL, custom_IDtype = NULL, custom_version = NULL) {
  db_list <- list()
  for (sps in species) {
    message("Species: ", sps)
    default_IDtypes <- list(
      "GO" = "entrez_id", "GO_BP" = "entrez_id", "GO_CC" = "entrez_id", "GO_MF" = "entrez_id", "KEGG" = "entrez_id",
      "WikiPathway" = "entrez_id", "Reactome" = "entrez_id", "CORUM" = "symbol",
      "MP" = "symbol", "DO" = "symbol", "HPO" = "symbol", "PFAM" = "entrez_id", "Chromosome" = "entrez_id",
      "GeneType" = "entrez_id", "Enzyme" = "entrez_id", "TF" = "symbol",
      "CSPA" = "symbol", "Surfaceome" = "symbol", "SPRomeDB" = "entrez_id", "VerSeDa" = "symbol",
      "TFLink" = "symbol", "hTFtarget" = "symbol", "TRRUST" = "symbol", "JASPAR" = "symbol", "ENCODE" = "symbol",
      "MSigDB" = c("symbol", "ensembl_id"), "CellTalk" = "symbol", "CellChat" = "symbol"
    )
    if (!is.null(custom_TERM2GENE)) {
      if (length(db) > 1) {
        stop("When building a custom database, the length of 'db' must be 1.")
      }
      if (is.null(custom_IDtype) || is.null(custom_species) || is.null(custom_version)) {
        stop("When building a custom database, 'custom_IDtype', 'custom_species' and 'custom_version' must be provided.")
      }
      custom_IDtype <- match.arg(custom_IDtype, choices = c("symbol", "entrez_id", "ensembl_id"))
      default_IDtypes[[db]] <- custom_IDtype
    }
    if (isFALSE(db_update) && is.null(custom_TERM2GENE)) {
      for (term in db) {
        # Try to load cached database, if already generated.
        dbinfo <- ListDB(species = sps, db = term)
        if (nrow(dbinfo) > 0 && !is.null(dbinfo)) {
          if (db_version == "latest") {
            pathname <- dbinfo[order(dbinfo[["timestamp"]], decreasing = TRUE)[1], "file"]
          } else {
            pathname <- dbinfo[grep(db_version, dbinfo[["db_version"]], fixed = TRUE)[1], "file"]
            if (is.na(pathname)) {
              warning("There is no ", db_version, " version of the database. Use the latest version.", immediate. = TRUE)
              pathname <- dbinfo[order(dbinfo[["timestamp"]], decreasing = TRUE)[1], "file"]
            }
          }
          if (!is.na(pathname)) {
            header <- readCacheHeader(pathname)
            cached_version <- strsplit(header[["comment"]], "\\|")[[1]][1]
            message("Loading cached db: ", term, " version:", cached_version, " created:", header[["timestamp"]])
            db_loaded <- loadCache(pathname = pathname)
            Sys.sleep(0.5)
            db_list[[sps]][[term]] <- db_loaded
          }
        }
      }
    }

    db_species <- setNames(object = rep(sps, length(db)), nm = db)

    sp <- unlist(strsplit(sps, split = "_"))
    org_sp <- paste0("org.", paste0(substring(sp, 1, 1), collapse = ""), ".eg.db")
    org_key <- "ENTREZID"
    if (sps == "Arabidopsis_thaliana") {
      biomart <- "plants_mart"
      org_sp <- "org.At.tair.db"
      org_key <- "TAIR"
      default_IDtypes[c("GO", "GO_BP", "GO_CC", "GO_MF", "PFAM", "Chromosome", "GeneType", "Enzyme")] <- "tair_locus"
    }
    if (sps == "Saccharomyces_cerevisiae") {
      org_sp <- "org.Sc.sgd.db"
      # org_key <- "SGD"
      # default_IDtypes[c("GO", "GO_BP", "GO_CC", "GO_MF", "PFAM", "Chromosome", "GeneType", "Enzyme")] <- "sgd_gene"
    }

    ## Prepare ----------------------------------------------------------------------
    if (any(!sps %in% names(db_list)) || any(!db %in% names(db_list[[sps]]))) {
      orgdb_dependent <- c("GO", "GO_BP", "GO_CC", "GO_MF", "PFAM", "Chromosome", "GeneType", "Enzyme")
      if (any(orgdb_dependent %in% db)) {
        status <- tryCatch(
          {
            check_R(c(org_sp, "GO.db", "GOSemSim"))
          },
          error = identity
        )
        if (inherits(status, "error")) {
          warning("Annotation package ", org_sp, " does not exist.", immediate. = TRUE)
          if (isTRUE(convert_species)) {
            warning("Use the human annotation to create the ", paste0(intersect(db, orgdb_dependent), collapse = ","), " database for ", sps, immediate. = TRUE)
            org_sp <- "org.Hs.eg.db"
            db_species[intersect(db, orgdb_dependent)] <- "Homo_sapiens"
          } else {
            stop("Stop the preparation.")
          }
        }
        suppressPackageStartupMessages(require(org_sp, character.only = TRUE, quietly = TRUE))
        orgdb <- get(org_sp)
      }
      if ("PFAM" %in% db) {
        check_R("PFAM.db")
      }
      if ("Reactome" %in% db) {
        check_R("reactome.db")
      }
      # if ("MeSH" %in% db) {
      #   check_R(c("AHMeSHDbs", "MeSHDbi", "MeSH.db", "AnnotationHub"))
      # }

      if (is.null(custom_TERM2GENE)) {
        ## GO ---------------------------------------------------------------------------------
        if (any(db %in% c("GO", "GO_BP", "GO_CC", "GO_MF")) && any(!intersect(db, c("GO", "GO_BP", "GO_CC", "GO_MF")) %in% names(db_list[[sps]]))) {
          terms <- db[db %in% c("GO", "GO_BP", "GO_CC", "GO_MF")]
          bg <- suppressMessages(select(orgdb, keys = keys(orgdb), columns = c("GOALL", org_key)))
          bg <- unique(bg[!is.na(bg[["GOALL"]]), c("GOALL", "ONTOLOGYALL", org_key), drop = FALSE])
          bg2 <- suppressMessages(select(GO.db::GO.db, keys = keys(GO.db::GO.db), columns = c("GOID", "TERM")))
          bg <- merge(x = bg, by.x = "GOALL", y = bg2, by.y = "GOID", all.x = TRUE)
          for (subterm in terms) {
            message("Preparing database: ", subterm)
            if (subterm == "GO") {
              TERM2GENE <- bg[, c("GOALL", org_key)]
              TERM2NAME <- bg[, c("GOALL", "TERM")]
              colnames(TERM2GENE) <- c("Term", default_IDtypes[[subterm]])
              colnames(TERM2NAME) <- c("Term", "Name")
              TERM2NAME[["ONTOLOGY"]] <- bg[["ONTOLOGYALL"]]
              semData <- NULL
            } else {
              simpleterm <- unlist(strsplit(subterm, split = "_"))[2]
              TERM2GENE <- bg[which(bg[["ONTOLOGYALL"]] %in% simpleterm), c("GOALL", org_key)]
              TERM2NAME <- bg[which(bg[["ONTOLOGYALL"]] %in% simpleterm), c("GOALL", "TERM")]
              colnames(TERM2GENE) <- c("Term", default_IDtypes[[subterm]])
              colnames(TERM2NAME) <- c("Term", "Name")
              TERM2NAME[["ONTOLOGY"]] <- simpleterm
              semData <- suppressMessages(GOSemSim::godata(orgdb, ont = simpleterm))
            }
            TERM2GENE <- na.omit(unique(TERM2GENE))
            TERM2NAME <- na.omit(unique(TERM2NAME))
            version <- packageVersion(org_sp)
            db_list[[db_species[subterm]]][[subterm]][["TERM2GENE"]] <- TERM2GENE
            db_list[[db_species[subterm]]][[subterm]][["TERM2NAME"]] <- TERM2NAME
            db_list[[db_species[subterm]]][[subterm]][["semData"]] <- semData
            db_list[[db_species[subterm]]][[subterm]][["version"]] <- version
            if (sps == db_species[subterm]) {
              saveCache(db_list[[db_species[subterm]]][[subterm]],
                key = list(version, as.character(db_species[subterm]), subterm),
                comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species[subterm], "-", subterm)
              )
            }
          }
        }

        ## KEGG ---------------------------------------------------------------------------
        if (any(db == "KEGG") && (!"KEGG" %in% names(db_list[[sps]]))) {
          message("Preparing database: KEGG")
          check_R("httr")
          orgs <- kegg_get("https://rest.kegg.jp/list/organism")
          kegg_sp <- orgs[grep(gsub(pattern = "_", replacement = " ", x = sps), orgs[, 3]), 2]
          if (length(kegg_sp) == 0) {
            warning("Failed to prepare the KEGG database for ", db_species["KEGG"], immediate. = TRUE)
            if (isTRUE(convert_species) && db_species["KEGG"] != "Homo_sapiens") {
              warning("Use the human annotation to create the KEGG database for ", sps, immediate. = TRUE)
              db_species["KEGG"] <- "Homo_sapiens"
              kegg_sp <- "hsa"
              return(NULL)
            } else {
              stop("Stop the preparation.")
            }
          }
          kegg_db <- "pathway"

          kegg_pathwaygene_url <- paste0("https://rest.kegg.jp/link/", kegg_sp, "/", kegg_db, collapse = "")
          TERM2GENE <- kegg_get(kegg_pathwaygene_url)
          colnames(TERM2GENE) <- c("Pathway", "KEGG_ID")
          kegg_geneconversion_url <- paste0("https://rest.kegg.jp/conv/ncbi-geneid/", kegg_sp)
          GENECONV <- kegg_get(kegg_geneconversion_url)
          colnames(GENECONV) <- c("KEGG_ID", "ENTREZID")
          TERM2GENE <- merge(x = TERM2GENE, y = GENECONV, by = "KEGG_ID", all.x = TRUE)
          TERM2GENE[, "Pathway"] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, "Pathway"])
          TERM2GENE[, "ENTREZID"] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2GENE[, "ENTREZID"])
          TERM2GENE <- TERM2GENE[, c("Pathway", "ENTREZID")]

          kegg_pathwayname_url <- paste0("https://rest.kegg.jp/list/", kegg_db, "/", kegg_sp, collapse = "")
          TERM2NAME <- kegg_get(kegg_pathwayname_url)
          colnames(TERM2NAME) <- c("Pathway", "Name")
          TERM2NAME[, "Pathway"] <- gsub(pattern = "[^:]+:", replacement = "", x = TERM2NAME[, "Pathway"])
          TERM2NAME[, "Name"] <- gsub(pattern = paste0(" - ", paste0(unlist(strsplit(db_species["KEGG"], split = "_")), collapse = " "), ".*$"), replacement = "", x = TERM2NAME[, "Name"])
          TERM2NAME <- TERM2NAME[TERM2NAME[, "Pathway"] %in% TERM2GENE[, "Pathway"], , drop = FALSE]

          colnames(TERM2GENE) <- c("Term", default_IDtypes[["KEGG"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          # kegg_info <- readLines("https://rest.kegg.jp/info/hsa")
          kegg_info <- strsplit(httr::content(httr::GET(paste0("https://rest.kegg.jp/info/", kegg_sp))), split = "\n")[[1]]
          version <- gsub(".*(?=Release)", replacement = "", x = kegg_info[grepl("Release", x = kegg_info)], perl = TRUE)
          db_list[[db_species["KEGG"]]][["KEGG"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["KEGG"]]][["KEGG"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["KEGG"]]][["KEGG"]][["version"]] <- version
          if (sps == db_species["KEGG"]) {
            saveCache(db_list[[db_species["KEGG"]]][["KEGG"]],
              key = list(version, as.character(db_species["KEGG"]), "KEGG"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["KEGG"], "-KEGG")
            )
          }
        }

        ## WikiPathway ---------------------------------------------------------------------------
        if (any(db == "WikiPathway") && (!"WikiPathway" %in% names(db_list[[sps]]))) {
          message("Preparing database: WikiPathway")
          tempdir <- tempdir()
          gmt_files <- list.files(tempdir)[grep(".gmt", x = list.files(tempdir))]
          if (length(gmt_files) > 0) {
            file.remove(paste0(tempdir, "/", gmt_files))
          }
          temp <- tempfile()
          download(url = "https://wikipathways-data.wmcloud.org/current/gmt", destfile = temp)
          lines <- paste0(readLines(temp, warn = FALSE), collapse = " ")
          gmtfiles <- unlist(regmatches(lines, m = gregexpr("(?<=>)wikipathways-\\S+\\.gmt\\b", lines, perl = TRUE)))
          wiki_sp <- sps
          gmtfile <- gmtfiles[grep(wiki_sp, gmtfiles, fixed = TRUE)]
          if (length(gmtfile) == 0) {
            warning("Failed to prepare the WikiPathway database for ", db_species["WikiPathway"], immediate. = TRUE)
            if (isTRUE(convert_species) && db_species["WikiPathway"] != "Homo_sapiens") {
              warning("Use the human annotation to create the WikiPathway database for ", sps, immediate. = TRUE)
              db_species["WikiPathway"] <- "Homo_sapiens"
              wiki_sp <- "Homo_sapiens"
              gmtfile <- gmtfiles[grep(wiki_sp, gmtfiles, fixed = TRUE)]
            } else {
              stop("Stop the preparation.")
            }
          }
          version <- strsplit(gmtfile, split = "-")[[1]][[2]]
          download(url = paste0("https://wikipathways-data.wmcloud.org/current/gmt/", gmtfile), destfile = temp)
          wiki_gmt <- read.gmt(temp)
          unlink(temp)
          wiki_gmt <- apply(wiki_gmt, 1, function(x) {
            wikiid <- strsplit(x[["term"]], split = "%")[[1]][3]
            wikiterm <- strsplit(x[["term"]], split = "%")[[1]][1]
            gmt <- x[["gene"]]
            data.frame(v0 = wikiid, v1 = gmt, v2 = wikiterm, stringsAsFactors = FALSE)
          })
          bg <- do.call(rbind.data.frame, wiki_gmt)
          TERM2GENE <- bg[, c(1, 2)]
          TERM2NAME <- bg[, c(1, 3)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["WikiPathway"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["WikiPathway"]]][["WikiPathway"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["WikiPathway"]]][["WikiPathway"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["WikiPathway"]]][["WikiPathway"]][["version"]] <- version
          if (sps == db_species["WikiPathway"]) {
            saveCache(db_list[[db_species["WikiPathway"]]][["WikiPathway"]],
              key = list(version, as.character(db_species["WikiPathway"]), "WikiPathway"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["WikiPathway"], "-WikiPathway")
            )
          }
        }

        ## Pathwaycommons ---------------------------------------------------------------------------
        # check_R("paxtoolsr")

        ## Reactome ---------------------------------------------------------------------------
        if (any(db == "Reactome") && (!"Reactome" %in% names(db_list[[sps]]))) {
          message("Preparing database: Reactome")
          reactome_sp <- gsub(pattern = "_", replacement = " ", x = sps)
          df_all <- suppressMessages(select(reactome.db::reactome.db, keys = keys(reactome.db::reactome.db), columns = c("PATHID", "PATHNAME")))
          df <- df_all[grepl(pattern = paste0("^", reactome_sp, ": "), x = df_all$PATHNAME), , drop = FALSE]
          if (nrow(df) == 0) {
            if (isTRUE(convert_species) && db_species["Reactome"] != "Homo_sapiens") {
              warning("Use the human annotation to create the Reactome database for ", sps, immediate. = TRUE)
              db_species["Reactome"] <- "Homo_sapiens"
              reactome_sp <- gsub(pattern = "_", replacement = " ", x = "Homo_sapiens")
              df <- df_all[grepl(pattern = paste0("^", reactome_sp, ": "), x = df_all$PATHNAME), , drop = FALSE]
            } else {
              stop("Stop the preparation.")
            }
          }
          df <- na.omit(df)
          df$PATHNAME <- gsub(x = df$PATHNAME, pattern = paste0("^", reactome_sp, ": "), replacement = "", perl = TRUE)
          TERM2GENE <- df[, c(2, 1)]
          TERM2NAME <- df[, c(2, 3)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["Reactome"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- packageVersion("reactome.db")
          db_list[[db_species["Reactome"]]][["Reactome"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["Reactome"]]][["Reactome"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["Reactome"]]][["Reactome"]][["version"]] <- version
          if (sps == db_species["Reactome"]) {
            saveCache(db_list[[db_species["Reactome"]]][["Reactome"]],
              key = list(version, as.character(db_species["Reactome"]), "Reactome"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["Reactome"], "-Reactome")
            )
          }
        }

        ## CORUM ---------------------------------------------------------------------------
        if (any(db == "CORUM") && (!"CORUM" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the CORUM database for ", sps, immediate. = TRUE)
              db_species["CORUM"] <- "Homo_sapiens"
            } else {
              warning("CORUM database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: corum")
          url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/corum/gene_set_library_crisp.gmt.gz"
          temp <- tempfile(fileext = ".gz")
          download(url = url, destfile = temp)
          R.utils::gunzip(temp)
          TERM2GENE <- read.gmt(gsub(".gz", "", temp))
          version <- "Harmonizome 3.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["CORUM"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["CORUM"]]][["CORUM"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["CORUM"]]][["CORUM"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["CORUM"]]][["CORUM"]][["version"]] <- version
          if (sps == db_species["CORUM"]) {
            saveCache(db_list[[db_species["CORUM"]]][["CORUM"]],
              key = list(version, as.character(db_species["CORUM"]), "CORUM"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["CORUM"], "-CORUM")
            )
          }
        }

        ## DGI ---------------------------------------------------------------------------
        # if (any(db == "DGI") && (!"DGI" %in% names(db_list[[sps]]))) {
        #   if (sps != "Homo_sapiens") {
        #     if (isTRUE(convert_species)) {
        #       warning("Use the human annotation to create the DGI database for ", sps, immediate. = TRUE)
        #       db_species["DGI"] <- "Homo_sapiens"
        #     } else {
        #       warning("DGI database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
        #       stop("Stop the preparation.")
        #     }
        #   }
        #   message("Preparing database: DGI")
        #   temp <- tempfile()
        #   download(url = "https://www.dgidb.org/downloads", destfile = temp)
        #   lines <- readLines(temp, warn = FALSE)
        #   lines <- lines[grep("data/monthly_tsvs/.*interactions.tsv", lines, perl = TRUE)]
        #   lines <- gsub("(<td><a href=\\\")|(\\\">interactions.tsv</a></td>)", "", lines)
        #   version <- strsplit(lines[1], split = "/")[[1]][3]
        #   download(url = paste0("https://www.dgidb.org/", lines[1]), destfile = temp)
        #   dgi <- read.table(temp, header = TRUE, sep = "\t", fill = TRUE, quote = "")
        #   unlink(temp)
        #   dgi <- dgi[!is.na(dgi[["entrez_id"]]), , drop = FALSE]
        #   dgi <- dgi[, c("entrez_id", "drug_claim_name")]
        #   dgi[, "drug_claim_name"] <- toupper(dgi[, "drug_claim_name"])
        #   TERM2GENE <- dgi[, c("drug_claim_name", "entrez_id")]
        #   TERM2NAME <- dgi[, c("drug_claim_name", "drug_claim_name")]
        #   colnames(TERM2GENE) <- c("Term", default_IDtypes["DGI"])
        #   colnames(TERM2NAME) <- c("Term", "Name")
        #   TERM2GENE <- na.omit(unique(TERM2GENE))
        #   TERM2NAME <- na.omit(unique(TERM2NAME))
        #   db_list[[db_species["DGI"]]][["DGI"]][["TERM2GENE"]] <- TERM2GENE
        #   db_list[[db_species["DGI"]]][["DGI"]][["TERM2NAME"]] <- TERM2NAME
        #   db_list[[db_species["DGI"]]][["DGI"]][["version"]] <- version
        #   if (sps == db_species["DGI"]) {
        #     saveCache(db_list[[db_species["DGI"]]][["DGI"]],
        #       key = list(version, as.character(db_species["DGI"]), "DGI"),
        #       comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["DGI"], "-DGI")
        #     )
        #   }
        # }

        ## MP ---------------------------------------------------------------------------
        if (any(db == "MP") && (!"MP" %in% names(db_list[[sps]]))) {
          if (sps != "Mus_musculus") {
            if (isTRUE(convert_species)) {
              warning("Use the mouse annotation to create the MP database for ", sps, immediate. = TRUE)
              db_species["MP"] <- "Mus_musculus"
            } else {
              warning("MP database only support Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: MP")
          temp <- tempfile()
          download(url = "http://www.informatics.jax.org/downloads/reports/", destfile = temp)
          version <- readLines(temp, warn = FALSE)
          version <- version[grep("MGI_PhenoGenoMP.rpt", version)]
          version <- strsplit(version, split = "  </td><td align=\"right\">")[[1]][2]
          download(url = "http://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt", destfile = temp)
          mp_name <- read.table(temp, header = FALSE, sep = "\t", fill = TRUE, quote = "")
          rownames(mp_name) <- mp_name[, 1]
          download(url = "http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt", destfile = temp)
          gene_id <- read.table(temp, header = FALSE, row.names = NULL, sep = "\t", fill = TRUE, quote = "")
          gene_id <- gene_id[, 1:15]
          colnames(gene_id) <- gene_id[1, ]
          gene_id <- gene_id[gene_id[, 2] %in% c("Gene", "Pseudogene"), , drop = FALSE]
          rownames(gene_id) <- gene_id[, 1]

          # download(url = "http://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt", destfile = temp) # 43.0 MB
          # mp_gene <- read.table(temp, header = FALSE, sep = "\t", fill = TRUE, quote = "")
          # mp_gene[["symbol"]]<- gene_id[mp_gene[["V6"]], "3. marker symbol"]
          # mp_gene[["MP"]] <- mp_name[mp_gene[, "V4"], 2]
          # TERM2GENE <- mp_gene[, c("V4", "symbol")]
          # TERM2NAME <- mp_gene[, c("V4", "MP")]

          download(url = "http://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt", destfile = temp) # 32.4 MB
          mp_gene <- read.table(temp, header = FALSE, sep = "\t", fill = TRUE, quote = "")
          mp_gene[["symbol"]] <- gene_id[mp_gene[["V7"]], "3. marker symbol"]
          mp_gene[["MP"]] <- mp_name[mp_gene[, "V5"], 2]
          TERM2GENE <- mp_gene[, c("V5", "symbol")]
          TERM2NAME <- mp_gene[, c("V5", "MP")]

          colnames(TERM2GENE) <- c("Term", default_IDtypes[["MP"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["MP"]]][["MP"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["MP"]]][["MP"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["MP"]]][["MP"]][["version"]] <- version
          if (sps == db_species["MP"]) {
            saveCache(db_list[[db_species["MP"]]][["MP"]],
              key = list(version, as.character(db_species["MP"]), "MP"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["MP"], "-MP")
            )
          }
        }

        ## DO ---------------------------------------------------------------------------
        if (any(db == "DO") && (!"DO" %in% names(db_list[[sps]]))) {
          message("Preparing database: DO")
          temp <- tempfile(fileext = ".tsv.gz")
          download(url = "https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_COMBINED.tsv.gz", destfile = temp)
          R.utils::gunzip(temp)
          do_all <- read.table(gsub(".gz", "", temp), header = TRUE, sep = "\t", fill = TRUE, quote = "")
          version <- gsub(pattern = ".*Alliance Database Version: ", replacement = "", x = grep("Alliance Database Version", readLines(gsub(".gz", "", temp), warn = FALSE), perl = TRUE, value = TRUE))
          unlink(temp)
          do_sp <- gsub(pattern = "_", replacement = " ", x = sps)
          do_df <- do_all[do_all[["DBobjectType"]] == "gene" & do_all[["SpeciesName"]] == do_sp, , drop = FALSE]
          if (nrow(do_df) == 0) {
            if (isTRUE(convert_species) && db_species["DO"] != "Homo_sapiens") {
              warning("Use the human annotation to create the DO database for ", sps, immediate. = TRUE)
              db_species["DO"] <- "Homo_sapiens"
              do_sp <- gsub(pattern = "_", replacement = " ", x = "Homo_sapiens")
              do_df <- do_all[do_all[["DBobjectType"]] == "gene" & do_all[["SpeciesName"]] == do_sp, , drop = FALSE]
            } else {
              stop("Stop the preparation.")
            }
          }
          TERM2GENE <- do_df[, c("DOID", "DBObjectSymbol")]
          TERM2NAME <- do_df[, c("DOID", "DOtermName")]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["DO"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["DO"]]][["DO"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["DO"]]][["DO"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["DO"]]][["DO"]][["version"]] <- version
          if (sps == db_species["DO"]) {
            saveCache(db_list[[db_species["DO"]]][["DO"]],
              key = list(version, as.character(db_species["DO"]), "DO"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["DO"], "-DO")
            )
          }
        }

        ## HPO ---------------------------------------------------------------------------
        if (any(db == "HPO") && (!"HPO" %in% names(db_list[[sps]]))) {
          message("Preparing database: HPO")
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the HPO database for ", sps, immediate. = TRUE)
              db_species["HPO"] <- "Homo_sapiens"
            } else {
              warning("HPO database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          temp <- tempfile()
          download(url = "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases?per_page=1", destfile = temp)
          release <- readLines(temp, warn = FALSE)
          version <- regmatches(release, m = regexpr("(?<=tag_name\\\":\\\")\\S+(?=\\\",\\\"target_commitish)", release, perl = T))

          download(url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt", destfile = temp)
          hpo <- read.table(temp, header = TRUE, sep = "\t", fill = TRUE, quote = "")
          unlink(temp)

          TERM2GENE <- hpo[, c("hpo_id", "gene_symbol")]
          TERM2NAME <- hpo[, c("hpo_id", "hpo_name")]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["HPO"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["HPO"]]][["HPO"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["HPO"]]][["HPO"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["HPO"]]][["HPO"]][["version"]] <- version
          if (sps == db_species["HPO"]) {
            saveCache(db_list[[db_species["HPO"]]][["HPO"]],
              key = list(version, as.character(db_species["HPO"]), "HPO"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["HPO"], "-HPO")
            )
          }
        }

        ## PFAM ---------------------------------------------------------------------------
        if (any(db == "PFAM") && (!"PFAM" %in% names(db_list[[sps]]))) {
          message("Preparing database: PFAM")
          if (!"PFAM" %in% columns(orgdb)) {
            warning("PFAM is not in the orgdb: ", orgdb, " . Skip the preparation.", immediate. = TRUE)
          } else {
            bg <- suppressMessages(select(orgdb, keys = keys(orgdb), columns = c("PFAM", org_key)))
            bg <- unique(bg[!is.na(bg$PFAM), c("PFAM", org_key), drop = FALSE])
            bg2 <- as.data.frame(PFAM.db::PFAMDE2AC[mappedkeys(PFAM.db::PFAMDE2AC)])
            rownames(bg2) <- bg2[["ac"]]
            bg[["PFAM_name"]] <- bg2[bg$PFAM, "de"]
            bg[is.na(bg[["PFAM_name"]]), "PFAM_name"] <- bg[is.na(bg[["PFAM_name"]]), "PFAM"]
            TERM2GENE <- bg[, c("PFAM", org_key)]
            TERM2NAME <- bg[, c("PFAM", "PFAM_name")]
            colnames(TERM2GENE) <- c("Term", default_IDtypes[["PFAM"]])
            colnames(TERM2NAME) <- c("Term", "Name")
            TERM2GENE <- na.omit(unique(TERM2GENE))
            TERM2NAME <- na.omit(unique(TERM2NAME))
            version <- packageVersion(org_sp)
            db_list[[db_species["PFAM"]]][["PFAM"]][["TERM2GENE"]] <- TERM2GENE
            db_list[[db_species["PFAM"]]][["PFAM"]][["TERM2NAME"]] <- TERM2NAME
            db_list[[db_species["PFAM"]]][["PFAM"]][["version"]] <- version
            if (sps == db_species["PFAM"]) {
              saveCache(db_list[[db_species["PFAM"]]][["PFAM"]],
                key = list(version, as.character(db_species["PFAM"]), "PFAM"),
                comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["PFAM"], "-PFAM")
              )
            }
          }
        }

        ## Chromosome ---------------------------------------------------------------------------
        if (any(db == "Chromosome") && (!"Chromosome" %in% names(db_list[[sps]]))) {
          message("Preparing database: Chromosome")
          orgdbCHR <- get(paste0(gsub(pattern = ".db", "", org_sp), "CHR"))
          chr <- as.data.frame(orgdbCHR[mappedkeys(orgdbCHR)])
          chr[, 2] <- paste0("chr", chr[, 2])
          TERM2GENE <- chr[, c(2, 1)]
          TERM2NAME <- chr[, c(2, 2)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["Chromosome"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- packageVersion(org_sp)
          db_list[[db_species["Chromosome"]]][["Chromosome"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["Chromosome"]]][["Chromosome"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["Chromosome"]]][["Chromosome"]][["version"]] <- version
          if (sps == db_species["Chromosome"]) {
            saveCache(db_list[[db_species["Chromosome"]]][["Chromosome"]],
              key = list(version, as.character(db_species["Chromosome"]), "Chromosome"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["Chromosome"], "-Chromosome")
            )
          }
        }

        ## GeneType ---------------------------------------------------------------------------
        if (any(db == "GeneType") && (!"GeneType" %in% names(db_list[[sps]]))) {
          message("Preparing database: GeneType")
          if (!"GENETYPE" %in% columns(orgdb)) {
            warning("GENETYPE is not in the orgdb: ", org_sp, " . Skip the preparation.", immediate. = TRUE)
          } else {
            bg <- suppressMessages(select(orgdb, keys = keys(orgdb), columns = c("GENETYPE", org_key)))
            TERM2GENE <- bg[, c("GENETYPE", org_key)]
            TERM2NAME <- bg[, c("GENETYPE", "GENETYPE")]
            colnames(TERM2GENE) <- c("Term", default_IDtypes[["GeneType"]])
            colnames(TERM2NAME) <- c("Term", "Name")
            TERM2GENE <- na.omit(unique(TERM2GENE))
            TERM2NAME <- na.omit(unique(TERM2NAME))
            version <- packageVersion(org_sp)
            db_list[[db_species["GeneType"]]][["GeneType"]][["TERM2GENE"]] <- TERM2GENE
            db_list[[db_species["GeneType"]]][["GeneType"]][["TERM2NAME"]] <- TERM2NAME
            db_list[[db_species["GeneType"]]][["GeneType"]][["version"]] <- version
            if (sps == db_species["GeneType"]) {
              saveCache(db_list[[db_species["GeneType"]]][["GeneType"]],
                key = list(version, as.character(db_species["GeneType"]), "GeneType"),
                comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["GeneType"], "-GeneType")
              )
            }
          }
        }

        ## Enzyme ---------------------------------------------------------------------------
        if (any(db == "Enzyme") && (!"Enzyme" %in% names(db_list[[sps]]))) {
          message("Preparing database: Enzyme")
          if (!"ENZYME" %in% columns(orgdb)) {
            warning("ENZYME is not in the orgdb: ", orgdb, " . Skip the preparation.", immediate. = TRUE)
          } else {
            bg <- suppressMessages(select(orgdb, keys = keys(orgdb), columns = c("ENZYME", org_key)))
            bg1 <- bg2 <- na.omit(bg)
            bg1[, "ENZYME"] <- sapply(strsplit(bg1[, "ENZYME"], "\\."), function(x) paste0(head(x, 1), collapse = "."))
            bg2[, "ENZYME"] <- sapply(strsplit(bg2[, "ENZYME"], "\\."), function(x) paste0(head(x, 2), collapse = "."))
            bg <- unique(rbind(bg1, bg2))
            bg[, "ENZYME"] <- gsub(pattern = "\\.-$", "", x = bg[, 2])
            bg[, "ENZYME"] <- paste0("ec:", bg[, "ENZYME"])
            temp <- tempfile()
            download(url = "https://ftp.expasy.org/databases/enzyme/enzclass.txt", destfile = temp)
            enzyme <- read.table(temp, header = FALSE, sep = "\t", fill = TRUE, quote = "")
            enzyme <- enzyme[grep("-.-", enzyme[, 1], fixed = TRUE), , drop = FALSE]
            enzyme <- do.call(rbind, strsplit(enzyme[, 1], split = ". -.-  "))
            enzyme[, 1] <- paste0("ec:", gsub(pattern = "( )|(. -)", replacement = "", enzyme[, 1]))
            enzyme[, 2] <- gsub(pattern = "(^ )|(\\.$)", replacement = "", enzyme[, 2])
            rownames(enzyme) <- enzyme[, 1]
            for (i in seq_len(nrow(enzyme))) {
              if (grepl(".", enzyme[i, 1], fixed = TRUE)) {
                enzyme[i, 2] <- paste0(enzyme[strsplit(enzyme[i, 1], ".", fixed = TRUE)[[1]][1], 2], "(", enzyme[i, 2], ")")
              }
            }
            unlink(temp)
            bg[, "Name"] <- enzyme[bg[, "ENZYME"], 2]
            TERM2GENE <- bg[, c("ENZYME", org_key)]
            TERM2NAME <- bg[, c("ENZYME", "Name")]
            colnames(TERM2GENE) <- c("Term", default_IDtypes[["Enzyme"]])
            colnames(TERM2NAME) <- c("Term", "Name")
            TERM2GENE <- na.omit(unique(TERM2GENE))
            TERM2NAME <- na.omit(unique(TERM2NAME))
            version <- packageVersion(org_sp)
            db_list[[db_species["Enzyme"]]][["Enzyme"]][["TERM2GENE"]] <- TERM2GENE
            db_list[[db_species["Enzyme"]]][["Enzyme"]][["TERM2NAME"]] <- TERM2NAME
            db_list[[db_species["Enzyme"]]][["Enzyme"]][["version"]] <- version
            if (sps == db_species["Enzyme"]) {
              saveCache(db_list[[db_species["Enzyme"]]][["Enzyme"]],
                key = list(version, as.character(db_species["Enzyme"]), "Enzyme"),
                comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["Enzyme"], "-Enzyme")
              )
            }
          }
        }

        ## TF ---------------------------------------------------------------------------
        if (any(db == "TF") && (!"TF" %in% names(db_list[[sps]]))) {
          message("Preparing database: TF")

          # AnimalTFDB4
          status <- tryCatch(
            {
              temp <- tempfile()
              url <- paste0("http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/", sps, "_TF")
              download(url = url, destfile = temp)
              tf <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
              url <- paste0("http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/Cof_list_final/", sps, "_Cof")
              download(url = url, destfile = temp)
              tfco <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
              if (!"Symbol" %in% colnames(tf)) {
                if (isTRUE(convert_species) && db_species["TF"] != "Homo_sapiens") {
                  warning("Use the human annotation to create the TF database for ", sps, immediate. = TRUE)
                  db_species["TF"] <- "Homo_sapiens"
                  url <- paste0("http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/Homo_sapiens_TF")
                  download(url = url, destfile = temp)
                  tf <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
                  url <- paste0("http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/Cof_list_final/Homo_sapiens_Cof")
                  download(url = url, destfile = temp)
                  tfco <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
                } else {
                  stop("Stop the preparation.")
                }
              }
              unlink(temp)
              version <- "AnimalTFDB4"
            },
            error = identity
          )

          # AnimalTFDB3
          if (inherits(status, "error")) {
            temp <- tempfile()
            url <- paste0("https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/", sps, "_TF")
            download(url = url, destfile = temp)
            tf <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
            url <- paste0("https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/", sps, "_TF_cofactors")
            download(url = url, destfile = temp)
            tfco <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
            if (!"Symbol" %in% colnames(tf)) {
              if (isTRUE(convert_species) && db_species["TF"] != "Homo_sapiens") {
                warning("Use the human annotation to create the TF database for ", sps, immediate. = TRUE)
                db_species["TF"] <- "Homo_sapiens"
                url <- paste0("https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/Homo_sapiens_TF")
                download(url = url, destfile = temp)
                tf <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
                url <- paste0("https://raw.githubusercontent.com/GuoBioinfoLab/AnimalTFDB3/master/AnimalTFDB3/static/AnimalTFDB3/download/Homo_sapiens_TF_cofactors")
                download(url = url, destfile = temp)
                tfco <- read.table(temp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
              } else {
                stop("Stop the preparation.")
              }
            }
            unlink(temp)
            version <- "AnimalTFDB3"
          }

          TERM2GENE <- rbind(data.frame("Term" = "TF", "symbol" = tf[["Symbol"]]), data.frame("Term" = "TF cofactor", "symbol" = tfco[["Symbol"]]))
          TERM2NAME <- data.frame("Term" = c("TF", "TF cofactor"), "Name" = c("TF", "TF cofactor"))
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["TF"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["TF"]]][["TF"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["TF"]]][["TF"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["TF"]]][["TF"]][["version"]] <- version
          if (sps == db_species["TF"]) {
            saveCache(db_list[[db_species["TF"]]][["TF"]],
              key = list(version, as.character(db_species["TF"]), "TF"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["TF"], "-TF")
            )
          }
        }

        ## CSPA ---------------------------------------------------------------------------
        if (any(db == "CSPA") && (!"CSPA" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the CSPA database for ", sps, immediate. = TRUE)
              db_species["CSPA"] <- "Homo_sapiens"
            } else {
              warning("CSPA database only support Homo_sapiens and Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          check_R("openxlsx")
          message("Preparing database: CSPA")
          temp <- tempfile(fileext = ".xlsx")
          url <- "https://wlab.ethz.ch/cspa/data/S1_File.xlsx"
          download(url = url, destfile = temp, mode = ifelse(.Platform$OS.type == "windows", "wb", "w"))
          surfacepro <- openxlsx::read.xlsx(temp, sheet = 1)
          unlink(temp)
          surfacepro <- surfacepro[surfacepro[["organism"]] == switch(db_species["CSPA"],
            "Homo_sapiens" = "Human",
            "Mus_musculus" = "Mouse"
          ), , drop = FALSE]
          TERM2GENE <- data.frame("Term" = "SurfaceProtein", "symbol" = surfacepro[["ENTREZ.gene.symbol"]])
          TERM2NAME <- data.frame("Term" = "SurfaceProtein", "Name" = "SurfaceProtein")
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["CSPA"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- "CSPA"
          db_list[[db_species["CSPA"]]][["CSPA"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["CSPA"]]][["CSPA"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["CSPA"]]][["CSPA"]][["version"]] <- version
          if (sps == db_species["CSPA"]) {
            saveCache(db_list[[db_species["CSPA"]]][["CSPA"]],
              key = list(version, as.character(db_species["CSPA"]), "CSPA"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["CSPA"], "-CSPA")
            )
          }
        }

        ## Surfaceome ---------------------------------------------------------------------------
        if (any(db == "Surfaceome") && (!"Surfaceome" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the Surfaceome database for ", sps, immediate. = TRUE)
              db_species["Surfaceome"] <- "Homo_sapiens"
            } else {
              warning("Surfaceome database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          check_R("openxlsx")
          message("Preparing database: Surfaceome")
          temp <- tempfile(fileext = ".xlsx")
          url <- "http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"
          download(url = url, destfile = temp, mode = ifelse(.Platform$OS.type == "windows", "wb", "w"))
          surfaceome <- openxlsx::read.xlsx(temp, sheet = 2, colNames = TRUE, startRow = 2)
          unlink(temp)
          TERM2GENE <- data.frame("Term" = "SurfaceProtein", "symbol" = surfaceome[["UniProt.gene"]])
          TERM2NAME <- data.frame("Term" = "SurfaceProtein", "Name" = "SurfaceProtein")
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["Surfaceome"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- "Surfaceome"
          db_list[[db_species["Surfaceome"]]][["Surfaceome"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["Surfaceome"]]][["Surfaceome"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["Surfaceome"]]][["Surfaceome"]][["version"]] <- version
          if (sps == db_species["Surfaceome"]) {
            saveCache(db_list[[db_species["Surfaceome"]]][["Surfaceome"]],
              key = list(version, as.character(db_species["Surfaceome"]), "Surfaceome"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["Surfaceome"], "-Surfaceome")
            )
          }
        }

        ## SPRomeDB ---------------------------------------------------------------------------
        if (any(db == "SPRomeDB") && (!"SPRomeDB" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the SPRomeDB database for ", sps, immediate. = TRUE)
              db_species["SPRomeDB"] <- "Homo_sapiens"
            } else {
              warning("SPRomeDB database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: SPRomeDB")
          temp <- tempfile()
          url <- "http://119.3.41.228/SPRomeDB/files/download/secreted_proteins_SPRomeDB.csv"
          download(url = url, destfile = temp)
          spromedb <- read.csv(temp, header = TRUE)
          unlink(temp)
          TERM2GENE <- data.frame("Term" = "SecretoryProtein", "entrez_id" = unlist(strsplit(spromedb$Gene_ID, ";")))
          TERM2NAME <- data.frame("Term" = "SecretoryProtein", "Name" = "SecretoryProtein")
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["SPRomeDB"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- "SPRomeDB"
          db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]][["version"]] <- version
          if (sps == db_species["SPRomeDB"]) {
            saveCache(db_list[[db_species["SPRomeDB"]]][["SPRomeDB"]],
              key = list(version, as.character(db_species["SPRomeDB"]), "SPRomeDB"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["SPRomeDB"], "-SPRomeDB")
            )
          }
        }

        ## VerSeDa ---------------------------------------------------------------------------
        if (any(db == "VerSeDa") && (!"VerSeDa" %in% names(db_list[[sps]]))) {
          temp <- tempfile()
          download(url = "http://genomics.cicbiogune.es/VerSeDa/downloads.php", destfile = temp)
          verseda_sps <- readLines(temp)
          verseda_sps <- regmatches(verseda_sps, m = regexpr("(?<=Downloads/)\\S+(?=\\.zip)", verseda_sps, perl = TRUE))
          verseda_sps <- setdiff(verseda_sps, c("NonRefined", "NonRefined_Curated", "Refined", "Refined_Curated"))
          if (!tolower(sps) %in% verseda_sps) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the VerSeDa database for ", sps, immediate. = TRUE)
              db_species["VerSeDa"] <- "Homo_sapiens"
            } else {
              warning("VerSeDa database only support ", paste0(verseda_sps, collapse = ","), ". Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: VerSeDa")
          temp <- tempfile(fileext = ".zip")
          url <- paste0("http://genomics.cicbiogune.es/VerSeDa/Downloads/", tolower(db_species["VerSeDa"]), ".zip")
          download(url = url, destfile = temp)
          con <- unz(temp, paste0(tolower(db_species["VerSeDa"]), "/", tolower(db_species["VerSeDa"]), "_Refined.sequences"))
          verseda <- readLines(con)
          close(con)
          unlink(temp)
          verseda <- verseda[grep("^>", verseda)]
          verseda <- gsub("^>|\\.\\d+", "", verseda)
          verseda_id <- GeneConvert(geneID = verseda, geneID_from_IDtype = c("ensembl_peptide_id", "refseq_peptide", "refseq_peptide_predicted", "uniprot_isoform", "uniprotswissprot", "uniprotsptrembl"), geneID_to_IDtype = "symbol", species_from = db_species["VerSeDa"])
          TERM2GENE <- data.frame("Term" = "SecretoryProtein", "symbol" = unique(verseda_id$geneID_expand$symbol))
          TERM2NAME <- data.frame("Term" = "SecretoryProtein", "Name" = "SecretoryProtein")
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["VerSeDa"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          version <- "VerSeDa"
          db_list[[db_species["VerSeDa"]]][["VerSeDa"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["VerSeDa"]]][["VerSeDa"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["VerSeDa"]]][["VerSeDa"]][["version"]] <- version
          if (sps == db_species["VerSeDa"]) {
            saveCache(db_list[[db_species["VerSeDa"]]][["VerSeDa"]],
              key = list(version, as.character(db_species["VerSeDa"]), "VerSeDa"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["VerSeDa"], "-VerSeDa")
            )
          }
        }

        ## TFLink ---------------------------------------------------------------------------
        if (any(db == "TFLink") && (!"TFLink" %in% names(db_list[[sps]]))) {
          tflink_sp <- c("Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Danio_rerio", "Drosophila_melanogaster", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae")
          if (!sps %in% tflink_sp) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the TFLink database for ", sps, immediate. = TRUE)
              db_species["TFLink"] <- "Homo_sapiens"
            } else {
              warning("TFLink database only support ", paste0(tflink_sp, collapse = ","), ". Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: TFLink")
          url <- paste0("https://cdn.netbiol.org/tflink/download_files/TFLink_", db_species["TFLink"], "_interactions_All_GMT_proteinName_v1.0.gmt")
          temp <- tempfile()
          download(url = url, destfile = temp)
          TERM2GENE <- read.gmt(temp)
          version <- "v1.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["TFLink"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["TFLink"]]][["TFLink"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["TFLink"]]][["TFLink"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["TFLink"]]][["TFLink"]][["version"]] <- version
          if (sps == db_species["TFLink"]) {
            saveCache(db_list[[db_species["TFLink"]]][["TFLink"]],
              key = list(version, as.character(db_species["TFLink"]), "TFLink"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["TFLink"], "-TFLink")
            )
          }
        }

        ## hTFtarget ---------------------------------------------------------------------------
        if (any(db == "hTFtarget") && (!"hTFtarget" %in% names(db_list[[sps]]))) {
          if (!sps %in% "Homo_sapiens") {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the hTFtarget database for ", sps, immediate. = TRUE)
              db_species["hTFtarget"] <- "Homo_sapiens"
            } else {
              warning("hTFtarget database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: hTFtarget")
          url <- paste0("http://bioinfo.life.hust.edu.cn/static/hTFtarget/file_download/tf-target-infomation.txt")
          temp <- tempfile()
          download(url = url, destfile = temp)
          TERM2GENE <- read.table(temp, header = TRUE, fill = T, sep = "\t")
          version <- "v1.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["hTFtarget"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["hTFtarget"]]][["hTFtarget"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["hTFtarget"]]][["hTFtarget"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["hTFtarget"]]][["hTFtarget"]][["version"]] <- version
          if (sps == db_species["hTFtarget"]) {
            saveCache(db_list[[db_species["hTFtarget"]]][["hTFtarget"]],
              key = list(version, as.character(db_species["hTFtarget"]), "hTFtarget"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["hTFtarget"], "-hTFtarget")
            )
          }
        }

        ## TRRUST ---------------------------------------------------------------------------
        if (any(db == "TRRUST") && (!"TRRUST" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the TRRUST database for ", sps, immediate. = TRUE)
              db_species["TRRUST"] <- "Homo_sapiens"
            } else {
              warning("hTFtarget database only support Homo_sapiens and Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: TRRUST")
          url <- switch(db_species["TRRUST"],
            "Homo_sapiens" = "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.human.tsv",
            "Mus_musculus" = "https://raw.githubusercontent.com/bioinfonerd/Transcription-Factor-Databases/master/Ttrust_v2/trrust_rawdata.mouse.tsv.gz"
          )
          if (endsWith(url, "gz")) {
            temp <- tempfile(fileext = ".gz")
            download(url = url, destfile = temp)
            R.utils::gunzip(temp)
            TERM2GENE <- read.table(gsub(".gz$", "", temp), header = FALSE, fill = T, sep = "\t")[, 1:2]
          } else {
            temp <- tempfile()
            download(url = url, destfile = temp)
            TERM2GENE <- read.table(temp, header = FALSE, fill = T, sep = "\t")[, 1:2]
          }
          version <- "v2.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["TRRUST"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["TRRUST"]]][["TRRUST"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["TRRUST"]]][["TRRUST"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["TRRUST"]]][["TRRUST"]][["version"]] <- version
          if (sps == db_species["TRRUST"]) {
            saveCache(db_list[[db_species["TRRUST"]]][["TRRUST"]],
              key = list(version, as.character(db_species["TRRUST"]), "TRRUST"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["TRRUST"], "-TRRUST")
            )
          }
        }

        ## JASPAR ---------------------------------------------------------------------------
        if (any(db == "JASPAR") && (!"JASPAR" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the JASPAR database for ", sps, immediate. = TRUE)
              db_species["JASPAR"] <- "Homo_sapiens"
            } else {
              warning("JASPAR database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: JASPAR")
          url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/jasparpwm/gene_set_library_crisp.gmt.gz"
          temp <- tempfile(fileext = ".gz")
          download(url = url, destfile = temp)
          R.utils::gunzip(temp)
          TERM2GENE <- read.gmt(gsub(".gz", "", temp))
          version <- "Harmonizome 3.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["JASPAR"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["JASPAR"]]][["JASPAR"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["JASPAR"]]][["JASPAR"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["JASPAR"]]][["JASPAR"]][["version"]] <- version
          if (sps == db_species["JASPAR"]) {
            saveCache(db_list[[db_species["JASPAR"]]][["JASPAR"]],
              key = list(version, as.character(db_species["JASPAR"]), "JASPAR"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["JASPAR"], "-JASPAR")
            )
          }
        }

        ## ENCODE ---------------------------------------------------------------------------
        if (any(db == "ENCODE") && (!"ENCODE" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the ENCODE database for ", sps, immediate. = TRUE)
              db_species["ENCODE"] <- "Homo_sapiens"
            } else {
              warning("ENCODE database only support Homo_sapiens. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: ENCODE")
          url <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/encodetfppi/gene_set_library_crisp.gmt.gz"
          temp <- tempfile(fileext = ".gz")
          download(url = url, destfile = temp)
          R.utils::gunzip(temp)
          TERM2GENE <- read.gmt(gsub(".gz", "", temp))
          version <- "Harmonizome 3.0"
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["ENCODE"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["ENCODE"]]][["ENCODE"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["ENCODE"]]][["ENCODE"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["ENCODE"]]][["ENCODE"]][["version"]] <- version
          if (sps == db_species["ENCODE"]) {
            saveCache(db_list[[db_species["ENCODE"]]][["ENCODE"]],
              key = list(version, as.character(db_species["ENCODE"]), "ENCODE"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["ENCODE"], "-ENCODE")
            )
          }
        }

        ## MSigDB ---------------------------------------------------------------------------
        if (any(grepl("MSigDB", db)) && (!"MSigDB" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the MSigDB database for ", sps, immediate. = TRUE)
              db_species["MSigDB"] <- "Homo_sapiens"
            } else {
              warning("MSigDB database only support Homo_sapiens and Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: MSigDB")

          temp <- tempfile()
          download(url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/", destfile = temp)
          version <- readLines(temp)
          version <- version[grep("alt=\"\\[DIR\\]\"", version)]
          version <- version[grep(switch(db_species["MSigDB"],
            "Homo_sapiens" = "Hs",
            "Mus_musculus" = "Mm"
          ), version)]
          version <- version[length(version)]
          version <- regmatches(version, m = regexpr("(?<=href\\=\")\\S+(?=/\"\\>)", version, perl = TRUE))
          # version="2023.1.Hs" "2023.2.Hs"

          url <- paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/", version, "/msigdb.v", version, ".json")
          download(url = url, destfile = temp)
          lines <- paste0(readLines(temp, warn = FALSE), collapse = "")
          lines <- gsub("\"", "", lines)
          lines <- gsub("^\\{|\\}$", "", lines)
          terms <- strsplit(lines, "\\},")[[1]]
          term_list <- list()
          for (idx in seq_along(terms)) {
            term_content <- terms[[idx]]
            term_name <- trimws(gsub(":|_|\\{.*", " ", term_content))
            term_id <- trimws(regmatches(term_content, m = regexpr("(?<=systematicName:)\\S+?(?=,)", term_content, perl = TRUE)))
            term_gene <- trimws(regmatches(term_content, m = regexpr(pattern = "(?<=geneSymbols:\\[)\\S+?(?=\\],)", term_content, perl = TRUE)))
            term_collection <- trimws(regmatches(term_content, m = regexpr(pattern = "(?<=collection:)\\S+?(?=\\,)", term_content, perl = TRUE)))
            term_collection <- gsub(":.*", "", term_collection)
            term_list[[idx]] <- c(id = term_id, "name" = term_name, "gene" = term_gene, "collection" = term_collection)
          }
          df <- as.data.frame(do.call(rbind, term_list))
          df$gene <- strsplit(df$gene, split = ",")
          df <- unnest(df, cols = "gene")

          TERM2NAME <- df[, c(1, 2, 4)]
          TERM2GENE <- df[, c(1, 3)]
          colnames(TERM2NAME) <- c("Term", "Name", "Collection")
          colnames(TERM2GENE) <- c("Term", paste0(default_IDtypes[["MSigDB"]], collapse = "."))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          TERM2GENE <- na.omit(unique(TERM2GENE))

          db_list[[db_species["MSigDB"]]][["MSigDB"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["MSigDB"]]][["MSigDB"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["MSigDB"]]][["MSigDB"]][["version"]] <- version
          if (sps == db_species["MSigDB"]) {
            saveCache(db_list[[db_species["MSigDB"]]][["MSigDB"]],
              key = list(version, as.character(db_species["MSigDB"]), "MSigDB"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["MSigDB"], "-MSigDB")
            )
          }

          for (collection in unique(TERM2NAME[["Collection"]])) {
            db_species[paste0("MSigDB_", collection)] <- db_species["MSigDB"]
            default_IDtypes[[paste0("MSigDB_", collection)]] <- default_IDtypes[["MSigDB"]]
            TERM2NAME_sub <- TERM2NAME[TERM2NAME[["Collection"]] == collection, , drop = FALSE]
            TERM2GENE_sub <- TERM2GENE[TERM2GENE[["Term"]] %in% TERM2NAME_sub[["Term"]], , drop = FALSE]
            db_list[[db_species["MSigDB"]]][[paste0("MSigDB_", collection)]][["TERM2GENE"]] <- TERM2GENE_sub
            db_list[[db_species["MSigDB"]]][[paste0("MSigDB_", collection)]][["TERM2NAME"]] <- TERM2NAME_sub
            db_list[[db_species["MSigDB"]]][[paste0("MSigDB_", collection)]][["version"]] <- version
            if (sps == db_species["MSigDB"]) {
              saveCache(db_list[[db_species["MSigDB"]]][[paste0("MSigDB_", collection)]],
                key = list(version, as.character(db_species["MSigDB"]), paste0("MSigDB_", collection)),
                comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["MSigDB"], "-MSigDB_", collection)
              )
            }
          }
        }

        ## CellTalk ---------------------------------------------------------------------------
        if (any(db == "CellTalk") && (!"CellTalk" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the CellTalk database for ", sps, immediate. = TRUE)
              db_species["CellTalk"] <- "Homo_sapiens"
            } else {
              warning("CellTalk database only support Homo_sapiens and Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: CellTalk")
          url <- switch(db_species["CellTalk"],
            "Homo_sapiens" = "https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/human_lr_pair.rds",
            "Mus_musculus" = "https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/mouse_lr_pair.rds"
          )

          temp <- tempfile()
          download(url = url, destfile = temp)
          lr <- readRDS(temp)
          version <- "v1.0"

          lr[["ligand_gene_symbol2"]] <- paste0("ligand_", lr[["ligand_gene_symbol"]])
          lr[["receptor_gene_symbol2"]] <- paste0("receptor_", lr[["receptor_gene_symbol"]])
          TERM2GENE <- rbind(
            data.frame("Term" = lr[["ligand_gene_symbol2"]], "symbol" = lr[["receptor_gene_symbol"]]),
            data.frame("Term" = lr[["receptor_gene_symbol2"]], "symbol" = lr[["ligand_gene_symbol"]])
          )
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["CellTalk"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["CellTalk"]]][["CellTalk"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["CellTalk"]]][["CellTalk"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["CellTalk"]]][["CellTalk"]][["version"]] <- version
          if (sps == db_species["CellTalk"]) {
            saveCache(db_list[[db_species["CellTalk"]]][["CellTalk"]],
              key = list(version, as.character(db_species["CellTalk"]), "CellTalk"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["CellTalk"], "-CellTalk")
            )
          }
        }

        ## CellChat ---------------------------------------------------------------------------
        if (any(db == "CellChat") && (!"CellChat" %in% names(db_list[[sps]]))) {
          if (!sps %in% c("Homo_sapiens", "Mus_musculus")) {
            if (isTRUE(convert_species)) {
              warning("Use the human annotation to create the CellChat database for ", sps, immediate. = TRUE)
              db_species["CellChat"] <- "Homo_sapiens"
            } else {
              warning("CellChat database only support Homo_sapiens and Mus_musculus. Consider using convert_species=TRUE", immediate. = TRUE)
              stop("Stop the preparation.")
            }
          }
          message("Preparing database: CellChat")
          url <- paste0(
            "https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.",
            switch(db_species["CellChat"],
              "Homo_sapiens" = "human.rda",
              "Mus_musculus" = "mouse.rda",
              "Danio_rerio" = "zebrafish.rda"
            )
          )
          temp <- tempfile()
          download(url = url, destfile = temp)
          load(temp)
          lr <- get(paste0("CellChatDB.", switch(db_species["CellChat"],
            "Homo_sapiens" = "human",
            "Mus_musculus" = "mouse",
            "Danio_rerio" = "zebrafish"
          )))[["interaction"]]
          download(url = "https://raw.githubusercontent.com/sqjin/CellChat/master/DESCRIPTION", destfile = temp)
          version <- grep(pattern = "Version", x = readLines(temp), value = TRUE)
          version <- gsub(pattern = "(.*Version: )|(</td>)", replacement = "", x = version)
          unlink(temp)

          lr_list <- strsplit(lr$interaction_name, split = "_")
          lr[["ligand_gene_symbol"]] <- paste0("ligand_", sapply(lr_list, function(x) x[[1]]))
          lr[["receptor_list"]] <- lapply(lr_list, function(x) paste0("receptor_", x[2:length(x)]))
          lr <- unnest(data = lr, cols = "receptor_list", keep_empty = FALSE)
          TERM2GENE <- rbind(
            data.frame("Term" = lr[["ligand_gene_symbol"]], "symbol" = gsub(pattern = "receptor_", replacement = "", lr[["receptor_list"]])),
            data.frame("Term" = lr[["receptor_list"]], "symbol" = gsub(pattern = "ligand_", replacement = "", lr[["ligand_gene_symbol"]]))
          )

          if (db_species["CellChat"] == "Homo_sapiens") {
            TERM2GENE[["symbol"]] <- toupper(TERM2GENE[["symbol"]])
          } else if (db_species["CellChat"] == "Mus_musculus") {
            TERM2GENE[["symbol"]] <- capitalize(TERM2GENE[["symbol"]], force_tolower = TRUE)
          } else if (db_species["CellChat"] == "Danio_rerio") {
            TERM2GENE[["symbol"]] <- tolower(TERM2GENE[["symbol"]])
          }
          TERM2NAME <- TERM2GENE[, c(1, 1)]
          colnames(TERM2GENE) <- c("Term", default_IDtypes[["CellChat"]])
          colnames(TERM2NAME) <- c("Term", "Name")
          TERM2GENE <- na.omit(unique(TERM2GENE))
          TERM2NAME <- na.omit(unique(TERM2NAME))
          db_list[[db_species["CellChat"]]][["CellChat"]][["TERM2GENE"]] <- TERM2GENE
          db_list[[db_species["CellChat"]]][["CellChat"]][["TERM2NAME"]] <- TERM2NAME
          db_list[[db_species["CellChat"]]][["CellChat"]][["version"]] <- version
          if (sps == db_species["CellChat"]) {
            saveCache(db_list[[db_species["CellChat"]]][["CellChat"]],
              key = list(version, as.character(db_species["CellChat"]), "CellChat"),
              comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", db_species["CellChat"], "-CellChat")
            )
          }
        }
      } else {
        ## Custom ---------------------------------------------------------------------------
        db_species[db] <- custom_species
        if (sps != custom_species) {
          if (isTRUE(convert_species)) {
            warning("Use the ", custom_species, " annotation to create the ", db, " database for ", sps, immediate. = TRUE)
          } else {
            warning(db, " database only support ", custom_species, ". Consider using convert_species=TRUE", immediate. = TRUE)
            stop("Stop the preparation.")
          }
        }
        TERM2GENE <- custom_TERM2GENE
        colnames(TERM2GENE) <- c("Term", custom_IDtype)
        if (is.null(custom_TERM2NAME)) {
          TERM2NAME <- TERM2GENE[, c(1, 1)]
        } else {
          TERM2NAME <- custom_TERM2NAME
        }
        colnames(TERM2NAME) <- c("Term", "Name")

        TERM2GENE <- na.omit(unique(TERM2GENE))
        TERM2NAME <- na.omit(unique(TERM2NAME))
        db_list[[db_species[db]]][[db]][["TERM2GENE"]] <- TERM2GENE
        db_list[[db_species[db]]][[db]][["TERM2NAME"]] <- TERM2NAME
        db_list[[db_species[db]]][[db]][["version"]] <- custom_version
        if (sps == db_species[db]) {
          saveCache(db_list[[db_species[db]]][[db]],
            key = list(custom_version, as.character(db_species[db]), db),
            comment = paste0(custom_version, " nterm:", length(TERM2NAME[[1]]), "|", db_species[db], "-", db)
          )
        }
      }

      ## MeSH ---------------------------------------------------------------------------
      # if (any(db == "MeSH") && (!"MeSH" %in% names(db_list[[sps]]))) {
      #   message("Preparing database: MeSH")
      #   # dir.create("~/.cache/R/AnnotationHub",recursive = TRUE,showWarnings = FALSE)
      #   ### A (Anatomy);B (Organisms);C (Diseases);D (Chemicals and Drugs);
      #   ### E (Analytical Diagnostic and Therapeutic Techniques and Equipment);F (Psychiatry and Psychology);
      #   ### G (Phenomena and Processes);H (Disciplines and Occupations);
      #   ### I (Anthropology, Education, Sociology and Social Phenomena);J (Technology and Food and Beverages);
      #   ### K (Humanities);L (Information Science);M (Persons);N (Health Care);
      #   ### V (Publication Type);Z (Geographical Locations)
      #
      # library(AnnotationHub)
      # ah <- AnnotationHub()
      # qr <- query(ah, paste(c("MeSHDb for", sp), collapse = " "))
      # if (length(qr) == 0) {
      #   stop("no MeSH records found for ", sp)
      # }
      # qr <- qr[length(qr)]
      # db <- qr[[1]]
      #
      # library("meshr")
      # library("MeSHDbi")
      # MeSH.db <- MeSHDbi::MeSHDb(db)
      # geneid <- keys(MeSH.db, keytype = "GENEID")
      #
      # RSQLite::dbConnect(db)
      # library("MeSHDbi")
      # meshdb <- MeSHDbi::MeSHDb(db)
      # MeSHDbi::dbconn(meshdb)
      # MeSHDbi::dbschema(meshdb)
      # MeSHDbi::dbconn(meshdb)
      #
      #   select(file, keys = keys(file), columns = columns(file))
      #
      #
      #   bg <- AnnotationDbi::select(meshdb, keys = keys(meshdb, keytype = "GENEID"), keytype = "GENEID", columns = c("GENEID", "MESHID", "MESHCATEGORY"))
      #   # bg <-  bg[which(bg$GENEID %in% keys(orgdb)),]
      #   # bg_all <- AnnotationDbi::select(MeSH.db, keys=keys(MeSH.db,keytype = "MESHID"),columns = c("MESHID","MESHTERM"),keytype = "MESHID")
      #   # saveRDS(bg_all,"./MeSHID2Term.rds")
      #   bg_all <- readRDS("./MeSHID2Term.rds")
      #   bg <- merge(x = bg, by.x = "MESHID", y = bg_all, by.y = "MESHID", all.x = TRUE)
      #   bg[which(is.na(bg$MESHTERM)), "MESHTERM"] <- bg[which(is.na(bg$MESHTERM)), "MESHID"]
      #   assign(paste0(sps, "_meshall"), bg)
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
      #   db_list[[sps]][["PFAM"]][["TERM2GENE"]] <- TERM2GENE
      #   db_list[[sps]][["PFAM"]][["TERM2NAME"]] <- TERM2NAME
      #   db_list[[sps]][["PFAM"]][["version"]] <- version
      #   saveCache(db_list[[sps]][["PFAM"]],
      #     key = list(sps, "PFAM"),
      #     comment = paste0(version, " nterm:", length(unique(TERM2NAME[[1]])))
      #   )
      # }
    }

    ### Convert species
    if (!all(db_species == sps)) {
      for (term in names(db_species[db_species != sps])) {
        message("Convert species for the database: ", term)
        sp_from <- db_species[term]
        db_info <- db_list[[sp_from]][[names(sp_from)]]
        TERM2GENE <- db_info[["TERM2GENE"]]
        TERM2NAME <- db_info[["TERM2NAME"]]
        if (is.na(default_IDtypes[term])) {
          IDtype <- colnames(TERM2GENE)[2]
        } else {
          IDtype <- default_IDtypes[[term]]
        }
        if (grepl("MSigDB_", term)) {
          TERM2GENE_map <- db_list[[sps]][["MSigDB"]][["TERM2GENE"]]
          TERM2GENE <- TERM2GENE_map[TERM2GENE_map[["Term"]] %in% TERM2GENE[["Term"]], , drop = FALSE]
          TERM2NAME <- TERM2NAME[TERM2NAME[["Term"]] %in% TERM2GENE[["Term"]], , drop = FALSE]
        } else {
          res <- GeneConvert(
            geneID = as.character(unique(TERM2GENE[, 2])),
            geneID_from_IDtype = IDtype,
            geneID_to_IDtype = "ensembl_id",
            species_from = sp_from,
            species_to = sps,
            Ensembl_version = Ensembl_version,
            mirror = mirror,
            biomart = biomart,
            max_tries = max_tries
          )
          if (is.null(res$geneID_res)) {
            warning("Failed to convert species for the database: ", term, immediate. = TRUE)
            next
          }
          map <- res$geneID_collapse
          TERM2GENE[["ensembl_id-converted"]] <- map[as.character(TERM2GENE[, 2]), "ensembl_id"]
          TERM2GENE <- unnest(TERM2GENE, cols = "ensembl_id-converted", keep_empty = FALSE)
          TERM2GENE <- TERM2GENE[, c("Term", "ensembl_id-converted")]
          colnames(TERM2GENE) <- c("Term", "ensembl_id")
          TERM2NAME <- TERM2NAME[TERM2NAME[["Term"]] %in% TERM2GENE[["Term"]], , drop = FALSE]
        }

        db_info[["TERM2GENE"]] <- unique(TERM2GENE)
        db_info[["TERM2NAME"]] <- unique(TERM2NAME)
        version <- paste0(db_info[["version"]], "(converted from ", sp_from, ")")
        db_info[["version"]] <- version
        db_list[[sps]][[term]] <- db_info
        default_IDtypes[[term]] <- "ensembl_id"
        ### save cache
        saveCache(db_list[[sps]][[term]],
          key = list(version, sps, term),
          comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", sps, "-", term)
        )
      }
    }

    ### Convert ID types
    for (term in names(db_list[[sps]])) {
      IDtypes <- db_IDtypes[!db_IDtypes %in% colnames(db_list[[sps]][[term]][["TERM2GENE"]])]
      if (length(IDtypes) > 0) {
        message("Convert ID types for the database: ", term)
        TERM2GENE <- db_list[[sps]][[term]][["TERM2GENE"]]
        TERM2NAME <- db_list[[sps]][[term]][["TERM2NAME"]]
        if (is.na(default_IDtypes[term])) {
          IDtype <- colnames(TERM2GENE)[2]
        } else {
          IDtype <- default_IDtypes[[term]]
        }
        if (grepl("MSigDB_", term)) {
          map <- db_list[[sps]][["MSigDB"]][["TERM2GENE"]][, -1, drop = FALSE]
          map <- aggregate(map, by = list(map[[1]]), FUN = function(x) list(unique(x)))
          rownames(map) <- map[, 1]
        } else {
          res <- GeneConvert(
            geneID = as.character(unique(TERM2GENE[, 2])),
            geneID_from_IDtype = IDtype,
            geneID_to_IDtype = IDtypes,
            species_from = sps,
            species_to = sps,
            Ensembl_version = Ensembl_version,
            mirror = mirror,
            biomart = biomart,
            max_tries = max_tries
          )
          if (is.null(res$geneID_res)) {
            warning("Failed to convert species for the database: ", term, immediate. = TRUE)
            next
          }
          map <- res$geneID_collapse
        }
        for (type in IDtypes) {
          TERM2GENE[[type]] <- map[as.character(TERM2GENE[, 2]), type]
          TERM2GENE <- unnest(TERM2GENE, cols = type, keep_empty = TRUE)
        }
        db_list[[sps]][[term]][["TERM2GENE"]] <- TERM2GENE
        ### save cache
        version <- db_list[[sps]][[term]][["version"]]
        saveCache(db_list[[sps]][[term]],
          key = list(version, sps, term),
          comment = paste0(version, " nterm:", length(TERM2NAME[[1]]), "|", sps, "-", term)
        )
      }
    }
  }
  return(db_list)
}

#' Perform the enrichment analysis (over-representation) on the genes
#'
#' @inheritParams GeneConvert
#' @param srt A Seurat object containing the results of differential expression analysis (RunDEtest).
#' If specified, the genes and groups will be extracted from the Seurat object automatically.
#' If not specified, the \code{geneID} and \code{geneID_groups} arguments must be provided.
#' @param group_by A character vector specifying the grouping variable in the Seurat object. This argument is only used if \code{srt} is specified.
#' @param test.use A character vector specifying the test to be used in differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param DE_threshold A character vector specifying the filter condition for differential expression analysis. This argument is only used if \code{srt} is specified.
#' @param geneID A character vector specifying the gene IDs.
#' @param geneID_groups A factor vector specifying the group labels for each gene.
#' @param geneID_exclude A character vector specifying the gene IDs to be excluded from the analysis.
#' @param IDtype A character vector specifying the type of gene IDs in the \code{srt} object or \code{geneID} argument. This argument is used to convert the gene IDs to a different type if \code{IDtype} is different from \code{result_IDtype}.
#' @param result_IDtype A character vector specifying the desired type of gene ID to be used in the output. This argument is used to convert the gene IDs from \code{IDtype} to \code{result_IDtype}.
#' @param species A character vector specifying the species for which the analysis is performed.
#' @param db A character vector specifying the name of the database to be used for enrichment analysis.
#' @param db_update A logical value indicating whether the gene annotation databases should be forcefully updated. If set to FALSE, the function will attempt to load the cached databases instead. Default is FALSE.
#' @param db_version A character vector specifying the version of the database to be used. This argument is ignored if \code{db_update} is \code{TRUE}. Default is "latest".
#' @param db_combine A logical value indicating whether to combine multiple databases into one. If TRUE, all database specified by \code{db} will be combined as one named "Combined".
#' @param convert_species A logical value indicating whether to use a species-converted database when the annotation is missing for the specified species. The default value is TRUE.
#' @param TERM2GENE A data frame specifying the gene-term mapping for a custom database. The first column should contain the term IDs, and the second column should contain the gene IDs.
#' @param TERM2NAME A data frame specifying the term-name mapping for a custom database. The first column should contain the term IDs, and the second column should contain the corresponding term names.
#' @param minGSSize A numeric value specifying the minimum size of a gene set to be considered in the enrichment analysis.
#' @param maxGSSize A numeric value specifying the maximum size of a gene set to be considered in the enrichment analysis.
#' @param unlimited_db A character vector specifying the names of databases that do not have size restrictions.
#' @param GO_simplify A logical value indicating whether to simplify the GO terms. If \code{TRUE}, additional results with simplified GO terms will be returned.
#' @param GO_simplify_cutoff A character vector specifying the filter condition for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_method A character vector specifying the method to be used for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param simplify_similarityCutoff A numeric value specifying the similarity cutoff for simplification of GO terms. This argument is only used if \code{GO_simplify} is \code{TRUE}.
#' @param BPPARAM A BiocParallelParam object specifying the parallel back-end to be used for parallel computation. Defaults to BiocParallel::bpparam().
#' @param seed The random seed for reproducibility. Defaults to 11.
#'
#' @returns
#' If input is a Seurat object, returns the modified Seurat object with the enrichment result stored in the tools slot.
#'
#' If input is a geneID vector with or without geneID_groups, return the enrichment result directly.
#'
#' Enrichment result is a list with the following component:
#' \itemize{
#'  \item{\code{enrichment}:}{ A data.frame containing all enrichment results.}
#'  \item{\code{results}:}{ A list of \code{enrichResult} objects from the DOSE package.}
#'  \item{\code{geneMap}:}{ A data.frame containing the ID mapping table for input gene IDs.}
#'  \item{\code{input}:}{ A data.frame containing the input gene IDs and gene ID groups.}
#'  \item{\code{DE_threshold}:}{ A specific threshold for differential expression analysis (only returned if input is a Seurat object).}
#' }
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}} \code{\link{EnrichmentPlot}} \code{\link{RunGSEA}} \code{\link{GSEAPlot}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunDEtest(pancreas_sub, group_by = "CellType")
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub, group_by = "CellType", DE_threshold = "p_val_adj < 0.05",
#'   db = "GO_BP", species = "Mus_musculus"
#' )
#' EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType", plot_type = "comparison")
#'
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub, group_by = "CellType", DE_threshold = "p_val_adj < 0.05",
#'   db = c("MSigDB", "MSigDB_MH"), species = "Mus_musculus"
#' )
#' EnrichmentPlot(pancreas_sub, db = "MSigDB", group_by = "CellType", plot_type = "comparison")
#' EnrichmentPlot(pancreas_sub, db = "MSigDB_MH", group_by = "CellType", plot_type = "comparison")
#'
#' # Remove redundant GO terms
#' pancreas_sub <- RunEnrichment(srt = pancreas_sub, group_by = "CellType", db = "GO_BP", GO_simplify = TRUE, species = "Mus_musculus")
#' EnrichmentPlot(pancreas_sub, db = "GO_BP_sim", group_by = "CellType", plot_type = "comparison")
#'
#' # Use a combined database
#' pancreas_sub <- RunEnrichment(
#'   srt = pancreas_sub, group_by = "CellType",
#'   db = c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
#'   db_combine = TRUE,
#'   species = "Mus_musculus"
#' )
#' EnrichmentPlot(pancreas_sub, db = "Combined", group_by = "CellType", plot_type = "comparison")
#'
#' # Or use "geneID" and "geneID_groups" as input to run enrichment
#' de_df <- dplyr::filter(pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox, p_val_adj < 0.05)
#' enrich_out <- RunEnrichment(geneID = de_df[["gene"]], geneID_groups = de_df[["group1"]], db = "GO_BP", species = "Mus_musculus")
#' EnrichmentPlot(res = enrich_out, db = "GO_BP", plot_type = "comparison")
#'
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers ipcid ipclock ipcunlock
#' @importFrom clusterProfiler enricher simplify
#' @export
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
RunSlingshot <- function(srt, group.by, reduction = NULL, dims = NULL, start = NULL, end = NULL, prefix = NULL,
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
  if (is.null(dims)) {
    dims <- 1:2
  }

  set.seed(seed)
  sl <- slingshot(
    data = as.data.frame(srt_sub[[reduction]]@cell.embeddings[, dims]),
    clusterLabels = as.character(srt_sub[[group.by, drop = TRUE]]),
    start.clus = start, end.clus = end, ...
  )

  srt@tools[[paste("Slingshot", group.by, reduction, sep = "_")]] <- sl
  df <- as.data.frame(slingPseudotime(sl))
  colnames(df) <- paste0(prefix, colnames(df))
  if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
      df <- max(df, na.rm = TRUE) - df
    }
  }
  srt <- AddMetaData(srt, metadata = df)
  srt <- AddMetaData(srt, metadata = slingBranchID(sl), col.name = paste0(prefix, "BranchID"))

  if (isTRUE(show_plot)) {
    if (ncol(srt[[reduction]]@cell.embeddings) == 2 || ncol(srt[[reduction]]@cell.embeddings) > 3) {
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 2, type = "lineages", col = "black")
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scp(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 3, col = 1:length(slingshot::SlingshotDataSet(sl)@lineages))
      p <- CellDimPlot(srt, group.by = group.by, reduction = reduction, dims = c(1, 2), lineages = colnames(df))
      print(p)
    } else if (ncol(srt[[reduction]]@cell.embeddings) == 3) {
      p <- CellDimPlot3D(srt, group.by = group.by, reduction = reduction, lineages = colnames(df))
      print(p)
    }
  }
  return(srt)
}

#' Run Monocle2 analysis
#'
#' Runs the Monocle2 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param slot The slot in the Seurat object to use for analysis. Default is "counts".
#' @param expressionFamily The distribution family to use for modeling gene expression. Default is "negbinomial.size".
#' @param features A vector of gene names or indices specifying the features to use in the analysis. Defaults to NULL, in which case features were determined by \code{feature_type}.
#' @param feature_type The type of features to use in the analysis. Possible values are "HVF" for highly variable features
#'   or "Disp" for features selected based on dispersion. Default is "HVF".
#' @param disp_filter A string specifying the filter to use when \code{feature_type} is "Disp". Default is
#'   "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit".
#' @param max_components The maximum number of dimensions to use for dimensionality reduction. Default is 2.
#' @param reduction_method The dimensionality reduction method to use. Possible values are "DDRTree" and "UMAP". Default is "DDRTree".
#' @param norm_method The normalization method to use. Possible values are "log" and "none". Default is "log".
#' @param residualModelFormulaStr A model formula specifying the effects to subtract. Default is NULL.
#' @param pseudo_expr Amount to increase expression values before dimensionality reduction. Default is 1.
#' @param root_state The state to use as the root of the trajectory. If NULL, will prompt for user input.
#' @param seed An integer specifying the random seed to use. Default is 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   pancreas_sub <- RunMonocle2(srt = pancreas_sub)
#'   names(pancreas_sub@tools$Monocle2)
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#'
#'   pancreas_sub <- RunMonocle2(
#'     srt = pancreas_sub,
#'     feature_type = "Disp", disp_filter = "mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#' }
#' @importFrom Seurat DefaultAssay CreateDimReducObject GetAssayData VariableFeatures FindVariableFeatures
#' @importFrom SeuratObject as.sparse
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @importFrom utils select.list
#' @export
RunMonocle2 <- function(srt, assay = NULL, slot = "counts", expressionFamily = "negbinomial.size",
                        features = NULL, feature_type = "HVF", disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit",
                        max_components = 2, reduction_method = "DDRTree", norm_method = "log", residualModelFormulaStr = NULL, pseudo_expr = 1,
                        root_state = NULL, seed = 11) {
  set.seed(seed)
  check_R(c("monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM"))

  if (!"package:DDRTree" %in% search()) {
    attachNamespace("DDRTree")
  }

  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- as.sparse(GetAssayData(srt, assay = assay, slot = slot))
  p_data <- srt@meta.data
  f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
  pd <- new("AnnotatedDataFrame", data = p_data)
  fd <- new("AnnotatedDataFrame", data = f_data)
  cds <- monocle::newCellDataSet(expr_matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily = do.call(get(expressionFamily, envir = getNamespace("VGAM")), args = list())
  )
  if (any(c("negbinomial", "negbinomial.size") %in% expressionFamily)) {
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- BiocGenerics::estimateDispersions(cds)
  }
  if (is.null(features)) {
    if (feature_type == "HVF") {
      features <- VariableFeatures(srt, assay = assay)
      if (length(features) == 0) {
        features <- VariableFeatures(FindVariableFeatures(srt, assay = assay), assay = assay)
      }
    }
    if (feature_type == "Disp") {
      features <- subset(monocle::dispersionTable(cds), eval(rlang::parse_expr(disp_filter)))$gene_id
    }
  }
  message("features number: ", length(features))
  cds <- monocle::setOrderingFilter(cds, features)
  p <- monocle::plot_ordering_genes(cds)
  print(p)

  cds <- monocle::reduceDimension(
    cds = cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr
  )
  cds <- orderCells(cds)

  embeddings <- t(cds@reducedDimS)
  colnames(embeddings) <- paste0(cds@dim_reduce_type, "_", 1:ncol(embeddings))
  srt[[cds@dim_reduce_type]] <- CreateDimReducObject(embeddings = embeddings, key = paste0(cds@dim_reduce_type, "_"), assay = assay)
  srt[["Monocle2_State"]] <- cds[["State"]]
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimS))
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimK))
  }
  edge_df <- as_data_frame(cds@minSpanningTree)
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  trajectory <- geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
  p <- CellDimPlot(srt, group.by = "Monocle2_State", reduction = reduction_method, label = TRUE, force = TRUE) +
    trajectory
  print(p)
  if (is.null(root_state)) {
    root_state <- select.list(sort(unique(cds[["State"]])), title = "Select the root state to order cells:")
    if (root_state == "" || length(root_state) == 0) {
      root_state <- NULL
    }
  }
  cds <- orderCells(cds, root_state = root_state)
  srt[["Monocle2_State"]] <- cds[["State"]]
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(cds = cds, features = features, trajectory = trajectory)

  p1 <- CellDimPlot(srt, group.by = "Monocle2_State", reduction = reduction_method, label = TRUE, force = TRUE) + trajectory
  p2 <- FeatureDimPlot(srt, features = "Monocle2_Pseudotime", reduction = reduction_method) + trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}

orderCells <- function(cds, root_state = NULL, num_paths = NULL, reverse = NULL) {
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths <- 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- as_matrix(stats::dist(adjusted_S))
    cellPairwiseDistances(cds) <- as_matrix(stats::dist(adjusted_S))
    gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::minimum.spanning.tree(gp)
    monocle::minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- monocle:::pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- monocle:::extract_good_branched_ordering(res$subtree, res$root, monocle::cellPairwiseDistances(cds), num_paths, FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    monocle::minSpanningTree(cds) <- igraph::as.undirected(order_list$cell_ordering_tree)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
    monocle::minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- igraph::as.undirected(order_list$cell_ordering_tree)
  } else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    K_old <- monocle::reducedDimK(cds)
    old_dp <- monocle::cellPairwiseDistances(cds)
    old_mst <- monocle::minSpanningTree(cds)
    old_A <- monocle::reducedDimA(cds)
    old_W <- monocle::reducedDimW(cds)
    cds <- project2MST(cds, monocle:::project_point_to_line_segment)
    monocle::minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(igraph::V(old_mst)$name == root_cell, arr.ind = TRUE)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- igraph::V(monocle::minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(igraph::degree(monocle::minSpanningTree(cds)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering_new_pseudotime[row.names(Biobase::pData(cds)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      cds[["State"]] <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
    cds@reducedDimK <- K_old
    cds@cellPairwiseDistances <- old_dp
    cds@minSpanningTree <- old_mst
    cds@reducedDimA <- old_A
    cds@reducedDimW <- old_W
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(monocle::minSpanningTree(cds))[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}
project2MST <- function(cds, Projection_Method) {
  dp_mst <- monocle::minSpanningTree(cds)
  Z <- monocle::reducedDimS(cds)
  Y <- monocle::reducedDimK(cds)
  cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as_matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(igraph::V(dp_mst)[suppressWarnings(nei(closest_vertex_names[i], mode = "all"))])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- monocle:::projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      if (!inherits(projection, "matrix")) {
        projection <- as_matrix(projection)
      }
      P[, i] <- projection[which(distance == min(distance))[1], ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- as_matrix(stats::dist(t(P)))
  min_dist <- min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  monocle::cellPairwiseDistances(cds) <- dp
  gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- igraph::minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}
extract_ddrtree_ordering <- function(cds, root_cell, verbose = TRUE) {
  dp <- monocle::cellPairwiseDistances(cds)
  dp_mst <- monocle::minSpanningTree(cds)
  curr_state <- 1
  res <- list(subtree = dp_mst, root = root_cell)
  states <- rep(1, ncol(dp))
  names(states) <- igraph::V(dp_mst)$name
  pseudotimes <- rep(0, ncol(dp))
  names(pseudotimes) <- igraph::V(dp_mst)$name
  parents <- rep(NA, ncol(dp))
  names(parents) <- igraph::V(dp_mst)$name
  mst_traversal <- igraph::graph.dfs(dp_mst,
    root = root_cell, mode = "all",
    unreachable = FALSE, father = TRUE
  )
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  for (i in 1:length(mst_traversal$order)) {
    curr_node <- mst_traversal$order[i]
    curr_node_name <- igraph::V(dp_mst)[curr_node]$name
    if (is.na(mst_traversal$father[curr_node]) == FALSE) {
      parent_node <- mst_traversal$father[curr_node]
      parent_node_name <- igraph::V(dp_mst)[parent_node]$name
      parent_node_pseudotime <- pseudotimes[parent_node_name]
      parent_node_state <- states[parent_node_name]
      curr_node_pseudotime <- parent_node_pseudotime +
        dp[curr_node_name, parent_node_name]
      if (igraph::degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
    } else {
      parent_node <- NA
      parent_node_name <- NA
      curr_node_pseudotime <- 0
    }
    curr_node_state <- curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  ordering_df <- data.frame(
    sample_name = names(states), cell_state = factor(states),
    pseudo_time = as.vector(pseudotimes), parent = parents
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

#' Run Monocle3 analysis
#'
#' Runs the Monocle3 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param slot The slot in the Seurat object to use for analysis. Default is "counts".
#' @param reduction The reduction used. Defaults to NULL, in which case the default reduction of the Seurat object is used.
#' @param clusters The cluster variable in the Seurat object to use for analysis. Defaults to NULL, in which case use Monocle clusters is used.
#' @param graph The name of the graph slot in the Seurat object to use for analysis. Defaults to NULL, in which case Monocle graph is used.
#' @param partition_qval The q-value threshold for partitioning cells. Defaults to 0.05.
#' @param k The number of nearest neighbors to consider for clustering. Defaults to 50.
#' @param cluster_method The clustering method to use. Defaults to "louvain".
#' @param num_iter The number of iterations for clustering. Defaults to 2.
#' @param resolution The resolution parameter for clustering. Defaults to NULL.
#' @param use_partition Whether to use partitions to learn disjoint graph in each partition. If not specified, user will be prompted for input. Defaults to NULL.
#' @param close_loop Whether to close loops in the graph. Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param root_cells The root cells to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param seed The random seed to use for reproducibility. Defaults to 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   # Use Monocle clusters to infer the trajectories
#'   pancreas_sub <- RunMonocle3(srt = pancreas_sub, reduction = "UMAP")
#'   names(pancreas_sub@tools$Monocle3)
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   milestones <- pancreas_sub@tools$Monocle3$milestones
#'
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory + milestones
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "UMAP", theme_use = "theme_blank") + trajectory
#'
#'   ## Select the lineage using monocle3::choose_graph_segments
#'   # cds <- pancreas_sub@tools$Monocle3$cds
#'   # cds_sub <- monocle3::choose_graph_segments(cds, starting_pr_node = NULL, ending_pr_nodes = NULL)
#'   # pancreas_sub$Lineages_1 <- NA
#'   # pancreas_sub$Lineages_1[colnames(cds_sub)] <- pancreas_sub$Monocle3_Pseudotime[colnames(cds_sub)]
#'   # CellDimPlot(pancreas_sub, group.by = "SubCellType", lineages = "Lineages_1", lineages_span = 0.1, theme_use = "theme_blank")
#'
#'   # Use Seurat clusters to infer the trajectories
#'   pancreas_sub <- Standard_SCP(pancreas_sub)
#'   CellDimPlot(pancreas_sub, group.by = c("Standardclusters", "CellType"), label = TRUE, theme_use = "theme_blank")
#'   pancreas_sub <- RunMonocle3(srt = pancreas_sub, clusters = "Standardclusters")
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "StandardUMAP2D", theme_use = "theme_blank") + trajectory
#'
#'   # Use custom graphs and cell clusters to infer the partitions and trajectories, respectively
#'   pancreas_sub <- Standard_SCP(pancreas_sub, cluster_resolution = 5)
#'   CellDimPlot(pancreas_sub, group.by = c("Standardclusters", "CellType"), label = TRUE)
#'   pancreas_sub <- RunMonocle3(
#'     srt = pancreas_sub,
#'     clusters = "Standardclusters", graph = "Standardpca_SNN"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_partitions", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle3_clusters", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
#'   FeatureDimPlot(pancreas_sub, features = "Monocle3_Pseudotime", reduction = "StandardUMAP2D", theme_use = "theme_blank") + trajectory
#' }
#' @importFrom SeuratObject as.sparse Embeddings Loadings Stdev
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color
#' @importFrom utils packageVersion select.list
#' @export
RunMonocle3 <- function(srt, assay = NULL, slot = "counts",
                        reduction = DefaultReduction(srt), clusters = NULL, graph = NULL, partition_qval = 0.05,
                        k = 50, cluster_method = "louvain", num_iter = 2, resolution = NULL,
                        use_partition = NULL, close_loop = TRUE,
                        root_pr_nodes = NULL, root_cells = NULL, seed = 11) {
  set.seed(seed)
  if (!requireNamespace("monocle3", quietly = TRUE) || packageVersion("monocle3") < package_version("1.2.0")) {
    check_R("cole-trapnell-lab/monocle3", force = TRUE)
  }
  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- as.sparse(GetAssayData(srt, assay = assay, slot = slot))
  p_data <- srt@meta.data
  f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
  cds <- monocle3::new_cell_data_set(
    expression_data = expr_matrix,
    cell_metadata = p_data,
    gene_metadata = f_data
  )
  if (!"Size_Factor" %in% colnames(cds@colData)) {
    size.factor <- paste0("nCount_", assay)
    if (size.factor %in% colnames(srt@meta.data)) {
      cds[["Size_Factor"]] <- cds[[size.factor, drop = TRUE]]
    }
  }
  reduction <- reduction %||% DefaultReduction(srt)
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- Embeddings(srt[[reduction]])
  loadings <- Loadings(object = srt[[reduction]])
  if (length(loadings) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["gene_loadings"]] <- loadings
  }
  stdev <- Stdev(object = srt[[reduction]])
  if (length(stdev) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["prop_var_expl"]] <- stdev
  }

  if (!is.null(clusters)) {
    if (!is.null(graph)) {
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = srt[[graph]],
        weighted = TRUE
      )
      cluster_result <- list(
        g = g,
        relations = NULL,
        distMatrix = "matrix",
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(as.factor(srt[[clusters, drop = TRUE]])),
          modularity = NA_real_
        )
      )
      if (length(unique(cluster_result$optim_res$membership)) > 1) {
        cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, cluster_result$optim_res, partition_qval)
        partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
        partitions <- as.factor(partitions)
      } else {
        partitions <- rep(1, ncol(srt))
      }
      names(partitions) <- colnames(cds)
      cds@clusters[["UMAP"]] <- list(
        cluster_result = cluster_result,
        partitions = partitions,
        clusters = as.factor(srt[[clusters, drop = TRUE]])
      )
      cds[["clusters"]] <- cds[[clusters]]
      cds <- monocle3:::add_citation(cds, "clusters")
      cds <- monocle3:::add_citation(cds, "partitions")
    } else {
      cds <- monocle3::cluster_cells(cds,
        reduction_method = "UMAP", partition_qval = partition_qval,
        k = k, cluster_method = cluster_method, num_iter = num_iter, resolution = resolution
      )
      cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters <- as.factor(srt[[clusters, drop = TRUE]])
    }
  } else {
    cds <- monocle3::cluster_cells(cds,
      reduction_method = "UMAP", partition_qval = partition_qval,
      k = k, cluster_method = cluster_method, num_iter = num_iter, resolution = resolution
    )
    cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
  }
  srt[["Monocle3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
  srt[["Monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions
  p1 <- CellDimPlot(srt, "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE)
  p2 <- CellDimPlot(srt, "Monocle3_clusters", reduction = reduction, label = FALSE, force = TRUE)
  print(p1)
  Sys.sleep(1)
  print(p2)
  if (is.null(use_partition)) {
    use_partition <- select.list(c(TRUE, FALSE), title = "Whether to use partitions to learn disjoint graph in each partition?")
    if (use_partition == "" || length(use_partition) == 0) {
      use_partition <- TRUE
    }
  }
  cds <- monocle3::learn_graph(cds = cds, use_partition = use_partition, close_loop = close_loop)

  reduced_dim_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  edge_df <- as_data_frame(cds@principal_graph[["UMAP"]])
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
  mst_leaf_nodes <- monocle3:::leaf_nodes(cds, "UMAP")
  mst_root_nodes <- monocle3:::root_nodes(cds, "UMAP")
  pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
  point_df <- data.frame(nodes = names(pps), x = reduced_dim_coords[pps, 1], y = reduced_dim_coords[pps, 2])
  point_df[, "is_branch"] <- names(pps) %in% names(mst_branch_nodes)
  trajectory <- list(
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
  )
  milestones <- list(
    geom_point(
      data = point_df[point_df[["is_branch"]] == FALSE, , drop = FALSE], aes(x = x, y = y),
      shape = 21, color = "white", fill = "black", size = 3, stroke = 1
    ),
    geom_point(
      data = point_df[point_df[["is_branch"]] == TRUE, , drop = FALSE], aes(x = x, y = y),
      shape = 21, color = "white", fill = "red", size = 3, stroke = 1
    ),
    new_scale_color(),
    geom_text_repel(
      data = point_df, aes(x = x, y = y, label = nodes, color = is_branch),
      fontface = "bold", min.segment.length = 0,
      point.size = 3, max.overlaps = 100,
      bg.color = "white", bg.r = 0.1, size = 3.5
    ),
    scale_color_manual(values = setNames(c("red", "black"), nm = c(TRUE, FALSE)))
  )
  p <- CellDimPlot(srt, group.by = "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE) +
    trajectory + milestones
  print(p)

  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    root_pr_nodes <- select.list(names(pps), title = "Select the root nodes to order cells, or leave blank for interactive selection:", multiple = TRUE)
    if (root_pr_nodes == "" || length(root_pr_nodes) == 0) {
      root_pr_nodes <- NULL
    }
  }
  cds <- monocle3::order_cells(cds, root_pr_nodes = root_pr_nodes, root_cells = root_cells)
  pseudotime <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime[is.infinite(pseudotime)] <- NA
  srt[["Monocle3_Pseudotime"]] <- pseudotime
  srt@tools$Monocle3 <- list(cds = cds, trajectory = trajectory, milestones = milestones)

  p1 <- CellDimPlot(srt, group.by = "Monocle3_partitions", reduction = reduction, label = FALSE, force = TRUE) +
    trajectory
  p2 <- FeatureDimPlot(srt, features = "Monocle3_Pseudotime", reduction = reduction) +
    theme(legend.position = "none") +
    trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}

#' RunDynamicFeatures
#'
#' Calculates dynamic features for lineages in a single-cell RNA-seq dataset.
#'
#' @param srt A Seurat object.
#' @param lineages A character vector specifying the lineage names for which dynamic features should be calculated.
#' @param features A character vector specifying the features (genes or metadata variables) for which dynamic features should be calculated. If NULL, n_candidates must be provided.
#' @param suffix A character vector specifying the suffix to append to the output slot names for each lineage. Defaults to the lineage names.
#' @param n_candidates An integer specifying the number of candidate features to select when features is NULL. Defaults to 1000.
#' @param minfreq An integer specifying the minimum frequency threshold for candidate features. Features with a frequency less than minfreq will be excluded. Defaults to 5.
#' @param family A character or character vector specifying the family of distributions to use for the generalized additive models (GAMs). If family is set to NULL, the appropriate family will be automatically determined based on the data. If length(family) is 1, the same family will be used for all features. Otherwise, family must have the same length as features.
#' @param slot A character vector specifying the slot in the Seurat object to use. Default is "counts".
#' @param assay A character vector specifying the assay in the Seurat object to use. Default is NULL.
#' @param libsize A numeric or numeric vector specifying the library size correction factors for each cell. If NULL, the library size correction factors will be calculated based on the expression matrix. If length(libsize) is 1, the same value will be used for all cells. Otherwise, libsize must have the same length as the number of cells in srt. Defaults to NULL.
#' @param BPPARAM A BiocParallelParam object specifying the parallelization parameters for computing dynamic features. Defaults to BiocParallel::bpparam().
#' @param seed An integer specifying the seed for random number generation. Defaults to 11.
#'
#' @return Returns the modified Seurat object with the calculated dynamic features stored in the tools slot.
#'
#' @seealso \code{\link{DynamicHeatmap}} \code{\link{DynamicPlot}} \code{\link{RunDynamicEnrichment}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
#' names(pancreas_sub@tools$DynamicFeatures_Lineage1)
#' head(pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures)
#' ht <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 6, reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   srt = pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Nnat", "Irx1"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers
#' @export
RunDynamicFeatures <- function(srt, lineages, features = NULL, suffix = lineages,
                               n_candidates = 1000, minfreq = 5,
                               family = NULL,
                               slot = "counts", assay = NULL, libsize = NULL,
                               BPPARAM = BiocParallel::bpparam(), seed = 11) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  assay <- assay %||% DefaultAssay(srt)

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start RunDynamicFeatures"))
  message("Workers: ", bpworkers(BPPARAM))

  check_R("mgcv")
  meta <- c()
  gene <- c()
  if (!is.null(features)) {
    gene <- features[features %in% rownames(srt[[assay]])]
    meta <- features[features %in% colnames(srt@meta.data)]
    isnum <- sapply(srt@meta.data[, meta, drop = FALSE], is.numeric)
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
    Y <- rbind(Y, t(srt@meta.data[, meta, drop = FALSE]))
  }

  features_list <- c()
  srt_sub_list <- list()
  for (l in lineages) {
    srt_sub <- subset(srt, cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])])
    if (is.null(features)) {
      if (is.null(n_candidates)) {
        stop("'features' or 'n_candidates' must provided at least one.")
      }
      HVF <- VariableFeatures(FindVariableFeatures(srt_sub, nfeatures = n_candidates, assay = assay), assay = assay)
      HVF_counts <- srt_sub[[assay]]@counts[HVF, , drop = FALSE]
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
      stop("'family' must be one character or a vector of the same length as features.")
    }
  }

  for (i in seq_along(lineages)) {
    l <- lineages[i]
    srt_sub <- srt_sub_list[[l]]
    t <- srt_sub[[l, drop = TRUE]]
    t <- t[is.finite(t)]
    t_ordered <- t[order(t)]
    Y_ordered <- as_matrix(Y[features, names(t_ordered), drop = FALSE])
    l_libsize <- Y_libsize[names(t_ordered)]
    raw_matrix <- as_matrix(cbind(data.frame(pseudotime = t_ordered), t(Y_ordered)))

    # df <- data.frame(x = rowMeans(Y_ordered), y = MatrixGenerics::rowVars(Y_ordered))
    # p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
    #   geom_point() +
    #   geom_abline(slope = 1, intercept = 0, color = "red") +
    #   labs(title = l, subtitle = "Var VS Mean") +
    #   theme_scp()
    # print(p)

    message("Calculate dynamic features for ", l, "...")
    system.time({
      gam_out <- bplapply(seq_len(nrow(Y_ordered)), function(n, Y_ordered, t_ordered, l_libsize, family) {
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
        peaktime <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
        valleytime <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])

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
          peaktime = peaktime, valleytime = valleytime,
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
  message(paste0("[", time_end, "] ", "RunDynamicFeatures done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  return(srt)
}

#' RunDynamicEnrichment
#'
#' This function calculates gene-set scores from the specified database (\code{db}) for each lineage using the specified scoring method (\code{score_method}).
#' It then treats these scores as expression values and uses them as input to the RunDynamicFeatures function to identify dynamically enriched terms along the lineage.
#'
#' @inheritParams RunEnrichment
#' @inheritParams DynamicHeatmap
#' @inheritParams CellScoring
#' @param score_method The method to use for scoring. Can be "Seurat", "AUCell", or "UCell". Defaults to "Seurat".
#'
#' @seealso \code{\link{RunDynamicFeatures}} \code{\link{DynamicHeatmap}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSlingshot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP")
#' pancreas_sub <- RunDynamicFeatures(pancreas_sub, lineages = "Lineage1", n_candidates = 200)
#' ht1 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   cell_annotation = "SubCellType",
#'   n_split = 4
#' )
#' ht1$plot
#'
#' pancreas_sub <- RunDynamicEnrichment(
#'   srt = pancreas_sub,
#'   lineages = "Lineage1",
#'   score_method = "AUCell",
#'   db = "GO_BP",
#'   species = "Mus_musculus"
#' )
#' ht2 <- DynamicHeatmap(
#'   srt = pancreas_sub,
#'   assay = "GO_BP",
#'   lineages = "Lineage1_GO_BP",
#'   cell_annotation = "SubCellType",
#'   n_split = 4,
#'   split_method = "kmeans-peaktime"
#' )
#' ht2$plot
#' @importFrom Seurat NormalizeData VariableFeatures FindVariableFeatures as.SingleCellExperiment AddMetaData
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply bpprogressbar<- bpRNGseed<- bpworkers
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @export
RunDynamicEnrichment <- function(srt, lineages,
                                 score_method = "AUCell",
                                 slot = "data", assay = NULL,
                                 min_expcells = 20, r.sq = 0.2, dev.expl = 0.2, padjust = 0.05,
                                 IDtype = "symbol", species = "Homo_sapiens",
                                 db = "GO_BP", db_update = FALSE, db_version = "latest", convert_species = TRUE,
                                 Ensembl_version = 103, mirror = NULL,
                                 TERM2GENE = NULL, TERM2NAME = NULL, minGSSize = 10, maxGSSize = 500,
                                 BPPARAM = BiocParallel::bpparam(), seed = 11) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed
  assay <- assay %||% DefaultAssay(srt)

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start RunDynamicFeatures"))
  message("Workers: ", bpworkers(BPPARAM))

  feature_union <- c()
  cell_union <- c()
  dynamic <- list()
  for (l in lineages) {
    if (!paste0("DynamicFeatures_", l) %in% names(srt@tools)) {
      stop(l, " info not found in the srt object. Should perform RunDynamicFeatures first!")
    }
    DynamicFeatures <- srt@tools[[paste0("DynamicFeatures_", l)]][["DynamicFeatures"]]
    DynamicFeatures <- DynamicFeatures[DynamicFeatures$exp_ncells > min_expcells & DynamicFeatures$r.sq > r.sq & DynamicFeatures$dev.expl > dev.expl & DynamicFeatures$padjust < padjust, , drop = FALSE]
    dynamic[[l]] <- DynamicFeatures
    feature_union <- c(feature_union, DynamicFeatures[, "features"])
    cell_union <- c(cell_union, rownames(srt@tools[[paste0("DynamicFeatures_", l)]][["raw_matrix"]]))
  }
  feature_union <- unique(feature_union)

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
  }

  for (i in seq_along(db)) {
    term <- db[i]
    TERM2GENE_tmp <- db_list[[species]][[term]][["TERM2GENE"]][, c("Term", IDtype)]
    TERM2NAME_tmp <- db_list[[species]][[term]][["TERM2NAME"]]
    dup <- duplicated(TERM2GENE_tmp)
    na <- rowSums(is.na(TERM2GENE_tmp)) > 0
    TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], , drop = FALSE]

    term_use <- unique(TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% feature_union, "Term"])
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, "Term"] %in% term_use, , drop = FALSE]
    TERM2NAME_tmp <- TERM2NAME_tmp[TERM2NAME_tmp[, "Term"] %in% term_use, , drop = FALSE]
    rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
    TERM2GENE_tmp <- TERM2GENE_tmp[TERM2GENE_tmp[, IDtype] %in% rownames(srt[[assay]]), , drop = FALSE]
    feature_list <- split(TERM2GENE_tmp[, IDtype], TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"])
    GSSize <- sapply(feature_list, length)
    feature_list <- feature_list[GSSize >= minGSSize & GSSize <= maxGSSize]

    srt <- CellScoring(
      srt = srt,
      features = feature_list,
      method = score_method,
      classification = FALSE,
      slot = slot,
      assay = assay,
      name = term,
      new_assay = TRUE,
      BPPARAM = BPPARAM
    )
    srt <- RunDynamicFeatures(
      srt = srt,
      lineages = lineages,
      features = rownames(srt[[term]]@counts),
      suffix = paste(lineages, term, sep = "_"),
      assay = term
    )
  }

  time_end <- Sys.time()
  message(paste0("[", time_end, "] ", "RunDynamicEnrichment done"))
  message("Elapsed time:", format(round(difftime(time_end, time_start), 2), format = "%Y-%m-%d %H:%M:%S"))

  return(srt)
}

#' @importFrom reticulate py_to_r
py_to_r_auto <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    x <- py_to_r(x)
  }
  return(x)
}

#' Convert a seurat object to an anndata object using reticulate
#'
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @param srt A Seurat object.
#' @param assay_X Assay to convert as the main data matrix (X) in the anndata object.
#' @param slot_X Slot name for assay_X in the Seurat object.
#' @param assay_layers Assays to convert as layers in the anndata object.
#' @param slot_layers Slot names for the assay_layers in the Seurat object.
#' @param convert_tools Logical indicating whether to convert the tool-specific data.
#' @param convert_misc Logical indicating whether to convert the miscellaneous data.
#' @param features Optional vector of features to include in the anndata object. Defaults to all features in assay_X.
#' @param verbose Logical indicating whether to print verbose messages during the conversion process.
#'
#' @return An anndata object.
#'
#' @return A \code{anndata} object.
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ### Or save as an h5ad file or a loom file
#' # adata$write_h5ad("pancreas_sub.h5ad")
#' # adata$write_loom("pancreas_sub.loom", write_obsm_varm = TRUE)
#'
#' @importFrom reticulate import np_array
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @export
srt_to_adata <- function(srt, features = NULL,
                         assay_X = "RNA", slot_X = "counts",
                         assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                         convert_tools = FALSE, convert_misc = FALSE, verbose = TRUE) {
  check_Python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }

  if (is.null(features)) {
    features <- rownames(srt[[assay_X]])
  }
  if (length(slot_layers) == 1) {
    slot_layers <- rep(slot_layers, length(assay_layers))
    names(slot_layers) <- assay_layers
  } else if (length(slot_layers) != length(assay_layers)) {
    stop("slot_layers must be one character or the same length of the assay_layers")
  }

  sc <- import("scanpy", convert = FALSE)
  np <- import("numpy", convert = FALSE)

  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }

  var <- srt[[assay_X]]@meta.features[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(as.character(var[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }
  if (length(VariableFeatures(srt, assay = assay_X) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in% VariableFeatures(srt, assay = assay_X)
  }

  # X <- t(as_matrix(slot(srt@assays[[assay_X]], slot_X)[features, , drop = FALSE]))
  X <- t(GetAssayData(srt, assay = assay_X, slot = slot_X)[features, , drop = FALSE])
  adata <- sc$AnnData(
    X = np_array(X, dtype = np$float32),
    obs = obs,
    var = cbind(data.frame(features = features), var)
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_X]) {
    if (assay %in% assay_layers) {
      layer <- t(GetAssayData(srt, assay = assay, slot = slot_layers[assay]))
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          stop(
            "The following features in the '", assay_X, "' assay can not be found in the '", assay, "' assay:\n  ",
            paste0(head(colnames(X)[!colnames(X) %in% colnames(layer)], 10), collapse = ","), "..."
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      if (isTRUE(verbose)) {
        message("Assay '", assay, "' is in the srt object but not converted.")
      }
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

  # varm_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   if (ncol(srt[[reduction]]@feature.loadings) > 0) {
  #     varm_list[[paste0(reduction, "_feature.loadings")]] <- srt[[reduction]]@feature.loadings
  #   }
  # }
  # if (length(varm_list) > 0) {
  #   adata$varm <- varm_list
  # }

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
    if (isTRUE(verbose)) {
      message("'misc' slot is not converted.")
    }
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'tools' slot is not converted.")
    }
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}

#' Convert an anndata object to a seurat object using reticulate
#'
#' @param adata a connected python anndata object.
#'
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(adata = adata, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' srt <- adata_to_srt(adata)
#' srt
#'
#' ### Or convert a h5ad file to Seurat object
#' # library(reticulate)
#' # check_Python("scanpy")
#' # sc <- import("scanpy")
#' # adata <- sc$read_h5ad("pancreas.h5ad")
#' # srt <- adata_to_srt(adata)
#' # srt
#' @importFrom Seurat CreateSeuratObject CreateAssayObject CreateDimReducObject AddMetaData
#' @importFrom SeuratObject as.Graph as.sparse
#' @importFrom reticulate iterate
#' @importFrom Matrix t
#' @export
adata_to_srt <- function(adata) {
  if (!inherits(adata, "python.builtin.object")) {
    stop("'adata' is not a python.builtin.object.")
  }
  sc <- import("scanpy", convert = TRUE)
  np <- import("numpy", convert = TRUE)

  x <- t(py_to_r_auto(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- as.sparse(x[1:nrow(x), , drop = FALSE])
  }
  rownames(x) <- py_to_r_auto(adata$var_names$values)
  colnames(x) <- py_to_r_auto(adata$obs_names$values)

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(py_to_r_auto(adata$obs))
    colnames(metadata) <- make.names(colnames(metadata))
  }

  srt <- CreateSeuratObject(counts = x, meta.data = metadata)

  if (inherits(adata$layers, "python.builtin.object")) {
    keys <- iterate(adata$layers$keys())
  } else {
    keys <- names(adata$layers)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      layer <- py_to_r_auto(adata$layers[[k]])
      if (!inherits(layer, c("Matrix", "matrix"))) {
        stop(paste0("The object in '", k, "' layers is not a matrix: ", paste0(class(adata$layers[[k]]), collapse = ",")))
      }
      layer <- t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- as.sparse(layer[1:nrow(layer), , drop = FALSE])
      }
      rownames(layer) <- py_to_r_auto(adata$var_names$values)
      colnames(layer) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateAssayObject(counts = layer)
    }
  }

  if (inherits(adata$obsm, "python.builtin.object")) {
    keys <- iterate(adata$obsm$keys())
  } else {
    keys <- names(adata$obsm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsm <- tryCatch(py_to_r_auto(adata$obsm[[k]]), error = identity)
      if (inherits(obsm, "error")) {
        warning("'obsm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- as_matrix(obsm)
      }
      k <- gsub(pattern = "^X_", replacement = "", x = py_to_r_auto(k))
      colnames(obsm) <- paste0(k, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateDimReducObject(embeddings = obsm, assay = "RNA", key = paste0(gsub(pattern = "_", replacement = "", x = k), "_"))
    }
  }

  if (inherits(adata$obsp, "python.builtin.object")) {
    keys <- iterate(adata$obsp$keys())
  } else {
    keys <- names(adata$obsp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsp <- tryCatch(py_to_r_auto(adata$obsp[[k]]), error = identity)
      if (inherits(obsp, "error")) {
        warning("'obsp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- as.sparse(obsp[1:nrow(obsp), , drop = FALSE])
      }
      colnames(obsp) <- py_to_r_auto(adata$obs_names$values)
      rownames(obsp) <- py_to_r_auto(adata$obs_names$values)
      obsp <- as.Graph(obsp[seq_len(nrow(obsp)), , drop = FALSE])
      DefaultAssay(object = obsp) <- "RNA"
      srt[[py_to_r_auto(k)]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- AddMetaData(srt[["RNA"]], metadata = as.data.frame(py_to_r_auto(adata$var)))
  }

  if (inherits(adata$varm, "python.builtin.object")) {
    keys <- iterate(adata$varm$keys())
  } else {
    keys <- names(adata$varm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varm <- tryCatch(py_to_r_auto(adata$varm[[k]]), error = identity)
      if (inherits(varm, "error")) {
        warning("'varm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- as_matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r_auto(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.loadings"]][[py_to_r_auto(k)]] <- varm
    }
  }

  if (inherits(adata$varp, "python.builtin.object")) {
    keys <- iterate(adata$varp$keys())
  } else {
    keys <- names(adata$varp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varp <- tryCatch(py_to_r_auto(adata$varp[[k]]), error = identity)
      if (inherits(varp, "error")) {
        warning("'varp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- as_matrix(varp)
      }
      colnames(varp) <- py_to_r_auto(adata$var_names$values)
      rownames(varp) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.graphs"]][[py_to_r_auto(k)]] <- varp
    }
  }

  if (inherits(adata$uns, "python.builtin.object")) {
    keys <- iterate(adata$uns$keys())
  } else {
    keys <- names(adata$uns)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      uns <- tryCatch(py_to_r_auto(adata$uns[[k]]), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r_auto(k)]] <- uns
      } else {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
    }
  }
  return(srt)
}

maxDepth <- function(x, depth = 0) {
  if (is.list(x)) {
    return(max(unlist(lapply(x, maxDepth, depth + 1))))
  } else {
    return(depth)
  }
}

#' @importFrom reticulate py_to_r
check_python_element <- function(x, depth = maxDepth(x)) {
  if (depth == 0 || !is.list(x) || !inherits(x, "python.builtin.object")) {
    if (inherits(x, "python.builtin.object")) {
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
      if (inherits(element, "python.builtin.object")) {
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
#' PAGA is a graph-based method used to infer cellular trajectories.
#' This function runs the PAGA analysis on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay_X Assay to convert as the main data matrix (X) in the anndata object.
#' @param slot_X Slot name for assay_X in the Seurat object.
#' @param assay_layers Assays to convert as layers in the anndata object.
#' @param slot_layers Slot names for the assay_layers in the Seurat object.
#' @param adata An anndata object.
#' @param group_by Variable to use for grouping cells in the Seurat object.
#' @param linear_reduction Linear reduction method to use, e.g., "PCA".
#' @param nonlinear_reduction Non-linear reduction method to use, e.g., "UMAP".
#' @param basis The basis to use for reduction, e.g., "UMAP".
#' @param n_pcs Number of principal components to use for linear reduction. Default is 30.
#' @param n_neighbors Number of neighbors to use for constructing the KNN graph. Default is 30.
#' @param use_rna_velocity Whether to use RNA velocity for PAGA analysis. Default is FALSE.
#' @param vkey The name of the RNA velocity data to use if \code{use_rna_velocity} is TRUE. Default is "stochastic".
#' @param embedded_with_PAGA Whether to embed data using PAGA layout. Default is FALSE.
#' @param paga_layout The layout for plotting PAGA graph. See See \href{https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.paga.html}{layout} param in scanpy.pl.paga function.
#' @param threshold The threshold for plotting PAGA graph. Edges for weights below this threshold will not draw.
#' @param point_size The point size for plotting.
#' @param infer_pseudotime Whether to infer pseudotime.
#' @param root_group The group to use as the root for pseudotime inference.
#' @param root_cell The cell to use as the root for pseudotime inference.
#' @param n_dcs TThe number of diffusion components to use for pseudotime inference.
#' @param n_branchings Number of branchings to detect.
#' @param min_group_size The minimum size of a group (as a fraction of the total number of cells) to consider it as a potential branching point.
#' @param palette The palette to use for coloring cells.
#' @param palcolor A vector of colors to use as the palette.
#' @param show_plot Whether to show the PAGA plot.
#' @param dpi The DPI (dots per inch) for saving the PAGA plot.
#' @param save Whether to save the PAGA plots.
#' @param dirpath The directory to save the PAGA plots.
#' @param fileprefix The file prefix to use for the PAGA plots.
#' @param return_seurat Whether to return a Seurat object instead of an anndata object. Default is TRUE.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{PAGAPlot}} \code{\link{CellDimPlot}} \code{\link{RunSCVELO}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPAGA(srt = pancreas_sub, assay_X = "RNA", group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "draw_graph_fr")
#' PAGAPlot(pancreas_sub, reduction = "UMAP")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", paga = pancreas_sub@misc$paga)
#'
#' pancreas_sub <- RunPAGA(
#'   srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP",
#'   embedded_with_PAGA = TRUE, infer_pseudotime = TRUE, root_group = "Ductal"
#' )
#' head(pancreas_sub[[]])
#' names(pancreas_sub@reductions)
#' FeatureDimPlot(pancreas_sub, features = "dpt_pseudotime", reduction = "PAGAUMAP2D")
#' PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "PAGAUMAP2D", paga = pancreas_sub@misc$paga)
#'
#' @export
RunPAGA <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                    adata = NULL, group_by = NULL,
                    linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                    n_pcs = 30, n_neighbors = 30, use_rna_velocity = FALSE, vkey = "stochastic",
                    embedded_with_PAGA = FALSE, paga_layout = "fr", threshold = 0.1, point_size = 20,
                    infer_pseudotime = FALSE, root_group = NULL, root_cell = NULL, n_dcs = 10, n_branchings = 0, min_group_size = 0.01,
                    palette = "Paired", palcolor = NULL,
                    show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                    return_seurat = !is.null(srt)) {
  check_Python("scanpy")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
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
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$PAGA, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(paga)|(distances)|(connectivities)|(draw_graph)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run scVelo workflow
#'
#' scVelo is a scalable toolkit for RNA velocity analysis in single cells. This function runs scVelo workflow on a Seurat object.
#'
#' @inheritParams RunPAGA
#' @param mode Velocity estimation model to use, e.g., "stochastic".
#' @param fitting_by Method used to fit gene velocities, e.g., "stochastic".
#' @param magic_impute Flag indicating whether to perform magic imputation.
#' @param knn The number of nearest neighbors for magic.MAGIC.
#' @param t power to which the diffusion operator is powered for magic.MAGIC.
#' @param min_shared_counts Minimum number of counts (both unspliced and spliced) required for a gene.
#' @param n_pcs Number of principal components (PCs) used for velocity estimation.
#' @param n_neighbors Number of nearest neighbors used for velocity estimation.
#' @param stream_smooth Multiplication factor for scale in Gaussian kernel around grid point.
#' @param stream_density Controls the closeness of streamlines. When density = 2 (default), the domain is divided into a 60x60 grid, whereas density linearly scales this grid. Each cell in the grid can have, at most, one traversing streamline.
#' @param arrow_length Length of arrows.
#' @param arrow_size Size of arrows.
#' @param arrow_density Amount of velocities to show.
#' @param denoise Boolean flag indicating whether to denoise.
#' @param denoise_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param kinetics Boolean flag indicating whether to estimate RNA kinetics.
#' @param kinetics_topn Number of genes with highest likelihood selected to infer velocity directions.
#' @param calculate_velocity_genes Boolean flag indicating whether to calculate velocity genes.
#' @param top_n The number of top features to plot.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}} \code{\link{VelocityPlot}} \code{\link{CellDimPlot}} \code{\link{RunPAGA}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, assay_X = "RNA", group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' head(pancreas_sub[[]])
#' names(pancreas_sub@assays)
#'
#' FeatureDimPlot(pancreas_sub, c("stochastic_length", "stochastic_confidence"))
#' FeatureDimPlot(pancreas_sub, "stochastic_pseudotime")
#' VelocityPlot(pancreas_sub, reduction = "UMAP", plot_type = "stream")
#' CellDimPlot(pancreas_sub, group.by = "SubCellType", reduction = "UMAP", pt.size = NA, velocity = "stochastic")
#'
#' pancreas_sub <- Standard_SCP(pancreas_sub, normalization_method = "SCT", nonlinear_reduction = "tsne")
#' pancreas_sub <- RunSCVELO(srt = pancreas_sub, assay_X = "SCT", group_by = "SubCellType", linear_reduction = "Standardpca", nonlinear_reduction = "StandardTSNE2D")
#'
#' @export
RunSCVELO <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                      adata = NULL, group_by = NULL,
                      linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                      mode = "stochastic", fitting_by = "stochastic",
                      magic_impute = FALSE, knn = 5, t = 2,
                      min_shared_counts = 30, n_pcs = 30, n_neighbors = 30,
                      stream_smooth = NULL, stream_density = 2,
                      arrow_length = 5, arrow_size = 5, arrow_density = 0.5,
                      denoise = FALSE, denoise_topn = 3, kinetics = FALSE, kinetics_topn = 100,
                      calculate_velocity_genes = FALSE, top_n = 6, n_jobs = 1,
                      palette = "Paired", palcolor = NULL,
                      show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                      return_seurat = !is.null(srt)) {
  check_Python("scvelo")
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
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
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$SCVELO, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = paste0("(velocity)|(distances)|(connectivities)|(Ms)|(Mu)|(", paste(mode, collapse = ")|("), ")|(paga)"), overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run Palantir analysis
#'
#' @inheritParams RunPAGA
#' @param dm_n_components The number of diffusion components to calculate.
#' @param dm_alpha Normalization parameter for the diffusion operator.
#' @param dm_n_eigs Number of eigen vectors to use.
#' @param early_group Name of the group to start Palantir analysis from.
#' @param early_cell Name of the cell to start Palantir analysis from.
#' @param terminal_groups Character vector specifying terminal groups for Palantir analysis.
#' @param terminal_cells Character vector specifying terminal cells for Palantir analysis.
#' @param num_waypoints Number of waypoints to be included.
#' @param scale_components Should the cell fate probabilities be scaled for each component independently?
#' @param use_early_cell_as_start Should the starting cell for each terminal group be set as early_cell?
#' @param adjust_early_cell Whether to adjust the early cell to the cell with the minimum pseudotime value.
#' @param adjust_terminal_cells hether to adjust the terminal cells to the cells with the maximum pseudotime value for each terminal group.
#' @param max_iterations Maximum number of iterations for pseudotime convergence.
#' @param n_jobs The number of parallel jobs to run.
#'
#' @seealso \code{\link{srt_to_adata}}
#'
#' @examples
#' data("pancreas_sub")
#' pancreas_sub <- RunPalantir(
#'   srt = pancreas_sub, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP",
#'   early_group = "Ductal", use_early_cell_as_start = TRUE,
#'   terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
#' )
#' head(pancreas_sub[[]])
#' FeatureDimPlot(pancreas_sub, c("palantir_pseudotime", "palantir_diff_potential"))
#' FeatureDimPlot(pancreas_sub, paste0(c("Alpha", "Beta", "Delta", "Epsilon"), "_diff_potential"))
#'
#' @export
RunPalantir <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, group_by = NULL,
                        linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                        n_pcs = 30, n_neighbors = 30, dm_n_components = 10, dm_alpha = 0, dm_n_eigs = NULL,
                        early_group = NULL, early_cell = NULL, terminal_cells = NULL, terminal_groups = NULL,
                        num_waypoints = 1200, scale_components = TRUE, use_early_cell_as_start = TRUE,
                        adjust_early_cell = FALSE, adjust_terminal_cells = FALSE,
                        max_iterations = 25, n_jobs = 8,
                        point_size = 20, palette = "Paired", palcolor = NULL,
                        show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                        return_seurat = !is.null(srt)) {
  check_Python("palantir")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by) && any(!is.null(early_group), !is.null(terminal_groups))) {
    stop("'group_by' must be provided when early_group or terminal_groups provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
  }
  if (is.null(early_cell) && is.null(early_group)) {
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
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$Palantir, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(palantir)|(dm_kernel)|(_diff_potential)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

#' Run WOT analysis
#'
#' @inheritParams RunPAGA
#' @param time_field A character string specifying the column name in `adata.obs` or `srt@meta.data` that contains the time information.
#' @param growth_iters An integer specifying the number of growth iterations to perform during the OT Model computation. Default is 3.
#' @param tmap_out A character string specifying the path to store the computed transport maps.
#' @param time_from A numeric value specifying the starting time point for trajectory and fate analysis.
#' @param time_to A numeric value specifying the ending time point for trajectory and fate analysis. If not provided, only trajectory and fate analysis for the specified `time_from` will be performed.
#' @param get_coupling A logical value indicating whether to compute and store the coupling matrix between the specified `time_from` and `time_to`. Default is FALSE.
#' @param recalculate A logical value indicating whether to recalculate the transport maps even if they already exist at the specified `tmap_out` location. Default is FALSE.
#'
#' @seealso \code{\link{srt_to_adata}}
#' @export
RunWOT <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                   adata = NULL, group_by = NULL,
                   time_field = "Time", growth_iters = 3L, tmap_out = "tmaps/tmap_out",
                   time_from = NULL, time_to = NULL, get_coupling = FALSE, recalculate = FALSE,
                   palette = "Paired", palcolor = NULL,
                   show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                   return_seurat = !is.null(srt)) {
  check_Python("wot")
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(time_field)) {
    stop("'time_field' must be provided.")
  }
  if (is.null(time_from)) {
    stop("'time_from' must be provided.")
  }
  if (isTRUE(get_coupling) && is.null(time_to)) {
    warning("The 'get_coupling' paramter is only valid when 'time_to' is specified.")
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
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$WOT, args)

  if (isTRUE(return_seurat)) {
    srt_out <- adata_to_srt(adata)
    if (is.null(srt)) {
      return(srt_out)
    } else {
      srt_out1 <- SrtAppend(srt_raw = srt, srt_append = srt_out)
      srt_out2 <- SrtAppend(srt_raw = srt_out1, srt_append = srt_out, pattern = "(trajectory_)|(fates_)|(transition_)|(coupling_)", overwrite = TRUE, verbose = FALSE)
      return(srt_out2)
    }
  } else {
    return(adata)
  }
}

RunCellRank <- function(srt = NULL, assay_X = "RNA", slot_X = "counts", assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                        adata = NULL, group_by = NULL, n_jobs = 1,
                        linear_reduction = NULL, nonlinear_reduction = NULL, basis = NULL,
                        mode = "stochastic", fitting_by = "stochastic",
                        magic_impute = FALSE, knn = 5, t = 2,
                        min_shared_counts = 30, n_pcs = 30, n_neighbors = 30, approx = TRUE,
                        stream_smooth = NULL, stream_density = 2,
                        arrow_size = 5, arrow_length = 5, arrow_density = 0.5,
                        s_genes = NULL, g2m_genes = NULL, calculate_velocity_genes = FALSE,
                        denoise = FALSE, kinetics = FALSE, axis = "equal",
                        show_plot = TRUE, save = FALSE, dpi = 300, dirpath = "./", fileprefix = "",
                        return_seurat = !is.null(srt)) {
  check_Python("cellrank")
  if (isTRUE(magic_impute)) {
    check_Python("magic-impute")
  }
  if (all(is.null(srt), is.null(adata))) {
    stop("One of 'srt', 'adata' must be provided.")
  }
  if (is.null(group_by)) {
    stop("'group_by' must be provided.")
  }
  if (is.null(linear_reduction) && is.null(nonlinear_reduction)) {
    stop("'linear_reduction' or 'nonlinear_reduction' must be provided at least one.")
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
  args <- args[!names(args) %in% c("srt", "assay_X", "slot_X", "assay_layers", "slot_layers", "return_seurat", "palette", "palcolor")]

  if (!is.null(srt)) {
    args[["adata"]] <- srt_to_adata(
      srt = srt,
      assay_X = assay_X, slot_X = slot_X,
      assay_layers = assay_layers, slot_layers = slot_layers
    )
  }
  groups <- py_to_r_auto(args[["adata"]]$obs)[[group_by]]
  args[["palette"]] <- palette_scp(levels(groups) %||% unique(groups), palette = palette, palcolor = palcolor)

  SCP_analysis <- reticulate::import_from_path("SCP_analysis", path = system.file("python", package = "SCP", mustWork = TRUE), convert = TRUE)
  adata <- do.call(SCP_analysis$CellRank, args)

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
