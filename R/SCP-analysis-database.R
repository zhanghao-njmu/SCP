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
