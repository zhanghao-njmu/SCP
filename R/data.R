#' A subsetted version of mouse 'pancreas' datasets
#'
#' Mouse pancreatic endocrinogenesis dataset from \href{https://doi.org/10.1242/dev.173849}{Bastidas-Ponce et al. (2019)}. A total of 1000 cells were downsampled to form the \code{pancreas_sub} dataset.
#'
#' @format A \code{Seurat} object.
#' @concept data
#' @source \url{https://scvelo.readthedocs.io/scvelo.datasets.pancreas/} \url{https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(Seurat)
#'   library(reticulate)
#'   check_Python("scvelo")
#'   scv <- import("scvelo")
#'   adata <- scv$datasets$pancreas()
#'   pancreas <- adata_to_srt(adata)
#'   set.seed(11)
#'   pancreas_sub <- subset(pancreas, cells = sample(colnames(pancreas), size = 1000))
#'   pancreas_sub <- pancreas_sub[rowSums(pancreas_sub@assays$RNA@counts) > 0, ]
#'   pancreas_sub[["CellType"]] <- pancreas_sub[["clusters_coarse"]]
#'   pancreas_sub[["SubCellType"]] <- pancreas_sub[["clusters"]]
#'   pancreas_sub[["clusters_coarse"]] <- pancreas_sub[["clusters"]] <- NULL
#'   pancreas_sub[["Phase"]] <- ifelse(pancreas_sub$S_score > pancreas_sub$G2M_score, "S", "G2M")
#'   pancreas_sub[["Phase"]][apply(pancreas_sub[[]][, c("S_score", "G2M_score")], 1, max) < 0, ] <- "G1"
#'   pancreas_sub[["Phase", drop = TRUE]] <- factor(pancreas_sub[["Phase", drop = TRUE]], levels = c("G1", "S", "G2M"))
#'   pancreas_sub[["PCA"]] <- pancreas_sub[["X_pca"]]
#'   pancreas_sub[["UMAP"]] <- pancreas_sub[["X_umap"]]
#'   pancreas_sub[["X_umap"]] <- pancreas_sub[["X_pca"]] <- NULL
#'   VariableFeatures(pancreas_sub) <- rownames(pancreas_sub[["RNA"]])[which(pancreas_sub[["RNA"]]@meta.features$highly_variable_genes == "True")]
#'   # usethis::use_data(pancreas_sub, compress = "xz")
#' }
#' }
#' @name pancreas_sub
NULL

#' A subsetted version of human 'panc8' datasets
#'
#' Human pancreatic islet cell datasets produced across four technologies, CelSeq (GSE81076) CelSeq2 (GSE85241), Fluidigm C1 (GSE86469), and SMART-Seq2 (E-MTAB-5061) from \href{https://github.com/satijalab/seurat-data}{SeuratData} package. For each data set in `panc8`, 200 cells were downsampled to form the \code{panc8_sub} dataset.
#'
#' @format A \code{Seurat} object.
#' @concept data
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076} \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi} \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/} \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi} \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   data("pancreas_sub")
#'   if (!require("SeuratData", quietly = TRUE)) {
#'     devtools::install_github("satijalab/seurat-data")
#'   }
#'   library(SeuratData)
#'   library(Seurat)
#'   suppressWarnings(InstallData("panc8"))
#'   data("panc8")
#'   set.seed(11)
#'   cells_sub <- unlist(lapply(split(colnames(panc8), panc8$dataset), function(x) sample(x, size = 200)))
#'   panc8_sub <- subset(panc8, cells = cells_sub)
#'   panc8_sub <- panc8_sub[rowSums(panc8_sub@assays$RNA@counts) > 0, ]
#'   panc8_sub <- panc8_sub[toupper(rownames(panc8_sub)) %in% toupper(rownames(pancreas_sub)), ]
#'   panc8_sub <- UpdateSeuratObject(panc8_sub)
#'   # usethis::use_data(panc8_sub, compress = "xz")
#' }
#' }
#' @name panc8_sub
NULL

#' A subsetted version of 'ifnb' datasets
#'
#' Human PBMC control/IFNB-stimulated dataset
#'
#' @format A \code{Seurat} object.
#' @concept data
#' @source \url{https://www.nature.com/articles/nbt.4042}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   if (!require("SeuratData", quietly = TRUE)) {
#'     devtools::install_github("satijalab/seurat-data")
#'   }
#'   library(SeuratData)
#'   library(Seurat)
#'   suppressWarnings(InstallData("ifnb"))
#'   data("ifnb")
#'   set.seed(11)
#'   cells_sub <- unlist(lapply(split(colnames(ifnb), ifnb$stim), function(x) sample(x, size = 1000)))
#'   ifnb_sub <- subset(ifnb, cells = cells_sub)
#'   ifnb_sub <- ifnb_sub[rowSums(ifnb_sub@assays$RNA@counts) > 0, ]
#'   ifnb_sub <- UpdateSeuratObject(ifnb_sub)
#'   # usethis::use_data(ifnb_sub, compress = "xz")
#' }
#' }
#' @name ifnb_sub
NULL

#' Excluded words in keyword enrichment analysis and extraction
#'
#' The variable "words_excluded" represents the words that are excluded during keyword enrichment analysis or keyword extraction process. These mainly include words that are excessively redundant or of little value.
#' @concept data
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   words_excluded <- c(
#'     "the", "is", "and", "or", "a", "in", "on", "under", "between", "of", "through", "via", "along", "that", "for", "with", "within", "without",
#'     "cell", "cellular", "dna", "rna", "protein", "peptide", "amino", "acid", "development", "involved", "organization", "system",
#'     "regulation", "regulated", "positive", "negative", "response", "process", "processing", "small", "large", "change"
#'   )
#'   # usethis::use_data(words_excluded)
#' }
#' }
#' @name words_excluded
NULL

#' A list of palettes for use in data visualization
#'
#' @concept data
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_R(c("stringr", "RColorBrewer", "ggsci", "Redmonder", "rcartocolor", "nord", "viridis", "pals", "oompaBase", "dichromat", "jaredhuling/jcolors"))
#'   library(stringr)
#'   library(RColorBrewer)
#'   library(ggsci)
#'   library(Redmonder)
#'   library(rcartocolor)
#'   library(nord)
#'   library(viridis)
#'   library(pals)
#'   library(dichromat)
#'   library(jcolors)
#'   brewer.pal.info <- RColorBrewer::brewer.pal.info
#'   ggsci_db <- ggsci:::ggsci_db
#'   redmonder.pal.info <- Redmonder::redmonder.pal.info
#'   metacartocolors <- rcartocolor::metacartocolors
#'   rownames(metacartocolors) <- metacartocolors$Name
#'   nord_palettes <- nord::nord_palettes
#'   viridis_names <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
#'   viridis_palettes <- lapply(setNames(viridis_names, viridis_names), function(x) viridis::viridis(100, option = x))
#'   ocean_names <- names(pals:::syspals)[grep("ocean", names(pals:::syspals))]
#'   ocean_palettes <- pals:::syspals[ocean_names]
#'   dichromat_palettes <- dichromat::colorschemes
#'   jcolors_names <- paste0("jcolors-", c("default", "pal2", "pal3", "pal4", "pal5", "pal6", "pal7", "pal8", "pal9", "pal10", "pal11", "pal12", "rainbow"))
#'   custom_names <- c("jet", "simspec", "GdRd")
#'   custom_palettes <- list(
#'     oompaBase::jetColors(N = 100),
#'     c("#c22b86", "#f769a1", "#fcc5c1", "#253777", "#1d92c0", "#9ec9e1", "#015b33", "#42aa5e", "#d9f0a2", "#E66F00", "#f18c28", "#FFBB61"),
#'     c("gold", "red3")
#'   )
#'   names(custom_palettes) <- custom_names
#'
#'   palette_list <- list()
#'   all_colors <- c(
#'     rownames(brewer.pal.info), names(ggsci_db), rownames(redmonder.pal.info),
#'     rownames(metacartocolors), names(nord_palettes), names(viridis_palettes),
#'     ocean_names, names(dichromat_palettes), jcolors_names,
#'     custom_names
#'   )
#'   for (pal in all_colors) {
#'     if (!pal %in% all_colors) {
#'       stop(paste0("Invalid pal Must be one of ", paste0(all_colors, collapse = ",")))
#'     }
#'     if (pal %in% rownames(brewer.pal.info)) {
#'       pal_n <- brewer.pal.info[pal, "maxcolors"]
#'       pal_category <- brewer.pal.info[pal, "category"]
#'       if (pal_category == "div") {
#'         palcolor <- rev(brewer.pal(name = pal, n = pal_n))
#'       } else {
#'         if (pal == "Paired") {
#'           palcolor <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
#'         } else {
#'           palcolor <- brewer.pal(name = pal, n = pal_n)
#'         }
#'       }
#'       if (pal_category == "qual") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% names(ggsci_db)) {
#'       if (pal %in% c("d3", "uchicago", "material")) {
#'         for (subpal in names(ggsci_db[[pal]])) {
#'           palcolor <- ggsci_db[[pal]][[subpal]]
#'           if (pal == "material") {
#'             attr(palcolor, "type") <- "continuous"
#'           } else {
#'             attr(palcolor, "type") <- "discrete"
#'           }
#'           palette_list[[paste0(pal, "-", subpal)]] <- palcolor
#'         }
#'         next
#'       } else {
#'         palcolor <- ggsci_db[[pal]][[1]]
#'         if (pal == "gsea") {
#'           attr(palcolor, "type") <- "continuous"
#'         } else {
#'           attr(palcolor, "type") <- "discrete"
#'         }
#'       }
#'     } else if (pal %in% rownames(redmonder.pal.info)) {
#'       pal_n <- redmonder.pal.info[pal, "maxcolors"]
#'       pal_category <- redmonder.pal.info[pal, "category"]
#'       if (pal_category == "div") {
#'         palcolor <- rev(redmonder.pal(name = pal, n = pal_n))
#'       } else {
#'         palcolor <- redmonder.pal(name = pal, n = pal_n)
#'       }
#'       if (pal_category == "qual") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% rownames(metacartocolors)) {
#'       pal_n <- metacartocolors[pal, "Max_n"]
#'       palcolor <- carto_pal(name = pal, n = pal_n)
#'       if (pal_category == "qualitative") {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% names(nord_palettes)) {
#'       palcolor <- nord_palettes[[pal]]
#'       attr(palcolor, "type") <- "discrete"
#'     } else if (pal %in% names(viridis_palettes)) {
#'       palcolor <- viridis_palettes[[pal]]
#'       attr(palcolor, "type") <- "continuous"
#'     } else if (pal %in% names(ocean_palettes)) {
#'       palcolor <- ocean_palettes[[pal]]
#'       attr(palcolor, "type") <- "continuous"
#'     } else if (pal %in% names(dichromat_palettes)) {
#'       palcolor <- dichromat_palettes[[pal]]
#'       if (pal %in% c("Categorical.12", "SteppedSequential.5")) {
#'         attr(palcolor, "type") <- "discrete"
#'       } else {
#'         attr(palcolor, "type") <- "continuous"
#'       }
#'     } else if (pal %in% jcolors_names) {
#'       palcolor <- jcolors(palette = gsub("jcolors-", "", pal))
#'       if (pal %in% paste0("jcolors-", c("pal10", "pal11", "pal12", "rainbow"))) {
#'         attr(palcolor, "type") <- "continuous"
#'       } else {
#'         attr(palcolor, "type") <- "discrete"
#'       }
#'     } else if (pal %in% custom_names) {
#'       palcolor <- custom_palettes[[pal]]
#'       if (pal %in% c("jet")) {
#'         attr(palcolor, "type") <- "continuous"
#'       } else {
#'         attr(palcolor, "type") <- "discrete"
#'       }
#'     }
#'     palette_list[[pal]] <- palcolor
#'   }
#'   # usethis::use_data(palette_list)
#' }
#' }
#' @name palette_list
NULL

#' Reference datasets for cell type annotation in single-cell RNA data
#'
#' @concept data
#' @source \url{https://github.com/ggjlab/scHCL} \url{https://github.com/ggjlab/scMCA} \url{https://github.com/ggjlab/scZCL}
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_R(c("ggjlab/scZCL", "ggjlab/scHCL", "ggjlab/scMCA"))
#'   ref_scHCL <- Seurat::NormalizeData(scHCL::ref.expr)
#'   ref_scMCA <- Seurat::NormalizeData(scMCA::ref.expr)
#'   ref_scZCL <- Seurat::NormalizeData(scZCL::ref.expr)
#'   Encoding(colnames(ref_scHCL)) <- "latin1"
#'   colnames(ref_scHCL) <- iconv(colnames(ref_scHCL), "latin1", "UTF-8")
#'   Encoding(colnames(ref_scMCA)) <- "latin1"
#'   colnames(ref_scMCA) <- iconv(colnames(ref_scMCA), "latin1", "UTF-8")
#'   Encoding(colnames(ref_scZCL)) <- "latin1"
#'   colnames(ref_scZCL) <- iconv(colnames(ref_scZCL), "latin1", "UTF-8")
#'   # usethis::use_data(ref_scHCL, compress = "xz")
#'   # usethis::use_data(ref_scMCA, compress = "xz")
#'   # usethis::use_data(ref_scZCL, compress = "xz")
#' }
#' }
#' @usage ref_scHCL
#' @usage ref_scMCA
#' @usage ref_scZCL
#' @name ref_scHCL
#' @aliases ref_scHCL ref_scMCA ref_scZCL
NULL

#' Embryonic Development Database from LifeMap Discovery
#' @concept data
#' @source \url{https://discovery.lifemapsc.com/in-vivo-development/cellular}
#' @usage lifemap_cell
#' @usage lifemap_compartment
#' @usage lifemap_organ
#' @name lifemap_cell
#' @aliases lifemap_cell lifemap_compartment lifemap_organ
NULL
