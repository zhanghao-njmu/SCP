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
"pancreas_sub"

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
"panc8_sub"

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
"ifnb_sub"

# ref_scHCL <- scHCL::ref.expr
# ref_scMCA <- scMCA::ref.expr
# ref_scZCL <- scZCL::ref.expr
# colnames(ref_scHCL) <- stringi::stri_trans_general(colnames(ref_scHCL), "latin-ascii")
# colnames(ref_scMCA) <- stringi::stri_trans_general(colnames(ref_scMCA), "latin-ascii")
# colnames(ref_scZCL) <- stringi::stri_trans_general(colnames(ref_scZCL), "latin-ascii")
# ref_scHCL <- as.matrix(log1p(ref_scHCL))
# ref_scMCA <- as.matrix(log1p(ref_scMCA))
# ref_scZCL <- as.matrix(log1p(ref_scZCL))
# save(ref_scHCL, file = "data/ref_scHCL.rda",version = 2,compress = "xz")
# save(ref_scMCA, file = "data/ref_scMCA.rda",version = 2,compress = "xz")
# save(ref_scZCL, file = "data/ref_scZCL.rda",version = 2,compress = "xz")

# load("~/Git/pancreas.rda")
# pancreas_sub <- subset(pancreas, cells = sample(colnames(pancreas),size = 1000));
# CellDimPlot(pancreas_sub,reduction = "umap",group.by = "CellType")
# save(pancreas_sub, file = "data/pancreas1k.rda",version = 2,compress = "xz")

# lifemap_cell <- readRDS("data/cell_gene_list.rds") %>% bind_rows()
# lifemap_compartment <- readRDS("data/compartment_gene_list.rds") %>% bind_rows()
# lifemap_organ <- readRDS("data/organ_gene_list.rds") %>% bind_rows()
# save(lifemap_cell, file = "data/lifemap_cell.rda",version = 2,compress = "xz",compression_level = 9)
# save(lifemap_compartment, file = "data/lifemap_compartment.rda",version = 2,compress = "xz",compression_level = 9)
# save(lifemap_organ, file = "data/lifemap_organ.rda",version = 2,compress = "xz",compression_level = 9)
NULL
