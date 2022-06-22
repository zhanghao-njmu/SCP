# ref_scHCL <- scHCL::ref.expr
# ref_scMCA <- scMCA::ref.expr
# ref_scZCL <- scZCL::ref.expr
# ### textclean::replace_non_ascii(colnames(ref_scHCL))
# colnames(ref_scHCL) <- stringi::stri_trans_general(colnames(ref_scHCL), "latin-ascii")
# colnames(ref_scMCA) <- stringi::stri_trans_general(colnames(ref_scMCA), "latin-ascii")
# colnames(ref_scZCL) <- stringi::stri_trans_general(colnames(ref_scZCL), "latin-ascii")
# ref_scHCL <- as.matrix(log1p(ref_scHCL))
# ref_scMCA <- as.matrix(log1p(ref_scMCA))
# ref_scZCL <- as.matrix(log1p(ref_scZCL))
# save(ref_scHCL, file = "~/Git/SCP/data/ref_scHCL.rda",version = 2,compress = "xz",ascii = TRUE)
# save(ref_scMCA, file = "~/Git/SCP/data/ref_scMCA.rda",version = 2,compress = "xz",ascii = TRUE)
# save(ref_scZCL, file = "~/Git/SCP/data/ref_scZCL.rda",version = 2,compress = "xz",ascii = TRUE)

# load("~/Git/pancreas.rda")
# pancreas_sub <- subset(pancreas, cells = sample(colnames(pancreas),size = 1000));
### ClassDimPlot(pancreas_sub,reduction = "umap",group.by = "CellType")
### save(pancreas_sub, file = "~/Git/SCP/data/pancreas1k.rda",version = 2,compress = "xz",ascii = TRUE)



# library(stringr)
# library(RColorBrewer)
# library(ggsci)
# library(Redmonder)
# library(rcartocolor)
# library(nord)
# library(viridis)
# library(pals)
# brewer.pal.info <- RColorBrewer::brewer.pal.info
# ggsci_db <- ggsci:::ggsci_db
# redmonder.pal.info <- Redmonder::redmonder.pal.info
# metacartocolors <- rcartocolor::metacartocolors
# rownames(metacartocolors) <- metacartocolors$Name
# nord_palettes <- nord::nord_palettes
# viridis_names <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
# viridis_palettes <- lapply(setNames(viridis_names, viridis_names), function(x) viridis::viridis(100, option = x))
# ocean_names <- names(pals:::syspals)[str_detect(names(pals:::syspals), "ocean")]
# ocean_palettes <- pals:::syspals[ocean_names]
# custom_names <- c("jet","simspec")
# custom_palettes <- list(oompaBase::jetColors(N = 100),
#                      c("#c22b86","#f769a1","#fcc5c1","#253777","#1d92c0","#9ec9e1","#015b33","#42aa5e","#d9f0a2","#E66F00","#f18c28","#FFBB61"))
# names(custom_palettes) <- custom_names
#
# palette_list <- list()
# all_colors <- c(rownames(brewer.pal.info), names(ggsci_db), rownames(redmonder.pal.info),
#                 rownames(metacartocolors), names(nord_palettes), names(viridis_palettes),
#                 ocean_names,custom_names)
# for (pal in all_colors) {
#   if (!pal %in% all_colors) {
#     stop(paste0("Invalid pal Must be one of ", paste0(all_colors, collapse = ",")))
#   }
#   if (pal %in% rownames(brewer.pal.info)) {
#     pal_n <- brewer.pal.info[pal, "maxcolors"]
#     pal_category <- brewer.pal.info[pal, "category"]
#     if (pal_category == "div") {
#       pal_color <- rev(brewer.pal(name = pal, n = pal_n))
#     } else {
#       if (pal == "Paired") {
#         pal_color <- brewer.pal(12, "Paired")[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
#       } else {
#         pal_color <- brewer.pal(name = pal, n = pal_n)
#       }
#     }
#   } else if (pal %in% names(ggsci_db)) {
#     if (pal %in% c("d3", "uchicago", "material")) {
#       for (subpal in names(ggsci_db[[pal]])) {
#         pal_color <- ggsci_db[[pal]][[subpal]]
#         pal_n <- length(pal_color)
#         palette_list[[paste0(pal, "-", subpal)]][["pal_color"]] <- pal_color
#         palette_list[[paste0(pal, "-", subpal)]][["pal_n"]] <- pal_n
#       }
#       next
#     } else {
#       pal_color <- ggsci_db[[pal]][[1]]
#       pal_n <- length(pal_color)
#     }
#   } else if (pal %in% rownames(redmonder.pal.info)) {
#     pal_n <- redmonder.pal.info[pal, "maxcolors"]
#     pal_category <- redmonder.pal.info[pal, "category"]
#     if (pal_category == "div") {
#       pal_color <- rev(redmonder.pal(name = pal, n = pal_n))
#     } else {
#       pal_color <- redmonder.pal(name = pal, n = pal_n)
#     }
#   } else if (pal %in% rownames(metacartocolors)) {
#     pal_n <- metacartocolors[pal, "Max_n"]
#     pal_color <- carto_pal(name = pal, n = pal_n)
#   } else if (pal %in% names(nord_palettes)) {
#     pal_color <- nord_palettes[[pal]]
#     pal_n <- length(pal_color)
#   } else if (pal %in% names(viridis_palettes)) {
#     pal_color <- viridis_palettes[[pal]]
#     pal_n <- length(pal_color)
#   }else if (pal %in% names(ocean_palettes)) {
#     pal_color <- ocean_palettes[[pal]]
#     pal_n <- length(pal_color)
#   }else if (pal %in% custom_names) {
#     pal_color <- custom_palettes[[pal]]
#     pal_n <- length(pal_color)
#   }
#   palette_list[[pal]][["pal_color"]] <- pal_color
#   palette_list[[pal]][["pal_n"]] <- pal_n
# }
# save(palette_list, file = "~/Git/SCP/data/palette_list.rda", version = 2)

# lifemap_cell <- readRDS("/home/zhanghao/Git/SCP/data/cell_gene_list.rds") %>% bind_rows()
# lifemap_compartment <- readRDS("/home/zhanghao/Git/SCP/data/compartment_gene_list.rds") %>% bind_rows()
# lifemap_organ <- readRDS("/home/zhanghao/Git/SCP/data/organ_gene_list.rds") %>% bind_rows()
# save(lifemap_cell, file = "~/Git/SCP/data/lifemap_cell.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)
# save(lifemap_compartment, file = "~/Git/SCP/data/lifemap_compartment.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)
# save(lifemap_organ, file = "~/Git/SCP/data/lifemap_organ.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)

NULL
