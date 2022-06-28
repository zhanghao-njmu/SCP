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

# lifemap_cell <- readRDS("/home/zhanghao/Git/SCP/data/cell_gene_list.rds") %>% bind_rows()
# lifemap_compartment <- readRDS("/home/zhanghao/Git/SCP/data/compartment_gene_list.rds") %>% bind_rows()
# lifemap_organ <- readRDS("/home/zhanghao/Git/SCP/data/organ_gene_list.rds") %>% bind_rows()
# save(lifemap_cell, file = "~/Git/SCP/data/lifemap_cell.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)
# save(lifemap_compartment, file = "~/Git/SCP/data/lifemap_compartment.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)
# save(lifemap_organ, file = "~/Git/SCP/data/lifemap_organ.rda",version = 2,compress = "xz",compression_level = 9,ascii = TRUE)
NULL
