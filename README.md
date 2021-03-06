
# SCP: Single Cell Pipeline

<!-- badges: start -->

[![version](https://img.shields.io/github/r-package/v/zhanghao-njmu/SCP)](https://github.com/zhanghao-njmu/SCP)
[![R-CMD-check](https://github.com/zhanghao-njmu/SCP/workflows/R-CMD-check/badge.svg)](https://github.com/zhanghao-njmu/SCP/actions)
[![codecov](https://codecov.io/gh/zhanghao-njmu/SCP/branch/main/graph/badge.svg)](https://codecov.io/gh/zhanghao-njmu/SCP?branch=main)
[![codesize](https://img.shields.io/github/languages/code-size/zhanghao-njmu/SCP.svg)](https://github.com/zhanghao-njmu/SCP)
<!-- badges: end -->

The SCP package provides a comprehensive set of tools for single cell
data processing and downstream analysis.

The package includes facilities for:

-   Integrated single cell quality control methods.
-   One-stop single-cell pipeline embedded with multiple methods for
    normalization, feature reduction, and cell population
    identification.
-   One-stop single-cell pipeline embedded with multiple data
    integration methods.
-   Multiple single cell downstream analyses such as identification of
    differential features, enrichment analysis, GSEA analysis,
    identification of dynamic features,
    [PAGA](https://github.com/theislab/paga), [RNA
    velocity](https://github.com/theislab/scvelo), etc.
-   Multiple methods for automatic annotation of single-cell data and
    methods for projection between single-cell datasets.
-   Publication-quality visualization of multiple analysis results.
-   Fast deployment of single-cell data to a [shiny
    app](https://shiny.rstudio.com/).

The functions in the SCP package are all developed around the [Seurat
object](https://github.com/mojaveazure/seurat-object) and compatible
with other Seurat functions.

## Installation

You can install the development version of SCP from
[GitHub](https://github.com/zhanghao-njmu/SCP) with:

``` r
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("zhanghao-njmu/SCP")
```

## Example

### Load the Data

The analysis is based on a subsetted version of [mouse pancreas
data](https://doi.org/10.1242/dev.173849).

``` r
library(SCP)
data("pancreas1k")
ClassDimPlot(pancreas1k, group.by = c("CellType", "SubCellType"), reduction = "UMAP", theme_use = "theme_blank")
```

<img src="man/figures/README-library-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpDimPlot(pancreas1k,
  features = c("Sox9", "Neurog3", "Fev", "Rbp4"),
  reduction = "UMAP", theme_use = "theme_blank"
)
```

<img src="man/figures/README-library-2.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpDotPlot(pancreas1k,
  features = c(
    "Sox9", "Anxa2", "Bicc1", # Ductal
    "Neurog3", "Hes6", # EPs
    "Fev", "Neurod1", # Pre-endocrine
    "Rbp4", "Pyy", # Endocrine
    "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
  ),
  cell_split_by = c("CellType", "SubCellType")
)
```

<img src="man/figures/README-library-3.png" width="100%" style="display: block; margin: auto;" />

### CellQC

``` r
pancreas1k <- RunCellQC(pancreas1k)
ClassDimPlot(pancreas1k, group.by = "CellQC", reduction = "UMAP")
```

<img src="man/figures/README-RunCellQC-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassStatPlot(pancreas1k, stat.by = "CellQC", group.by = "CellType", label = TRUE)
```

<img src="man/figures/README-RunCellQC-2.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassStatPlot(pancreas1k, stat.by = c("db_qc", "outlier_qc", "umi_qc", "gene_qc",
    "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc"), plot_type = "upset",
    stat_level = "Fail")
```

<img src="man/figures/README-RunCellQC-3.png" width="100%" style="display: block; margin: auto;" />

### Standard SCP

``` r
pancreas1k <- Standard_SCP(pancreas1k)
ClassDimPlot(pancreas1k, group.by = c("CellType", "SubCellType"), reduction = "StandardUMAP2D",
    theme_use = "theme_blank")
```

<img src="man/figures/README-Standard_SCP-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassDimPlot3D(pancreas1k, group.by = "SubCellType")
```

### Integration SCP

Example data for integration is [panc8(eight human pancreas
datasets)](https://github.com/satijalab/seurat-data)

``` r
if (!require("SeuratData", quietly = TRUE)) {
    devtools::install_github("satijalab/seurat-data")
}
library(SeuratData)
suppressWarnings(InstallData("panc8"))
data("panc8")
cell_sub <- unlist(lapply(split(colnames(panc8), panc8$tech), function(x) sample(x,
    size = 500)))
panc8 <- subset(panc8, cells = cell_sub)
panc8 <- Integration_SCP(panc8, batch = "tech", integration_method = "Seurat")
ClassDimPlot(panc8, group.by = c("celltype", "tech"), reduction = "SeuratUMAP2D",
    theme_use = "theme_blank")
```

<img src="man/figures/README-Integration_SCP-1.png" width="100%" style="display: block; margin: auto;" />

### Cell projection between single-cell datasets

``` r
library(stringr)
panc8 <- RenameFeatures(panc8, newnames = make.unique(str_to_title(rownames(panc8))))
pancreas1k <- RunKNNMap(srt_query = pancreas1k, srt_ref = panc8, ref_umap = "SeuratUMAP2D")
ProjectionPlot(srt_query = pancreas1k, srt_ref = panc8, query_group = "SubCellType",
    ref_group = "celltype")
```

<img src="man/figures/README-RunKNNMap-1.png" width="100%" style="display: block; margin: auto;" />

### Cell annotation using bulk RNA-seq datasets

``` r
data("ref_scMCA")
pancreas1k <- RunKNNPredict(srt_query = pancreas1k, bulk_ref = ref_scMCA, filter_lowfreq = 20)
ClassDimPlot(pancreas1k, group.by = "knnpredict_classification", reduction = "UMAP",
    label = TRUE)
```

<img src="man/figures/README-RunKNNPredict-bulk-1.png" width="100%" style="display: block; margin: auto;" />

### Cell annotation using single-cell datasets

``` r
pancreas1k <- RunKNNPredict(srt_query = pancreas1k, srt_ref = panc8, ref_group = "celltype",
    filter_lowfreq = 20)
ClassDimPlot(pancreas1k, group.by = "knnpredict_classification", reduction = "UMAP",
    label = TRUE)
```

<img src="man/figures/README-RunKNNPredict-scrna-1.png" width="100%" style="display: block; margin: auto;" />

### PAGA analysis

``` r
pancreas1k <- RunPAGA(srt = pancreas1k, group_by = "SubCellType", liner_reduction = "PCA",
    nonliner_reduction = "UMAP", return_seurat = TRUE)
PAGAPlot(pancreas1k, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE)
```

<img src="man/figures/README-RunPAGA-1.png" width="100%" style="display: block; margin: auto;" />

### Velocity analysis

``` r
pancreas1k <- RunSCVELO(srt = pancreas1k, group_by = "SubCellType", liner_reduction = "PCA",
    nonliner_reduction = "UMAP", return_seurat = TRUE)
VelocityPlot(pancreas1k, reduction = "UMAP", group_by = "SubCellType")
```

<img src="man/figures/README-RunSCVELO-1.png" width="100%" style="display: block; margin: auto;" />

``` r
VelocityPlot(pancreas1k, reduction = "UMAP", plot_type = "stream")
```

<img src="man/figures/README-RunSCVELO-2.png" width="100%" style="display: block; margin: auto;" />

### Differential expression analysis

``` r
library(dplyr)
pancreas1k <- RunDEtest(pancreas1k, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
VolcanoPlot(pancreas1k, group_by = "CellType")
```

<img src="man/figures/README-RunDEtest-1.png" width="100%" style="display: block; margin: auto;" />

``` r
DEGs <- filter(pancreas1k@tools$DEtest_CellType$AllMarkers_wilcox, avg_log2FC > 1 &
    p_val_adj < 0.05)
ht <- ExpHeatmap(pancreas1k, features = DEGs$gene, feature_split = DEGs$group1, cell_split_by = "CellType",
    species = "Mus_musculus", anno_enrichmnet = TRUE, anno_features = TRUE, row_title_size = 0,
    height = 5, width = 7)
print(ht$plot)
```

<img src="man/figures/README-DEGsPlot-1.png" width="100%" style="display: block; margin: auto;" />

### Enrichment analysis(over-representation)

``` r
pancreas1k <- RunEnrichment(pancreas1k, group_by = "CellType", enrichment = "GO_BP",
    species = "Mus_musculus", DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
EnrichmentPlot(pancreas1k, group_by = "CellType", group_use = c("Ductal", "Endocrine"),
    plot_type = "bar")
```

<img src="man/figures/README-RunEnrichment-1.png" width="100%" style="display: block; margin: auto;" />

``` r
EnrichmentPlot(pancreas1k, group_by = "CellType", group_use = c("Ductal", "Endocrine"),
    plot_type = "wordcloud")
```

<img src="man/figures/README-RunEnrichment-2.png" width="100%" style="display: block; margin: auto;" />

### Enrichment analysis(GSEA)

``` r
pancreas1k <- RunGSEA(pancreas1k, group_by = "CellType", enrichment = "GO_BP", species = "Mus_musculus")
GSEAPlot(pancreas1k, group_by = "CellType", group_use = "Endocrine")
```

<img src="man/figures/README-RunGSEA-1.png" width="100%" style="display: block; margin: auto;" />

``` r
GSEAPlot(pancreas1k, group_by = "CellType", group_use = "Endocrine", geneSetID = "GO:0007186")
```

<img src="man/figures/README-RunGSEA-2.png" width="100%" style="display: block; margin: auto;" />

### Dynamic features analysis

``` r
pancreas1k <- RunSlingshot(pancreas1k, group.by = "SubCellType", reduction = "UMAP",
    show_plot = TRUE)
```

<img src="man/figures/README-RunDynamicFeatures-1.png" width="100%" style="display: block; margin: auto;" />

``` r
pancreas1k <- RunDynamicFeatures(pancreas1k, lineages = c("Lineage1", "Lineage2"),
    n_candidates = 200)
```

``` r
ht <- DynamicHeatmap(srt = pancreas1k, lineages = c("Lineage1", "Lineage2"), cell_annotation = "SubCellType",
    n_split = 5, reverse_ht = "Lineage1", species = "Mus_musculus", anno_enrichmnet = TRUE,
    anno_features = TRUE, height = 5, width = 7, use_raster = FALSE)
print(ht$plot)
```

<img src="man/figures/README-DynamicHeatmap-1.png" width="100%" style="display: block; margin: auto;" />

``` r
DynamicPlot(pancreas1k, lineages = c("Lineage1", "Lineage2"), features = c("Plk1",
    "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"), group.by = "SubCellType", compare_lineages = TRUE,
    compare_features = FALSE)
```

<img src="man/figures/README-DynamicPlot-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpVlnPlot(pancreas1k, features = c("Sox9", "Neurod2", "Isl1", "Rbp4"), group.by = "SubCellType",
    bg.by = "CellType", comparisons = list(c("Ductal", "Ngn3 low EP"), c("Ngn3 high EP",
        "Pre-endocrine"), c("Alpha", "Beta")), multiplegroup_comparisons = TRUE)
```

<img src="man/figures/README-ExpVlnPlot-1.png" width="100%" style="display: block; margin: auto;" />
