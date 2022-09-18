
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
-   Pipelines embedded with multiple methods for normalization, feature
    reduction, and cell population identification (standard Seurat
    workflow).
-   Pipelines embedded with multiple data integration methods, including
    Uncorrected, [Seurat](https://github.com/satijalab/seurat),
    [scVI](https://github.com/scverse/scvi-tools),
    [MNN](http://www.bioconductor.org/packages/release/bioc/html/batchelor.html),
    [fastMNN](http://www.bioconductor.org/packages/release/bioc/html/batchelor.html),
    [Harmony](https://github.com/immunogenomics/harmony),
    [Scanorama](https://github.com/brianhie/scanorama),
    [BBKNN](https://github.com/Teichlab/bbknn),
    [CSS](https://github.com/quadbiolab/simspec),
    [LIGER](https://github.com/welch-lab/liger),
    [Conos](https://github.com/kharchenkolab/conos).
-   Multiple single cell downstream analyses such as identification of
    differential features, enrichment analysis, GSEA analysis,
    identification of dynamic features,
    [PAGA](https://github.com/theislab/paga), [RNA
    velocity](https://github.com/theislab/scvelo), etc.
-   Multiple methods for automatic annotation of single-cell data and
    methods for projection between single-cell datasets.
-   High quality data visualization methods.
-   Fast deployment of single-cell data into a [shiny
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

### Requirement for python functions in SCP

To run functions such as `RunSCVELO` or `RunPAGA`, SCP requires python
3.7-3.9 to be installed in the environment.

Check the version of python in the terminal:

``` shell
python3 --version
```

or in the R environment:

``` r
if (!require("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
py <- Sys.which("python3")
reticulate:::python_version(py)
```

Then run `PrepareVirtualEnv` to create a standalone python virtual
environment for SCP and install the necessary packages.

``` r
PrepareVirtualEnv(python = py, pipy_mirror = "https://pypi.tuna.tsinghua.edu.cn/simple", remove_old = TRUE)
```

## Example

### Load the Data

The analysis is based on a subsetted version of [mouse pancreas
data](https://doi.org/10.1242/dev.173849).

``` r
library(SCP)
data("pancreas1k")
ClassDimPlot(
  srt = pancreas1k, group.by = c("CellType", "SubCellType"),
  reduction = "UMAP", theme_use = "theme_blank"
)
```

<img src="README/README-library-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpDimPlot(
  srt = pancreas1k, features = c("Sox9", "Neurog3", "Fev", "Rbp4"),
  reduction = "UMAP", theme_use = "theme_blank"
)
```

<img src="README/README-library-2.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpDotPlot(
  srt = pancreas1k,
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

<img src="README/README-library-3.png" width="100%" style="display: block; margin: auto;" />

### CellQC

``` r
pancreas1k <- RunCellQC(srt = pancreas1k)
ClassDimPlot(srt = pancreas1k, group.by = "CellQC", reduction = "UMAP")
```

<img src="README/README-RunCellQC-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassStatPlot(srt = pancreas1k, stat.by = "CellQC", group.by = "CellType", label = TRUE)
```

<img src="README/README-RunCellQC-2.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassStatPlot(
  srt = pancreas1k,
  stat.by = c(
    "db_qc", "outlier_qc", "umi_qc", "gene_qc",
    "mito_qc", "ribo_qc", "ribo_mito_ratio_qc", "species_qc"
  ),
  plot_type = "upset", stat_level = "Fail"
)
```

<img src="README/README-RunCellQC-3.png" width="100%" style="display: block; margin: auto;" />

### Standard pipeline in SCP

``` r
pancreas1k <- Standard_SCP(srt = pancreas1k)
ClassDimPlot(
  srt = pancreas1k, group.by = c("CellType", "SubCellType"),
  reduction = "StandardUMAP2D", theme_use = "theme_blank"
)
```

<img src="README/README-Standard_SCP-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassDimPlot3D(srt = pancreas1k, group.by = "SubCellType")
```

![ClassDimPlot3D](README/README-ClassDimPlot3D-1.png)

``` r
ExpDimPlot3D(srt = pancreas1k, features = c("Sox9", "Neurog3", "Fev", "Rbp4"))
```

![ExpDimPlot3D](README/README-ExpDimPlot3D-1.png)

### Integration pipeline in SCP

Example data for integration is [panc8(eight human pancreas
datasets)](https://github.com/satijalab/seurat-data)

``` r
if (!require("SeuratData", quietly = TRUE)) {
  devtools::install_github("satijalab/seurat-data")
}
library(SeuratData)
suppressWarnings(InstallData("panc8"))
data("panc8")
panc8 <- Integration_SCP(srtMerge = panc8, batch = "tech", integration_method = "Seurat")
panc8 <- Integration_SCP(srtMerge = panc8, batch = "tech", integration_method = "Harmony", nonliner_reduction = "pacmap")
ClassDimPlot(
  srt = panc8, group.by = c("celltype", "tech"), reduction = "SeuratUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)
```

<img src="README/README-Integration_SCP-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ClassDimPlot(
  srt = panc8, group.by = c("celltype", "tech"), reduction = "HarmonyPACMAP2D",
  title = "Harmony", theme_use = "theme_blank"
)
```

<img src="README/README-Integration_SCP-2.png" width="100%" style="display: block; margin: auto;" />

### Cell projection between single-cell datasets

``` r
library(stringr)
panc8 <- RenameFeatures(srt = panc8, newnames = make.unique(str_to_title(rownames(panc8))))
pancreas1k <- RunKNNMap(srt_query = pancreas1k, srt_ref = panc8, ref_umap = "SeuratUMAP2D")
ProjectionPlot(
  srt_query = pancreas1k, srt_ref = panc8,
  query_group = "SubCellType", ref_group = "celltype"
)
```

<img src="README/README-RunKNNMap-1.png" width="100%" style="display: block; margin: auto;" />

### Cell annotation using bulk RNA-seq datasets

``` r
data("ref_scMCA")
pancreas1k <- RunKNNPredict(srt_query = pancreas1k, bulk_ref = ref_scMCA, filter_lowfreq = 20)
ClassDimPlot(srt = pancreas1k, group.by = "knnpredict_classification", reduction = "UMAP", label = TRUE)
```

<img src="README/README-RunKNNPredict-bulk-1.png" width="100%" style="display: block; margin: auto;" />

### Cell annotation using single-cell datasets

``` r
pancreas1k <- RunKNNPredict(
  srt_query = pancreas1k, srt_ref = panc8,
  ref_group = "celltype", filter_lowfreq = 20
)
ClassDimPlot(srt = pancreas1k, group.by = "knnpredict_classification", reduction = "UMAP", label = TRUE)
```

<img src="README/README-RunKNNPredict-scrna-1.png" width="100%" style="display: block; margin: auto;" />

### PAGA analysis

``` r
pancreas1k <- RunPAGA(
  srt = pancreas1k, group_by = "SubCellType",
  liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE
)
PAGAPlot(srt = pancreas1k, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE)
```

<img src="README/README-RunPAGA-1.png" width="100%" style="display: block; margin: auto;" />

### Velocity analysis

``` r
pancreas1k <- RunSCVELO(
  srt = pancreas1k, group_by = "SubCellType",
  liner_reduction = "PCA", nonliner_reduction = "UMAP", return_seurat = TRUE
)
VelocityPlot(srt = pancreas1k, reduction = "UMAP", group_by = "SubCellType")
```

<img src="README/README-RunSCVELO-1.png" width="100%" style="display: block; margin: auto;" />

``` r
VelocityPlot(srt = pancreas1k, reduction = "UMAP", plot_type = "stream")
```

<img src="README/README-RunSCVELO-2.png" width="100%" style="display: block; margin: auto;" />

### Differential expression analysis

``` r
pancreas1k <- RunDEtest(srt = pancreas1k, group_by = "CellType", only.pos = FALSE, fc.threshold = 1)
VolcanoPlot(srt = pancreas1k, group_by = "CellType")
```

<img src="README/README-RunDEtest-1.png" width="100%" style="display: block; margin: auto;" />

``` r
DEGs <- pancreas1k@tools$DEtest_CellType$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
ht <- ExpHeatmap(
  srt = pancreas1k, features = DEGs$gene, feature_split = DEGs$group1, cell_split_by = "CellType",
  species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  row_title_size = 0, height = 5, width = 7
)
print(ht$plot)
```

<img src="README/README-DEGsPlot-1.png" width="100%" style="display: block; margin: auto;" />

### Enrichment analysis(over-representation)

``` r
pancreas1k <- RunEnrichment(
  srt = pancreas1k, group_by = "CellType", enrichment = "GO_BP", species = "Mus_musculus",
  DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = pancreas1k, group_by = "CellType", group_use = c("Ductal", "Endocrine"),
  plot_type = "bar"
)
```

<img src="README/README-RunEnrichment-1.png" width="100%" style="display: block; margin: auto;" />

``` r
EnrichmentPlot(
  srt = pancreas1k, group_by = "CellType", group_use = c("Ductal", "Endocrine"),
  plot_type = "wordcloud"
)
```

<img src="README/README-RunEnrichment-2.png" width="100%" style="display: block; margin: auto;" />

``` r
EnrichmentPlot(
  srt = pancreas1k, group_by = "CellType", group_use = c("Ductal", "Endocrine"),
  plot_type = "wordcloud", word_type = "feature"
)
```

<img src="README/README-RunEnrichment-3.png" width="100%" style="display: block; margin: auto;" />

### Enrichment analysis(GSEA)

``` r
pancreas1k <- RunGSEA(
  srt = pancreas1k, group_by = "CellType", enrichment = "GO_BP", species = "Mus_musculus",
  DE_threshold = "p_val_adj < 0.05"
)
GSEAPlot(srt = pancreas1k, group_by = "CellType", group_use = "Endocrine")
```

<img src="README/README-RunGSEA-1.png" width="100%" style="display: block; margin: auto;" />

``` r
GSEAPlot(srt = pancreas1k, group_by = "CellType", group_use = "Endocrine", geneSetID = "GO:0007186")
```

<img src="README/README-RunGSEA-2.png" width="100%" style="display: block; margin: auto;" />

### Dynamic features analysis

``` r
pancreas1k <- RunSlingshot(srt = pancreas1k, group.by = "SubCellType", reduction = "UMAP", show_plot = TRUE)
```

<img src="README/README-RunDynamicFeatures-1.png" width="100%" style="display: block; margin: auto;" />

``` r
pancreas1k <- RunDynamicFeatures(srt = pancreas1k, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
```

``` r
ht <- DynamicHeatmap(
  srt = pancreas1k, lineages = c("Lineage1", "Lineage2"), cell_annotation = "SubCellType",
  n_split = 5, reverse_ht = "Lineage1",
  species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  height = 5, width = 7, use_raster = FALSE
)
print(ht$plot)
```

<img src="README/README-DynamicHeatmap-1.png" width="100%" style="display: block; margin: auto;" />

``` r
DynamicPlot(
  srt = pancreas1k, lineages = c("Lineage1", "Lineage2"), group.by = "SubCellType",
  features = c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
  compare_lineages = TRUE, compare_features = FALSE
)
```

<img src="README/README-DynamicPlot-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ExpVlnPlot(
  srt = pancreas1k, group.by = "SubCellType", bg.by = "CellType",
  features = c("Sox9", "Neurod2", "Isl1", "Rbp4"),
  comparisons = list(
    c("Ductal", "Ngn3 low EP"),
    c("Ngn3 high EP", "Pre-endocrine"),
    c("Alpha", "Beta")
  ),
  multiplegroup_comparisons = TRUE
)
```

<img src="README/README-ExpVlnPlot-1.png" width="100%" style="display: block; margin: auto;" />
