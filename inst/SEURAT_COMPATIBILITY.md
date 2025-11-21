# Seurat V4/V5 Compatibility Guide

## Overview

SCP now supports both Seurat V4 (>= 4.2.0) and Seurat V5 (>= 5.0.0). This document explains the compatibility considerations and how to use SCP with different Seurat versions.

## Installation

### With Seurat V4

```r
# Install Seurat V4
install.packages("Seurat")

# Install SCP
devtools::install_github("zhanghao-njmu/SCP")
```

### With Seurat V5

```r
# Install Seurat V5
install.packages("Seurat")

# Or install from GitHub for latest V5 features
devtools::install_github("satijalab/seurat", ref = "main")

# Install SCP
devtools::install_github("zhanghao-njmu/SCP")
```

## Key Differences Between Seurat V4 and V5

### 1. Data Access (slots vs layers)

**Seurat V4:**
```r
# V4 uses slot parameter
counts <- GetAssayData(srt, slot = "counts")
data <- GetAssayData(srt, slot = "data")
scaled <- GetAssayData(srt, slot = "scale.data")
```

**Seurat V5:**
```r
# V5 uses layer parameter (but slot still works for backwards compatibility)
counts <- GetAssayData(srt, layer = "counts")
data <- GetAssayData(srt, layer = "data")
scaled <- GetAssayData(srt, layer = "scale.data")
```

**SCP automatically handles both:**
- SCP's internal functions detect your Seurat version
- All SCP functions work with both V4 and V5
- No code changes required from users

### 2. Integration Workflow

**Seurat V4:**
```r
# V4 integration returns an integrated assay
srt <- IntegrateData(anchorset = anchors)
DefaultAssay(srt) <- "integrated"
```

**Seurat V5:**
```r
# V5 integration can work with layers
srt <- IntegrateLayers(object = srt, method = HarmonyIntegration)
# OR use V4-style for compatibility
srt <- IntegrateData(anchorset = anchors)
```

**Using SCP:**
```r
# SCP's Integration_SCP works with both versions
srt_integrated <- Integration_SCP(
  srtMerge = srt,
  batch = "batch_column",
  integration_method = "Harmony"
)
```

### 3. Differential Expression

**Seurat V4:**
- Standard log-fold change calculation
- Pseudocount added at cell level

**Seurat V5:**
- Uses `presto` package when available (faster)
- Pseudocount added at group level
- **Important:** logFC values may differ from V4

**Impact on SCP:**
```r
# SCP's RunDEtest works with both versions
markers <- RunDEtest(
  srt = srt,
  group_by = "celltype",
  fc.threshold = 1.5
)

# Note: If you need V4-style results in V5, you can:
# 1. Uninstall presto package
# 2. Or interpret logFC values accordingly
```

### 4. SCTransform

**Seurat V4:**
- Default: `vst.flavor = "v1"`

**Seurat V5:**
- Default: `vst.flavor = "v2"` (recommended)
- v2 provides better performance

**Using SCP:**
```r
# SCP automatically uses appropriate defaults
srt <- Standard_SCP(
  srt = srt,
  normalization_method = "SCT"
)

# Or specify explicitly
srt <- Seurat::SCTransform(srt, vst.flavor = "v2")
```

## Common Issues and Solutions

### Issue 1: "layer/slot parameter not recognized"

**Cause:** Code written for one Seurat version running on another

**Solution:** SCP handles this automatically. If you're using direct Seurat functions:

```r
# Instead of:
data <- GetAssayData(srt, slot = "data")  # May fail in V5

# Use SCP's compatibility layer (internal):
# SCP functions automatically handle version differences
```

### Issue 2: Different logFC values

**Cause:** V5 uses different pseudocount strategy

**Solution:**
- This is expected behavior
- V5 logFC estimates may be higher
- Be consistent with Seurat version in your analysis
- Document which version you used

### Issue 3: Integration results differ

**Cause:** V5 has different integration algorithms

**Solution:**
- Use same integration method consistently
- SCP's `Integration_SCP` provides consistent interface
- Results may differ slightly between versions

## Checking Your Seurat Version

```r
# Check installed Seurat version
packageVersion("Seurat")

# Check if V5 features are available
if (packageVersion("Seurat") >= "5.0.0") {
  message("Using Seurat V5")
} else {
  message("Using Seurat V4")
}
```

## Testing Compatibility

SCP includes comprehensive compatibility tests:

```r
# Run compatibility tests
library(testthat)
library(SCP)

# Test Seurat compatibility
test_file("tests/testthat/test-seurat-v5-compatibility.R")

# Test ggplot2 compatibility
test_file("tests/testthat/test-ggplot2-compatibility.R")
```

## Migrating Between Versions

### From V4 to V5

Most code should work without changes. Be aware of:

1. **LogFC differences** in differential expression
2. **Integration** may produce slightly different results
3. **SCTransform** default changed to v2

```r
# If you need V4-compatible results in V5
srt <- SCTransform(srt, vst.flavor = "v1")
```

### From V5 to V4

Should work seamlessly, as V5 maintains backwards compatibility.

## Best Practices

1. **Document your Seurat version**
   ```r
   sessionInfo()  # Always include in your analysis
   ```

2. **Be consistent within a project**
   - Don't mix V4 and V5 results
   - Use same version for all samples

3. **Test with your data**
   ```r
   # Quick test
   data("pancreas_sub")
   srt <- RunCellQC(pancreas_sub)
   # Should work regardless of version
   ```

4. **Use SCP's high-level functions**
   - They handle version differences automatically
   - More stable across Seurat updates

## Getting Help

If you encounter compatibility issues:

1. Check your Seurat version: `packageVersion("Seurat")`
2. Update to latest SCP: `devtools::install_github("zhanghao-njmu/SCP")`
3. Run compatibility tests
4. Report issues at: https://github.com/zhanghao-njmu/SCP/issues

Include in your report:
- Seurat version
- SCP version
- Minimal reproducible example
- Error message

## Additional Resources

- [Seurat V5 Announcement](https://satijalab.org/seurat/articles/announcements.html)
- [Seurat V5 Documentation](https://satijalab.org/seurat/)
- [SCP Documentation](https://zhanghao-njmu.github.io/SCP/)

---

**Last Updated:** 2024-11
**SCP Version:** 0.5.6+
**Supported Seurat Versions:** 4.2.0 - 5.x.x
