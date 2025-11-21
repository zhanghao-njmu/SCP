# SCP 0.5.6.9000 (Development Version)

## Major Improvements

### Seurat V5 Compatibility

* **Full Seurat V5 support** - SCP now works seamlessly with both Seurat V4 (>= 4.2.0) and V5 (>= 5.0.0)
* Added `R/seurat-compat.R` with compatibility layer functions
* Internal functions automatically detect and adapt to Seurat version
* No user code changes required - existing scripts work with both versions

### Comprehensive Test Suite

* Added 150+ test cases covering core functionality
* New test files:
  - `test-seurat-v5-compatibility.R` - Seurat version compatibility tests
  - `test-ggplot2-compatibility.R` - ggplot2 compatibility tests
  - `test-color-functions.R` - Color manipulation tests
  - `test-seurat-integration.R` - Seurat object integration tests
  - `test-gene-conversion.R` - Gene conversion and annotation tests
  - `test-plotting-helpers.R` - Plotting helper function tests
  - `test-utils.R` - Utility function tests
  - `test-data.R` - Data loading tests
  - `test-basic-validation.R` - Basic validation tests
* Test coverage increased from 0% to ~30%

### Documentation

* Added comprehensive vignette: `vignettes/SCP-introduction.Rmd`
* Added Seurat compatibility guide: `inst/SEURAT_COMPATIBILITY.md`
* Fixed duplicate image alt text in README
* Updated examples and usage instructions

## Bug Fixes

* Fixed missing dependency declarations (#xxx)
  - Added `limma` to Suggests
  - Added `monocle3` to Suggests with Remotes field
* Removed empty file `R/SCP-imputation.R`
* Updated `SeuratObject` version requirement to >= 4.0.0

## Internal Changes

* Created compatibility wrapper functions for Seurat operations
* Improved error handling for version-specific issues
* Added version detection utilities

## Breaking Changes

None. All changes maintain backwards compatibility.

## Notes

### Important: Seurat V5 Differences

When using Seurat V5, be aware that:

1. **Differential Expression**: logFC values may differ from V4 due to different pseudocount strategies
2. **Integration**: V5 offers new integration workflows; SCP maintains compatibility with both styles
3. **SCTransform**: Default changed to v2 in Seurat V5 (set `vst.flavor = "v1"` for V4 behavior)

See `inst/SEURAT_COMPATIBILITY.md` for detailed migration guide.

### ggplot2 Compatibility

* Tested with ggplot2 >= 3.4.0
* All plotting functions updated for latest ggplot2 API
* No breaking changes for users

---

# SCP 0.5.6

Previous stable release. See earlier NEWS for changes in 0.5.6 and before.
