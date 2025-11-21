#' Seurat V4/V5 Compatibility Helpers
#'
#' These functions provide compatibility layer between Seurat V4 and V5
#' @keywords internal
NULL

#' Get Assay Data with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param slot Slot name (V4 style) or layer name (V5 style)
#' @param assay Assay name
#' @param layer Layer name (V5 style, optional)
#' @keywords internal
#' @noRd
compat_get_assay_data <- function(object, slot = NULL, layer = NULL, assay = NULL) {
  seurat_version <- packageVersion("Seurat")

  # Determine which parameter to use
  if (!is.null(layer)) {
    target <- layer
  } else if (!is.null(slot)) {
    target <- slot
  } else {
    target <- "data"
  }

  # Try V5 style first (layer parameter)
  if (seurat_version >= "5.0.0") {
    tryCatch({
      if (!is.null(assay)) {
        return(Seurat::GetAssayData(object, layer = target, assay = assay))
      } else {
        return(Seurat::GetAssayData(object, layer = target))
      }
    }, error = function(e) {
      # Fall back to V4 style if layer doesn't work
      if (!is.null(assay)) {
        return(Seurat::GetAssayData(object, slot = target, assay = assay))
      } else {
        return(Seurat::GetAssayData(object, slot = target))
      }
    })
  } else {
    # V4 style (slot parameter)
    if (!is.null(assay)) {
      return(Seurat::GetAssayData(object, slot = target, assay = assay))
    } else {
      return(Seurat::GetAssayData(object, slot = target))
    }
  }
}

#' Set Assay Data with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param new.data New data to set
#' @param slot Slot name (V4 style) or layer name (V5 style)
#' @param assay Assay name
#' @param layer Layer name (V5 style, optional)
#' @keywords internal
#' @noRd
compat_set_assay_data <- function(object, new.data, slot = NULL, layer = NULL, assay = NULL) {
  seurat_version <- packageVersion("Seurat")

  # Determine which parameter to use
  if (!is.null(layer)) {
    target <- layer
  } else if (!is.null(slot)) {
    target <- slot
  } else {
    target <- "data"
  }

  # Try V5 style first (layer parameter)
  if (seurat_version >= "5.0.0") {
    tryCatch({
      if (!is.null(assay)) {
        return(Seurat::SetAssayData(object, layer = target, new.data = new.data, assay = assay))
      } else {
        return(Seurat::SetAssayData(object, layer = target, new.data = new.data))
      }
    }, error = function(e) {
      # Fall back to V4 style
      if (!is.null(assay)) {
        return(Seurat::SetAssayData(object, slot = target, new.data = new.data, assay = assay))
      } else {
        return(Seurat::SetAssayData(object, slot = target, new.data = new.data))
      }
    })
  } else {
    # V4 style (slot parameter)
    if (!is.null(assay)) {
      return(Seurat::SetAssayData(object, slot = target, new.data = new.data, assay = assay))
    } else {
      return(Seurat::SetAssayData(object, slot = target, new.data = new.data))
    }
  }
}

#' Check Seurat version
#'
#' @return Character indicating Seurat major version
#' @keywords internal
#' @noRd
get_seurat_version <- function() {
  version <- packageVersion("Seurat")
  if (version >= "5.0.0") {
    return("v5")
  } else if (version >= "4.0.0") {
    return("v4")
  } else {
    return("v3")
  }
}

#' Normalize data with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param ... Additional arguments passed to NormalizeData
#' @keywords internal
#' @noRd
compat_normalize_data <- function(object, ...) {
  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    # V5 may have different defaults
    return(Seurat::NormalizeData(object, ...))
  } else {
    # V4 style
    return(Seurat::NormalizeData(object, ...))
  }
}

#' Scale data with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param ... Additional arguments passed to ScaleData
#' @keywords internal
#' @noRd
compat_scale_data <- function(object, ...) {
  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    return(Seurat::ScaleData(object, ...))
  } else {
    return(Seurat::ScaleData(object, ...))
  }
}

#' Find integration anchors with V4/V5 compatibility
#'
#' @param object.list List of Seurat objects
#' @param ... Additional arguments
#' @keywords internal
#' @noRd
compat_find_integration_anchors <- function(object.list, ...) {
  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    # V5 style - may need to handle layers differently
    return(Seurat::FindIntegrationAnchors(object.list = object.list, ...))
  } else {
    # V4 style
    return(Seurat::FindIntegrationAnchors(object.list = object.list, ...))
  }
}

#' Integrate data with V4/V5 compatibility
#'
#' @param anchorset Anchor set object
#' @param ... Additional arguments
#' @keywords internal
#' @noRd
compat_integrate_data <- function(anchorset, ...) {
  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    # V5 uses IntegrateLayers for some workflows
    # But IntegrateData should still work for backwards compatibility
    return(Seurat::IntegrateData(anchorset = anchorset, ...))
  } else {
    # V4 style
    return(Seurat::IntegrateData(anchorset = anchorset, ...))
  }
}

#' Run differential expression with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param ... Additional arguments passed to FindMarkers
#' @keywords internal
#' @noRd
compat_find_markers <- function(object, ...) {
  seurat_version <- packageVersion("Seurat")

  if (seurat_version >= "5.0.0") {
    # V5 uses presto by default and has different logFC calculation
    # Users should be aware that logFC values may differ
    return(Seurat::FindMarkers(object, ...))
  } else {
    # V4 style
    return(Seurat::FindMarkers(object, ...))
  }
}

#' SCTransform with V4/V5 compatibility
#'
#' @param object Seurat object
#' @param vst.flavor SCTransform flavor ("v1" or "v2")
#' @param ... Additional arguments
#' @keywords internal
#' @noRd
compat_sctransform <- function(object, vst.flavor = NULL, ...) {
  seurat_version <- packageVersion("Seurat")

  # Set default vst.flavor based on version if not specified
  if (is.null(vst.flavor)) {
    if (seurat_version >= "5.0.0") {
      vst.flavor <- "v2"  # V5 default
    } else {
      vst.flavor <- "v1"  # V4 default
    }
  }

  return(Seurat::SCTransform(object, vst.flavor = vst.flavor, ...))
}
