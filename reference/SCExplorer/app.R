# !/usr/bin/env Rscript
library(SCP)
initial_plot_dpi <- 100
initial_panel_dpi <- 300
initial_arrange <- "Row"
initial_ncol <- 3
initial_size <- 4
initial_coExp <- "No"
initial_theme2 <- "theme_scp"
initial_theme1 <- "theme_scp"
initial_palette2 <- "Spectral"
initial_palette1 <- "Paired"
initial_label <- "No"
initial_slot <- NULL
initial_assay <- NULL
initial_feature <- NULL
initial_metaname <- NULL
initial_reduction <- NULL
initial_dataset <- NULL
title <- "SCExplorer"
MetaFile <- "Meta.hdf5"
DataFile <- "Data.hdf5"
base_dir <- "/home/runner/work/SCP/SCP/docs/reference/SCExplorer"

check_R(c("shiny", "shinycssloaders"))
library(shiny)
DataFile_full <- paste0(base_dir, "/", DataFile)
MetaFile_full <- paste0(base_dir, "/", MetaFile)

data_group <- rhdf5::h5ls(DataFile_full)$group
meta_group <- rhdf5::h5ls(MetaFile_full)$group
group <- intersect(data_group, meta_group)
group <- group[group != "/"]
group <- as.character(sapply(group, function(x) substr(x, 2, nchar(x))))
if (length(group) == 0) {
  stop("Can not find the shared group names in the DataFile and the MetaFile. They may not correspond to the same project.")
}
if (is.null(initial_dataset)) {
  initial_dataset <- group[1]
}
if (substr(initial_dataset, 1, 1) == "/") {
  initial_dataset <- substr(initial_dataset, 2, nchar(initial_dataset))
}
if (!initial_dataset %in% group) {
  stop("Dataset ", group, " is not in the DataFile and the MetaFile")
}

assays <- unique(na.omit(sapply(strsplit(data_group[grep(initial_dataset, data_group)], "/"), function(x) x[3])))
slots <- unique(na.omit(sapply(strsplit(data_group[grep(initial_dataset, data_group)], "/"), function(x) x[4])))
if (is.null(initial_assay)) {
  initial_assay <- as.character(rhdf5::h5read(DataFile_full, name = paste0("/", initial_dataset, "/Default_assay")))
}
if (is.null(initial_slot)) {
  initial_slot <- ifelse("data" %in% slots, "data", slots[1])
}
if (!initial_assay %in% assays) {
  stop("initial_assay is not in the dataset ", initial_dataset, " in the DataFile")
}
if (!initial_slot %in% slots) {
  stop("initial_slot is not in the dataset ", initial_slot, " in the DataFile")
}

data <- HDF5Array::TENxMatrix(filepath = DataFile_full, group = paste0("/", initial_dataset, "/", initial_assay, "/", initial_slot))
all_features <- colnames(data)

meta_struc <- rhdf5::h5ls(MetaFile_full)
meta_features_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", initial_dataset, "/metadata.stat/asfeatures"))
meta_groups_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", initial_dataset, "/metadata.stat/asgroups"))
reduction_name <- meta_struc[meta_struc$group == paste0("/", initial_dataset, "/reductions"), "name"]
default_reduction <- as.character(rhdf5::h5read(MetaFile_full, name = paste0("/", initial_dataset, "/reductions.stat/Default_reduction")))

if (is.null(initial_reduction)) {
  initial_reduction <- default_reduction
}
if (is.null(initial_metaname)) {
  initial_metaname <- "orig.ident"
}
if (is.null(initial_feature)) {
  initial_feature <- meta_features_name[1]
}

initial_srt_tmp <- FetchH5(
  DataFile = DataFile_full, MetaFile = MetaFile_full, name = initial_dataset,
  features = initial_feature, slot = initial_slot, assay = initial_assay,
  metanames = initial_metaname, reduction = initial_reduction
)

initial_p1_dim <- ClassDimPlot(
  srt = initial_srt_tmp, group.by = initial_metaname, reduction = initial_reduction,
  label = ifelse(initial_label == "Yes", TRUE, FALSE), palette = initial_palette1, theme_use = initial_theme1,
  ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE), force = TRUE
)
initial_p1_dim <- panel_fix(initial_p1_dim, height = initial_size, raster = TRUE, dpi = initial_panel_dpi, verbose = FALSE)
initial_p2_dim <- ExpDimPlot(
  srt = initial_srt_tmp, features = initial_feature, reduction = initial_reduction, slot = "data",
  calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE), palette = initial_palette2, theme_use = initial_theme2,
  ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE)
)
initial_p2_dim <- panel_fix(initial_p2_dim, height = initial_size, raster = TRUE, dpi = initial_panel_dpi, verbose = FALSE)
initial_p2_vln <- ExpVlnPlot(
  srt = initial_srt_tmp, features = initial_feature, group.by = initial_metaname,
  calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE), palette = initial_palette2,
  ncol = initial_ncol, byrow = ifelse(initial_arrange == "Row", TRUE, FALSE), force = TRUE
)
initial_p2_vln <- panel_fix(initial_p2_vln, height = initial_size, raster = TRUE, dpi = initial_panel_dpi, verbose = FALSE)

initial_plot3d <- max(sapply(names(initial_srt_tmp@reductions), function(r) dim(initial_srt_tmp[[r]])[2])) >= 3
if (isTRUE(initial_plot3d)) {
  initial_p1_3d <- ClassDimPlot3D(
    srt = initial_srt_tmp, group.by = initial_metaname, reduction = initial_reduction, palette = initial_palette1,
    force = TRUE
  )
  initial_p2_3d <- ExpDimPlot3D(
    srt = initial_srt_tmp, features = initial_feature, reduction = initial_reduction,
    calculate_coexp = ifelse(initial_coExp == "Yes", TRUE, FALSE)
  )
} else {
  initial_p1_3d <- initial_p2_3d <- NULL
}

# ui ----------------------------------------------------------------------
ui <- fluidPage(
  navbarPage(
    title = title,
    tabPanel(
      "Cell classification view",
      sidebarPanel(
        width = 3,
        selectInput(
          inputId = "dataset1",
          label = "Select a dataset",
          choices = group,
          selected = substr(initial_dataset, 2, nchar(initial_dataset))
        ),
        selectInput(
          inputId = "reduction1",
          label = "Select a reduction",
          choices = reduction_name,
          selected = initial_reduction
        ),
        selectInput(
          inputId = "class1",
          label = "Select a categorical variable",
          choices = meta_groups_name,
          selected = initial_metaname
        ),
        selectInput(
          inputId = "split1",
          label = "Select a variable for splitting",
          choices = c("None", meta_groups_name),
          selected = "None"
        ),
        radioButtons(
          inputId = "label1",
          label = "Whether labeled",
          choices = c("Yes", "No"),
          selected = initial_label,
          inline = TRUE
        ),
        selectInput(
          inputId = "palette1",
          label = "Select a palette",
          choices = names(palette_list),
          selected = initial_palette1
        ),
        selectInput(
          inputId = "theme1",
          label = "Select a theme",
          choices = c("theme_scp", "theme_blank"),
          selected = initial_theme1
        ),
        numericInput(
          inputId = "size1",
          label = "Panel size",
          value = initial_size,
          min = 1,
          max = 10,
          step = 0.1,
          width = "150px"
        ),
        numericInput(
          inputId = "panel_dpi1",
          label = "Resolution of the panel",
          value = initial_panel_dpi,
          min = 50,
          max = 1000,
          step = 50,
          width = "150px"
        ),
        numericInput(
          inputId = "plot_dpi1",
          label = "Resolution of the plot",
          value = initial_plot_dpi,
          min = 50,
          max = 1000,
          step = 50,
          width = "150px"
        ),
        numericInput(
          inputId = "ncol1",
          label = "Number of columns",
          value = initial_ncol,
          min = 1,
          max = 100,
          step = 1,
          width = "150px"
        ),
        radioButtons(
          inputId = "arrange1",
          label = "Arrange by",
          choices = c("Row", "Column"),
          selected = initial_arrange,
          inline = TRUE
        ),
        actionButton(inputId = "submit1", label = "Submit", icon("play")),
        helpText("Click the submit button to update the plot displayed in the right panel.")
      ),
      mainPanel(
        width = 9,
        fluidPage(
          tabsetPanel(
            tabPanel(
              title = "2D plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  shinycssloaders::withSpinner(plotOutput("plot1", height = "100%", width = "100%"))
                )
              )
            ),
            tabPanel(
              title = "3D plot",
              column(
                width = 12, offset = 0, style = "padding:0px;margin:0%",
                div(
                  style = "overflow-x: auto;",
                  shinycssloaders::withSpinner(plotly::plotlyOutput("plot1_3d", height = "100%", width = "100%"))
                )
              )
            )
          )
        )
      )
    ),
    tabPanel(
      "Feature expression view",
      sidebarPanel(
        width = 3,
        selectInput(
          inputId = "dataset2",
          label = "Select a dataset",
          choices = group,
          selected = substr(initial_dataset, 2, nchar(initial_dataset))
        ),
        selectInput(
          inputId = "reduction2",
          label = "Select a reduction",
          choices = reduction_name,
          selected = initial_reduction
        ),
        selectInput(
          inputId = "assays2",
          label = "Select a assay",
          choices = assays,
          selected = initial_assay
        ),
        selectInput(
          inputId = "slots2",
          label = "Select a slot",
          choices = slots,
          selected = initial_slot
        ),
        selectizeInput(
          inputId = "features2",
          label = "Select genes",
          choices = NULL,
          selected = initial_feature,
          multiple = TRUE,
          options = list(maxOptions = 20, maxItems = 20)
        ),
        textAreaInput(
          inputId = "gene_area",
          label = "Input genes",
          height = "200px",
          placeholder = paste(sample(all_features, 4), collapse = "
")
        ),
        selectInput(
          inputId = "split2",
          label = "Select a variable for splitting",
          choices = c("None", meta_groups_name),
          selected = "None"
        ),
        selectInput(
          inputId = "group2",
          label = "Select a grouping variable",
          choices = c("None", meta_groups_name),
          selected = initial_metaname
        ),
        radioButtons(
          inputId = "coExp2",
          label = "Calculate co-expression?",
          choices = c("Yes", "No"),
          selected = initial_coExp,
          inline = TRUE
        ),
        selectInput(
          inputId = "palette2",
          label = "Select a palette",
          choices = names(palette_list),
          selected = initial_palette2
        ),
        selectInput(
          inputId = "theme2",
          label = "Select a theme",
          choices = c("theme_scp", "theme_blank"),
          selected = initial_theme2
        ),
        numericInput(
          inputId = "size2",
          label = "Panel size",
          value = initial_size,
          min = 1,
          max = 10,
          step = 0.1,
          width = "150px"
        ),
        numericInput(
          inputId = "panel_dpi2",
          label = "Resolution of the panel",
          value = initial_panel_dpi,
          min = 50,
          max = 1000,
          step = 50,
          width = "150px"
        ),
        numericInput(
          inputId = "plot_dpi2",
          label = "Resolution of the plot",
          value = initial_plot_dpi,
          min = 50,
          max = 1000,
          step = 50,
          width = "150px"
        ),
        numericInput(
          inputId = "ncol2",
          label = "Number of columns",
          value = initial_ncol,
          min = 1,
          max = 100,
          step = 1,
          width = "150px"
        ),
        radioButtons(
          inputId = "arrange2",
          label = "Arrange by",
          choices = c("Row", "Column"),
          selected = initial_arrange,
          inline = TRUE
        ),
        actionButton(inputId = "submit2", label = "Submit", icon("play")),
        helpText("Click the submit button to update the plot displayed in the right panel.")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel(
            title = "2D plot",
            column(
              width = 12, offset = 0, style = "padding:0px;margin:0%",
              div(
                style = "overflow-x: auto;",
                shinycssloaders::withSpinner(plotOutput("plot2", height = "100%", width = "100%"))
              )
            )
          ),
          tabPanel(
            title = "3D plot",
            column(
              width = 12, offset = 0, style = "padding:0px;margin:0%",
              div(
                style = "overflow-x: auto;",
                shinycssloaders::withSpinner(plotly::plotlyOutput("plot2_3d", height = "100%", width = "100%"))
              )
            )
          ),
          tabPanel(
            title = "Violin plot",
            column(
              width = 12, offset = 0, style = "padding:0px;margin:0%",
              div(
                style = "overflow-x: auto;",
                shinycssloaders::withSpinner(plotOutput("plot2_vln", height = "100%", width = "100%"))
              )
            )
          )
        )
      )
    )
  )
)

# server ------------------------------------------------------------------
server <- function(input, output, session) {

  # initial  ----------------------------------------------------------------
  updateSelectizeInput(session, "features2", choices = c(meta_features_name, all_features), selected = initial_feature, server = TRUE)

  output$plot1 <- renderPlot(
    {
      initial_p1_dim
    },
    width = attr(initial_p1_dim, "size")$width * initial_plot_dpi,
    height = attr(initial_p1_dim, "size")$height * initial_plot_dpi,
    res = initial_plot_dpi
  )

  output$plot1_3d <- plotly::renderPlotly({
    initial_p1_3d
  })

  output$plot2 <- renderPlot(
    {
      initial_p2_dim
    },
    width = attr(initial_p2_dim, "size")$width * initial_plot_dpi,
    height = attr(initial_p2_dim, "size")$height * initial_plot_dpi,
    res = initial_plot_dpi
  )

  output$plot2_3d <- plotly::renderPlotly({
    initial_p2_3d
  })

  output$plot2_vln <- renderPlot(
    {
      initial_p2_vln
    },
    width = attr(initial_p2_vln, "size")$width * initial_plot_dpi,
    height = attr(initial_p2_vln, "size")$height * initial_plot_dpi,
    res = initial_plot_dpi
  )

  # change dataset  ----------------------------------------------------------------
  observeEvent(input$dataset1, {
    meta_groups_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset1, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset1, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset1, "/reductions.stat/Default_reduction")))
    updateSelectInput(session, "reduction1", choices = reduction_name, selected = default_reduction)
    updateSelectInput(session, "class1", choices = meta_groups_name, selected = "orig.ident")
    updateSelectInput(session, "split1", choices = c("None", meta_groups_name), selected = "None")
  })

  observeEvent(input$dataset2, {
    assays <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[3])))
    slots <- unique(na.omit(sapply(strsplit(data_group[grep(input$dataset2, data_group)], "/"), function(x) x[4])))
    default_assay <- as.character(rhdf5::h5read(DataFile_full, name = paste0("/", input$dataset2, "/Default_assay")))
    default_slot <- ifelse("data" %in% slots, "data", slots[1])
    updateSelectInput(session, "assays2", choices = assays, selected = default_assay)
    updateSelectInput(session, "slots2", choices = slots, selected = default_slot)

    data <- HDF5Array::TENxMatrix(filepath = DataFile_full, group = paste0("/", input$dataset2, "/", default_assay, "/", default_slot))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset2, "/metadata.stat/asfeatures"))
    meta_groups_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset2, "/metadata.stat/asgroups"))
    reduction_name <- meta_struc[meta_struc$group == paste0("/", input$dataset2, "/reductions"), "name"]
    default_reduction <- as.character(rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset2, "/reductions.stat/Default_reduction")))
    updateSelectInput(session, "reduction2", choices = reduction_name, selected = default_reduction)
    updateSelectizeInput(session, "features2",
      choices = c(meta_features_name, all_features), selected = meta_features_name[1],
      options = list(maxOptions = 20, maxItems = 20), server = TRUE
    )
    updateSelectInput(session, "split2", choices = c("None", meta_groups_name), selected = "None")
    updateSelectInput(session, "group2", choices = meta_groups_name, selected = "orig.ident")
  })

  # submit1  ----------------------------------------------------------------
  observeEvent(input$submit1, ignoreInit = FALSE, {
    if (input$split1 == "None") {
      split1 <- NULL
    } else {
      split1 <- input$split1
    }

    # message("DataFile_full:", DataFile_full)
    # message("MetaFile_full:", MetaFile_full)
    # message("input$dataset1:", input$dataset1)
    # message("initial_slot:", initial_slot)
    # message("c(input$class1, split1):", c(input$class1, split1))
    # message("input$reduction1:", input$reduction1)

    srt_tmp <- FetchH5(
      DataFile = DataFile_full, MetaFile = MetaFile_full, name = input$dataset1,
      features = NULL, slot = initial_slot, assay = NULL,
      metanames = c(input$class1, split1), reduction = input$reduction1
    )

    p1_dim <- ClassDimPlot(
      srt = srt_tmp, group.by = input$class1, split.by = split1, reduction = input$reduction1,
      label = ifelse(input$label1 == "Yes", TRUE, FALSE), palette = input$palette1, theme_use = input$theme1,
      ncol = input$ncol1, byrow = ifelse(input$arrange1 == "Row", TRUE, FALSE), force = TRUE
    )
    p1_dim <- panel_fix(p1_dim, height = input$size1, raster = TRUE, dpi = input$panel_dpi1, verbose = FALSE)

    output$plot1 <- renderPlot(
      {
        p1_dim
      },
      width = attr(p1_dim, "size")$width * input$plot_dpi1,
      height = attr(p1_dim, "size")$height * input$plot_dpi1,
      res = input$plot_dpi1
    )

    plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
    if (isTRUE(plot3d)) {
      p1_3d <- ClassDimPlot3D(
        srt = srt_tmp, group.by = input$class1, reduction = input$reduction1, palette = input$palette1,
        force = TRUE
      )
    } else {
      p1_3d <- NULL
    }

    output$plot1_3d <- plotly::renderPlotly({
      p1_3d
    })
  })

  # submit2  ----------------------------------------------------------------
  observeEvent(input$submit2, ignoreInit = FALSE, {
    data <- HDF5Array::TENxMatrix(filepath = DataFile_full, group = paste0("/", input$dataset2, "/", initial_assay, "/", initial_slot))
    all_features <- colnames(data)
    meta_features_name <- rhdf5::h5read(MetaFile_full, name = paste0("/", input$dataset2, "/metadata.stat/asfeatures"))

    if (input$split2 == "None") {
      split2 <- NULL
    } else {
      split2 <- input$split2
    }
    input_features <- input$features2
    if (is.null(input_features)) {
      message("input gene is null")
      input_features <- meta_features_name[1]
    }
    gene_area <- gsub(x = unlist(strsplit(input$gene_area, "(\
)|(\
)", perl = TRUE)), pattern = " ", replacement = "")
    input_features <- c(as.character(input_features), as.character(gene_area))
    input_features <- unique(input_features[input_features %in% c(all_features, meta_features_name)])

    # message("DataFile_full:", DataFile_full)
    # message("MetaFile_full:", MetaFile_full)
    # message("input$dataset2:", input$dataset2)
    # message("input_features:", paste0(input_features,collapse = ","))
    # message("initial_slot:", initial_slot)
    # message("initial_assay:", initial_assay)
    # message("c(input$group2, split2):", c(input$group2, split2))
    # message("input$reduction2:", input$reduction2)

    srt_tmp <- FetchH5(
      DataFile = DataFile_full, MetaFile = MetaFile_full, name = input$dataset2,
      features = input_features, slot = input$slots2, assay = input$assays2,
      metanames = c(input$group2, split2), reduction = input$reduction2
    )

    p2_dim <- ExpDimPlot(
      srt = srt_tmp, features = input_features, split.by = split2, reduction = input$reduction2, slot = "data",
      calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE), palette = input$palette2, theme_use = input$theme2,
      ncol = input$ncol2, byrow = ifelse(input$arrange2 == "Row", TRUE, FALSE)
    )
    print(p2_dim)
    p2_dim <- panel_fix(p2_dim, height = input$size2, raster = TRUE, dpi = input$panel_dpi2, verbose = FALSE)

    output$plot2 <- renderPlot(
      {
        p2_dim
      },
      width = attr(p2_dim, "size")$width * input$plot_dpi2,
      height = attr(p2_dim, "size")$height * input$plot_dpi2,
      res = input$plot_dpi2
    )

    p2_vln <- ExpVlnPlot(
      srt = srt_tmp, features = input_features, group.by = input$group2, split.by = split2,
      calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE), palette = input$palette2,
      ncol = input$ncol2, byrow = ifelse(input$arrange2 == "Row", TRUE, FALSE), force = TRUE
    )
    p2_vln <- panel_fix(p2_vln, height = input$size2, raster = TRUE, dpi = input$panel_dpi2, verbose = FALSE)

    output$plot2_vln <- renderPlot(
      {
        p2_vln
      },
      width = attr(p2_vln, "size")$width * input$plot_dpi2,
      height = attr(p2_vln, "size")$height * input$plot_dpi2,
      res = input$plot_dpi2
    )

    plot3d <- max(sapply(names(srt_tmp@reductions), function(r) dim(srt_tmp[[r]])[2])) >= 3
    if (isTRUE(plot3d)) {
      p2_3d <- ExpDimPlot3D(
        srt = srt_tmp, features = input_features, reduction = input$reduction2,
        calculate_coexp = ifelse(input$coExp2 == "Yes", TRUE, FALSE)
      )
    } else {
      p2_3d <- NULL
    }

    output$plot2_3d <- plotly::renderPlotly({
      p2_3d
    })
  })
}

shinyApp(ui = ui, server = server)
