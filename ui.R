library(shiny)
library(leaflet)
library(shinydashboard)
library(shinyjs)
library(ggplot2)

# Read the species data
species_data <- read.csv("IWJV/data/spp_list_04sep2024.csv")
species_choices <- setNames(paste0("Species ", 1:35), species_data$common.name)

ui <- dashboardPage(
  dashboardHeader(title = "Management Scenario Tool"),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      fileInput("shapefile", "Upload Shapefile (.zip containing .shp, .shx, .dbf, .prj)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".shx", ".dbf", ".prj")),
      selectInput("species", "Choose Species", choices = species_choices),
      sliderInput("tree_canopy_cover", "Tree Canopy Cover (%)", min = -100, max = 100, value = 0, step = 5),
      sliderInput("qmd_rmrs", "Quadratic Mean Diameter (QMD) RMRS (%)", min = -100, max = 100, value = 0, step = 5),
      sliderInput("tpa_dead", "Trees Per Acre (TPA) Dead (%)", min = -100, max = 100, value = 0, step = 5),
      sliderInput("tpa_live", "Trees Per Acre (TPA) Live (%)", min = -100, max = 100, value = 0, step = 5),
      actionButton("runModel", "Run Model", class = "btn-primary", icon = icon("play"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    fluidRow(
      box(
        title = "Abundance Map",
        status = "primary",
        solidHeader = TRUE,
        width = 12,
        leafletOutput("abundanceMap"),
        downloadButton("downloadMap", "Download Map", class = "btn-info")
      )
    ),
    fluidRow(
      box(
        title = "Abundance Statistics",
        status = "warning",
        solidHeader = TRUE,
        width = 6,
        tableOutput("abundanceStats"),
        downloadButton("downloadStats", "Download Statistics", class = "btn-info")
      ),
      box(
        title = "Species Abundance",
        status = "success",
        solidHeader = TRUE,
        width = 6,
        plotOutput("speciesAbundancePlot"),
        downloadButton("downloadPlot", "Download Plot", class = "btn-info")
      )
    )
  )
)
