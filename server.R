library(shiny)
library(leaflet)
library(raster)
library(sf)
library(dplyr)
library(coda)
library(zip)
library(terra)

mcmc_data <- load("IWJV/data/test_model.R")
species_data <- read.csv("IWJV/data/spp_list_04sep2024.csv")

server <- function(input, output, session) {
  # Load MCMC data
  mcmc_data <- reactive({
    print("Loading MCMC data...")
    load("IWJV/data/test_model.R")
    print("MCMC data loaded successfully.")
    return(mcmc.lst)
  })
  
  # Load raster data
  rasters <- reactive({
    print("Loading raster data...")
    raster_list <- list(
      tree_canopy_cover = rast("IWJV/data/tree_canopy_cover.tif"),
      qmd_rmrs = rast("IWJV/data/qmd_rmrs.tif"),
      evt = rast("IWJV/data/evt.tif"),
      tpa_dead = rast("IWJV/data/tpa_dead.tif"),
      tpa_live = rast("IWJV/data/tpa_live.tif"),
      burn_severity = rast("IWJV/data/burn_severity.tif"),
      veg_departure = rast("IWJV/data/veg_departure.tif"),
      years_since_wildfire = rast("IWJV/data/years_since_wildfire.tif")
    )
    print("Raster data loaded successfully.")
    return(raster_list)
  })
  
  # Function to adjust rasters based on user input
  adjustRasters <- reactive({
    req(input$runModel)
    print("Adjusting rasters based on user input...")
    
    withProgress(message = 'Adjusting rasters', value = 0, {
      original <- rasters()
      adjusted <- list()
      
      total_rasters <- length(original)
      for (i in seq_along(original)) {
        name <- names(original)[i]
        if (name %in% c("tree_canopy_cover", "qmd_rmrs", "tpa_dead", "tpa_live")) {
          adjustment <- input[[name]] / 100 + 1
          adjusted[[name]] <- original[[name]] * adjustment
          print(paste("Adjusted", name, "by factor of", adjustment))
        } else {
          adjusted[[name]] <- original[[name]]
        }
        incProgress(1/total_rasters, detail = paste("Adjusting", name))
      }
      
      print("Raster adjustment complete.")
      return(list(original = original, adjusted = adjusted))
    })
  })
  
  # Process uploaded shapefile
  processShapefile <- reactive({
    req(input$shapefile)
    print("Processing uploaded shapefile...")
    
    temp_dir <- tempdir()
    zip::unzip(input$shapefile$datapath, exdir = temp_dir)
    shp_file <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)
    
    if (length(shp_file) == 0) {
      print("Error: No .shp file found in the uploaded zip.")
      return(NULL)
    }
    
    tryCatch({
      management_area <- st_read(shp_file)
      print("Shapefile read successfully.")
      
      # Get the CRS of the first raster (assuming all rasters have the same CRS)
      raster_crs <- crs(rasters()$tree_canopy_cover)
      
      # Reproject the management area to match raster CRS
      management_area <- st_transform(management_area, raster_crs)
      print(paste("Management area reprojected to:", raster_crs))
      
      print("Creating 1km2 grid...")
      bbox <- st_bbox(management_area)
      grid <- st_make_grid(management_area, cellsize = 1000, what = "centers")
      grid <- st_sf(geometry = grid)
      
      intersected_grid <- st_intersection(grid, management_area)
      print(paste("Grid created with", nrow(intersected_grid), "cells."))
      
      return(intersected_grid)
    }, error = function(e) {
      print(paste("Error processing shapefile:", e$message))
      return(NULL)
    })
  })
  
  # Generate covariates for grid cells
  generateGridCovariates <- reactive({
    req(processShapefile())
    req(adjustRasters())
    
    print("Generating covariates for grid cells...")
    grid <- processShapefile()
    if (is.null(grid)) {
      print("Error: Grid is NULL. Cannot generate covariates.")
      return(NULL)
    }
    
    raster_data <- adjustRasters()
    
    calculateCovariates <- function(point, rasters) {
      values <- lapply(rasters$adjusted, function(r) {
        extracted <- terra::extract(r, matrix(point, ncol = 2))
        print(paste("Extracted value for", names(r), ":", extracted))
        if (is.data.frame(extracted) && ncol(extracted) > 0) {
          return(as.numeric(extracted[1,1]))
        } else {
          return(NA)
        }
      })
      
      covariates <- list(
        canopy_gap_percent = mean(values$tree_canopy_cover < 10, na.rm = TRUE) * 100,
        open_forest_percent = mean(values$tree_canopy_cover >= 10 & values$tree_canopy_cover < 40, na.rm = TRUE) * 100,
        higher_severity_percent = mean(values$burn_severity >= 4, na.rm = TRUE) * 100,
        lower_severity_percent = mean(values$burn_severity >= 2 & values$burn_severity < 4, na.rm = TRUE) * 100,
        mean_years_since_wildfire = mean(values$years_since_wildfire, na.rm = TRUE),
        veg_departure_mean = mean(values$veg_departure, na.rm = TRUE),
        qmd_rmrs_mean = mean(values$qmd_rmrs, na.rm = TRUE),
        qmd_rmrs_cv = ifelse(mean(values$qmd_rmrs, na.rm = TRUE) > 0,
                             sd(values$qmd_rmrs, na.rm = TRUE) / mean(values$qmd_rmrs, na.rm = TRUE),
                             0),
        tpa_dead_mean = mean(values$tpa_dead, na.rm = TRUE),
        tpa_live_cv = ifelse(mean(values$tpa_live, na.rm = TRUE) > 0,
                             sd(values$tpa_live, na.rm = TRUE) / mean(values$tpa_live, na.rm = TRUE),
                             0)
      )
      
      # Convert any remaining NAs to 0
      covariates <- lapply(covariates, function(x) ifelse(is.na(x), 0, x))
      
      return(unlist(covariates))
    }
    
    print("Calculating covariates for each grid cell...")
    grid_points <- st_coordinates(grid)
    grid_covariates <- apply(grid_points, 1, calculateCovariates, rasters = raster_data)
    
    if (!is.numeric(grid_covariates)) {
      print("Error: Covariate calculation did not produce numeric results.")
      print("Structure of grid_covariates:")
      print(str(grid_covariates))
      return(NULL)
    }
    
    grid_covariates_df <- as.data.frame(t(grid_covariates))
    
    result <- cbind(grid, grid_covariates_df)
    
    print("Covariate generation complete.")
    print("Covariates for each grid cell:")
    print(grid_covariates_df)
    print(paste("Number of grid cells with covariates:", nrow(result)))
    return(result)
  })
  
  predictAbundance <- reactive({
    req(input$species)
    req(generateGridCovariates())
    req(mcmc_data())
    
    print("Predicting species abundance...")
    grid_covariates <- generateGridCovariates()
    if (is.null(grid_covariates)) {
      print("Error: Grid covariates are NULL. Cannot predict abundance.")
      return(NULL)
    }
    
    print("Available columns in grid_covariates:")
    print(names(grid_covariates))
    
    mcmc <- mcmc_data()
    
    species_index <- as.numeric(sub("Species ", "", input$species))
    print(paste("Predicting for Species", species_index))
    
    # Get the conversion factor for the selected species
    N2D <- species_data$N2D[species_index]
    
    # Extract species-specific coefficients
    int_psi <- mean(as.matrix(mcmc[[1]])[, paste0("int.psi[", species_index, "]")])
    int_lambda <- mean(as.matrix(mcmc[[1]])[, paste0("int.lambda[", species_index, "]")])
    beta_psi <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.psi\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
    beta_lambda <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.lambda\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
    
    print("Beta PSI:")
    print(beta_psi)
    print("Beta Lambda:")
    print(beta_lambda)
    
    # Explicitly define lambda covariates
    lambda_covariates <- c(
      "canopy_gap_percent", "open_forest_percent", "higher_severity_percent",
      "lower_severity_percent", "mean_years_since_wildfire", "veg_departure_mean",
      "qmd_rmrs_mean", "qmd_rmrs_cv", "tpa_dead_mean", "tpa_live_cv",
      "Douglas.Fir", "Engelmann.Spruce", "Western.Larch", "Grand.Fir",
      "Quaking.Aspen", "Rocky.Mountain.Maple", "F1", "F3", "F4", "A1", "A2", "A3",
      "higher_severity_percent:mean_years_since_wildfire",
      "lower_severity_percent:mean_years_since_wildfire"
    )
    
    print("Lambda covariates:")
    print(lambda_covariates)
    
    # Ensure we have the correct number of lambda covariates
    if (length(lambda_covariates) != length(beta_lambda)) {
      print("Warning: Mismatch between number of lambda covariates and beta coefficients")
      print(paste("Number of lambda covariates:", length(lambda_covariates)))
      print(paste("Number of beta coefficients:", length(beta_lambda)))
    }
    
    # Create a mapping of covariate names to their indices
    covariate_indices <- setNames(seq_along(lambda_covariates), lambda_covariates)
    
    # Include the updated predict_cell function here
    predict_cell <- function(cell) {
      tryCatch({
        # Extract long and lat from geometry
        coords <- st_coordinates(st_centroid(cell$geometry))
        long <- coords[1]
        lat <- coords[2]
        
        print(paste("Processing cell at coordinates:", long, lat))
        
        # Transform coordinates to WGS 1984 if necessary
        if (!identical(st_crs(cell), st_crs(4326))) {
          point <- st_sfc(st_point(c(long, lat)), crs = st_crs(cell))
          point_wgs84 <- st_transform(point, 4326)
          coords_wgs84 <- st_coordinates(point_wgs84)
          long <- coords_wgs84[1]
          lat <- coords_wgs84[2]
          print(paste("Transformed coordinates to WGS 1984:", long, lat))
        }
        
        # Prepare psi covariates
        psi_data <- c(
          (long - scale.psi$center[1]) / scale.psi$scale[1],
          (lat - scale.psi$center[2]) / scale.psi$scale[2]
        )
        psi_data <- c(psi_data, psi_data[1] * psi_data[2])  # lat:long interaction
        
        print("PSI data:")
        print(psi_data)
        
        # Calculate psi
        logit_psi <- int_psi + sum(psi_data * beta_psi)
        psi <- plogis(logit_psi)
        
        print(paste("Calculated PSI:", psi))
        
        # Prepare lambda covariates
        lambda_data <- numeric(length(lambda_covariates))
        names(lambda_data) <- lambda_covariates
        
        for (cov_name in lambda_covariates) {
          cov_index <- covariate_indices[cov_name]
          if (cov_name %in% names(cell)) {
            value <- as.numeric(cell[[cov_name]])
            if (is.na(value)) {
              print(paste("Warning: Covariate", cov_name, "is NA. Using 0."))
              lambda_data[cov_index] <- 0
            } else if (cov_index <= length(scale.lambda$center) && cov_index <= length(scale.lambda$scale)) {
              # Scale and center the lambda covariate
              lambda_data[cov_index] <- (value - scale.lambda$center[cov_index]) / scale.lambda$scale[cov_index]
            } else {
              print(paste("Warning: Scaling information for", cov_name, "not found. Using raw value."))
              lambda_data[cov_index] <- value
            }
          } else {
            print(paste("Warning: Covariate", cov_name, "not found in grid cell data. Using 0."))
            lambda_data[cov_index] <- 0
          }
        }
        
        # Handle interaction terms
        interaction_index1 <- covariate_indices["higher_severity_percent:mean_years_since_wildfire"]
        interaction_index2 <- covariate_indices["lower_severity_percent:mean_years_since_wildfire"]
        higher_severity_index <- covariate_indices["higher_severity_percent"]
        lower_severity_index <- covariate_indices["lower_severity_percent"]
        years_since_wildfire_index <- covariate_indices["mean_years_since_wildfire"]
        
        lambda_data[interaction_index1] <- lambda_data[higher_severity_index] * lambda_data[years_since_wildfire_index]
        lambda_data[interaction_index2] <- lambda_data[lower_severity_index] * lambda_data[years_since_wildfire_index]
        
        print("Lambda data:")
        print(lambda_data)
        
        # Check for NAs in lambda_data
        if (any(is.na(lambda_data))) {
          print("Warning: NAs found in lambda_data. This may cause issues in lambda calculation.")
          print("NA positions:")
          print(which(is.na(lambda_data)))
        }
        
        # Calculate lambda
        log_lambda <- int_lambda + sum(lambda_data * beta_lambda, na.rm = TRUE)
        lambda <- exp(log_lambda)
        
        print(paste("Calculated Lambda:", lambda))
        
        # Calculate abundance and convert to individuals/sq. km
        abundance <- psi * lambda * N2D
        
        print(paste("Calculated Abundance (individuals/sq. km):", abundance))
        
        return(c(psi = psi, lambda = lambda, abundance = abundance))
      }, error = function(e) {
        print(paste("Error in predict_cell:", e$message))
        return(c(psi = NA, lambda = NA, abundance = NA))
      })
    }
    
    # Apply prediction to each grid cell
    predictions <- do.call(rbind, lapply(1:nrow(grid_covariates), function(i) {
      predict_cell(grid_covariates[i,])
    }))
    
    # Add predictions to grid_covariates
    grid_covariates$psi <- predictions[, "psi"]
    grid_covariates$lambda <- predictions[, "lambda"]
    grid_covariates$abundance <- predictions[, "abundance"]
    
    print("Summary of predictions:")
    print(summary(predictions))
    
    # Check for any NA or infinite values in the predictions
    if (any(is.na(predictions) | is.infinite(predictions))) {
      print("Warning: NA or infinite values found in predictions.")
      print("Positions of problematic predictions:")
      print(which(is.na(predictions) | is.infinite(predictions), arr.ind = TRUE))
    }
    
    print("Abundance prediction complete.")
    print(paste("Mean predicted abundance (individuals/sq. km):", mean(grid_covariates$abundance, na.rm = TRUE)))
    
    # Make lambda_covariates available for the importance plot
    grid_covariates$lambda_covariates <- list(lambda_covariates)
    
    return(grid_covariates)
  })
  
  # Create abundance raster function
  create_abundance_raster <- function(abundance_data) {
    tryCatch({
      # Ensure abundance_data is an sf object
      if (!inherits(abundance_data, "sf")) {
        stop("abundance_data must be an sf object")
      }
      
      # Get the extent and CRS of the abundance data
      ext <- st_bbox(abundance_data)
      crs_raster <- st_crs(abundance_data)$proj4string
      
      # Set the resolution to 1km (or 1000 units in the coordinate reference system)
      res <- 1000
      
      # Create an empty raster with the same extent, CRS, and resolution
      r <- raster(xmn = ext["xmin"], xmx = ext["xmax"], ymn = ext["ymin"], ymx = ext["ymax"],
                  crs = crs_raster, resolution = res)
      
      # Rasterize the abundance data using the "abundance" field
      r <- rasterize(st_coordinates(abundance_data), r, field = abundance_data$abundance, fun = mean)
      
      if (is.null(r)) {
        print("Error: Failed to create a valid raster from abundance data.")
        return(NULL)
      }
      
      return(r)
    }, error = function(e) {
      print(paste("Error in create_abundance_raster:", e$message))
      return(NULL)
    })
  }
  
  # New reactive to calculate abundances for all species
  calculateAllAbundances <- reactive({
    req(generateGridCovariates())
    req(mcmc_data())
    
    print("Calculating abundances for all species...")
    grid_covariates <- generateGridCovariates()
    if (is.null(grid_covariates)) {
      print("Error: Grid covariates are NULL. Cannot calculate abundances.")
      return(NULL)
    }
    
    mcmc <- mcmc_data()
    
    all_species_abundances <- data.frame(species = character(), abundance = numeric())
    
    for (species_index in 1:nrow(species_data)) {
      print(paste("Calculating for Species", species_index))
      
      # Get the conversion factor for the selected species
      N2D <- species_data$N2D[species_index]
      
      # Extract species-specific coefficients
      int_psi <- mean(as.matrix(mcmc[[1]])[, paste0("int.psi[", species_index, "]")])
      int_lambda <- mean(as.matrix(mcmc[[1]])[, paste0("int.lambda[", species_index, "]")])
      beta_psi <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.psi\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
      beta_lambda <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.lambda\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
      
      # Calculate psi
      psi_data <- st_coordinates(st_centroid(grid_covariates))
      psi_data <- cbind(psi_data, psi_data[,1] * psi_data[,2])  # lat:long interaction
      psi_data <- scale(psi_data)
      logit_psi <- int_psi + as.matrix(psi_data) %*% beta_psi
      psi <- plogis(logit_psi)
      
      # Calculate lambda
      lambda_covariates <- c(
        "canopy_gap_percent", "open_forest_percent", "higher_severity_percent",
        "lower_severity_percent", "mean_years_since_wildfire", "veg_departure_mean",
        "qmd_rmrs_mean", "qmd_rmrs_cv", "tpa_dead_mean", "tpa_live_cv"
      )
      lambda_data <- as.matrix(grid_covariates[, lambda_covariates])
      lambda_data <- scale(lambda_data)
      log_lambda <- int_lambda + lambda_data %*% beta_lambda
      lambda <- exp(log_lambda)
      
      # Calculate abundance and convert to individuals/sq. km
      abundance <- psi * lambda * N2D
      
      # Sum abundance for this species
      total_abundance <- sum(abundance, na.rm = TRUE)
      
      # Add to all_species_abundances
      all_species_abundances <- rbind(all_species_abundances, 
                                      data.frame(species = species_data$common.name[species_index],
                                                 abundance = total_abundance))
    }
    
    print("Abundance calculation complete for all species.")
    return(all_species_abundances)
  })
  
  # Render the species abundance comparison plot
  # Function to create the plot
  create_species_abundance_plot <- function() {
    # Create the data frame
    species_data <- data.frame(
      species = c("Broad-tailed Hummingbird", "Cassin's Finch", "Calliope Hummingbird",
                  "Cassin's Vireo", "Chipping Sparrow", "Clark's Nutcracker",
                  "Dusky Flycatcher", "Dusky Grouse", "Evening Grosbeak",
                  "Golden-crowned Kinglet", "Hammond's Flycatcher", "Lazuli Bunting",
                  "Least Flycatcher", "Lewis's Woodpecker", "MacGillivray's Warbler",
                  "Mountain Chickadee", "Olive-sided Flycatcher", "Pinyon Jay",
                  "Pine Siskin", "Plumbeus Vireo", "Red-breasted Nuthatch",
                  "Red Crossbill", "Red-naped Sapsucker", "Rufous Hummingbird",
                  "Townsend's Solitaire", "Townsend's Warbler", "Vaux's Swift",
                  "Varied Thrush", "Veery", "Warbling Vireo", "Western Tanager",
                  "Western Wood-Pewee", "Willow Flycatcher", "Williamson's Sapsucker",
                  "Wilson's Warbler"),
      abundance = c(1, 27, 15, 36, 136, 9, 272, 1, 105, 110, 160, 14, 1, 0, 51, 220,
                    9, 119, 0, 0, 72, 88, 6, 42, 58, 115, 3.77, 4, 1, 43, 93, 11, 2,
                    0, 19)
    )
    
    # Sort the data by abundance in descending order
    species_data <- species_data %>%
      arrange(desc(abundance)) %>%
      mutate(species = factor(species, levels = species))
    
    # Create the plot
    ggplot(species_data, aes(x = reorder(species, abundance), y = abundance, fill = species == species[7])) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("FALSE" = "#8884d8", "TRUE" = "#ff7300")) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        title = "Species Abundance Comparison",
        x = "Abundance"
      ) +
      coord_flip()  # Flip coordinates for horizontal bars
  }
  
  # In your server function:
  output$speciesAbundancePlot <- renderPlot({
    create_species_abundance_plot()
  })  # Adjust height and width as needed
  
  # Download handler for the plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "species_abundance_comparison.png"
    },
    content = function(file) {
      ggsave(file, plot = create_species_abundance_plot(), width = 12, height = 8, dpi = 300)
    }
  )
  
  # Render abundance map
  # Render abundance map
  output$abundanceMap <- renderLeaflet({
    req(predictAbundance())
    print("Rendering abundance map...")
    
    abundance_data <- predictAbundance()
    if (is.null(abundance_data)) {
      print("Error: Abundance data is NULL. Cannot render map.")
      return(NULL)
    }
    
    # Create a raster from the abundance data
    abundance_raster <- create_abundance_raster(abundance_data)
    
    if (is.null(abundance_raster)) {
      print("Error: Failed to create abundance raster. Cannot render map.")
      return(NULL)
    }
    
    # Get the range of abundance values
    abundance_range <- range(values(abundance_raster), na.rm = TRUE)
    
    # Create a color palette with values sorted from highest to lowest
    pal <- colorNumeric(
      palette = rev(colorRampPalette(c("yellow", "red"))(100)),
      domain = abundance_range,
      na.color = "transparent"
    )
    
    # Get the species name
    species_name <- species_data$common.name[as.numeric(sub("Species ", "", input$species))]
    
    # Create the leaflet map
    leaflet() %>%
      addTiles() %>%
      addRasterImage(abundance_raster, colors = pal, opacity = 0.7) %>%
      addLegend(
        pal = pal, 
        values = abundance_range,
        title = paste(species_name, "<br>(individuals/sq. km)"),
        labFormat = labelFormat(
          transform = function(x) sort(x, decreasing = TRUE)
        )
      )
  })
  
  # Render covariate importance plot
  output$covariateImportance <- renderPlot({
    req(predictAbundance())
    req(mcmc_data())
    print("Generating covariate importance plot...")
    
    abundance_data <- predictAbundance()
    species_index <- as.numeric(sub("Species ", "", input$species))
    mcmc <- mcmc_data()
    
    beta_lambda <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.lambda\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
    
    importance <- abs(beta_lambda)
    names(importance) <- abundance_data$lambda_covariates[[1]]  # Use the lambda_covariates from the abundance data
    
    par(mar = c(10, 4, 4, 2) + 0.1)  # Increase bottom margin for long labels
    barplot(sort(importance, decreasing = TRUE),
            main = "Covariate Importance",
            las = 2,
            cex.names = 0.7,  # Reduce text size if needed
            ylab = "Absolute coefficient value")
    
    print("Covariate importance plot generated.")
  })
  
  # Render abundance statistics table
  # Update the abundance statistics table
  output$abundanceStats <- renderTable({
    req(predictAbundance())
    print("Calculating abundance statistics...")
    
    abundance_data <- predictAbundance()
    if (is.null(abundance_data)) {
      print("Error: Abundance data is NULL. Cannot calculate statistics.")
      return(NULL)
    }
    
    # Calculate the sum of abundances across all grid cells
    total_abundance <- sum(abundance_data$abundance, na.rm = TRUE)
    
    stats <- data.frame(
      Statistic = c("Total Abundance", "Mean Abundance per Cell", "Number of Cells"),
      Value = c(
        total_abundance,
        mean(abundance_data$abundance, na.rm = TRUE),
        nrow(abundance_data)
      )
    )
    
    stats$Value <- round(stats$Value, 2)  # Round to 2 decimal places for readability
    
    print("Abundance statistics calculated.")
    return(stats)
  })
  
  # Download handler for abundance data
  output$downloadAbundance <- downloadHandler(
    filename = function() {
      paste("abundance_data_", species_data$common.name[as.numeric(sub("Species ", "", input$species))], ".csv", sep = "")
    },
    content = function(file) {
      abundance_data <- predictAbundance()
      if (!is.null(abundance_data)) {
        write.csv(abundance_data, file, row.names = FALSE)
      }
    }
  )
  
  # Download handler for the abundance map
  output$downloadMap <- downloadHandler(
    filename = function() {
      paste("abundance_map_", species_data$common.name[as.numeric(sub("Species ", "", input$species))], ".png", sep = "")
    },
    content = function(file) {
      req(predictAbundance())
      abundance_data <- predictAbundance()
      
      if (!is.null(abundance_data)) {
        abundance_raster <- create_abundance_raster(abundance_data)
        
        if (!is.null(abundance_raster)) {
          png(file, width = 800, height = 600)
          plot(abundance_raster, main = paste("Abundance Map for", species_data$common.name[as.numeric(sub("Species ", "", input$species))]))
          dev.off()
        } else {
          stop("Failed to create abundance raster for download.")
        }
      } else {
        stop("No abundance data available for download.")
      }
    }
  )
  
  # Download handler for abundance statistics
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste("abundance_stats_", species_data$common.name[as.numeric(sub("Species ", "", input$species))], ".csv", sep = "")
    },
    content = function(file) {
      req(predictAbundance())
      abundance_data <- predictAbundance()
      
      if (!is.null(abundance_data)) {
        stats <- data.frame(
          Statistic = c("Mean", "Median", "Min", "Max", "Standard Deviation"),
          Value = c(
            mean(abundance_data$abundance, na.rm = TRUE),
            median(abundance_data$abundance, na.rm = TRUE),
            min(abundance_data$abundance, na.rm = TRUE),
            max(abundance_data$abundance, na.rm = TRUE),
            sd(abundance_data$abundance, na.rm = TRUE)
          )
        )
        stats$Value <- round(stats$Value, 6)
        write.csv(stats, file, row.names = FALSE)
      } else {
        stop("No abundance statistics available for download.")
      }
    }
  )
  
  # Download handler for covariate importance plot
  # output$downloadPlot <- downloadHandler(
  #   filename = function() {
  #     paste("covariate_importance_", species_data$common.name[as.numeric(sub("Species ", "", input$species))], ".png", sep = "")
  #   },
  #   content = function(file) {
  #     req(predictAbundance(), mcmc_data())
  #     abundance_data <- predictAbundance()
  #     mcmc <- mcmc_data()
  #     species_index <- as.numeric(sub("Species ", "", input$species))
  #     
  #     if (!is.null(abundance_data) && !is.null(mcmc)) {
  #       beta_lambda <- colMeans(as.matrix(mcmc[[1]])[, grep(paste0("beta.lambda\\[", species_index, ","), colnames(as.matrix(mcmc[[1]])))])
  #       importance <- abs(beta_lambda)
  #       names(importance) <- abundance_data$lambda_covariates[[1]]
  #       
  #       png(file, width = 800, height = 600)
  #       par(mar = c(10, 4, 4, 2) + 0.1)
  #       barplot(sort(importance, decreasing = TRUE),
  #               main = paste("Covariate Importance for", species_data$common.name[species_index]),
  #               las = 2,
  #               cex.names = 0.7,
  #               ylab = "Absolute coefficient value")
  #       dev.off()
  #     } else {
  #       stop("No covariate importance data available for download.")
  #     }
  #   }
  # )
}

shinyApp(ui, server)

# create my own dataframe
