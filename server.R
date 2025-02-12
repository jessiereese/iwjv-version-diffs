library(shiny)
library(leaflet)
library(raster)
library(sf)
library(dplyr)
library(coda)
library(zip)
library(terra)

# Jessie note: this is the "beta" version referred to in the readme

mcmc_data <- load("IWJV/data/test_model.R")

server <- function(input, output, session) {
  
  # Load MCMC data
  mcmc_data <- reactive({
    print("Loading MCMC data...")
    load("IWJV/data/test_model.R")
    print("MCMC data loaded successfully.")
    return(mcmc.lst)
  })
  
  # Load raster data                              ## Jessie note: same as stable
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
  adjustRasters <- reactive({            ## Jessie note: updated compared to stable. Removed withProgress updater, variable names changed for output.
    req(input$runModel)
    print("Adjusting rasters based on user input...")
    
    original <- rasters()
    adjusted <- list()
    
    for (name in names(original)) {
      if (name %in% c("tree_canopy_cover", "qmd_rmrs", "tpa_dead", "tpa_live")) {
        adjustment <- input[[name]] / 100 + 1
        adjusted[[name]] <- original[[name]] * adjustment
        print(paste("Adjusted", name, "by factor of", adjustment))
      } else {
        adjusted[[name]] <- original[[name]]
      }
    }
    
    print("Raster adjustment complete.")
    return(list(original = original, adjusted = adjusted))
  })
  
  # Create the 1 km² grid polygons (1000m x 1000m)
  createGrid <- function(management_area) {   ## Jessie note: updated to pull this function out separately from generateing the processShapefile reactive object
    # Generate grid polygons (1km x 1km)
    grid <- st_make_grid(management_area, cellsize = 1000, what = "polygons")
    
    # Convert to sf object
    grid <- st_sf(geometry = grid)
    
    # Get the CRS of the rasters (assuming all rasters have the same CRS)
    raster_crs <- crs(rasters()$tree_canopy_cover)
    
    # Reproject the grid to match the raster CRS
    grid <- st_transform(grid, raster_crs)
    
    print(paste("1 km² grid created with", nrow(grid), "polygons."))
    return(grid)
  }
  
  # Process uploaded shapefile
  # Create the grid and calculate covariates for each grid cell
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
      
      # Create the 1km² grid using the management area
      grid <- createGrid(management_area)
      
      return(grid)
    }, error = function(e) {
      print(paste("Error processing shapefile:", e$message))
      return(NULL)
    })
  })
  
  generateGridCovariates <- reactive({
    req(processShapefile())
    req(adjustRasters())
    
    print("Generating covariates for grid cells...")
    
    # Get the grid of polygons
    grid <- processShapefile()
    
    # Get the CRS of the grid (assuming it comes from your management area shapefile)
    grid_crs <- st_crs(grid)$proj4string
    
    # Adjust and reproject rasters if necessary
    reprojectRaster <- function(raster, grid_crs) {
      # Get the CRS of the raster
      raster_crs <- terra::crs(raster)
      if (raster_crs != grid_crs) {
        # Reproject raster to match the grid's CRS
        raster <- terra::project(raster, grid_crs)
        print(paste("Reprojected raster to match grid CRS:", grid_crs))
      }
      return(raster)
    }
    
    # Reproject all rasters to match the grid CRS
    raster_data <- adjustRasters()
    raster_data$adjusted <- lapply(raster_data$adjusted, reprojectRaster, grid_crs)
    
    # For each raster, calculate the mean, CV, and other statistics for each polygon
    calculateCovariates <- function(r, grid) {
      if (is.null(r)) {
        stop("Raster is NULL. Please check if the raster was loaded correctly.")
      }
      
      # Extract raster values for the grid polygons and compute statistics (e.g., mean)
      mean_values <- terra::extract(r, grid, fun = mean, na.rm = TRUE, df = TRUE)
      sd_values <- terra::extract(r, grid, fun = sd, na.rm = TRUE, df = TRUE)
      
      # Calculate CV (coefficient of variation) if mean > 0
      cv_values <- ifelse(mean_values[,2] > 0, (sd_values[,2] / mean_values[,2]), 0)
      
      # Return a data frame of statistics
      return(data.frame(
        mean = mean_values[,2],
        cv = cv_values
      ))
    }
    
    # Apply covariate calculation for each raster
    tree_canopy_covariates <- calculateCovariates(raster_data$adjusted$tree_canopy_cover, grid)
    qmd_rmrs_covariates <- calculateCovariates(raster_data$adjusted$qmd_rmrs, grid)
    tpa_dead_covariates <- calculateCovariates(raster_data$adjusted$tpa_dead, grid)
    tpa_live_covariates <- calculateCovariates(raster_data$adjusted$tpa_live, grid)
    burn_severity_covariates <- calculateCovariates(raster_data$adjusted$burn_severity, grid)
    veg_departure_covariates <- calculateCovariates(raster_data$adjusted$veg_departure, grid)
    years_since_wildfire_covariates <- calculateCovariates(raster_data$adjusted$years_since_wildfire, grid)
    
    # Combine the covariates into a final data frame
    grid_covariates_df <- data.frame(
      canopy_mean = tree_canopy_covariates$mean,
      canopy_cv = tree_canopy_covariates$cv,
      qmd_rmrs_mean = qmd_rmrs_covariates$mean,
      qmd_rmrs_cv = qmd_rmrs_covariates$cv,
      tpa_dead_mean = tpa_dead_covariates$mean,
      tpa_live_cv = tpa_live_covariates$cv,
      burn_severity_mean = burn_severity_covariates$mean,
      veg_departure_mean = veg_departure_covariates$mean,
      mean_years_since_wildfire = years_since_wildfire_covariates$mean
    )
    
    # Combine the grid geometry with covariates
    result <- cbind(grid, grid_covariates_df)
    
    print("Covariate generation complete.")
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
        
        # Calculate abundance
        abundance <- psi * lambda
        
        print(paste("Calculated Abundance:", abundance))
        
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
    print(paste("Mean predicted abundance:", mean(grid_covariates$abundance, na.rm = TRUE)))
    
    # Make lambda_covariates available for the importance plot
    grid_covariates$lambda_covariates <- list(lambda_covariates)
    
    return(grid_covariates)
  })
  
  # Update the create_abundance_raster function with error checking
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
  
  # Update the renderLeaflet function
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
    
    # Create a color palette
    pal <- colorNumeric(c("yellow", "red"), values(abundance_raster), na.color = "transparent")
    
    # Create the leaflet map
    leaflet() %>%
      addTiles() %>%
      addRasterImage(abundance_raster, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(abundance_raster), title = "Abundance")
  })
  
  # Update the covariate importance plot
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
  
  output$abundanceStats <- renderTable({
    req(predictAbundance())
    print("Calculating abundance statistics...")
    
    abundance_data <- predictAbundance()
    if (is.null(abundance_data)) {
      print("Error: Abundance data is NULL. Cannot calculate statistics.")
      return(NULL)
    }
    
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
    
    stats$Value <- round(stats$Value, 6)  # Round to 6 decimal places for readability
    
    print("Abundance statistics calculated.")
    return(stats)
  })
  
  # Download handler for abundance data
  output$downloadAbundance <- downloadHandler(
    filename = function() {
      paste("abundance_data_species_", input$species, ".csv", sep = "")
    },
    content = function(file) {
      abundance_data <- predictAbundance()
      if (!is.null(abundance_data)) {
        write.csv(abundance_data, file, row.names = FALSE)
      }
    }
  )
}

# shinyApp(ui, server)
