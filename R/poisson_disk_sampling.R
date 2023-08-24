# Fallback implementation of Poisson Disk Sampling (Bridson's algorithm).

.create_grid <- function(.polygon, .min_dist, .planar_crs = NA) {
  
  # default to Web Mercator:
  if (is.na(.planar_crs)) {
    .planar_crs <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  }
  
  polygon_planar <- sf::st_transform(.polygon, crs = .planar_crs)
  
  # Create grid:
  cell_size <- .min_dist / sqrt(2)
  grid <- sf::st_make_grid(polygon_planar, cellsize = rep(cell_size, 2))
  grid <- sf::st_sf(grid)
  grid$id <- seq_len(nrow(grid))
  grid$contains_point <- FALSE
  
  grid[, c("id", "contains_point", "grid")]
}

# Finding cell neighbors:

.find_neighbors <- function(.id, .grid) {
  # Find direct neighbors:
  cell <- .grid[.grid$id == .id, ]
  direct_neighbors <- sf::st_intersects(cell, .grid)[[1]]

  # Find indirect neighbors:
  direct_neighbors |>
    lapply(function(x) sf::st_intersects(.grid[.grid$id == x, ], .grid)[[1]]) |>
    do.call(c, args = _) |>
    unique() |>
    (function(v) v[v != .id])()
}

.make_bagel <- function(.point, .min_dist) {
  double_radius <- sf::st_buffer(.point, .min_dist * 2)
  radius <- sf::st_buffer(.point, .min_dist)
  
  sf::st_difference(double_radius, radius)
}

.new_counter <- function() {
  i <- 0
  
  function() {
    i <<- i + 1
    i
  }
}

#' Generate quasi-random points using Poisson Disk Sampling
#'
#' Generate quasi-random points inside the bounding box of an `sf` POLYGON or
#' MULTIPOLYGON object. This function is heavily inspired by Bridson's algorithm
#' for fast n-dimensional Poisson Disk Sampling. 
#' 
#' @param polygon Polygon to use as base for sampling.
#' @param min_dist Minimum distance between points (in units of the object's crs; 
#' check using `sf::st_crs()`)
#' @param planar_crs Intermediate crs to use for point generation; defaults to 
#' Web Mercator.
#' @param k Number of points to sample (i.e. try to place) for each active point/
#' iteration; defaults to 30.
#' @param cores Number of CPU cores/threads to allocate for the parallelizable 
#' sequences. By default, parallelization occurs across all available cores.
#' @return A `sf`-object of type POINT.
#' @examples
#' germany <- rnaturalearth::ne_countries(country = "Germany", returnclass = "sf")
#' test <- poisson_disk_sampling(germany, min_dist = 100000)
#'
#' germany |> 
#'   ggplot() +
#'   geom_sf() +
#'   geom_sf(test, mapping = aes())
#' @export
poisson_disk_sampling <- function(polygon, min_dist, planar_crs = NA, k = 30, cores = parallel::detectCores()) {
  origin_crs <- sf::st_crs(polygon)
  
  future::plan(future::multisession, workers = cores)
  
  # Setting up:
  # Generating grid structure & a lookup for which cells are "neighbors"
  # so we don't have to check separately for every point later on:
  message("Generating grid...")
  grid <- .create_grid(.polygon = polygon, .min_dist = min_dist, .planar_crs = planar_crs)
  box <- sf::st_union(grid)
  
  message("Generating lookup of neighboring cells...")
  neighbors <- purrr::map(grid$id, function(x) .find_neighbors(x, grid), .progress = TRUE)
  
  message("Entering routine...\n")
  
  # here: create active list, starting with one random point:
  active_list <- sf::st_sample(grid, 1)
  point <- active_list # Initialize collection of valid points with active list
  valid_points <- sf::st_sf(point)
  
  # Keep track of which grid cells already contain points:
  containing_cell <- sf::st_intersects(active_list, grid)[[1]]
  valid_points$grid_id <- containing_cell
  
  grid$contains_point[grid$id == containing_cell] <- TRUE
  
  counter <- .new_counter()
  
  suppressWarnings({
    while (length(active_list) > 0) {
      index <- sample(seq_along(active_list), 1)
      active_point <- active_list[index]
      
      # Generate candidates:
      bagel <- .make_bagel(active_point, min_dist)
      bagel <- sf::st_intersection(bagel, box) # in case the bagel exceeds the grid
      candidates <- sf::st_sample(bagel, k)
      
      candidates <- sf::st_sf(candidates)
      candidates$grid_id <- NA_real_
      
      for (i in seq_len(nrow(candidates))) {
        
        containing_cell <- sf::st_intersects(candidates[i, ], grid)[[1]]
        
        if (grid[containing_cell, ]$contains_point) {
          next
        }
        
        neighborhood <- neighbors[[containing_cell]]
        
        if (any(grid[neighborhood, ]$contains_point)) {
          neighboring_pts <- valid_points[valid_points$grid_id %in% neighborhood, ]
          
          if (any(as.numeric(sf::st_distance(candidates[i, ], neighboring_pts)) < min_dist)) {
            next
          }
        }
        
        # It should already be appended here to circumvent this:
        
        grid$contains_point[grid$id == containing_cell] <- TRUE
        candidates$grid_id[i] <- containing_cell
        
        to_append <- candidates[i, ]
        names(to_append)[1] <- "point"
        sf::st_geometry(to_append) <- "point"
        
        valid_points <- rbind(valid_points, to_append)
      }
      
      valid_candidates <- candidates[!is.na(candidates$grid_id), ]
      names(valid_candidates)[1] <- "point"
      sf::st_geometry(valid_candidates) <- "point"
      
      # ...maybe put all of this ^ into a function since we are not modifying
      # existing data structures until later...?
      
      # Drop active point from active list:
      active_list <- c(active_list, valid_candidates$point)
      active_list <- active_list[-index]
      # # Append valid candidates to valid_points:
      # valid_points <- rbind(valid_points, valid_candidates)
      
      cat(
        "Iteration no. ", counter(), " complete.\n",
        "Valid candidates: ", nrow(valid_candidates), " points.\n",
        "Active List: ", length(active_list), " points.\n\n",
        "Total number of points generated: ", nrow(valid_points), "\n\n\n",
        sep = ""
      )
    }
  })
  
  valid_points
}
