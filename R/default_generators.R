# Default generators

#' Generate random points using Simple Sequential Inhibition
#' 
#' Wrapper for `spatstat.random::rSSI()` to generate quasi-random points within an `sf`-
#' POLYGON or MULTIPOLYGON (similar to `sf::st_sample()`). For some MULTIPOLYGONS,
#' you may need to turn spherical geometry off (`sf::sf_use_s2(FALSE)`)
#' 
#' @param polygon Boundary polygon
#' @param .n How many points to place (defaults to `Inf`, meaning as many points as
#' possible will be fit)
#' @param .min_dist Minimum distance between points (in units of `polygon`'s CRS). 
#' @param .planar_crs Planar CRS to use for coercion to `spatstat`-object; defaults
#' to Web Mercator.
#' @param ... Additional arguments to pass to `spatstat.random::rSSI()`
#' @return A `sf`-object of type POINT.
#' @examples
#' library(ggplot2)
#'
#' germany <- rnaturalearth::ne_countries(country = "Germany", returnclass = "sf")["admin"]
#' germany_pts <- st_random_points(germany, .min_dist = 40000)
#'
#' germany |> 
#'   ggplot() +
#'   geom_sf() +
#'   geom_sf(germany_pts, mapping = aes()) +
#'   theme_void()
#' @export
st_random_points <- function(polygon, .n = Inf, .min_dist, .planar_crs = NA, ...) {
  origin_crs <- sf::st_crs(polygon)
  web_mercator <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  if (is.na(.planar_crs)) {
    .planar_crs <- web_mercator
  }
  
  # Transform to planar:
  polygon_planar <- 
    polygon |> 
    sf::st_union() |> 
    sf::st_transform(.planar_crs) |> 
    spatstat.geom::as.owin()
  
  # Random points:
  pts <- spatstat.random::rSSI(
    r = .min_dist,
    n = .n,
    win = polygon_planar,
    ...
  )
  
  # Coerce back into `sf`-obj:
  pts_sf <- sf::st_as_sf(pts)
  pts_sf <- pts_sf[pts_sf$label != "window", ]
  sf::st_crs(pts_sf) <- .planar_crs
  
  # Transform back to origin CRS:
  sf::st_transform(pts_sf, crs = origin_crs)
}

#' Random voronoi-tiles within polygon
#'
#' Generate quasi-random voronoi-tiles within a `sf` POLYGON or MULTIPOLYGON. For 
#' some MULTIPOLYGONS, you may need to turn spherical geometry off 
#' (`sf::sf_use_s2(FALSE)`)
#'
#' @param polygon Boundary polygon (a POLYGON or MULTIPOLYGON)
#' @param min_dist Minimum distance between tiles' centroids (in units of `polygon`'s CRS)
#' @param n Number of tiles to generate (defaults to `Inf`, meaning as many as can be
#' fit without violating `minimum_distance` will be fit).
#' @param planar_crs Planar projection to use for generating centroids; defaults to 
#' Web Mercator.
#' @param ... Additional arguments to pass to `spatstat.random::rSSI()`, which handles the 
#' quasi-random generation of centroids for the tiles.
#' @return A `sf`-object of type MULTIPOLYGON.
#' @examples
#' library(ggplot2)
#'
#' germany <- rnaturalearth::ne_countries(country = "Germany", returnclass = "sf")["admin"]
#'
#' germany_vor <- st_random_voronoi(germany, min_dist = 40000)
#' ggplot(germany_vor) +
#'   geom_sf() +
#'   theme_void()
#' @export
st_random_voronoi <- function(polygon, min_dist, n = Inf, planar_crs = NA, ...) {
  
  message("Generating centroids...\n")
  pts <- st_random_points(polygon, .n = n, .min_dist = min_dist, .planar_crs = planar_crs, ...)
  
  message("Tessellating...\n")
  geom <- 
    pts |> 
    sf::st_union() |> 
    sf::st_voronoi() |> 
    sf::st_collection_extract() |> 
    sf::st_intersection(polygon)
  
  geom <- sf::st_sf(geom)
  geom[["id"]] <- seq_len(nrow(geom))
  
  geom[, c("id", "geom")]
}

germany <- rnaturalearth::ne_countries(country = "Germany", returnclass = "sf")

germany  |> 
  sf::st_crs()  |> 
  unclass()  |> 
  _$wkt  |> 
  stringr::str_extract('LENGTHUNIT\\[(.*)\\]')
