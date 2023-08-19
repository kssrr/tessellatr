# `tessellatr` - Generate quasi-random points or tiles inside `sf`-POLYGONS

Simple wrappers for `spatstat.random::rSSI()` to generate quasi-random points or tiles inside `sf`-POLYGON or MULTIPOLYGON objects. Sometimes, regular random sampling is not enough, and you may want a more uniform sample distribution, where no two points are too close together:

<p align="center"><img src="https://github.com/kssrr/tessellatr/assets/121236725/b880c106-2d45-41bf-bc1a-7bad76a16b47" alt="Random vs. Quasi-Random" width="400"></p>

There is a convenient implementation in `spatstat.random::rSSI()`, but coercion between `sf` and `spatstat`-objects is messy, which is what these wrappers handle for you. 

## Installation

The package is not on CRAN, but you can install it from Github using

```
devtools::install_github("kssrr/tessellatr")
```

Alternatively, you can just copy the functions from `R/default_generators.R` into your session.

## Example usage

First, we need a POLYGON or MULTIPOLYGON object. Both `st_random_points()` and `st_random_voronoi()` require you to specify a minimum distance between the points for sampling; the unit here is the unit of the CRS of the POLYGON (for example, meters for WGS84):

```{r}
library(ggplot2)

germany <- rnaturalearth::ne_countries(country = "Germany", returnclass = "sf")
pts <- tessellatr::st_random_points(germany, .min_dist = 20000) # distance in m. for this CRS

germany |> 
  ggplot() +
  geom_sf() +
  geom_sf(pts, mapping = aes()) +
  theme_void()
```

<p align="center"><img src="https://github.com/kssrr/tessellatr/assets/121236725/71a1d871-8380-4a80-904e-71dca130fcb2" alt="Random Points" width="400"></p>

```{r}
vor <- tessellatr::st_random_voronoi(germany, 20000)

vor |> 
  ggplot() +
  geom_sf() +
  theme_void()
```

<p align="center"><img src="https://github.com/kssrr/tessellatr/assets/121236725/2bca6da5-2cdd-4c56-b8cc-e679876dea02" alt="Random Voronoi Polygons" width="400"></p>
