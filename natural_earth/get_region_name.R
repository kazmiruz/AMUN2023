### Location Identification ###

# install.packages("rnaturalearthdata") # Install if not already installed
# install.packages("sf") # Install if not already installed
library(rnaturalearth)
library(sf)
library(dplyr)
library(purrr)

# Load and transform datasets to WGS84
marine_polys <- st_read("natural_earth/ne_10m_geography_marine_polys/ne_10m_geography_marine_polys.shp") %>% st_transform(crs = 4326) %>%
  mutate(name = ifelse(str_detect(name, "River"), paste("The", name), name))
rivers <- st_read("natural_earth/ne_10m_rivers_lake_centerlines_scale_rank/ne_10m_rivers_lake_centerlines_scale_rank.shp") %>%
  st_transform(crs = 4326)
provinces <- ne_states(returnclass = "sf") %>% st_transform(crs = 4326)

# Function to identify regions with buffers and default to "Canada" if no match is found
get_region_name <- function(latitude, longitude) {
  
  # Create a point in WGS84 for the given latitude and longitude
  point <- st_sfc(st_point(c(longitude, latitude)), crs = 4326) %>% st_sf(geometry = .)
  
  # Check for marine region with 2 km buffer
  marine_region <- st_join(st_buffer(st_transform(point, 32618), dist = 2000) %>% st_transform(4326), marine_polys, join = st_intersects)
  if (nrow(marine_region) > 0 && !is.na(marine_region$name[1])) {
    return(paste("Canada:", marine_region$name[1]))
  }

  # Check for river with 0.5 km buffer
  river_match <- st_join(st_buffer(st_transform(point, 32618), dist = 500) %>% st_transform(4326), rivers, join = st_intersects)
  if (nrow(river_match) > 0 && !is.na(river_match$name[1])) {
    return(paste("Canada: The", river_match$name[1]))
  }

  # Check for province with 1 km buffer
  province_region <- st_join(st_buffer(st_transform(point, 32618), dist = 1000) %>% st_transform(4326), provinces, join = st_intersects)
  if (nrow(province_region) > 0 && !is.na(province_region$name[1])) {
    return(paste("Canada:", province_region$name[1]))
  }
  
  # Default to "Canada" if no match
  return("Canada")
  
}

sf_use_s2(FALSE) # Disable S2 Geometry Processing
