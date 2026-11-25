# NBN species forecasting 2026: extra cell covariates
# Lucy Somekh

# Builds a table with one row per 5 km grid cell in the two target regions, describing where recorders are likely to go and the habitat each species needs
# Output: extra_data/cell_covariates.csv
#   protected_share: share of the cell that is a nature reserve or SSSI
#   population: number of people living in the cell
#   woodland_pct, heath_pct, freshwater_pct, coastal_pct, saltwater_pct, urban_pct: % of the cell covered by each type of land cover
#   river_km: length of river in the cell

# Why: the model forecasts where a species will be RECORDED, so testing some assumptions here:
# Reserves & SSSIs attrac surveys, more people nearby means more casual sightings, & habitat says where a species could exist in a cell with no records yet.
# These were all tested in 03_forecast_pipeline.R but only river_km (for otter) passed the 6-of-8 rule, so it's the only one used

library(sf)
library(terra)

#### Settings ####
# Raw downloads from websites, saved in extra_data/raw/

# WorldPop UK population 2020, 1 km, UN-adjusted (https://hub.worldpop.org/geodata/summary?id=39457)
# Natural England Open Data: SSSIs, Local Nature Reserves and National Nature Reserves (England)
# NatureScot Open Data: SSSIs (Scotland)
# UKCEH Land Cover Map 2025, 1 km percentage target class (https://catalogue.ceh.ac.uk/documents/c71ad6f1-200c-4fdf-b946-a2e3923794ee)
# OS Open Rivers (https://osdatahub.os.uk/data/downloads/open/OpenRivers)
raw_dir  <- "extra_data/raw"
out_file <- "extra_data/cell_covariates.csv"

protected_files <- c("SSSI_England_8759627969800523937.gpkg",
                     "Local_Nature_Reserves_England_2264923962887655393.gpkg",
                     "National_Nature_Reserves_England_-3260325413303632251.gpkg",
                     "Sites_of_Special_Scientific_Interest_8447690424398606363.gpkg")  # Scotland
population_file <- "gbr_ppp_2020_1km_Aggregated_UNadj.tif"
landcover_file  <- "gb2025lcm1km_percentage_target.tif"
rivers_file     <- "oprvrs_gb.gpkg"

# Find a downloaded file anywhere inside extra_data/raw/
find_file <- function(name) {
  pattern <- paste0("^", gsub(".", "\\.", name, fixed = TRUE), "$")   # exact file name
  list.files(raw_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)[1]
}


#### Rebuild the 5 km grid cells as squares ####
# Only given the centre of each cell, not its outline. Since the grid lines up with the British National Grid (EPSG:27700) I worked out the pattern from the cell coordinates:
# column = (grid_id + 131) %% 132,  row = (grid_id + 131) %/% 132
# centre easting = 2384 + 5000 * column,  centre northing = 2853 + 5000 * row
cells <- read.csv("data/env_predictors/gridded/gridden_env_2025.csv")[, c("grid_id", "Region")]
cells <- cells[cells$Region %in% c("South East", "Scotland"), ]

cx <- 2384 + 5000 * ((cells$grid_id + 131) %%  132)
cy <- 2853 + 5000 * ((cells$grid_id + 131) %/% 132)

# 5 km square around each centre
squares <- lapply(seq_along(cx), function(i) {
  st_polygon(list(cbind(cx[i] + c(-2500, 2500, 2500, -2500, -2500),
                        cy[i] + c(-2500, -2500, 2500, 2500, -2500))))
})
cells_sf <- st_sf(cells, geometry = st_sfc(squares, crs = 27700))


#### Protected areas: share of each cell covered ####
# Put everything in British National Grid so areas are in square metres
# st_make_valid repairs small errors in some downloaded shapes.
protected <- lapply(protected_files, function(f) {
  st_geometry(st_make_valid(st_transform(st_read(find_file(f), quiet = TRUE), 27700)))
})

# Sites overlap (e.g. a National Nature Reserve inside an SSSI), so merge
# them into one shape first, otherwise the same ground is counted twice
protected <- st_union(do.call(c, protected))

pieces  <- st_intersection(cells_sf, protected)
covered <- aggregate(area ~ grid_id, FUN = sum,
                     data = data.frame(grid_id = pieces$grid_id, area = as.numeric(st_area(pieces))))

cells$protected_share <- 0                                     # no protected land -> 0
cells$protected_share[match(covered$grid_id, cells$grid_id)] <- covered$area / 25e6   # cell = 25 million m^2


#### Population ####
# WorldPop is a ~1 km grid in latitude/longitude, so I put the squares into latitude/longitude and add up the people in each
# exact = TRUE gives the fraction of each pixel inside the cell, so pixels on an edge are shared
pop <- rast(find_file(population_file))
px  <- extract(pop, vect(st_transform(cells_sf, 4326)), exact = TRUE)
per_cell <- tapply(px[[2]] * px$fraction, px$ID, sum, na.rm = TRUE)

cells$population <- 0                                          # all-sea cells -> 0
cells$population[as.integer(names(per_cell))] <- per_cell


#### Land cover: % of each cell in each habitat ####
# 21 layers, one per land-cover class, each the % of a 1 km square
# Layer order: 1 broadleaved woodland, 2 coniferous woodland, 3 arable, 4 improved grassland, 5 neutral grassland, 6 calcareous grassland, 7 acid grassland, 8 fen/marsh/swamp, 9 heather, 10 heather grassland, 11 bog, 12 inland rock, 13 saltwater, 14 freshwater, 15-18 supralittoral/littoral rock & sediment, 19 saltmarsh, 20 urban, 21 suburban.
# Each cell gets the average of the 1 km squares in it, weighted by how much of each square is inside the cell
lcm <- rast(find_file(landcover_file))
px  <- extract(lcm, vect(cells_sf), exact = TRUE)
pct <- as.matrix(px[, 2:22])

cell_mean <- function(classes) {
  v  <- rowSums(pct[, classes, drop = FALSE])
  ok <- !is.na(v)                                              # sea pixels are NA
  num <- tapply(v[ok] * px$fraction[ok], px$ID[ok], sum)
  den <- tapply(px$fraction[ok],         px$ID[ok], sum)
  out <- rep(0, nrow(cells))
  out[as.integer(names(num))] <- num / den
  out
}

cells$woodland_pct   <- cell_mean(1:2)      # red squirrel
cells$heath_pct      <- cell_mean(9:10)     # adder
cells$freshwater_pct <- cell_mean(14)       # otter, great crested newt, osprey
cells$coastal_pct    <- cell_mean(15:19)    # puffin, white-tailed eagle
cells$saltwater_pct  <- cell_mean(13)       # puffin, white-tailed eagle
cells$urban_pct      <- cell_mean(20:21)    # pipistrelle


#### Rivers: km of river in each cell ####
rivers <- st_read(find_file(rivers_file), layer = "watercourse_link", quiet = TRUE)
rivers <- st_geometry(st_transform(rivers, 27700))
rivers <- rivers[lengths(st_intersects(rivers, st_union(cells_sf))) > 0]   # just ones near our cells (faster)

bits <- st_intersection(cells_sf, rivers)
river_len <- aggregate(km ~ grid_id, FUN = sum,
                       data = data.frame(grid_id = bits$grid_id, km = as.numeric(st_length(bits)) / 1000))

cells$river_km <- 0
cells$river_km[match(river_len$grid_id, cells$grid_id)] <- river_len$km


#### Save ####
write.csv(cells, out_file, row.names = FALSE)

#### END ####
