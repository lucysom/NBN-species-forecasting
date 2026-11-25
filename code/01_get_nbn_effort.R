# NBN species forecasting 2026: recording effort from the NBN Atlas
# Lucy Somekh

# Downloads the number of records per 1 km square per year (2007-2025) from the NBN Atlas, for all species & for a few species groups, in the two target regions.  
# Then adds the 1 km squares up into the 5 km grid cells.
# Output: extra_data/nbn_effort_by_cell.csv (grid_id, year, group, records)

# Why: a cell where lots of people record birds, moths, bats etc. is a cell where a given species is more likely to be recorded too
# Records of just the 9 species would give thinner picture of recording effort

# Downloaded on 18 Nov 2026
# N.b. The NBN Atlas keeps receiving records, so I asusme thatrunning this again later will give slightly different counts
# the file I used for the forecasts is saved in extra_data/ to make sure its reproducable


#### Settings ####
years <- 2007:2025

# Groups to count ("all" = every species), as NBN Atlas filters
groups <- c(
  all        = "",
  birds      = 'class:"Aves"',
  mammals    = 'class:"Mammalia"',
  bats       = 'order:"Chiroptera"',
  moths      = 'order:"Lepidoptera"',
  amphibians = 'class:"Amphibia"',
  reptiles   = 'class:"Reptilia"'
)

# Latitude/longitude boxes a little larger than each target region, so
# edge cells are fully covered
boxes <- list(
  south_east = c(lat_min = 50.45, lat_max = 52.30, lon_min = -2.05, lon_max = 1.55),
  scotland   = c(lat_min = 54.55, lat_max = 60.90, lon_min = -8.70, lon_max = -0.70)
)

raw_dir  <- "extra_data/nbn_effort_raw"
out_file <- "extra_data/nbn_effort_by_cell.csv"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)


#### Download record counts per 1 km square ####

# One request = one region, one group, one year
get_counts <- function(box, year, group_filter) {
  fq <- c(paste0("year:", year),
          '-occurrence_status:"absent"',
          paste0("latitude:[",  box["lat_min"], " TO ", box["lat_max"], "]"),
          paste0("longitude:[", box["lon_min"], " TO ", box["lon_max"], "]"))
  if (group_filter != "") fq <- c(fq, group_filter)
  url <- paste0("https://records-ws.nbnatlas.org/occurrences/facets/download",
                "?q=", URLencode("*:*", reserved = TRUE),
                paste0("&fq=", sapply(fq, URLencode, reserved = TRUE), collapse = ""),
                "&facets=grid_ref_1000&count=true&lookup=false")
  d <- read.csv(url, stringsAsFactors = FALSE)
  names(d) <- c("gridref", "records")
  d
}

for (region in names(boxes)) {
  for (g in names(groups)) {
    for (y in years) {
      f <- file.path(raw_dir, paste0(region, "_", g, "_", y, ".csv"))
      if (file.exists(f)) next                 # already downloaded
      write.csv(get_counts(boxes[[region]], y, groups[[g]]), f, row.names = FALSE)
      Sys.sleep(1)                             # pause between requests
    }
  }
}


#### Convert 1 km grid references to 5 km grid_id ####
# Convert to British National Grid metres, take the centre of the 1 km square, and find which 5 km cell it falls in
# A 1 km square that crosses a cell edge goes to the cell containing its centre

gridref_to_grid_id <- function(gr) {
  alphabet <- strsplit("ABCDEFGHJKLMNOPQRSTUVWXYZ", "")[[1]]   # OS grid letters skip "I"
  ok <- grepl("^[A-HJ-Z]{2}[0-9]{4}$", gr)                    # proper 1 km refs only
  g  <- gr[ok]
  l1 <- match(substr(g, 1, 1), alphabet) - 1
  l2 <- match(substr(g, 2, 2), alphabet) - 1
  easting  <- (((l1 - 2) %% 5) * 5 + (l2 %% 5)) * 100000 +
              as.numeric(substr(g, 3, 4)) * 1000 + 500
  northing <- ((19 - (l1 %/% 5) * 5) - (l2 %/% 5)) * 100000 +
              as.numeric(substr(g, 5, 6)) * 1000 + 500
  id <- rep(NA, length(gr))
  id[ok] <- ((northing - 353) %/% 5000) * 132 + (easting + 116) %/% 5000 - 131
  id
}


#### Add up into 5 km cells & save ####
cells <- read.csv("data/env_predictors/gridded/gridden_env_2025.csv")
target_cells <- cells$grid_id[cells$Region %in% c("South East", "Scotland")]

pieces <- list()
for (f in list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE)) {
  parts <- strsplit(sub("\\.csv$", "", basename(f)), "_")[[1]]   # e.g. south, east, birds, 2020
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$grid_id <- gridref_to_grid_id(d$gridref)
  d <- d[d$grid_id %in% target_cells, ]
  if (nrow(d) == 0) next
  s <- aggregate(records ~ grid_id, data = d, FUN = sum)
  s$year  <- as.integer(parts[length(parts)])
  s$group <- parts[length(parts) - 1]
  pieces[[f]] <- s
}

effort <- do.call(rbind, pieces)[, c("grid_id", "year", "group", "records")]
write.csv(effort, out_file, row.names = FALSE)

#### END ####
