# NBN species forecasting 2026
# Lucy Somekh

# Overview: model based on the assumption that effectively what is being tested isn't is this cell suitable for this species?
# BUT: will this species be RECORDED in this cell in 2026?

# So the model is based on:
# (1) whether the species was recorded in that cell recently (recent presence)
# (2) whether it was recorded in neighbouring cells (to account for spread)
# (3) how much recording happens in that cell at all (to account for effort)
# (4) for otter only: length of river in the cell (the only habitat measure that passed the test below)

# For each species I fit a logistic regression on the above
# Tested by "forecasting" each past year 2018-2025 using only earlier data, then used to forecast 2026

# Two versions of the forecast:
#   A - the model as it is
#   B - A adjusted for reporting delay
#       Explained in submission report but, records reach the NBN Atlas months or years late, so if the 2026 answer key is extracted soon after the end of 2026 some records won't be there yet.  
#       To account for this, B multiplies each probability by the share of records that usually arrive quickly
# B is my actual submission but I have included A (in submission/alternative_no_delay_adjustment/) in case I've broken some rule I haven't spotted here by accounting for 2026 records...
# Sorry if that's the case! Hopefully can use A instead if needed?

# I've included some things I tested and left out (all made no consistent difference once a cell's history is known):
# climate, land use, elevation, soil, nature reserves/SSSIs, population, land-cover habitats
# Also tested dropping 2020 (COVID) from training but this made the model slightly worse, so I kept it
# These tests can basically be ignored but including for posterity

# Rule for adding a predictor: it needs an ecological reason AND must improve the score in at least 6 of the 8 backtest years

# Needs the two tables made by 01_get_nbn_effort.R & 02_prepare_cell_covariates.R (both saved in extra_data/)


#### Settings ####
species_regions <- c(
  Euplagia_quadripunctaria  = "South East",
  Fratercula_arctica        = "South East",
  Haliaeetus_albicilla      = "South East",
  Lutra_lutra               = "South East",
  Pandion_haliaetus         = "Scotland",
  Pipistrellus_pipistrellus = "South East",
  Sciurus_vulgaris          = "South East",
  Triturus_cristatus        = "South East",
  Vipera_berus              = "South East"
)

# The group of related species used to measure recording effort for each species (because people who record moths are the ones likely to record a Jersey tiger, bat groups record pipistrelles, etc.)
effort_group <- c(
  Euplagia_quadripunctaria  = "moths",
  Fratercula_arctica        = "birds",
  Haliaeetus_albicilla      = "birds",
  Lutra_lutra               = "mammals",
  Pandion_haliaetus         = "birds",
  Pipistrellus_pipistrellus = "bats",
  Sciurus_vulgaris          = "mammals",
  Triturus_cristatus        = "amphibians",
  Vipera_berus              = "reptiles"
)

# Species-specific extra predictors (only one that passed the 6-of-8 rule)
extra_predictors <- list(
  Lutra_lutra = "river_km"     # for otters, better in all 8 years
)

# Species mostly recorded in summer
# For these I can't tell fast-reporting datasets from slow ones because records only from January-June 2026, so they get the median delay factor of the others
summer_species <- c("Euplagia_quadripunctaria", "Fratercula_arctica", "Pandion_haliaetus")

forecast_year   <- 2026
backtest_years  <- 2018:2025  # years to test the model
train_window    <- 8          # how many past years to train on
min_presences   <- 20         # fewer presences than this in training -> use base rate
                              # (too little data to fit a regression reliably; puffin & osprey)
prob_floor      <- 0.001      # keep regression predictions away from 0 & 1: the deviance punishes a confident wrong answer enormously
prob_ceiling    <- 0.99
early_2026_prob <- 0.90       # probability for cells that already have a 2026 record

grid_width  <- 132            # the 5 km grid is 132 columns wide:
grid_offset <- 131            # grid_id = row * 132 + column - 131

species_dir     <- "data/NBN-atlas/species"
env_file        <- "data/env_predictors/gridded/gridden_env_2025.csv"  # used only for the list of cells
effort_file     <- "extra_data/nbn_effort_by_cell.csv"
covariates_file <- "extra_data/cell_covariates.csv"
out_dir         <- "submission"
alt_dir         <- "submission/alternative_no_delay_adjustment"


#### Load data ####

# The list of all 5 km grid cells, with their region & coordinates
cells_all <- read.csv(env_file)[, c("grid_id", "Region",
                                    "Longitude..WGS84.", "Latitude..WGS84.")]

# 110 grid cells have no region assigned.  No species records fall in them & they aren't part of any target region, so dropping them
cells_all <- cells_all[!is.na(cells_all$Region) & cells_all$Region != "", ]

# Processed occurrence records for each species
# Using the contests' processed files rather than the raw downloads, assuming the 2026 answer key will likely be built the same way
occ_list <- list()
for (sp in names(species_regions)) {
  f <- file.path(species_dir, sp, paste0(tolower(sp), "_processed.csv"))
  d <- read.csv(f)
  occ_list[[sp]] <- d[, c("Occurrence.ID", "grid_id", "Year", "Region")]
}

# Recorder effort (1): how many records of any of the 9 species each cell has in each year
# Trying to account for how often each cell is visited by recorders
all_records  <- do.call(rbind, occ_list)
effort_count <- table(paste(all_records$grid_id, all_records$Year))
effort_count <- setNames(as.numeric(effort_count), names(effort_count))

# Recorder effort (2): records of ALL species, & of each species group, per cell per year, downloaded from the NBN Atlas by get_nbn_effort.R
# Fuller picture of where recording happens than just the 9 species
ext <- read.csv(effort_file)
ext_count <- setNames(ext$records, paste(ext$group, ext$grid_id, ext$year))

# Region totals per group per year
# To turn a cell's count into its SHARE of the region's records: total recording grows almost every year (and the latest year is incomplete), so a share compares like with like
ext$Region <- cells_all$Region[match(ext$grid_id, cells_all$grid_id)]
reg_tot <- aggregate(records ~ Region + group + year, data = ext, FUN = sum)
reg_tot <- setNames(reg_tot$records, paste(reg_tot$Region, reg_tot$group, reg_tot$year))

# Fixed facts about each cell (made by prepare_cell_covariates.R)
# Only river_km is used (for otter)
covariates <- read.csv(covariates_file)

#### Helper functions ####
# Was the species recorded in each of these cells in year `yr`? (1 = yes)
# Turns every (cell, year) pair into a text key like "1938 2024"
recorded_in <- function(occ_keys, ids, yr) {
  as.integer(paste(ids, yr) %in% occ_keys)
}

# Number of records of the 9 species in each cell in year `yr` (0 if none)
effort_in <- function(ids, yr) {
  v <- effort_count[paste(ids, yr)]
  v[is.na(v)] <- 0
  as.numeric(v)
}

# A cell's share of its region's records (for one species group), averaged over the last 3 years, on a log scale
# The +0.5 stops cells w/ no records getting log(0)
effort_share <- function(group, region, ids, T) {
  shares <- sapply((T - 1):(T - 3), function(t) {
    n <- ext_count[paste(group, ids, t)]
    n[is.na(n)] <- 0
    (as.numeric(n) + 0.5) / reg_tot[[paste(region, group, t)]]
  })
  log(rowMeans(shares))
}

# Build table for one species, one region & one target year T
# One row per grid cell
# Every predictor uses only years BEFORE T
make_features <- function(sp, occ_keys, cells, T) {
  id     <- cells$grid_id
  region <- species_regions[[sp]]

  # The outcome: recorded in year T?  (For T = 2026 this is unknown so just a column of zeroes, not used)
  y <- recorded_in(occ_keys, id, T)

  # History in this cell
  last  <- recorded_in(occ_keys, id, T - 1)                    # last year (0/1)
  prev2 <- recorded_in(occ_keys, id, T - 2) +
           recorded_in(occ_keys, id, T - 3)                    # the 2 years before (0-2)
  older_mat <- sapply((T - 10):(T - 4), function(t) recorded_in(occ_keys, id, t))
  older <- as.integer(rowSums(older_mat) > 0)                  # any record 4-10 years ago (0/1)

  # Neighbours: how many of the 8 surrounding cells had a record in the last 3 years (0-8).  Captures spread into new cells
  row <- (id + grid_offset) %/% grid_width
  col <- (id + grid_offset) %%  grid_width
  nb  <- 0
  for (dr in -1:1) {
    for (dc in -1:1) {
      if (dr == 0 && dc == 0) next                             # skip the cell itself
      h <- (row + dr) * grid_width + (col + dc) - grid_offset  # neighbour's grid_id
      nb_recent <- recorded_in(occ_keys, h, T - 1) +
                   recorded_in(occ_keys, h, T - 2) +
                   recorded_in(occ_keys, h, T - 3)
      nb <- nb + as.integer(nb_recent > 0)
    }
  }

  # Recorder effort over the last 3 years
  effort    <- log1p(effort_in(id, T - 1) + effort_in(id, T - 2) + effort_in(id, T - 3))  # 9 contest species, log scale
  e_all_rel <- effort_share("all", region, id, T)                # all species, share of region
  e_grp_rel <- effort_share(effort_group[[sp]], region, id, T)   # this species' group, share of region

  # Length of river in the cell (used for otters only)
  river_km <- covariates$river_km[match(id, covariates$grid_id)]

  data.frame(grid_id = id, year = T, y = y, last = last, prev2 = prev2,
             older = older, nb = nb, effort = effort,
             e_all_rel = e_all_rel, e_grp_rel = e_grp_rel, river_km = river_km)
}

# Fit the model on a training table & predict for a new table
fit_and_predict <- function(sp, train, newdata) {
  # Base rate = share of cell-years occupied in the last 5 training years w/ a small adjustment (+0.5, +1) so it's never exactly 0
  recent  <- train[train$year >= max(train$year) - 4, ]
  base    <- (sum(recent$y) + 0.5) / (nrow(recent) + 1)

  if (sum(train$y) < min_presences) {
    # Too few presences to fit a regression (e.g. puffin in the South East): predict the same small base rate everywhere
    p <- rep(base, nrow(newdata))
  } else {
    # Logistic regression: the probability of a record depends on history, neighbours, effort (& river length for otters)
    predictors <- c("last", "prev2", "older", "nb", "effort", "e_all_rel", "e_grp_rel",
                    extra_predictors[[sp]])
    model <- suppressWarnings(
      glm(reformulate(predictors, response = "y"), family = binomial, data = train)
    )
    p <- predict(model, newdata = newdata, type = "response")

    # Keep probabilities away from exactly 0 & 1
    # Not applied to the base rate above: for puffin the honest probability is below 0.001)
    p <- pmin(pmax(p, prob_floor), prob_ceiling)
  }

  p
}

# Mean binomial deviance per grid cell (lower is better)
# deviance = -2 * mean( y*log(p) + (1-y)*log(1-p) )
# Following the example notebook: sum over cells instead of averaging
binomial_deviance <- function(y, p) {
  -2 * mean(y * log(p) + (1 - y) * log(1 - p))
}

# Turn all record IDs into text so they can be matched between files
id_text <- function(x) {
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE) else as.character(x)
}

# Assign each record to the nearest grid-cell centre
# Records more than 4 km from any cell centre are dropped
nearest_cell <- function(lon, lat) {
  km_x <- function(lon, lat) lon * 111.32 * cos(lat * pi / 180)
  cx <- km_x(cells_all$Longitude..WGS84., cells_all$Latitude..WGS84.)
  cy <- cells_all$Latitude..WGS84. * 110.57
  sapply(seq_along(lon), function(i) {
    d2 <- (cx - km_x(lon[i], lat[i]))^2 + (cy - lat[i] * 110.57)^2
    if (min(d2) > 4^2) NA else cells_all$grid_id[which.min(d2)]
  })
}


#### Reporting delay & early-2026 records (from the raw downloads) ####
# The raw NBN downloads were taken on 4 Aug 2026
# For each species I use them to identify:
# 1. Fast datasets: datasets that had already sent in 2026 records by August 2026.  These report quickly
#
# 2. Delay factor: in 2021-2023, what share of the species' recorded cells had at least one record from a fast dataset?
#    That's roughly the share of cells that I EXPECT would already be in the answer key if it were extracted soon after the end of 2026.  
#    Used for model B.
#    For summer species the Jan-Jun 2026 records can't show which datasets are fast, so they get the median of the others.
#
# 3. Early-2026 records: a few verified records already made in 2026 (Jan-Jun).
#    Contests' original processing left them out, but a cell that already has a 2026 record will most likely count as occupied, so I raise its probability
#    Follow same quality filters as original processing (code/env-data-processing.R)

delay_factor <- c()
fast_keys    <- list()   # cell-years recorded by fast datasets (for "early extraction" check)
early_2026   <- list()

for (sp in names(species_regions)) {
  raw_file <- list.files(file.path(species_dir, sp), pattern = "^records-.*\\.csv$",
                         full.names = TRUE)[1]
  raw <- read.csv(raw_file)

  # 1. Fast datasets
  fast_sets <- unique(raw$Dataset.name[which(raw$Year == forecast_year)])

  # Link each processed record to the dataset it came from
  proc <- occ_list[[sp]]
  proc$dataset <- raw$Dataset.name[match(id_text(proc$Occurrence.ID), id_text(raw$Occurrence.ID))]
  proc$fast <- proc$dataset %in% fast_sets
  fast_keys[[sp]] <- unique(paste(proc$grid_id[proc$fast], proc$Year[proc$fast]))

  # 2. Delay factor (needs at least 30 recorded cell-years to be worth estimating)
  r <- proc[proc$Region == species_regions[[sp]] & proc$Year %in% 2021:2023, ]
  if (nrow(r) > 0 && !(sp %in% summer_species)) {
    cell_years <- aggregate(fast ~ grid_id + Year, data = r, FUN = any)
    if (nrow(cell_years) >= 30) delay_factor[sp] <- mean(cell_years$fast)
  }

  # 3. Early-2026 records
  keep <- raw$Year == forecast_year &
          !is.na(raw$Coordinate.Uncertainty..m.) &
          raw$Coordinate.Uncertainty..m. <= 5000 &
          raw$Identification.verification.status %in%
            c("Accepted", "Accepted - considered correct", "Accepted - correct") &
          !is.na(raw$Latitude..WGS84.)
  raw <- raw[which(keep), ]
  early_2026[[sp]] <- if (nrow(raw) == 0) integer(0) else
    unique(na.omit(nearest_cell(raw$Longitude..WGS84., raw$Latitude..WGS84.)))
  rm(raw)
}

# Species without their own factor get the median of the others
missing <- setdiff(names(species_regions), names(delay_factor))
delay_factor[missing] <- median(delay_factor)
delay_factor <- delay_factor[names(species_regions)]
cat("Delay factors (share of cells reported quickly):\n")
print(round(delay_factor, 2))


#### Backtest: forecast each year 2018-2025 using only earlier data & score it ####
# Scored 2 ways:
# (1) against the data as it is now (complete for these past years)
# (2) against an "early extraction" version that only counts records from fast datasets
# (N.b. Summer species are scored on the full data both times, as fast datasets can't be identified)

backtest      <- data.frame()
feature_store <- list()   # keep each species' feature tables for the final forecast

for (sp in names(species_regions)) {
  region   <- species_regions[[sp]]
  cells    <- cells_all[cells_all$Region == region, ]
  occ_keys <- unique(paste(occ_list[[sp]]$grid_id, occ_list[[sp]]$Year))

  # Build the feature table for every year
  years_needed <- (min(backtest_years) - train_window):forecast_year
  feats <- lapply(years_needed, function(t) make_features(sp, occ_keys, cells, t))
  names(feats) <- years_needed
  feature_store[[sp]] <- feats

  for (T in backtest_years) {
    train <- do.call(rbind, feats[as.character((T - train_window):(T - 1))])
    test  <- feats[[as.character(T)]]

    # Base rate: same probability in every cell (for comparison)
    recent <- train[train$year >= T - 5, ]
    p_base <- rep((sum(recent$y) + 0.5) / (nrow(recent) + 1), nrow(test))

    # Model A, and model B = A x delay factor
    p_A <- fit_and_predict(sp, train, test)
    p_B <- p_A * delay_factor[[sp]]

    # "Early extraction" version of the truth
    y_early <- if (sp %in% summer_species) test$y else recorded_in(fast_keys[[sp]], test$grid_id, T)

    backtest <- rbind(backtest, data.frame(
      species = sp, year = T, cells = nrow(test), presences = sum(test$y),
      deviance_base_rate = binomial_deviance(test$y, p_base),
      deviance_A         = binomial_deviance(test$y, p_A),
      deviance_B         = binomial_deviance(test$y, p_B),
      early_deviance_A   = binomial_deviance(y_early, p_A),
      early_deviance_B   = binomial_deviance(y_early, p_B)
    ))
  }
  cat("Backtest done:", sp, "\n")
}

# Average over the 8 backtest years for each species, then over species
# (each species counts equally)
score_cols  <- c("presences", "deviance_base_rate", "deviance_A", "deviance_B",
                 "early_deviance_A", "early_deviance_B")
summary_tab <- aggregate(backtest[, score_cols], by = list(species = backtest$species), FUN = mean)
mean_row    <- data.frame(species = "MEAN (all species equal)", t(colMeans(summary_tab[, score_cols])))
mean_row$presences <- NA
summary_tab <- rbind(summary_tab, mean_row)
print(summary_tab, digits = 3)

dir.create(alt_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(backtest,    file.path(out_dir, "backtest_deviance_by_year.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "backtest_deviance_summary.csv"), row.names = FALSE)


#### Final forecasts for 2026 ####
# Same model but trained on 2018-2025
for (sp in names(species_regions)) {
  region <- species_regions[[sp]]
  feats  <- feature_store[[sp]]
  train  <- do.call(rbind, feats[as.character((forecast_year - train_window):(forecast_year - 1))])
  new    <- feats[[as.character(forecast_year)]]

  p_A <- fit_and_predict(sp, train, new)
  p_B <- p_A * delay_factor[[sp]]

  # Cells that already have a 2026 record get at least early_2026_prob
  already <- new$grid_id %in% early_2026[[sp]]
  p_A[already] <- pmax(p_A[already], early_2026_prob)
  p_B[already] <- pmax(p_B[already], early_2026_prob)

  cells <- cells_all[match(new$grid_id, cells_all$grid_id), ]
  file_name <- paste0("forecast_", tolower(sp), "_", forecast_year, "_",
                      gsub(" ", "_", region), ".csv")
  for (version in c("B", "A")) {
    out <- data.frame(
      grid_id                     = new$grid_id,
      Region                      = region,
      Longitude..WGS84.           = cells$Longitude..WGS84.,
      Latitude..WGS84.            = cells$Latitude..WGS84.,
      occurrence_probability_2026 = if (version == "B") p_B else p_A,
      forecast_year               = forecast_year
    )
    folder <- if (version == "B") out_dir else alt_dir
    write.csv(out, file.path(folder, file_name), row.names = FALSE, quote = FALSE)
  }

  cat(sprintf("%-27s %4d cells | expected occupied cells A %6.1f, B %6.1f | cells with early-2026 record %d\n",
              sp, nrow(new), sum(p_A), sum(p_B), sum(already)))
}

#### END ####
