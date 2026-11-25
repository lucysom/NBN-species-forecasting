# Introduction {#introduction .unnumbered}

This report outlines a workflow for forecasting where nine species will
be recorded in 2026, using recording effort and each cell's recent
record history, rather than habitat suitability.

Having tested various habitat/niche-based models, it appears that
whether each species is *recorded* in each 5 km cell depends as much (if
not more) on whether someone looks, and whether the record reaches the
NBN Atlas, as on habitat. So I forecast the recording process: recent
records in the cell, records in neighbouring cells (spread), and how
much recording happens there (effort).

I produce two versions of the forecast: **Model A**, the model as
described below, and **Model B**, the same model adjusted for delays in
records reaching the NBN Atlas (Section 3). Model B is my submitted
forecast.

# Load Data

I begin by loading the processed occurrence datasets for the nine
species.

    for (sp in names(species_regions)) {
      f <- file.path(species_dir, sp, paste0(tolower(sp), "_processed.csv"))
      d <- read.csv(f)
      occ_list[[sp]] <- d[, c("Occurrence.ID", "grid_id", "Year", "Region")]
    }

The raw downloads (4 August 2026) are used only for a few verified
January--June 2026 records and to identify which dataset each record
came from (Section 3). I also built a table of fixed cell
characteristics (`02_prepare_cell_covariates.R`): nature reserve and
SSSI cover (Natural England, NatureScot), population (WorldPop 2020),
land cover (UKCEH LCM 2025) and river length (OS Open Rivers). Sources
listed in `data_citations.csv`.

# Prepare effort data

Records of the nine contest species alone give a thin picture of
recording effort, so I also download record counts per 1 km square per
year (2007--2025) from the NBN Atlas, for all species and for birds,
mammals, bats, moths, amphibians and reptiles, and sum them into the
5 km grid (downloaded 18 November 2026).

Total recording grows almost every year and the latest year is
incomplete, so each count becomes the cell's share of the region's
records that year, averaged over the previous three years:

    effort_share <- function(group, region, ids, T) {
      shares <- sapply((T - 1):(T - 3), function(t) {
        n <- ext_count[paste(group, ids, t)]
        n[is.na(n)] <- 0
        (as.numeric(n) + 0.5) / reg_tot[[paste(region, group, t)]]
      })
      log(rowMeans(shares))
    }

# Model Structure

Every cell in the target region is included every year, so cells without
a record are real "not recorded" cells rather than randomly sampled
pseudo-absences.

The model uses seven predictors (eight for otters), all from years
before the forecast year:

- **History (3):** whether the species was recorded in the cell last
  year (yes/no); in how many of the two years before that (0--2); and
  whether it was recorded at any time 4--10 years ago (yes/no)

- **Neighbours (1):** how many of the 8 surrounding cells had a record
  in the last 3 years (0--8)

- **Effort (3):** records of the nine contest species in the cell over
  the last 3 years (log scale); the cell's share of the region's records
  for all species; and its share for the species' own group (e.g. moths
  for Jersey tiger)

- **River length (1, but otters only)**

A predictor was only added if it had an ecological reason and improved
the score in at least 6 of the 8 backtest years. Climate, land use,
elevation, soil, reserves/SSSIs, population and land-cover habitats all
failed this test: they change little between years, and a cell's record
history already reflects both suitable habitat and active recorders.
River length for otters passed in all 8 years.

**Reporting delay (Model B).** Records can reach the NBN Atlas months or
years after they are made, so if the 2026 answer key is extracted soon
after the end of 2026, some 2026 records will be missing. Model A
ignores this. I built Model B based on the assumption that the contest's
2026 data will be less complete than earlier years' data are now, and
that allowing for this will improve the forecast for the purposes of the
contest. For each species I estimated the share of recorded cells that
are usually reported quickly: I treated datasets that had already
submitted 2026 records by August 2026 as "fast", and calculated the
share of recorded cells in 2021--2023 that came from them. Model B
multiplies each Model A probability by this share: 0.42 for white-tailed
eagle (most of its records arrive late), 0.74--0.93 for the others, and
the median (0.78) for Jersey tiger (its flight season starts after the
supplied 2026 records end) and for puffin and osprey (too few records in
their target regions to estimate a factor). Model B still scored
slightly better than Model A (Table 1) so, even if the assumption is
wrong and the answer key is complete, little is lost.

    p_A <- fit_and_predict(sp, train, test)
    p_B <- p_A * delay_factor[[sp]]

# Model Training

A logistic regression is fitted per species on the previous 8 years
(2018--2025 for the 2026 forecast). Puffin and osprey have fewer than 20
presences in their regions, so they receive a flat base rate.
Predictions are kept between 0.001 and 0.99. Cells with verified
January--June 2026 records receive at least 0.9.

    model <- glm(reformulate(predictors, response = "y"),
                 family = binomial, data = train)
    p <- predict(model, newdata = newdata, type = "response")
    p <- pmin(pmax(p, prob_floor), prob_ceiling)

# Model Validation

I forecast each year 2018--2025 using only earlier data (a rolling
backtest) and scored each with mean binomial deviance per cell. The
total for model B (mean $\times$ number of cells) is given for
comparison with a summed deviance.

    binomial_deviance <- function(y, p) {
      -2 * mean(y * log(p) + (1 - y) * log(1 - p))
    }

  Species                               Region (cells)       Base rate     Model A     Model B   Model B, total
  ------------------------------------- ------------------ ----------- ----------- ----------- ----------------
  *Euplagia quadripunctaria*            South East (762)         0.443       0.310       0.298            227.0
  *Fratercula arctica*                  South East (762)         0.006       0.006       0.006              4.5
  *Haliaeetus albicilla*                South East (762)         0.137       0.121       0.128             97.8
  *Lutra lutra*                         South East (762)         0.485       0.401       0.399            304.1
  *Pandion haliaetus*                   Scotland (3,064)         0.011       0.011       0.011             33.0
  *Pipistrellus pipistrellus*           South East (762)         1.110       0.944       0.960            731.6
  *Sciurus vulgaris*                    South East (762)         0.194       0.058       0.058             43.9
  *Triturus cristatus*                  South East (762)         0.765       0.717       0.673            512.8
  *Vipera berus*                        South East (762)         0.458       0.299       0.295            224.9
  **Mean (species weighted equally)**                        **0.401**   **0.319**   **0.314**        **242.2**

  : Mean binomial deviance per cell, 2018--2025 backtest. Model A: no
  delay adjustment; Model B: adjusted for reporting delay.

"Base rate" gives every cell the same probability; Model B is 22%
better. Scored only on records from fast datasets (an "early extraction"
answer key), mean deviance is 0.296 for A and 0.275 for B. This partly
favours B by design, so is a scenario rather than a test. **I submit
Model B**; Model A is in `submission/alternative_no_delay_adjustment/`.

I re-ran the example approach (gradient boosting with 3:1
pseudo-absences) and scored it the same way for 2023--2025 (eight
species; it cannot be trained for osprey). Its mean deviance was 4.79,
against 0.31 here. Its probabilities are far too high (averaging 0.99
for otter, recorded in about 6% of cells).

# Reproducing {#reproducing .unnumbered}

All code is in R, in `code/`, and is run from the top folder of the
repository:

1.  `01_get_nbn_effort.R` downloads the effort counts (output saved in
    `extra_data/`; rerun would give slightly different counts as new
    records arrive)

2.  `02_prepare_cell_covariates.R` builds the cell covariates from the
    downloads listed in its header (output saved in `extra_data/`)

3.  `03_forecast_pipeline.R` runs the backtest and writes both
    forecasts.
