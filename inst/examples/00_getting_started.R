## ============================================================================
##  pollmap - getting started
##
##  A self-contained tour of the package. Everything below runs on small
##  synthetic landscapes and on the data bundled inside the package, so you do
##  NOT need to download the case-study GIS data to try it.
##
##  In RStudio: File > Open File... then click "Source" (Ctrl+Shift+S).
##  Or run it section by section.
## ============================================================================

## ---- 1. Install -----------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("ehsanrahimi666/-EhsanRahimi-pollmap-")

library(pollmap)
library(terra)

## ---- 2. What ships with the package ---------------------------------------
extdata <- system.file("extdata", package = "pollmap")
list.files(extdata)

# 26 wild bee species, with nesting substrate, seasonal activity, foraging
# distance and relative abundance (the last derived from GBIF record counts)
guilds <- poll_guilds(read.csv(file.path(extdata, "swiss_guild_table.csv")))
guilds
head(as.data.frame(guilds)[, c("SPECIES", "alpha", "relative_abundance")], 5)

# nesting and floral scores for 32 CORINE land-cover classes
bio <- poll_biophysical(read.csv(file.path(extdata, "clc_biophysical_table.csv")))
bio

# the occurrence records used in the paper (GBIF DOI 10.15468/dl.zhsrsb)
length(list.files(file.path(extdata, "occurrences")))

## ---- 3. Build a small CORINE-coded landscape ------------------------------
set.seed(1)
lulc <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 6000,
             ymin = 0, ymax = 6000, crs = "EPSG:2056")
values(lulc) <- 211L                 # arable land everywhere
lulc[8:22,  8:22]  <- 311L           # broad-leaved forest block
lulc[38:52, 6:20]  <- 231L           # pasture
lulc[10:20, 42:54] <- 112L           # discontinuous urban
lulc[30:36, 30:50] <- 243L           # agriculture with natural vegetation
names(lulc) <- "clc"
plot(lulc, main = "Synthetic CORINE-coded landscape")

## ---- 4. The original Lonsdorf model ---------------------------------------
lc <- poll_landcover(data.frame(
  lucode  = bio$lucode,
  nesting = pmax(bio$nesting_cavity_availability_index,
                 bio$nesting_ground_availability_index),
  floral  = rowMeans(bio[, c("floral_resources_spring_index",
                             "floral_resources_summer_index")])))
g <- poll_guild(alpha = 500, name = "generic_wild_bee")

m_lons <- poll_lonsdorf(lulc, lc, g)
m_lons
plot(m_lons)

## ---- 5. The seasonal, multi-guild InVEST model (all 26 species) -----------
m_kennedy <- poll_kennedy(lulc, bio, guilds)
m_kennedy
plot(m_kennedy)

## ---- 6. Compare five model formulations on identical inputs ---------------
scores <- data.frame(
  lucode = bio$lucode,
  availability = rowMeans(bio[, c("nesting_ground_availability_index",
                                  "floral_resources_spring_index",
                                  "floral_resources_summer_index")]))

supplies <- list(
  Lonsdorf = poll_lonsdorf(lulc, lc, g, layers = "supply")$map,
  Kennedy  = m_kennedy$map[["supply"]],
  Modified = poll_modified(lulc, lc, g, layers = "supply")$map,
  ESTIMAP  = poll_estimap(lulc, scores, g)$map,
  CPF      = poll_cpf(lulc, lc, g, layers = "supply")$map)

cmp <- poll_compare(supplies)
round(cmp$correlation, 3)            # who agrees with whom
plot(c(cmp$consensus, cmp$cv), main = c("consensus", "disagreement (CV)"))

## ---- 7. Demand, and the supply-demand balance -----------------------------
crop <- lulc * 0
crop[30:36, 30:50] <- 2310L          # fruit  (dependence 0.65)
crop[38:52, 6:20]  <- 1110L          # wheat  (dependence 0)
dep <- poll_crop_dependence(
  read.csv(file.path(extdata, "clms_crop_dependence.csv")),
  crop = "crop", dependence = "dependence")

demand <- poll_demand(crop, dep)
bal <- poll_balance(m_kennedy$map[["supply"]], demand)
plot(bal)                            # balance + a Boolean deficit layer

## ---- 8. Scenario: what if the forest is cleared? --------------------------
scen <- poll_scenario(poll_lonsdorf, lulc, lc, g, layers = "supply",
                      from = 311L, to = 211L)
plot(scen$delta, main = "Change in supply after clearing forest")

## ---- 9. Uncertainty in the resource scores --------------------------------
k <- 12
priors <- data.frame(lucode = lc$lucode,
                     nesting_shape1 = pmax(0.1, lc$nesting * k),
                     nesting_shape2 = pmax(0.1, (1 - lc$nesting) * k),
                     floral_shape1  = pmax(0.1, lc$floral * k),
                     floral_shape2  = pmax(0.1, (1 - lc$floral) * k))
draw <- beta_landcover(priors)
unc <- poll_sensitivity(function() poll_lonsdorf(lulc, draw(), g, layers = "supply"),
                        n = 20)
plot(unc$map[["cv"]], main = "Parameter uncertainty (CV of 20 draws)")

## ---- 10. Thinning real occurrence records ---------------------------------
occ <- read.csv(file.path(extdata, "occurrences", "Bombus_pascuorum.csv"))
nrow(occ)
occ_thin <- poll_thin(occ, dist = 10)          # one record per 10 km cell
nrow(occ_thin)
attr(occ_thin, "thin")

## ---- 11. A distribution model, fitted and honestly evaluated --------------
preds <- rast(nrows = 40, ncols = 40, nlyr = 2,
              xmin = 0, xmax = 400, ymin = 0, ymax = 400)
names(preds) <- c("temperature", "precipitation")
values(preds) <- cbind(runif(1600), runif(1600))

xy  <- xyFromCell(preds, which(values(preds[["temperature"]]) > 0.6))
occ_sim <- data.frame(Longitude = xy[, 1], Latitude = xy[, 2])

fit <- poll_sdm(occ_sim, preds, method = "glm", n_background = 300, seed = 1)
fit
suit <- poll_sdm_project(fit, preds)
plot(suit, main = "Modelled habitat suitability")

# k-fold cross-validation with AUC and TSS
ev <- poll_sdm_eval(occ_sim, preds, k = 3, block = FALSE,
                    n_background = 300, seed = 1)
ev

## ---- 12. Reporting: ordinal classes and area, structure by class ----------
cl <- poll_classify(m_kennedy$map[["supply"]], n = 5)
cl$area
plot(cl$map, main = "Pollination supply, five classes")

poll_metrics(m_kennedy$map[["supply"]], lulc)$class_summary

## ---- 13. Neutral landscapes for controlled experiments --------------------
land <- poll_simulate(nrow = 80, ncol = 80, p_habitat = 0.3, autocorr = 3,
                      seed = 1)
plot(land, main = "Neutral landscape (30% habitat, clustered)")

cat("\nDone. See help(package = 'pollmap') for the full function list,\n",
    "and the vignette: vignette('pollmap-intro').\n", sep = "")
