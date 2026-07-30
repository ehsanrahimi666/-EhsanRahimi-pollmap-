# pollmap

<!-- badges: start -->
<!-- badges: end -->

**pollmap** is an R package for spatially modelling wild-bee pollination in
agricultural landscapes. It provides a raster-native, reproducible framework
that unifies:

- the **distance-decay family** of land-use models — Lonsdorf *et al.* (2009),
  the seasonal/guild formulation used in InVEST (following Kennedy *et al.* 2013
  and Olsson *et al.* 2015), and the quality-weighted, edge-sensitive variant of
  Rahimi *et al.* (2021);
- **species-distribution** pathways that couple pollinator occurrence models to
  pollination-service maps; and
- **pollinator-dependence demand** mapping and supply–demand matching.

The package targets ecological analysis (relative indices of supply, abundance
and demand) rather than economic valuation.

## Status

Feature-complete (v0.1.0), 73 unit tests, validated on the InVEST Willamette and
Swiss (CORINE / WorldClim / GBIF) datasets.

**Land-use supply** (one shared engine): `poll_kernel`, `poll_lonsdorf`,
`poll_kennedy` (InVEST), `poll_modified` (PollMap), `poll_estimap`, `poll_mce`,
`poll_cpf`.
**SDM pathway**: `poll_sdm`, `poll_sdm_project`, `poll_sdm_supply`.
**Demand**: `poll_crop_dependence`, `poll_demand`, `poll_balance`.
**Analysis**: `poll_scenario`, `poll_sensitivity` (+ `beta_landcover`),
`poll_validate`, `poll_compare`.
**Landscape / utilities**: `poll_simulate`, `poll_fragment`, `poll_metrics`,
`poll_classify`.
**Inputs**: `poll_guild`, `poll_landcover`, `poll_biophysical`, `poll_guilds`.

See the `pollmap-intro` vignette for worked examples and `NEWS.md` for the
roadmap.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ehsanrahimi666/-EhsanRahimi-pollmap-")
```

You will also need [`terra`](https://rspatial.github.io/terra/).

## Quick start

```r
library(pollmap)
library(terra)

# A 20 x 20 landscape: farmland (code 1) with a forest nesting strip (code 2)
lulc <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 2000, ymin = 0, ymax = 2000)
values(lulc) <- 1L
lulc[, 1:3] <- 2L

# Land-cover parameters: nesting and floral indices in [0, 1]
lc <- poll_landcover(
  data.frame(lucode  = c(1L, 2L),
             nesting = c(0.0, 1.0),
             floral  = c(0.8, 0.1))
)

# A pollinator guild with a 300 m mean foraging range
g <- poll_guild(alpha = 300, name = "wild_bee")

# Run the model
m <- poll_lonsdorf(lulc, lc, g)
m
summary(m)
plot(m)
```

## References

Lonsdorf, E., Kremen, C., Ricketts, T., Winfree, R., Williams, N. & Greenleaf,
S. (2009) Modelling pollination services across agricultural landscapes.
*Annals of Botany*, 103, 1589–1600.

Rahimi, E., Barghjelveh, S. & Dong, P. (2021) Using the Lonsdorf model for
estimating habitat loss and fragmentation effects on pollination service.
*Ecological Processes*, 10, 22.
