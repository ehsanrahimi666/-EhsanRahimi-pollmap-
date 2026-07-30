# pollmap 0.1.0

First feature-complete version: all v1 modules implemented, tested (73 unit
tests), and validated on the InVEST Willamette and Swiss (CORINE / WorldClim /
GBIF) datasets.

## Land-use supply models (share one foraging engine)

* `poll_kernel()` — distance-weighted foraging convolution (engine).
* `poll_lonsdorf()` — original 2009 two-step model.
* `poll_kennedy()` — full InVEST: substrates x seasons x guilds + central-place
  foraging.
* `poll_modified()` — PollMap: quality-weighted, flight-range model.
* `poll_estimap()` — relative potential with abiotic activity multiplier.
* `poll_mce()` — weighted linear combination (+ AHP weights, consistency ratio).
* `poll_cpf()` — central-place-foraging (linear net-gain) variant.

## Species-distribution pathway

* `poll_sdm()`, `poll_sdm_project()`, `poll_sdm_supply()`, `predict.poll_model()`
  — ensemble occurrence models (glm/gam/rf; flexsdm-ready) coupled to a community
  supply surface.

## Demand

* `poll_crop_dependence()`, `poll_demand()`, `poll_balance()` — crop
  pollinator-dependence demand and supply-demand matching.

## Cross-cutting analysis

* `poll_scenario()`, `poll_sensitivity()` (+ `beta_landcover()`),
  `poll_validate()`, `poll_compare()`.

## Landscape and utilities

* `poll_simulate()`, `poll_fragment()`, `poll_metrics()`, `poll_classify()`.

## Inputs and data

* Input constructors: `poll_guild()`, `poll_landcover()`, `poll_biophysical()`,
  `poll_guilds()`.
* Bundled literature-based default tables: CORINE biophysical, Swiss guilds,
  CLMS crop dependence.
* S3 result class `poll_map` with `print`/`summary`/`plot`.

## Roadmap

* FFT fast-path for `poll_kernel()` (with a focal/FFT equivalence test).
* Full flexsdm ensemble integration and joint SDMs.
* Plant-pollinator network module; remotely sensed floral-resource ingestion.
