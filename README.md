
<!-- README.md is generated from README.Rmd. Please edit that file -->
<img src="cell_migration_logo.png" width="50%" alt="CellMigRation">

An R package for tracking cells and analyzing their trajectories.


Installation
============

You can install the current version of CellMigRation from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ocbe-uio/CellMigRation")
```

Example
=======

Example usage:

``` r
library(CellMigRation)

data(TrajectoryDataset)
df <- CellMig(TrajectoryDataset[1:100, ])
preproc <- rmPreProcessing(df, PixelSize = 1.24, TimeInterval = 100)
#> This dataset contains:  1 Cells 
#> The minimum number of steps:  100 
#> The maximum number of steps:  100 
#> Number of cells with a total number of steps less than  100 steps : 0 
#> All the tracks are adjusted to have only  100  steps
```

Shiny application
=================

Some CellMigRation features are also available as a standalone web application powered by an R package called Shiny. **The Shiny app is still under construction**, so it is not feature-complete and is not guaranteed to perform as expected. Until release, we advise you to use the R package directly.

To access the CellMigRation Shiny app, visit <https://ocbe.shinyapps.io/CellMigRation/>.

Statistics for developers
=========================

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ocbe-uio/CellMigRation.svg?branch=master)](https://travis-ci.org/ocbe-uio/CellMigRation) <!-- badges: end -->
