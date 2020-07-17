
<!--
###############################################################################
######## DO NOT EDIT THIS FILE DIRECTLY. PLEASE READ THE COMMENT BELOW ########
###############################################################################

README.md is generated from README.Rmd. Please edit the README.Rmd file and
regenerate README.md by running `rmarkdown::render("README.Rmd)` in R.

###############################################################################
# DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE COMMENT ABOVE #
###############################################################################
-->
cellmigRation
=============

<img src="cell_migration_logo.png" width="50%" alt="cellmigRation">

An R package for tracking cells and analyzing their trajectories.

Installation
============

cellmigRation is under active development and a stable version is yet to be released. However, you can install the current development version of cellmigRation from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ocbe-uio/cellmigRation")
```

Example
=======

Example usage:

``` r
library(cellmigRation)

data(TrajectoryDataset)
df <- CellMig(TrajectoryDataset[1:100, ])
preproc <- rmPreProcessing(df, PixelSize = 1.24, TimeInterval = 10)
#> This dataset contains: 1 cell(s) in total
#> This dataset contains: 1 cell(s) with more than three steps in their tracks
#> The minimum number of steps:  100 
#> The maximum number of steps:  100 
#> Number of cells with a total number of steps less than  100 steps : 0 
#> All the tracks are adjusted to have only  100  steps
```

Shiny application
=================

Some cellmigRation features are also available as a standalone web application powered by an R package called Shiny. **The Shiny app is still under construction**, so it is not feature-complete and is not guaranteed to perform as expected. Until release, we advise you to use the R package directly.

To access the cellmigRation Shiny app, visit <https://ocbe.shinyapps.io/cellmigRation/>.

Statistics for developers
=========================

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ocbe-uio/cellmigRation.svg?branch=master)](https://travis-ci.org/ocbe-uio/cellmigRation) <!-- badges: end -->
