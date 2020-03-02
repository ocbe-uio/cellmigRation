
<!-- README.md is generated from README.Rmd. Please edit that file -->
CellMigRation
=============

![CellMigRation](CellMigRationLogo.png)

<!-- badges: start -->
[![Build Status](https://travis-ci.org/ocbe-uio/CellMigRation.svg?branch=master)](https://travis-ci.org/ocbe-uio/CellMigRation) <!-- badges: end -->

An R package for tracking cells and analyzing their trajectories.

Installation
------------

You can install the current version of CellMigRation from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ocbe-uio/CellMigRation")
```

Example
-------

Example usage:

``` r
library(CellMigRation)

data(Trajectory_dataset)
df <- CellMig(Trajectory_dataset[1:100, ])
preproc <- rmPreProcessing(df, PixelSize = 1.24, TimeInterval = 100)
#> This dataset contains:  1 Cells 
#> The minimum number of steps:  100 
#> The maximum number of steps:  100 
#> Number of cells with a total number of steps less than  100 steps : 0 
#> All the tracks are adjusted to have only  100  steps
```
