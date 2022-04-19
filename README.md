
<!--
###############################################################################
######## DO NOT EDIT THIS FILE DIRECTLY. PLEASE READ THE COMMENT BELOW ########
###############################################################################

README.md is generated from README.Rmd. Please edit the README.Rmd file and
regenerate README.md by running the following in R:

rmarkdown::render('README.Rmd', output_format = 'github_document')

If you have GNU Make installed, you can also simply run "make" and Make will
automatically update README.md if it sees changes in README.Rmd.

###############################################################################
# DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE COMMENT ABOVE #
###############################################################################
-->

# cellmigRation

<img src="cell_migration_logo.png" width="50%" alt="cellmigRation">

An R package for tracking cells and analyzing their trajectories.

<!--
###############################################################################
#  DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE TOP COMMENT  #
###############################################################################
-->

# Installation

cellmigRation is under active development and a stable version is yet to
be released. However, you can install the current development version of
cellmigRation from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ocbe-uio/cellmigRation")
```

<!--
###############################################################################
#  DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE TOP COMMENT  #
###############################################################################
-->

# Example

Example usage:

``` r
library(cellmigRation)
#> Lastar nødvendig pakke: foreach

data(TrajectoryDataset)
df <- CellMig(TrajectoryDataset[1:100, ])
preproc <- rmPreProcessing(df, PixelSize = 1.24, TimeInterval = 100)
#> This dataset contains: 1 cell(s) in total
#> This dataset contains: 1 cell(s) with more than three steps in their tracks
#> The minimum number of steps:  100 
#> The maximum number of steps:  100 
#> Number of cells with a total number of steps less than  100 steps : 0 
#> All the tracks are adjusted to have only  100  steps
```

<!--
###############################################################################
#  DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE TOP COMMENT  #
###############################################################################
-->

# Shiny application

Some cellmigRation features are also available as a standalone web
application powered by an R package called Shiny.

You can run the Shiny app locally on your computer. To do so, download a
copy of this repository and either run `make runshiny` or `R -e
"shiny::runApp('Shinyapp', port=3029)"` from its root directory.

<!--
###############################################################################
#  DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE TOP COMMENT  #
###############################################################################
-->

# Statistics for developers

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/ocbe-uio/cellmigRation.svg?branch=master)](https://travis-ci.org/ocbe-uio/cellmigRation)
<!-- badges: end -->

<!--
###############################################################################
#  DO NOT EDIT README.md. YOU WILL LOSE YOUR CHANGES. PLEASE SEE TOP COMMENT  #
###############################################################################
-->
