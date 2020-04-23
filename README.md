
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
-----------------

Some CellMigRation features are also available as a standalone web application powered by an R package called Shiny. Once published, the app will run from a remote server which doesn't require the end user to have an R installation. As this app is still under construction, though, it should be run from a local installation of R.

Here are two ways to run the Shiny app. Both assume the user has the `shiny` package installed.

### Running the GitHub-hosted application from a local R installation

1.  From an R terminal, execute the following command:

    ``` r
    shiny::runGitHub("CellMigRation", "ocbe-uio", "shiny", "Shinyapp")
    ```

### Running a 100% locally-hosted application

1.  Download the CellMigRation repository
2.  From the root folder of the (uncompressed) file, run the following in R:

    ``` r
    shiny::runApp("Shinyapp")
    ```

Step 1 can be skipped if you already have a clone of CellMigRation in your machine.

### Running and disconnecting from the Shiny application

Once the Shiny app is called, a browser window should open with the running application. After you are done using it, close the browser tab and press <kbd>Ctrl</kbd>+<kbd>C</kbd> in your R terminal to quit the application and regain control of the R prompt.
