#' @title Track And Analyze Cell Migrations
#' @description Track Fluorescent Cells and Analyze their movements.
#' Compute Migration Statistics and advanced metrics to understand
#' motility and movement characteristics of a population of cells.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no};
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @keywords internal
"_PACKAGE"




#' @title Sample Stack of Fluorescent Cells
#' @description Sample Stack of Fluorescent Cells to be used for
#' computing
#' cell tracks and stats
#' @usage data(TrackCellsDataset)
#' @format a \code{trackedCells} object including 10 stacks
#'
#' @examples
#' data(TrackCellsDataset)
#'
#' @keywords internal
"TrackCellsDataset"



#' Trajectories of 147 cells
#'
#' A dataset containing the coordinates and the ID of 147 cells from
#' wound scratch migration experiment
#'
#' @usage data(WSADataset)
#' @format A data frame with 11970 rows and 4 columns
#' @examples
#' data(WSADataset)
#'
#' @keywords internal
"WSADataset"


#' Trajectories of 350 cells
#'
#' A dataset containing the coordinates and the ID of 350 cells from
#' a dense random migration experiment
#'
#' @usage data(TrajectoryDataset)
#' @format A data frame with 50216 rows and 4 columns
#' @examples
#' data(TrajectoryDataset)
#'
"TrajectoryDataset"



#' @title Intermediates and Results from Cell Tracking Analyses
#' @description Intermediates and Results from Cell Tracking Analyses,
#' used to build the package vignette.
#' @usage data(VignBuilderDataset)
#' @format a list including 21 elements
#'
#' @examples
#' data(VignBuilderDataset)
#'
#' @keywords internal
"VignBuilderDataset"

