##
##
## ~~ All classes ~~
#
#

#' The trackedCells Class.
#'
#' An S4 class to represent a set of cells whose movements were tracked
#' over time.
#'
#' @slot images is a list of imported images
#' @slot proc_images is a list of processed images
#' @slot ops is a list keeping track of the operations executed on the object
#' @slot optimized is a list including results of the params
#' auto-optimization (optional)
#' @slot centroids is a list of detected centroids
#' @slot positions is a data.frame of cell positions across stacks
#' @slot tracks is a numeric matrix of cell tracks
#' @slot params is a list of parameters used for the analysis
#' @slot stats is a list of stats computed for the cell tracks
#' @slot metadata is a list including labels about the image, and
#' the experiment
#'
#' @author Damiano Fantini \email{damiano.fantini@@gmail.com}
#'
#' @name trackedCells-class
#' @rdname trackedCells-class
#' @exportClass trackedCells
#' @export
#' @return An S4-class object
trackedCells <- setClass(
  #Name of the class
  "trackedCells",

  #define slots
  slots = list(
    images="list",
    proc_images="list",
    ops = "list",
    optimized = "list",
    centroids = "list",
    positions = "data.frame",
    tracks = "matrix",
    params = "list",
    stats = "list",
    metadata = "list"
  ),

  # Make a function to check compatibility
  validity = function(object)
  {
    if(length(object@images) > 1){
      return(TRUE)
    } else {
      return("Malformed Data")
    }
  }
)



#' The CellMig Class.
#'
#' The CellMig class represents objects storing all information for both
#' random migration (RM) and wound scratch assay (WSA). It comprises 14 slots.
#'
#' @slot trajdata The raw trajectory data matrix organized into
#' four columns: cell ID, X coordinates, Y coordinates and Track number,
#' which is the track's path order.
#' @slot adjDS A data frame of the trajectory data passed
#' from the WSAprep function.
#' @slot cellpos A binary vector showing on which side of the wound cells
#' are located. "0" refers to a cell located above the wound whereas "1"
#' refers to a cell located below the wound.
#' @slot parE A numeric vector contains estimations for the
#' imageH, woundH, upperE and lowerE.
#' @slot preprocessedDS list object of data frames, each data frame
#' shows the trajectories of a single cell.
#' @slot DRtable A data frame of the results of running the
#' DiRatio() function.
#' @slot MSDtable A data frame of the results of running the MSD() function.
#' @slot PerAanSpeedtable A data frame of the results of running
#' the PerAndSpeed() function.
#' @slot DACtable A data frame of the results of running
#' the DiAutoCor() function.
#' @slot VACtable A data frame of the results of running
#' the VeAutoCor() function.
#' @slot ForMigtable A data frame of the results of running
#' the ForwardMigration() function.
#' @slot FMItable A data frame of the results of running the FMI() function.
#' @slot results A data frame of all the results.
#' @slot parCor A data frame for Parameters Correlation.
#' @slot meta A list including experiment name, meta data and other information.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#' @name CellMig-class
#' @rdname CellMig-class
#' @exportClass CellMig
#' @export
#' @return An S4-class object
CellMig <- setClass(
  # Name of your class
  "CellMig",

  slots = list(trajdata = "data.frame",
               adjDS= "data.frame",
               cellpos = "vector",
               parE= "vector",
               preprocessedDS= "list",
               DRtable= "data.frame",
               MSDtable= "data.frame",
               PerAanSpeedtable= "data.frame",
               DACtable= "data.frame",
               VACtable= "data.frame",
               ForMigtable= "data.frame",
               FMItable= "data.frame",
               results="data.frame",
               parCor="matrix",
               meta="list"),

  # Make a function to check compatibility
  validity = function(object)
  {
    msg <- NULL
    if ( nrow(object@trajdata) < 2 ){
      msg <- "input data must have more than one row"
    } else if ( ncol(object@trajdata) < 4 ){
      msg <- "input data must have four columns"
    }
    if (is.null(msg)) {
      return(TRUE)
    } else {
      return(msg)
    }
  })




