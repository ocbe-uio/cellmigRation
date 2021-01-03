##
##
## ~~ All Methods ~~
#
# 

# ~~~
# 1) trackedCells-methods
# ~~~

#' Constructor method of the trackedCells Class.
#'
#' @name trackedCells
#' @rdname trackedCells-class
#' @param .Object the trackedCells object being built
#' @param x imported TIFF image data
#' @aliases initialize,trackedCells-method
#' @import methods
#' 
setMethod("initialize", "trackedCells",
          function(.Object, x) {
            .Object <- callNextMethod(.Object)
            
            # check data and args
            chk <- (sum(names(x) %in% c("images", "dim", "attributes")) == 3)
            
            if(!chk)
              stop("Malformed Data")
            
            # Initialize
            P0 <- data.frame(row = 1, col = 1, tau = 1)
            P0 <- P0[-1, ]
            
            T0 <- matrix(NA, ncol = 4, nrow = 0)
            
            O0 <- list(images = 1,
                       optimized_params = 0,
                       custom_params = 0,
                       track = 0,
                       stats = 0
            )
            
            #Assign
            .Object@images <- x
            .Object@proc_images <- list()
            .Object@ops <- O0
            .Object@optimized <- list()
            .Object@centroids <- list()
            .Object@positions <- P0
            .Object@tracks <- T0
            .Object@params <- list()
            .Object@stats <- list()
            .Object@metadata <- list(tiff_file = NA,
                                     experiment = NA,
                                     condition = NA,
                                     replicate = NA)
            
            # return
            .Object
            
          })


setMethod("show", signature(object = "trackedCells"),
          function(object) {
            LNS <- list(
              " ~~~ An S4 trackedCells object ~~~",
              "",
              paste0("      + Num of images: ", object@images$dim$NumberImages),
              paste0("      + Optimized Params: ", ifelse(object@ops$optimized_params == 1, "Yes", "No")),
              paste0("      + Run w/ Custom Params: ", ifelse(object@ops$custom_params == 1, "Yes", "No")),
              paste0("      + Cells Tracked: ", ifelse(object@ops$track == 1, "Yes", "No")),
              paste0("      + Stats Computed: ", ifelse(object@ops$stats == 1, "Yes", "No")),
              "")
            
            for(lni in LNS){
              cat(lni, sep = "\n")
            }
          })




## ~~~~~~~~~~~~~~~~~~

#' Method getCellImages
#' 
#' Retrieve Images from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getCellImages
#' @exportMethod getCellImages
setGeneric("getCellImages", function(x) {
  standardGeneric("getCellImages")
})

#' @rdname getCellImages
#' @aliases getCellImages,trackedCells-method
setMethod("getCellImages", "trackedCells",
          function(x){
            x@images
          })


#' Method getProcessedImages
#' 
#' Retrieve Processed Images from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getProcessedImages
#' @exportMethod getProcessedImages
setGeneric("getProcessedImages", function(x) {
  standardGeneric("getProcessedImages")
})

#' @rdname getProcessedImages
#' @aliases getProcessedImages,trackedCells-method
setMethod("getProcessedImages", "trackedCells",
          function(x){
            x@proc_images
          })


#' Method getImageCentroids
#' 
#' Retrieve Image Centroids from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getImageCentroids
#' @exportMethod getImageCentroids
setGeneric("getImageCentroids", function(x) {
  standardGeneric("getImageCentroids")
})

#' @rdname getImageCentroids
#' @aliases getImageCentroids,trackedCells-method
setMethod("getImageCentroids", "trackedCells",
          function(x){
            x@centroids
          })




#' Method getProcessingStatus
#' 
#' Retrieve Processing Status from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getProcessingStatus
#' @exportMethod getProcessingStatus
setGeneric("getProcessingStatus", function(x) {
  standardGeneric("getProcessingStatus")
})

#' @rdname getProcessingStatus
#' @aliases getProcessingStatus,trackedCells-method
setMethod("getProcessingStatus", "trackedCells",
          function(x){
            x@ops
          })

#' Method getCellTracks
#' 
#' Retrieve Cell Tracks from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getCellTracks
#' @exportMethod getCellTracks
setGeneric("getCellTracks", function(x) {
  standardGeneric("getCellTracks")
})

#' @rdname getCellTracks
#' @aliases getCellTracks,trackedCells-method
setMethod("getCellTracks", "trackedCells",
          function(x){
            x@tracks
          })


#' Method getCellTrackMeta
#' 
#' Retrieve Metadata from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getCellTrackMeta
#' @exportMethod getCellTrackMeta
setGeneric("getCellTrackMeta", function(x) {
  standardGeneric("getCellTrackMeta")
})

#' @rdname getCellTrackMeta
#' @aliases getCellTrackMeta,trackedCells-method
setMethod("getCellTrackMeta", "trackedCells",
          function(x){
            x@metadata
          })


#' Method getCellTrackStats
#' 
#' Retrieve Stats from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getCellTrackStats
#' @exportMethod getCellTrackStats
setGeneric("getCellTrackStats", function(x) {
  standardGeneric("getCellTrackStats")
})

#' @rdname getCellTrackStats
#' @aliases getCellTrackStats,trackedCells-method
setMethod("getCellTrackStats", "trackedCells",
          function(x){
            x@stats
          })




#' Method getOptimizedParams
#' 
#' Retrieve Optimized Params from a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' 
#' @rdname getOptimizedParams
#' @exportMethod getOptimizedParams
setGeneric("getOptimizedParams", function(x) {
  standardGeneric("getOptimizedParams")
})

#' @rdname getOptimizedParams
#' @aliases getOptimizedParams,trackedCells-method
setMethod("getOptimizedParams", "trackedCells",
          function(x){
            x@optimized
          })







## ~~~~~~~~~~~~~~~~~~


#' Method setTrackedCellsMeta
#' 
#' Set Metadata of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param meta a list including all metadata
#' 
#' @rdname setTrackedCellsMeta
#' @exportMethod setTrackedCellsMeta
setGeneric("setTrackedCellsMeta", function(x, meta) {
  standardGeneric("setTrackedCellsMeta")
})


#' @rdname setTrackedCellsMeta
#' @aliases setTrackedCellsMeta,trackedCells,list-method
setMethod("setTrackedCellsMeta", 
          signature(x="trackedCells", 
                    meta = "list"),
          function(x, meta){
            x@metadata <- meta
            return(x)
          })


#' Method setProcessingStatus
#' 
#' Set Operation Status of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param slot string pointing to the slot to be updated
#' @param value numeric value to be written
#' 
#' @rdname setProcessingStatus
#' @exportMethod setProcessingStatus
setGeneric("setProcessingStatus", function(x, slot, value) {
  standardGeneric("setProcessingStatus")
})


#' @rdname setProcessingStatus
#' @aliases setProcessingStatus,trackedCells,character,numeric-method
setMethod("setProcessingStatus", 
          signature(x="trackedCells", 
                    slot = "character", 
                    value = "numeric"),
          function(x, slot, value){
            x@ops[[slot]] <- value
            return(x)
          })


#' Method setTrackingStats
#' 
#' Set Tracking Statistics of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param population population-level statistics
#' @param cells cell-level statistics
#' 
#' @rdname setTrackingStats
#' @exportMethod setTrackingStats
setGeneric("setTrackingStats", function(x, population, cells) {
  standardGeneric("setTrackingStats")
})


#' @rdname setTrackingStats
#' @aliases setTrackingStats,trackedCells,ANY,ANY-method
setMethod("setTrackingStats", 
          signature(x="trackedCells", 
                    population = "ANY", 
                    cells = "ANY"),
          function(x, population, cells){
            x@stats <- list(population = population, cells = cells)
            return(x)
          })





#' Method setOptimizedParams
#' 
#' Set Optimized Params of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param auto_params automatically selected parameters
#' @param results optimization analysis results
#' 
#' @rdname setOptimizedParams
#' @exportMethod setOptimizedParams
setGeneric("setOptimizedParams", function(x, auto_params, results) {
  standardGeneric("setOptimizedParams")
})


#' @rdname setOptimizedParams
#' @aliases setOptimizedParams,trackedCells,ANY,ANY-method
setMethod("setOptimizedParams", 
          signature(x="trackedCells", 
                    auto_params = "ANY", 
                    results = "ANY"),
          function(x, auto_params, results){
            x@optimized <- list(auto_params = auto_params, results = results)
            return(x)
          })




#' Method setProcessedImages
#' 
#' Set Processed Images of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param procImages a list including all metadata
#' 
#' @rdname setProcessedImages
#' @exportMethod setProcessedImages
setGeneric("setProcessedImages", function(x, procImages) {
  standardGeneric("setProcessedImages")
})


#' @rdname setProcessedImages
#' @aliases setProcessedImages,trackedCells,list-method
setMethod("setProcessedImages", 
          signature(x="trackedCells", 
                    procImages = "list"),
          function(x, procImages){
            x@proc_images <- procImages
            return(x)
          })


#' Method setTrackedCentroids
#' 
#' Set Centroids of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param centroids a list including all metadata
#' 
#' @rdname setTrackedCentroids
#' @exportMethod setTrackedCentroids
setGeneric("setTrackedCentroids", function(x, centroids) {
  standardGeneric("setTrackedCentroids")
})


#' @rdname setTrackedCentroids
#' @aliases setTrackedCentroids,trackedCells,list-method
setMethod("setTrackedCentroids", 
          signature(x="trackedCells", 
                    centroids = "list"),
          function(x, centroids){
            x@centroids <- centroids
            return(x)
          })


#' Method setTrackedPositions
#' 
#' Set positions of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param positions a data.frame including all positions
#' 
#' @rdname setTrackedPositions
#' @exportMethod setTrackedPositions
setGeneric("setTrackedPositions", function(x, positions) {
  standardGeneric("setTrackedPositions")
})


#' @rdname setTrackedPositions
#' @aliases setTrackedPositions,trackedCells,data.frame-method
setMethod("setTrackedPositions", 
          signature(x="trackedCells", 
                    positions = "data.frame"),
          function(x, positions){
            x@positions <- positions
            return(x)
          })




#' Method setCellTracks
#' 
#' Set Tracks of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param tracks a matrix including all cell tracks
#' 
#' @rdname setCellTracks
#' @exportMethod setCellTracks
setGeneric("setCellTracks", function(x, tracks) {
  standardGeneric("setCellTracks")
})


#' @rdname setCellTracks
#' @aliases setCellTracks,trackedCells,matrix-method
setMethod("setCellTracks", 
          signature(x="trackedCells", 
                    tracks = "matrix"),
          function(x, tracks){
            x@tracks <- tracks
            return(x)
          })




#' Method setAnalyticParams
#' 
#' Set Analytic Params of a \code{trackedCells} object.
#' 
#' @param x a \code{trackedCells}-class object
#' @param params a list including all params
#' 
#' @rdname setAnalyticParams
#' @exportMethod setAnalyticParams
setGeneric("setAnalyticParams", function(x, params) {
  standardGeneric("setAnalyticParams")
})


#' @rdname setAnalyticParams
#' @aliases setAnalyticParams,trackedCells,list-method
setMethod("setAnalyticParams", 
          signature(x="trackedCells", 
                    params = "list"),
          function(x, params){
            x@params <- params
            return(x)
          })




















## ~~~~~~~~~~~~~~~~~~




# ~~~
# 2) CellMig-methods
# ~~~

#' Constructor method of the CellMig Class.
#'
#' @name CellMig
#' @rdname CellMig-class
#' @param .Object the CellMig object being built
#' @param trajdata data frame including trajectory data
#' @aliases initialize,CellMig-method
#' @import methods
#' 
setMethod("initialize",
          signature = "CellMig",
          definition = function(.Object, trajdata){
            .Object <- callNextMethod(.Object)
            .Object@trajdata <- trajdata
            
            return(.Object)
          }
)


#' Wrapper function CellMig.
#'
#' @name CellMig
#' @rdname CellMig-class
#' @param ... arguments to pass to the CellMig constructor
#' @importFrom methods new
#' @export
CellMig <- function(...) new("CellMig", ...) 


setMethod("show", signature(object = "CellMig"),
          function(object) {
            LNS <- list(
              " ~~~ An S4 CellMig object ~~~",
              "")
            
            for(lni in LNS){
              cat(lni, sep = "\n")
            }
          })


#' Method getCellMigSlot
#' 
#' Get Data from a slot in a \code{CellMig} object.
#' 
#' @param x a \code{CellMig}-class object
#' @param slot string pointing to the slot to be retrieved
#' 
#' @rdname getCellMigSlot
#' @exportMethod getCellMigSlot
setGeneric("getCellMigSlot", function(x, slot) {
  standardGeneric("getCellMigSlot")
})


#' @rdname getCellMigSlot
#' @aliases getCellMigSlot,CellMig,character-method
setMethod("getCellMigSlot", 
          signature("CellMig", "character"),
          function(x, slot){
            objList <- list(trajdata = x@trajdata, 
                            adjDS = x@adjDS,
                            cellpos = x@cellpos, 
                            parE = x@parE, 
                            preprocessedDS = x@preprocessedDS, 
                            DRtable = x@DRtable, 
                            MSDtable = x@MSDtable, 
                            PerAanSpeedtable = x@PerAanSpeedtable, 
                            DACtable = x@DACtable, 
                            VACtable = x@VACtable, 
                            ForMigtable = x@ForMigtable, 
                            FMItable = x@FMItable, 
                            results = x@results, 
                            parCor = x@parCor, 
                            meta = x@meta)
            
            y <- objList[[slot[1] ]]
            return(y)
          })




#' Method setCellMigSlot
#' 
#' Set Data of a slot in a \code{CellMig} object.
#' 
#' @param x a \code{CellMig}-class object
#' @param slot string pointing to the slot to be updated
#' @param value ANY value to be written
#' 
#' @rdname setCellMigSlot
#' @exportMethod setCellMigSlot
setGeneric("setCellMigSlot", function(x, slot, value) {
  standardGeneric("setCellMigSlot")
})





#' @rdname setCellMigSlot
#' @aliases setCellMigSlot,CellMig,character,ANY-method
setMethod("setCellMigSlot", 
          signature(x="CellMig", 
                    slot = "character", value = "ANY"),
          function(x, slot, value){
            objList <- list(trajdata = x@trajdata, 
                            adjDS = x@adjDS,
                            cellpos = x@cellpos, 
                            parE = x@parE, 
                            preprocessedDS = x@preprocessedDS, 
                            DRtable = x@DRtable, 
                            MSDtable = x@MSDtable, 
                            PerAanSpeedtable = x@PerAanSpeedtable, 
                            DACtable = x@DACtable, 
                            VACtable = x@VACtable, 
                            ForMigtable = x@ForMigtable, 
                            FMItable = x@FMItable, 
                            results = x@results, 
                            parCor = x@parCor, 
                            meta = x@meta)
            
            objList[[ as.character(slot)[1] ]] <- value
            
            x@trajdata <- objList$trajdata
            x@adjDS <- objList$adjDS
            x@cellpos <- objList$cellpos
            x@parE <- objList$parE
            x@preprocessedDS <- objList$preprocessedDS
            x@DRtable <- objList$DRtable
            x@MSDtable <- objList$MSDtable
            x@PerAanSpeedtable <- objList$PerAanSpeedtable
            x@DACtable <- objList$DACtable
            x@VACtable <- objList$VACtable
            x@ForMigtable <- objList$ForMigtable
            x@FMItable <- objList$FMItable
            x@results <- objList$results
            x@parCor <- objList$parCor
            x@meta <- objList$meta
            
            return(x)
          })




#' Method setExpName
#' 
#' Set Experiment Name of a \code{CellMig} object.
#' 
#' @param x a \code{CellMig}-class object
#' @param ExpName string corresponding to the ExpName
#' 
#' @rdname setExpName
#' @exportMethod setExpName
setGeneric("setExpName", function(x, ExpName) {
  standardGeneric("setExpName")
})


#' @rdname setExpName
#' @aliases setExpName,CellMig,character-method
setMethod("setExpName", 
          signature(x="CellMig", 
                    ExpName = "character"),
          function(x, ExpName){
            
            x@meta$expName <- ExpName
            return(x)
          })



