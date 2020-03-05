##
##
## ~~ All Methods ~~
#
# 

# ~~~
# 1) trackedCells-methods
# ~~~

#' @import methods
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





# ~~~
# 2) CellMig-methods
# ~~~
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

