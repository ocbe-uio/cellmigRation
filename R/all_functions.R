##
##
## ~~ All f(x) ~~
#
#

#' Return the Next Odd Integer
#'
#' Returns the smallest odd number bigger than the number(s) provided as the argument
#'
#' @param x a vector of class numeric
#'
#' @return a vector of class integer
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::NextOdd(2:5)
#'
#' @keywords internal
NextOdd <- function(x) {
  y <- base::floor(x) + 1
  y <- base::ifelse(y %% 2 == 0, y + 1, y)
  return(y)
}


#' Shift Array Circularly
#'
#' Circularly shift the elements in an array by a user-defined number of positions.
#' This emulates the behavior of the corresponding Matlab Circhsift function.
#'
#' @param x a character, numeric, or logical vector with at least n + 1 elements
#' @param n an integer corresponding to the number of positions for the shift
#'
#' @return a vector corresponding to x (same size, same class), whose elements have been shifted
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::circshift(1:10, -2)
#'
#' @keywords internal
circshift <- function(x, n = 1) {

  n <- as.integer(n[1])
  nn <- abs(n)

  if(is.vector(x) || is.list(x)) {
    len <- length(x)
  } else if (is.data.frame(x) || is.matrix(x)) {
    len <- nrow(x)
  }

  if (len == 1)
    return(x)

  if (nn >= len)
    stop("Bad n!")

  nu_idx <- 1:len
  if(n > 0) {
    nu_idx <- c((length(nu_idx) - nn + 1):length(nu_idx), 1:(length(nu_idx) - nn))
  } else if (n < 0) {
    nu_idx <- c((nn + 1):length(nu_idx), 1:nn)
  }

  if(is.vector(x) || is.list(x)) {
    y <- x[nu_idx]
  } else if (is.data.frame(x) || is.matrix(x)) {
    y <- x[nu_idx,]
  }
  return(y)
}


#' Add Dimension to a Molten Data Frame
#'
#' Creates a new (molten) data matrix where all elements of y ar added to each row of x.
#' Each row in x is recycled for each element in y. Elements in y are added as
#' the first column in the returned matrix.
#'
#' @param x a matrix or data.frame with at least 1 row and 1 column.
#' @param y a vector with elements that will be added to x
#'
#' @return a matrix with an extra column as compared to x
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::AddDimension(x = cbind(1:4, 4:1), y = c(9, 7))
#'
#' @keywords internal
AddDimension <- function(x, y) {
  w0 <- lapply(1:nrow(x), function(j1) {
    tmp <- x[j1, ]
    w1 <- lapply(1:length(y), function(j2) {
      c(y[j2], tmp)
    })
    do.call(rbind, w1)
  })
  out <- do.call(rbind, w0)
  out <- as.matrix(out)
  rownames(out) <- NULL
  colnames(out) <- NULL
  return(out)
}


#' Make Hypercube
#'
#' Creates a Molten Hypercube with a user-defined number of dimensions.
#' The values supplied by the user
#' are used to fill each dimension. All possible combination of values are included in
#' the resulting hyper cube.
#'
#' @param vals vector of values used to fill the hyper cube
#' @param dims integer indicating the number of dimensions. The resulting molden
#' data frame will have a number of
#' columns equal to dims
#'
#' @return Matrix corresponding to a molten hyper cube. The number of columns is equal to dims;
#' the number of rows is equal to length(vals) ^ dims
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::MakeHypercube(1:3, 3)
#'
#' @keywords internal
MakeHypercube <- function(vals, dims) {
  xi <- as.matrix(cbind(vals))
  yi <- vals
  if (dims > 1){
    for(i in 1:(dims-1)) {
      xi <- AddDimension(x = xi, y = yi)
    }
  } else {
    return(NULL)
  }
  return(xi)
}


#' Clean And Reformat a Numeric Matrix
#'
#' Convert any matrix-lie object to a numeric Matrix, and coerces all the elements to integer.
#' Row names and column names are removed.
#'
#' @param x matrix or data.frame including numeric data (or data that can be coerced to integer)
#'
#' @return numeric matrix with all its elements coerced to integer
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' tmp <- data.frame(A = 1:4, B=c(3.1, 2.8, 3.3, 9.1), C = FALSE)
#' cellmigRation:::matfix(tmp)
#'
#' @keywords internal
matfix <- function(x) {
  xx <- as.data.frame(x)
  for ( j in 1:ncol(xx)) {
    xx[, j] <- as.integer(xx[, j])
  }
  colnames(xx) <- NULL
  rownames(xx) <- NULL
  return(as.matrix(xx))
}


#' Linear Convolution of a Numeric Matrix
#'
#' Performs a linear convoltion of a Numeric Matrix, using a user-supplied linear kernel.
#' The convolution can be executed in a column-wise fashion by setting the col.wise argument to TRUE.
#' Alternatively, the convolution is performed in a row-wise fashion.
#'
#' @param x numeric matrix that will be used as input for the convoltion;
#' this matrix typically corresponds to an image where signal (high values) indicates the
#' presence of a cell or a cell-like particle
#' @param krnl numeric vector corresponding to th kernel that will be used for the convolution.
#' Briefly, the kernel includes the weights that will be used to compute a weighted sum at each
#' position of the input numeric matrix
#' @param col.wise logical; shall the linear convolution be performed in a column-wise or
#' row-wise fashion
#'
#' @return Linearly convoluted numeric matrix. The resulting matrix has the same
#' dimensions of the inut matrix
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' graphics::par(mfrow = c(1, 2))
#' tmp <- sapply(1:12, function(i) { (6 + abs(i - 6)) * c(1:10, 10:1) })
#' cnv.tmp <- cellmigRation:::LinearConv2(tmp, c(-3, 0, 3))
#' graphics::image(tmp); graphics::image(cnv.tmp)
#' @importFrom graphics par image
#'
#' @keywords internal
LinearConv2 <- function(x, krnl, col.wise = TRUE)
{

  # Adjust based on col.wise
  if (col.wise) {
    xx <- t(x)
  } else {
    xx <- x
  }

  # Enlarge x based on kernel size
  ncl <- ncol(xx)
  tmp.i <- sapply(1:floor(length(krnl)/2), function(w) {xx[,1]})
  tmp.f <- sapply(1:floor(length(krnl)/2), function(w) {xx[,ncl]})
  X <- cbind(tmp.i, xx, tmp.f)

  # Proceed with convolution
  Y <- do.call(rbind, lapply(1:nrow(X), function(ri) {
    sapply(1:(ncol(X) - length(krnl) + 1), function(ci) {
      xcoord <- ci:(ci+length(krnl)-1)
      tmp <- X[ri, xcoord]
      as.numeric(rbind(tmp) %*% cbind(krnl))
    })
  }))

  if (col.wise)
    Y <- t(Y)

  return(Y)
}



#' Visualize a matrix image
#'
#' Shows an image representation of a numeric matrix. Typically, this is a non-negative numeric matrix,
#' where signal (high values) corresponds to the presence of cells, or cell-like particles
#'
#' @param img_mtx numeric matrix corresponding to a image
#' @param col character vector corresponding to a valid color palette
#' @param ... additional arguments will be passed to graphics::image()
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @keywords cellTracker
#'
#' @examples
#' x <- sapply(1:20, function(i) {runif(n = 20, min = 0, max = 10)})
#' cellmigRation:::VisualizeImg(x)
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image
#'
#' @keywords internal
#' @export
VisualizeImg <- function(img_mtx, col = NULL, ...)
{

  if(is.list(img_mtx)) {
    img_mtx <- img_mtx[[1]]
  }

  if (is.null(col)) {
    col <- grDevices::colorRampPalette(c("white", "blue4"))(100)
  }

  if(!is.matrix(img_mtx))
    stop("The IMG is not a matrix!")

  m <- nrow(img_mtx)
  n <- ncol(img_mtx)
  xx <- t(img_mtx[m:1, ])
  graphics::image(xx, col = col, ...)
}




#' Visualize Cells in an Image Stack
#'
#' Visualize objects that were identified as cells in a given image stack
#'
#' @param tc_obj a \code{trackedCells} object
#' @param stack index of the image stack to use
#' @param pnt.cex cex of the points drawn around cells
#' @param txt.cex cex of the text used to annotate cells
#' @param offset offset value for the annotation
#' @param main string used for the plot title, can be NULL= NULL
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("The following example will take few min to complete")
#' \dontrun{
#'   x <- cellmigRation::TrackCellsDataset
#'   x <- CellTracker(x, lnoise = 5, diameter = 16, threshold = 10)
#'   VisualizeStackCentroids(tc_obj = x, stack = 2)
#' }
#'
#'
#' @export
VisualizeStackCentroids <- function(tc_obj, stack = 1,
                                    pnt.cex = 1.2, txt.cex = 0.9,
                                    offset = 0.18, main = NULL) {

  b <- tc_obj@proc_images$images[[stack]]
  cnt <- tc_obj@centroids[[stack]]

  if(is.null(main)){
    main <- paste0("Stack num. ", stack)
  }

  VisualizeImg(img_mtx = b, las = 1, main = main)
  VisualizeCntr(centroids = cnt, width_px = ncol(b), height_px = nrow(b),
                pnt.cex = pnt.cex, txt.cex = txt.cex, offset = offset)

}

#' Visualize Centroids
#'
#' Annotates centroids over an image
#'
#' @param centroids centroid data.frame
#' @param width_px width of the image in pixels
#' @param height_px height of the image in pixels
#' @param pnt.cex cex of the point (circle) drawn around each cell
#' @param txt.cex cex of the text used to annotate the image
#' @param offset offset for the text annotations
#' @param col color of the points, e.g. "red2"
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x1 <- data.frame(row = c(50, 80, 20, 65, 99),
#'                  col = c(15, 25, 50, 65, 86))
#' plot(2, 2, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", las = 2)
#' cellmigRation:::VisualizeCntr(x1, width_px = 100, height_px = 100)
#'
#'
#' @importFrom graphics points text
#'
#' @keywords internal
VisualizeCntr <- function(centroids, width_px, height_px, pnt.cex = 1.2,
                          txt.cex = 0.9, offset = 0.18, col = "red2")
{
  cnt <- centroids
  points(x = ((cnt$col - 1) / (width_px - 1)),
         y = 1-((cnt$row - 1) / (height_px - 1)),
         cex = pnt.cex, col = col)

  points(x = ((cnt$col - 1) / (width_px - 1)),
         y = 1-((cnt$row - 1) / (height_px - 1)),
         cex = 1.2, col = col)

  text(x = ((cnt$col - 1) / (width_px - 1)),
       y = 1-((cnt$row - 1) / (height_px - 1)),
       labels = 1:nrow(cnt), font = 4,
       cex = txt.cex, col = col,
       pos = 4, offset = offset)

  #return()
}



#' Visualize Cell Tracks originating at an Image Stack
#'
#' Visualize Cell Tracks that originated at an Image Stack of interest
#'
#' @param tc_obj a trackedCells object
#' @param stack index of the stack
#' @param pnt.cex cex of the point drawn around each cell
#' @param lwd width of the lines visualizing cell tracks
#' @param col color of the points and the tracks, e.g.: "red2"
#' @param col.untracked color of the points that were not tracked further, e.g.: "gray45"
#' @param main string used as plot title, can be NULL
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("The following example will take few min to complete")
#' \dontrun{
#'   x <- cellmigRation::TrackCellsDataset
#'   x <- CellTracker(x, lnoise = 5, diameter = 16, threshold = 10)
#'   visualizeCellTracks(tc_obj = x, stack = 2)
#' }
#'
#'
#' @importFrom stats setNames
#'
#' @export
visualizeCellTracks <- function(tc_obj, stack = 1,
                                pnt.cex = 1.2, lwd = 1.6,
                                col = "red2", col.untracked = "gray45",
                                main = NULL) {

  if (is.null(main)) {
    main <- paste0("Tracks of Cells in Stack num. ", stack)
  }

  # Retrieve anr show image

  b <- tc_obj@proc_images$images[[stack]]
  VisualizeImg(img_mtx = b, las = 1, main = main)

  # Rerieve tracks / centroids
  #cnt <- tracked_cells$centroids[[stack]]

  cnt <- tc_obj@tracks
  cnt <- cnt[cnt[, 3]>= stack, ]
  cids_stack <- cnt[cnt[, 3] == stack, 4]

  cnt_plt <- cnt[cnt[, 4] %in% cids_stack, ]

  id_2pls <- unique(cnt_plt[duplicated(cnt_plt[,4]), 4])
  id_1shr <- unique(cnt_plt[!cnt_plt[, 4] %in% id_2pls, 4])

  if(length(id_1shr) > 0) {
    tmp_cnt <- cnt_plt[cnt_plt[, 4] %in% id_1shr, ]

    tmp_cnt <- tmp_cnt[tmp_cnt[, 3] == stack, ]
    tmp_cnt <- setNames(as.data.frame(tmp_cnt),
                        nm = colnames(tc_obj@centroids[[1]]))

    VisualizeCntr(centroids = tmp_cnt, width_px = ncol(b), height_px = nrow(b),
                  pnt.cex = pnt.cex, txt.cex = 0.00001, offset = 0.1, col = col.untracked)

  }

  if(length(id_2pls) > 0) {
    tmp_cnt <- cnt_plt[cnt_plt[, 4] %in% id_2pls, ]
    tmp_cnt <- setNames(as.data.frame(tmp_cnt),
                        nm = colnames(tc_obj@centroids[[1]]))

    # Use dedicated f(x)
    visualizeTrcks(tracks = tmp_cnt, width_px = ncol(b), height_px = nrow(b),
                   i.slice = stack, pnt.cex = pnt.cex, lwd = lwd, col = col)

  }

  # DOne , no return needed
  # return()
}

#' Visualize Cell Tracks
#'
#' Annotates an image with cell centroids by adding cell ROIs and drawing cell tracks
#'
#' @param tracks cell tracks
#' @param width_px width in pixels
#' @param height_px height in pixels
#' @param i.slice index of the stack slice to use
#' @param pnt.cex cex for the points (circles) drawn around the cells
#' @param lwd lwd of cell tracks
#' @param col color used for the cell tracks, .g. "red"
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/celltracker/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x1 <- data.frame(c(10, 30, 25, 55, 43, 39, 75, 72),
#'                  c(22, 28, 35, 24, 31, 39, 65, 73),
#'                  c( 1,  2,  3,  4,  5,  6,  7,  8),
#'                  c( 1,  1,  1,  1,  1,  1,  1,  1))
#' plot(2, 2, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", las = 2)
#' cellmigRation:::visualizeTrcks(x1, width_px = 100, height_px = 100)
#'
#'
#' @importFrom graphics points lines
#'
#' @keywords internal
visualizeTrcks <- function(tracks, width_px, height_px, i.slice = 1, pnt.cex = 1.2, lwd = 1.2, col = "red")
{
  allcnt <- setNames(as.data.frame(tracks),
                     nm = c("row", "col", "slice", "cell"))
  cnt <- allcnt[allcnt$slice == i.slice, ]
  cnt <- cnt[order(cnt$cell),]


  all_cells_slice <- cnt$cell
  cKeep <- sapply(all_cells_slice, function(jj) {
    sum(allcnt$cell == jj) > 1
  })

  # Cell outile

  graphics::points(x = ((cnt$col - 1) / (width_px - 1)),
                   y = 1-((cnt$row - 1) / (height_px - 1)),
                   cex = pnt.cex, col = ifelse(cKeep, col, "gray75"))

  for(j in sort(unique(cnt$cell))){
    TMP <- allcnt[allcnt$cell == j,]
    graphics::lines(x = ((TMP$col - 1) / (width_px - 1)),
                    y = 1-((TMP$row - 1) / (height_px - 1)),
                    lwd = lwd, col = col)

  }
  #return()
}




#' Import Image from TIFF
#'
#' Import a .tif stack containing fluorescently labeled point particles to be tracked
#'
#' @param tiff_file path to a TIFF file to be read in
#' @param experiment string, a label to describe the experiment (optional). Can be NULL
#' @param condition string, a label to describe the experimental condition (optional). Can be NULL
#' @param replicate string, a label to identify the replicate (optional). Can be NULL
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' \dontrun{
#' # Let `path/to/tiff_file.tiff` be the path to tiff file we want to import
#' x <- LoadTiff(tiff_file = "path/to/tiff_file.tiff")
#' }
#'
#' @importFrom tiff readTIFF
#'
#' @export
LoadTiff <- function(tiff_file, experiment = NULL, condition = NULL, replicate = NULL)
{
  myIMG <- suppressMessages(
    suppressWarnings(
      tiff::readTIFF(source = tiff_file,
                     native = FALSE,
                     all = TRUE,
                     info = TRUE,
                     as.is = TRUE)
    )
  )

  if (!is.list(myIMG))
    myIMG <- list(myIMG)

  if (is.null(experiment)) {
    experiment <- NA
  } else {
    experiment <- tryCatch(as.character(experiment[1]), error = function(e) NA)
  }

  if (is.null(replicate)) {
    replicate <- NA
  } else {
    replicate <- tryCatch(as.character(replicate[1]), error = function(e) NA)
  }

  if (is.null(condition)) {
    condition <- NA
  } else {
    condition <- tryCatch(as.character(condition[1]), error = function(e) NA)
  }

  # num of images
  NumberImages <- length(myIMG)

  # m means width... should be n of cols
  mImage <- ncol(myIMG[[1]])

  # n means height... should be n of rows
  nImage <- nrow(myIMG[[1]])

  # Get image INFO
  InfoImage <- try({lapply(myIMG, attributes)}, silent = TRUE)

  # Get image matrices
  FinalImage <- try({lapply(myIMG, function(x) {
    sapply(1:ncol(x), function(ii) {as.numeric(x[,ii])})
  })}, silent = TRUE)

  #return(list(images = FinalImage,
  #            dim = list(NumberImages  = NumberImages , width_m = mImage, height_n = nImage),
  #            attributes = InfoImage))
  img_list <-   list(images = FinalImage,
                     dim = list(NumberImages  = NumberImages , width_m = mImage, height_n = nImage),
                     attributes = InfoImage)

  Y <- new(Class = "trackedCells", img_list)

  # Attach labels
  Y@metadata <- list(tiff_file = sub("^.*[/]([^/]+$)", "\\1", tiff_file),
                     experiment = experiment,
                     condition = condition,
                     replicate = replicate)

  return(Y)
}



#' Validate Centroids
#'
#' Validate parameters used to identify cells in a image stack. A figure containing
#' current image frame with identified particles labeled with circles and numerical tags is generated.
#' This function is included for consistency and compatibility reasons
#' with the original fastTracks software (Matlab). Also, consider using
#' VisualizeStackCentroids() or visualizeCellTracks() instead.
#'
#' @param stack stack of images to be evaluated
#' @param slice index of the frame within the stack to be evaluated
#' @param lobject integer, length in pixels somewhat larger than a typical object (cell)
#' @param threshold the minimum brightness of a pixel that might be local maxima. NOTE:
#' Make it big and the code runs faster but you might miss some particles.
#' Make it small and you'll get everything and it'll be slow.
#' @param pnt.cex cex of the circle drawn around each cell
#' @param txt.cex cex of the text used for annotating cells
#' @param offset offset used for annotating cells
#'
#' @return data.frame of centroid positions
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("the following example may take up to one minute to run")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' y <- cellmigRation:::CentroidValidation(stack = x@@images, slice = 1,
#'                                         lobject = 6, threshold = 20)
#' y[1:10,]
#' }
#'
#' @importFrom graphics box axis points text
#'
#' @keywords internal
CentroidValidation <- function(stack, slice, lobject, threshold,
                               pnt.cex = 1.2, txt.cex = 0.85, offset = 0.18)
{
  a <- stack$images[[slice]]
  b <- bpass(image_array = a, lnoise = 1, lobject = lobject, threshold = threshold)
  pk = pkfnd(im = b, th = threshold, sz = lobject+1)
  cnt = cntrd(im = b, mx = pk, sz = lobject+1)

  VisualizeImg(b, axes = FALSE)
  graphics::box()
  my_xax <- ncol(b)
  my_yax <- nrow(b)

  my_xax <- unique(c(seq(1, my_xax, by = 100), my_xax))
  my_yax <- unique(c(seq(1, my_yax, by = 100), my_yax))

  axis(side = 1,at = ((my_xax - 1) / max(my_xax)), labels = my_xax)
  axis(side = 2,at = 1-((my_yax - 1) / max(my_yax)), labels = my_yax, las = 1)


  points(x = ((cnt$col - 1) / (ncol(b) - 1)),
         y = 1-((cnt$row - 1) / (nrow(b) - 1)),
         cex = pnt.cex, col = "red2")

  points(x = ((cnt$col - 1) / (ncol(b) - 1)),
         y = 1-((cnt$row - 1) / (nrow(b) - 1)),
         cex = 1.2, col = "red2")

  text(x = ((cnt$col - 1) / (ncol(b) - 1)),
       y = 1-((cnt$row - 1) / (nrow(b) - 1)),
       labels = 1:nrow(cnt), font = 4,
       cex = txt.cex, col = "red2",
       pos = 4, offset = offset)

  return(cnt)
}



#' Perform a bandpass by convolving with an appropriate kernel
#'
#' Implements a real-space bandpass filter that suppresses pixel noise and long-wavelength
#' image variations while retaining information of a characteristic size.
#' First, a lowpassed image is produced by convolving the original with a gaussian.
#' Next, a second lowpassed image is produced by convolving the original with a
#' boxcar function. By subtracting the boxcar version from the gaussian version,
#' we are using the boxcar version to perform a highpass.
#' This code 'bpass.pro' is copyright 1997, John C. Crocker and
#' David G. Grier.  It should be considered 'freeware'- and may be
#' distributed freely in its original form when properly attributed.
#'
#' @param image_array Numeric matrix corresponding to the image to be filtered
#' @param lnoise Characteristic lengthscale of noise in pixels.
#' @param lobject Integer length in pixels somewhat larger than a typical object
#' @param threshold By default, after the convolution, any negative pixels are reset to 0.
#' Threshold changes the threshhold for setting pixels to 0. Positive values may be useful for removing
#' stray noise or small particles.
#'
#' @return Numeric matrix corresponding to the filtered image
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cellmigRation::TrackCellsDataset@@images$images[[1]][380:480, 300:400]
#' y0 <- cellmigRation:::bpass(x0, lnoise = 5, lobject = 16, threshold = 10)
#' par(mfrow = c(1, 2))
#' image(x0); title("original")
#' image(y0); title("after bpass")
#'
#' @keywords internal
bpass <- function(image_array, lnoise, lobject = NULL, threshold)
{

  cstm_normalize <- function(x) { x/sum(x) }

  # Make kernel (linear)
  gaussian_kernel <- cstm_normalize(exp(-(seq(-2.5, 2.5, length.out = ((10 * lnoise) + 1))^2)))

  if (!is.null(lobject))
    boxcar_kernel <- cstm_normalize(rep(1, times = (2 * lobject) + 1))


  gconv <- LinearConv2(t(image_array), gaussian_kernel)
  gconv <- LinearConv2(t(gconv), gaussian_kernel)

  #VisualizeImg(image_array, col = colorRampPalette(c("white", "red2"))(100))
  #VisualizeImg(gconv)

  if (!is.null(lobject)) {
    bconv <- LinearConv2(t(image_array), boxcar_kernel)
    bconv <- LinearConv2(t(bconv), boxcar_kernel)
    filtered <- gconv - bconv
  } else {
    filtered <- gconv
  }

  # Zero out the values on the edges to signal that they're not useful.
  lzero <- max(lobject, ceiling(5*lnoise))

  filtered[1:(round(lzero)),] <- 0
  filtered[(nrow(filtered) - round(lzero) + 1):nrow(filtered),] <- 0

  filtered[, 1:(round(lzero))] <- 0
  filtered[, (ncol(filtered) - round(lzero) + 1):ncol(filtered)] <- 0

  # Zero all values below threshold
  filtered[filtered < threshold] <- 0

  return(filtered)
}



#' Calculates Centroids
#'
#' Calculates the centroid of bright spots to sub-pixel accuracy. Inspired by Grier &
#' Crocker's feature for IDL, but greatly simplified and optimized for matlab, and then
#' further ported to R. CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#'
#' @param im numeric matrix corresponding to the image to process
#' @param mx location of local maxima to pixel-levels accuracy
#' @param sz diameter of the window over which to average to calculate the centroid. should be big enough.
#' @param interactive numeric; if set to 1 (or any positive number), an image showing the
#' computed centroids will be visualized
#'
#' @return a data.frame with 4 columns, containing, x, y, brightness, and
#' the square of the radius of gyration for each cell.
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("the following example may take up to one minute to run")
#' \dontrun{
#' x0 <- cellmigRation::TrackCellsDataset@@images$images[[1]][100:500, 200:600]
#' b <- cellmigRation:::bpass(image_array = x0,
#'                            lnoise = 1,
#'                            lobject = 15,
#'                            threshold = 10)
#' pk <- cellmigRation:::pkfnd(b, 10, 17)
#' cnt <- cellmigRation:::cntrd(im = b, mx = pk, sz = 17)
#' cnt[1:10,]
#' }
#'
#' @importFrom graphics box axis title
#'
#' @keywords internal
cntrd <- function(im, mx, sz, interactive = NULL)
{
  # check interactive
  if(is.null(interactive))
    interactive <- 0

  # check sz
  if ((sz/2) == (floor(sz/2))) {
    sz <- sz + 1
    message("sz must be odd, like bpass")
    message(paste0("sz set to ", sz))
  }

  # check mx
  if (is.null(mx) || (!is.data.frame(mx)) || nrow(mx) == 0) {
    message('there were no positions inputted into cntrd. check your pkfnd theshold')
    return(NULL)
  }

  # Compute
  r <- (sz+1)/2

  # Create mask - window around trial location over which to calculate the centroid
  m <- 2*r
  x <- 0:(m-1)
  cent <- (m-1)/2
  x2 <- (x-cent) ^ 2
  dst <- do.call(rbind, lapply(1:m, function(i){
    sqrt((i-1-cent)^2+x2)
  }))

  ind <- dst < r

  msk <- sapply(1:ncol(ind), function(j) {as.numeric(ind[,j])})
  dst2 <- msk * (dst^2)
  ndst2 <- sum(dst2, na.rm = TRUE)

  nr <- nrow(im)
  nc <- ncol(im)

  # remove all potential locations within distance sz from edges of image
  ind <- mx$col > 1.5 * sz & mx$col < nc - 1.5*sz
  mx <- mx[ind, ]
  ind <- mx$row > (1.5*sz) & mx$row < nr - 1.5*sz
  mx <- mx[ind, ]

  nmx <- nrow(mx)

  # inside of the window, assign an x and y coordinate for each pixel
  xl <- do.call(rbind, lapply(1:(2*r), function(j) {
    (1:(2*r))
  }))
  yl <- t(xl)

  #loop through all of the candidate positions
  pts <- list()
  for (i in 1:nmx) {
    #create a small working array around each candidate location, and apply the window function
    tmp <- msk * im[(mx$row[i] - floor(r) + 1):(mx$row[i] + floor(r)),
                    (mx$col[i] - floor(r) + 1):(mx$col[i] + floor(r))]

    #calculate the total brightness
    norm <- sum(tmp, na.rm = TRUE)

    #calculate the weigthed average x location
    xavg <- sum(tmp * xl) / norm

    #calculate the weighted average y location
    yavg <- sum(tmp * yl) / norm

    #calculate the radius of gyration^2
    #rg=(sum(sum(tmp.*dst2))/ndst2);
    rg <- sum(tmp * dst2)/norm

    #concatenate it up
    pts[[length(pts) +1 ]] <- data.frame(row = mx$row[i]+yavg-r,
                                         col = mx$col[i] + xavg - r,
                                         norm = norm ,
                                         rg = rg)
    if (as.numeric(interactive) > 0) {
      VisualizeImg(img_mtx = tmp, axes = FALSE)
      box()
      axis(side = 1, at = seq(0, 1, length.out = ncol(tmp)),
           labels = (mx$row[i] - floor(r) + 1):(mx$row[i] + floor(r)))
      axis(side = 2, at = seq(0, 1, length.out = nrow(tmp)),
           labels = (mx$col[i] + floor(r)):(mx$col[i] - floor(r) + 1), las = 1)
      title(main = paste0("Cell number #", i), ylab = "y_pixel",
            xlab = "x_pixel", font = 2, cex = 0.9)

      # Wait for user input from keyboard
      readline("Press Enter for Next Cell...")
    }
  }
  pts <- do.call(rbind, pts)
  rownames(pts) <- NULL
  return(pts)
}





#' Find Signal Peaks
#'
#' Finds local maxima in an image to pixel level accuracy. This provides a rough guess
#' of particle centers to be used by cntrd(). Inspired by the lmx subroutine of
#' Grier and Crocker's. CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#'
#' @param im image to process, particle should be bright spots on dark background
#' with little noise ofen an bandpass filtered brightfield image
#' @param th the minimum brightness of a pixel that might be local maxima.
#' NOTE: Make it big and the code runs faster but you might miss some particles.
#' Make it small and you'll get everything and it'll be slow.
#' @param sz if your data is noisy, (e.g. a single particle has multiple local maxima),
#' then set this optional keyword to a value slightly larger than the diameter of your blob.
#' If multiple peaks are found withing a radius of sz/2 then the code will keep only the brightest.
#' Also gets rid of all peaks within sz of boundary
#'
#' @return a numeric data.frame with two columns, with the coordinates of local maxima
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("the following example may take up to one minute to run")
#' \dontrun{
#' x0 <- cellmigRation::TrackCellsDataset@@images$images[[1]][100:500, 200:600]
#' b <- cellmigRation:::bpass(image_array = x0,
#'                            lnoise = 1,
#'                            lobject = 15,
#'                            threshold = 10)
#' pk <- cellmigRation:::pkfnd(b, 10, 17)
#' cnt <- cellmigRation:::cntrd(im = b, mx = pk, sz = 17)
#' cnt[1:10,]
#' }
#'
#' @keywords internal
pkfnd <- function(im, th, sz=NULL)
{

  # # nested f(x)
  #my_melt <- function(data, varnames = NULL) {
  #
  #  out.list <- list()
  #  for(ci in 1:ncol(data)) {
  #    for (ri in 1:nrow(data)) {
  #      out.list[[length(out.list) + 1]] <-
  #        data.frame(Var1 = ri, Var2 = ci, value = data[ri, ci],
  #                   stringsAsFactors = FALSE)
  #    }
  #  }
  #  OUT <- do.call(rbind, out.list)
  #
  #  if (!is.null(varnames)) {
  #    colnames(OUT)[1:2] <- varnames
  #  }
  #  return(OUT)
  #}

  # find all the pixels above threshold
  ind <- im > th
  nr <- nrow(im)
  nc <- ncol(im)

  # melt to have a list of points above threshold
  ind2 <- reshape2::melt(data = ind, varnames = c("row", "col"))
  ind2 <- ind2[ind2$value, ]

  # check each pixel above threshold to see if it's brighter than it's neighbors
  # THERE'S GOT TO BE A FASTER WAY OF DOING THIS.  I'M CHECKING SOME MULTIPLE TIMES,
  # BUT THIS DOESN'T SEEM THAT SLOW COMPARED TO THE OTHER ROUTINES, ANYWAY.
  keep <- list()
  for(i in 1:nrow(ind2)) {
    ri <- ind2$row[i]
    ci <- ind2$col[i]

    if (ri>1 & ri<nr & ci>1 & ci<nc) {
      z1 <- im[ri, ci]
      z2 <- as.numeric(im[(ri-1):(ri+1), (ci-1):(ci+1)])

      if (sum(z1 < z2, na.rm = TRUE) == 0) {
        keep[[length(keep) + 1]] <- i
      }
    }
  }

  # Next step
  npks <- length(keep)
  if(npks > 0) {
    keep <- do.call(c, keep)
    mx <- ind2[keep, ]
  } else {
    return(NULL)
  }

  # if size is specified, then get ride of pks within size of boundary (i.e., a margin from image edges)
  if (!is.null(sz) & npks>0) {
    # throw out all pks within sz of boundary;
    keep <- mx$row > sz & mx$row < (nr - sz + 1) & mx$col > sz & mx$col < (nc - sz + 1)
    mx<-mx[keep,]
  }

  # prevent from finding peaks within size of each other
  npks <- nrow(mx)
  if (!is.null(sz) & npks > 1) {
    # CREATE AN IMAGE WITH ONLY PEAKS
    mask <- matrix(FALSE, nrow = nrow(im), ncol = ncol(im))
    for(i in 1:nrow(mx)) {
      mask[mx$row[i], mx$col[i]] <- TRUE
    }

    tmp <- matrix(0, nrow=nrow(im), ncol = ncol(im))
    tmp[mask] <- im[mask]

    # LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK THE BRIGHTEST
    for (i in 1:nrow(mx)) {
      astep <- floor(sz/2)
      roi <- tmp[(mx$row[i] - astep):(mx$row[i] + astep), (mx$col[i] - astep):(mx$col[i] + astep)]
      imax <- which.max(roi)
      chkrow <- imax %% nrow(roi)
      myrow <- ifelse(chkrow == 0, nrow(roi), chkrow)
      mycol <- ifelse(chkrow == 0, floor(imax/nrow(roi)), floor(imax/nrow(roi)) + 1)
      mv <- roi[myrow, mycol]

      tmp[(mx$row[i] - astep):(mx$row[i] + astep), (mx$col[i] - astep):(mx$col[i] + astep)] <- 0
      tmp[(mx$row[i] - astep + myrow - 1), (mx$col[i] - astep + mycol - 1)] <- mv
    }

    ind <- tmp > th
    nr <- nrow(tmp)
    nc <- ncol(tmp)

    # melt to have a list of points above threshold
    ind.f <- reshape2::melt(data = ind, varnames = c("row", "col"))
    ind.f <- ind.f[ind.f$value, 1:2]
    rownames(ind.f) <- NULL

    return(ind.f)
  } else {
    return(NULL)
  }
}


#' Build a Centroid Array
#'
#' Create an array containing centroid data for particles identified in each frame
#' of the imported TIFF image stack
#'
#' @param stack 3D matrix loaded to workspace from .tif stack
#' @param lobject Integer length in pixels somewhat larger than a typical object
#' @param threshold the minimum brightness of a pixel that might be local maxima
#'
#' @return data.frame of centroids, with 4 columns corresponding to x-position of centroid,
#' y-postion of centroid, brightness, and square of the radius of gyration
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("the following example may take up to several minutes to run")
#' \dontrun{
#' x0 <- cellmigRation::TrackCellsDataset@@images
#' y0 <- cellmigRation:::CentroidArray(x0, 16, 10)
#' class(y0)
#' y0[[1]][1:10,]
#' }
#'
#'
#' @keywords internal
CentroidArray <- function(stack, lobject, threshold)
{
  #determine the number of slices within the stack
  m <- stack$dim$width_m
  n <- stack$dim$height_n
  p <- stack$dim$NumberImages

  centroid <- list()
  for(i in 1:p) {
    a <-  stack$images[[i]]
    b <- bpass(image_array = a,
               lnoise = 1,
               lobject = lobject,
               threshold = quantile(a, 0.25)) # maybe set to 0 or to threshold

    pk <- pkfnd(b, threshold, lobject+1)
    cnt <- cntrd(im = b, mx = pk, sz = lobject + 1)

    if(is.null(cnt) || nrow(cnt) < 1) {
      message(paste0('No centroids detectd in frame ', i, '...
                     \nCheck nuclei validation settings for this frame.'))
    }
    centroid[[length(centroid) + 1]] <- cnt
  }
  return(centroid)
}


#' Detect Linear Paricle Diameters
#'
#' Estimates the diameters of particles in a numeric or logical vector
#'
#'
#' @param x numeric or logical vector
#'
#' @return data.frame including two columns: MPOS indicates the centroid position of a particle,
#' and LEN indicates the diameter size
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::DetectRadii(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1))
#'
#' @keywords internal
DetectRadii <- function(x) {

  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]

  if(length(table(x)) > 2) {
    my.mean <- mean(x, na.rm = TRUE)
    x[x >= my.mean] <- 1
    x[x < my.mean] <- 0
  }

  radii <- list()
  xx <- which(x == 1)

  if (length(xx) > 1) {
    LN <- 1
    p0 <- xx[1]
    p1 <- xx[1]

    for (j in 2:length(xx)) {
      if (xx[j] == (xx[(j-1)] + 1)) {
        LN <- LN + 1
        p1 <- xx[j]
        if (j == length(xx)) {
          yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
          radii[[length(radii) + 1]]  <- yy
        }
      } else {
        yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
        radii[[length(radii) + 1]]  <- yy
        p0 <- xx[j]
        p1 <- xx[j]
        LN <- 1
      }
    }

  } else if (length(xx) == 1) {

    yy <- data.frame(MPOS = xx, LEN = 1)
    radii[[length(radii) + 1]]  <- yy
  }


  if (length(radii) > 0 ) {
    radii <- do.call(rbind, radii)
  } else {
    radii <- NULL
  }

  return(radii)
}


#' Detect Paricle Diameters in a Numeric matrix
#'
#' Estimates the diameters of particles in a numeric matrix
#'
#'
#' @param x numeric matrix corresponding to a digital image
#' @param px.margin integer, number of pixels used as margin while searching/filtering for neighboring particles
#' @param min.px.diam integer, minimum diameter of a particle (cell).
#' Particles with a diameter smaller than min.px.diam are discarded
#' @param quantile.val numeric, must be bigger than 0 and smaller than 1.
#' Quantile for discriminating signal and background; only pixels with intensity higher than the corresponding
#' quantile will count as signal while estimating particle diameters
#' @param plot logial, shall a histogram of the distribution of diameters be shown
#'
#' @return list including summary stats and data about the particles found in the image
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' a <- cbind(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1),
#'            c(1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
#'            c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'            c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
#'            c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0))
#' graphics::image(a)
#' b <- EstimateDiameterRange(a, min.px.diam = 2)
#' print(b$estim.cell.num)
#' print(b$raw)
#'
#' @importFrom graphics image hist
#'
#' @export
EstimateDiameterRange <- function(x, px.margin = 2,
                                  min.px.diam = 5,
                                  quantile.val = 0.99,
                                  plot = TRUE) {

  QNTS <- as.numeric(quantile(x, probs = quantile.val[1]))

  # Adjust if quantile.val is too low (few cells)
  tmp.xx <- as.numeric(x)
  max.sig <- max(tmp.xx, na.rm = TRUE)
  min.sig <- min(tmp.xx, na.rm = TRUE)
  if (QNTS == min.sig && max.sig > min.sig) {
    QNTS <- mean(c(min.sig, min(tmp.xx[tmp.xx > min.sig], na.rm = TRUE)), na.rm = TRUE)
  }

  B <- x
  B[B < QNTS] <- 0
  B[B >= QNTS] <- 1

  rdds <- do.call(rbind, lapply(1:ncol(B), function(ii) {

    out <- DetectRadii(B[,ii])
    if (!is.null(out)) {
      data.frame(RPOS = out$MPOS, CPOS = ii, LEN = out$LEN)
    }
  }))
  rdds$KEEP <- TRUE

  for (j in 1:nrow(rdds)) {
    if (rdds$KEEP[j]){

      tdm <- ( 2 * px.margin)  + rdds$LEN[j]
      ROWmin <- rdds$RPOS[j] - (0.5 * tdm)
      ROWmax <- rdds$RPOS[j] + (0.5 * tdm)
      COLmin <- rdds$CPOS[j] - (0.5 * tdm)
      COLmax <- rdds$CPOS[j] + (0.5 * tdm)

      keep <- rdds$RPOS >= ROWmin & rdds$RPOS <= ROWmax &
        rdds$CPOS >= COLmin & rdds$CPOS <= COLmax & rdds$KEEP
      keep <- which(keep)
      keep <- keep[keep != j]
      if (length(keep) > 0) {
        curVal <- rdds$LEN[j]
        allValz <- rdds$LEN[keep]

        if (sum(curVal > allValz) == length(allValz)) {
          rdds$KEEP[keep] <- FALSE
        } else {
          rdds$KEEP[j] <- FALSE
        }
      }
    }
  }

  FINL <- rdds[rdds$KEEP,]
  FINL <- FINL[FINL[, "LEN"] >= min.px.diam, ]

  yy <- list(estim.cell.num = sum(FINL$KEEP),
             q50.diam = median(FINL$LEN, na.rm = TRUE),
             q75.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.75)),
             q90.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.90)),
             q95.diam = as.numeric(quantile(FINL$LEN, na.rm = TRUE, probs = 0.95)),
             raw = FINL)

  if (plot) {
    try(hist(FINL$LEN, breaks = seq(min(FINL$LEN, na.rm = TRUE),
                                    max(FINL$LEN, na.rm = TRUE), length.out = 20),
             xlab = "Particle Diameter", las = 1, main = "Diam. Distribution",
             col = "aquamarine3"), silent = TRUE); box()
  }

  return(yy)
}


#' Track cells
#'
#' Constructs n-dimensional trajectories from a scrambled list of particle
#' coordinates determined at discrete times (e.g. in consecutive image frames)
#'
#'
#' @param xyzs an array listing the scrambled coordinates and data of the
#' different particles at different times
#' @param maxdisp an estimate of the maximum distance that a particle would
#' move in a single time interval
#' @param params a list containing a few tracking parameters that are needed
#' for the analysis
#'
#'
#' @return data.frame including cell tracks data
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- data.frame(row = c(1, 30, 50, 5, 35, 55, 6, 56, 7, 58),
#'                  col = c(1, 30, 50, 5, 35, 55, 6, 56, 7, 58),
#'                  tau = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4))
#' cellmigRation:::track(x0, maxdisp = 10, params = NULL)
#'
#'
#' @keywords internal
track <- function(xyzs, maxdisp, params)
{
  #% ; maxdisp: an estimate of the maximum distance that a particle
  #% ;     would move in a single time interval.(see Restrictions)
  #%  OPTIONAL INPUT:
  #%   param:  a structure containing a few tracking parameters that are
  #%       needed for many applications.  If param is not included in the
  #%       function call, then default values are used.  If you set one value
  #%       make sure you set them all:
  #% ;         param.mem: this is the number of time steps that a particle can be
  #% ;             'lost' and then recovered again.  If the particle reappears
  #% ;             after this number of frames has elapsed, it will be
  #% ;             tracked as a new particle. The default setting is zero.
  #% ;             this is useful if particles occasionally 'drop out' of
  #% ;             the data.
  #% ;         param.dim: if the user would like to unscramble non-coordinate data
  #% ;             for the particles (e.g. apparent radius of gyration for
  #% ;             the particle images), then positionlist should
  #% ;             contain the position data in positionlist(0:param.dim-1,*)
  #% ;             and the extra data in positionlist(param.dim:d-1,*). It is then
  #% ;             necessary to set dim equal to the dimensionality of the
  #% ;             coordinate data to so that the track knows to ignore the
  #% ;             non-coordinate data in the construction of the
  #% ;             trajectories. The default value is two.
  #% ;         param.good: set this keyword to eliminate all trajectories with
  #% ;             fewer than param.good valid positions.  This is useful
  #% ;             for eliminating very short, mostly 'lost' trajectories
  #% ;             due to blinking 'noise' particles in the data stream.
  #%;          param.quiet: set this keyword to 1 if you don't want any text
  #% ; OUTPUTS:
  #  % ; result:  a list containing the original data rows sorted
  #% ;     into a series of trajectories.  To the original input
  #% ;     data structure there is appended an additional column
  #% ;     containing a unique 'id number' for each identified
  #% ;     particle trajectory.  The result array is sorted so
  #% ;     rows with corresponding id numbers are in contiguous
  #% ;     blocks, with the time variable a monotonically
  #% ;     increasing function inside each block.  For example:
  #  % ;
  #% ;
  #% ;     NB: for t=1 in the example above, one particle temporarily
  #% ;     vanished.  As a result, the trajectory id=1 has one time
  #% ;     missing, i.e. particle loss can cause time gaps to occur
  #% ;     in the corresponding trajectory list. In contrast:
  #  % ;
  #% ;     IDL> res = track(pos,5)
  #% ;
  #% ;     track will return the result 'res'
  #% ;         (x)      (y)      (t)          (id)
  #% ;     res = 15.1000      22.6000      0.00000      0.00000
  #% ;                   3.60000      5.00000      0.00000      1.00000
  #% ;               4.10000      5.50000      1.00000      1.00000
  #% ;               6.20000      4.30000      2.00000      1.00000
  #% ;               15.9000      20.7000      2.00000      2.00000
  #% ;
  #% ;     where the reappeared 'particle' will be labelled as new
  #% ;     rather than as a continuation of an old particle since
  #% ;     mem=0.  It is up to the user to decide what setting of
  #% ;     'mem' will yeild the highest fidelity .
  #% ;
  #% ; SIDE EFFECTS:
  #  % ; Produces informational messages.  Can be memory intensive for
  #% ; extremely large data sets.
  #% ; RESTRICTIONS:
  #  % ; maxdisp should be set to a value somewhat less than the mean
  #% ; spacing between the particles. As maxdisp approaches the mean
  #% ; spacing the runtime will increase significantly. The function
  #% ; will produce an error message: "Excessive Combinatorics!" if
  #% ; the run time would be too long, and the user should respond
  #% ; by re-executing the function with a smaller value of maxdisp.
  #% ; Obviously, if the particles being tracked are frequently moving
  #% ; as much as their mean separation in a single time step, this
  #% ; function will not return acceptable trajectories.
  #% ; PROCEDURE:
  #  % ; Given the positions for n particles at time t(i), and m possible
  #% ; new positions at time t(i+1), this function considers all possible
  #% ; identifications of the n old positions with the m new positions,
  #% ; and chooses that identification which results in the minimal total
  #% ; squared displacement. Those identifications which don't associate
  #% ; a new position within maxdisp of an old position ( particle loss )
  #% ; penalize the total squared displacement by maxdisp^2. For non-
  #% ; interacting Brownian particles with the same diffusivity, this
  #% ; algorithm will produce the most probable set of identifications
  #% ; ( provided maxdisp >> RMS displacement between frames ).
  #% ; In practice it works reasonably well for systems with oscillatory,
  #% ; ballistic, correlated and random hopping motion, so long as single
  #% ; time step displacements are reasonably small.  NB: multidimensional
  #% ; functionality is intended to facilitate tracking when additional
  #% ; information regarding target identity is available (e.g. size or
  #% ; color).  At present, this information should be rescaled by the
  #% ; user to have a comparable or smaller (measurement) variance than
  #% ; the spatial displacements.
  #% ;

  # Initialize and stuff
  warn_log <- list()

  warn_message <- function(warn_log, quiet = FALSE) {

    warn_cycles <- NULL

    if (is.list(warn_log) && length(warn_log) > 0) {

      warn_cycles <- sort(unique(do.call(c, warn_log)))

      if (!quiet) {
        message(paste0("Difficult combinatorics encountered while processing slide(s): ",
                       paste(warn_cycles, collapse = ", "), "."))
      }
    }
    return(warn_cycles)
  }



  dd <- ncol(xyzs)

  # use default parameters if none given
  # if nargin==2
  # default values
  memory_b <- 0
  goodenough <- 0
  dim <- dd - 1
  quiet <- FALSE
  force_exec <- FALSE


  if(!is.null(params)) {
    if(is.list(params)) {

      if(!is.null(params$memory_b) && is.numeric(params$memory_b)) {
        memory_b <- params$memory_b
      }

      if(!is.null(params$goodenough) && is.numeric(params$goodenough)) {
        goodenough <- params$goodenough
      }

      if(!is.null(params$dim) && is.numeric(params$dim)) {
        dim <- params$dim
      }

      if(!is.null(params$quiet) && is.logical(params$quiet)) {
        quiet <- params$quiet
      }


      if(!is.null(params$force_exec) && is.logical(params$force_exec)) {
        force_exec <- params$force_exec
      }

    }
  }

  # % checking the input time vector
  # THis should be monotonically not-decreasing and not identical
  tau <- xyzs[, dd]
  st <- tau[2:length(tau)] - tau[1:(length(tau) - 1)]

  if (sum(st < 0) > 0) {
    message("", appendLF = TRUE)
    message("The time vector (tau) is not ordered")
    return(NULL)
  }

  if (length(unique(tau)) == 1) {
    message("", appendLF = TRUE)
    message('All positions are at the same time... go back!')
    return(NULL)
  }

  #--remove if useless
  info <- 1
  w <- which(st > 0)
  z <- length(w)
  z <- z + 1

  # % partitioning the data with unique times
  # the first two lines were skipped in the original file
  # they are included here for completeness
  #res = unq(t);
  # implanting unq directly

  indices <- which(tau - circshift(tau, -1) != 0)
  count <- length(indices)

  if (count > 0) {
    res <- indices
  } else{
    res = length(tau)-1
  }

  res <- c(1, res, length(tau))
  ngood <- res[2] - res[1] + 1
  eyes <- 1:ngood
  pos <- xyzs[eyes, 1:dim]
  istart <- 2
  n <- ngood;

  zspan <- 50;
  if (n > 200) {
    zspan <- 20
  }

  if (n > 500){
    zspan <- 10
  }

  # initialize a matrix with -1
  resx <- matrix((-1), nrow = zspan, ncol = n)

  # initialize a second matrix with -1
  bigresx <- matrix((-1), nrow = z, ncol = n)
  mem <- matrix(0, nrow = n, ncol = 1)
  #%  whos resx
  #%  whos bigresx
  uniqid <- 1:n;
  maxid <- n;

  # initialize olis
  #olist <- data.frame(x= 0, y = 0)
  #olist <- c(0,0)
  olist <- list()

  if (goodenough > 0) {
    dumphash <- matrix(0, nrow = n, ncol = 1)
    nvalid <- matrix(1, nrow = n, ncol = 1)
  }

  #%  whos eyes;
  resx[1,] <- eyes

  #% setting up constants
  maxdisq <- maxdisp^2

  #% (Little) John calls this the setup for "fancy code" ???
  # Robin replies: Fancy? Where? You got to be kidding, man!!!

  notnsqrd <- (sqrt(n*ngood) > 200) && (dim < 7)

  if (notnsqrd) {
    #%;   construct the vertices of a 3x3x3... d-dimensional hypercube
    numbs <- 0:2
    cube <- MakeHypercube(vals = numbs, dims = dim)

    #%   calculate a blocksize which may be greater than maxdisp, but which
    #%   keeps nblocks reasonably small.
    volume <- 1
    for (d in 0:(dim-1)) {
      minn <- min(xyzs[w, (d+1)])
      maxx = max(xyzs[w, (d+1)])
      volume <- volume * (maxx-minn)
    }

    # volume;
    blocksize <- max( c(maxdisp,((volume)/(20*ngood))^(1.0/dim)) )
  }


  ### %   Start the main loop over the frames.
  for (i in istart:z){

    #message(paste0("i=", i))
    ispan <- ((i-1) %% zspan) + 1
    # %disp(ispan)
    # % get new particle positions
    m <- res[(i+1)] - res[(i)]
    # res[i]
    eyes <- 1:m
    eyes <- eyes + res[i]

    if (m > 0) {
      xyi <- xyzs[eyes, 1:dim]
      found <- matrix(0, nrow = m, ncol = 1)

      # % THE TRIVIAL BOND CODE BEGINS
      if (notnsqrd) {

        # %Use the raster metric code to do trivial bonds

        #% construct "s", a one dimensional parameterization of the space
        #% which consists of the d-dimensional raster scan of the volume.)

        abi <- matfix(xyi/blocksize)
        abpos <- matfix(pos/blocksize)
        si <- matrix(0, nrow = m, ncol = 1)
        spos <- matrix(0, nrow = n, ncol = 1)
        dimm <- matrix(0, nrow=dim, ncol=1)
        coff <- 1

        for (j in 1:dim){
          minn <- min(c(as.numeric(abi[,j]),
                        as.numeric(abpos[, j])), na.rm = TRUE)
          maxx <- max(c(as.numeric(abi[,j]),
                        as.numeric(abpos[,j])), na.rm = TRUE)
          abi[, j] <- abi[, j] - minn
          abpos[,j] <- abpos[,j] - minn
          dimm[j,1] <- maxx-minn + 1
          si <- si + abi[,j] * coff
          spos <- spos + abpos[,j]*coff
          coff <- dimm[j,1]*coff
        }
        nblocks <- coff
        #% trim down (intersect) the hypercube if its too big to fit in the
        #% particle volume. (i.e. if dimm(j) lt 3)

        cub <- cube
        deg <- which(dimm[,1] < 3)
        if (length(deg) > 0) {
          for (j in 0:(length(deg)-1)){
            cub <- cub[which(cub[, deg[j+1]] < dimm[deg[j+1],1]) ,]
          }
        }

        # % calculate the "s" coordinates of hypercube (with a corner @ the origin)
        scube <- matrix(0, nrow = nrow(cub), ncol=1)
        coff <- 1
        for (j in 1:dim){
          scube <- scube + (cub[,j] * coff)
          coff <- coff*dimm[j, 1]
        }

        # % shift the hypercube "s" coordinates to be centered around the origin
        coff <- 1
        for (j in 1:dim){
          if (dimm[j, 1] > 3) {
            scube <- scube - coff
          }
          coff <- dimm[j, 1] * coff
        }
        scube <- (scube + nblocks) %% nblocks

        # get the sorting for the particles by their "s" positions.
        ed <- sort(si)
        isort <- order(si)

        #% make a hash table which will allow us to know which new particles
        #% are at a given si.
        strt <- matrix((-1), nrow = nblocks, ncol = 1)
        fnsh <- matrix(0, nrow = nblocks, ncol = 1)
        h <- which(si == 0)
        lh <- length(h)
        if (lh > 0) {
          si[h] <- 1
        }

        for (j in 1:m){
          if (strt[si[isort[j]], 1] == (-1)){
            strt[si[isort[j]],1] <- j
            fnsh[si[isort[j]], 1] <- j
          } else {
            fnsh[si[isort[j]], 1] <- j
          }
        }
        if (lh > 0) {
          si[h] <- 0
        }


        coltot <- matrix(0, nrow = m, ncol = 1)
        rowtot <- matrix(0, nrow = n, ncol = 1)
        which1 <- matrix(0, nrow = n, ncol = 1)

        for (j in 1:n){

          map <- matfix(-1)

          scub_spos <- scube + spos[j];
          s <- scub_spos %% nblocks
          whzero <- which(s == 0 )
          if (length(whzero > 0)){
            nfk <- which(s !=0 )
            s <- s[nfk]
          }

          w <- which(strt[s, 1] != (-1))

          ngood <- length(w)
          ltmax <- 0
          if (ngood != 0){

            s <- s[w]
            for (k in 1:ngood){
              map = c(map, isort[strt[s[k]]:fnsh[s[k]]])
            }
            map <- map[2:length(map)]
            #%                     if length(map) == 2
            #%                         if (map(1) - map(2)) == 0
            #%                             map = unique(map);
            #%                          end
            #%                     end
            #%   map = map(umap);
            #%end
            #% find those trival bonds
            distq <- matrix(0, nrow=length(map), ncol=1)
            for (d in 1:dim){
              distq <- distq + (xyi[map,d] - pos[j,d])^2
            }
            ltmax <- distq < maxdisq

            rowtot[j, 1] <- sum(ltmax)

            if (rowtot[j] >= 1){
              w <- which(ltmax == 1)
              coltot[map[w], 1] <- coltot[ map[w], 1] +1
              which1[j, 1] <- map[ w[1]]
            }
          }
        }


        ntrk <- matfix(n - sum(rowtot == 0))

        w <- which(rowtot == 1)
        ngood <- length(w)


        if (ngood != 0) {
          #ww <- which(coltot( which1[w] ) == 1);
          ww <- which(coltot[which1[w]] == 1)

          ngood <- length(ww)
          if (ngood != 0){
            # %disp(size(w(ww)))
            resx[ispan, w[ww]] <- eyes[which1[w[ww]]]
            found[which1[w[ww]]] <- 1
            rowtot[w[ww]] = 0;
            coltot[which1[w[ww]]] <- 0
          }
        }

        labely <- which(rowtot > 0)
        ngood <- length(labely)
        if (ngood != 0){
          labelx <- which(coltot > 0)

          nontrivial <- 1
        } else {
          nontrivial <- 0
        }

      } else {

        # % THE TRIVIAL BOND {else} block CODE BEGINS
        #%   or: Use simple N^2 time routine to calculate trivial bonds

        #% let's try a nice, loopless way!
        #% don't bother tracking perm. lost guys.
        wh <- which(pos[,1] >= 0)
        ntrack <- length(wh)
        if (ntrack == 0){
          message('There are no valid particles to track!')
          break
        }

        # yma initialization was added
        xmat <- matrix(0, nrow = ntrack, ncol = m)
        ymat <- matrix(0, nrow = ntrack, ncol = m)

        count <- 0
        for (kk in 1:ntrack) {
          for (ll in 1:m) {
            xmat[kk,ll] <- count
            count <- count+1
          }
        }
        count <- 0

        # if there are not enough cols or rows, add them and set to 0
        if (nrow(ymat) < m) {
          TMP <- matrix(0, nrow = (m - nrow(ymat)), ncol = ncol(ymat))
          ymat <- rbind(ymat, TMP)
        }
        if (ncol(ymat) < ntrack) {
          TMP <- matrix(0, nrow = nrow(ymat), ncol = (ntrack - ncol(ymat)))
          ymat <- cbind(ymat, TMP)
        }

        for (kk in 1:m) {
          for (ll in 1:ntrack) {
            ymat[kk,ll] <- count
            count <- count+1
          }
        }

        xmat <- (xmat %% m) + 1
        ymat <- t((ymat %% ntrack) +1)
        lenxn <- nrow(xmat)
        lenxm <- ncol(xmat)
        #%            whos ymat
        #%            whos xmat
        #%            disp(m)

        for (d in 1:dim) {
          x <- xyi[,d]
          y <- pos[wh,d]

          xm <- sapply(1:ncol(xmat), function(jj) {
            tcljj <- xmat[, jj]
            x[tcljj]
          })

          #ym <- y[ymat[1:lenxn, 1:lenxm]]
          tmpymat <- ymat[1:lenxn, 1:lenxm]
          ym <- sapply(1:ncol(tmpymat), function(jj) {
            tcljj <- tmpymat[, jj]
            y[tcljj]
          })

          if (nrow(xm) != nrow(ym) || ncol(xm) != ncol(ym)) {
            xm <- t(xm)
          }

          if (d == 1) {
            dq <- (xm -ym)^2
            #%dq = (x(xmat)-y(ymat(1:lenxn,1:lenxm))).^2;
          } else {
            dq <- dq + (xm-ym)^2
            #%dq = dq + (x(xmat)-y(ymat(1:lenxn,1:lenxm)) ).^2;
          }
        }

        ltmax <- 1 * (dq < maxdisq)

        #% figure out which trivial bonds go with which

        rowtot <- matrix(0, nrow = n, ncol = 1)
        rowtot[wh, 1] <- apply(ltmax, 1, sum)

        if (ntrack > 1) {
          coltot <- apply(ltmax, 2, sum, na.rm = TRUE)
        } else {
          coltot <- ltmax
        }
        which1 <- matrix(0, nrow = n, ncol = 1)
        for (j in 1:ntrack) {
          mx  <- max(ltmax[j, ], na.rm = TRUE)
          w <- which.max(ltmax[j, ])
          which1[wh[j]] <- w
        }

        ntrk <- matfix( n - sum(rowtot == 0))
        w <- which( rowtot == 1)
        ngood <- length(w)
        if (ngood != 0) {
          ww <- which(coltot[which1[w]] == 1)
          ngood <- length(ww)
          if (ngood != 0) {
            resx[ ispan, w[ww] ] <- eyes[ which1[w[ww]]]
            found[which1[ w[ww]]] <- 1
            rowtot[w[ww]] <- 0
            coltot[which1[w[ww]]] <- 0
          }
        }

        labely <- which(rowtot > 0)
        ngood <- length(labely)

        if (ngood != 0) {
          labelx <- which(coltot > 0)
          nontrivial <- 1
        } else {
          nontrivial <- 0
        }
      }

      # %THE TRIVIAL BOND CODE ENDS

      if (nontrivial == 1){

        xdim <- length(labelx)
        ydim <- length(labely)

        #%  make a list of the non-trivial bonds

        bonds <- list()
        bondlen <- list()

        for (j in 1:ydim) {
          distq <- matrix(0, nrow = xdim, ncol = 1)

          for (d in 1:dim) {
            #%distq
            distq <- distq + cbind((xyi[labelx,d] - pos[labely[j],d])^2)
            #%distq
          }

          w <- which(distq < maxdisq) - 1
          ngood <- length(w)
          newb <- rbind(w, rep(j, times = ngood))

          bonds[[(length(bonds) + 1)]] <- t(newb)
          bondlen[[(length(bondlen) + 1)]] = distq[w + 1]
        }

        bonds <- do.call(rbind, bonds)
        bondlen <- do.call(c, bondlen)

        numbonds <- length(bonds[,1])
        mbonds <- bonds;
        #max([xdim,ydim]);


        if (max(c(xdim,ydim)) < 4){
          nclust <- 1
          maxsz <- 0
          mxsz <- xdim
          mysz <- ydim
          bmap <- matrix((-1), nrow = length(bonds[,1]) + 1, 1)

        } else {

          #  %   THE SUBNETWORK CODE BEGINS
          lista <- matrix(0, nrow = numbonds, ncol = 1)
          listb <- matrix(0, nrow = numbonds, ncol = 1)
          nclust <- 0
          maxsz <- 0
          thru <- xdim

          while (thru != 0) {
            #%  the following code extracts connected
            #%   sub-networks of the non-trivial
            #%   bonds.  NB: lista/b can have redundant entries due to
            #%   multiple-connected subnetworks

            w <- which(bonds[, 2] >= 0)

            lista[1] = bonds[w[1],2]
            listb[1] = bonds[w[1],1]
            bonds[w[1], ] <- (-1) * (nclust+1)
            # bonds;
            adda <- 1
            addb <- 1
            donea <- 0
            doneb <- 0
            if ((donea != adda) || (doneb != addb)){
              true <- FALSE
            } else {
              true <- TRUE
            }

            while (!true){

              if (donea != adda) {
                w <- which(bonds[,2] == lista[donea+1])
                ngood <- length(w)
                if (ngood != 0) {
                  listb[(addb+1):(addb+ngood),1] <- bonds[w,1]
                  bonds[w,] <- (-1)*(nclust+1)
                  addb <- addb+ngood;
                }
                donea <- donea + 1
              }
              if (doneb != addb){
                w <- which(bonds[,1] == listb[doneb+1])
                ngood <- length(w);
                if (ngood != 0) {
                  lista[(adda+1):(adda+ngood),1] <- bonds[w,2]
                  bonds[w,] <- (-1)*(nclust+1)
                  adda <- adda+ngood;
                }
                doneb <- doneb + 1
              }
              if ((donea != adda) || (doneb != addb)){
                true <- FALSE
              } else {
                true = TRUE
              }
            }

            pp <- sort(listb[1:doneb])
            pqx <- order(listb[1:doneb])
            #%unx =  unq(listb(1:doneb),pqx);
            #%implanting unq directly
            arr <- listb[1:doneb]
            q <- arr[pqx]
            indices <- which(q != circshift(q,-1))
            count <- length(indices)
            if (count > 0){
              unx <- pqx[indices]
            } else {
              unx <- length(q) -1
            }

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            xsz <- length(unx)

            pp <- sort(lista[1:donea])
            pqy <- order(lista[1:donea])

            #%uny =  unq(lista(1:donea),pqy);
            #%implanting unq directly
            arr <- lista[1:donea]
            q <- arr[pqy]
            indices <- which(q != circshift(q,-1))
            count <- length(indices)
            if (count > 0){
              uny <- pqy[indices]
            } else {
              uny <- length(q) -1
            }

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ysz <- length(uny)
            if ((xsz*ysz) > maxsz){
              maxsz <- xsz*ysz
              mxsz <- xsz
              mysz <- ysz
            }

            thru <- thru - xsz
            nclust <- nclust + 1
          }
          bmap <- bonds[,2]
        }

        #% THE SUBNETWORK CODE ENDS
        #% put verbose in for Jaci

        ## Adjusting nclust
        all_clusts <- unique(abs(bmap))
        nclust <- length(all_clusts)

        #%   THE PERMUTATION CODE BEGINS
        for (nc in 1:nclust){

          #message(paste0("nc=", nc))
          w <- which(bmap == (-1)*(nc))

          nbonds <- length(w)
          bonds <- mbonds[w,]
          lensq <- bondlen[w]

          pq <- sort(bonds[,1])
          st <- order(bonds[,1])
          #%un = unq(bonds(:,1),st);
          #%implanting unq directly
          arr <- bonds[,1]
          q <- arr[st]
          indices <- which(q != circshift(q,-1))
          count <- length(indices)
          if (count > 0) {
            un <- st[indices]
          } else {
            un <- length(q) - 1
          }

          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          uold <- bonds[un,1]
          nold <- length(uold)

          #%un = unq(bonds(:,2));
          #%implanting unq directly
          indices <- which(bonds[, 2] != circshift(bonds[, 2], -1))
          count <- length(indices)
          if (count > 0){
            un <- indices
          } else {
            un <- length(bonds[,2]) -1
          }

          # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          unew <- bonds[un,2]
          nnew <- length(unew)

          if (nnew > 5){
            rnsteps <- 1
            for (ii in 1:nnew){
              rnsteps <- rnsteps * length(which(bonds[,2] == unew[ii]))
              if (rnsteps >= 50000 && rnsteps < 200000){

                warn_log[[length(warn_log) + 1]]  <- i
                #message('Warning: difficult combinatorics encountered')
              } else if (rnsteps >= 200000 && !force_exec){

                #message(paste0("i=", i, "... nc=", nc))
                #message('Excessive Combinitorics LOOK WHAT YOU HAVE DONE TO ME!!!')
                #close(uniquevar);

                warn_message(warn_log = warn_log, quiet = quiet)
                message(paste0("Excessive Combinitorics encountered while processing slide ", i,
                               ". Quitting now... Try using a smaller maxdisp."))

                return(NULL)

              } else if (rnsteps < 5000000 && force_exec) {

                warn_log[[length(warn_log) + 1]]  <- i

              } else if (rnsteps >= 200000) {

                warn_message(warn_log = warn_log, quiet = quiet)
                message(paste0("Excessive Combinitorics encountered while processing slide ", i,
                               ". Quitting now... Try using a smaller maxdisp."))

                return(NULL)
              }
            }
          }

          st <- rep(0, times = nnew)
          fi <- rep(0, times = nnew)

          h <- rep(0, times = nbonds)
          ok <- rep(1, times = nold)
          nlost <- (nnew - nold) > 0


          for (ii in 1:nold) {
            h[which(bonds[,1] == uold[ii])] <- ii
          }
          st[1] <- 1
          fi[nnew] <- nbonds; ##----------------------% check this later
          if (nnew > 1){
            sb <- bonds[, 2]
            sbr <- circshift(sb,1)
            sbl <- circshift(sb,-1)
            st[2:length(st)] <- which(sb[2:length(sb)] != sbr[2:length(sbr)]) + 1
            fi[1:(nnew-1)] <- which(sb[1:(nbonds-1)] != sbl[1:(nbonds-1)])
          }
          #%                if i-1 == 13
          #%                    hi
          #%                end


          checkflag <- 0
          while (checkflag != 2){

            pt <- st - 1
            lost <- matrix(0, nrow = nnew, ncol = 1)
            who <- 0
            losttot <- 0
            mndisq <- nnew*maxdisq


            while (who != (-1)){

              if (pt[(who+1)] != fi[(who+1)]){


                w <- which(ok[h[(pt[(who+1)]+1):(fi[(who+1)])]]!=0) ###---------------% check this -1
                ngood <- length(w)
                if (ngood > 0){
                  if (pt[(who+1)] != st[(who+1)]-1) {
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  pt[(who+1)] <- pt[(who+1)] + w[1]
                  ok[h[pt[(who+1)]]] <- 0
                  if (who == (nnew - 1)){
                    ww <- which(lost == 0)
                    dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)

                    if (dsq < mndisq){
                      minbonds <- pt[ww]
                      mndisq <- dsq
                    }
                  } else {
                    who <- who+1
                  }
                } else {
                  if (!lost[(who+1)] & (losttot != nlost)){
                    lost[(who+1)] <- 1
                    losttot <- losttot + 1
                    if (pt[(who+1)] != st[(who+1)] - 1){
                      ok[h[pt[(who+1)]]] <- 1
                    }
                    if (who == (nnew-1)){
                      ww <- which(lost == 0)
                      dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)
                      if (dsq < mndisq){
                        minbonds <- pt[ww]
                        mndisq <- dsq
                      }
                    } else {
                      who <- who + 1
                    }

                  } else {
                    if (pt[(who+1)] != (st[(who+1)] - 1)){
                      ok[h[pt[(who+1)]]] <- 1
                    }
                    pt[(who+1)] <- st[(who+1)] - 1
                    if (lost[(who+1)]){
                      lost[(who+1)] <- 0
                      losttot <- losttot -1
                    }
                    who <- who - 1
                  }
                }

              } else {
                if (!lost[(who+1)] && (losttot != nlost)){
                  lost[(who+1)] <- 1
                  losttot <- losttot + 1
                  if (pt[(who+1)] != st[(who+1)]-1){
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  if (who == (nnew - 1)) {
                    ww <- which(lost == 0)
                    dsq <- sum(lensq[pt[ww]]) + (losttot * maxdisq)

                    if (dsq < mndisq){
                      minbonds <- pt[ww]
                      mndisq <- dsq
                    }
                  } else {
                    who <- who + 1
                  }
                } else {
                  if (pt[(who+1)] != st[(who+1)] - 1){
                    ok[h[pt[(who+1)]]] <- 1
                  }
                  pt[(who+1)] <- st[(who+1)] - 1
                  if (lost[(who+1)]){
                    lost[(who+1)] <- 0
                    losttot <- losttot - 1
                  }
                  who <- who -1
                }
              }
            }

            checkflag <- checkflag + 1
            if (checkflag == 1){
              plost <- min(c(matfix(mndisq/maxdisq) , (nnew -1)))
              if (plost > nlost){
                nlost <- plost
              } else {
                checkflag <- 2
              }
            }
          }
          #%   update resx using the minimum bond configuration

          resx[ispan, labely[bonds[minbonds, 2]]] <- eyes[labelx[(bonds[minbonds,1] + 1)]]
          found[labelx[(bonds[minbonds,1] + 1)], 1] <- 1

        }

        #%   THE PERMUTATION CODE ENDS
      }

      w <- which(resx[ispan,] >= 0)
      nww <- length(w)

      if (nww > 0){
        pos[w,] <- xyzs[resx[ispan,w], (1:dim)]
        if (goodenough > 0){
          nvalid[w] <- nvalid[w] + 1
        }
      }  #-----------------------  %go back and add goodenough keyword thing
      newguys <- which(found == 0)
      nnew <- length(newguys)

      if (nnew > 0) {             ##% & another keyword to workout inipos
        newarr <- matrix(-1, nrow = zspan, ncol = nnew)

        # cbind?
        resx <- cbind(resx, newarr)

        resx[ispan, ((n+1):ncol(resx))] <- eyes[newguys]
        pos <- rbind(pos, xyzs[eyes[newguys],(1:dim)])
        nmem <- matrix(0, nrow = nnew, ncol = 1)
        mem <- c(mem, nmem)
        nun <- 1:nnew
        uniqid <- c(uniqid, ((nun) + maxid))
        maxid <- maxid + nnew
        if (goodenough > 0){
          dumphash <- c(dumphash, t(matrix(0, nrow = 1, ncol = nnew)))
          nvalid <- c(nvalid, t(matrix(1, nrow = 1, ncol = nnew)))
        }

        #% put in goodenough
        n <- n + nnew

      }

    } else {
      #' Warning- No positions found for t='
      message("@@", appendLF = FALSE)
    }

    w <- which(resx[ispan,] != (-1))
    nok <- length(w)
    if (nok != 0){
      mem[w] <- 0
    }

    #---------------------------------------------------
    mem <- mem + (0 + (cbind(resx[ispan,]) == -1))
    wlost <- which(mem == memory_b+1)
    nlost <- length(wlost)

    if (nlost > 0){
      pos[wlost, ] <- (-maxdisp)
      if (goodenough > 0){
        wdump <- which(nvalid[wlost] < goodenough)
        ndump <- length(wdump);
        if (ndump > 0){
          dumphash[wlost[wdump]] <- 1
        }
      }
      #% put in goodenough keyword stuff if
    }

    if ((ispan == zspan) | (i == z)){
      nold <- length(bigresx[1,])
      nnew <- n - nold;
      if (nnew > 0){

        newarr <- matrix(-1, nrow = z, ncol = nnew)

        ## bigresx <- c(bigresx, newarr)
        bigresx <- cbind(bigresx, newarr)
      }
      if (goodenough > 0){
        if ((sum(dumphash)) > 0){
          wkeep <- which(dumphash == 0)
          nkeep <- length(wkeep)
          resx <- resx[ ,wkeep]
          bigresx <- bigresx[, wkeep]
          pos <- pos[wkeep, ]
          mem <- mem[wkeep]
          uniqid <- uniqid[wkeep]
          nvalid <- nvalid[wkeep]
          n <- nkeep
          dumphash <- matrix(0, nrow = nkeep, ncol = 1)
        }
      }

      #% again goodenough keyword
      if (!quiet) {

        message(paste0(i, ' of ' , z, ' done. Tracking ', ntrk, ' particles. ', n, ' tracks total.'))

      }

      if (!is.matrix(bigresx) || nrow(resx) > nrow(bigresx)) {
        bigresx <- rbind(bigresx)
        bigresx <- rbind(bigresx,
                         matrix(-1, nrow = (nrow(resx) - nrow(bigresx)),
                                ncol = ncol(bigresx)))
      }

      bigresx[(i-(ispan)+1):i,]  <- resx[1:ispan,]
      resx <- matrix((-1), nrow = zspan, ncol = n)

      wpull <- which(pos[ ,1] == (-1 * maxdisp))
      npull <- length(wpull)

      if (npull > 0){
        lillist <- list()
        for (ipull in 1:npull){
          wpull2 <- which(bigresx[, wpull[ipull]] != (-1))
          npull2 <- length(wpull2)


          thing = cbind(bigresx[wpull2,wpull[ipull]],
                        rep(x = uniqid[wpull[ipull]], times = npull2))

          #thing <- c(bigresx[wpull2, wpull[ipull]], ),zeros(npull2,1)+uniqid(wpull(ipull))];
          #lillist = [lillist;thing];
          lillist[[length(lillist) + 1]] <- thing

        }
        olist[[length(olist) + 1]] <- do.call(rbind, lillist)

      }



      wkeep <- which(pos[, 1] >= 0)
      nkeep <- length(wkeep)
      if (nkeep == 0) {
        message ('Were going to crash now, no particles....')
      }
      resx <- resx[,wkeep]
      bigresx <- bigresx[, wkeep]
      pos <- pos[wkeep, ]
      mem <- mem[wkeep]
      uniqid <- uniqid[wkeep]
      n <- nkeep;
      dumphash <- matrix(0, nrow = nkeep, ncol =1)
      if (goodenough > 0){
        nvalid <- nvalid[wkeep]
      }
    }
    #waitbar(i / z)
  }

  if (goodenough > 0){
    nvalid <- apply(bigresx >= 0 , 2, sum)
    wkeep <- which(nvalid >= goodenough)
    nkeep <- length(wkeep)
    if (nkeep == 0){
      for (i in 1:10){
        message('You are not going any further, check your params and data')
      }
      message('the code broke at line 1995')
      return()
    }
    if (nkeep < n){
      bigresx <- bigresx[, wkeep]
      n <- nkeep
      uniqid <- uniqid[wkeep]
      pos <- pos[wkeep, ]
    }
  }

  wpull <- which(pos[, 1] != ((-2) * maxdisp))
  npull <- length(wpull);
  if (npull > 0) {
    lillist <- list()
    for (ipull in 1:npull){
      wpull2 <- which(bigresx[, wpull[ipull]] != (-1))
      npull2 <- length(wpull2)
      thing <- cbind(bigresx[wpull2, wpull[ipull]],
                     rep(uniqid[wpull[ipull]], times = npull2))
      lillist[[length(lillist) + 1]] <- thing
    }

    olist[[length(olist) + 1]] <- do.call(rbind, lillist)
  }

  olist <- do.call(rbind, olist)
  #%bigresx = 0;
  #%resx = 0;

  nolist <- nrow(olist)
  res <- matrix(0, nrow = nolist, ncol = (dd+1))
  for (j in 1:dd){
    res[, j] <- xyzs[olist[, 1], j]
  }
  res[, dd+1] <- olist[,2]

  #% this is uberize included for simplicity of a single monolithic code

  ndat <- ncol(res)
  newtracks <- res


  #%u=unq(newtracks(:,ndat));

  #% inserting unq
  indices <- which(newtracks[, ndat] != circshift(newtracks[, ndat], -1))
  count <- length(indices)
  if (count > 0){
    u <- indices
  } else {
    u <- nrow(newtracks) - 1
  }


  ntracks <- length(u)
  u <- c(0, u)
  for (i in 2: (ntracks + 1)){
    newtracks[(u[(i-1)]+1):u[i], ndat] = (i - 1)
  }

  #% end of uberize code
  warn_message(warn_log = warn_log, quiet = quiet)
  return(newtracks)
}


#' Compute Cell Migration Statistics
#'
#' Calculate the statistics from X/Y positional data obtained from cell tracks
#'
#' @param tracks data.frame with cell tracks information
#' @param interval_time integer, time interval between two successive frames were taken
#' @param pixel_micron integer, image resolution, i.e. number of pixels per micron
#'
#' @return list of stats calculated for the cell tracks. Info include variables of speed,
#' distance, euclidean displacement, persistence, angular displacement,
#' yFMI, xFMI, y-displacement, x-displacement and frames
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cbind(
#'   c(30, 35, 1, 5, 6, 7, 50, 55, 56, 58),
#'   c(29, 37, 2, 7, 4, 9, 40, 50, 59, 49),
#'   c( 1,  2,  1, 2, 3, 4, 1, 2, 3, 4),
#'   c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
#' cellmigRation:::MigrationStats(x0, 10, 10)
#'
#'
#' @keywords internal
MigrationStats <- function(tracks, interval_time, pixel_micron) {

  #
  speed <- list()
  distance <- list()
  frames <- list()
  euclid <- list()
  persistence <- list()
  initial <- list()
  final <- list()
  deltaX = list()
  deltaY = list()

  # Keep track
  kept <- list()

  cell_number <- sort(unique(tracks[,4]))

  for (i in cell_number){
    data <- tracks[tracks[,4] == i, ]

    if (!"matrix" %in% class(data) ) {
      next
    } else (
      kept[[length(kept) + 1]] <- i
    )

    # obtain X-axis and Y-axis positional data for the specified track
    X <- data[,1]
    Y <- data[,2]

    x1 <- X[1]
    xEnd <- X[length(X)]
    y1 <- Y[1]
    yEnd <- Y[length(Y)]

    initial[[length(initial) + 1]] <- c(x=x1, y=y1)
    final[[length(final) + 1]] <- c(x=xEnd, y=yEnd)
    delX <- as.numeric(xEnd - x1)
    delY <- as.numeric(yEnd - y1)
    deltaX[[length(deltaX) + 1]] <- delX
    deltaY[[length(deltaY) + 1]] <- delY
    # calculate euclidean distance (vector displacement of the cell)
    E <- as.numeric(sqrt((delX)^2 + (delY)^2))
    euclid[[length(euclid) + 1]] <- E

    # add subsequent displacements of the cell
    cumulative_displacements <- as.numeric(cumsum(sqrt(diff(X)^2 + diff(Y)^2)) )

    # sum of the displacements between each cell centroid for the given
    # track

    distance[[length(distance) + 1]] <- max(cumulative_displacements)

    # calculate cell persistence
    persistence[[length(persistence) + 1]] = E/max(cumulative_displacements)

    # total number of frames that cell centroid was tracked ( can be
    # greater than number of frames where centroid was identified given the
    # param.mem parameter
    # total number of time intervals through which the cell has been tracked
    totalframes = data[nrow(data), 3] - data[1, 3]

    # sum of all individual displacemnts divided by the time that cell
    # centroid was tracked
    ds_dt <- max(cumulative_displacements)/(totalframes*interval_time)
    speed[[length(speed) + 1]] <- ds_dt

    frames[[length(frames) + 1]] <- totalframes
  }

  # Expand resulting lists
  speed <- do.call(c, speed)
  distance <- do.call(c, distance)
  frames <- do.call(c, frames)
  euclid <- do.call(c, euclid)
  persistence <- do.call(c, persistence)
  initial <- do.call(rbind, initial)
  final <- do.call(rbind, final)
  kept <- do.call(c, kept)
  deltaX <- do.call(c, deltaX)
  deltaY <- do.call(c, deltaY)


  # calculate angular displacement of cells trajectory
  arccos <- deltaX/euclid
  theta <- acos(arccos)

  for (j in 1:length(arccos)){
    if  (arccos[j] < 0 & deltaY[j] > 0) {
      theta[j] <- 2*pi - theta[j]
    }

    if (arccos[j] > 0 & deltaY[j] > 0) {
      theta[j] <- 2*pi - theta[j]
    }
  }

  #theta = theta.*((2*pi)/360);
  deltaY <-  deltaY * (-1)
  yfmi <- deltaY / distance
  xfmi <- deltaX / distance
  deltaY <- deltaY * pixel_micron
  deltaX <- deltaX * pixel_micron
  speed <- speed * pixel_micron;
  distance <- distance * pixel_micron
  euclid <- euclid * pixel_micron

  OUT <- list(
    speed = speed,
    distance = distance,
    frames = frames,
    euclid = euclid,
    persistence = persistence,
    initial = initial,
    final = final,
    yfmi = yfmi,
    xfmi = xfmi,
    deltaX = deltaX,
    deltaY = deltaY,
    kept = kept,
    theta = theta
  )
  return(OUT)
}






#' Compute Tracks Stats
#'
#' Wrapper for the MigrationStats() function. It computes statistics for a
#' trackedCells object where cells have already been tracked.
#'
#' @param tc_obj a \code{trackedCells} object
#' @param time_between_frames integer, time interval between two successive frames were taken
#' @param resolution_pixel_per_micron integer, image resolution, i.e. number of pixels per micron
#'
#' @return a trackedCells object, including cell track statistics
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- CellTracker(tc_obj = x, lnoise = 6, diameter = 16, threshold = 10)
#' x <- ComputeTracksStats(x, time_between_frames = 10, resolution_pixel_per_micron = 20)
#' getCellsStats(x)[1:10,]
#' }
#'
#' @export
ComputeTracksStats <- function(tc_obj, time_between_frames, resolution_pixel_per_micron)
{

  if(tc_obj@ops$track == 0)
    stop("You need to run CellTracker() before computing stats")

  # RETRIEVE
  my_tracks <- tc_obj@tracks

  # DO
  handles <- MigrationStats(tracks = my_tracks, interval_time = time_between_frames,
                            pixel_micron = resolution_pixel_per_micron);


  sz <- length(handles$speed)
  handles$cell_number <- 1:sz

  cell_stats <- data.frame(Cell_Number = handles$cell_number,
                           Speed = handles$speed,
                           Distance = handles$distance,
                           Displacement = handles$euclid,
                           Persistence = handles$persistence,
                           Degrees = handles$theta,
                           YFMI = handles$yfmi,
                           XFMI = handles$xfmi,
                           Y_displacement = handles$deltaY,
                           X_displacement = handles$deltaX,
                           Frames = handles$frames,
                           stringsAsFactors = FALSE)

  # to organize population stats
  my_colz <- c("Speed", "Distance", "Displacement", "Persistence", "YFMI",
               "XFMI", "Y_displacement", "X_displacement")

  my_rows <- lapply(my_colz, function(cl) {
    tmp <- cell_stats[, cl]
    data.frame(mean = mean(tmp, na.rm = TRUE),
               SD = sd(tmp, na.rm = TRUE),
               median = median(tmp, na.rm = TRUE),
               min = min(tmp, na.rm = TRUE),
               max = max(tmp, na.rm = TRUE))
  })
  my_rows <- do.call(rbind, my_rows)
  rownames(my_rows) <- my_colz

  # compute sum of cos and sin of angles
  r <- sum(exp(1i*handles$theta))

  # obtain mean angle
  meanTheta <- Arg(r)
  degrees <- meanTheta/pi*180
  my_rows <- rbind(my_rows,
                   Angle = data.frame(mean = sz, SD = degrees, median = NA, min = NA, max = NA))
  rownames(my_rows)[c(2,3)] <- c("Total_displacement", "Euclidean_displacement")
  pop_stats <- my_rows

  # Attach, return
  tc_obj@ops$stats <- 1
  tc_obj@stats <- list(population = pop_stats,
                       cells = cell_stats)

  return(tc_obj)
}


#' Optimize Detection Params
#'
#' Optimize Detection Parameters for running a cell tracking job
#'
#' @details The lnoise param is used to guide a lowpass blurring operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param is used for a background correction following
#' the initial image convolution
#' \itemize{
#'
#'   \item \strong{lnoise}: Characteristic lengthscale of noise in pixels.
#'   Additive noise averaged over this length should vanish. May assume any positive floating value.
#'   May be also set to 0, in which case only the highpass "background subtraction" operation is performed.
#'
#'   \item \strong{lobject} Integer length in pixels somewhat larger than a typical object.
#'   Can also be set to 0, in which case only the lowpass "blurring" operation defined by lnoise is done
#'   without the background subtraction defined by lobject
#'
#'   \item \strong{threshold} Numeric. By default, after the convolution, any negative pixels are reset
#'   to 0.  Threshold changes the threshhold for setting pixels to 0.  Positive values may be useful
#'   for removing stray noise or small particles.
#'
#' }
#'
#'
#'
#' @param tc_obj a \code{trackedCells} object
#' @param lnoise_range numeric vector of lnoise values to be used in the optimization step. Can be NULL
#' @param min.px.diam integer, minimum diameter of a particle (cell).
#' Particles with a diameter smaller than min.px.diam are discarded
#' @param diameter_range numeric vector of diameter values to be used in the optimization step. Can be NULL
#' @param threshold_range numeric vector of threshold values to be used in the optimization step. Can be NULL
#' @param target_cell_num integer, the expected (optimal) number of cells to be detected in each frame
#' @param threads integer, number of cores to use for parallelization
#' @param quantile.val numeric, argument passed to EstimateDiameterRange(). If NULL, it is defaulted to 0.99
#' @param px.margin numeric, argument passed to EstimateDiameterRange(). If NULL, it ia defaulted to 2
#' @param plot if `TRUE`, plots results in the end
#' @param verbose shall information about the progress of the operation be printed to screen/console
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(tc_obj = x)
#' getOptimizedParams(x)
#' }
#'
#'
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom graphics par
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @import foreach
#'
#'
#' @export
OptimizeParams <- function(tc_obj, lnoise_range = NULL, min.px.diam = 5,
                           diameter_range = NULL, threshold_range = NULL,
                           target_cell_num = NULL, threads = 1,
                           quantile.val = NULL, px.margin= NULL,
                           plot=FALSE, verbose = FALSE)

{
  # do
  stack_img <- tc_obj@images

  # Nested f(x)
  all_combos <- function(...){
    xx <- list(...)
    zz <- names(xx)

    # Init
    out <- data.frame(xx[[1]], stringsAsFactors = FALSE)
    colnames(out) <- zz[1]

    # Keep attaching
    for (j in 2:length(xx)) {
      TMP <- xx[[j]]
      nuOUT <- list()
      for (z in TMP){
        tmpout <- cbind(out, data.frame(TMP = z, stringsAsFactors = FALSE))
        nuOUT[[length(nuOUT) + 1]] <- tmpout
      }
      out <- do.call(rbind, nuOUT)
      colnames(out)[j] <- zz[j]
    }
    return(out)
  }

  ## ----- debugging -----
  #bpass = cellmigRation:::bpass
  #pkfnd = cellmigRation:::pkfnd
  #VisualizeImg = cellmigRation:::VisualizeImg
  #cntrd = cellmigRation:::cntrd
  #NextOdd = cellmigRation:::NextOdd
  #VisualizeCntr = cellmigRation:::VisualizeCntr
  #track = cellmigRation:::track
  ## ----- endo of debugging -----

  if (!verbose) {
    tryCatch(sink(file = "/dev/null", type = "message"), error = function(e) {NULL})
    tryCatch(sink(file = "/dev/null", type = "output"), error = function(e) {NULL})

    on.exit(expr = {
      tryCatch(sink(file = NULL, type = "message"), error = function(e) {NULL});
      tryCatch(sink(file = NULL, type = "output"), error = function(e) {NULL})})
  }


  # select mid signal image
  imgSums <- sapply(stack_img$images, sum, na.rm = TRUE)
  med.i <- ifelse(length(imgSums) %% 2 == 0, length(imgSums) / 2, (0.5 * length(imgSums) + 0.5))
  r.i <- order(imgSums)[med.i]

  tmp_img <- stack_img$images[[r.i]]

  # Define param ranges
  if (is.null(px.margin)) {
    px.margin <- 2
  }
  if (is.null(quantile.val)) {
    quantile.val <- 0.99
  }

  estRDI <- tryCatch({
    EstimateDiameterRange(x = tmp_img, px.margin = px.margin,
                          min.px.diam = min.px.diam,
                          quantile.val = quantile.val, plot = FALSE)},
    error = function(e) NULL)

  # diam range
  if(is.null(diameter_range) && !is.null(estRDI)) {
    diameter_range <- c(floor(estRDI$q75.diam - 1), ceiling(1.25 * as.numeric(estRDI$q95.diam)))
    diameter_range[diameter_range < min.px.diam] <- min.px.diam
    diameter_range <- unique(as.integer(diameter_range))
    diameter_range <- unique(as.integer(
      seq(min(diameter_range), max(diameter_range), length.out = 3)))

  } else if (is.null(diameter_range)) {
    diameter_range <- c(10, 30, 90)
  }

  # num cell
  if(is.null(target_cell_num) && !is.null(estRDI)) {
    target_cell_num <- estRDI$estim.cell.num
  } else if (is.null(target_cell_num)) {
    target_cell_num <- 100
  }

  # Define param ranges
  if(is.null(lnoise_range))
    lnoise_range <- unique(as.integer(seq(from = min.px.diam,
                                          to = quantile(estRDI$raw$LEN, probs = 0.25),
                                          length.out = 3)))

  if(is.null(threshold_range)) {
    threshold_range <- seq(max(0, (min(tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)], na.rm = TRUE) - 1)),
                           (1 + quantile(tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)], probs = 0.75)),
                           length.out = 4)
    threshold_range <- unique(as.integer(threshold_range))
  }

  # Al params
  all_params <- all_combos(image.i = med.i,
                           lnoise = lnoise_range,
                           diameter = diameter_range,
                           threshold = threshold_range)

  all_params <- all_params[all_params$diameter > (4 + all_params$lnoise), ]
  rownames(all_params) <- NULL

  if(nrow(all_params) < 4) {
    message("There is a problem with the param ranges that were submitted")
    message("Please, try again with different param ranges")
    return(tc_obj)
  }

  if(nrow(all_params) > 20) {

    all_params$TMPdiff <- all_params$diameter - all_params$lnoise
    all_params <- all_params[order(all_params$TMPdiff, decreasing = TRUE),]
    rownames(all_params) <- NULL
    all_params <- all_params[1:20,]
  }


  # Verbose
  if (verbose) {

    message(paste0("Testing ", nrow(all_params), " combination(s) of params."), appendLF = TRUE)
    message("This may take some time.", appendLF = TRUE)

    message("Processing ", appendLF = FALSE)
  }

  ##
  ## Parallelize please
  j <- NULL

  # how many cores can we use?
  num_parallelCores <- threads
  debugging <- TRUE

  max.cores <- parallel::detectCores()
  max.cores <- max.cores - 1
  max.cores <- ifelse(max.cores < 1, 1, max.cores)
  my.test <- 1 <= num_parallelCores & num_parallelCores <= max.cores
  use.cores <- ifelse(my.test, num_parallelCores, max.cores)

  # Adjust if NA
  if (is.na(use.cores)) {
    use.cores <- 1
  }

  # cores = 1, do not parallelize
  if (use.cores == 1) {

    # Initialize collector (list)
    all_results <- list()

    for (i in 1:nrow(all_params)){

      # Verbose
      if (verbose)
        message(".", appendLF = FALSE)

      #VisualizeImg(tmp_img)
      b <- bpass(image_array = tmp_img,
                 lnoise = all_params$lnoise[i],
                 lobject = all_params$diameter[i],
                 threshold = all_params$threshold[i])
      tmpOUT <- list(img = b)

      tryCatch({
        pk <- suppressMessages(
          pkfnd(im = b,
                th = all_params$threshold[i],
                sz = NextOdd(all_params$diameter[i])))

        cnt <- suppressMessages(
          cntrd(im = b, mx = pk,
                sz = NextOdd(all_params$diameter[i])))

        tmpOUT[["count"]] <- nrow(cnt)

      }, error = function(e) {
        tmpOUT[["count"]] <- 0

      })
      all_results[[i]] <- tmpOUT
    }

    # cores > 1, DO parallelize!
  } else {

    if (verbose) {
      cl <- suppressMessages(parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }

    suppressMessages(doParallel::registerDoParallel(cl))

    # Nothing to export! "tmp_img", "all_params" automatically exported
    #stuffToExp <- c("tmp_img", "all_params")
    stuffToExp <- c()
    suppressMessages(parallel::clusterExport(cl, stuffToExp))

    ## %dopar%
    all_results <-
      tryCatch(foreach::foreach(j = (1:nrow(all_params)),
                                .verbose = verbose,
                                .packages = "cellmigRation") %dopar% {

                                  # Verbose
                                  message(".", appendLF = FALSE)

                                  #VisualizeImg(tmp_img)
                                  b <- bpass(image_array = tmp_img,
                                             lnoise = all_params$lnoise[j],
                                             lobject = all_params$diameter[j],
                                             threshold = all_params$threshold[j])
                                  tmpOUT <- list(img = b)

                                  tryCatch({
                                    pk <- suppressMessages(
                                      pkfnd(im = b,
                                            th = all_params$threshold[j],
                                            sz = NextOdd(all_params$diameter[j])))

                                    cnt <- suppressMessages(
                                      cntrd(im = b, mx = pk,
                                            sz = NextOdd(all_params$diameter[j])))

                                    tmpOUT[["count"]] <- nrow(cnt)

                                  }, error = function(e) {
                                    tmpOUT[["count"]] <- 0

                                  })
                                  tmpOUT
                                }, error = (function(e) {
                                  print(e)
                                  try(parallel::stopCluster(cl), silent = TRUE)
                                  return(NULL)
                                }))
    message("Done!", appendLF = TRUE)
    try({suppressWarnings(parallel::stopCluster(cl))}, silent = TRUE)
  }

  # Attach counts
  all_params$counts <- do.call(c, lapply(all_results, function(x) {
    tmp <- x$count
    ifelse(is.null(tmp), 0, tmp)}))
  all_params$i <- 1:nrow(all_params)

  # Return
  all_params$diff100 <- abs(target_cell_num - all_params$counts)
  ord_params <- all_params[order(all_params$diff100), ]
  ret.i <- head(ord_params$i, n = 9)
  best_params <- list()

  top.i <- 1
  for (ri in ret.i) {

    if (top.i == 1) {
      best_params[["lnoise"]] <- ord_params$lnoise[ord_params$i == ri]
      best_params[["diameter"]] <- ord_params$diameter[ord_params$i == ri]
      best_params[["threshold"]] <- ord_params$threshold[ord_params$i == ri]
    }

    myLAB <- paste0("Pick #", top.i, "; Cell_count=", ord_params$counts[ord_params$i == ri], "\n")
    myLAB <- paste0(myLAB, "lnoise=", ord_params$lnoise[ord_params$i == ri], "; ",
                    "diameter=", ord_params$diameter[ord_params$i == ri], "; ",
                    "threshold=", ord_params$threshold[ord_params$i == ri])

    if (plot) {
      curPAR <- par(no.readonly = TRUE)
      par(mfrow = c(3, 3))
      on.exit(expr = {par(curPAR)})
      VisualizeImg(img_mtx = all_results[[ri]]$img, main = myLAB)
    }

    top.i <- top.i + 1
  }

  # Extract_all_img
  allIMG <- lapply(all_results, function(x) {x$img})

  #return(list(auto_params = best_params,
  #            results = all_params,
  #            images = allIMG))


  tc_obj@ops$optimized_params <- 1
  tc_obj@optimized <- list(auto_params = best_params,
                           results = all_params)

  return(tc_obj)
}


#' Compute Cell Tracks
#'
#' Analyze Stacks, detect cells in each frame, and analyze cell tracks over time
#'
#' @details The lnoise param is used to guide a lowpass blurring operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param is used for a background correction following
#' the initial image convolution
#' \itemize{
#'
#'   \item \strong{lnoise}: Characteristic lengthscale of noise in pixels.
#'   Additive noise averaged over this length should vanish. May assume any positive floating value.
#'   May be also set to 0, in which case only the highpass "background subtraction" operation is performed.
#'
#'   \item \strong{lobject} Integer length in pixels somewhat larger than a typical object.
#'   Can also be set to 0, in which case only the lowpass "blurring" operation defined by lnoise is done
#'   without the background subtraction defined by lobject
#'
#'   \item \strong{threshold} Numeric. By default, after the convolution, any negative pixels are reset
#'   to 0.  Threshold changes the threshhold for setting pixels to 0.  Positive values may be useful
#'   for removing stray noise or small particles.
#'
#' }
#'
#'
#' @param tc_obj a \code{trackedCells} object.
#' @param import_optiParam_from a \code{trackedCells} object (optional) used to
#' import optimized parameters; can be NULL.
#' @param min_frames_per_cell numeric, minimum number of consecutive frames in which
#' a cell shall be found in order to retain that cell in the final cell tracks data.frame. Defaults to 1.
#' @param lnoise numeric, lnoise parameter; can be NULL if OptimizeParams() has already been run
#' @param diameter numeric, diameter parameter; can be NULL if OptimizeParams() has already been run
#' @param threshold numeric, threshold parameter; can be NULL if OptimizeParams() has already been run
#' @param maxDisp numeric,  maximum displacement of a cell per time interval.
#' When many cells are detected in each frame, small maxDisp values should be used.
#' @param memory_b numeric, memory_b parameter as used in the original track.m function.
#' In the current R implementation, only the value memory_b=0 is accepted
#' @param goodenough numeric, goodenough parameter as used in the original track.m function.
#' In the current R implementation, only the value goodenough=0 is accepted
#' @param threads integer, number of cores to use for parallelization
#' @param show_plots logical, shall cells detected in each frame of the image stack be visualized
#' @param verbose logical, shall info about the progress of the cell tracking job be printed
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(x)
#' x <- CellTracker(x)
#' getTracks(x)[1:10,]
#' }
#'
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#'
#' @export
CellTracker <- function(tc_obj, import_optiParam_from = NULL,
                        min_frames_per_cell = 1,
                        lnoise = NULL, diameter = NULL,
                        threshold = NULL, maxDisp = NULL,
                        memory_b = 0, goodenough = 0,
                        threads = 1, show_plots = FALSE,
                        verbose = FALSE)
{
  # get stuff
  stack_img <- tc_obj@images
  optimal_params <- tc_obj@optimized

  if (length(optimal_params) > 0) {
    my.lnoise <- optimal_params$auto_params$lnoise
    my.diameter <- optimal_params$auto_params$diameter
    my.threshold <- optimal_params$auto_params$threshold
  } else {
    my.lnoise <- NA
    my.diameter <- NA
    my.threshold <- NA
  }

  # if a `import_optiParam_from` is specified, and
  if (!is.null(import_optiParam_from)) {
    chk.tco <- 0
    if ("trackedCells" %in% class(import_optiParam_from)) {
      if(import_optiParam_from@ops$optimized_params == 1) {

        # import
        my.lnoise <- import_optiParam_from@optimized$auto_params$lnoise
        my.diameter <- import_optiParam_from@optimized$auto_params$diameter
        my.threshold <- import_optiParam_from@optimized$auto_params$threshold
        chk.tco <- 1

      }
    }
    if (chk.tco == 0) {
      message("Could not import optimal params from the provided object.")
      message("Is it a 'trackedCells'-class object? Did you run `OptimizeParams()`?")
    }
  }

  custom_params_flag <- 0

  if(!is.null(lnoise) && is.numeric(lnoise)){
    my.lnoise <- lnoise[1]
    custom_params_flag <- 1
  }

  if(!is.null(diameter) && is.numeric(diameter)){
    my.diameter <- diameter[1]
    custom_params_flag <- 1
  }

  if(!is.null(threshold) && is.numeric(threshold)){
    my.threshold <- threshold[1]
    custom_params_flag <- 1
  }

  # Max Disp
  if(!is.null(maxDisp) && is.numeric(maxDisp)) {
    maxDisp <- maxDisp[1]
  } else {
    maxDisp <- NULL
  }

  # Other params
  quiet <- verbose
  force_exec <- FALSE
  j <- NULL

  # At this moment, let's play safe!
  # Impose the following
  if (memory_b != 0)
    message("Currently, only memory_b=0 is supported... Resetting.")
  memory_b <- 0

  if (goodenough != 0)
    message("Currently, only goodenough=0 is supported... Resetting.")
  goodenough <- 0

  # In the end, my params are:
  lnoise <- my.lnoise
  diameter <- my.diameter
  threshold <- my.threshold

  track_params <- list(lnoise = lnoise,
                       diameter = diameter,
                       threshold = threshold,
                       maxDisp = 0,
                       goodenough = goodenough,
                       memory_b = memory_b,
                       force_exec = force_exec,
                       quiet = quiet,
                       verbose = verbose,
                       show_plots = show_plots)

  # Final check
  if (sum(sapply(track_params, is.na)) > 0) {
    message("Make sure to set all params for the analysis, or run OptimizeParams()")
    return(tc_obj)
  }

  if(!is.null(maxDisp))
    track_params$maxDisp <- maxDisp

  ## ----- debugging -----
  #bpass = cellmigRation:::bpass
  #pkfnd = cellmigRation:::pkfnd
  #VisualizeImg = cellmigRation:::VisualizeImg
  #cntrd = cellmigRation:::cntrd
  #NextOdd = cellmigRation:::NextOdd
  #VisualizeCntr = cellmigRation:::VisualizeCntr
  #track = cellmigRation:::track
  ## ----- endo of debugging -----

  # Load stack
  stack <- stack_img

  InfoImage <- stack$attributes[[1]]
  mImage <- stack$dim$width_m
  nImage <- stack$dim$height_n
  NumberImages <- stack$dim$NumberImages
  FinalImage <- stack$images

  ## ----------- Evaluate centroids ---------------------

  # locate centroids, via CentroidArray
  if (verbose)
    message("Computing centroid positions", appendLF = FALSE)

  ##
  ## Parallelize please
  if (!verbose) {
    tryCatch(sink(file = "/dev/null", type = "message"), error = function(e) {NULL})
    tryCatch(sink(file = "/dev/null", type = "output"), error = function(e) {NULL})

    on.exit(expr = {
      tryCatch(sink(file = NULL, type = "message"), error = function(e) {NULL});
      tryCatch(sink(file = NULL, type = "output"), error = function(e) {NULL})})
  }


  # how many cores can we use?
  num_parallelCores <- threads
  debugging <- TRUE

  max.cores <- parallel::detectCores()
  max.cores <- max.cores - 1
  max.cores <- ifelse(max.cores < 1, 1, max.cores)
  my.test <- 1 <= num_parallelCores & num_parallelCores <= max.cores
  use.cores <- ifelse(my.test, num_parallelCores, max.cores)

  # fix for NA cores
  if (is.na(use.cores)) {use.cores <- 1}

  # cores = 1, do not parallelize
  if (use.cores == 1) {

    # Init collectors
    all_centroids <- list()
    all_b <- list()

    for (i in 1:NumberImages) {

      if (verbose)
        message(".", appendLF = FALSE)

      # generate an 1xP array with each column containing centroid output for
      # individual frames
      a <- FinalImage[[i]]
      b <- tryCatch({bpass(image_array = a, lnoise = lnoise,
                           lobject = diameter, threshold = threshold)},
                    error = function(e) {NULL})
      pk <- tryCatch({pkfnd(im = b, th = threshold, sz = NextOdd(diameter))},
                     error = function(e) {NULL})
      cnt <- tryCatch({cntrd(im = b, mx = pk, sz = NextOdd(diameter))},
                      error = function(e) {NULL})
      #b <- bpass(image_array = a, lnoise = lnoise, lobject = diameter, threshold = threshold)
      #pk <- pkfnd(im = b, th = threshold, sz = NextOdd(diameter))
      #cnt <- cntrd(im = b, mx = pk, sz = NextOdd(diameter))

      if (show_plots) {
        VisualizeImg(img_mtx = b, las = 1, main = paste0("Stack num. ", i))
        VisualizeCntr(centroids = cnt, width_px = ncol(b), height_px = nrow(b))
      }

      # determine that frame s has at least 1 valid centroid
      if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
        all_centroids[[length(all_centroids) + 1]] <- cnt
        all_b[[length(all_b) + 1]] <- b

      } else {
        message(paste0('No centroids detectd in frame ', i, ' in the current stack'))
        message('Please, check nuclei validation settings for this image stack.')
      }
    }
    if (verbose)
      message("", appendLF = TRUE)

  } else {

    if (verbose) {
      cl <- suppressMessages(parallel::makeCluster(use.cores, outfile = ""))
    } else {
      cl <- suppressMessages(parallel::makeCluster(use.cores))
    }

    suppressMessages(doParallel::registerDoParallel(cl))
    # Nothing to export! ""FinalImage", "all_params" automatically exported
    #stuffToExp <- c("FinalImage", "all_params")
    stuffToExp <- c()
    suppressMessages(parallel::clusterExport(cl, stuffToExp))

    ## %dopar%
    all_results <-
      tryCatch(foreach::foreach(j = (1:NumberImages),
                                .verbose = verbose,
                                .packages = "cellmigRation") %dopar% {

                                  # Verbose
                                  message(".", appendLF = FALSE)

                                  # generate an 1xP array with each column containing centroid output for
                                  # individual frames
                                  a <- FinalImage[[j]]
                                  #b <- bpass(image_array = a, lnoise = lnoise, lobject = diameter, threshold = threshold)
                                  #pk <- pkfnd(im = b, th = threshold, sz = NextOdd(diameter))
                                  #cnt <- cntrd(im = b, mx = pk, sz = NextOdd(diameter))
                                  b <- tryCatch({bpass(image_array = a, lnoise = lnoise,
                                                       lobject = diameter, threshold = threshold)},
                                                error = function(e) {NULL})
                                  pk <- tryCatch({pkfnd(im = b, th = threshold, sz = NextOdd(diameter))},
                                                 error = function(e) {NULL})
                                  cnt <- tryCatch({cntrd(im = b, mx = pk, sz = NextOdd(diameter))},
                                                  error = function(e) {NULL})

                                  # determine that frame s has at least 1 valid centroid
                                  if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
                                    tmpOUT <- list(cnt = cnt, b = b, j = j)

                                  } else {
                                    #message(paste0("No centroids detectd in frame ",
                                    #               i, " in the current stack"))
                                    #message("Please, check nuclei validation settings for this image stack.")
                                    errCNT <- data.frame(row = 1, col = 1, norm = 1, rg = 1)
                                    tmpOUT <- list(cnt = errCNT[-1, ], b = b, j = j)

                                  }
                                  tmpOUT

                                }, error = (function(e) {
                                  print(e)
                                  try(parallel::stopCluster(cl), silent = TRUE)
                                  return(NULL)
                                }))

    message("Done!", appendLF = TRUE)
    try({suppressWarnings(parallel::stopCluster(cl))}, silent = TRUE)

    re.idx <- order(do.call(c, lapply(all_results, function(x) {x$j})))
    all_results <- all_results[re.idx]
    skpd.frames <- list()
    all_centroids <- lapply(all_results, function(x) {x$cnt})
    all_b <- lapply(all_results, function(x) {x$b})

    # Visualize if needed
    if (show_plots) {
      for (ii in 1:length(all_results)){
        bii <- all_results[[ii]]$b
        cntii <- all_results[[ii]]$cnt
        VisualizeImg(img_mtx = bii, las = 1, main = paste0("Stack num. ", ii))
        VisualizeCntr(centroids = cntii, width_px = ncol(bii), height_px = nrow(bii))
        bii<-NULL; cntii<-NULL
      }
    }
  }

  # Position list (reformated centroid data for track.m input)
  OUT_centroids <- all_centroids
  # Remove columns that contain brightness and sqare of radius of gyration
  # this is the equivalent of position.m function
  # Also, create a frame(tau)label for each array of centroid data
  #for (ti in 1:length(all_centroids)) {
  #
  #    # retain only positional data by removing columns 3 and 4
  #    all_centroids[[ti]] <- all_centroids[[ti]] [, 1:2]
  #    all_centroids[[ti]]$tau <- ti
  #  }
  med.celln.perfr <- tryCatch({median(do.call(c, lapply(all_centroids, nrow)), na.rm = TRUE)},
                              error = function(e) {100})
  max.celln <- med.celln.perfr * 1.55
  min.celln <- med.celln.perfr * 0.45
  keepXX <- do.call(c, lapply(all_centroids, function(xx){
    nrow(xx) >= min.celln && nrow(xx) <= max.celln}))
  all_centroids <- all_centroids[keepXX]

  # Updated to avoid NO cells frames
  all_centroids2 <- list()
  for (ti in 1:length(all_centroids)) {
    ttmp <- all_centroids[[ti]]
    if (nrow(ttmp) > 0) {
      ttmp <- ttmp[, 1:2]
      ttmp$tau <- (length(all_centroids2) + 1)
      all_centroids2[[length(all_centroids2) + 1]] <- ttmp
    }
  }

  # create a matrix that contains centroid data in sequential order by frame(tau)
  #pos <- do.call(rbind, all_centroids)
  pos <- do.call(rbind, all_centroids2)

  tracks2 <- NULL
  if (!is.null(maxDisp)) {
    tracks2 <- tryCatch(track(xyzs = pos, maxdisp = maxDisp, params = track_params), error = function(e) NULL)
  }

  if (is.null(tracks2)) {
    tmp.Area <- tc_obj@images$dim$width_m * tc_obj@images$dim$height_n
    max.disp <- as.integer(as.numeric(sqrt(tmp.Area)) / 5)
    allDisp <- seq(from=max.disp, to = 5, by = -3)
    jj0 <- 1
    while(is.null(tracks2) && jj0 < length(allDisp) ) {
      tracks2 <- tryCatch(
        suppressMessages(track(xyzs = pos, maxdisp = allDisp[jj0], params = track_params)), error = function(e) NULL)
      jj0 <- jj0 + 1
    }
    if (!is.null(tracks2)) {
      track_params$maxDisp <- allDisp[jj0]
      message(paste0("The following maxDisp value was used for this analysis: ", allDisp[jj0]))

    } else {
      message("a reasonable MaxDisp value couldn't be found! Sorry!")
      return(NULL)

    }
  }

  #tracks <- track(xyzs = pos, maxdisp = maxDisp, params = track_params)
  tracks <- tracks2

  # init num of cells
  init.cell.n <- length(unique(tracks[,4]))

  if (!is.null(min_frames_per_cell) &&
      is.numeric(min_frames_per_cell) &&
      length(min_frames_per_cell) == 1 &&
      min_frames_per_cell > 1) {

    all.cids <- table(tracks[, 4])
    all.cids <- data.frame(cell.id = as.numeric(names(all.cids)),
                           count = as.numeric(all.cids))


    all.cids <- all.cids[all.cids$count >= min_frames_per_cell, ]
    keep.cid <- all.cids$cell.id

    tracks <- tracks[ tracks[,4] %in% unique(keep.cid), ]

  } else {
    min_frames_per_cell <- 1
  }
  track_params$min_frames_per_cell <- min_frames_per_cell

  # end num of cells
  end.cell.n <- length(unique(tracks[,4]))

  # message
  message(paste0("Tot num of cells detected in the image stack: ", init.cell.n, "; Cells retained after filtering: ", end.cell.n))

  ### generate tracks
  #tracks <- track(xyzs = pos, maxdisp = maxDisp, params = track_params)

  # pack and return
  #OUT <- list(images = all_b,
  #            centroids = OUT_centroids,
  #            positions = pos,
  #            tracks = tracks,
  #            params = track_params)

  tc_obj@proc_images <- list(images = all_b)
  tc_obj@centroids <- OUT_centroids
  tc_obj@positions <- pos
  tc_obj@tracks <- tracks
  tc_obj@params <- track_params

  tc_obj@ops$track <- 1
  tc_obj@ops$custom_params <- custom_params_flag

  return(tc_obj)
}


#' Get Track Data
#'
#' Extract Track Data from a trackedCells object
#'
#' @param tc_obj a \code{trackedCells} object
#' @param attach_meta logical, shall metaData be attached to tracks
#'
#' @return a data.frame including cell tracks data
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(x)
#' x <- CellTracker(x)
#' getTracks(x)[1:10,]
#' }
#'
#'
#' @export
getTracks <- function(tc_obj, attach_meta = FALSE)
{
  tmp <- tc_obj@tracks
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c("Y", "X", "frame.ID", "cell.ID")
  TMP <- tmp[, c("frame.ID", "X", "Y", "cell.ID")]
  # Adjust as per S request
  TMP <- TMP[, c(4, 2, 3, 1)]
  rownames(TMP) <- NULL
  if(attach_meta && nrow(TMP) > 0) {

    TMP$tiff_file = tc_obj@metadata$tiff_file
    TMP$experiment = tc_obj@metadata$experiment
    TMP$condition = tc_obj@metadata$condition
    TMP$replicate = tc_obj@metadata$replicate
  } else if (nrow(TMP) < 1) {
    return (NULL)
  }

  return(TMP)
}


#' Get Cell population stats
#'
#' Extract cell population statistics from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return data.frame including cell population stats
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(x)
#' x <- CellTracker(x)
#' x <- ComputeTracksStats(x, 10, 10)
#' getPopulationStats(x)
#' }
#'
#'
#' @export
getPopulationStats <- function(tc_obj)
{
  if (tc_obj@ops$stats == 1) {
    return(tc_obj@stats$population)
  } else {
    message("Stats have not been computed yet. Please, run `ComputeTracksStats()`. Thanks.")
  }
}



#' Get Cell migration stats
#'
#' Extract cell migration statistics from a trackedCells object
#'
#' @param tc_obj a \code{trackedCells} object
#'
#' @return data.frame including cell migration stats
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(x)
#' x <- CellTracker(x)
#' x <- ComputeTracksStats(x, 10, 10)
#' getCellsStats(x)[1:10,]
#' }
#'
#'
#'
#' @export
getCellsStats <- function(tc_obj)
{
  if (tc_obj@ops$stats == 1) {
    return(tc_obj@stats$cells)
  } else {
    message("Stats have not been computed yet. Please, run `ComputeTracksStats()`. Thanks.")
  }
}


#' Get MetaData
#'
#' Extract MetaData from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return a list including four items: tiff filename, experiment name, condition label, and replicate ID.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cellmigRation::TrackCellsDataset
#' getCellsMeta(x0)
#'
#' @export
getCellsMeta <- function(tc_obj)
{
  tmp <- tc_obj@metadata
  return(tmp)
}


#' Get Image Stacks
#'
#' Extract Images Stacks from a trackedCells object
#'
#' @param tc_obj a \code{trackedCells} object
#'
#' @return a list including stack images (formatted as numeric matrices)
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cellmigRation::TrackCellsDataset
#' y0 <- getImageStacks(x0)
#' graphics::image(y0[[1]][200:400, 100:300])
#'
#' @export
getImageStacks <- function(tc_obj)
{
  tmp <- tc_obj@images$images
  return(tmp)
}



#' Get Auto Optimized Parameters
#'
#' Extract Parameters that were automatically optimized
#'
#' @param tc_obj a \code{trackedCells} object
#'
#' @return a list including optimized parameter values (lnoise, diameter, and threshold)
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' message("this example may take up to several mins to complete")
#' \dontrun{
#' x <- cellmigRation::TrackCellsDataset
#' x <- OptimizeParams(x)
#' getOptimizedParams(x)
#' }
#'
#'
#' @export
getOptimizedParams <- function(tc_obj)
{
  chk1 <- tc_obj@ops$optimized_params
  y <- list()
  if (!is.null(chk1) && chk1 == 1) {
    y <- tc_obj@optimized$auto_params
  } else {
    message("There were NO optimized params to return")
  }
  return(y)
}


#' Set MetaData
#'
#' Write/Replace MetaData of a trackedCells object
#'
#' @param tc_obj a \code{trackedCells} object
#' @param experiment string, a label to describe the experiment (optional). Can be NULL
#' @param condition string, a label to describe the experimental condition (optional). Can be NULL
#' @param replicate string, a label to identify the replicate (optional). Can be NULL
#'
#' @return a list including three items: experiment name, condition label, and replicate ID.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cellmigRation::TrackCellsDataset
#' x0 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DMSO")
#' getCellsMeta(x0)
#'
#'
#' @export
setCellsMeta <- function(tc_obj, experiment = NULL,
                         condition = NULL, replicate = NULL)
{

  if (is.null(experiment)) {
    experiment <- NA
  } else {
    experiment <- tryCatch(as.character(experiment[1]), error = function(e) NA)
  }

  if (is.null(replicate)) {
    replicate <- NA
  } else {
    replicate <- tryCatch(as.character(replicate[1]), error = function(e) NA)
  }

  if (is.null(condition)) {
    condition <- NA
  } else {
    condition <- tryCatch(as.character(condition[1]), error = function(e) NA)
  }

  FILENM <- tc_obj@metadata$tiff_file

  tmp <- list(tiff_file = FILENM,
              experiment = experiment,
              condition = condition,
              replicate = replicate)
  tc_obj@metadata <- tmp
  return(tc_obj)
}


#' Aggregate trackedCells Objects
#'
#' Aggregate two or more trackedCells-class objects together. Input objects must carry information of cell
#' tracks (otherwise an error will be raised). All tracks form the different experiments/images are returned in a
#' large data.frame. A new unique ID is assigned to specifically identify each cell track from each image/experiment.
#'
#' @param x a \code{trackedCells}-class object where cells have already been tracked
#' @param ... one or more trackedCells-class object(s) where cells have already been tracked
#' @param meta_id_field string, can take one of the following values, c("tiff_file", "experiment",
#' "condition", "replicate"). Indicates the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return An aggregate data.frame including all cells that were tracked over two or more images/experiments.
#' The data.frame includes the following columns: "new.ID", "frame.ID", "X", "Y", "cell.ID", "tiff_name",
#' "experiment", "condition", "replicate". The "new.ID" uniquely identifies a cell in a given image/experiment.
#'
#' @details each trackedCells-class object passed to this function requires a unique identifier (such as a unique
#' tiff_file name). Any of the metadata columns can be used as unique ID for an image/experiment. The function
#' will raise an error if non-unique identifiers are found across the input objects.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom graphics par image
#'
#' @export
aggregateTrackedCells <- function(x, ...,
                                  meta_id_field = c("tiff_file", "experiment",
                                                    "condition", "replicate"))
{
  # Inner fxs
  check_trobj <- function(xx) {
    RT <- FALSE
    if (!is.null(xx)) {
      if ("trackedCells" %in% class(xx)) {
        if (nrow(xx@tracks) > 0) {
          RT <- TRUE
        }
      }
    }
    return(RT)
  }

  compute_mult <- function(xx) {
    zz <- nchar(xx) + 2
    out <- (10 ^ zz)
    return(out)
  }

  meta_id_field <- match.arg(arg = meta_id_field,
                             choices = c("tiff_file", "experiment",
                                         "condition", "replicate"),
                             several.ok = FALSE)
  y <- list(...)
  test1 <- FALSE
  if (length(y) > 0) {
    test1 <- sum(do.call(c, lapply(y, check_trobj))) == length(y)
  }

  # first check
  stopifnot(check_trobj(x), test1)

  big.list <- list(x)
  for(yi in y) {
    big.list[[length(big.list) + 1]] <- yi
  }

  # Chek names are different
  all_ids <- do.call(c, lapply(big.list, function(xx) {xx@metadata[[meta_id_field]] }))
  all_ids <- as.character(all_ids)
  unq_ids <- unique(all_ids)

  # second check
  stopifnot(length(unq_ids) == length(all_ids))

  # Adjust
  my_tracks <- lapply(big.list, getTracks, attach_meta = TRUE)

  my_tracks <- do.call(rbind, my_tracks)
  my_tracks[,"new.ID"] <- factor(my_tracks[,meta_id_field], levels = unq_ids)
  my_tracks[,"new.ID"] <- as.numeric(my_tracks[,"new.ID"])
  my.mult <- compute_mult(max(my_tracks[, "cell.ID"], na.rm = TRUE))
  my_tracks[,"new.ID"] <- (my.mult * my_tracks[,"new.ID"]) + my_tracks[, "cell.ID"]

  # Adjust as per S request
  #keep.colz <- c('new.ID', 'frame.ID', 'X', 'Y', 'cell.ID', 'tiff_file', 'experiment', 'condition', 'replicate')
  keep.colz <- c('new.ID', 'X', 'Y', 'frame.ID', 'cell.ID', 'tiff_file', 'experiment', 'condition', 'replicate')
  out <- my_tracks[, keep.colz]
  rownames(out) <- NULL

  return(out)
}


#' Filter an Aggregated Table of Cell Tracks
#'
#' Filter an Aggregated Table (data.frame) of cell tracks (from multiple images/experiments) and
#' retain cell tracks from images/experiments of interest
#'
#' @param x data.frame, is an aggregated Table of Cell Tracks. Must include the following columns:
#' "new.ID", "frame.ID", "X", "Y", "cell.ID", "tiff_name", "experiment", "condition", "replicate"
#' @param id_list character vector, indicates the IDs (such as tiff_filenames) to be retained in the
#' output data.frame
#' @param meta_id_field string, can take one of the following values, c("tiff_file", "experiment",
#' "condition", "replicate"). Indicates the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return data.frame, a filtered aggregated Table of Cell Tracks
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' A <- data.frame(new.ID = 1:10,
#'                 frame.ID = 10:1,
#'                 X = sample(1:100, size = 10),
#'                 Y = sample(1:100, size = 10),
#'                 cell.ID = c(rep(1, 5), rep(2, 5)),
#'                 tiff_file= c(rep("ii", 3), rep("jj", 5), rep('kk', 2)))
#' FilterTrackedCells(A, id_list = c("ii", "kk"))
#'
#' @export
FilterTrackedCells <- function(x, id_list,
                               meta_id_field = c("tiff_file", "experiment",
                                                 "condition", "replicate")) {

  meta_id_field <- match.arg(arg = meta_id_field,
                             choices = c("tiff_file", "experiment",
                                         "condition", "replicate"),
                             several.ok = FALSE)

  REQd <- c("new.ID", "frame.ID", "X", "Y", "cell.ID", meta_id_field)
  CHK1 <- sum(REQd %in% colnames(x)) == length(REQd)

  stopifnot(CHK1)

  xx <- x[x[, meta_id_field] %in% id_list, ]
  return(xx)
}

#
## --- !!! --- DF part ends here --- !!! ---
#

#' @title Data preprocessing for random migration (RM)
#'
#' @description This function allows preprocessing of the trajectory data from random
#' migration (RM) experiments.
#'
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack. Default is 10 min.
#' @param PixelSize A numeric value of the physical size of a pixel. Default is 1.24.
#' @param FrameN A numeric value of the number of frames. Default is NULL
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @return An CellMig class object with preprocessed data.
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD, FrameN=100)
#' rmTD
#' @export
rmPreProcessing = function(object, PixelSize=1.24,
                           TimeInterval=10, FrameN=NULL) {
  msg <- NULL
  if ( ! is.data.frame(object@trajdata)){
    msg <- c(msg, "input data must be data.frame")
  }
  if ( ! is.numeric(PixelSize)) stop( "PixelSize has to be a positive number" ) else if ( PixelSize<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  object@adjDS <- object@trajdata
  df<-object@adjDS
  df<-df[,1:3]                                        # Removing the unnecessary columns
  spl<-split(df,df[,1])
  cat("This dataset contains:", length(spl), "cell(s) in total\n")

  tbd<-c()
  LENspl<-length(spl)
  for (i in 1: LENspl){
	  if (length(spl[[i]][,1])<4){                     # Removing cells with less than 4 tracked frames
		  tbd<-c(tbd,i)
	  }
  }
  spl[tbd]<-NULL
  df<-do.call(rbind.data.frame, spl)

  L<-length(df[,1])
  df[,4:26]<-rep(0,L)
  df[,27]<-rep(NA,L)                                  # to be used for migration type
  colnames(df)<-c("ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P","Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY","New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F","Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type")
  ID_split <- split(df, df$ID)                        #Splitting the data frame based on the ID
  cat("This dataset contains:", length(ID_split), "cell(s) with more than three steps in their tracks\n")

  for(j in 1:length(ID_split)){                        # Having the ID =group order
    ID_split[[j]][1]=j
  }

  for(j in 1:length(ID_split)){                        # adjusting x and y (starting from 0 & being multiplied by H)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(2:MM, function(i){
      ID_split[[j]][i,4]=PixelSize*((ID_split[[j]][i,2])-( ID_split[[j]][1,2]))      # x2-x1
      ID_split[[j]][i,5]=PixelSize*(( ID_split[[j]][1,3])-(ID_split[[j]][i,3]))      # y2-y1
      return(ID_split[[j]][i,4:5])
    }))
    ID_split[[j]][2:MM,4:5] <- as.data.frame(res)
    ID_split[[j]][,4:5] <- lapply(ID_split[[j]][,4:5], as.numeric)
  }



  for(j in 1:length(ID_split)){                    # removing the old x and y
    ID_split[[j]]=ID_split[[j]][-2]            # removing x column
    ID_split[[j]]=ID_split[[j]][-2]            # removing the y column [-2] is used because x column is gone.
  }



  for(j in 1:length(ID_split)){                    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(1:MM, function(i){
      ID_split[[j]][i,4]=(ID_split[[j]][i+1,2])-( ID_split[[j]][i,2])                                         # creating values for dx
      ID_split[[j]][,4][is.na(ID_split[[j]][,4])] <- 0
      ID_split[[j]][i,5]= ( ID_split[[j]][i+1,3])-(ID_split[[j]][i,3])                                        # creating values for dy
      ID_split[[j]][,5][is.na(ID_split[[j]][,5])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,6]=sqrt((ID_split[[j]][i,4])^2 + (ID_split[[j]][i,5])^2)                                # creating values for dis
      ID_split[[j]][i,7]= acos((ID_split[[j]][i,4])/(ID_split[[j]][i,6]))
      ID_split[[j]][,7][is.na(ID_split[[j]][,7])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,11]=((ID_split[[j]][i,6])/TimeInterval)^2                                               # creating values for Square Speed
      return(ID_split[[j]][i,c(4:7, 11)])
    }))

    ID_split[[j]][1:MM,c(4:7, 11)] <- as.data.frame(res)
    ID_split[[j]][,c(4:7, 11)] <- lapply(ID_split[[j]][,c(4:7, 11)], as.numeric)
  }

  for(j in 1:length(ID_split)){
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res1 <- t(sapply(1:MM, function(i){
      ID_split[[j]][,12]=cumsum(ID_split[[j]][,6])                                                            # creating values for cumsum
      ID_split[[j]][i,13]= sqrt(((ID_split[[j]][i+1,2])^2)+((ID_split[[j]][i+1,3])^2))/(ID_split[[j]][i,12])  # creating values for cumulative directionality ratio
      return(ID_split[[j]][i,12:13])
    }))
    ID_split[[j]][1:MM,12:13] <- as.data.frame(res1)
    ID_split[[j]][,12:13] <- lapply(ID_split[[j]][,12:13], as.numeric)
  }

  for(j in 1:length(ID_split)){              # creating values for  rel.ang.P  (step to the previous)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM1<-MM-1
    res <- sapply(1:MM1, function(i){

      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]>=0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]<0) ){
        ID_split[[j]][i,8]= abs(ID_split[[j]][i+1,7])+abs(ID_split[[j]][i,7])
      }
      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]<0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]>=0) ){
        ID_split[[j]][i,8]=ID_split[[j]][i+1,7]-ID_split[[j]][i,7]
      }
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])<= (-pi), 2*pi+(ID_split[[j]][i,8]),(ID_split[[j]][i,8]))    # adjusting the rel.ang
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])> pi,(ID_split[[j]][i,8])-2*pi,(ID_split[[j]][i,8]))
      return(ID_split[[j]][i, 8])
    })
    ID_split[[j]][1:MM1, 8] <- as.data.frame(res)
  }

  cosine.P<-data.frame()
  for(j in 1:length(ID_split)){              # creating values for  cosine.P  based on rel.ang.P
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      ID_split[[j]][i,9]<-cos(ID_split[[j]][i,8])
      return(ID_split[[j]][i,9])
    })
    ID_split[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-ID_split[[j]][,9]
  }


  for(j in 1:length(ID_split)){              # Computing persistence time   (based on rel.ang.P)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      if(abs(ID_split[[j]][i,8])<=1.5707963268){
        ID_split[[j]][i,10]= TimeInterval
      }
      if(abs(ID_split[[j]][i,8])>1.5707963267){
        ID_split[[j]][i,10]= 0
      }
      return(ID_split[[j]][i,10])
    })
    ID_split[[j]][1:MM, 10] <- as.data.frame(res)
  }


  for(j in 1:length(ID_split)){              # Computing Acceleration
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM2<-MM-2
    res <- sapply(1:MM2, function(i){
      ID_split[[j]][i,24]= (sqrt(ID_split[[j]][i+1,11])-sqrt(ID_split[[j]][i,11]))/TimeInterval
      return(ID_split[[j]][i,24])
    })
    ID_split[[j]][1:MM2, 24] <- as.data.frame(res)
  }


  size<-c()
  for(j in 1:length(ID_split)){
    size[j]<-length(ID_split[[j]][,1])
  }
  S<-summary(size)
  names(S)<-NULL

  if (is.null(FrameN)){
    IncompleteTracks<-subset(size,S<S[6])
    for(j in 1:length(ID_split)){
      ID_split[[j]]<-ID_split[[j]][1:S[1],]
      ID_split[[j]][S[1],4:25]=0
      ID_split[[j]][S[1],10]=TimeInterval
      ID_split[[j]][1,10]=0
    }

    cat("The minimum number of steps: ",S[1],"\n")
    cat("The maximum number of steps: ",S[6],"\n")
    cat("Number of cells with a total number of steps less than ",S[6],"steps",":",length(IncompleteTracks),"\n")
    cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
    PreprocessedData<-ID_split
    object@preprocessedDS<-PreprocessedData

  }else{
    if ( ! is.numeric(FrameN)) stop( "FrameN has to be a positive number" ) else if ( FrameN<= 0 ) stop( "FrameN has to be a positive number" )
    if (FrameN>S[6]) stop( paste0("No cells have ",FrameN, " steps in their tracks"))
    okTracks<-c()
    for(j in 1:length(ID_split)){
      if (length(ID_split[[j]][,1])>=FrameN){
        okTracks<-c(okTracks,j)
      }
    }

    ID_split=ID_split[okTracks]
    for(j in 1:length(ID_split)){
      ID_split[[j]]<-ID_split[[j]][1:FrameN,]
      ID_split[[j]][1:FrameN,1]<-j
      ID_split[[j]][FrameN,4:25]=0
      ID_split[[j]][FrameN,10]=TimeInterval
      ID_split[[j]][1,10]=0
    }
    cat("The desired number of steps: ",FrameN,"\n")
    cat("The maximum number of steps: ",S[6],"\n")
    cat("Only: ", length(ID_split), " cells were selected","\n")

    cat("All the tracks of the selected cells are adjusted to have only ",FrameN," steps","\n")
    PreprocessedData<-ID_split
    object@preprocessedDS<-PreprocessedData

  }
  return(object)
}


#' @title Data preprocessing for wound scratch assay (WSA).
#'
#' @description This function allows filtering of cells and preprocessing of the trajectory data from wound scratch assay (WSA) experiments.
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param PixelSize A numeric value of the physical size of a pixel.
#' @param FrameN A numeric value of the number of frames. Default is NULL
#' @param imageH A numeric value of the image height.
#' @param woundH A numeric value of the image height.
#' @param upperE A numeric value of the upper edge of the wound.
#' @param lowerE A numeric value of the lower edge of the wound.
#' @param mar A numeric value of the margin to be used to narrow the clearing zone inside the zone.
#' @param clearW A logical vector that allows removing the cells within the wound. Default is TRUE.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @return An CellMig class object with filtered, annotated and preprocessed data.
#' @examples
#'
#' data(WSADataset)
#' wasDF=WSADataset[1:1000,]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' @export
wsaPreProcessing = function(object, PixelSize=1.24,
                            TimeInterval=10, FrameN=NULL,
                            imageH=1500, woundH=600,
                            upperE=400, lowerE=1000,
                            mar=75, clearW=TRUE) {
  msg <- NULL
  if ( ! is.data.frame(object@trajdata)){
    msg <- c(msg, "input data must be data.frame")
  }
  if ( ! is.numeric(PixelSize)) stop( "PixelSize has to be a positive number" ) else if ( PixelSize<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(imageH)) stop( "imageH has to be a positive number" ) else if ( imageH<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(woundH)) stop( "woundH has to be a positive number" ) else if ( woundH<= 0 ) stop( "woundH has to be a positive number" )
  if ( ! is.numeric(upperE) ) stop( "upperE has to be a positive number" ) else if ( upperE<= 0 ) stop( "upperE has to be a positive number" )
  if ( ! is.numeric(lowerE) ) stop( "lowerE has to be a positive number" ) else if ( lowerE<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(mar)) stop( "mar has to be a positive number" ) else if ( mar<= 0 ) stop( "mar has to be a positive number" )
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  if (clearW == TRUE){
    dff<-object@trajdata
    dff<-dff[,1:4]                                        # Removing the unnecessary columns
    splitFORu<-split(dff,dff[,1])
    for (i in 1:length(splitFORu)){
      if ((splitFORu[[i]][1,3]>= (upperE + mar) & splitFORu[[i]][1,3]<= (lowerE -mar)) & (splitFORu[[i]][1,4]<=20 )){   # to remove cells within the wound
         splitFORu[[i]]<-NA
      }
    }
    deletedCells<-splitFORu[is.na(splitFORu)]   ########## To get the deleted cells
    cat(paste0(length(deletedCells)," cells were inside the wound and they have been removed"),"\n")
    Ftable<-splitFORu[!is.na(splitFORu)]
    Ftable<-do.call(rbind.data.frame, Ftable)    ########## to convert the list to a data frame
    rownames(Ftable)<-NULL
    object@adjDS <- Ftable
  }else{
    object@adjDS <- object@trajdata
  }

  ####### to set the cells orientation

  hh<- upperE + ((lowerE-upperE)/2)
  CellOr<-c()
  finaltable<- object@adjDS
  finaltable<-split(finaltable,finaltable[,1])
  for (i in 1:length(finaltable)){
    if (finaltable[[i]][1,3]<=hh){
      CellOr[i]=0
    }else{
      CellOr[i]=1
    }
  }
  object@cellpos<-CellOr

  df<-object@adjDS
  df<-df[,1:3]                                        # Removing the unnecessary columns
  spl<-split(df,df[,1])
  cat("This dataset contains:", length(spl), "cell(s) in total\n")

  tbd<-c()
  LENspl<-length(spl)
  for (i in 1: LENspl){
	  if (length(spl[[i]][,1])<4){                     # Removing cells with less than 4 tracked frames
		  tbd<-c(tbd,i)
	  }
  }
  spl[tbd]<-NULL
  df<-do.call(rbind.data.frame, spl)

  L<-length(df[,1])
  df[,4:26]<-rep(0,L)
  df[,27]<-rep(NA,L)                                  # to be used for migration type
  colnames(df)<-c("ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P","Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY","New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F","Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type")
  ID_split <- split(df, df$ID)                        #Splitting the data frame based on the ID
  cat("This dataset contains: ",length(ID_split),"Cells with more than three steps in their tracks","\n")

  for(j in 1:length(ID_split)){                        # Having the ID =group order
    ID_split[[j]][1]=j
  }

  for(j in 1:length(ID_split)){                        # adjusting x and y (starting from 0 & being multiplied by H)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(2:MM, function(i){
      ID_split[[j]][i,4]=PixelSize*((ID_split[[j]][i,2])-( ID_split[[j]][1,2]))      # x2-x1
      ID_split[[j]][i,5]=PixelSize*(( ID_split[[j]][1,3])-(ID_split[[j]][i,3]))      # y2-y1
      return(ID_split[[j]][i,4:5])
    }))
    ID_split[[j]][2:MM,4:5] <- as.data.frame(res)
    ID_split[[j]][,4:5] <- lapply(ID_split[[j]][,4:5], as.numeric)
  }



  for(j in 1:length(ID_split)){                    # removing the old x and y
    ID_split[[j]]=ID_split[[j]][-2]            # removing x column
    ID_split[[j]]=ID_split[[j]][-2]            # removing the y column [-2] is used because x column is gone.
  }



  for(j in 1:length(ID_split)){                    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(1:MM, function(i){
      ID_split[[j]][i,4]=(ID_split[[j]][i+1,2])-( ID_split[[j]][i,2])                                         # creating values for dx
      ID_split[[j]][,4][is.na(ID_split[[j]][,4])] <- 0
      ID_split[[j]][i,5]= ( ID_split[[j]][i+1,3])-(ID_split[[j]][i,3])                                        # creating values for dy
      ID_split[[j]][,5][is.na(ID_split[[j]][,5])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,6]=sqrt((ID_split[[j]][i,4])^2 + (ID_split[[j]][i,5])^2)                                # creating values for dis
      ID_split[[j]][i,7]= acos((ID_split[[j]][i,4])/(ID_split[[j]][i,6]))
      ID_split[[j]][,7][is.na(ID_split[[j]][,7])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,11]=((ID_split[[j]][i,6])/TimeInterval)^2                                               # creating values for Square Speed
      return(ID_split[[j]][i,c(4:7, 11)])
    }))

    ID_split[[j]][1:MM,c(4:7, 11)] <- as.data.frame(res)
    ID_split[[j]][,c(4:7, 11)] <- lapply(ID_split[[j]][,c(4:7, 11)], as.numeric)
  }

  for(j in 1:length(ID_split)){
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res1 <- t(sapply(1:MM, function(i){
      ID_split[[j]][,12]=cumsum(ID_split[[j]][,6])                                                            # creating values for cumsum
      ID_split[[j]][i,13]= sqrt(((ID_split[[j]][i+1,2])^2)+((ID_split[[j]][i+1,3])^2))/(ID_split[[j]][i,12])  # creating values for cumulative directionality ratio
      return(ID_split[[j]][i,12:13])
    }))
    ID_split[[j]][1:MM,12:13] <- as.data.frame(res1)
    ID_split[[j]][,12:13] <- lapply(ID_split[[j]][,12:13], as.numeric)
  }

  for(j in 1:length(ID_split)){              # creating values for  rel.ang.P  (step to the previous)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM1<-MM-1
    res <- sapply(1:MM1, function(i){

      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]>=0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]<0) ){
        ID_split[[j]][i,8]= abs(ID_split[[j]][i+1,7])+abs(ID_split[[j]][i,7])
      }
      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]<0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]>=0) ){
        ID_split[[j]][i,8]=ID_split[[j]][i+1,7]-ID_split[[j]][i,7]
      }
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])<= (-pi), 2*pi+(ID_split[[j]][i,8]),(ID_split[[j]][i,8]))    # adjusting the rel.ang
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])> pi,(ID_split[[j]][i,8])-2*pi,(ID_split[[j]][i,8]))
      return(ID_split[[j]][i, 8])
    })
    ID_split[[j]][1:MM1, 8] <- as.data.frame(res)
  }

  cosine.P<-data.frame()
  for(j in 1:length(ID_split)){              # creating values for  cosine.P  based on rel.ang.P
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      ID_split[[j]][i,9]<-cos(ID_split[[j]][i,8])
      return(ID_split[[j]][i,9])
    })
    ID_split[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-ID_split[[j]][,9]
  }


  for(j in 1:length(ID_split)){              # Computing persistence time   (based on rel.ang.P)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      if(abs(ID_split[[j]][i,8])<=1.5707963268){
        ID_split[[j]][i,10]= TimeInterval
      }
      if(abs(ID_split[[j]][i,8])>1.5707963267){
        ID_split[[j]][i,10]= 0
      }
      return(ID_split[[j]][i,10])
    })
    ID_split[[j]][1:MM, 10] <- as.data.frame(res)
  }


  for(j in 1:length(ID_split)){              # Computing Acceleration
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM2<-MM-2
    res <- sapply(1:MM2, function(i){
      ID_split[[j]][i,24]= (sqrt(ID_split[[j]][i+1,11])-sqrt(ID_split[[j]][i,11]))/TimeInterval
      return(ID_split[[j]][i,24])
    })
    ID_split[[j]][1:MM2, 24] <- as.data.frame(res)
  }
  size<-c()
  for(j in 1:length(ID_split)){
    size[j]<-length(ID_split[[j]][,1])
  }
  S<-summary(size)
  names(S)<-NULL

  if (is.null(FrameN)){
    IncompleteTracks<-subset(size,S<S[6])
    for(j in 1:length(ID_split)){
      ID_split[[j]]<-ID_split[[j]][1:S[1],]
      ID_split[[j]][S[1],4:25]=0
      ID_split[[j]][S[1],10]=TimeInterval
      ID_split[[j]][1,10]=0
    }

    cat("The minimum number of steps: ",S[1],"\n")
    cat("The maximum number of steps: ",S[6],"\n")
    cat("Number of cells with a total number of steps less than ",S[6],"steps",":",length(IncompleteTracks),"\n")
    cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
    PreprocessedData<-ID_split
    object@preprocessedDS<-PreprocessedData

  }else{
    if ( ! is.numeric(FrameN)) stop( "FrameN has to be a positive number" ) else if ( FrameN<= 0 ) stop( "FrameN has to be a positive number" )
    if (FrameN>S[6]) stop( paste0("No cells have ",FrameN, " steps in their tracks"))
    okTracks<-c()
    for(j in 1:length(ID_split)){
      if (length(ID_split[[j]][,1])>=FrameN){
        okTracks<-c(okTracks,j)
      }
    }

    ID_split=ID_split[okTracks]
    for(j in 1:length(ID_split)){
      ID_split[[j]]<-ID_split[[j]][1:FrameN,]
      ID_split[[j]][1:FrameN,1]<-j
      ID_split[[j]][FrameN,4:25]=0
      ID_split[[j]][FrameN,10]=TimeInterval
      ID_split[[j]][1,10]=0
    }
    cat("The desired number of steps: ",FrameN,"\n")
    cat("The maximum number of steps: ",S[6],"\n")
    cat("Only: ", length(ID_split), " cells were selected","\n")

    cat("All the tracks of the selected cells are adjusted to have only ",FrameN," steps","\n")
    PreprocessedData<-ID_split
    object@preprocessedDS<-PreprocessedData

  }
  return(object)
}

#' @title A 2D rose-plot
#'
#' @description Plotting the trajectory data of all cells.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param Type has to be one of the following: c("p", "l", "b", "o")
#' "p": Points
#' "l": Lines
#' "b": Both
#' "o": Both "overplotted"
#' @param FixedField logical(1) Allows generating a plot with fixed field 800um x 800um. Default is TRUE.
#' @param export if `TRUE` (default), exports plot to JPG file
#' @return A 2D rose-plot showing the tracks of all cells.
#' @details  The visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' plotAllTracks(rmTD, ExpName="Test",Type="l", FixedField=TRUE,export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points lines
#' @export
plotAllTracks= function(object, ExpName="ExpName", Type="l", FixedField=TRUE, export=FALSE) {

  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }

  Len<-length(Object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <- grDevices::rainbow(1023)
    colo2 <- grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }

  if ( FixedField == TRUE){
	graphics::plot(Object[[1]][1:Step,2], Object[[1]][1:Step,3],
                 type=Type, xlab="X (um)", ylab="Y (um)",
                 col=color[1], las=1, xlim=c(-400,400),cex.lab=0.7,
                 ylim=c(-400,400), main=ExpName)
  	for(n in 2:Len){
    		points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(-500,500)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(-500,500)
  	graphics::lines(x, y, type='l', col="black")
	if (export) grDevices::jpeg(paste0(ExpName,"_All_tracks_plot.jpg"),width = 4, height = 4, units = 'in', res = 300)
  		graphics::plot(Object[[1]][1:Step,2], Object[[1]][1:Step,3], type=Type,
                 xlab="X (um)", ylab="Y (um)", col=color[1],
                 las=1, xlim=c(-400,400), ylim=c(-400,400), main=ExpName)
  	for(n in 2:Len){
    		graphics::points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(-500,500)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(-500,500)
  	graphics::lines(x, y, type='l', col="black")
  }else{
	MinX<-c()
  	MaxX<-c()
  	MinY<-c()
  	MaxY<-c()
  	for(j in 1:Len){
    		minX=min(Object[[j]][1:Step,2])
    		minY=min(Object[[j]][1:Step,3])
    		maxX=max(Object[[j]][1:Step,2])
    		maxY=max(Object[[j]][1:Step,3])
    		MinX[j]<-c(minX)
    		MaxX[j]<-c(maxX)
    		MinY[j]<-c(minY)
    		MaxY[j]<-c(maxY)
  	}
 	RangeX=c(MinX,MaxX)
  	RangeY=c(MinY,MaxY)
  	graphics::plot(Object[[1]][1:Step,2], Object[[1]][1:Step,3],
                 type=Type, xlab="X (um)", ylab="Y (um)",
                 col=color[1], las=1, xlim=range(RangeX),cex.lab=0.7,
                 ylim=range(RangeY), main=ExpName)
  	for(n in 2:Len){
    		points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(min(RangeX)-100,max(RangeX)+100)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(min(RangeY)-100,max(RangeY)+100)
  	graphics::lines(x, y, type='l', col="black")

  	if (export) grDevices::jpeg(paste0(ExpName,"_All_tracks_plot.jpg"),width = 4, height = 4, units = 'in', res = 300)
  	graphics::plot(Object[[1]][1:Step,2], Object[[1]][1:Step,3], type=Type,
                 xlab="X (um)", ylab="Y (um)", col=color[1],
                 las=1, xlim=range(RangeX), ylim=range(RangeY), main=ExpName)
  	for(n in 2:Len){
    		graphics::points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(min(RangeX)-100,max(RangeX)+100)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(min(RangeY)-100,max(RangeY)+100)
  	graphics::lines(x, y, type='l', col="black")
    }
    if (export) {
    grDevices::dev.off()
    cat("The plot is saved in your directory [use getwd()]","\n")
  }
}

#' @title A 2D rose-plot of sample cells
#'
#' @description Plotting the trajectory data of some cells.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param Type has to be one of the following: c("p", "l", "b", "o")
#' @param export if `TRUE` (default), exports plot to JPG file
#' "p": Points
#' "l": Lines
#' "b": Both
#' "o": Both "overplotted"
#' @param FixedField logical(1) Allows generating a plot with fixed field 800um x 800um. Default is TRUE.
#' @param celNum A numeric value showing the desired number of cells to be plotted.
#' @return A 2D rose-plot showing the tracks of sample cells selected randomly based on the desired number of cells selected by the user.
#' @details  The visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' plotSampleTracks(rmTD, ExpName="Test",Type="l", FixedField=TRUE, celNum=5, export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points lines
#' @export
plotSampleTracks= function(object, ExpName="ExpName", Type="l", celNum=35,FixedField=TRUE,export=FALSE) {
 if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")

  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  Len<-length(Object)
  if ( ! is.numeric(celNum) ) stop( "celNum has to be a positive number" ) else if ( celNum > Len ) stop( "The cellNum should be less than the total number of cells" )
  Step<-length(Object[[1]][,1])

  color <-c()
  if (Len> 1023){
     	colnum= Len-1023
 	color1 <- grDevices::rainbow(1023)
  	colo2 <- grDevices::rainbow(colnum)
  	color=c(color1 ,colo2)
  }else{
  	color <- grDevices::rainbow(Len)
  }
  OBJ<-c(1:Len)
  cells=sample(OBJ,celNum)
  cells=sort(cells)
  cat(paste0("The plot contains the following cells: "),"\n")
  cat(cells,"\n")

  if ( FixedField == TRUE){
	graphics::plot(Object[[cells[1]]][1:Step,2], Object[[cells[1]]][1:Step,3],
                 type=Type, xlab="X (um)", ylab="Y (um)",
                 col=color[1], las=1, xlim=c(-400,400),cex.lab=0.7,
                 ylim=c(-400,400), main=ExpName)
      CELLS<-cells[-1]
  	for(n in CELLS){
    		points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(-500,500)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(-500,500)
  	graphics::lines(x, y, type='l', col="black")
	if (export) grDevices::jpeg(paste0(ExpName,"_Sample_tracks_plot.jpg"),width = 4, height = 4, units = 'in', res = 300)
  		graphics::plot(Object[[cells[1]]][1:Step,2], Object[[cells[1]]][1:Step,3],type=Type,
                 xlab="X (um)", ylab="Y (um)", col=color[1],
                 las=1, xlim=c(-400,400), ylim=c(-400,400), main=ExpName)
      CELLS<-cells[-1]
  	for(n in CELLS){
    		graphics::points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(-500,500)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(-500,500)
  	graphics::lines(x, y, type='l', col="black")
  }else{
	MinX<-c()
  	MaxX<-c()
  	MinY<-c()
  	MaxY<-c()
  	for(j in cells){
    		minX=min(Object[[j]][1:Step,2])
    		minY=min(Object[[j]][1:Step,3])
    		maxX=max(Object[[j]][1:Step,2])
    		maxY=max(Object[[j]][1:Step,3])
    		MinX<-c(MinX,minX)
    		MaxX<-c(MaxX,maxX)
    		MinY<-c(MinY,minY)
    		MaxY<-c(MaxY,maxY)
  	}
 	RangeX=c(MinX,MaxX)
  	RangeY=c(MinY,MaxY)
  	graphics::plot(Object[[cells[1]]][1:Step,2], Object[[cells[1]]][1:Step,3],
                 type=Type, xlab="X (um)", ylab="Y (um)",
                 col=color[1], las=1, xlim=range(RangeX),cex.lab=0.7,
                 ylim=range(RangeY), main=ExpName)
      CELLS<-cells[-1]
 	for(n in CELLS){
    		points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(min(RangeX)-100,max(RangeX)+100)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(min(RangeY)-100,max(RangeY)+100)
  	graphics::lines(x, y, type='l', col="black")

  	if (export) grDevices::jpeg(paste0(ExpName,"_Sample_tracks_plot.jpg"),width = 4, height = 4, units = 'in', res = 300)
  	graphics::plot(Object[[cells[1]]][1:Step,2], Object[[cells[1]]][1:Step,3], type=Type,
                 xlab="X (um)", ylab="Y (um)", col=color[1],
                 las=1, xlim=range(RangeX), ylim=range(RangeY), main=ExpName)
      CELLS<-cells[-1]
  	for(n in CELLS){
    		graphics::points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    		end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    		graphics::points(end,pch=16,col=color[n], cex = 0.6)
  	}
  	x=c(min(RangeX)-100,max(RangeX)+100)
  	y=c(0,0)
  	graphics::lines(x, y, type='l', col="black")
  	x=c(0,0)
  	y=c(min(RangeY)-100,max(RangeY)+100)
  	graphics::lines(x, y, type='l', col="black")
    }
    if (export) {
    grDevices::dev.off()
    cat("The plot is saved in your directory [use getwd()]","\n")
  }
}


#' @title A 3D rose-plot of all cells
#'
#' @description Plotting the trajectory data of all cells in 3D.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#'
#'
#' @return A 3D rose-plot showing the tracks of all cells.
#' @details  The 3D visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
#'
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' \dontrun{
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' plot3DAllTracks(rmTD, VS=3, size=2)
#' }
#'
#' @importFrom grDevices rainbow
#' @importFrom rgl plot3d
#'
#' @export
plot3DAllTracks= function(object, VS=3, size=2) {

  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }

  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )

  plotTable<-data.frame()
  Len<-length(Object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <- grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }
  col=c()
  coll=c()
  for(i in 1:Len){
    firstNum=((i*Step)-Step)+1
    RowNum=c(firstNum:(i*Step))
    plotTable[RowNum,1]= Object[[i]][1:Step,2]
    plotTable[RowNum,2]=Object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
  }
  rgl::plot3d(plotTable, col=coll, type="p", size=size, axes=FALSE,xlab=" ", ylab=" ",zlab=" ")
}


#' @title A 3D rose-plot
#' @description Plotting the trajectory data of particular cells in 3D.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#' @param cells A numeric vector containing the cell's numbers to be plotted.
#'
#' @return A 3D rose-plot showing the tracks of particular cells.
#' @details  The 3D visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' plot3DTracks(rmTD, VS=3, size=2, cells=1:5)
#'
#' @importFrom grDevices rainbow
#' @importFrom rgl plot3d
#' @importFrom stats complete.cases
#'
#' @export
plot3DTracks= function(object, VS=3, size=2, cells) {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )

  plotTable<-data.frame()
  Len<-length(Object)
  cat(paste0("The plot contains ",length(cells), " Cells"),"\n")
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-grDevices::rainbow(Len)
  }
  cells=sort(cells)
  col=c()
  coll=c()
  for(i in cells){
    firstNum=((i*Step)-Step)+1
    RowNum=c(firstNum:(i*Step))
    plotTable[RowNum,1]= Object[[i]][1:Step,2]
    plotTable[RowNum,2]=Object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
    NewplotTable<-plotTable[stats::complete.cases(plotTable),]

  }
  rgl::plot3d(NewplotTable, col=coll, type="p", size=size, axes=FALSE,xlab=" ", ylab=" ",zlab=" ")
}


#' @title A graphical display of the track of each cell.
#' @description Plotting the trajectory data of each cell.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param Type has to be one of the following: [p, l, b, o]
#' "p": Points
#' "l": Lines
#' "b": Both
#' "o": Both "overplotted"
#' @param FixedField logical(1) Allows generating individual plots with fixed field. Default is TRUE.
#' @param export if `TRUE` (default), exports plot to JPG file
#'
#' @return 2D rose-plots of the cells' track Separately.
#' @details  The visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' PlotTracksSeparately(rmTD, ExpName="Test",
#'                      Type="b", FixedField=FALSE, export = FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points title
#'
#'@export
PlotTracksSeparately= function(object, ExpName="ExpName",
                               Type="l", FixedField=TRUE, export=FALSE) {

  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  Len<-length(Object)
  cat(paste0(Len," plots will be generated in a folder called:",ExpName,"_Tracks","\n"))
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <- grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-grDevices::rainbow(Len)
  }
  if (export) dir.create(paste0(ExpName,"_Tracks"))
  d=getwd()
  if (export) setwd(paste0(d,"/",paste0(ExpName,"_Tracks")))

  if ( FixedField == TRUE){
    MinX<-c()
    MaxX<-c()
    MinY<-c()
    MaxY<-c()
    for(j in 1:Len){
      minX=min(Object[[j]][1:Step,2])
      minY=min(Object[[j]][1:Step,3])
      maxX=max(Object[[j]][1:Step,2])
      maxY=max(Object[[j]][1:Step,3])
      MinX[j]<-c(minX)
      MaxX[j]<-c(maxX)
      MinY[j]<-c(minY)
      MaxY[j]<-c(maxY)
    }
    RangeX=c(MinX,MaxX)
    RangeY=c(MinY,MaxY)
    for(n in 1:Len){
      if (export) grDevices::jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      graphics::plot(Object[[n]][1:Step,2], Object[[n]][1:Step,3],
                     type=Type, xlab="x (um)", ylab="y (um)",
                     col=color[n], las=1, xlim=range(RangeX), ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      graphics::lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      graphics::lines(x, y, type='l', col="black")
      end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
      graphics::points(end,pch=16,col=color[n], cex = 1)
      graphics::title(main=paste0("Cell Number  ", n),col.main="black")
      if (export) grDevices::dev.off()
    }
  }else{
    for(n in 1:Len){
      RangeX= Object[[n]][1:Step,2]
      RangeY= Object[[n]][1:Step,3]
      if (export) grDevices::jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      graphics::plot(Object[[n]][1:Step,2],Object[[n]][1:Step,3],type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,xlim=range(RangeX),ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      graphics::lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      graphics::lines(x, y, type='l', col="black")
      end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
      graphics::points(end,pch=16,col=color[n], cex = 1)
      graphics::title(main=paste0("Cell Number  ", n),col.main="black")
      if (export) grDevices::dev.off()
    }

  }
  setwd(d)
}

#' @title Persistence and Speed
#' @description The PerAndSpeed() generates data and plots for persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param PtSplot A logical vector that allows generating individual plots of persistence time vs speed per cell. Default is TRUE.
#' @param AllPtSplot  A logical vector that allows generating a plot of persistence time vs speed for all cells. Default is TRUE.
#' @param ApSplot A logical vector that allows generating individual plots of angular persistence vs speed per cell. Default is TRUE.
#' @param AllApSplot A logical vector that allows generating a plot of angular persistence vs speed of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#'
#' @return An CellMig class object with a data frame and plots. The data frame is stored in the PerAanSpeedtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- PerAndSpeed(rmTD,TimeInterval=10,ExpName="ExpName", export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off chull
#' @importFrom stats cor.test lm median sd
#' @importFrom graphics plot title abline legend
#' @importFrom sp Polygon
#' @importFrom matrixStats rowMedians
#' @importFrom vioplot vioplot
#' @importFrom utils write.csv
#'
#'
#'@export
PerAndSpeed= function(object, TimeInterval=10,
                      ExpName="ExpName", PtSplot=TRUE,
                      AllPtSplot=TRUE, ApSplot=TRUE,
                      AllApSplot=TRUE, export=FALSE) {

  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  d=getwd()
  if (export) {
    dir.create(paste0(ExpName,"-PerResults"))
    setwd(paste0(d,"/",paste0(ExpName,"-PerResults")))
  }

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }
  PerResultsTable<-data.frame()                             # creating a table to store the persistence results
  for(j in 1:length(Object)){                               # creating values (NA and persistence time )for  persistence    (step to the original)
    MM<-length(Object[[j]][,1])
    MM2<-MM-1
    Ptime<-Object[[j]][1:MM2,10]
    MeanPerTime<-round(mean(Ptime),digits=2)             # computing mean persistence time
    PerTimLen<-Ptime                                     # computing the number of persistence steps to be used in computing the persistence ratio
    PerTimLen[PerTimLen==0]<-NA
    PerTimLen<-PerTimLen[!is.na(PerTimLen)]
    PerTimLen<-length(PerTimLen)
    PerRatio<-round(PerTimLen/MM2, digits=2)              # computing persistence ratio
    DD.Ptime<-Object[[j]][1:MM,10]                        # computing the direction deviating time
    DD.Ptime[DD.Ptime==0]<-1                              # replacing the 0 with 5 and the 5 with 0
    DD.Ptime[DD.Ptime==TimeInterval]<-0
    DD.Ptime[DD.Ptime==1]<-TimeInterval
    MeanDD.PerTime<-round(mean(DD.Ptime),digits=2)

    PerResultsTable[1,j]<-j
    PerResultsTable[2,j]<-MeanPerTime
    PerResultsTable[3,j]<-MeanDD.PerTime
    PerResultsTable[4,j]<-PerRatio
  }

  VelPerTable<-data.frame()               # creating a table to store the mean velocity with correspondence with the persistence time
  for(j in 1:length(Object)){
    MM=length(Object[[j]][,1])
    Ptime<-Object[[j]][1:MM,10]
    Ptime0<-c(0,Ptime)                      #adding a "0" value in the beginning
    Ptime0[Ptime0==TimeInterval]<-NA        #replacing the "TimeInterval" values with NAs
    Ptime0TF<- !is.na(Ptime0)
    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval <- rowN[Ptime0TF]

    fillIdx <- cumsum(Ptime0TF)
    s<-Tval[fillIdx]                     #replacing the NAs with the previous true value
    justTrueVal<-Tval
    s[Ptime0TF]=0                        #replacing the true values with 0
    finalID<-s[-1]                       # removing the first added value

    Vel<-Object[[j]][1:MM,11]
    Per<-Object[[j]][1:MM,10]
    Per[Per==0]<-NA
    tab<-data.frame(Vel,Per,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]         # to exclude the last row since it has no velocity
    tabb<-split(tab, finalIDD)     #to avoid taking the last row

    meanVEL<-c()
    Per.time<-c()
    res1 <- sapply(1:length(tabb), function(i){
      meanVEL= round(sqrt(mean(tabb[[i]][,1])),digits=2)
      return(meanVEL)
    })
    meanVEL= res1
    meanVEL=meanVEL[-1]

    res2 <- sapply(1:length(tabb), function(i){
      Per.time=sum(tabb[[i]][,2])
      return(Per.time)
    })
    Per.time= res2
    Per.time=Per.time[-1]
    w<-which.max(Per.time)
    PerResultsTable[5,j]<-Per.time[w]


    Ptime<-Object[[j]][1:MM,10]           # computing the "meanVel.for0"
    Ptime00<-c(1,Ptime)                   # adding a "1" value in the beginning  this value should not be 0 because we will not catch 0 persistence if it is in the beginning.
    Ptime00[Ptime00==0]<-NA               # replacing the "0" values with NAs
    Ptime00TF<- !is.na(Ptime00)
    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval0 <- rowN[Ptime00TF]
    #length(Tval0)

    TTT<-c()
    res3 <- sapply(1:(length (Tval0)-1), function(i){
      TTT<- (Tval0[i+1]- Tval0[i])-1
      return(TTT)
    })
    TTT=res3
    TTT[TTT==0]<-NA
    F.ID.for.0.per<-TTT[!is.na(TTT)]
    Final.ID.for.0.per<-c()
    for (i in 1:length(F.ID.for.0.per)){
      Final.ID.for.0.per<-c(Final.ID.for.0.per,rep(i,(F.ID.for.0.per[i])))
    }

    Vel<-Object[[j]][1:MM,11]
    Per<-Object[[j]][1:MM,10]
    Per[Per==0]<-NA
    tab<-data.frame()
    tab<-data.frame(Vel,Per,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]
    tabb<-split(tab, finalIDD)     #to avoid taking the last row
    t<-tabb[[1]]

    tabb1<-cbind(t,Final.ID.for.0.per)
    tabbb<-split(tabb1, Final.ID.for.0.per)
    meanVel.for0<-c()

    res4 <- sapply(1:length(tabbb), function(i){
      meanVel.for0<-round(sqrt(mean(tabbb[[i]][,1])),digits=2)
      return(meanVel.for0)
    })

    meanVel.for0=res4
    zerooPrep<-rep(0,length(meanVel.for0))
    Per.time1<-c(zerooPrep,Per.time)
    meanVEL1<-c(meanVel.for0,meanVEL)
    colNumP<-j+j-1
    colNumV<-j+j
    rowNum<-c(1:length(meanVEL1))
    VelPerTable[rowNum,colNumP]<-Per.time1
    VelPerTable[rowNum,colNumV]<-meanVEL1
    PT<-VelPerTable[,j+j]
    c<-stats::cor.test( ~ VelPerTable[,j+j-1]+ PT, method = "spearman",exact=FALSE)                 #testing the correlation
    cc<-unlist(c[4])
    ccPV<-round(cc, digits = 3)
    PerResultsTable[6,j]<-ccPV               # Speed vs  persistence time  Spearman correlation


    if ( PtSplot == TRUE){
      if (export) {
        grDevices::jpeg(paste0(ExpName," Persist Time vs Speed",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      }
      graphics::plot(VelPerTable[,j+j],VelPerTable[,j+j-1],
                     pch=16,type="p",ylab="Persistence Time (min)",
                     xlab=" Mean Speed during persistence time  (um/min)",col=color[j],las=1)
      reg<-stats::lm(PT~VelPerTable[,j+j-1])
      #abline(reg,untf=FALSE,col="black")
      graphics::title(main=paste0("Cell Number  ", j,"   Speed vs Persistence Time"),cex.main =0.7 ,sub=paste0("Spearman's rank correlation coefficient = ",ccPV),col.sub="red")
      if (export) grDevices::dev.off()

    }

  }

  ## All cells (Persistence times  vs Speed)
  allper<-VelPerTable[,1]
  allvel<-VelPerTable[,2]
  for(j in 1:length(Object)){
    allper<-c(allper,VelPerTable[,j+j-1])
    allvel<-c(allvel,VelPerTable[,j+j])
  }
  allper<-allper[!is.na(allper)]
  allvel<-allvel[!is.na(allvel)]
  allvel<-(allvel/TimeInterval)*60
  all.per.vel.table<-data.frame(allper,allvel)
  reg<-stats::lm(allper~allvel)
  c<-stats::cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)                 #testing the correlation
  cc<-unlist(c[4])
  ccP<-round(cc, digits = 3)
  PerResultsTable[6,(length(Object)+1)]<-ccP                                          # Speed vs  persistence time  Spearman correlation

  if ( AllPtSplot == TRUE){
    if (export) {
      grDevices::jpeg(paste0(ExpName,"_Persist_Time_vs_Speed-All_Cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
    }
    graphics::plot(allvel,allper,type="p",pch=16,ylab="Persist_Time (min)",xlab=" Mean Speed during persistence time (um/h)",col="black",las=1)
    graphics::abline(reg,untf=FALSE,col="red")
    graphics::title("Speed vs Persist Time (All cells)",cex.main = 0.7,sub=paste0("Spearman's rank correlation coefficient = ",ccP),col.sub="red")
    if (export) grDevices::dev.off()

  }

  for(j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    MM<-length(Object[[j]][,1])
    MM2<-MM-1
    Root.Median.Square.Speed<-round(sqrt(stats::median(Object[[j]][1:MM2,11])),digits = 3)*60
    PerResultsTable[7,j]<-Root.Median.Square.Speed

    wma<-which.max(sqrt(Object[[j]][,11]))
    wmi<-which.min(sqrt(Object[[j]][1:MM2,11]))

    PerResultsTable[8,j]<-round(sqrt(Object[[j]][wma,11]),digits=3)* 60
    PerResultsTable[9,j]<-round(sqrt(Object[[j]][wmi,11]),digits=3)* 60

    mean.cosineP<-round(mean(Object[[j]][1:(MM2-1),9],na.rm = TRUE),digits = 3)
    PerResultsTable[10,j]<-mean.cosineP
    s<-stats::cor.test( ~ sqrt(Object[[j]][1:MM2,11])+ Object[[j]][1:MM2,9],
                        method = "spearman",exact=FALSE)                 #testing the correlation
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    PerResultsTable[11,j]<-VEvsCOSP

    data<-Object[[j]][1:MM2,2:3]
    data1<-Object[[j]][1:round(MM2/4),2:3]
    data2<-Object[[j]][round(MM2/4):round(MM2/2),2:3]
    data3<-Object[[j]][round(MM2/2):round(MM2*3/4),2:3]
    data4<-Object[[j]][round(MM2*3/4):MM2,2:3]

    ch  <- grDevices::chull(data)
    ch1 <- grDevices::chull(data1)
    ch2 <- grDevices::chull(data2)
    ch3 <- grDevices::chull(data3)
    ch4 <- grDevices::chull(data4)

    coords  <- data[c(ch, ch[1]), ]  # closed polygon
    coords1 <- data1[c(ch1, ch1[1]), ]  # closed polygon
    coords2 <- data2[c(ch2, ch2[1]), ]  # closed polygon
    coords3 <- data3[c(ch3, ch3[1]), ]  # closed polygon
    coords4 <- data4[c(ch4, ch4[1]), ]  # closed polygon

    p = sp::Polygon(coords)
    p1 = sp::Polygon(coords1)
    p2 = sp::Polygon(coords2)
    p3 = sp::Polygon(coords3)
    p4 = sp::Polygon(coords4)

    EmptyArea=abs(p@area -(p1@area + p2@area + p3@area + p4@area))
    SegmentedCA<-p1@area + p2@area + p3@area + p4@area

    PerResultsTable[12,j]=round(p@area,digits=3)
    PerResultsTable[13,j]=round(SegmentedCA,digits=3)
    PerResultsTable[14,j]=round(EmptyArea,digits=3)


    PerResultsTable[15,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3)
    PerResultsTable[16,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3) - abs(round(sum(Object[[j]][,8])/6.28,digits=3))


    if ( ApSplot == TRUE){
      if (export) {
        grDevices::jpeg(
          paste0(ExpName,"_Angular_Persistence_vs_Speed",j,".jpg",width = 4, height = 4, units = 'in', res = 300)
        )
      }
      Speed=(sqrt(Object[[j]][1:MM2,11])/TimeInterval)*60
      graphics::plot(Speed,Object[[j]][1:MM2,9],pch=16,type="p",ylab="Angular Persistence (cosine)",xlab=" Instantaneous Speed (um/h)",col=color[j],las=1, xlim=c(0,3))
      reg<-stats::lm(Object[[j]][1:MM2,9]~Speed)
      graphics::abline(reg,untf=FALSE,col="black")
      graphics::title(main=paste0("Cell Number  ", j," Instantaneous Speeds vs Angular Persistence "),cex.main = 0.7,sub=paste0("Spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
      if (export) grDevices::dev.off()
    }

    PoCos<- subset(Object[[j]][1:(MM2-1),9],Object[[j]][1:(MM2-1),9]>0)
    NeCos<- subset(Object[[j]][1:(MM2-1),9],Object[[j]][1:(MM2-1),9]<=0)

    PerResultsTable[17,j]<-round(stats::median(PoCos,na.rm = TRUE),digits = 3)
    PerResultsTable[18,j]<-round(stats::median(NeCos,na.rm = TRUE),digits = 3)

    AvSp<-(Object[[j]][1:length(Object[[j]][,1])-1,6]/TimeInterval)*60

    PerResultsTable[19,j]<-round(stats::median(AvSp,na.rm = TRUE),digits = 3)
    s=summary(AvSp)
    PerResultsTable[20,j]<-round((s[5]-s[2])/s[3],digits=4)
    PerResultsTable[21,j]<-round(mean(AvSp),digits = 3)
    PerResultsTable[22,j]<-round(stats::sd(AvSp),digits = 3)
  }

  #(all cells)  persistence vs inst.speed
  MM<-length(Object[[1]][,1])
  MM2<-MM-1
  cosine.P<-data.frame()
  for(j in 1:length(Object)){              # creating values for  cosine.P  based on rel.ang.P
    M<- Object[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      Object[[j]][i,9]<-cos(Object[[j]][i,8])
      return(Object[[j]][i,9])
    })
    Object[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-Object[[j]][,9]
  }

  RM<-round(matrixStats::rowMedians(as.matrix(cosine.P[1:MM2,]),na.rm = TRUE),digits=3)
  Speed<-data.frame()
  for (j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    Speed[1:MM2,j]<-round(sqrt(Object[[j]][1:MM2,11]),digits = 3)
  }
  RowmeanSpeed<-round(matrixStats::rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-stats::cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  VEvsCOSP<-round(ss, digits = 3)
  PerResultsTable[11,(length(Object)+1)]<-VEvsCOSP

  if ( AllApSplot == TRUE){
    if (export) {
      grDevices::jpeg(
        paste0(
          ExpName,
          " All_Cells_Average_Angular_Persistence_vs_Average_Speed.jpg"
        ),width = 4, height = 4, units = 'in', res = 300)
    }
    MS<-max(RowmeanSpeed)*60
    graphics::plot(RowmeanSpeed*60,RM,pch=16,type="p",ylab="Average Angular Persistence (cosine)",xlab=" Average Instantaneous Speed (um/h)",col="black",las=1,xlim=c(0,MS))
    NewSpeed=RowmeanSpeed*60
    reg<-stats::lm(RM~NewSpeed)
    graphics::abline(reg,untf=FALSE,col="red")
    graphics::title(main=paste0("All Cells Instantaneous Speed vs Angular Persistence "),cex.main = 0.7,
                    sub=paste0("Spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
    if (export) grDevices::dev.off()
  }

  rownames(PerResultsTable)<-c("Cell Number","Mean Persist Time (min)","Mean persist Deviating Time (min)","Persistence Ratio",
                               "Maximum Persistence period (min)","Persistence Time vs Speed (SCC)","RMSS (um per h)","Maximum Speed (um per h)","Minimum Speed (um per h)",
                               "Mean Angular Persistence (cosine)","Instantaneous Speed vs Angular Persistence (SCC)","Covered Area (um2)","Segmented Covered Area (um2)","Empty Area (um2)","Number of complete rotations",
                               "Number of canceled rotations","Mean Persistence Angle (cosine)","Mean Deviating Angle (cosine)","Median Speed","Speed QBCV",
                               "Mean Speed (um per h)","Speed standard deviation (um per h)")



  RM1<-round(matrixStats::rowMedians(as.matrix(PerResultsTable),na.rm = TRUE),digits=3)
  PerResultsTable[c(2:5,7:10,12:22),(length(Object)+1)]<-RM1[c(2:5,7:10,12:22)]

  RMSS<-as.numeric(PerResultsTable[7,1:length(PerResultsTable[1,])-1])
  if (export) grDevices::jpeg(paste0(ExpName,"_RMSS_profile_of_all_cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  cells<-c(1:(length(PerResultsTable[1,])-1))
  MS<-max(RMSS)
  graphics::plot(cells,RMSS,pch=16,type="o",ylab = 'RMSS(um/h)',xlab = 'Cells',las=1,ylim = c(0, MS))
  graphics::title(main="RMSS Profile of all cells",cex.main = 1)
  graphics::abline(h=stats::median(RMSS[which(!is.na(RMSS))]),col="red")
  graphics::abline(h=mean(RMSS[which(!is.na(RMSS))]),col="blue")
  graphics::legend(1, y=200, legend=c("Mean RMSSs","Median RMSSs"), col=c("blue","red"),lty=1, cex=0.8)
  if (export) grDevices::dev.off()

  if (export) {
    grDevices::jpeg(paste0(ExpName,"_RMSS_ViolinPlot_of_all_cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  }
  graphics::plot(1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '', ylab = 'RMSS(um/h)', xaxt = 'n',las=1)
  graphics::title("RMSS of all cells",cex.main = 1)
  vioplot::vioplot(RMSS, at = 1, add = TRUE, col = "gray")
  if (export) grDevices::dev.off()


  SPEED<-as.numeric(PerResultsTable[19,1:length(PerResultsTable[1,])-1])
  if (export) grDevices::jpeg(paste0(ExpName,"_Speed_profile_of_all_cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  cells<-c(1:(length(PerResultsTable[1,])-1))
  MS<-max(SPEED)
  graphics::plot(cells,SPEED,pch=16,type="o",ylab = 'Speed(um/h)',xlab = 'Cells',las=1,ylim = c(0, MS))
  graphics::title(main="Speed Profile of all cells",cex.main = 0.7)
  graphics::abline(h=stats::median(SPEED[which(!is.na(SPEED))]),col="red")
  graphics::abline(h=mean(SPEED[which(!is.na(SPEED))]),col="blue")
  graphics::legend(1, y=200, legend=c("Mean Speed","Median Speed"), col=c("blue","red"),lty=1, cex=0.8)
  if (export) grDevices::dev.off()

  if (export) grDevices::jpeg(paste0(ExpName,"_Speed_ViolinPlot_of_all_cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  graphics::plot(1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '', ylab = 'Speed(um/h)', xaxt = 'n',las=1)
  graphics::title("Speed of all cells",cex.main = 1)
  vioplot::vioplot(SPEED, at = 1, add = TRUE, col = "gray")
  if (export) grDevices::dev.off()

  if (export) grDevices::jpeg(paste0(ExpName,"_Instantaneous_Speed_VS_Persistence Ratio_all_cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  SPEED<-as.numeric(PerResultsTable[21,1:length(PerResultsTable[1,])-1])
  PerR<- as.numeric(PerResultsTable[4,1:length(PerResultsTable[1,])-1])
  MS<-max(SPEED)
  graphics::plot(SPEED,PerR,pch=16,type="p",xlab="Average Instantaneous Speed (um/h)",ylab=" Persistence Ratio",col="black",las=1,xlim=c(0,MS),ylim=c(0,1))
  reg<-stats::lm(PerR~SPEED)
  graphics::abline(reg,untf=FALSE,col="red")
  RowmeanSpeed<-round(matrixStats::rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-cor.test( ~ SPEED+ PerR, method = "spearman",exact=FALSE)                 # testing the correlation
  ss<-unlist(s[4])
  SCC<-round(ss, digits = 3)
  graphics::title(main=paste0("All Cells Average Speed vs Persistence Ratio "),cex.main =0.7,sub=paste0("Spearman's rank correlation coefficient = ",SCC),col.sub="red")
  if (export) grDevices::dev.off()

  PerResultsTable[1,(length(Object)+1)]<-"All Cells"
  object@PerAanSpeedtable <-PerResultsTable
  setwd(d)
  if (export) {
    utils::write.csv(
      PerResultsTable,
      file = paste0(ExpName,"-PerResultsTable.csv")
    )
    cat(
      "Results are saved as: ",
      paste0(ExpName,"-PerResultsTable.csv" ),
      "in your directory [use getwd()]\n"
    )
  }
  return(object)
}



#' @title Directionality Table
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straight line length between the start point and the endpoint of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return An CellMig class object with a data frame stored in the DRtable slot. It contains nine rows: "Cell Number","Directionality Ratio","Mean Cumulative Directionality Ratio","Stable Directionality Ratio", "Number of returns","Min CumDR","Location of Min CumDR, Steps with less CumDR than DR","Directional Persistence".
#'
#' @details  Directionality Ratio and Directional persistence
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- DiRatio(rmTD, ExpName="Test", export=FALSE)
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
DiRatio = function(object,TimeInterval=10,ExpName="ExpName", export=FALSE) {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  d=getwd()
  if (export) {
    dir.create(paste0(ExpName,"-DRResults"))
    setwd(paste0(d,"/",paste0(ExpName,"-DRResults")))
  }
  Step<-length(Object[[1]][,1])
  DRResultsTable<-data.frame()
  DIR.RATIO<-c()
  for(j in 1:length(Object)){                             # calculating the cumsum of distance for each cell
    MM<-Step
    MM2<-MM-1
    end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])       # finding the cordinates of the final point in the track.
    start<-cbind(0,0)
    final.dis=SpatialTools::dist2(start, end)                               # calculating the final distance
    total.dis=sum(Object[[j]][1:MM2,6])                     # calculating the total distance
    Dir.Ratio=round(final.dis/total.dis ,digits = 3)          # calculating the Dir.Ratio
    mean.Dir.Ratio<-round(mean(Object[[j]][1:MM2,13],na.rm = TRUE) ,digits = 3)
    StableDR<- (1-(mean.Dir.Ratio-Dir.Ratio))* mean.Dir.Ratio
    StableDR<-round(StableDR,digits = 3)

    DRResultsTable[1,j]<- j
    DRResultsTable[2,j]<-Dir.Ratio
    DRResultsTable[3,j]<-mean.Dir.Ratio
    DRResultsTable[4,j]<- StableDR

  }


  for(j in 1:length(Object)){                                   #### Adding min CumDR  and number of angles greater than .75
    MM<-Step
    MM2<-MM-1
    p1<-Object[[j]][1:MM2,9]
    returns<-subset(p1,p1<(-0.87))                            # greater than 150 degrees
    DRResultsTable[5,j]<-length(returns)

    p2<-Object[[j]][1:MM2,13]
    w<-which.min(p2)
    DRResultsTable[6,j]<-round(p2[w], digits=3)
    DRResultsTable[7,j]<-w
    DR<-as.numeric(DRResultsTable[2,j])
    lessThanMINcumdr<-subset(p2,p2<DR)                        # number of the steps that have a cumdr less than the final dr
    DRResultsTable[8,j]<-length(lessThanMINcumdr)

    Ptime<-Object[[j]][1:MM2,10]
    PerTimLen<-Ptime                                     # computing the number of persistence steps to be used in computing the persistence ratio
    PerTimLen[PerTimLen==0]<-NA
    PerTimLen<-PerTimLen[!is.na(PerTimLen)]
    PerTimLen<-length(PerTimLen)
    PerRatio<-round(PerTimLen/MM2, digits=2)
    DRResultsTable[9,j]<- PerRatio + DRResultsTable[4,j]

  }
  rownames(DRResultsTable)<-c("Cell Number","Directionality Ratio","Mean Cumulative Directionality Ratio","Stable Directionality Ratio",
                              "Number of returns","Min CumDR",paste0("Location of Min CumDR (out of ",Step-1,")"),"Steps with less CumDR than DR","Directional Persistence")
  RM1<-round(matrixStats::rowMedians(as.matrix(DRResultsTable),na.rm = TRUE),digits=3)
  DRResultsTable[,(length(Object)+1)]<-RM1
  DRResultsTable[1,(length(Object)+1)]<-"All Cells"
  object@DRtable<-DRResultsTable
  setwd(d)

  if (export) {
    utils::write.csv(
      DRResultsTable,
      file = paste0(ExpName,"-DRResultsTable.csv")
    )
    cat(
      "Results are saved as: ",
      paste0(ExpName,"-DRResultsTable.xlsx" ),
      "in your directory [use getwd()]",
      "\n"
    )
  }
  return(object)
}


#' @title Directionality Ratio plots
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straightline length between the start point and the endpoint of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param export if `TRUE` (default), exports plot to JPG file
#'
#' @return Directionality Ratio plots
#'
#' @details  Directionality Ratio
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom grDevices rainbow jpeg dev.off rgb
#' @importFrom graphics plot axis title lines polygon
#' @importFrom matrixStats rowMedians rowSds
#'
#' @export
DiRatio.Plot = function(object,TimeInterval=10,ExpName=ExpName, export=FALSE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  Step<-length(Object[[1]][,1])
  color <-c()
  Len<-length(Object)
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-grDevices::rainbow(Len)
  }
  d=getwd()

  if (file.exists(paste0(d,"/",paste0(ExpName,"-DRResults")))){
    setwd(paste0(d,"/",paste0(ExpName,"-DRResults")))

  } else {
    stop("Please run DiRatio() first and export the results")
  }

  DIR.RATIO.AllCells<-data.frame()
  for(j in 1:length(Object)){
    MM<-Step
    MM2<-MM-1
    Time<-c(1:MM2)
    t<-60/TimeInterval
    xM<-MM2*TimeInterval
    xMM<-round(xM/60)
    xMMM<-c(0:xMM)
    xMMM1<-(c(0:xMM)*(60/TimeInterval))
    Xaxis<-c(1:MM2)
    if (export) grDevices::jpeg(paste0(ExpName,"-D.R.plot-",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
    p<-graphics::plot(Time,Object[[j]][1:MM2,13], type="l",col=color[j],xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1)
    graphics::axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
    graphics::title(main=paste0("Cell Number  ", j),col.main=color[j])
    if (export) grDevices::dev.off()
    DIR.RATIO.AllCells[1:MM2,j]<-Object[[j]][1:MM2,13]
  }
  mycolblue <- grDevices::rgb(0, 0, 255, maxColorValue = 255, alpha = 100, names = "blue")    #transparent color
  mean.DIR.RATIO.AllCells<-matrixStats::rowMedians(as.matrix(DIR.RATIO.AllCells),na.rm = TRUE)
  SD.DIR.RATIO.AllCells<-matrixStats::rowSds(as.matrix(DIR.RATIO.AllCells),na.rm = TRUE)    # is a function in matrixStats
  meanSDp<-mean.DIR.RATIO.AllCells+SD.DIR.RATIO.AllCells
  meanSDp=ifelse(meanSDp>=1,1,meanSDp)
  meanSDn<-mean.DIR.RATIO.AllCells-SD.DIR.RATIO.AllCells
  meanSDpn<-c(meanSDp,meanSDn)

  if (export) {
    grDevices::jpeg(paste0(ExpName,"directionality ratio for all cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
  }
  p<-graphics::plot(Time,mean.DIR.RATIO.AllCells, type="l",col="black",xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1)
  graphics::axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
  graphics::lines(Time,meanSDp, col="black")
  graphics::lines(Time,meanSDn, col="black")
  graphics::polygon(c(Time, rev(Time)), c(meanSDp, rev(meanSDn)),col = mycolblue , border = NA)
  graphics::title(main="Directionality Ratio - All Cells",col.main="black")
  if (export) grDevices::dev.off()
  setwd(d)
  cat("Plots are saved in a folder in your directory [use getwd()]","\n")
}



#' @title Mean Square Displacement
#' @description The MSD function automatically computes the mean square displacements across several sequential time intervals. MSD parameters are used to assess the area explored by cells over time.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param ffLAG A numeric value to be used to get the number of lags for the  Furth formula fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param SlopePlot A logical vector that allows generating individual plots showing the slope of the mean square displacement of the movement of individual cells. Default is TRUE.
#' @param AllSlopesPlot A logical vector that allows generating a plot showing the slope of the mean square displacement of the movement of all cells. Default is TRUE.
#' @param FurthPlot A logical vector that allows generating individual plots fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method per cell. Default is TRUE.
#' @param AllFurthPlot A logical vector that allows generating a plot fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method for all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#'
#' @return An CellMig class object with a data frame and plots. The data frame is stored in the MSDtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF <- TrajectoryDataset[1:600, ]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- MSD(rmTD, sLAG=0.25, ffLAG=0.25, export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom stats lm coef
#' @importFrom graphics par plot title abline lines
#' @importFrom matrixStats rowMedians rowSds
#' @importFrom FME modCost modFit
#' @importFrom utils write.csv
#'
#' @export
MSD <- function(object, TimeInterval=10,
               ExpName="ExpName",
               ExpDir=tempdir(),
               sLAG=0.25, ffLAG=0.25,
               SlopePlot=TRUE, AllSlopesPlot=TRUE,
               FurthPlot=TRUE, AllFurthPlot=TRUE, export=FALSE) {
  # ============================================================================
  # Validation
  # ============================================================================
  if (!is.numeric(TimeInterval) | TimeInterval<= 0) {
    stop( "TimeInterval has to be a positive number")
  }
  if (!is.numeric(sLAG) | sLAG<= 0) {
    stop( "sLAG has to be a positive number")
  }
  if (!is.numeric(ffLAG) | ffLAG<= 0) {
    stop( "ffLAG has to be a positive number")
  }
  # ============================================================================
  # Operations
  # ============================================================================
  Object <- object@preprocessedDS
  if (export) {
    SavePath <- paste0(ExpDir, "/", ExpName,"-MSDResults")
    if (!dir.exists(SavePath)) dir.create(SavePath)
  }
  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }

  MSDResultsTable<-data.frame()
  MSD.table<-data.frame() # creating a table that has all the MSDs to be able to compute the mean and sd
  for(j in 1:length(Object)){
    meanSD<-c()
    LAG<-round(Step*sLAG) # number of lags is based on sLAG
    for(lag in 1:LAG){
      res <- sapply(1:Step, function(i){
        Object[[j]][i,22]=(((Object[[j]][i+lag,2]- Object[[j]][i,2])^2)+((Object[[j]][i+lag,3]- Object[[j]][i,3])^2))
        return(Object[[j]][i,22])
      })
      Object[[j]][1:Step, 22] <- as.data.frame(res)
      Object[[j]][,22][is.na(Object[[j]][,22])] <- 0
      meanSD[lag]<-mean(Object[[j]][1:(Step-LAG),22])
      Object[[j]][,22]=0
    }
    MSDResultsTable[1,j]<-j
    MSDResultsTable[2,j]<-round(meanSD[1],digits=3)
    MSD.table[1:LAG,j]<-meanSD
    Yaxis<-assign(paste0("MSD.Cell.",j),meanSD)
    Xaxis<-c(1:LAG)
    NewrowMeans<-meanSD[1:(round(LAG*sLAG))]                                    # plotting the regression line based on sLAG
    NewXaxis<-Xaxis[1:(round(LAG*sLAG))]
    reg<-stats::lm(NewrowMeans~ NewXaxis)
    reg1<-round(stats::coef(stats::lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
    MSDResultsTable[3,j]<-reg1
    # --------------------------------------------------------------------------
    # Generating separate slope plots
    # --------------------------------------------------------------------------
    if (SlopePlot) {
      if (export) {
        grDevices::jpeg(
          paste0(SavePath, "/", ExpName, "-MSD.plot.Cell", j, ".jpg"),
          width = 4, height = 4, units = 'in', res = 300
        )
      }
      xn <- expression(paste("MSD (um"^2, ")"))
      graphics::par(mar=c(5.1, 5.9, 4.1, 1.1), mgp=c(3.5, 0.5, 0), las=0)
      graphics::plot(Xaxis,Yaxis, type="p",col=color[j],xlab="Lag",ylab=xn,pch=19,las=1,log="xy",cex=2)
      graphics::title(main=paste0("Cell Number  ", j," -  MSD Slope = ",reg1),col.main="black")
      graphics::abline(reg,untf=TRUE,col="black")
      if (export) grDevices::dev.off()
    }
  }
  # ============================================================================
  # Generating combined slopes plot
  # ============================================================================
  if (AllSlopesPlot) {
    RM1<-matrixStats::rowMedians(as.matrix(MSD.table),na.rm = TRUE)
    RSD<-matrixStats::rowSds(as.matrix(MSD.table),na.rm = TRUE)
    xn <- expression(paste("MSD (um"^2, ")"))
    LAG<-round(Step*sLAG)
    Xaxis<-c(1:LAG)
    if (export) {
      jpeg(
        paste0(SavePath, "/", ExpName, "-MSD.plot All Cells.jpg"),
        width = 4, height = 4, units = 'in', res = 300
      )
    }
    xn <- expression(paste("MSD (um"^2, ")"))
    graphics::par(mar=c(5.1, 5.5, 4.1, 0.9), mgp=c(3.5, .5, 0), las=0)
    graphics::plot(Xaxis,RM1, type="p",col="black",xlab="Lag",ylab=xn,pch=19,las=1,cex=1.5,log="xy")
    NewrowMeans<-RM1[1:(round(LAG*sLAG)+1)]                                  # best fit based on sLAG *sLAG
    NewXaxis<-Xaxis[1:(round(LAG*sLAG)+1)]
    reg<-stats::lm(NewrowMeans~ NewXaxis)
    reg1<-round(stats::coef(stats::lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
    graphics::abline(reg,untf=TRUE,col="red")
    graphics::title(main=paste0("All Cells -  MSD Slope = ",reg1),col.main="black")
    if (export) grDevices::dev.off()
  }
  for (j in 1: length(MSD.table[1,])){                      # Fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method
    LAG<-round(Step*ffLAG)
    y=MSD.table[1:LAG,j]
    t <- c(1:length(y))
    Data<-c(t,y)
    Data <- matrix(ncol =  2, byrow = FALSE, data =Data)
    colnames(Data) <- c("time","y")
    parms <- c(D=1,P=1)
    parms["D"] <- 1
    parms["P"] <- 1
    Cost <- function(Par) {
      D <- Par[1]
      P <- Par[2]
      out <- cbind(time = t, y = (D*4*(t-P*(1-(exp(-t/P))))))
      return(FME::modCost(obs = Data, model = out))
    }
    Fit<-FME::modFit(p = c(D = 1, P = 1), lower = c(0,0),upper=c(100,100),f = Cost,method = "Nelder-Mead")
    Summary<-summary(Fit)
    MSDResultsTable[4,j]<-round(Fit$par[1],digits=3)
    MSDResultsTable[5,j]<-round(Fit$par[2],digits=3)
    MSDResultsTable[6,j]<-round(Summary$par[1,4],digits=3)
    MSDResultsTable[7,j]<-round(Summary$par[2,4],digits=3)
    # --------------------------------------------------------------------------
    # Generating individual Furth plots
    # --------------------------------------------------------------------------
    if (FurthPlot){
      if (export) {
        grDevices::jpeg(
          paste0(SavePath, "/", ExpName, "-MSD N-M bestfit Cell", j, ".jpg"),
          width = 4, height = 4, units = 'in', res = 300
        )
      }
      graphics::plot(Data, pch = 16,col=color[j], cex = 1.5, xlab = "Lags", ylab = "MSD")
      x<-seq(0,LAG,1)
      Model <- function(p, x) return(data.frame(x = x, y = p[1]*4*(x- p[2]*(1-(exp(-x/p[2]))))))
      graphics::lines(Model(Fit$par, x),col="black")
      graphics::title(main=paste0("Cell Number  ", j,"    Nelder-Mead best-fit"),col.main="black",
                      sub=paste0("D = ",round(Fit$par[1],digits=3),
                                 "          P = ", round(Fit$par[2],digits=3)),col.sub="red")
      if (export) grDevices::dev.off()
    }
  }

  RM1<-matrixStats::rowMedians(as.matrix(MSD.table),na.rm = TRUE)
  MSDResultsTable[1,(length(Object)+1)]<-"All Cells"
  MSDResultsTable[2,(length(Object)+1)]<-round(RM1[1],digits=3)

  RM<-c(0,RM1)
  Xaxis<-c(1:round(Step*ffLAG))
  NewrowMeans<-RM1[1:(round(LAG*ffLAG)+1)]                                  # best fit based on ffLAG *ffLAG
  NewXaxis<-Xaxis[1:(round(LAG*ffLAG)+1)]

  reg<-stats::lm(NewrowMeans~ NewXaxis)
  reg1<-round(stats::coef(stats::lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
  MSDResultsTable[3,(length(Object)+1)]<-reg1
  y=RM1[1:round(Step*ffLAG)]
  t <- c(1:length(y))
  Data<-c(t,y)
  Data <- matrix(ncol =  2, byrow = FALSE, data =Data)
  colnames(Data) <- c("time","y")
  parms <- c(D=10,P=1)
  parms["D"] <- 1
  parms["P"] <- 1
  Cost <- function(Par) {
    D <- Par[1]
    P <- Par[2]
    out <- cbind(time = t, y = (D*4*(t-P*(1-(exp(-t/P))))))
    return(FME::modCost(obs = Data, model = out))
  }
  Fit<-FME::modFit(p = c(D = 1, P = 1), lower = c(0, 0),upper=c(100,100),f = Cost,method = "Nelder-Mead")
  MSDResultsTable[4,(length(Object)+1)]<-round(Fit$par[1],digits=3)
  MSDResultsTable[5,(length(Object)+1)]<-round(Fit$par[2],digits=3)
  MSDResultsTable[6,(length(Object)+1)]<-round(Summary$par[1,4],digits=3)
  MSDResultsTable[7,(length(Object)+1)]<-round(Summary$par[2,4],digits=3)
  # ============================================================================
  # Generating consolidated Furth plot
  # ============================================================================
  if (AllFurthPlot) {
    if (export) {
      grDevices::jpeg(
        paste0(SavePath, "/", ExpName, "-MSD N-M bestfit All Cells.jpg"),
        width = 4, height = 4, units = 'in', res = 300
      )
    }
    graphics::plot(Data, pch = 16,col="black", cex = 1.5, xlab = "Lags", ylab = "MSD")
    x<-seq(0,LAG,1)
    Model <- function(p, x) return(data.frame(x = x, y = p[1]*4*(x- p[2]*(1-(exp(-x/p[2]))))))
    graphics::lines(Model(Fit$par, x),col="red")
    graphics::title(main=paste0("All Cells Nelder-Mead best-fit"),col.main="black",
          sub=paste0("D = ",round(Fit$par[1],digits=3),
                     "          P = ", round(Fit$par[2],digits=3)),col.sub="red")
    if (export) grDevices::dev.off()
  }

  rownames(MSDResultsTable)<-c("Cell Number","MSD (lag=1)", "MSD slope", "N-M best fit (Furth) [D]","N-M best fit (Furth) [P]","The significance of fitting D","The significance of fitting P")
  object@MSDtable<-MSDResultsTable
  if (export) {
    utils::write.csv(
      MSDResultsTable,
      file = paste0(SavePath, "/", ExpName,"-MSDResultsTable.csv")
    )
    message("Results were saved to ", ExpDir)
  }
  return(object)
}



#' @title Direction AutoCorrelation
#'
#' @description The DiAutoCor function automatically compute the angular persistence across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param sPLOT A logical vector that allows generating individual plots showing the angular persistence across several sequantial time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot showing the angular persistence across several sequantial time intervals of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return An CellMig class Object with a data frame and plots. The data frame, which contains six rows: "Cell Number", "Angular Persistence", "Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)", "Stable Direction AutoCorrelation through the track" and "Difference between Mean DA and Intercept DA".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#'
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- DiAutoCor(
#'    rmTD, TimeInterval=10, ExpName="ExpName", sLAG=0.25, sPLOT=FALSE,
#'    aPLOT=FALSE, export=FALSE
#' )
#'
#' @importFrom grDevices rainbow
#' @importFrom stats lm predict median
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics plot lines title abline legend
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
DiAutoCor= function(object, TimeInterval=10,
                    ExpName="ExpName", sLAG=0.25,
                    sPLOT=TRUE, aPLOT=TRUE, export=FALSE) {

  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  Object<-object@preprocessedDS

  if ( ! is.list(Object) ){
    stop("Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  if (export) {
    dir.create(paste0(ExpName,"-DIAutoCorResults"))
    setwd(paste0(d,"/",paste0(ExpName,"-DIAutoCorResults")))
  }

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }

  DA.ResultsTable<-data.frame()
  DA.table<-data.frame()

  for(j in 1:length(Object)){
    cos.diff<-c()
    LAG<-round(Step*sLAG)     #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:Step, function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2]    # newdx
        Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3]    # newdy
        return(Object[[j]][i,14:15])
      }))
      Object[[j]][1:Step,14:15] <- as.data.frame(res)


      Object[[j]][,14:15] <- lapply(Object[[j]][,14:15], as.numeric)
      Object[[j]][,14][is.na(Object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      Object[[j]][,15][is.na(Object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:Step, function(i){
        Object[[j]][i,16]=acos((Object[[j]][i+lag,2]-Object[[j]][i,2])/sqrt((Object[[j]][i+lag,2]-Object[[j]][i,2])^2 +(Object[[j]][i+lag,3]-Object[[j]][i,3])^2)) # to find the abs angle
        return(Object[[j]][i,16])
      })
      Object[[j]][1:(Step),16] <- res1
      Object[[j]][,16][is.na(Object[[j]][,16])] <- 0                                # to remove NA and replace it with 0


      res2 <- sapply(1:(Step - lag), function(i){
        if((Object[[j]][i+1,15]<0) && (Object[[j]][i,15]>=0)||(Object[[j]][i+1,15]>=0) && (Object[[j]][i,15]<0)){
          Object[[j]][i,17]= abs(Object[[j]][i+1,16])+abs(Object[[j]][i,16])
        }
        if((Object[[j]][i+1,15]<0) && (Object[[j]][i,15]<0)||(Object[[j]][i+1,15]>=0) && (Object[[j]][i,15]>=0) ){
          Object[[j]][i,17]=Object[[j]][i+1,16]-Object[[j]][i,16]
        }

        return(Object[[j]][i,17])
      })
      Object[[j]][1:(Step - lag),17] <- res2

      res3 <- t(sapply(1:(Step - lag), function(i){
        Object[[j]][i,17]<-ifelse((Object[[j]][i,17])<= (-pi), 2*pi+(Object[[j]][i,17]),(Object[[j]][i,17]))    # adjusting the ang.diff
        Object[[j]][i,17]<-ifelse((Object[[j]][i,17])>= pi,(Object[[j]][i,17])-2*pi,(Object[[j]][i,17]))
        Object[[j]][i,18]<-cos(Object[[j]][i,17])
        return(Object[[j]][i,17:18])
      }))
      Object[[j]][1:(Step - lag),17:18] <- as.data.frame(res3)
      Object[[j]][,17:18] <- lapply(Object[[j]][,17:18], as.numeric)

      for(i in 1:LAG){
        cos.diff[lag]<-mean(Object[[j]][1:(Step-lag)-1,18], na.rm=TRUE)  # computing the cosine mean
      }
    }
    DA.table[1:LAG,j]<-cos.diff
    assign(paste0("DA.Cell.",j),cos.diff)
    DA.ResultsTable[1,j]<-j
    DA.ResultsTable[2,j]<-round(cos.diff[1],digits=3)
    lags<-c(1:length(DA.table[,1]))
    lags2<- lags^2
    quadratic.model<-c()
    quadratic.m <-stats::lm(DA.table[,j]~ lags + lags2)
    c<-quadratic.m
    cc<-unlist(c)
    quadratic.model[j]<-cc[1]
    DA.ResultsTable[3,j]<-round(unlist(cc[1]),digits=3)
    ccc<-unlist(cc[1])
    DA.ResultsTable[4,j]<-round(mean(DA.table[1:LAG,j]),digits=3)
    DA.ResultsTable[5,j]<-round((1-(mean(DA.table[1:LAG,j])-unlist(cc[1])))* mean(DA.table[1:LAG,j]),digits=3)
    DA.ResultsTable[6,j]<-round(mean(DA.table[1:LAG,j])-unlist(cc[1]),digits=3)

    timevalues <- seq(1, length(lags), 1)
    predictedcounts <- stats::predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

    if ( sPLOT == TRUE){
      Xaxis<-c(1:LAG)
      Yaxis<-cos.diff
      if (export) {
        grDevices::jpeg(
          paste0(ExpName,"Direction Autocorrelation.plot.Cell",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      }
      graphics::plot(Xaxis,Yaxis, type="o",ylim=c(-1,1),xlim=c(0,lag),col=color[j],xlab="Lag",ylab="Cosine",pch=19,las=1,cex=1.2)
      xx<-c(0,1)
      yy<-c(1,cos.diff[1])
      graphics::lines(xx,yy, type='l',col=color[j])
      graphics::lines(timevalues, predictedcounts, col = "black", lwd = 3)
      graphics::title(main=paste0("Cell Number  ", j, "   DA quadratic model"),col.main="darkgreen",
            sub=paste0(" Intercept of DA quadratic model = ",round(ccc, digits=3)),col.sub="red")
      if (export) grDevices::dev.off()
    }

    Object[[j]][,15:18]=0
  }
  RM1<-matrixStats::rowMedians(as.matrix(DA.table),na.rm = TRUE)
  DA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
  DA.ResultsTable[2,(length(Object)+1)]<-round(RM1[1],digits=3)
  lags<-c(1:length(DA.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-stats::lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  DA.ResultsTable[3,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
  DA.ResultsTable[4,(length(Object)+1)]<-round(stats::median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=3)
  DA.ResultsTable[5,(length(Object)+1)]<-round((1-(stats::median(as.numeric(DA.ResultsTable[4,1:length(Object)]))-unlist(cc[1])))* stats::median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=3)
  DA.ResultsTable[6,(length(Object)+1)]<-round(median(as.numeric(DA.ResultsTable[4,1:length(Object)]))-unlist(cc[1]),digits=3)

  ccc<-unlist(cc[1])
  timevalues <- seq(1, length(lags), 1)
  predictedcounts <- stats::predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

  if ( aPLOT == TRUE){
    Xaxis<-c(1:LAG)
    Yaxis<-RM1

    if (export) {
      grDevices::jpeg(
        paste0(ExpName,"-Direction Autocorrelation All Cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
    }
    graphics::plot(Xaxis,Yaxis, type="o",ylim=c(-1,1),xlim=c(0,lag),col="black",xlab="Lag",ylab="Cosine",pch=19,las=1,cex=1.2)
    xx<-c(0,1)
    yy<-c(1,RM1[1])
    graphics::lines(xx,yy, type='l',col="black")
    graphics::lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
    MDA<-round(stats::median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=2)
    graphics::abline(h=MDA,col="blue",lwd = 2)
    graphics::title(main=paste0("All Cells - DA quadratic model"),col.main="darkgreen",
          sub=paste0(" Intercept of DA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    graphics::legend(1, y=-0.82, legend=c("Mean Direction AutoCorrelation","Quadratic model"), col=c("blue","darkgreen"),lty=1, cex=0.8)
    if (export) grDevices::dev.off()
  }

  rownames(DA.ResultsTable)<-c("Cell Number","Angular Persistence","Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)","Stable Direction AutoCorrelation through the track",
                               "Difference between Mean DA and Intercept DA" )
  object@DACtable<-DA.ResultsTable
  setwd(d)
  if (export) {
    utils::write.csv(
      DA.ResultsTable,
      file = paste0(ExpName,"-DA.ResultsTable.csv")
    )
    cat("Results are saved in your directory [use getwd()]","\n")
  }
  return(object)
}


#' @title Velocity AutoCorrelation
#'
#' @description The VeAutoCor function automatically compute the changes in both speed and direction across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param sPLOT A logical vector that allows generating individual plots showing the velocity across several sequantial time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot showing the velocity across several sequantial time intervals of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return Plots and a data frame, which contains six rows: "Cell Number", .
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:1000,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- VeAutoCor(
#'    rmTD, TimeInterval=10, ExpName="ExpName", sLAG=0.25, sPLOT=FALSE,
#'    aPLOT=FALSE, export=FALSE
#' )
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom stats lm predict median
#' @importFrom graphics plot lines title abline
#' @importFrom utils write.csv
#'
#' @export
VeAutoCor = function(object, TimeInterval=10,
                     ExpName="ExpName", sLAG=0.25,
                     sPLOT=TRUE, aPLOT=TRUE, export=FALSE) {

  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  if (export) {
    dir.create(paste0(ExpName,"-VeAutoCorResults"))
    setwd(paste0(d,"/",paste0(ExpName,"-VeAutoCorResults")))
  }

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }

  VA.ResultsTable<-data.frame()
  VAC.table<-data.frame()     # creating a table that has all the VAC to be able to compute the mean
  VAC.first.value<-c()        # to save the fist VAC non-normalized value for all cells
  VAC.second.value<-c()        # to save the second VAC non-normalized value for all cells

  for(j in 1:length(Object)){
    meanVAC<-c()
    LAG<-round(Step*sLAG)              #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:(Step - 1), function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2]    # newdx
        Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3]    # newdy
        return(Object[[j]][i,14:15])
      }))
      Object[[j]][1:(Step -1),14:15] <- as.data.frame(res)
      Object[[j]][,14:15] <- lapply(Object[[j]][,14:15], as.numeric)
      Object[[j]][,14][is.na(Object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      Object[[j]][,15][is.na(Object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:(Step - lag), function(i){                                               # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,23]=((Object[[j]][i,14]* Object[[j]][i+lag,14])+ (Object[[j]][i,15]* Object[[j]][i+lag,15]))/ ((lag*TimeInterval)^2)
        return(Object[[j]][i,23])
      })
      Object[[j]][1:(Step - lag),23] <- res1
      meanVAC[lag]<-mean(Object[[j]][1:(Step - lag),23])

    }
    VAC.first.value[j]<-meanVAC[1]

    NORMmeanVAC<-meanVAC/meanVAC[1]
    VAC.table[1:LAG,j]<-meanVAC
    assign(paste0("VAC.Cell.",j),meanVAC)
    VAC.second.value[j]<-NORMmeanVAC[2]

    VA.ResultsTable[1,j]<-j
    VA.ResultsTable[2,j]<-round(meanVAC[1],digits=3)        # VA (lag =1)
    VA.ResultsTable[3,j]<-round(NORMmeanVAC[2],digits=3)    # VA (lag =2)

    lags<-c(1:length(VAC.table[,1]))
    lags2<- lags^2
    quadratic.model<-c()
    quadratic.m <-stats::lm(VAC.table[,j]~ lags + lags2)
    c<-quadratic.m
    cc<-unlist(c)
    quadratic.model[j]<-cc[1]

    ccc<-unlist(cc[1])
    VA.ResultsTable[4,j]<-round(ccc,digits=3)
    VA.ResultsTable[5,j]<-round(mean(meanVAC),digits=3)      # mean VA (all lags)
    timevalues <- seq(1, length(lags), 1)
    predictedcounts <- stats::predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

    if (sPLOT == TRUE){
      Xaxis<-c(1:LAG)
      Yaxis<-meanVAC
      if (export) {
        grDevices::jpeg(
          paste0(ExpName,"Velocity Autocorrelation.plot.Cell",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      }
      graphics::plot(Xaxis,Yaxis, type="o",ylim=range(meanVAC),xlim=c(0,lag),col=color[j],xlab="Lag",ylab="Velocity  Autocorrelation",pch=19,las=1,cex=1.2)
      graphics::lines(timevalues, predictedcounts, col = "black", lwd = 3)
      graphics::title(main=paste0("Cell Number  ", j, "   VA quadratic model"),col.main="darkgreen",
            sub=paste0(" Intercept of VA quadratic model = ",round(ccc, digits=3)),col.sub="red")
      if (export) grDevices::dev.off()
    }

    Object[[j]][,15:16]=0
  }

  RM1<-matrixStats::rowMedians(as.matrix(VAC.table),na.rm = TRUE)
  VA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
  VA.ResultsTable[2,(length(Object)+1)]<-round(median(VAC.first.value),digits=3)
  VA.ResultsTable[3,(length(Object)+1)]<-round(median(VAC.second.value),digits=3)

  lags<-c(1:length(VAC.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  VA.ResultsTable[4,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
  VA.ResultsTable[5,(length(Object)+1)]<-round(stats::median(as.numeric(VA.ResultsTable[5,1:length(Object)])),digits=3)

  ccc<-unlist(cc[1])
  timevalues <- seq(1, length(lags), 1)
  predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

  if ( aPLOT == TRUE){
    Xaxis<-c(1:LAG)
    Yaxis<-RM1
    if (export) {
      grDevices::jpeg(paste0(ExpName,"-Velocity Autocorrelation All Cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
    }
    graphics::plot(Xaxis,Yaxis, type="o",ylim=range(RM1),xlim=c(0,lag),col="black",xlab="Lag",ylab="Velocity  Autocorrelation",pch=19,las=1,cex=1.2)
    graphics::lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
    MVA<-round(stats::median(as.numeric(VA.ResultsTable[5,1:length(Object)])),digits=2)
    graphics::abline(h=MVA,col="blue",lwd = 2)
    graphics::title(main=paste0("All Cells - VA quadratic model"),col.main="darkgreen",
          sub=paste0(" Intercept of VA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    if (export) grDevices::dev.off()
  }

  rownames(VA.ResultsTable)<-c("Cell Number","Velocity AutoCorrelation (lag=1)","2nd normalized Velocity AutoCorrelation",
			       "Intercept of VA quadratic model","Mean Velocity AutoCorrelation (all lags)")
  object@VACtable<-VA.ResultsTable
  setwd(d)
  if (export) {
    utils::write.csv(
      VA.ResultsTable,
      file = paste0(ExpName,"-VA.ResultsTable.csv")
    )
    cat("Results are saved in your directory [use getwd()]\n")
  }
  return(object)
}


#' @title Forward Migration
#'
#' @description The ForwardMigration function automatically generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param sfptPLOT A logical vector that allows generating individual plots of persistence time vs speed per cell. Default is TRUE.
#' @param afptPLOT  A logical vector that allows generating a plot of persistence time vs speed for all cells. Default is TRUE.
#' @param sfpPLOT A logical vector that allows generating individual plots of angular persistence vs speed per cell. Default is TRUE.
#' @param afpPLOT A logical vector that allows generating a plot of angular persistence vs speed of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return  An CellMig class Object with a data frame and plots. The data frame is stored in the ForMigtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF=WSADataset[1:1000,]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' wsaTD <-ForwardMigration(
#'    wsaTD, TimeInterval=10, ExpName="ExpName", export=FALSE
#' )
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom stats cor.test lm
#' @importFrom graphics plot title abline
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
ForwardMigration <- function(
    object,
    TimeInterval = 10,
    ExpName      = "ExpName",
    sfptPLOT     = TRUE,
    afptPLOT     = TRUE,
    sfpPLOT      = TRUE,
    afpPLOT      = TRUE,
    export       = FALSE
){
    if (!is.numeric(TimeInterval)) {
        stop("TimeInterval has to be a positive number")
    } else if (TimeInterval <= 0) {
        stop("TimeInterval has to be a positive number")
    }
    Object <- object@preprocessedDS
    UPorDO <- object@cellpos
    msg    <- NULL
    if (!is.list(Object)) {
        msg <- c(
            msg, "Input data must be a list.",
            "Please run the PreProcessing step first either ",
            "rmPreProcessing() or wsaPreProcessing()"
        )
    }
    d     <- getwd()
    Len   <- length(Object)
    Step  <- length(Object[[1]][, 1])
    color <- c()
    if (Len > 1023) {
        colnum <- Len-1023
        color1 <- grDevices::rainbow(1023)
        colo2  <- grDevices::rainbow(colnum)
        color  <- c(color1, colo2)
    } else {
        color <- grDevices::rainbow(Len)
    }
    if (export) {
      dir.create(paste0(ExpName, "-ForwardMigrationResults"))
      setwd(paste0(d, "/", paste0(ExpName,"-ForwardMigrationResults")))
    }

    #if ( length(UPorDO) != length(Object))
    # stop("UPorDO needs to be a vector with same number of elements as cells")

    # Defining if the cell is ubove (1) the wound or below (0) it
    for (j in 1:length(Object)) {
        Object[[j]][, 25] <- UPorDO[j]
    }

    # creating values for  rel.ang.F  (step to the original)
    for(j in 1:length(Object)){
        MM  <- Step
        MM1 <- MM - 1
        res <- sapply(1:MM1, function(i) {

        if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){
            Object[[j]][i,19]= 1.5707963268 - abs(Object[[j]][i,7])
        }
        if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
            Object[[j]][i,19]= abs(Object[[j]][i,7])+1.5707963268

        }
        Object[[j]][i,19]<-ifelse((Object[[j]][i,19])<= (-pi), 2*pi+(Object[[j]][i,19]),(Object[[j]][i,19]))    # adjusting the rel.ang
        Object[[j]][i,19]<-ifelse((Object[[j]][i,19])>= pi,(Object[[j]][i,19])-2*pi,(Object[[j]][i,19]))
        return(Object[[j]][i, 19])
        })
        Object[[j]][1:MM1, 19] <- as.data.frame(res)
    }


  cosine.FP<-data.frame()
  for(j in 1:length(Object)){              # creating values for  cosine based on rel.ang.F
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM, function(i){
      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){   # upper cell going up or lower cell going down
        Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
      }

      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
        Object[[j]][i,20]<-cos(Object[[j]][i,19])
      }
      return(Object[[j]][i,20])
    })

    Object[[j]][1:MM, 20] <- as.data.frame(res)
    cosine.FP[1:MM,j]<-Object[[j]][,20]
  }

  for(j in 1:length(Object)){             ## Forward Pesrsistence time, FP deviating time and FP ratio #
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM1, function(i){
      if(abs(Object[[j]][i,19])<1.5707963268){
        Object[[j]][i,21]= TimeInterval
      }
      if(abs(Object[[j]][i,19])>=1.5707963267){
        Object[[j]][i,21]= 0
      }
      return(Object[[j]][i,21])
    })
    Object[[j]][1:MM1,21] <- as.data.frame(res)
    Object[[j]][MM1,21] <- TimeInterval

  }


  FMResultsTable<-data.frame()                                  # creating a table to store the forward migration results
  for(j in 1:length(Object)){                                   # creating values (NA and forward persistence time )for  forward persistence    (step to the forward movement)
    MM<-Step
    F.P.time<-Object[[j]][1:MM,21]
    MeanF.P.time<-round(mean(F.P.time),digits=2)             # computing mean  Forward Pesrsistence time
    F.P.time.Len<-F.P.time                                   # computing the number of FP steps to be used in computing the FP ratio
    F.P.time.Len[F.P.time.Len==0]<-NA
    F.P.time.Len<-F.P.time.Len[!is.na(F.P.time.Len)]
    F.P.time.Len<-length(F.P.time.Len)
    F.P.Ratio<- round(F.P.time.Len/MM, digits=2)               # computing FP ratio

    DD.F.P.time<-Object[[j]][1:MM,21]                        # computing the FP deviating time
    DD.F.P.time[DD.F.P.time==0]<-1
    DD.F.P.time[DD.F.P.time==TimeInterval]<-0
    DD.F.P.time[DD.F.P.time==1]<-TimeInterval
    MeanDD.F.P.time<-round(mean(DD.F.P.time),digits=2)

    FMResultsTable[1,j]<-j
    FMResultsTable[2,j]<-MeanF.P.time
    FMResultsTable[3,j]<-MeanDD.F.P.time
    FMResultsTable[4,j]<-F.P.Ratio
  }


  VelFPTable<-data.frame()               # creating a table to store the mean velocity with correspondence with the FP time
  for(j in 1:length(Object)){
    MM=Step          ###########computing the "finalID" to be used for the splitting
    FPtime<-Object[[j]][1:MM,21]
    FPtime0<-c(0,FPtime)                   #adding a "0" value in the beginning
    FPtime0[FPtime0==TimeInterval]<-NA                #replacing the "5" values with NAs
    FPtime0TF<- !is.na(FPtime0)
    MM1<-Step+1
    rowN<-c(1:MM1)
    Tval <- rowN[FPtime0TF]


    fillIdx <- cumsum(FPtime0TF)
    s<-Tval[fillIdx]                     #replacing the NAs with the previous true value
    justTrueVal<-Tval
    s[FPtime0TF]=0                       #replacing the true values with 0
    finalID<-s[-1]                       # removing the first added value

    Vel<-Object[[j]][1:MM,11]
    FP<-Object[[j]][1:MM,21]
    FP[FP==0]<-NA
    tab<-data.frame(Vel,FP,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]               # to exclude the last row since it has no velocity
    tabb<-split(tab, finalIDD)           # to avoid taking the last row

    meanVEL<-c()
    FP.time<-c()

    res1 <- sapply(1:length(tabb), function(i){
      meanVEL= round(sqrt(mean(tabb[[i]][,1])),digits=2)
      return(meanVEL)
    })
    meanVEL= res1
    meanVEL=meanVEL[-1]

    res2 <- sapply(1:length(tabb), function(i){
      FP.time=sum(tabb[[i]][,2])
      return(FP.time)
    })
    FP.time= res2
    FP.time=FP.time[-1]
    w<-which.max(FP.time)
    FMResultsTable[5,j]<-FP.time[w]



    FPtime<-Object[[j]][1:MM,21]        # computing the "meanVel.for0"
    FPtime00<-c(1,FPtime)                   #adding a "1" value in the beginning  this value should not be 0 because we will not catch 0 persistence if it is in the beginning.
    FPtime00[FPtime00==0]<-NA                #replacing the "0" values with NAs
    FPtime00TF<- !is.na(FPtime00)

    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval0 <- rowN[FPtime00TF]
    TTT<-c()
    res3 <- sapply(1:(length (Tval0)-1), function(i){
      TTT<- (Tval0[i+1]- Tval0[i])-1
      return(TTT)
    })

    TTT=res3
    TTT[TTT==0]<-NA
    F.ID.for.0.FP<-TTT[!is.na(TTT)]
    Final.ID.for.0.FP<-c()
    for (i in 1:length(F.ID.for.0.FP)){
      Final.ID.for.0.FP<-c(Final.ID.for.0.FP,rep(i,(F.ID.for.0.FP[i])))
    }

    Vel<-Object[[j]][1:MM,11]
    FP<-Object[[j]][1:MM,21]
    FP[FP==0]<-NA

    tab<-data.frame()
    tab<-data.frame(Vel,FP,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]
    tabb<-split(tab, finalIDD)     #to avoid taking the last row
    t<-tabb[[1]]
    tabb1<-cbind(t,Final.ID.for.0.FP)
    tabbb<-split(tabb1, Final.ID.for.0.FP)
    meanVel.for0<-c()

    res4 <- sapply(1:length(tabbb), function(i){
      meanVel.for0<-round(sqrt(mean(tabbb[[i]][,1])),digits=2)
      return(meanVel.for0)
    })
    meanVel.for0=res4

    zerooFP<-rep(0,length(meanVel.for0))
    FP.time1<-c(zerooFP,FP.time)
    meanVEL1<-c(meanVel.for0,meanVEL)
    colNumP<-j+j-1
    colNumV<-j+j
    rowNum<-c(1:length(meanVEL1))
    VelFPTable[rowNum,colNumP]<-FP.time1
    VelFPTable[rowNum,colNumV]<-meanVEL1
    c<-stats::cor.test( ~ VelFPTable[,j+j-1]+ VelFPTable[,j+j], method = "spearman",exact=FALSE)             #testing the correlation
    cc<-unlist(c[4])
    ccPV<-round(cc, digits = 3)
    FMResultsTable[6,j]<-ccPV               # Speed vs  persistence time  Spearman correlation



    if ( sfptPLOT == TRUE){
      if (export) grDevices::jpeg(paste0(ExpName," FP Time vs Speed",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      graphics::plot(VelFPTable[,j+j],VelFPTable[,j+j-1],pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Mean Speed during FP time (um/min)",col=color[j],las=1)
      reg<-stats::lm(VelFPTable[,j+j]~VelFPTable[,j+j-1])
      #abline(reg,untf=FALSE,col="red")
      graphics::title(main=paste0("Cell Number  ", j,"   Speed vs Forward Persistence Time"),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccPV),col.sub="red")
      if (export) grDevices::dev.off()

    }

  }

  ## All cells (FP times  vs Speed)
  allper<-VelFPTable[,1]
  allvel<-VelFPTable[,2]
  for(j in 1:length(Object)){
    allper<-c(allper,VelFPTable[,j+j-1])
    allvel<-c(allvel,VelFPTable[,j+j])
  }
  allper<-allper[!is.na(allper)]
  allvel<-allvel[!is.na(allvel)]
  all.per.vel.table<-data.frame(allper,allvel)
  reg<-stats::lm(allper~allvel)
  c<-stats::cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)                 #testing the correlation
  cc<-unlist(c[4])
  ccP<-round(cc, digits = 3)
  FMResultsTable[6,(length(Object)+1)]<-ccP                                          # Speed vs  persistence time  Spearman correlation

  if ( afptPLOT == TRUE){
    if (export) {
      grDevices::jpeg(paste0(ExpName," FP Time vs Speed - All Cells.jpg"),width = 4, height = 4, units = 'in', res = 300)
    }
    graphics::plot(allvel,allper,type="p",pch=16,ylab="FP Time (min)",xlab=" Mean Speed during FP time (um/min)",col="black",las=1)
    graphics::abline(reg,untf=FALSE,col="red")
    graphics::title("Speed vs FP Time (All cells)",cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccP),col.sub="red")
    if (export) grDevices::dev.off()

  }

  for(j in 1:length(Object)){    # calculating the Mean.Square.speed for each cell
    MM<-Step
    MM2<-MM-1
    Root.Mean.Square.Speed<-round(sqrt(mean(Object[[j]][1:MM2,11])),digits = 3)
    FMResultsTable[7,j]<-Root.Mean.Square.Speed

    mean.cosineFP<-round(mean(Object[[j]][1:MM2,20],na.rm = TRUE),digits = 3)
    FMResultsTable[8,j]<-mean.cosineFP
    s<-stats::cor.test( ~ sqrt(Object[[j]][1:MM2,11])+ Object[[j]][1:MM2,20], method = "spearman",exact=FALSE)                 #testing the correlation
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    FMResultsTable[9,j]<-VEvsCOSP

    if ( sfpPLOT == TRUE){
      if (export) grDevices::jpeg(paste0(ExpName," FP vs Speed",j,".jpg"),width = 4, height = 4, units = 'in', res = 300)
      graphics::plot(sqrt(Object[[j]][1:MM2,11]),Object[[j]][1:MM2,20],pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Instantaneous Speed (um/min)",col="black",las=1)
      reg<-stats::lm(Object[[j]][1:MM2,20]~sqrt(Object[[j]][1:MM2,11]))
      graphics::abline(reg,untf=FALSE,col="red")
      graphics::title(main=paste0("Cell Number  ", j,"   Speed vs Forward Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
      if (export) grDevices::dev.off()

    }

  }
  RM<-round(matrixStats::rowMedians(as.matrix(cosine.FP[1:(Step-1),]),na.rm = TRUE),digits=3)
  Speed<-data.frame()
  for (j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    MM<-Step
    MM2<-MM-1
    Speed[1:MM2,j]<-round(sqrt(Object[[j]][1:MM2,11]),digits = 3)
  }
  RowmeanSpeed<-round(matrixStats::rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-stats::cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  VEvsCOSP<-round(ss, digits = 3)
  FMResultsTable[9,(length(Object)+1)]<-VEvsCOSP

  if ( afpPLOT == TRUE){
    if (export) grDevices::jpeg(paste0(ExpName," All Cells FP vs Speed.jpg"),width = 4, height = 4, units = 'in', res = 300)
    graphics::plot(RowmeanSpeed,RM,pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Instantaneous Speed (um/min)",col="black",las=1)
    reg<-stats::lm(RM~RowmeanSpeed)
    graphics::abline(reg,untf=FALSE,col="red")
    graphics::title(main=paste0("All Cells Speed vs Forward Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
    if (export) grDevices::dev.off()
  }
  RM1<-round(matrixStats::rowMedians(as.matrix(FMResultsTable),na.rm = TRUE),digits=3)
  FMResultsTable[c(2:5,7:8),(length(Object)+1)]<-RM1[c(2:5,7:8)]
  FMResultsTable[1,(length(Object)+1)]<-"All Cells"
  rownames(FMResultsTable)<-c("Cell Number","Mean Forward Persist Time (min)","Mean Forward Persist Deviating Time (min)","Forward Persistence Ratio",
                              "Maximum Forward Persistence period","Forward Persistence Time vs Speed (SCC)","RMSS (um per min)","Mean Forward Angular Persistence (mean cos.F)","Instantaneous Speed vs Forward Persistence (SCC)")
  FMResultsTable<-FMResultsTable[-7,]
  object@ForMigtable=FMResultsTable

  setwd(d)
  if (export) {
    utils::write.csv(
      FMResultsTable,
      file = paste0(ExpName,"-FMResultsTable.csv")
    )
    cat(
      "Results are saved as: ",
      paste0(ExpName,"-FMResultsTable.csv" ),
      "in your directory [use getwd()]\n"
    )
  }
  return(object)
}



#' @title Forward Migration Index
#'
#' @description The FMI function automatically generates data for the forward migration index
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return  An CellMig class Object with a data frame. The data frame is stored in the FMItable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF=WSADataset[1:1000,]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' wsaTD <-FMI(wsaTD,TimeInterval=10,ExpName="ExpName", export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom SpatialTools dist2
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
FMI= function(object, TimeInterval=10, ExpName="ExpName", export=FALSE){

  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  Object<-object@preprocessedDS
  UPorDO<-object@cellpos
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-grDevices::rainbow(1023)
    colo2 <-grDevices::rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- grDevices::rainbow(Len)
  }
  for(j in 1:length(Object)){                             # Defining if the cell is ubove (1) the wound or below (0) it
    Object[[j]][,25]<-UPorDO[j]
  }

  FMIResultsTable<-data.frame()
  for(j in 1:length(Object)){                             # calculating the cumsum of distance for each cell
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM, function(i){
      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){   # upper cell going up or lower cell going down
        Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
      }

      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
        Object[[j]][i,20]<-cos(Object[[j]][i,19])
      }
      return(Object[[j]][i,20])
    })

    Object[[j]][1:MM, 20] <- as.data.frame(res)
    end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])       # finding the cordinates of the final point in the track.
    start<-cbind(0,0)
    final.dis=SpatialTools::dist2(start, end)
    alpha<-acos(Object[[j]][MM,2]/final.dis)
    bita<- 1.5707963268 - alpha
    FMI<- cos(bita)
    if(Object[[j]][1,25]==0 & Object[[j]][MM,3]>0){
      FMI<-(-FMI)
    }

    if(Object[[j]][1,25]==1 & Object[[j]][MM,3]<0){
      FMI<-(-FMI)
    }

    y=abs(Object[[j]][MM,3])
    cumDis=Object[[j]][MM1,12]
    FMIy=round(y/cumDis,digits=3)
    MMM<-round(MM/2)
    mid<-cbind(Object[[j]][MMM,2],Object[[j]][MMM,3])       # finding the cordinates of the final point in the mid track.
    start<-cbind(0,0)
    mid.dis=SpatialTools::dist2(start, mid)
    alphaM<-acos(Object[[j]][MMM,2]/mid.dis)
    bitaM<- 1.5707963268 - alphaM
    MTFMI<- cos(bitaM)
    if(Object[[j]][1,25]==0 & Object[[j]][MMM,3]>0){
      MTFMI<-(-MTFMI)
    }

    if(Object[[j]][1,25]==1 & Object[[j]][MMM,3]<0){
      MTFMI<-(-MTFMI)
    }

    p1<-Object[[j]][,20]
    returns<-subset(p1,p1<(-0.87))      # greater than 150 degrees

    yy=abs(Object[[j]][MMM,3])
    cumMDis=Object[[j]][round(MM1/2),12]
    MTFMIy=round(yy/cumMDis,digits=3)

    FMIResultsTable[1,j]<-j
    FMIResultsTable[2,j]<-round(FMI,digits=3)
    FMIResultsTable[3,j]<-round(FMIy,digits=3)
    FMIResultsTable[4,j]<-round(MTFMI,digits=3)
    FMIResultsTable[5,j]<-round(MTFMIy,digits=3)
    FMIResultsTable[6,j]<-round(abs(Object[[j]][MM,3]),digits=2)
    FMIResultsTable[7,j]<-round(length(returns))
  }
  for(j in 1:length(Object)){                            # creating values for  rel.ang.F  (step to the original)
    MM<-Step
    if((Object[[j]][1,25]==0) && (Object[[j]][MM,3]>0) || (Object[[j]][1,25]==1) && (Object[[j]][MM,3]<0)){
      FMIResultsTable[6,j]<- (-1) * FMIResultsTable[6,j]
    } else {
      FMIResultsTable[6,j]<- FMIResultsTable[6,j]
    }

  }
  RM1<-round(matrixStats::rowMedians(as.matrix(FMIResultsTable),na.rm = TRUE),digits=3)
  FMIResultsTable[,(length(Object)+1)]<-RM1
  FMIResultsTable[6,(length(Object)+1)]<-round(FMIResultsTable[6,(length(Object)+1)])
  FMIResultsTable[1,(length(Object)+1)]<-"All Cells"

  rownames(FMIResultsTable)<-c("Cell Number","FMI","FMIy", "MTFMI","MTFMIy","Deepness (um)","Number of backwards")

  if (export) {
    utils::write.csv(
      FMIResultsTable,
      file = paste0(ExpName,"-FMIResultsTable.csv")
    )
    cat(
      "Results are saved as: ",
      paste0(ExpName,"-FMIResultsTable.csv" ),
      "in your directory [use getwd()]\n"
    )
  }
  object@FMItable=FMIResultsTable

  return(object)
}



#' @title Final Results
#'
#' @description The FinRes function automatically generates a data frame that contains all the results.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ParCor A logical vector that allows generating a correlation table. Default is TRUE.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param ExpDir Directory to export the results to (if `export = TRUE`)
#' @param export if `TRUE` (default), exports function output to CSV file
#' @return  A data frame that contains all the results.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' \dontrun{
#' data(WSADataset)
#' wasDF <- WSADataset[1:1000, ]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' wsaTD <-FMI(wsaTD,TimeInterval=10,ExpName="ExpName.FNRS")
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,ExpName="ExpName.FNRS")
#' wsaTD <-FinRes(wsaTD,ExpName="ExpName",ParCor=FALSE, export=FALSE)
#' }
#'
#' @importFrom utils write.csv
#' @importFrom Hmisc rcorr
#'
#' @export
FinRes <- function(
  object,
  ExpName="ExpName",
  ExpDir=tempdir(),
  ParCor=TRUE,
  export=FALSE
) {
  # ============================================================================
  # Writing results
  # ============================================================================
  if (length(object@results) > 0) object@results <- data.frame() # reset results
  juxtaposeResults <- function(slt, obj=object) {
    # rbinds slt to object@results, removing the first row of slt
    new_results <- slot(obj, slt)
    tot_results <- slot(obj, "results")
    if (length(new_results) > 0) {
      new_results_1st_row           <- new_results[1, ]
      new_results_no_1st_row        <- new_results[-1, ]
      if (length(names(tot_results)) > 0) {
        names(new_results_no_1st_row) <- names(tot_results)
      }
      tot_results <- rbind(tot_results, new_results_no_1st_row)
      names(tot_results) <- new_results_1st_row
    }
    return(tot_results)
  }
  partial_result_slots <- c(
    "DRtable", "MSDtable", "PerAanSpeedtable", "DACtable", "VACtable",
    "ForMigtable", "FMItable"
  )
  for (i in partial_result_slots) object@results <- juxtaposeResults(i)

  # ============================================================================
  # Exporting results
  # ============================================================================
  if (export) {
    fileNamePath <- paste0(ExpDir, "/", ExpName, "-Final_Results.csv")
    utils::write.csv(object@results, file = fileNamePath)
    message("The table with the final results is saved to ", fileNamePath)
  }
  # ============================================================================
  # Calculating correlation table
  # ============================================================================
  if (ParCor) {
    R <- object@results
    R[, ] <- lapply(R[, ], function(x) as.numeric(gsub(",", ".", x)))
    Parameters.Correlation <- Hmisc::rcorr(t(R), type="spearman")
    object@parCor <- Parameters.Correlation$r
    if (export) {
      fileNamePath <- paste0(ExpDir, "/", ExpName,"-Parameters.Correlation.csv")
      utils::write.csv(Parameters.Correlation$r, file = fileNamePath)
      message("Parameters Correlation table is saved to ", fileNamePath)
    }
  }
  message("\nThese are the parameters in your final results:")
  print(rownames(object@results))
  # ============================================================================
  # Returning whole object
  # ============================================================================
  return(object)
}


#' @title PCA
#'
#' @description The CellMigPCA function automatically generates Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data.
#' @param parameters A numeric vector contains the parameters to be included in the Principal Component Analysis. These numbers can be obtained from the outcome of the FinRes() function.
#'
#' @return  PCA Graph of cells and PCA Graph of variables.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' \dontrun{
#' data(WSADataset)
#' wasDF=WSADataset[1:1000,]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' wsaTD <-FMI(wsaTD,TimeInterval=10,ExpName="ExpName")
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,ExpName="ExpName")
#' wsaTD <-FinRes(wsaTD,ExpName="ExpName",ParCor=FALSE)
#' PCAplot<-CellMigPCA(wsaTD,parameters=c(1,2))
#' }
#'
#' @importFrom FactoMineR PCA
#'
#' @export
CellMigPCA = function(object, ExpName="ExpName",
                      parameters=c(1,2,3)){

  if (!is.list(object) & !is(object, "CellMig")) {
    stop(
      "Input data must be a list. Please run the PreProcessing step first, ",
      "either rmPreProcessing() or wsaPreProcessing()"
    )
  }
  if ( length(object@results[,1])<1 ){
    stop("There are no results stored. Please run trajectory analysis first")
  }

  if ( length(parameters)<2){
    stop("At least two parameters are required to run the PCA")
  }
  df1<- object@results
  df1=df1[,-length(object@results[1,])]    #### excluding the last column since it is the avarage of all the cells
  tt<-t(df1)
  tt1=tt[,parameters]
  res <- FactoMineR::PCA(tt1)
}


#' @title PCA Clusters
#'
#' @description The CellMigPCAclust function automatically generates clusters based on the Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data.
#' @param parameters A numeric vector contains the parameters to be included in the Principal Component Analysis. These numbers can be obtained from the outcome of the FinRes() function.
#' @param export if `TRUE` (default), exports function output to CSV file
#'
#' @return  PCA Graph of cells and PCA Graph of variables.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' \dontrun{
#' data(WSADataset)
#' wasDF=WSADataset[1:1000,]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=95)
#' wsaTD <-FMI(wsaTD,TimeInterval=10,ExpName="ExpName")
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,ExpName="ExpName")
#' wsaTD <-FinRes(wsaTD,ExpName="ExpName",ParCor=FALSE)
#' PCAclust<-CellMigPCAclust(wsaTD,parameters=c(1,2))
#' }
#'
#' @importFrom FactoMineR PCA
#'
#' @export
CellMigPCAclust = function(object, ExpName="ExpName",
                           parameters=c(1,2,3), export=FALSE){

  if (!is.list(object) & !is(object, "CellMig")) {
    stop(
      "Input data must be a list. Please run the PreProcessing step first, ",
      "either rmPreProcessing() or wsaPreProcessing()"
    )
  }
  if ( length(object@results[,1])<1 ){
    stop("There are no results stored. Please run trajectory analysis first")
  }

  if ( length(parameters)<2){
    stop("At least two parameters are required to run the PCA")
  }
  df1<- object@results
  df1=df1[,-length(object@results[1,])]    #### excluding the last column since it is the avarage of all the cells
  tt<-t(df1)
  tt1=tt[,parameters]
  res <- FactoMineR::PCA(tt1)
  res.hcpc <- FactoMineR::HCPC(res)
  results<-res.hcpc$data
  results<-results[order(results$clust),]
  print(results)
  if (export) {
    utils::write.csv(
      results,
      file = paste0(ExpName,"-Clusters.csv")
    )
    cat(
      "The table of the clusters is saved as: ",
      paste0(ExpName,"-Clusters.csv"),
      " in your directory [use getwd()]\n"
    )
  }

}
