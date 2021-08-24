##
##
## ~~ All f(x) ~~
#
#

#' Return the Next Odd Integer
#'
#' Returns the smallest odd number bigger than the number(s) provided
#' as
#' the argument
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
#' cellmigRation:::NextOdd(seq(2,5,by=1))
#'
#' @keywords internal
NextOdd <- function(x) {
    y <- base::floor(x) + 1
    y <- base::ifelse(y %% 2 == 0, y + 1, y)
    return(y)
}


#' Shift Array Circularly
#'
#' Circularly shift the elements in an array by a user-defined number
#' of positions. This emulates the behavior of the corresponding
#' Matlab Circhsift function.
#'
#' @param x a character, numeric, or logical vector with at least
#' n + 1 elements
#' @param n an integer corresponding to the number of positions for
#' the shift
#'
#' @return a vector corresponding to x (same size, same class),
#' whose elements have been shifted
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::circshift(seq(1,10,by=1), -2)
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

    nu_idx <- seq(1,len,by=1)
    if(n > 0) {
        nu_idx <- c(
            seq((length(nu_idx) - nn + 1),length(nu_idx),by=1),
            seq(1,(length(nu_idx) - nn),by=1))
    } else if (n < 0) {
        nu_idx <- c(seq((nn + 1),length(nu_idx),by=1),
                    seq(1,nn,by=1))
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
#' Creates a new (molten) data matrix where all elements of y
#' are added to each row of x. Each row in x is recycled
#' for each element in y. Elements in y are added as
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
#' cellmigRation:::AddDimension(
#'     x = cbind(seq(1,4,by=1), seq(4,1,by=-1)),
#'     y = c(9, 7))
#'
#' @keywords internal
AddDimension <- function(x, y) {
    w0 <- lapply(seq(1,nrow(x),by=1), function(j1) {
        tmp <- x[j1, ]
        w1 <- lapply(seq(1,length(y),by=1), function(j2) {
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
#' Creates a Molten Hypercube with a user-defined number of
#' dimensions.
#' The values supplied by the user are used to fill each dimension.
#' All possible combination of values are included in
#' the resulting hyper cube.
#'
#' @param vals vector of values used to fill the hyper cube
#' @param dims integer indicating the number of dimensions.
#' The resulting molden data frame will have a number of columns
#' equal to dims
#'
#' @return Matrix corresponding to a molten hyper cube.
#' The number of columns is equal to dims;
#' the number of rows is equal to length(vals) ^ dims
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::MakeHypercube(seq(1,3,by=1), 3)
#'
#' @keywords internal
MakeHypercube <- function(vals, dims) {
    xi <- as.matrix(cbind(vals))
    yi <- vals
    if (dims > 1){
        for(i in seq(1,(dims-1),by=1)) {
            xi <- AddDimension(x = xi, y = yi)
        }
    } else {
        return(NULL)
    }
    return(xi)
}


#' Clean And Reformat a Numeric Matrix
#'
#' Convert any matrix-lie object to a numeric Matrix,
#' and coerces all the elements to integer.
#' Row names and column names are removed.
#'
#' @param x matrix or data.frame including numeric data
#' (or data that can be coerced to integer)
#'
#' @return numeric matrix with all its elements coerced to integer
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' tmp <- data.frame(A = c(1,2,3,4), B=c(3.1, 2.8, 3.3, 9.1), C = FALSE)
#' cellmigRation:::matfix(tmp)
#'
#' @keywords internal
matfix <- function(x) {
    xx <- as.data.frame(x)
    for ( j in seq(1,ncol(xx),by=1)) {
        xx[, j] <- as.integer(xx[, j])
    }
    colnames(xx) <- NULL
    rownames(xx) <- NULL
    return(as.matrix(xx))
}


#' Linear Convolution of a Numeric Matrix
#'
#' Performs a linear convolution of a Numeric Matrix, using a
#' user-supplied
#' linear kernel. The convolution can be executed in a column-wise
#' fashion
#' by setting the col.wise argument to TRUE. Alternatively,
#' the convolution is performed in a row-wise fashion.
#'
#' @param x numeric matrix that will be used as input for the
#' convolution;
#' this matrix typically corresponds to an image where signal
#' (high values) indicates the presence of a cell or a cell-like
#' particle
#' @param krnl numeric vector corresponding to the kernel
#' that will be used for the convolution. Briefly, the kernel
#' includes the weights that will be used to compute a weighted sum
#' at each position of the input numeric matrix
#' @param col.wise logical; shall the linear convolution be performed
#' in a column-wise or row-wise fashion?
#'
#' @return Linearly convoluted numeric matrix. The resulting matrix
#' has the same dimensions of the inut matrix
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' graphics::par(mfrow = c(1, 2))
#' tmp <- vapply(
#'     seq_len(12),
#'     function(i) {(6 + abs(i - 6)) * c(seq(1,10,by=1), seq(10,1,by=-1))},
#'     FUN.VALUE = numeric(20)
#' )
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
    tmp.i <- vapply(
        X = seq_len(floor(length(krnl) / 2)),
        FUN = function(w) {xx[, 1]},
        FUN.VALUE = numeric(length=nrow(xx))
    )
    tmp.f <- vapply(
        X = seq_len(floor(length(krnl) / 2)),
        FUN = function(w) {xx[, ncl]},
        FUN.VALUE = numeric(length=nrow(xx))
    )
    X <- cbind(tmp.i, xx, tmp.f)

    # Proceed with convolution
    Y <- do.call(rbind, lapply(seq_len(nrow(X)), function(ri) {
        vapply(seq_len(ncol(X) - length(krnl) + 1), function(ci) {
            xcoord <- seq(from = ci, to = (ci+length(krnl)-1), by = 1L)
            tmp <- X[ri, xcoord]
            as.numeric(rbind(tmp) %*% cbind(krnl))
        }, FUN.VALUE = numeric(1))
    }))

    if (col.wise)
        Y <- t(Y)

    return(Y)
}



#' Visualize a matrix image
#'
#' Shows an image representation of a numeric matrix. Typically,
#' this is a non-negative numeric matrix, where signal (high values)
#' corresponds to the presence of cells, or cell-like particles.
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
#' x <- vapply(
#'     seq_len(20),
#'     function(i) {runif(n = 20, min = 0, max = 10)},
#'     FUN.VALUE = numeric(20)
#' )
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
    xx <- t(img_mtx[seq(m,1,by=-1), ])
    graphics::image(xx, col = col, ...)
}




#' Visualize Cells in an Image Stack
#'
#' Visualize objects that were identified as cells in a given image
#' stack
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
#' # Representative output
#' x <- get(data(TrackCellsDataset))
#' VisualizeStackCentroids(tc_obj = x, stack = 2, pnt.cex = 5, offset = 1.3)
#'
#'
#' @export
VisualizeStackCentroids <- function(tc_obj, stack = 1,
                                    pnt.cex = 1.2, txt.cex = 0.9,
                                    offset = 0.18, main = NULL) {

    b <- getProcessedImages(tc_obj)$images[[stack]]
    cnt <- getImageCentroids(tc_obj)[[stack]]

    if(is.null(main)){
        main <- paste0("Stack num. ", stack)
    }

    VisualizeImg(img_mtx = b, las = 1, main = main)
    VisualizeCntr(
        centroids = cnt, width_px = ncol(b), height_px = nrow(b),
        pnt.cex = pnt.cex, txt.cex = txt.cex, offset = offset
    )

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
#'
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x1 <- data.frame(
#'     row = c(50, 80, 20, 65, 99),
#'     col = c(15, 25, 50, 65, 86))
#' plot(2, 2, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", las = 2)
#' cellmigRation:::VisualizeCntr(x1, width_px = 100, height_px = 100)
#'
#'
#' @importFrom graphics points text
#'
#' @keywords internal
VisualizeCntr <- function(
    centroids, width_px, height_px, pnt.cex = 1.2,
    txt.cex = 0.9, offset = 0.18, col = "red2")
{
    cnt <- centroids
    points(
        x = ((cnt$col - 1) / (width_px - 1)),
        y = 1-((cnt$row - 1) / (height_px - 1)),
        cex = pnt.cex, col = col
    )

    text(
        x = ((cnt$col - 1) / (width_px - 1)),
        y = 1-((cnt$row - 1) / (height_px - 1)),
        labels = seq(1,nrow(cnt),by=1), font = 4,
        cex = txt.cex, col = col,
        pos = 4, offset = offset
    )

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
#' @param col.untracked color of the points that were not tracked
#' further, e.g.: "gray45"
#' @param main string used as plot title, can be NULL
#'
#' @return None
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x <- get(data(TrackCellsDataset))
#' visualizeCellTracks(tc_obj = x, stack = 2)
#'
#'
#' @importFrom stats setNames
#'
#' @export
visualizeCellTracks <- function(
    tc_obj, stack = 1, pnt.cex = 1.2, lwd = 1.6,
    col = "red2", col.untracked = "gray45", main = NULL
) {

    if (is.null(main)) {
        main <- paste0("Tracks of Cells in Stack num. ", stack)
    }

    # Retrieve anr show image
    b <- getProcessedImages(tc_obj)$images[[stack]]
    VisualizeImg(img_mtx = b, las = 1, main = main)

    # Rerieve tracks / centroids
    #cnt <- tracked_cells$centroids[[stack]]
    cnt <- getCellTracks(tc_obj)
    cnt <- cnt[cnt[, 3]>= stack, ]
    cids_stack <- cnt[cnt[, 3] == stack, 4]
    cnt_plt <- cnt[cnt[, 4] %in% cids_stack, ]
    id_2pls <- unique(cnt_plt[duplicated(cnt_plt[,4]), 4])
    id_1shr <- unique(cnt_plt[!cnt_plt[, 4] %in% id_2pls, 4])

    if(length(id_1shr) > 0) {
        tmp_cnt <- cnt_plt[cnt_plt[, 4] %in% id_1shr, ]
        tmp_cnt <- tmp_cnt[tmp_cnt[, 3] == stack, ]
        tmp_cnt <- setNames(
            as.data.frame(tmp_cnt),
            nm = colnames(getImageCentroids(tc_obj)[[1]])
        )

        VisualizeCntr(
            centroids = tmp_cnt, width_px = ncol(b), height_px = nrow(b),
            pnt.cex = pnt.cex, txt.cex = 0.00001, offset = 0.1,
            col = col.untracked
        )

    }

    if(length(id_2pls) > 0) {
        tmp_cnt <- cnt_plt[cnt_plt[, 4] %in% id_2pls, ]
        tmp_cnt <- setNames(
            as.data.frame(tmp_cnt),
            nm = colnames(getImageCentroids(tc_obj)[[1]])
        )

        # Use dedicated f(x)
        visualizeTrcks(
            tracks = tmp_cnt, width_px = ncol(b), height_px = nrow(b),
            i.slice = stack, pnt.cex = pnt.cex, lwd = lwd, col = col
        )

    }

    # DOne , no return needed
    # return()
}

#' Visualize Cell Tracks
#'
#' Annotates an image with cell centroids by adding cell ROIs
#' and drawing cell tracks
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
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
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
visualizeTrcks <- function(
    tracks, width_px, height_px, i.slice = 1,
    pnt.cex = 1.2, lwd = 1.2, col = "red")
{
    allcnt <- setNames(
        as.data.frame(tracks),
        nm = c("row", "col", "slice", "cell")
    )
    cnt <- allcnt[allcnt$slice == i.slice, ]
    cnt <- cnt[order(cnt$cell),]


    all_cells_slice <- cnt$cell
    cKeep <- vapply(all_cells_slice, function(jj) {
        sum(allcnt$cell == jj) > 1
    }, FUN.VALUE = logical(1))

    # Cell outile

    graphics::points(
        x = ((cnt$col - 1) / (width_px - 1)),
        y = 1-((cnt$row - 1) / (height_px - 1)),
        cex = pnt.cex, col = ifelse(cKeep, col, "gray75")
    )

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
#' Import a .tif stack containing fluorescently labeled point particles
#' to be tracked
#'
#' @param tiff_file path to a TIFF file to be read in
#' @param experiment string, a label to describe the experiment
#' (optional)
#' @param condition string, a label to describe the experimental
#' condition
#' @param replicate string, a label to identify the replicate
#' (optional)
#'
#' @note `experiment`, `condition` and `replicate` are optional
#' arguments
#' and can be NULL.
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' # Let `path/to/tiff_file.tiff` be the path to tiff file we want to
#' # import. If an error is thrown, NULL is returned.
#' x <- LoadTiff(tiff_file = "path/to/tiff_file.tiff")
#'
#' @importFrom tiff readTIFF
#'
#' @export
LoadTiff <- function(
    tiff_file, experiment = NULL, condition = NULL, replicate = NULL){
    errTR <- tryCatch({
        myIMG <- suppressMessages(
            suppressWarnings(
                tiff::readTIFF(
                    source = tiff_file, native = FALSE, all = TRUE,
                    info = TRUE, as.is = TRUE)))

        if (!is.list(myIMG))
            myIMG <- list(myIMG)

        if (is.null(experiment)) {
            experiment <- NA
        } else {
            experiment <- tryCatch(
                as.character(experiment[1]), error = function(e) NA
            )
        }

        if (is.null(replicate)) {
            replicate <- NA
        } else {
            replicate <- tryCatch(
                as.character(replicate[1]), error = function(e) NA
            )
        }

        if (is.null(condition)) {
            condition <- NA
        } else {
            condition <- tryCatch(
                as.character(condition[1]), error = function(e) NA
            )
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
            vapply(
                seq_len(ncol(x)),
                function(ii) {as.numeric(x[,ii])},
                FUN.VALUE = numeric(nrow(x))
            )
        })}, silent = TRUE)

        img_list <- list(
            images = FinalImage,
            dim = list(
                NumberImages = NumberImages, width_m = mImage, height_n=nImage
            ),
            attributes = InfoImage
        )

        Y <- new(Class = "trackedCells", img_list)

        # Attach labels
        my.meta <- list(
            tiff_file = sub("^.*[/]([^/]+$)", "\\1", tiff_file),
            experiment = experiment,
            condition = condition,
            replicate = replicate
        )
        Y <- setTrackedCellsMeta(x = Y, meta = my.meta)
        NULL
    }, error = function(e) { e })

    if (!is.null(errTR)) {
        message(errTR)
        return(NULL)
    }
    return(Y)
}



#' Validate Centroids
#'
#' Validate parameters used to identify cells in a image stack.
#' A figure containing current image frame with identified particles
#' labeled with circles and numerical tags is generated. This function
#' is included for consistency and compatibility reasons
#' with the original fastTracks software (Matlab). Also, consider
#' using
#' VisualizeStackCentroids() or visualizeCellTracks() instead.
#'
#' @param stack stack of images to be evaluated
#' @param slice index of the frame within the stack to be evaluated
#' @param lobject integer, length in pixels somewhat larger than
#' a typical object (cell)
#' @param threshold the minimum brightness of a pixel that might
#' be local maxima. NOTE:
#' Make it big and the code runs faster but you might miss some
#' particles.
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
#' x <- get(data(TrackCellsDataset))
#' x <- getCellImages(x)
#' x$images[[1]] <- x$images[[1]][seq(110,160,by=1), seq(100,160,by=1)]
#' cellmigRation:::CentroidValidation(x, slice = 1, lobject =10, threshold = 5)
#'
#' @importFrom graphics box axis points text
#'
#' @keywords internal
CentroidValidation <- function(
    stack, slice, lobject, threshold,
    pnt.cex = 1.2, txt.cex = 0.85, offset = 0.18
) {
    a <- stack$images[[slice]]
    b <- bpass(
        image_array = a, lnoise = 1, lobject = lobject, threshold = threshold
    )
    pk = pkfnd(im = b, th = threshold, sz = lobject+1)
    cnt = cntrd(im = b, mx = pk, sz = lobject+1)

    VisualizeImg(b, axes = FALSE)
    graphics::box()
    my_xax <- ncol(b)
    my_yax <- nrow(b)

    my_xax <- unique(c(seq(1, my_xax, by = 100), my_xax))
    my_yax <- unique(c(seq(1, my_yax, by = 100), my_yax))

    axis(side = 1,at = ((my_xax - 1) / max(my_xax)), labels = my_xax)
    axis(side = 2,at = 1-((my_yax - 1) / max(my_yax)), labels = my_yax, las=1)

    points(
        x = ((cnt$col - 1) / (ncol(b) - 1)),
        y = 1-((cnt$row - 1) / (nrow(b) - 1)),
        cex = pnt.cex, col = "red2"
    )

    points(
        x = ((cnt$col - 1) / (ncol(b) - 1)),
        y = 1-((cnt$row - 1) / (nrow(b) - 1)),
        cex = 1.2, col = "red2")

    text(
        x = ((cnt$col - 1) / (ncol(b) - 1)),
        y = 1-((cnt$row - 1) / (nrow(b) - 1)),
        labels = seq(1,nrow(cnt),by=1), font = 4,
        cex = txt.cex, col = "red2",
        pos = 4, offset = offset)

    return(cnt)
}



#' Perform a bandpass by convolving with an appropriate kernel
#'
#' Implements a real-space bandpass filter that suppresses pixel noise
#' and
#' long-wavelength
#' image variations while retaining information of a characteristic
#' size.
#' First, a lowpassed image is produced by convolving the original
#' with
#' a gaussian.
#' Next, a second lowpassed image is produced by convolving the
#' original
#' with a
#' boxcar function. By subtracting the boxcar version from
#' the gaussian version,
#' we are using the boxcar version to perform a highpass.
#' This code 'bpass.pro' is copyright 1997, John C. Crocker and
#' David G. Grier.    It should be considered 'freeware'- and may be
#' distributed freely in its original form when properly attributed.
#'
#' @param image_array Numeric matrix corresponding to the image to be
#' filtered
#' @param lnoise Characteristic lengthscale of noise in pixels.
#' @param lobject Integer length in pixels somewhat larger than
#' a typical object
#' @param threshold By default, after the convolution, any negative
#' pixels
#' are reset to 0. Threshold changes the threshhold for setting
#' pixels
#' to 0. Positive values may be useful for removing
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
#' x0 <- cbind(
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 1, 1, 4, 2, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 1, 2, 6, 5, 3, 0, 0, 0, 1, 0, 0),
#'     c(0, 0, 0, 0, 5, 5, 6, 8, 6, 1, 0, 0, 6, 2, 0),
#'     c(0, 0, 2, 5, 8, 7, 3, 5, 1, 0, 0, 0, 6, 2, 0),
#'     c(0, 0, 1, 5, 8, 7, 4, 5, 2, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 5, 8, 7, 4, 5, 2, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 1, 4, 5, 2, 4, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 2, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0),
#'     c(9, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1),
#'     c(2, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1),
#'     c(0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
#' y0 <- cellmigRation:::bpass(x0, lnoise = 1, lobject = 5, threshold = 1)
#' par(mfrow = c(1, 2))
#' image(x0); title("original")
#' image(y0); title("after bpass")
#'
#' @keywords internal
bpass <- function(image_array, lnoise, lobject = NULL, threshold)
{

    cstm_normalize <- function(x) { x/sum(x) }

    # Make kernel (linear)
    gaussian_kernel <- cstm_normalize(
        exp(-(seq(-2.5, 2.5, length.out = ((10 * lnoise) + 1))^2))
    )

    if (!is.null(lobject)) {
        boxcar_kernel <- cstm_normalize(rep(1, times = (2 * lobject) + 1))
    }


    gconv <- LinearConv2(t(image_array), gaussian_kernel)
    gconv <- LinearConv2(t(gconv), gaussian_kernel)


    if (!is.null(lobject)) {
        bconv <- LinearConv2(t(image_array), boxcar_kernel)
        bconv <- LinearConv2(t(bconv), boxcar_kernel)
        filtered <- gconv - bconv
    } else {
        filtered <- gconv
    }

    # Zero out the values on the edges to signal that they're not useful.
    lzero <- max(lobject, ceiling(5*lnoise))

    filtered[seq(1, (round(lzero)), by = 1),] <- 0
    my.sel.range <- seq((nrow(filtered)-round(lzero)+1), nrow(filtered), by=1)
    filtered[my.sel.range, ] <- 0

    filtered[, seq(1, (round(lzero)), by = 1)] <- 0
    my.sel.range <- seq((ncol(filtered)-round(lzero)+1), ncol(filtered), by=1)
    filtered[, my.sel.range] <- 0

    # Zero all values below threshold
    filtered[filtered < threshold] <- 0

    return(filtered)
}



#' Calculates Centroids
#'
#' Calculates the centroid of bright spots to sub-pixel accuracy.
#' Inspired by Grier &
#' Crocker's feature for IDL, but greatly simplified and optimized for
#' MATLAB, and then
#' further ported to R.
#' CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#'
#' @param im numeric matrix corresponding to the image to process
#' @param mx location of local maxima to pixel-levels accuracy
#' @param sz diameter of the window over which to average to calculate
#' the centroid. should be big enough.
#' @param interactive numeric; if set to 1 (or any positive number),
#' an image showing the
#' computed centroids will be visualized
#'
#' @return a data.frame with 4 columns, containing, x, y, brightness,
#' and
#' the square of the radius of gyration for each cell.
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- getCellImages(get(data(TrackCellsDataset)))
#' x0 <- x0$images[[1]][seq(80,150,by=1), seq(80,150,by=1)]
#' b <- cellmigRation:::bpass(image_array = x0, lnoise = 2,
#'                            lobject = 15, threshold = 1)
#' pk <- cellmigRation:::pkfnd(b, th = 2, sz = 5)
#' cnt <- cellmigRation:::cntrd(im = b, mx = pk, sz = 5)
#' cnt
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
        message(
            'there were no positions inputted into cntrd. ',
            'check your pkfnd theshold'
        )
        return(NULL)
    }

    # Compute
    r <- (sz+1)/2

    # Create mask window around trial location over which to calculate
    # the centroid
    m <- 2*r
    x <- seq(0, (m-1), by = 1)
    cent <- (m-1)/2
    x2 <- (x-cent) ^ 2
    dst <- do.call(rbind, lapply(seq(1,m,by=1), function(i){
        sqrt((i-1-cent)^2+x2)
    }))

    ind <- dst < r
    msk <- ind * 1 # convert ind to numeric
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
    xl <- do.call(rbind, lapply(seq(1,(2*r),by = 1), function(j) {
        (seq(1, (2*r), by = 1))
    }))
    yl <- t(xl)

    #loop through all of the candidate positions
    pts <- list()
    for (i in seq(1, nmx, by = 1)) {
        # create a small working array around each candidate location,
        # and apply the window function
        tmp <- msk * im[
            seq((mx$row[i] - floor(r) + 1),(mx$row[i] + floor(r)),by=1),
            seq((mx$col[i] - floor(r) + 1),(mx$col[i] + floor(r)),by=1)
        ]

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
        pts[[length(pts) +1 ]] <- data.frame(
            row = mx$row[i]+yavg-r,
            col = mx$col[i] + xavg - r,
            norm = norm ,
            rg = rg
        )
        if (as.numeric(interactive) > 0) {
            VisualizeImg(img_mtx = tmp, axes = FALSE)
            box()
            axis(
                side = 1, at = seq(0, 1, length.out = ncol(tmp)),
                labels = seq((mx$row[i]-floor(r)+1),(mx$row[i]+floor(r)),by=1))
            axis(
                side = 2, at = seq(0, 1, length.out = nrow(tmp)),las = 1,
                labels = seq((mx$col[i]+floor(r)),(mx$col[i]-floor(r)+1),by=1))
            title(
                main = paste0("Cell number #", i), ylab = "y_pixel",
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
#' Finds local maxima in an image to pixel level accuracy.
#' This provides a rough guess
#' of particle centers to be used by cntrd(). Inspired by
#' the lmx subroutine of
#' Grier and Crocker's.
#' CREATED: Eric R. Dufresne, Yale University, Feb 4 2005.
#'
#' @param im image to process, particle should be bright spots on
#' dark background
#' with little noise ofen an bandpass filtered brightfield image
#' @param th the minimum brightness of a pixel that might be local
#' maxima.
#' NOTE: Make it big and the code runs faster but you might miss some
#' particles.
#' Make it small and you'll get everything and it'll be slow.
#' @param sz if your data is noisy, (e.g. a single particle has
#' multiple
#' local maxima),
#' then set this optional keyword to a value slightly larger than
#' the diameter of your blob.
#' If multiple peaks are found withing a radius of sz/2 then the code
#' will keep only the brightest.
#' Also gets rid of all peaks within sz of boundary
#'
#' @return a numeric data.frame with two columns,
#' with the coordinates of local maxima
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- getCellImages(get(data(TrackCellsDataset)))
#' x0 <- x0$images[[1]][seq(80,150,by=1), seq(80,150,by=1)]
#' b <- cellmigRation:::bpass(image_array = x0, lnoise = 2,
#'                            lobject = 15, threshold = 1)
#' pk <- cellmigRation:::pkfnd(b, th = 2, sz = 5)
#' pk
#'
#' @keywords internal
pkfnd <- function(im, th, sz=NULL)
{
    # find all the pixels above threshold
    ind <- im > th
    nr <- nrow(im)
    nc <- ncol(im)

    # melt to have a list of points above threshold
    ind2 <- reshape2::melt(data = ind, varnames = c("row", "col"))
    ind2 <- ind2[ind2$value, ]

    # check each pixel above threshold to see if it's brighter than
    # its neighbors
    # THERE'S GOT TO BE A FASTER WAY OF DOING THIS.
    # I'M CHECKING SOME MULTIPLE TIMES, BUT THIS DOESN'T SEEM THAT SLOW
    # COMPARED TO THE OTHER ROUTINES, ANYWAY.
    keep <- list()
    for(i in seq(1, nrow(ind2), by = 1)) {
        ri <- ind2$row[i]
        ci <- ind2$col[i]

        if (ri>1 & ri<nr & ci>1 & ci<nc) {
            z1 <- im[ri, ci]
            z2 <- as.numeric(
                im[seq((ri-1),(ri+1),by=1),
                    seq((ci-1),(ci+1), by = 1)])

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

    # if size is specified, then get ride of pks within size of boundary
    # (i.e., a margin from image edges)
    if (!is.null(sz) & npks>0) {
        # throw out all pks within sz of boundary;
        keep <-
            mx$row > sz & mx$row < (nr - sz + 1) &
            mx$col > sz & mx$col < (nc - sz + 1)
        mx<-mx[keep,]
    }

    # prevent from finding peaks within size of each other
    npks <- nrow(mx)
    if (!is.null(sz) & npks > 1) {
        # CREATE AN IMAGE WITH ONLY PEAKS
        mask <- matrix(FALSE, nrow = nrow(im), ncol = ncol(im))
        for(i in seq(1,nrow(mx), by=1)) {
            mask[mx$row[i], mx$col[i]] <- TRUE
        }

        tmp <- matrix(0, nrow=nrow(im), ncol = ncol(im))
        tmp[mask] <- im[mask]

        # LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK THE BRIGHTEST
        for (i in seq(1,nrow(mx), by=1)) {
            astep <- floor(sz/2)
            roi <- tmp[
                seq((mx$row[i] - astep), (mx$row[i] + astep), by=1),
                seq((mx$col[i] - astep), (mx$col[i] + astep), by=1)
            ]
            imax <- which.max(roi)
            chkrow <- imax %% nrow(roi)
            myrow <- ifelse(chkrow == 0, nrow(roi), chkrow)
            mycol <- ifelse(
                chkrow == 0, floor(imax/nrow(roi)), floor(imax/nrow(roi)) + 1
            )
            mv <- roi[myrow, mycol]

            tmp[
                seq((mx$row[i] - astep), (mx$row[i] + astep), by=1),
                seq((mx$col[i] - astep), (mx$col[i] + astep), by=1)
            ] <- 0
            tmp[
                (mx$row[i] - astep + myrow - 1),
                (mx$col[i] - astep + mycol - 1)
            ] <- mv
        }

        ind <- tmp > th
        nr <- nrow(tmp)
        nc <- ncol(tmp)

        # melt to have a list of points above threshold
        ind.f <- reshape2::melt(data = ind, varnames = c("row", "col"))
        ind.f <- ind.f[ind.f$value, c(1,2)]
        rownames(ind.f) <- NULL

        return(ind.f)
    } else {
        return(NULL)
    }
}


#' Build a Centroid Array
#'
#' Create an array containing centroid data for particles identified
#' in each frame of the imported TIFF image stack
#'
#' @param stack 3D matrix loaded to workspace from .tif stack
#' @param lobject Integer length in pixels somewhat larger
#' than a typical object
#' @param threshold the minimum brightness of a pixel that might be
#' local maxima
#' @param dryrun logical, shall the execution be skipped
#'
#' @return data.frame of centroids, with 4 columns corresponding to
#' x-position of centroid,
#' y-postion of centroid, brightness, and square of the radius
#' of gyration
#'
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' # by default, the dryrun argument is set to FALSE
#' df <- get(data(TrackCellsDataset))
#' x0 <- getCellImages(df)
#' y0 <- cellmigRation:::CentroidArray(x0, 16, 10, TRUE)
#' y0
#'
#'
#' @keywords internal
CentroidArray <- function(stack, lobject, threshold, dryrun=FALSE)
{
    #determine the number of slices within the stack
    m <- stack$dim$width_m
    n <- stack$dim$height_n
    p <- stack$dim$NumberImages

    centroid <- list()

    if (dryrun) {
        Y <- list(
            data.frame(
                row=c(29.2, 131.8), col=c(26.1, 27.3),
                norm = c(2909, 1021), rg= c(26.1, 27.1)),
            data.frame(
                row=c(24.5, 133.1), col=c(29.1, 32.3),
                norm=c(3221, 1551), rg=c(25.1, 29.1)))
        return(Y)
    }

    for(i in seq(1,p,by=1)) {
        a <- stack$images[[i]]
        b <- bpass(
            image_array = a,
            lnoise = 1,
            lobject = lobject,
            threshold = quantile(a, 0.25)
        ) # maybe set to 0 or to threshold

        pk <- pkfnd(b, threshold, lobject+1)
        cnt <- cntrd(im = b, mx = pk, sz = lobject + 1)

        if(is.null(cnt) || nrow(cnt) < 1) {
            message(
                paste0(
                    'No centroids detectd in frame ', i,
                    '...\nCheck nuclei validation settings for this frame.'
                )
            )
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
#' @return data.frame including two columns: MPOS indicates
#' the centroid position of a particle,
#' and LEN indicates the diameter size
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::DetectRadii(c(0,0,1,1,0,1,1,1,1,0,0, 1,0,0,1,1))
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

        for (j in seq(2,length(xx),by=1)) {
            if (xx[j] == (xx[(j-1)] + 1)) {
                LN <- LN + 1
                p1 <- xx[j]
                if (j == length(xx)) {
                    yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
                    radii[[length(radii) + 1]]    <- yy
                }
            } else {
                yy <- data.frame(MPOS = mean(c(p0, p1)), LEN = LN)
                radii[[length(radii) + 1]]    <- yy
                p0 <- xx[j]
                p1 <- xx[j]
                LN <- 1
            }
        }

    } else if (length(xx) == 1) {

        yy <- data.frame(MPOS = xx, LEN = 1)
        radii[[length(radii) + 1]]    <- yy
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
#' @param px.margin integer, number of pixels used as margin
#' while searching/filtering for neighboring particles
#' @param min.px.diam integer, minimum diameter of a particle (cell).
#' Particles with a diameter smaller than min.px.diam are discarded
#' @param quantile.val numeric, must be bigger than 0 and smaller
#' than 1.
#' Quantile for discriminating signal and background;
#' only pixels with intensity higher than the corresponding
#' quantile will count as signal while estimating particle diameters
#' @param plot logial, shall a histogram of the distribution
#' of diameters be shown
#'
#' @return list including summary stats and data about the particles
#' found in the image
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
EstimateDiameterRange <- function(
    x, px.margin = 2, min.px.diam = 5,
    quantile.val = 0.99, plot = TRUE)
{
    QNTS <- as.numeric(quantile(x, probs = quantile.val[1]))

    # Adjust if quantile.val is too low (few cells)
    tmp.xx <- as.numeric(x)
    max.sig <- max(tmp.xx, na.rm = TRUE)
    min.sig <- min(tmp.xx, na.rm = TRUE)
    if (QNTS == min.sig && max.sig > min.sig) {
        QNTS <- mean(
            c(min.sig, min(tmp.xx[tmp.xx > min.sig], na.rm = TRUE)),
            na.rm = TRUE
        )
    }

    B <- x
    B[B < QNTS] <- 0
    B[B >= QNTS] <- 1

    rdds <- do.call(rbind, lapply(seq(1,ncol(B),by=1), function(ii) {

        out <- DetectRadii(B[,ii])
        if (!is.null(out)) {
            data.frame(RPOS = out$MPOS, CPOS = ii, LEN = out$LEN)
        }
    }))
    rdds$KEEP <- TRUE

    for (j in seq(1,nrow(rdds))) {
        if (rdds$KEEP[j]){

            tdm <- ( 2 * px.margin)    + rdds$LEN[j]
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

    yy <- list(
        estim.cell.num = sum(FINL$KEEP),
        q50.diam = median(FINL$LEN, na.rm = TRUE),
        q75.diam = as.numeric(
            quantile(FINL$LEN, na.rm = TRUE, probs = 0.75)
        ),
        q90.diam = as.numeric(
            quantile(FINL$LEN, na.rm = TRUE, probs = 0.90)
        ),
        q95.diam = as.numeric(
            quantile(FINL$LEN, na.rm = TRUE, probs = 0.95)
        ),
        raw = FINL
    )

    if (plot) {
        try(hist(FINL$LEN, breaks = seq(
            min(FINL$LEN, na.rm = TRUE),
            max(FINL$LEN, na.rm = TRUE), length.out = 20
        ),
        xlab = "Particle Diameter", las = 1, main = "Diam. Distribution",
        col = "aquamarine3"), silent = TRUE); box()
    }

    return(yy)
}


#' Initialize Tracking parameters
#'
#' Initialize parameter variables used for the tracking
#'
#' @param dd numeric, value of the dd param
#' @param params a list containing a few tracking parameters that are
#' needed for the analysis
#'
#' @return a list including parsed arguments
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#'
#' @examples
#' cellmigRation:::initializeTrackParams(0, NULL)
#'
#'
#' @keywords internal
initializeTrackParams <- function(dd = 0, params){

    y <- list()

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
    y$memory_b <- memory_b
    y$goodenough <- goodenough
    y$dim <- dim
    y$quiet <- quiet
    y$force_exec <- force_exec

    return(y)
}


#' Print Warning Messages
#'
#' Print warning messages to the console when an issue is encountered
#'
#' @param warn_log list, including the warning logs
#' @param quiet logical, shall the warning be printed
#'
#' @return log of warning messages is returned
#'
#'
#' @examples
#' cellmigRation:::warnMessage(list("Hello world"), FALSE)
#'
#'
#' @keywords internal
warnMessage <- function(warn_log, quiet = FALSE) {

    warn_cycles <- NULL
    if (is.list(warn_log) && length(warn_log) > 0) {

        warn_cycles <- sort(unique(do.call(c, warn_log)))

        if (!quiet) {
            message(
                paste0(
                    "Difficult combinatorics encountered ",
                    "while processing slide(s): ",
                    paste(warn_cycles, collapse = ", "), "."
                )
            )
        }
    }
    return(warn_cycles)
}




#' Track Hypercube Build
#'
#' Build an hypercube used for the tracking step
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param env environment, including all objects used for the tracking
#'
#' @return NULL is returned while objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#' @examples
#' cellmigRation:::trackHypercubeBuild(data.frame(1), new.env())
#'
#'
#' @keywords internal
trackHypercubeBuild <- function(xyzs, env) {

    stopifnot(is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env

    tryCatch({
        suppressWarnings({
            # construct the vertices of a 3x3x3... d-dimensional hypercube
            eee$numbs <- c(0,1,2)
            eee$cube <- MakeHypercube(vals = eee$numbs, dims = eee$dim)

            # calculate a blocksize which may be greater than maxdisp,
            # but which keeps nblocks reasonably small.
            eee$volume <- 1
            for (d in seq(0,(eee$dim-1), by=1)) {
                eee$minn <- min(xyzs[eee$w, (d+1)])
                eee$maxx = max(xyzs[eee$w, (d+1)])
                eee$volume <- eee$volume * (eee$maxx-eee$minn)
            }

            # volume;
            eee$blocksize <-
                max(c(
                    eee$maxdisp,
                    ((eee$volume)/(20*eee$ngood))^(1.0/eee$dim)))
        })
    }, error = function(e) {
        message(
            "An issue may have occurred with the trackHypercubeBuild function")
    })
    return(NULL)
}







#' Internal Permutation
#'
#' Perform Internal Permutation as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned while objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#' @examples
#' cellmigRation:::internalPermutation(data.frame(1), 1, new.env())
#'
#'
#' @keywords internal
internalPermutation <- function(xyzs, maxdisp, env){

    stopifnot(
        is.data.frame(xyzs), is.numeric(maxdisp), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the internalPermutation function"

    tryCatch({
        if (eee$pt[(eee$who+1)] != eee$fi[(eee$who+1)]){

            eee$w <- which(
                eee$ok[eee$h[
                    seq((eee$pt[(eee$who+1)]+1),
                        (eee$fi[(eee$who+1)]),by=1)]]!=0
            ) # check -1
            eee$ngood <- length(eee$w)
            if (eee$ngood > 0){
                if (eee$pt[(eee$who+1)] != eee$st[(eee$who+1)]-1) {
                    eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 1
                }
                eee$pt[(eee$who+1)] <- eee$pt[(eee$who+1)] + eee$w[1]
                eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 0
                if (eee$who == (eee$nnew - 1)){
                    eee$ww <- which(eee$lost == 0)
                    eee$dsq <- sum(eee$lensq[eee$pt[eee$ww]]) +
                        (eee$losttot * eee$maxdisq)

                    if (eee$dsq < eee$mndisq){
                        eee$minbonds <- eee$pt[eee$ww]
                        eee$mndisq <- eee$dsq
                    }
                } else {
                    eee$who <- eee$who+1
                }
            } else {
                if(!eee$lost[(eee$who+1)] & (eee$losttot != eee$nlost)){
                    eee$lost[(eee$who+1)] <- 1
                    eee$losttot <- eee$losttot + 1
                    if (eee$pt[(eee$who+1)] != eee$st[(eee$who+1)] - 1){
                        eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 1
                    }
                    if (eee$who == (eee$nnew-1)){
                        eee$ww <- which(eee$lost == 0)
                        eee$dsq <- sum(eee$lensq[eee$pt[eee$ww]]) +
                            (eee$losttot * eee$maxdisq)
                        if (eee$dsq < eee$mndisq){
                            eee$minbonds <- eee$pt[eee$ww]
                            eee$mndisq <- eee$dsq
                        }
                    } else {
                        eee$who <- eee$who + 1
                    }

                } else {
                    if (eee$pt[(eee$who+1)] != (eee$st[(eee$who+1)] - 1)){
                        eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 1
                    }
                    eee$pt[(eee$who+1)] <- eee$st[(eee$who+1)] - 1
                    if (eee$lost[(eee$who+1)]){
                        eee$lost[(eee$who+1)] <- 0
                        eee$losttot <- eee$losttot -1
                    }
                    eee$who <- eee$who - 1
                }
            }

        } else {
            if (!eee$lost[(eee$who+1)] && (eee$losttot != eee$nlost)){
                eee$lost[(eee$who+1)] <- 1
                eee$losttot <- eee$losttot + 1
                if (eee$pt[(eee$who+1)] != eee$st[(eee$who+1)]-1){
                    eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 1
                }
                if (eee$who == (eee$nnew - 1)) {
                    eee$ww <- which(eee$lost == 0)
                    eee$dsq <- sum(eee$lensq[eee$pt[eee$ww]]) +
                        (eee$losttot * eee$maxdisq)

                    if (eee$dsq < eee$mndisq){
                        eee$minbonds <- eee$pt[eee$ww]
                        eee$mndisq <- eee$dsq
                    }
                } else {
                    eee$who <- eee$who + 1
                }
            } else {
                if (eee$pt[(eee$who+1)] != eee$st[(eee$who+1)] - 1){
                    eee$ok[eee$h[eee$pt[(eee$who+1)]]] <- 1
                }
                eee$pt[(eee$who+1)] <- eee$st[(eee$who+1)] - 1
                if (eee$lost[(eee$who+1)]){
                    eee$lost[(eee$who+1)] <- 0
                    eee$losttot <- eee$losttot - 1
                }
                eee$who <- eee$who -1
            }
        }
    }, error = function(e) {
        message(MZ3)
    })
    return(FALSE)
}








#' Run Tracking Permutation
#'
#' Perform Internal Permutation as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param nc numeric, value of the nc parameter
#' @param env environment, including all objects used for the tracking
#' @param i integer, index of the current cycle
#'
#' @return FALSE is returned while objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::runTrackingPermutation(data.frame(1), 1, 1, 1, new.env())
#'
#'
#' @keywords internal
runTrackingPermutation <- function(xyzs, maxdisp, nc, i, env) {

    stopifnot(
        is.numeric(nc), is.numeric(maxdisp),
        is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    eee$nc <- nc
    MZ1 <- "Excessive Combinitorics encountered while processing slide "
    MZ2 <- ". Quitting now... Try using a smaller maxdisp."
    MZ3 <- "An issue may have occurred with the runTrackingPermutation function"
    out.value <- FALSE

    tryCatch({

        eee$w <- which(eee$bmap == (-1)*(eee$nc))
        eee$nbonds <- length(eee$w)
        eee$bonds <- eee$mbonds[eee$w,]
        eee$lensq <- eee$bondlen[eee$w]

        eee$pq <- sort(eee$bonds[,1])
        eee$st <- order(eee$bonds[,1])
        #%implanting unq directly
        eee$arr <- eee$bonds[,1]
        eee$q <- eee$arr[eee$st]
        eee$indices <- which(eee$q != circshift(eee$q,-1))
        eee$count <- length(eee$indices)
        if (eee$count > 0) {
            eee$un <- eee$st[eee$indices]
        } else {
            eee$un <- length(eee$q) - 1
        }
        eee$uold <- eee$bonds[eee$un,1]
        eee$nold <- length(eee$uold)

        #%implanting unq directly
        eee$indices <- which(eee$bonds[, 2] != circshift(eee$bonds[, 2], -1))
        eee$count <- length(eee$indices)
        if (eee$count > 0){
            eee$un <- eee$indices
        } else {
            eee$un <- length(eee$bonds[,2]) -1
        }

        eee$unew <- eee$bonds[eee$un,2]
        eee$nnew <- length(eee$unew)

        if (eee$nnew > 5){
            eee$rnsteps <- 1
            for (ii in seq(1,eee$nnew,by=1)){

                eee$rnsteps <- eee$rnsteps * length(
                    which(eee$bonds[,2] == eee$unew[ii])
                )
                if (eee$rnsteps >= 50000 && eee$rnsteps < 200000){

                    eee$warn_log[[length(eee$warn_log) + 1]]  <- i
                } else if (eee$rnsteps >= 200000 && !eee$force_exec){

                    warnMessage(warn_log = eee$warn_log, quiet = eee$quiet)
                    message(paste0(MZ1, i, MZ2))
                    out.value <- NULL
                    return(NULL)

                } else if (eee$rnsteps < 5000000 && eee$force_exec) {
                    eee$warn_log[[length(eee$warn_log) + 1]] <- i
                } else if (eee$rnsteps >= 200000) {

                    warnMessage(warn_log = eee$warn_log, quiet = eee$quiet)
                    message(paste0(MZ1, i, MZ2))
                    out.value <- NULL
                    return(NULL)
                }
            }
        }

        eee$st <- rep(0, times = eee$nnew)
        eee$fi <- rep(0, times = eee$nnew)
        eee$h <- rep(0, times = eee$nbonds)
        eee$ok <- rep(1, times = eee$nold)
        eee$nlost <- (eee$nnew - eee$nold) > 0

        for (ii in seq(1,eee$nold,by=1)) {
            eee$h[which(eee$bonds[,1] == eee$uold[ii])] <- ii
        }
        eee$st[1] <- 1
        eee$fi[eee$nnew] <- eee$nbonds; ## check this later
        if (eee$nnew > 1){
            eee$sb <- eee$bonds[, 2]
            eee$sbr <- circshift(eee$sb,1)
            eee$sbl <- circshift(eee$sb,-1)
            eee$st[seq(2,length(eee$st),by=1)] <- which(
                eee$sb[seq(2,length(eee$sb),by=1)] !=
                    eee$sbr[seq(2,length(eee$sbr),by=1)]) + 1
            eee$fi[seq(1,(eee$nnew-1),by=1)] <- which(
                eee$sb[seq(1,(eee$nbonds-1),by=1)] !=
                    eee$sbl[seq(1,(eee$nbonds-1),by=1)]
            )
        }

        eee$checkflag <- 0
        while (eee$checkflag != 2){

            eee$pt <- eee$st - 1
            eee$lost <- matrix(0, nrow = eee$nnew, ncol = 1)
            eee$who <- 0
            eee$losttot <- 0
            eee$mndisq <- eee$nnew*eee$maxdisq

            while (eee$who != (-1)){
                internalPermutation(xyzs=xyzs, maxdisp=maxdisp, env=eee)
            }
            eee$checkflag <- eee$checkflag + 1
            if (eee$checkflag == 1){
                eee$plost <- min(
                    c(matfix(eee$mndisq/eee$maxdisq), (eee$nnew -1)))
                if (eee$plost > eee$nlost){
                    eee$nlost <- eee$plost
                } else {
                    eee$checkflag <- 2
                }
            }
        }
        eee$resx[eee$ispan, eee$labely[eee$bonds[eee$minbonds, 2]]] <-
            eee$eyes[eee$labelx[(eee$bonds[eee$minbonds,1] + 1)]]
        eee$found[eee$labelx[(eee$bonds[eee$minbonds,1] + 1)], 1] <- 1

    }, error = function(e) { message(e); message(MZ3) })
    return(out.value)
}





#' Subnetwork Tracking
#'
#' Perform Internal Subnetwork Tracking as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned while objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::subNetworkTracking(data.frame(1), 1, new.env())
#'
#' @keywords internal
subNetworkTracking <- function(xyzs, maxdisp, env){

    stopifnot(is.numeric(maxdisp), is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the subNetworkTracking function"

    tryCatch({
        suppressWarnings({

            # THE SUBNETWORK CODE BEGINS
            eee$lista <- matrix(0, nrow = eee$numbonds, ncol = 1)
            eee$listb <- matrix(0, nrow = eee$numbonds, ncol = 1)
            eee$nclust <- 0
            eee$maxsz <- 0
            eee$thru <- eee$xdim

            while (eee$thru != 0) {
                #%the following code extracts connected
                #% sub-networks of the non-trivial
                #% bonds. NB: lista/b can have redundant entries due to
                #% multiple-connected subnetworks
                eee$w <- which(eee$bonds[, 2] >= 0)
                eee$lista[1] = eee$bonds[eee$w[1],2]
                eee$listb[1] = eee$bonds[eee$w[1],1]
                eee$bonds[eee$w[1], ] <- (-1) * (eee$nclust+1)
                # bonds;
                eee$adda <- 1
                eee$addb <- 1
                eee$donea <- 0
                eee$doneb <- 0
                if ((eee$donea != eee$adda) || (eee$doneb != eee$addb)){
                    eee$true <- FALSE
                } else {
                    eee$true <- TRUE
                }

                while (!eee$true){
                    if (eee$donea != eee$adda) {
                        eee$w <- which(eee$bonds[,2] == eee$lista[eee$donea+1])
                        eee$ngood <- length(eee$w)
                        if (eee$ngood != 0) {
                            mSQ1 <- seq((eee$addb+1),(eee$addb+eee$ngood),by=1)
                            eee$listb[mSQ1, 1] <- eee$bonds[eee$w,1]
                            eee$bonds[eee$w,] <- (-1)*(eee$nclust+1)
                            eee$addb <- eee$addb+eee$ngood;
                        }
                        eee$donea <- eee$donea + 1
                    }
                    if (eee$doneb != eee$addb){
                        eee$w <- which(eee$bonds[,1] == eee$listb[eee$doneb+1])
                        eee$ngood <- length(eee$w);
                        if (eee$ngood != 0) {
                            mSQ2 <- seq((eee$adda+1),(eee$adda+eee$ngood),by=1)
                            eee$lista[mSQ2,1] <- eee$bonds[eee$w,2]
                            eee$bonds[eee$w,] <- (-1)*(eee$nclust+1)
                            eee$adda <- eee$adda+eee$ngood;
                        }
                        eee$doneb <- eee$doneb + 1
                    }
                    if ((eee$donea != eee$adda) || (eee$doneb != eee$addb)){
                        eee$true <- FALSE
                    } else {
                        eee$true = TRUE
                    }
                }

                eee$pp <- sort(eee$listb[seq(1,eee$doneb,by=1)])
                eee$pqx <- order(eee$listb[seq(1,eee$doneb,by=1)])
                #%unx =    unq(listb(1:doneb),pqx);
                #%implanting unq directly
                eee$arr <- eee$listb[seq(1,eee$doneb,by=1)]
                eee$q <- eee$arr[eee$pqx]
                eee$indices <- which(eee$q != circshift(eee$q,-1))
                eee$count <- length(eee$indices)
                if (eee$count > 0){
                    eee$unx <- eee$pqx[eee$indices]
                } else {
                    eee$unx <- length(eee$q) -1
                }

                eee$xsz <- length(eee$unx)
                eee$pp <- sort(eee$lista[seq(1,eee$donea,by=1)])
                eee$pqy <- order(eee$lista[seq(1,eee$donea,by=1)])

                #%implanting unq directly
                eee$arr <- eee$lista[seq(1,eee$donea,by=1)]
                eee$q <- eee$arr[eee$pqy]
                eee$indices <- which(eee$q != circshift(eee$q,-1))
                eee$count <- length(eee$indices)
                if (eee$count > 0){
                    eee$uny <- eee$pqy[eee$indices]
                } else {
                    eee$uny <- length(eee$q) -1
                }

                eee$ysz <- length(eee$uny)
                if ((eee$xsz*eee$ysz) > eee$maxsz){
                    eee$maxsz <- eee$xsz*eee$ysz
                    eee$mxsz <- eee$xsz
                    eee$mysz <- eee$ysz
                }

                eee$thru <- eee$thru - eee$xsz
                eee$nclust <- eee$nclust + 1
            }
            eee$bmap <- eee$bonds[,2]
        })
    }, error = function(e) {message(MZ3)})
    return(NULL)
}


#' Trivial Bond Tracking
#'
#' Perform Trivial Bond Tracking as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned when particles are correctly tracked; TRUE is
#' returned when no particles are found to be tracked; objects in env are
#' updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::trivialBondTracking(data.frame(1), new.env())
#'
#' @keywords internal
trivialBondTracking <- function(xyzs, env) {

    stopifnot(is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the trivialBondTracking function"

    tryCatch({
        suppressWarnings({
            #  THE TRIVIAL BOND {else} block CODE BEGINS
            #  or: Use simple N^2 time routine to calculate trivial bonds
            eee$wh <- which(eee$pos[,1] >= 0)
            eee$ntrack <- length(eee$wh)
            if (eee$ntrack == 0){
                message('There are no valid particles to track!')
                return(TRUE)
            }

            # yma initialization was added
            eee$xmat <- matrix(0, nrow = eee$ntrack, ncol = eee$m)
            eee$ymat <- matrix(0, nrow = eee$ntrack, ncol = eee$m)

            eee$count <- 0
            for (kk in seq(1,eee$ntrack,by=1)) {
                for (ll in seq(1,eee$m,by=1)) {
                    eee$xmat[kk,ll] <- eee$count
                    eee$count <- eee$count+1
                }
            }
            eee$count <- 0

            # if there are not enough cols or rows, add them and set to 0
            if (nrow(eee$ymat) < eee$m) {
                eee$TMP <- matrix(
                    0, nrow = (eee$m - nrow(eee$ymat)), ncol = ncol(eee$ymat))
                eee$ymat <- rbind(eee$ymat, eee$TMP)
            }
            if (ncol(eee$ymat) < eee$ntrack) {
                eee$TMP <- matrix(
                    0, nrow = nrow(eee$ymat),
                    ncol = (eee$ntrack - ncol(eee$ymat))
                )
                eee$ymat <- cbind(eee$ymat, eee$TMP)
            }

            for (kk in seq(1,eee$m,by=1)) {
                for (ll in seq(1,eee$ntrack,by=1)) {
                    eee$ymat[kk,ll] <- eee$count
                    eee$count <- eee$count+1
                }
            }

            eee$xmat <- (eee$xmat %% eee$m) + 1
            eee$ymat <- t((eee$ymat %% eee$ntrack) +1)
            eee$lenxn <- nrow(eee$xmat)
            eee$lenxm <- ncol(eee$xmat)

            for (d in seq_len(eee$dim)) {
                eee$x <- eee$xyi[,d]
                eee$y <- eee$pos[eee$wh,d]

                eee$xm <- vapply(seq_len(ncol(eee$xmat)), function(jj) {
                    tcljj <- eee$xmat[, jj]
                    eee$x[tcljj]
                }, FUN.VALUE = numeric(nrow(eee$xmat)))

                #ym <- y[ymat[1:lenxn, 1:lenxm]]
                eee$tmpymat <- eee$ymat[seq_len(eee$lenxn), seq_len(eee$lenxm)]
                eee$ym <- vapply(seq_len(ncol(eee$tmpymat)), function(jj) {
                    tcljj <- eee$tmpymat[, jj]
                    eee$y[tcljj]
                }, FUN.VALUE = numeric(nrow(eee$tmpymat)))

                if (nrow(eee$xm) != nrow(eee$ym) ||
                    ncol(eee$xm) != ncol(eee$ym)) {
                    eee$xm <- t(eee$xm)
                }

                if (d == 1) {
                    eee$dq <- (eee$xm - eee$ym)^2
                    #%dq = (x(xmat)-y(ymat(1:lenxn,1:lenxm))).^2;
                } else {
                    eee$dq <- eee$dq + (eee$xm-eee$ym)^2
                    #%dq = dq + (x(xmat)-y(ymat(1:lenxn,1:lenxm)) ).^2;
                }
            }

            eee$ltmax <- 1 * (eee$dq < eee$maxdisq)

            #% figure out which trivial bonds go with which
            eee$rowtot <- matrix(0, nrow = eee$n, ncol = 1)
            eee$rowtot[eee$wh, 1] <- apply(eee$ltmax, 1, sum)

            if (eee$ntrack > 1) {
                eee$coltot <- apply(eee$ltmax, 2, sum, na.rm = TRUE)
            } else {
                eee$coltot <- eee$ltmax
            }
            eee$which1 <- matrix(0, nrow = eee$n, ncol = 1)
            for (j in seq(1,eee$ntrack,by=1)) {
                eee$mx    <- max(eee$ltmax[j, ], na.rm = TRUE)
                eee$w <- which.max(eee$ltmax[j, ])
                eee$which1[eee$wh[j]] <- eee$w
            }

            eee$ntrk <- matfix( eee$n - sum(eee$rowtot == 0))
            eee$w <- which( eee$rowtot == 1)
            eee$ngood <- length(eee$w)
            if (eee$ngood != 0) {
                eee$ww <- which(eee$coltot[eee$which1[eee$w]] == 1)
                eee$ngood <- length(eee$ww)
                if (eee$ngood != 0) {
                    eee$resx[ eee$ispan, eee$w[eee$ww] ] <-
                        eee$eyes[ eee$which1[eee$w[eee$ww]]]
                    eee$found[eee$which1[ eee$w[eee$ww]]] <- 1
                    eee$rowtot[eee$w[eee$ww]] <- 0
                    eee$coltot[eee$which1[eee$w[eee$ww]]] <- 0
                }
            }

            eee$labely <- which(eee$rowtot > 0)
            eee$ngood <- length(eee$labely)

            if (eee$ngood != 0) {
                eee$labelx <- which(eee$coltot > 0)
                eee$nontrivial <- 1
            } else {
                eee$nontrivial <- 0
            }
        })
    }, error = function(e) {message(MZ3)})
    return(FALSE)
}


#' Inner Bond Raster
#'
#' Perform Inner Bond Raster as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param i numeric, index of the iteration cycle
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned; objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::innerBondRaster(data.frame(1), 1, 1, new.env())
#'
#'
#' @keywords internal
innerBondRaster  <- function(xyzs, maxdisp, j=1, env) {

    stopifnot(
        is.numeric(maxdisp), is.numeric(j),
        is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the innerBondRaster function"

    tryCatch({
        suppressWarnings({

            eee$map <- matfix(-1)
            eee$scub_spos <- eee$scube + eee$spos[j];
            eee$s <- eee$scub_spos %% eee$nblocks
            eee$whzero <- which(eee$s == 0 )
            if (length(eee$whzero > 0)){
                eee$nfk <- which(eee$s !=0 )
                eee$s <- eee$s[eee$nfk]
            }

            eee$w <- which(eee$strt[eee$s, 1] != (-1))
            eee$ngood <- eee$length(eee$w)
            eee$ltmax <- 0
            if (eee$ngood != 0){

                eee$s <- eee$s[eee$w]
                for (k in seq(1,eee$ngood,by=1)){
                    eee$map = c(eee$map,
                                eee$isort[ seq(
                                    eee$strt[eee$s[k]],
                                    eee$fnsh[eee$s[k]],by=1) ])}
                eee$map <- eee$map[seq(2,length(eee$map),by=1)]
                eee$distq <- matrix(0, nrow=length(eee$map), ncol=1)
                for (d in seq(1,eee$dim,by=1)){
                    eee$distq <- eee$distq + (
                        eee$xyi[eee$map,d] - eee$pos[j,d])^2
                }
                eee$ltmax <- eee$distq < eee$maxdisq

                rowtot[j, 1] <- sum(eee$ltmax)

                if (eee$rowtot[j] >= 1){
                    eee$w <- which(eee$ltmax == 1)
                    eee$coltot[eee$map[eee$w], 1] <-
                        eee$coltot[ eee$map[eee$w], 1] +1
                    which1[j, 1] <- eee$map[ eee$w[1]]
                }
            }
        })
    }, error = function(e) {message(MZ3)})
    return(FALSE)
}


#' Trivial Bond Raster
#'
#' Perform Trivial Bond Raster as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned; objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::trivialBondRaster(data.frame(1), new.env())
#'
#' @keywords internal
trivialBondRaster <- function(xyzs, env) {

    stopifnot(is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the trivialBondRaster function"

    tryCatch({
        suppressWarnings({
            # %Use the raster metric code to do trivial bonds
            eee$abi <- matfix(eee$xyi/eee$blocksize)
            eee$abpos <- matfix(eee$pos/eee$blocksize)
            eee$si <- matrix(0, nrow = eee$m, ncol = 1)
            eee$spos <- matrix(0, nrow = eee$n, ncol = 1)
            eee$dimm <- matrix(0, nrow=eee$dim, ncol=1)
            eee$coff <- 1

            for (j in seq(1,eee$dim,by=1)){
                eee$minn <- min(c(
                    as.numeric(eee$abi[,j]),
                    as.numeric(eee$abpos[, j])),na.rm = TRUE)
                eee$maxx <- max(c(
                    as.numeric(eee$abi[,j]),
                    as.numeric(eee$abpos[,j])),na.rm = TRUE)
                eee$abi[, j] <- eee$abi[, j] - eee$minn
                eee$abpos[,j] <- eee$abpos[,j] - eee$minn
                eee$dimm[j,1] <- eee$maxx-eee$minn + 1
                eee$si <- eee$si + eee$abi[,j] * eee$coff
                eee$spos <- eee$spos + eee$abpos[,j]*eee$coff
                eee$coff <- eee$dimm[j,1]*eee$coff
            }

            eee$nblocks <- eee$coff
            #% trim down (intersect) the hypercube if its too big to fit
            # in the particle volume. (i.e. if dimm(j) lt 3)
            eee$cub <- eee$cube
            eee$deg <- which(eee$dimm[,1] < 3)
            if (length(eee$deg) > 0) {
                for (j in seq(0,(length(eee$deg)-1), by=1)){
                    #cub <- cub[which(cub[, deg[j+1]] < dimm[deg[j+1],1]) ,]
                    abc01 <- eee$cub[, eee$deg[j+1]]
                    abc02 <- eee$dimm[eee$deg[j+1],1]
                    eee$cub <- eee$cub[which(abc01 < abc02) ,]
                }
            }

            # % calculate the "s" coordinates of hypercube
            # (with a corner @ the origin)
            eee$scube <- matrix(0, nrow = nrow(eee$cub), ncol=1)
            eee$coff <- 1
            for (j in seq(1,dim,by=1)){
                eee$scube <- eee$scube + (eee$cub[,j] * eee$coff)
                eee$coff <- eee$coff*eee$dimm[j, 1]
            }

            # % shift the hypercube "s" coordinates to be centered
            # around the origin
            eee$coff <- 1
            for (j in seq(1,eee$dim,by=1)){
                if (eee$dimm[j, 1] > 3) { eee$scube <- eee$scube - eee$coff }
                eee$coff <- eee$dimm[j, 1] * eee$coff
            }
            eee$scube <- (eee$scube + eee$nblocks) %% eee$nblocks

            # get the sorting for the particles by their "s" positions.
            eee$ed <- sort(eee$si)
            eee$isort <- order(eee$si)

            #% make a hash table which will allow us to know which
            # new particles are at a given si.
            eee$strt <- matrix((-1), nrow = eee$nblocks, ncol = 1)
            eee$fnsh <- matrix(0, nrow = eee$nblocks, ncol = 1)
            eee$h <- which(eee$si == 0)
            eee$lh <- length(eee$h)
            if (eee$lh > 0) {
                eee$si[eee$h] <- 1
            }

            for (j in seq(1,eee$m,by=1)){
                if (eee$strt[eee$si[eee$isort[j]], 1] == (-1)){
                    eee$strt[eee$si[eee$isort[j]],1] <- j
                    eee$fnsh[eee$si[eee$isort[j]], 1] <- j
                } else {
                    eee$fnsh[eee$si[eee$isort[j]], 1] <- j
                }
            }
            if (eee$lh > 0) {
                eee$si[eee$h] <- 0
            }

            eee$coltot <- matrix(0, nrow = eee$m, ncol = 1)
            eee$rowtot <- matrix(0, nrow = eee$n, ncol = 1)
            eee$which1 <- matrix(0, nrow = eee$n, ncol = 1)
            for (j in seq(1,eee$n,by=1)){
                zzz <- innerBondRaster(j)
            }

            eee$ntrk <- matfix(eee$n - sum(eee$rowtot == 0))
            eee$w <- which(eee$rowtot == 1)
            eee$ngood <- length(eee$w)

            if (eee$ngood != 0) {
                #ww <- which(coltot( which1[w] ) == 1);
                eee$ww <- which(eee$coltot[eee$which1[eee$w]] == 1)

                eee$ngood <- length(eee$ww)
                if (eee$ngood != 0){
                    # %disp(size(w(ww)))
                    eee$resx[eee$ispan, eee$w[eee$ww]] <-
                        eee$eyes[eee$which1[eee$w[eee$ww]]]
                    eee$found[eee$which1[eee$w[eee$ww]]] <- 1
                    eee$rowtot[eee$w[eee$ww]] = 0;
                    eee$coltot[eee$which1[eee$w[eee$ww]]] <- 0
                }
            }

            eee$labely <- which(eee$rowtot > 0)
            eee$ngood <- length(eee$labely)
            if (eee$ngood != 0){
                eee$labelx <- which(eee$coltot > 0)

                eee$nontrivial <- 1
            } else {
                eee$nontrivial <- 0
            }
        })
    }, error = function(e) {message(MZ3)})
    return(FALSE)
}


#' Non-Trivial Bond Tracking
#'
#' Perform Non-Trivial Bond Tracking as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param env environment, including all objects used for the tracking
#' @param i integer, index of the current cycle
#'
#' @return FALSE is returned; objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::nontrivialBondTracking(data.frame(1), 1, 1, new.env())
#'
#' @keywords internal
nontrivialBondTracking <- function(xyzs, maxdisp, i, env){

    stopifnot(
        is.numeric(maxdisp), is.data.frame(xyzs), is.environment(env))

    # link to environment
    ret.out <- FALSE
    eee <- env
    MZ3 <- "An issue may have occurred with the nontrivialBondTracking function"

    CHK1 <- tryCatch({
        suppressWarnings({
            eee$xdim <- length(eee$labelx)
            eee$ydim <- length(eee$labely)

            # make a list of the non-trivial bonds
            eee$bonds <- list()
            eee$bondlen <- list()

            for (j in seq(1,eee$ydim,by=1)) {
                eee$distq <- matrix(0, nrow = eee$xdim, ncol = 1)

                for (d in seq(1,eee$dim,by=1)) {
                    #%distq
                    eee$distq <- eee$distq + cbind(
                        (eee$xyi[eee$labelx,d] - eee$pos[eee$labely[j],d])^2
                    )
                }

                eee$w <- which(eee$distq < eee$maxdisq) - 1
                eee$ngood <- length(eee$w)
                eee$newb <- rbind(eee$w, rep(j, times = eee$ngood))

                eee$bonds[[(length(eee$bonds) + 1)]] <- t(eee$newb)
                eee$bondlen[[(length(eee$bondlen) + 1)]] <-
                    eee$distq[eee$w + 1]
            }

            eee$bonds <- do.call(rbind, eee$bonds)
            eee$bondlen <- do.call(c, eee$bondlen)
            eee$numbonds <- length(eee$bonds[,1])
            eee$mbonds <- eee$bonds;

            if (max(c(eee$xdim,eee$ydim)) < 4){
                eee$nclust <- 1
                eee$maxsz <- 0
                eee$mxsz <- eee$xdim
                eee$mysz <- eee$ydim
                eee$bmap <- matrix((-1), nrow = length(eee$bonds[,1])+1, 1)

            } else {
                subNetworkTracking(xyzs, maxdisp, env=eee)
            }

            ## Adjusting nclust
            eee$all_clusts <- unique(abs(eee$bmap))
            eee$nclust <- length(eee$all_clusts)

            # PERMUTATION CODE
            perm.sum <- list()
            for (nc in seq(1,eee$nclust,by=1)){
                perm.sum[[nc]] <-
                    runTrackingPermutation(
                        xyzs, maxdisp = maxdisp, nc = nc, i = i, env = eee)
            }
            if (length(perm.sum) != eee$nclust) {
                ret.out <- TRUE
            } else if (sum(do.call(c, perm.sum)) != 0) {
                ret.out <- TRUE
            }

        })
        0
    }, error = function(e) {message(MZ3); 1})

    if (CHK1 != 0) ret.out <- TRUE
    return(ret.out)
}


#' Tracking Slide Wrap Up
#'
#' Perform Tracking Slide Wrap Up as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param env environment, including all objects used for the tracking
#' @param i integer, index of the current cycle
#'
#' @return FALSE is returned; objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#'
#' @examples
#' cellmigRation:::trackSlideWrapUp(data.frame(1), 1, 1, new.env())
#'
#' @keywords internal
trackSlideWrapUp <- function(xyzs, maxdisp, i, env){

    stopifnot(
        is.numeric(maxdisp), is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the trackSlideWrapUp function"

    tryCatch({
        suppressWarnings({

            eee$nold <- length(eee$bigresx[1,])
            eee$nnew <- eee$n - eee$nold;
            if (eee$nnew > 0){

                eee$newarr <- matrix(-1, nrow = eee$z, ncol = eee$nnew)
                ## bigresx <- c(bigresx, newarr)
                eee$bigresx <- cbind(eee$bigresx, eee$newarr)
            }
            if (eee$goodenough > 0){
                if ((sum(eee$dumphash)) > 0){
                    eee$wkeep <- which(eee$dumphash == 0)
                    eee$nkeep <- length(eee$wkeep)
                    eee$resx <- eee$resx[ ,eee$wkeep]
                    eee$bigresx <- eee$bigresx[, eee$wkeep]
                    eee$pos <- eee$pos[eee$wkeep, ]
                    eee$mem <- eee$mem[eee$wkeep]
                    eee$uniqid <- eee$uniqid[eee$wkeep]
                    eee$nvalid <- eee$nvalid[eee$wkeep]
                    eee$n <- eee$nkeep
                    eee$dumphash <- matrix(0, nrow = eee$nkeep, ncol = 1)
                }
            }

            #% again goodenough keyword
            if (!eee$quiet) {
                MZ1 <- paste(
                    i, 'of' , eee$z, 'done. Tracking', eee$ntrk,
                    'particles.', eee$n, 'tracks total.')
                message(MZ1)
            }

            if (!is.matrix(eee$bigresx) || nrow(eee$resx) > nrow(eee$bigresx)){
                eee$bigresx <- rbind(eee$bigresx)
                eee$bigresx <- rbind(
                    eee$bigresx,
                    matrix(
                        -1, nrow = (nrow(eee$resx) - nrow(eee$bigresx)),
                        ncol = ncol(eee$bigresx)))
            }

            eee$bigresx[seq((i-(eee$ispan)+1),i,by=1),] <-
                eee$resx[seq(1,eee$ispan,by=1),]
            eee$resx <- matrix((-1), nrow = eee$zspan, ncol = eee$n)

            eee$wpull <- which(eee$pos[ ,1] == (-1 * maxdisp))
            eee$npull <- length(eee$wpull)

            if (eee$npull > 0){
                eee$lillist <- list()
                for (ipull in seq(1,eee$npull,by=1)){
                    eee$wpull2 <- which(eee$bigresx[, eee$wpull[ipull]] != (-1))
                    eee$npull2 <- length(eee$wpull2)

                    eee$thing = cbind(
                        eee$bigresx[eee$wpull2,eee$wpull[ipull]],
                        rep(x = eee$uniqid[eee$wpull[ipull]],
                            times = eee$npull2)
                    )
                    eee$lillist[[length(eee$lillist) + 1]] <- eee$thing
                }
                eee$olist[[length(eee$olist) + 1]] <-
                    do.call(rbind, eee$lillist)
            }

            eee$wkeep <- which(eee$pos[, 1] >= 0)
            eee$nkeep <- length(eee$wkeep)
            if (eee$nkeep == 0) {
                message ('Were going to crash now, no particles....')
            }
            eee$resx <- eee$resx[,eee$wkeep]
            eee$bigresx <- eee$bigresx[, eee$wkeep]
            eee$pos <- eee$pos[eee$wkeep, ]
            eee$mem <- eee$mem[eee$wkeep]
            eee$uniqid <- eee$uniqid[eee$wkeep]
            eee$n <- eee$nkeep;
            eee$dumphash <- matrix(0, nrow = eee$nkeep, ncol =1)
            if (eee$goodenough > 0){
                eee$nvalid <- eee$nvalid[eee$wkeep]
            }
        })
    }, error = function(e) {message(MZ3)})
    return(FALSE)
}


#' Tracking Slide Processing
#'
#' Frame-by-frame Slide Processing as part of the Cell Tracking Processing
#'
#' @param xyzs data.frame, including input cell centroid positions
#' @param maxdisp numeric, value of maximum cell dispersion in pixels
#' @param env environment, including all objects used for the tracking
#'
#' @return FALSE is returned as long as cells are detected and tracked;
#' TRUE is returned when no further cells to track are found;
#' objects in env are updated
#'
#' @details a message is printed if an issue (typically arising by a
#' non-suitable environment being passed as the env argument)  is detected.
#' See the example below.
#'
#' @examples
#' cellmigRation:::trackSlideProcessing(data.frame(1), 1, 1, new.env())
#'
#' @keywords internal
trackSlideProcessing <- function(xyzs, maxdisp, i, env) {

    stopifnot(
        is.numeric(i), is.numeric(maxdisp),
        is.data.frame(xyzs), is.environment(env))

    # link to environment
    eee <- env
    MZ3 <- "An issue may have occurred with the trackSlideProcessing function"
    ret.out <- FALSE

    tryCatch({
        suppressWarnings({
            eee$ispan <- ((i-1) %% eee$zspan) + 1
            eee$m <- eee$res[(i+1)] - eee$res[(i)]
            eee$eyes <- seq(1, eee$m, by=1)
            eee$eyes <- eee$eyes + eee$res[i]

            if (eee$m > 0) {
                eee$xyi <- xyzs[eee$eyes, seq(1,eee$dim,by=1)]
                eee$found <- matrix(0, nrow = eee$m, ncol = 1)

                # % THE TRIVIAL BOND CODE BEGINS
                if (eee$notnsqrd) {
                    zzzz <- trivialBondRaster(xyzs, env = eee)

                } else {
                    # TRIVIAL BOND CODE - ELSE
                    tbtrk <- trivialBondTracking(xyzs, env = eee)
                    if (tbtrk) { ret.out <- TRUE; return(TRUE) }
                }

                if (eee$nontrivial == 1){
                    ntbtrk <- nontrivialBondTracking(xyzs, maxdisp, i, eee)
                    if (ntbtrk) { ret.out <- TRUE }
                }
                eee$w <- which(eee$resx[eee$ispan,] >= 0)
                eee$nww <- length(eee$w)

                if (eee$nww > 0){
                    eee$pos[eee$w,] <-
                        xyzs[eee$resx[eee$ispan,eee$w], seq(1,eee$dim,by=1)]
                    if (eee$goodenough > 0){
                        eee$nvalid[eee$w] <- eee$nvalid[eee$w] + 1
                    }
                }
                eee$newguys <- which(eee$found == 0)
                eee$nnew <- length(eee$newguys)

                if (eee$nnew > 0) { ##% & another keyword to workout inipos
                    eee$newarr <- matrix(-1, nrow = eee$zspan, ncol = eee$nnew)
                    eee$resx <- cbind(eee$resx, eee$newarr)
                    eee$resx[eee$ispan, (seq((eee$n+1),ncol(eee$resx),by=1))]<-
                        eee$eyes[eee$newguys]
                    eee$pos <- rbind(eee$pos, xyzs[
                        eee$eyes[eee$newguys], (seq(1,eee$dim,by=1))])
                    eee$nmem <- matrix(0, nrow = eee$nnew, ncol = 1)
                    eee$mem <- c(eee$mem, eee$nmem)
                    eee$nun <- seq(1,eee$nnew,by=1)
                    eee$uniqid <- c(eee$uniqid, ((eee$nun) + eee$maxid))
                    eee$maxid <- eee$maxid + eee$nnew
                    if (eee$goodenough > 0){
                        eee$dumphash <- c(
                            eee$dumphash, t(matrix(0,nrow=1,ncol=eee$nnew)))
                        eee$nvalid <- c(
                            eee$nvalid, t(matrix(1,nrow=1,ncol=eee$nnew)))
                    }
                    eee$n <- eee$n + eee$nnew
                }

            } else {
                #' Warning- No positions found for t='
                message("@@", appendLF = FALSE)
            }

            eee$w <- which(eee$resx[eee$ispan,] != (-1))
            eee$nok <- length(eee$w)
            if (eee$nok != 0){
                eee$mem[eee$w] <- 0
            }

            eee$mem <- eee$mem + (0 + (cbind(eee$resx[eee$ispan,]) == -1))
            eee$wlost <- which(eee$mem == eee$memory_b+1)
            eee$nlost <- length(eee$wlost)

            if (eee$nlost > 0){
                eee$pos[eee$wlost, ] <- -(maxdisp)
                if (eee$goodenough > 0){
                    eee$wdump <- which(eee$nvalid[eee$wlost] < eee$goodenough)
                    eee$ndump <- length(eee$wdump);
                    if (eee$ndump > 0){
                        eee$dumphash[eee$wlost[eee$wdump]] <- 1
                    }
                }
            }

            if ((eee$ispan == eee$zspan) | (i == eee$z)){
                zzz <- trackSlideWrapUp(xyzs, maxdisp, i, env = eee)
            }
        })

    }, error = function(e) {
        message(MZ3)
    })
    return(ret.out)
}


#' Track cells
#'
#' Constructs n-dimensional trajectories from a scrambled list of
#' particle coordinates determined at discrete times (e.g. in consecutive
#' image frames)
#'
#'
#' @param xyzs an array listing the xy coordinates and data of
#' the different particles at different times
#' @param maxdisp an estimate of the maximum distance that a particle
#' would move in a single time interval
#' @param params a list containing a few tracking parameters that are
#' needed for the analysis
#'
#'
#' @return data.frame including cell tracks data
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
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
    # Initialize and stuff
    eee <- new.env(parent = emptyenv())
    eee$warn_log <- list()
    eee$dd <- ncol(xyzs)
    MZ4 <- "An issue may have occurred with the track function"
    errTrack <- FALSE
    tryCatch({
        suppressWarnings({

            eee$tmp.param <- initializeTrackParams(eee$dd, params)
            eee$memory_b <- eee$tmp.param$memory_b;
            eee$goodenough <- eee$tmp.param$goodenough
            eee$dim <- eee$tmp.param$dim;
            eee$quiet <- eee$tmp.param$quiet
            eee$force_exec <- eee$tmp.param$force_exec

            # Checking the input time vector
            # This should be monotonically not-decreasing and not identical
            eee$tau <- xyzs[, eee$dd]
            eee$st <- eee$tau[seq(2,length(eee$tau), by = 1)] -
                eee$tau[seq(1,(length(eee$tau) - 1), by = 1)]

            if (sum(eee$st < 0) > 0) {
                message("", appendLF = TRUE)
                message("The time vector (tau) is not ordered")
                return(NULL)
            }
            if (length(unique(eee$tau)) == 1) {
                message("", appendLF = TRUE)
                message('All positions are at the same time... go back!')
                return(NULL)
            }

            eee$info <- 1
            eee$w <- which(eee$st > 0)
            eee$z <- length(eee$w)
            eee$z <- eee$z + 1

            # implanting unq directly
            eee$indices <- which(eee$tau - circshift(eee$tau, -1) != 0)
            eee$count <- length(eee$indices)

            if (eee$count > 0) {
                eee$res <- eee$indices
            } else{
                eee$res = length(eee$tau)-1
            }

            eee$res <- c(1, eee$res, length(eee$tau))
            eee$ngood <- eee$res[2] - eee$res[1] + 1
            eee$eyes <- seq(1,eee$ngood,by=1)
            eee$pos <- xyzs[eee$eyes, seq(1,eee$dim, by = 1)]
            eee$istart <- 2
            eee$n <- eee$ngood;
            eee$zspan <- 50;
            if (eee$n > 200)
                eee$zspan <- 20
            if (eee$n > 500)
                eee$zspan <- 10

            # initialize 2 matrices with -1
            eee$resx <- matrix((-1), nrow = eee$zspan, ncol = eee$n)
            eee$bigresx <- matrix((-1), nrow = eee$z, ncol = eee$n)
            eee$mem <- matrix(0, nrow = eee$n, ncol = 1)
            eee$uniqid <- seq(1,eee$n,by=1);
            eee$maxid <- eee$n;

            # initialize olis
            eee$olist <- list()

            if (eee$goodenough > 0) {
                eee$dumphash <- matrix(0, nrow = eee$n, ncol = 1)
                eee$nvalid <- matrix(1, nrow = eee$n, ncol = 1)
            }
            eee$resx[1,] <- eee$eyes

            #% setting up constants
            eee$maxdisq <- maxdisp^2
            eee$notnsqrd <- (sqrt(eee$n*eee$ngood) > 200) && (eee$dim < 7)

            if (eee$notnsqrd)
                trackHypercubeBuild(xyzs, eee)

            # Start the main loop over the frames.
            for (i in seq(eee$istart,eee$z,by=1)){
                tspck <- trackSlideProcessing(
                    xyzs = xyzs, maxdisp = maxdisp, i = i, env = eee)
                if (tspck) { errTrack <- TRUE; break }
            }

            if (eee$goodenough > 0){
                eee$nvalid <- apply(eee$bigresx >= 0 , 2, sum)
                eee$wkeep <- which(eee$nvalid >= eee$goodenough)
                eee$nkeep <- length(eee$wkeep)
                if (eee$nkeep == 0){
                    message(MZ4)
                    return()
                }
                if (eee$nkeep < eee$n){
                    eee$bigresx <- eee$bigresx[, eee$wkeep]
                    eee$n <- eee$nkeep
                    eee$uniqid <- eee$uniqid[eee$wkeep]
                    eee$pos <- eee$pos[eee$wkeep, ]
                }
            }

            eee$wpull <- which(eee$pos[, 1] != ((-2) * maxdisp))
            eee$npull <- length(eee$wpull);
            if (eee$npull > 0) {
                eee$lillist <- list()
                for (ipull in seq(1,eee$npull,by=1)){
                    eee$wpull2 <- which(eee$bigresx[, eee$wpull[ipull]] != (-1))
                    eee$npull2 <- length(eee$wpull2)
                    eee$thing <- cbind(
                        eee$bigresx[eee$wpull2, eee$wpull[ipull]],
                        rep(eee$uniqid[eee$wpull[ipull]], times = eee$npull2)
                    )
                    eee$lillist[[length(eee$lillist) + 1]] <- eee$thing
                }
                eee$olist[[length(eee$olist)+1]] <- do.call(rbind, eee$lillist)
            }
            eee$olist <- do.call(rbind, eee$olist)
            eee$nolist <- nrow(eee$olist)
            eee$res <- matrix(0, nrow = eee$nolist, ncol = (eee$dd+1))

            for (j in seq(1, eee$dd, by=1)){
                eee$res[, j] <- xyzs[eee$olist[, 1], j]
            }
            eee$res[, eee$dd+1] <- eee$olist[,2]
            eee$ndat <- ncol(eee$res)
            eee$newtracks <- eee$res

            eee$indices <- which(
                eee$newtracks[, eee$ndat] !=
                    circshift(eee$newtracks[, eee$ndat], -1))
            eee$count <- length(eee$indices)
            if (eee$count > 0){
                eee$u <- eee$indices
            } else {
                eee$u <- nrow(eee$newtracks) - 1
            }

            eee$ntracks <- length(eee$u)
            eee$u <- c(0, eee$u)
            for (i in seq(2, (eee$ntracks + 1),by=1)){
                eee$newtracks[seq(
                    (eee$u[(i-1)]+1), eee$u[i], by=1), eee$ndat] = (i - 1)
            }
            warnMessage(warn_log = eee$warn_log, quiet = eee$quiet)
        })
    }, error = function(e) {
        message("An issue may have occurred with the track function")
        errTrack <- TRUE
    })

    if (errTrack)
        return(NULL)

    return(eee$newtracks)
}


#' Compute Cell Migration Statistics
#'
#' Calculate the statistics from X/Y positional data obtained from
#' cell tracks
#'
#' @param tracks data.frame with cell tracks information
#' @param interval_time integer, time interval between two successive
#' frames were taken
#' @param pixel_micron integer, image resolution, i.e. number of
#' pixels
#' per micron
#'
#' @return list of stats calculated for the cell tracks.
#' Info include variables of speed,
#' distance, euclidean displacement, persistence, angular
#' displacement,
#' yFMI, xFMI, y-displacement, x-displacement and frames
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- cbind(c(30, 35, 1, 5, 6, 7, 50, 55, 56, 58),
#'             c(29, 37, 2, 7, 4, 9, 40, 50, 59, 49),
#'             c( 1,  2, 1, 2, 3, 4,  1,  2,  3,  4),
#'             c( 1,  1, 2, 2, 2, 2,  3,  3,  3,  3))
#' cellmigRation:::MigrationStats(x0, 10, 10)
#'
#'
#' @keywords internal
MigrationStats <- function(tracks, interval_time, pixel_micron) {

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
        cumulative_displacements <- as.numeric(
            cumsum(sqrt(diff(X)^2 + diff(Y)^2))
        )

        # sum of the displacements between centroids for the given track
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

    for (j in seq(1,length(arccos),by=1)){
        if    (arccos[j] < 0 & deltaY[j] > 0) {
            theta[j] <- 2*pi - theta[j]
        }

        if (arccos[j] > 0 & deltaY[j] > 0) {
            theta[j] <- 2*pi - theta[j]
        }
    }

    #theta = theta.*((2*pi)/360);
    deltaY <-    deltaY * (-1)
    yfmi <- deltaY / distance
    xfmi <- deltaX / distance
    deltaY <- deltaY * pixel_micron
    deltaX <- deltaX * pixel_micron
    speed <- speed * pixel_micron;
    distance <- distance * pixel_micron
    euclid <- euclid * pixel_micron

    OUT <- list(
        speed = speed, distance = distance,
        frames = frames, euclid = euclid,
        persistence = persistence,
        initial = initial, final = final, yfmi = yfmi,
        xfmi = xfmi, deltaX = deltaX,
        deltaY = deltaY, kept = kept, theta = theta)
    return(OUT)
}



#' Compute Tracks Stats
#'
#' Wrapper for the MigrationStats() function. It computes statistics
#' for a
#' trackedCells object where cells have already been tracked.
#'
#' @param tc_obj a \code{trackedCells} object
#' @param time_between_frames integer, time interval between two
#' successive frames were taken
#' @param resolution_pixel_per_micron integer, image resolution,
#' i.e. number of pixels per micron
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
#' x <- get(data(TrackCellsDataset))
#' x <- ComputeTracksStats(x, time_between_frames = 10,
#'                         resolution_pixel_per_micron = 20)
#' getCellsStats(x)
#'
#'
#' @export
ComputeTracksStats <- function(
    tc_obj, time_between_frames, resolution_pixel_per_micron)
{

    if(getProcessingStatus(tc_obj)$track == 0)
        stop("You need to run CellTracker() before computing stats")

    # RETRIEVE
    my_tracks <- getCellTracks(tc_obj)

    # DO
    handles <- MigrationStats(
        tracks = my_tracks, interval_time = time_between_frames,
        pixel_micron = resolution_pixel_per_micron
    )


    sz <- length(handles$speed)
    handles$cell_number <- seq(1,sz,by=1)

    cell_stats <- data.frame(
        Cell_Number = handles$cell_number, Speed = handles$speed,
        Distance = handles$distance, Displacement = handles$euclid,
        Persistence = handles$persistence,  Degrees = handles$theta,
        YFMI = handles$yfmi, XFMI = handles$xfmi,
        Y_displacement = handles$deltaY, X_displacement = handles$deltaX,
        Frames = handles$frames, stringsAsFactors = FALSE )

    # to organize population stats
    my_colz <- c(
        "Speed", "Distance", "Displacement", "Persistence", "YFMI",
        "XFMI", "Y_displacement", "X_displacement"
    )

    my_rows <- lapply(my_colz, function(cl) {
        tmp <- cell_stats[, cl]
        data.frame(
            mean = mean(tmp, na.rm = TRUE),
            SD = sd(tmp, na.rm = TRUE),
            median = median(tmp, na.rm = TRUE),
            min = min(tmp, na.rm = TRUE),
            max = max(tmp, na.rm = TRUE)
        )
    })
    my_rows <- do.call(rbind, my_rows)
    rownames(my_rows) <- my_colz

    # compute sum of cos and sin of angles
    r <- sum(exp(1i*handles$theta))

    # obtain mean angle
    meanTheta <- Arg(r)
    degrees <- meanTheta/pi*180
    my_rows <- rbind(
        my_rows,
        Angle = data.frame(
            mean = sz, SD = degrees, median = NA, min = NA, max = NA
        )
    )
    rownames(my_rows)[c(2,3)] <- c(
        "Total_displacement",
        "Euclidean_displacement"
    )
    pop_stats <- my_rows
    tc_obj <- setProcessingStatus(tc_obj, slot = "stats", value = 1)
    tc_obj <- setTrackingStats(tc_obj, pop_stats, cell_stats)

    return(tc_obj)
}


#' Generate All Combinations
#'
#' Generate All Combinations as part of the Optimization Parameter process
#'
#' @param ... a series of arguments where each argument is a vector of values
#' to be combined together
#'
#' @return a data frame of combined parameters to be tested
#'
#' @details This is an internal function supporting the Optimization Parameter
#' steps
#'
#'
#' @examples
#' cellmigRation:::GenAllCombos(A=c(1,2,3), B = 10, C = c("x", "y", "z"))
#'
#' @keywords internal
GenAllCombos <- function(...){
    xx <- list(...)
    zz <- names(xx)

    # Init
    out <- data.frame(xx[[1]], stringsAsFactors = FALSE)
    colnames(out) <- zz[1]

    # Keep attaching
    for (j in seq(2,length(xx),by=1)) {
        TMP <- xx[[j]]
        nuOUT <- list()
        for (z in TMP){
            tmpout <- cbind(
                out, data.frame(TMP = z, stringsAsFactors = FALSE)
            )
            nuOUT[[length(nuOUT) + 1]] <- tmpout
        }
        out <- do.call(rbind, nuOUT)
        colnames(out)[j] <- zz[j]
    }
    return(out)
}



#' Prepare For Parameter Optimization
#'
#' Pre-processing as part of the Optimization Parameter process
#'
#' @param stack_img input image
#' @param lnoise_range numeric, lnoise values to test
#' @param min.px.diam numeric, minimum number of pixels of a cell signal
#' @param diameter_range numeric, numeric values for diameter
#' @param threshold_range numeric, numeric values for threshold
#' @param target_cell_num numeric, target number of cells, can be NULL
#' @param quantile.val numeric, the quantile prob used to suppress noise
#' @param px.margin numeric, the frame margin size
#' @param plot logical, shall a plot be printed
#' @param verbose logical, shall info be printed to console
#'
#' @return a data frame of combined parameters to be tested
#'
#' @details This is an internal function supporting the Optimization Parameter
#' steps
#'
#' @importFrom stats quantile
#'
#' @examples
#' x <- get(data(TrackCellsDataset))
#' x <- cellmigRation::getCellImages(x = x)
#' y <- cellmigRation:::Prep4OptimizeParams(stack_img = x)
#' y$success
#'
#' @keywords internal
Prep4OptimizeParams <- function(
    stack_img, lnoise_range = NULL, min.px.diam = 5,
    diameter_range = NULL, threshold_range = NULL,
    target_cell_num = NULL, quantile.val = NULL,
    px.margin= NULL, plot=FALSE, verbose = FALSE)

{

    tryCatch({suppressWarnings({

        # select mid signal image
        imgSums <- do.call(c, lapply(stack_img$images, FUN = sum, na.rm = TRUE))
        med.i <- ifelse(
            length(imgSums) %% 2 == 0,
            length(imgSums) / 2, (0.5 * length(imgSums) + 0.5))
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
            EstimateDiameterRange(
                x = tmp_img, px.margin = px.margin, min.px.diam = min.px.diam,
                quantile.val = quantile.val, plot = FALSE)},
            error = function(e) NULL)

        # diam range
        if(is.null(diameter_range) && !is.null(estRDI)) {
            diameter_range <- c(
                floor(estRDI$q75.diam - 1),
                ceiling(1.25 * as.numeric(estRDI$q95.diam)) )
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
            lnoise_range <- unique(
                as.integer(seq(
                    from = min.px.diam,
                    to = quantile(estRDI$raw$LEN, probs = 0.25),
                    length.out = 3)))

        if(is.null(threshold_range)) {
            threshold_range <-
                seq(max(0, (min(
                    tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)],
                    na.rm = TRUE) - 1)),
                    1 + quantile(
                        tmp_img[tmp_img > min(tmp_img, na.rm = TRUE)],
                        probs = 0.75), length.out = 4)
            threshold_range <- unique(as.integer(threshold_range))
        }

        # Al params
        all_params <- GenAllCombos(
            image.i = med.i, lnoise = lnoise_range, diameter = diameter_range,
            threshold = threshold_range)
        all_params <- all_params[all_params$diameter>(4 + all_params$lnoise), ]
        rownames(all_params) <- NULL

        if(nrow(all_params) < 3) {
            message("There is a problem with the params that were submitted")
            message("A minimum number of 3 combinations are expected!")
            message("Please, try again with different param ranges")
            return(list(success = FALSE, tmp_img = NULL, all_params = NULL))
        }

        if(nrow(all_params) > 20) {
            all_params$TMPdiff <- all_params$diameter - all_params$lnoise
            re.ii <- order(all_params$TMPdiff, decreasing = TRUE)
            all_params <- all_params[re.ii,]
            rownames(all_params) <- NULL
            all_params <- all_params[seq(1,20,by=1),]
        }

        # Verbose
        if (verbose) {
            message(paste0(
                "Testing ", nrow(all_params),
                " combination(s) of params."), appendLF = TRUE)
            message("This may take some time.", appendLF = TRUE)
            message("Processing ", appendLF = FALSE)
        }
    })}, error = function(e) {return(FALSE)})
    return(list(
        success = TRUE, tmp_img = tmp_img, all_params = all_params,
        target_cell_num = target_cell_num))
}



#' Non Paralle Parameter Optimization
#'
#' Non Parallel Optimization as part of the Optimization Parameter process
#'
#' @param tmp_img numeric matrix, corresponding to the
#' @param all_params data.frame, including all parameter combinations to test
#' @param verbose logical, shall info be printed to console
#'
#' @return a list including test results; an empty list
#' is returned if an error is encountered.
#'
#' @details This is an internal function supporting the Optimization Parameter
#' process
#'
#' @examples
#' cellmigRation:::NonParallel4OptimizeParams(matrix(1), data.frame(1), TRUE)
#'
#' @keywords internal
NonParallel4OptimizeParams <- function(
    tmp_img = tmp_img, all_params = all_params, verbose = verbose)

{
    # Initialize collector (list)
    all_results <- list()

    tryCatch({suppressWarnings({

        for (i in seq(1,nrow(all_params),by=1)){

            # Verbose
            if (verbose)
                message(".", appendLF = FALSE)

            #VisualizeImg(tmp_img)
            b <- bpass(
                image_array = tmp_img,
                lnoise = all_params$lnoise[i],
                lobject = all_params$diameter[i],
                threshold = all_params$threshold[i]
            )
            tmpOUT <- list(img = b, count = 0)

            tryCatch({
                pk <- suppressMessages(
                    pkfnd(
                        im = b, th = all_params$threshold[i],
                        sz = NextOdd(all_params$diameter[i])))

                cnt <- suppressMessages(cntrd(
                    im = b, mx = pk, sz = NextOdd(all_params$diameter[i])))

                tmpOUT[["count"]] <- nrow(cnt)

            }, error = function(e) { NULL })
            all_results[[i]] <- tmpOUT
        }
    })}, error = function(e) {NULL})

    return(all_results)
}



#' Paralle Parameter Optimization
#'
#' Parallel Optimization as part of the Optimization Parameter process
#'
#' @param tmp_img numeric matrix, corresponding to the
#' @param all_params data.frame, including all parameter combinations to test
#' @param use.cores numeric, number of cores to use for the parallelization
#' @param verbose logical, shall info be printed to console
#'
#' @return a list including test results; an empty list
#' is returned if an error is encountered.
#'
#' @details This is an internal function supporting the Optimization
#' Parameter process
#'
#' @examples
#' cellmigRation:::Parallel4OptimizeParams(matrix(1), data.frame(1), 1, TRUE)
#'
#' @importFrom parallel makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @import foreach
#'
#'
#'
#' @keywords internal
Parallel4OptimizeParams <- function(
    tmp_img = tmp_img, all_params = all_params,
    use.cores = use.cores, verbose = verbose)
{
    all_results <- list()
    j <- NULL

    if (use.cores == 1) {
        if (verbose) {
            message(paste0(
                "if a single thread is used, please run ",
                "the NonParallel4OptimizeParams() function."))
        }
        return(list())
    }

    tryCatch({suppressWarnings({

        if (verbose) {
            cl <- suppressMessages(
                parallel::makeCluster(use.cores, outfile = "")
            )
        } else {
            cl <- suppressMessages(parallel::makeCluster(use.cores))
        }

        suppressMessages(doParallel::registerDoParallel(cl))

        # Nothing to export! "tmp_img", "all_params" automatically exported
        #stuffToExp <- c("tmp_img", "all_params")
        stuffToExp <- c()
        suppressMessages(parallel::clusterExport(cl, stuffToExp))

        ## %dopar%
        all_results <- tryCatch(
            foreach::foreach(
                j = seq(1,nrow(all_params),by=1),
                .verbose = verbose,
                .packages = "cellmigRation"
            ) %dopar% {
                # Verbose
                message(".", appendLF = FALSE)

                #VisualizeImg(tmp_img)
                b <- bpass(
                    image_array = tmp_img, lnoise = all_params$lnoise[j],
                    lobject = all_params$diameter[j],
                    threshold = all_params$threshold[j])
                tmpOUT <- list(img = b)

                ChkErr <- tryCatch({
                    pk <- suppressMessages(
                        pkfnd(
                            im = b, th = all_params$threshold[j],
                            sz = NextOdd(all_params$diameter[j])))
                    cnt <- suppressMessages(cntrd(
                        im = b, mx = pk, sz = NextOdd(all_params$diameter[j])))
                    tmpOUT[["count"]] <- nrow(cnt)
                    0
                }, error = function(e) { 1 })
                if (ChkErr == 1) {tmpOUT[["count"]] <- 0}
                tmpOUT
            },
            error = (function(e) {
                print(e)
                try(parallel::stopCluster(cl), silent = TRUE)
                return(NULL)
            })
        )
        message("Done!", appendLF = TRUE)

    })}, error = function(e) {NULL})

    try({suppressWarnings(parallel::stopCluster(cl))}, silent = TRUE)

    if(is.null(all_results))
        all_results <- list()

    return(all_results)
}


#' Main Loop to Parameter Optimization
#'
#' Main Loop as part of the Optimization Parameter process
#'
#' @param tmp_img numeric matrix, corresponding to the
#' @param all_params data.frame, including all parameter combinations to test
#' @param threads numeric, number of cores to use for the parallelization
#' @param verbose logical, shall info be printed to console
#'
#' @return a list including test results; an empty list
#' is returned if an error is encountered.
#'
#' @details This is an internal function supporting the Optimization
#' Parameter process
#'
#' @examples
#' cellmigRation:::OptimizeParamsMainLoop(matrix(1), data.frame(1), 1, FALSE)
#'
#' @importFrom parallel detectCores
#'
#' @keywords internal
OptimizeParamsMainLoop <- function(
    tmp_img = tmp_img, all_params = all_params,
    threads = threads, verbose = verbose)
{
    # Initialize j
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

    if (use.cores == 1) {
        all_results <-
            NonParallel4OptimizeParams(
                tmp_img = tmp_img, verbose = verbose, all_params = all_params)
    } else {

        all_results <-
            Parallel4OptimizeParams(
                tmp_img = tmp_img, all_params = all_params,
                use.cores = use.cores, verbose = verbose)
    }
    return(all_results)
}



#' Sinking Output as part of Parameter Optimization
#'
#' Tool for Sinking Output as part of the Optimization Parameter process
#'
#' @param ty character, "M" for messages, anything else for output
#' @param fl filename to sink the output to
#'
#' @return value returned by the sink() function
#'
#' @details This is an internal function supporting the Optimization
#' Parameter process
#'
#' @examples
#' print(1)
#' cellmigRation:::sinkAway("O")
#' print(2)
#' cellmigRation:::sinkAway("O", NULL)
#' print(3)
#'
#' @keywords internal
sinkAway <- function(ty = "M", fl = "/dev/null") {

    if (ty == "M") {
        tyy <- "message"
    } else {
        tyy <- "output"
    }

    z <- tryCatch(suppressWarnings(
        sink(file = fl, type = tyy)),
        error = function(e) {NULL})
    return(z)
}


#' Finalize Output in Parameter Optimization
#'
#' Finalize Output as part of the Optimization Parameter process
#'
#' @param all_results list, including all intermediates
#' @param all_params data.frame, including all parameter combinations to test
#' @param target_cell_num numeric, number of expected cells
#' @param plot logical, shall a series of plots be generated
#'
#' @return a list including test results; an empty list
#' is returned if an error is encountered.
#'
#' @details This is an internal function supporting the Optimization
#' Parameter process
#'
#' @importFrom utils head
#' @importFrom graphics par
#'
#' @examples
#' cellmigRation:::FianlizeOptiParams(list(1), data.frame(1), 10, FALSE)
#'
#' @keywords internal
FianlizeOptiParams <- function(all_results, all_params, target_cell_num, plot)
{
    OUT <- list()
    tryCatch({suppressWarnings({
        if (plot) {
            curPAR <- par(no.readonly = TRUE)
            par(mfrow = c(3, 3))
            on.exit(expr = {par(curPAR)})
        }

        all_params$counts <- do.call(c, lapply(all_results, function(x) {
            tmp <- x$count
            ifelse(is.null(tmp), 0, tmp)}))

        all_params$i <- seq(1,nrow(all_params),by=1)

        # Return Prep
        all_params$diff100 <- abs(target_cell_num - all_params$counts)
        ord_params <- all_params[order(all_params$diff100), ]
        ret.i <- utils::head(ord_params$i, n = 9)
        best_params <- list()

        top.i <- 1
        for (ri in ret.i) {

            if (top.i == 1) {
                best_params[["lnoise"]] <-
                    ord_params$lnoise[ord_params$i == ri]
                best_params[["diameter"]] <-
                    ord_params$diameter[ord_params$i == ri]
                best_params[["threshold"]] <-
                    ord_params$threshold[ord_params$i==ri]
            }

            myLAB <- paste0(
                "Pick #", top.i, "; Cell_count=",
                ord_params$counts[ord_params$i == ri], "\n")
            myLAB <- paste0(
                myLAB, "lnoise=",
                ord_params$lnoise[ord_params$i == ri],
                "; ", "diameter=",
                ord_params$diameter[ord_params$i == ri], "; ",
                "threshold=",
                ord_params$threshold[ord_params$i == ri])

            if (plot) {
                sinkAway(ty = "O", fl = NULL)
                VisualizeImg(img_mtx = all_results[[ri]]$img, main = myLAB)
            }

            top.i <- top.i + 1
        }

        OUT <- list(all_params=all_params, best_params= best_params)
    })}, error = function(e) {NULL})

    return(OUT)
}






#' Optimize Detection Params
#'
#' Optimize Detection Parameters for running a cell tracking job
#'
#' @details The lnoise param is used to guide a lowpass blurring
#' operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param
#' is used for a background correction following
#' the initial image convolution
#' \itemize{
#'
#' \item \strong{lnoise}: Characteristic lengthscale of noise in
#' pixels.
#' Additive noise averaged over this length should vanish. May
#' assume
#' any positive floating value.
#' May be also set to 0, in which case only the highpass
#' "background subtraction" operation is performed.
#'
#' \item \strong{lobject} Integer length in pixels somewhat larger
#' than a typical object.
#' Can also be set to 0, in which case only the lowpass "blurring"
#' operation defined by lnoise is done
#' without the background subtraction defined by lobject
#'
#' \item \strong{threshold} Numeric. By default, after the
#' convolution,
#' any negative pixels are reset
#' to 0.    Threshold changes the threshhold for setting pixels to 0.
#' Positive values may be useful
#' for removing stray noise or small particles.
#'
#' }
#'
#'
#'
#' @param tc_obj a \code{trackedCells} object
#' @param lnoise_range numeric vector of lnoise values to be used
#' in the optimization step. Can be NULL
#' @param min.px.diam integer, minimum diameter of a particle (cell).
#' Particles with a diameter smaller than min.px.diam are discarded
#' @param diameter_range numeric vector of diameter values to be used
#' in the optimization step. Can be NULL
#' @param threshold_range numeric vector of threshold values to be
#' used
#' in the optimization step. Can be NULL
#' @param target_cell_num integer, the expected (optimal) number of
#' cells
#' to be detected in each frame
#' @param threads integer, number of cores to use for parallelization
#' @param quantile.val numeric, argument passed to
#' EstimateDiameterRange(). If NULL, it is defaulted to 0.99
#' @param px.margin numeric, argument passed to
#' EstimateDiameterRange(). If NULL, it ia defaulted to 2
#' @param plot if `TRUE`, plots results in the end
#' @param verbose shall information about the progress of the
#' operation be printed to screen/console
#' @param dryrun shall a dryrun be performed
#'
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x <- get(data(TrackCellsDataset))
#' x <- OptimizeParams(tc_obj = x, lnoise_range = c(5,7,10),
#'                     diameter_range = c(12,14,18),
#'                     threshold_range = c(4,7), dryrun = TRUE)
#' getOptimizedParameters(x)
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
OptimizeParams <- function(
    tc_obj, lnoise_range = NULL, min.px.diam = 5,
    diameter_range = NULL, threshold_range = NULL,
    target_cell_num = NULL, threads = 1,
    quantile.val = NULL, px.margin= NULL,
    plot=FALSE, verbose = FALSE, dryrun = FALSE)

{
    stack_img <- getCellImages(tc_obj)

    if (dryrun && getProcessingStatus(tc_obj)$optimized_params == 1)
        return(tc_obj)

    if (plot) {
        curPAR <- par(no.readonly = TRUE)
        par(mfrow = c(3, 3))
        on.exit(expr = {par(curPAR)})
    }

    if (!verbose) {
        if(file.exists("/dev/null")){
            zzz <- sinkAway(ty = "M", fl = "/dev/null")
            zzz <- sinkAway(ty = "O", fl = "/dev/null")

            on.exit(expr = {sinkAway(ty = "M", fl = NULL);
                sinkAway(ty = "O", fl = NULL)})

        } else {
            # create files
            zzz <- tryCatch(suppressWarnings(
                file.create("tmp.log.mssg.txt")), error = function(e) {NULL})
            zzz <- tryCatch(suppressWarnings(
                file.create("tmp.log.outp.txt")), error = function(e) {NULL})
            zzz <- sinkAway(ty = "M", fl = "tmp.log.mssg.txt")
            zzz <- sinkAway(ty = "O", fl = "tmp.log.outp.txt")

            # on exit, do...
            on.exit(expr = {

                # delete temp files
                zzz <- sinkAway(ty = "M", fl = NULL)
                zzz <- sinkAway(ty = "O", fl = NULL)

                # delete temp files
                zzz <- tryCatch(suppressWarnings(
                    file.remove("tmp.log.mssg.txt")),
                    error = function(e) {NULL})

                zzz <- tryCatch(suppressWarnings(
                    file.remove("tmp.log.outp.txt")),
                    error = function(e) {NULL})
            })
        }
    }

    prepOP <- Prep4OptimizeParams(
        stack_img = stack_img, lnoise_range = lnoise_range,
        min.px.diam = min.px.diam, diameter_range = diameter_range,
        threshold_range = threshold_range, target_cell_num = target_cell_num,
        quantile.val = quantile.val, px.margin = px.margin,
        plot = plot, verbose = verbose)

    if (!prepOP$success)
        return(tc_obj)

    tmp_img <- prepOP$tmp_img
    all_params <- prepOP$all_params
    target_cell_num <- prepOP$target_cell_num

    all_results <- OptimizeParamsMainLoop(
        tmp_img = tmp_img, all_params = all_params,
        threads = threads, verbose = verbose)

    # Prep for return
    yy <- FianlizeOptiParams(
        all_results=all_results, all_params = all_params,
        target_cell_num=target_cell_num, plot=plot)

    tc_obj <- setProcessingStatus(tc_obj, "optimized_params", 1)
    tc_obj <- setOptimizedParams(tc_obj, yy$best_params, yy$all_params)
    return(tc_obj)
}




#' Validate Tracking Arguments
#'
#' Tool for Validate Tracking Arguments
#'
#' @param import_optiParam_from a trackedCells object, defaults to NULL
#' @param lnoise numeric, lnoise value
#' @param diameter numeric, lnoise value
#' @param threshold numeric, lnoise value
#' @param maxDisp numeric, lnoise value
#' @param memory_b numeric, should be 0
#' @param goodenough numeric, should be 0
#' @param show_plots logical, shall plots be shown
#' @param verbose logical, shall info be printed to console
#'
#' @return list of processed data
#'
#' @details This is an internal function supporting the CellTracker function.
#'
#' @examples
#' cellmigRation:::ValidateTrackingArgs()
#'
#'
#'
#' @keywords internal
ValidateTrackingArgs <- function(
    import_optiParam_from = NULL, lnoise = NULL, diameter = NULL,
    threshold = NULL, maxDisp = NULL, memory_b = 0, goodenough = 0,
    show_plots = FALSE, verbose = FALSE)
{

    track_params <- list()
    MS1 <- "Could not import optimal params from the provided object."
    MS2 <- "Is it a 'trackedCells' object? Did you run `OptimizeParams()`?"
    MS3 <- "Please, set all params for the analysis or run OptimizeParams()"

    tryCatch({suppressWarnings({

        # if a `import_optiParam_from` is specified, and
        if (!is.null(import_optiParam_from)) {
            chk.tco <- 0
            if ("trackedCells" %in% class(import_optiParam_from)) {
                TMP <- getProcessingStatus(import_optiParam_from)
                if(TMP$optimized_params == 1) {

                    # import
                    TMP2 <- getOptimizedParams(import_optiParam_from)
                    my.lnoise <- TMP2$auto_params$lnoise
                    my.diameter <- TMP2$auto_params$diameter
                    my.threshold <- TMP2$auto_params$threshold
                    chk.tco <- 1
                }
            }
            if (chk.tco == 0) {
                message(MS1)
                message(MS2)
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

        track_params <- list(
            lnoise = lnoise, diameter = diameter,
            threshold = threshold, maxDisp = 0,
            goodenough = goodenough,
            memory_b = memory_b, force_exec = force_exec,
            quiet = quiet, verbose = verbose, show_plots = show_plots
        )

        # Final check
        if (any(is.na(track_params))) {
            message(MS3)
            return(list())
        }

        if(!is.null(maxDisp))
            track_params$maxDisp <- maxDisp

    })}, error = function(e) { NULL })

    return(list(params=track_params, custom_params_flag=custom_params_flag))
}


#' Single Core Tracking Loop
#'
#' Tool for Single Core Tracking Loop
#'
#' @param FinalImage a list of numeric matrices (images)
#' @param min_frames_per_cell numeric, minimum number of frames
#' @param track_params a list of tracking parameters
#'
#' @return list of processed data
#'
#' @details This is an internal function supporting the CellTracker function.
#'
#' @examples
#' cellmigRation:::NonParallelTrackLoop(list(), 1, list())
#'
#' @keywords internal
NonParallelTrackLoop <- function(FinalImage, min_frames_per_cell, track_params)
{

    all_centroids <- list()
    all_b <- list()

    tryCatch({

        NumberImages <- length(FinalImage)
        lnoise <- track_params$lnoise
        diameter <- track_params$diameter
        threshold<- track_params$threshold
        verbose <- track_params$verbose
        show_plots <- track_params$show_plots

        # Init collectors
        all_centroids <- list()
        all_b <- list()

        for (i in seq(1,NumberImages,by=1)) {

            if (verbose)
                message(".", appendLF = FALSE)

            # generate an 1xP array with each column containing centroid output
            # for
            # individual frames
            a <- FinalImage[[i]]
            b <- tryCatch({bpass(
                image_array = a, lnoise = lnoise,
                lobject = diameter, threshold = threshold)},
                error = function(e) {NULL})
            pk <- tryCatch({pkfnd(
                im=b, th=threshold, sz=NextOdd(diameter))},
                error = function(e) {NULL})
            cnt <- tryCatch({cntrd(
                im = b, mx = pk, sz = NextOdd(diameter))},
                error = function(e) {NULL})

            if (show_plots) {
                VisualizeImg(img_mtx=b, las=1, main=paste0("Stack num. ", i))
                VisualizeCntr(
                    centroids = cnt, width_px = ncol(b), height_px = nrow(b))
            }

            # determine that frame s has at least 1 valid centroid
            if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
                all_centroids[[length(all_centroids) + 1]] <- cnt
                all_b[[length(all_b) + 1]] <- b

            } else {
                message(paste(
                    'No centroids detectd in frame', i, 'in the current stack'))
                message(
                    'Please, check nuclei validation settings ',
                    'for this image stack.')
            }
        }
        if (verbose)
            message("", appendLF = TRUE)
    }, error = function(e) {NULL})

    OUT <- list(centroids = all_centroids, all_b = all_b)
    return(OUT)
}


#' Multi Core Tracking Loop
#'
#' Tool for Multi Core Tracking Loop
#'
#' @param FinalImage a list of numeric matrices (images)
#' @param use.cores numeric, number of cores to use
#' @param min_frames_per_cell numeric, minimum number of frames
#' @param track_params a list of tracking parameters
#'
#' @return list of processed data
#'
#' @details This is an internal function supporting the CellTracker function.
#'
#' @importFrom parallel detectCores makeCluster clusterExport
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' cellmigRation:::ParallelTrackLoop(list(), 1, 1, list())
#'
#'
#'
#' @keywords internal
ParallelTrackLoop <- function(
    FinalImage, use.cores, min_frames_per_cell, track_params)
{

    all_centroids <- list()
    all_b <- list()
    j <- NULL

    tryCatch({

        NumberImages <- length(FinalImage)
        lnoise <- track_params$lnoise
        diameter <- track_params$diameter
        threshold<- track_params$threshold
        verbose <- track_params$verbose
        show_plots <- track_params$show_plots

        if (verbose) {
            cl <- suppressMessages(parallel::makeCluster(use.cores, outfile=""))
        } else {
            cl <- suppressMessages(parallel::makeCluster(use.cores))
        }
        suppressMessages(doParallel::registerDoParallel(cl))

        # Nothing to export! ""FinalImage", "all_params" automatically exported
        # stuffToExp <- c("FinalImage", "all_params")
        stuffToExp <- c()
        suppressMessages(parallel::clusterExport(cl, stuffToExp))

        ## %dopar%
        all_results <-
            tryCatch(foreach::foreach(
                j = (seq(1,NumberImages,by=1)),
                .verbose = verbose, .packages = "cellmigRation"
            ) %dopar% {

                # Verbose
                message(".", appendLF = FALSE)

                # generate an 1xP array with each column containing centroid
                # output
                # for individual frames
                a <- FinalImage[[j]]
                b <- tryCatch({
                    bpass(
                        image_array = a, lnoise = lnoise,
                        lobject = diameter, threshold = threshold)
                }, error = function(e) {NULL} )
                pk <- tryCatch({pkfnd(
                    im = b, th = threshold, sz = NextOdd(diameter))},
                    error = function(e) {NULL} )
                cnt <- tryCatch({cntrd(
                    im = b, mx = pk, sz = NextOdd(diameter))},
                    error = function(e) {NULL})

                # determine that frame s has at least 1 valid centroid
                if(! is.null(cnt) && is.data.frame(cnt) && nrow(cnt) > 0) {
                    tmpOUT <- list(cnt = cnt, b = b, j = j)
                } else {
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
            for (ii in seq(1,length(all_results),by=1)){
                bii <- all_results[[ii]]$b
                cntii <- all_results[[ii]]$cnt
                VisualizeImg(
                    img_mtx = bii, las = 1, main = paste0("Stack num. ", ii) )
                VisualizeCntr(
                    centroids = cntii, width_px = ncol(bii),
                    height_px = nrow(bii) )
                bii<-NULL; cntii<-NULL;
            }
        }
    }, error = function(e) {NULL})
    OUT <- list(centroids = all_centroids, all_b = all_b)
    return(OUT)
}


#' Post Processing Tracking Routine
#'
#' Tool for Post Processing Tracking Routine
#'
#' @param tc_obj a trackedCells object
#' @param all_centroids list, frame by frame centroids
#' @param track_params a list of tracking parameters
#'
#' @return list of processed data
#'
#' @details This is an internal function supporting the CellTracker function.
#'
#' @examples
#' cellmigRation:::PostProcessTracking(list(1), list(2), list(A=3))
#'
#'
#'
#' @keywords internal
PostProcessTracking <- function(tc_obj, all_centroids, track_params) {

    tracks2 <- list()
    all_centroids2 <- list()
    MS1 <- "The following maxDisp value was used for this analysis: "

    errTrk <- FALSE
    errTrk <- tryCatch({

        # Position list (reformated centroid data for track.m input)
        OUT_centroids <- all_centroids

        # Remove columns that contain brightness and square of radius
        # of gyration; this is the equivalent of position.m function
        med.celln.perfr <- tryCatch({median(do.call(c, lapply(
            all_centroids, nrow)), na.rm = TRUE)},
            error = function(e) {100})
        max.celln <- med.celln.perfr * 1.55
        min.celln <- med.celln.perfr * 0.45
        keepXX <- do.call(c, lapply(all_centroids, function(xx){
            nrow(xx) >= min.celln && nrow(xx) <= max.celln}))
        all_centroids <- all_centroids[keepXX]

        # Updated to avoid NO cells frames
        all_centroids2 <- list()
        for (ti in seq(1,length(all_centroids),by=1)) {
            ttmp <- all_centroids[[ti]]
            if (nrow(ttmp) > 0) {
                ttmp <- ttmp[, c(1,2)]
                ttmp$tau <- (length(all_centroids2) + 1)
                all_centroids2[[length(all_centroids2) + 1]] <- ttmp
            }
        }

        # create a matrix that contains centroid data in sequential
        # order by frame(tau)
        pos <- do.call(rbind, all_centroids2)
        if (!is.null(track_params$maxDisp)) {
            if (is.finite(track_params$maxDisp)) {
                if (track_params$maxDisp == 0) {
                    track_params$maxDisp <-  NULL
                }
            }
        }

        tracks2 <- NULL
        if (!is.null(track_params$maxDisp)) {
            tracks2 <- tryCatch(
                track(
                    xyzs = pos, maxdisp = track_params$maxDisp,
                    params = track_params),
                error = function(e) NULL
            )
        }

        if (is.null(tracks2)) {
            tracks2 <- NULL
            tmp.Area <- getCellImages(tc_obj)$dim$width_m *
                getCellImages(tc_obj)$dim$height_n
            max.disp <- as.integer(as.numeric(sqrt(tmp.Area)) / 5)
            allDisp <- seq(from=max.disp, to = 5, by = -3)
            jj0 <- 1
            while(is.null(tracks2) && jj0 < length(allDisp) ) {
                message(allDisp[jj0])
                tracks2 <- tryCatch(suppressMessages(
                    track(
                        xyzs = pos, maxdisp = allDisp[jj0],
                        params = track_params)), error = function(e) { NULL })
                if (!is.null(tracks2)) {
                    if(nrow(tracks2) < 1) { tracks2 <- NULL }
                }
                jj0 <- jj0 + 1
            }
            if (!is.null(tracks2)) {
                track_params$maxDisp <- allDisp[jj0]
                message(paste0(MS1, allDisp[jj0]))

            } else {
                message("a reasonable MaxDisp value couldn't be found! Sorry!")
                return(NULL)
            }
        }
    }, error = function(e) {TRUE})

    if (is.null(errTrk))
        errTrk <- FALSE

    if(errTrk)
        return(NULL)

    return(list(tracks = tracks2, centroids = all_centroids2,
                track_params=track_params,
                OUT_centroids = OUT_centroids,
                pos = pos))
}



#' Cell Tracker Main Loop
#'
#' Tool for Cell Tracker Main Loop
#'
#' @param FinalImage list of numeric matrices
#' @param threads numeric, number of cores to use
#' @param tc_obj trackingCell object
#' @param min_frames_per_cell numeric, minimum number of frames per cell
#' @param track_params a list of tracking parameters
#'
#' @return list of processed data
#'
#' @details This is an internal function supporting the CellTracker function.
#'
#' @examples
#' cellmigRation:::CellTrackerMainLoop(list(1), 1, 1, list(1))
#'
#'
#' @keywords internal
CellTrackerMainLoop <- function(
    FinalImage, threads, tc_obj, min_frames_per_cell, track_params) {

    # how many cores can we use?
    num_parallelCores <- threads
    debugging <- TRUE
    OUT <- list()

    tryCatch({
        max.cores <- parallel::detectCores()
        max.cores <- max.cores - 1
        max.cores <- ifelse(max.cores < 1, 1, max.cores)
        my.test <- 1 <= num_parallelCores & num_parallelCores <= max.cores
        use.cores <- ifelse(my.test, num_parallelCores, max.cores)

        # fix for NA cores
        if (is.na(use.cores)) {use.cores <- 1}

        # cores = 1, do not parallelize
        if (use.cores == 1) {
            y0 <- NonParallelTrackLoop(
                FinalImage = FinalImage, min_frames_per_cell, track_params)
        } else {
            y0 <- ParallelTrackLoop(
                FinalImage=FinalImage, track_params=track_params,
                min_frames_per_cell=min_frames_per_cell, use.cores=use.cores)
        }

        y1 <- suppressMessages(
            PostProcessTracking(
                tc_obj=tc_obj, all_centroids = y0$centroids,
                track_params = track_params))
        track_params <- y1$track_params
        tracks <- y1$tracks

        # Finalize
        # init num of cells
        init.cell.n <- length(unique(tracks[,4]))

        if (!is.null(min_frames_per_cell) &&
            is.numeric(min_frames_per_cell) &&
            length(min_frames_per_cell) == 1 &&
            min_frames_per_cell > 1) {

            all.cids <- table(tracks[, 4])
            all.cids <- data.frame(
                cell.id = as.numeric(names(all.cids)),
                count = as.numeric(all.cids)
            )

            all.cids <- all.cids[all.cids$count >= min_frames_per_cell, ]
            keep.cid <- all.cids$cell.id
            tracks <- tracks[ tracks[,4] %in% unique(keep.cid), ]

        } else {
            min_frames_per_cell <- 1
        }
        track_params$min_frames_per_cell <- min_frames_per_cell
        OUT <- list(
            tracks = tracks, track_params = track_params,
            init.cell.n = init.cell.n, all_b = y0$all_b,
            OUT_centroids = y1$OUT_centroids, pos = y1$pos)
    }, error = function(e) {NULL})

    return(OUT)
}



#' Compute Cell Tracks
#'
#' Analyze Stacks, detect cells in each frame, and analyze cell tracks
#' over time
#'
#' @details The lnoise param is used to guide a lowpass blurring
#' operation, while the lobject param is used
#' to guide a highpass background subtraction. The threshold param
#' is used for a background correction following
#' the initial image convolution
#' \itemize{
#'
#' \item \strong{lnoise}: Characteristic lengthscale of noise in
#' pixels.
#' Additive noise averaged over this length should vanish.
#' May assume any positive floating value.
#' May be also set to 0, in which case only the highpass
#' "background subtraction" operation is performed.
#'
#' \item \strong{lobject} Integer length in pixels somewhat larger
#' than a typical object.
#' Can also be set to 0, in which case only the lowpass "blurring"
#' operation defined by lnoise is done
#' without the background subtraction defined by lobject
#'
#' \item \strong{threshold} Numeric. By default,
#' after the convolution, any negative pixels are reset
#' to 0.    Threshold changes the threshhold for setting pixels to 0.
#' Positive values may be useful
#' for removing stray noise or small particles.
#'
#' }
#'
#'
#' @param tc_obj a \code{trackedCells} object.
#' @param import_optiParam_from a \code{trackedCells} object
#' (optional)
#' used to
#' import optimized parameters; can be NULL.
#' @param min_frames_per_cell numeric, minimum number of consecutive
#' frames in which
#' a cell shall be found in order to retain that cell in the final
#' cell tracks data.frame. Defaults to 1.
#' @param lnoise numeric, lnoise parameter; can be NULL if
#' OptimizeParams() has already been run
#' @param diameter numeric, diameter parameter; can be NULL if
#' OptimizeParams() has already been run
#' @param threshold numeric, threshold parameter; can be NULL
#' if OptimizeParams() has already been run
#' @param maxDisp numeric,    maximum displacement of a cell per
#' time interval.
#' When many cells are detected in each frame, small maxDisp values
#' should be used.
#' @param memory_b numeric, memory_b parameter as used in
#' the original track.m function.
#' In the current R implementation, only the value memory_b=0 is
#' accepted
#' @param goodenough numeric, goodenough parameter as used in
#' the original track.m function.
#' In the current R implementation, only the value goodenough=0
#' is accepted
#' @param threads integer, number of cores to use for parallelization
#' @param show_plots logical, shall cells detected in each frame of
#' the image stack be visualized
#' @param verbose logical, shall info about the progress of the cell
#' tracking job be printed
#' @param dryrun logical, shall a dryrun be performed
#'
#' @return a trackedCells object
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x <- get(data(TrackCellsDataset))
#' x <- CellTracker(x, dryrun=TRUE)
#' getTracks(x)[seq(1,12,by=1),]
#'
#' @importFrom parallel detectCores makeCluster clusterExport
#' stopCluster
#'
#' @export
CellTracker <- function(
    tc_obj, import_optiParam_from = NULL, min_frames_per_cell = 1,
    lnoise = NULL, diameter = NULL, threshold = NULL, maxDisp = NULL,
    memory_b = 0, goodenough = 0, threads = 1, show_plots = FALSE,
    verbose = FALSE, dryrun=FALSE)
{
    if (dryrun && getProcessingStatus(tc_obj)$track == 1)
        return(tc_obj)

    # get stuff
    stack_img <- getCellImages(tc_obj)
    optimal_params <- getOptimizedParams(tc_obj)

    if (length(optimal_params) > 0) {
        my.lnoise <- optimal_params$auto_params$lnoise
        my.diameter <- optimal_params$auto_params$diameter
        my.threshold <- optimal_params$auto_params$threshold
    } else {
        my.lnoise <- NULL
        my.diameter <- NULL
        my.threshold <- NULL
    }
    CPm <- 0
    if (!is.null(lnoise)) {my.lnoise <- lnoise; CPm <- 1}
    if (!is.null(diameter)) {my.diameter <- diameter; CPm <- 1}
    if (!is.null(threshold)) {my.threshold <- threshold; CPm <- 1}

    track_params <- ValidateTrackingArgs(
        import_optiParam_from = import_optiParam_from,
        lnoise = my.lnoise, diameter = my.diameter, threshold = my.threshold,
        maxDisp = maxDisp, memory_b = memory_b, goodenough = goodenough,
        show_plots = show_plots,  verbose = verbose)

    custom_params_flag <- CPm
    track_params <- track_params$params

    # Load stack
    stack <- stack_img

    InfoImage <- stack$attributes[[1]]
    mImage <- stack$dim$width_m
    nImage <- stack$dim$height_n
    NumberImages <- stack$dim$NumberImages
    FinalImage <- stack$images

    # locate centroids, via CentroidArray
    if (verbose)
        message("Computing centroid positions", appendLF = FALSE)

    ## Parallelize please
    if (!verbose) {

        if(file.exists("/dev/null")){
            zzz <- sinkAway(ty = "M", fl = "/dev/null")
            zzz <- sinkAway(ty = "O", fl = "/dev/null")
            on.exit(expr = {
                sinkAway(ty = "O", fl = NULL);
                sinkAway(ty = "M", fl = NULL)})
        } else {
            # create files
            zzz <- tryCatch(suppressWarnings(
                file.create("tmp.log.mssg.txt")),
                error = function(e) {NULL})
            zzz <- tryCatch(suppressWarnings(
                file.create("tmp.log.outp.txt")),
                error = function(e) {NULL})

            zzz <- sinkAway(ty = "M", fl = "tmp.log.mssg.txt")
            zzz <- sinkAway(ty = "O", fl = "tmp.log.outp.txt")

            # on exit, do...
            on.exit(expr = {
                sinkAway(ty = "M", fl = NULL)
                sinkAway(ty = "O", fl = NULL)
                zzz <- tryCatch(suppressWarnings(
                    file.remove("tmp.log.mssg.txt")),
                    error = function(e) {NULL})
                zzz <- tryCatch(suppressWarnings(
                    file.remove("tmp.log.outp.txt")),
                    error = function(e) {NULL})
            })
        }
    }

    ctrks <- CellTrackerMainLoop(
        FinalImage = FinalImage, threads = threads,
        tc_obj = tc_obj, min_frames_per_cell = min_frames_per_cell,
        track_params = track_params)

    # end num of cells
    end.cell.n <- length(unique(ctrks$tracks[,4]))
    sinkAway(ty = "M", fl = NULL)
    MS3 <- paste0(
        "Total number of unique cells across images: ", ctrks$init.cell.n,
        "; Cells retained after filtering: ", end.cell.n )
    message(MS3)

    # pack and return
    tc_obj <- setProcessedImages(tc_obj, list(images = ctrks$all_b))
    tc_obj <- setTrackedCentroids(tc_obj, ctrks$OUT_centroids)
    tc_obj <- setTrackedPositions(tc_obj, ctrks$pos)
    tc_obj <- setCellTracks(tc_obj, ctrks$tracks)
    tc_obj <- setAnalyticParams(tc_obj, ctrks$track_params)
    tc_obj <- setProcessingStatus(tc_obj, slot = "track", 1)
    tc_obj <- setProcessingStatus(
        tc_obj, slot = "custom_params", custom_params_flag)
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
#' x <- get(data(TrackCellsDataset))
#' getTracks(x)[seq(1,10,by=1),]
#'
#'
#' @export
getTracks <- function(tc_obj, attach_meta = FALSE)
{
    tmp <- getCellTracks(tc_obj)
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- c("Y", "X", "frame.ID", "cell.ID")
    TMP <- tmp[, c("frame.ID", "X", "Y", "cell.ID")]
    # Adjust as per S request
    TMP <- TMP[, c(4, 2, 3, 1)]
    rownames(TMP) <- NULL
    if(attach_meta && nrow(TMP) > 0) {

        TMP$tiff_file = getCellTrackMeta(tc_obj)$tiff_file
        TMP$experiment = getCellTrackMeta(tc_obj)$experiment
        TMP$condition = getCellTrackMeta(tc_obj)$condition
        TMP$replicate = getCellTrackMeta(tc_obj)$replicate
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
#' x <- get(data(TrackCellsDataset))
#' getPopulationStats(x)
#'
#'
#' @export
getPopulationStats <- function(tc_obj)
{
    if (getProcessingStatus(tc_obj)$stats == 1) {
        return(getCellTrackStats(tc_obj)$population)
    } else {
        message(
            "Stats have not been computed yet. ",
            "Please, run `ComputeTracksStats()`. Thanks."
        )
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
#' x <- get(data(TrackCellsDataset))
#' getCellsStats(x)
#'
#'
#'
#' @export
getCellsStats <- function(tc_obj)
{
    if (getProcessingStatus(tc_obj)$stats == 1) {
        return(getCellTrackStats(tc_obj)$cells)
    } else {
        message(
            "Stats have not been computed yet. ",
            "Please, run `ComputeTracksStats()`. Thanks."
        )
    }
}


#' Get MetaData
#'
#' Extract MetaData from a trackedCells object
#'
#' @param tc_obj a trackedCells object
#'
#' @return a list including four items:
#' tiff filename, experiment name, condition label, and replicate ID.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- get(data(TrackCellsDataset))
#' getCellsMeta(x0)
#'
#' @export
getCellsMeta <- function(tc_obj)
{
    tmp <- getCellTrackMeta(tc_obj)
    return(tmp)
}


#' Get Image Stacks
#'
#' Extract Images Stacks from a trackedCells object
#'
#' @param tc_obj a \code{trackedCells} object
#'
#' @return a list including stack images (formatted as numeric
#' matrices)
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- get(data(TrackCellsDataset))
#' y0 <- getImageStacks(x0)
#' graphics::image(y0[[1]])
#'
#' @export
getImageStacks <- function(tc_obj)
{
    tmp <- getCellImages(tc_obj)$images
    return(tmp)
}



#' Get Auto Optimized Parameters
#'
#' Extract Parameters that were automatically optimized
#'
#' @param tc_obj a \code{trackedCells} object
#'
#' @return a list including optimized parameter values
#' (lnoise, diameter, and threshold)
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x <- get(data(TrackCellsDataset))
#' getOptimizedParameters(x)
#'
#'
#' @export
getOptimizedParameters <- function(tc_obj)
{
    chk1 <- getProcessingStatus(tc_obj)$optimized_params
    y <- list()
    if (!is.null(chk1) && chk1 == 1) {
        y <- getOptimizedParams(tc_obj)$auto_params
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
#' @param experiment string, a label to describe the experiment
#' (optional). Can be NULL
#' @param condition string, a label to describe the experimental
#' condition (optional). Can be NULL
#' @param replicate string, a label to identify the replicate
#' (optional). Can be NULL
#'
#' @return a list including three items: experiment name, condition label,
#' and replicate ID.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#' \url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
#'
#' @examples
#' x0 <- get(data(TrackCellsDataset))
#' x0 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DMSO")
#' getCellsMeta(x0)
#'
#'
#' @export
setCellsMeta <- function(
    tc_obj, experiment = NULL,condition = NULL, replicate = NULL
) {

    if (is.null(experiment)) {
        experiment <- NA
    } else {
        experiment <- tryCatch(
            as.character(experiment[1]), error = function(e) NA
        )
    }

    if (is.null(replicate)) {
        replicate <- NA
    } else {
        replicate <- tryCatch(
            as.character(replicate[1]), error = function(e) NA
        )
    }

    if (is.null(condition)) {
        condition <- NA
    } else {
        condition <- tryCatch(
            as.character(condition[1]), error = function(e) NA
        )
    }

    FILENM <- getCellTrackMeta(tc_obj)$tiff_file

    tmp <- list(
        tiff_file = FILENM, experiment = experiment,
        condition = condition, replicate = replicate)

    tc_obj <- setTrackedCellsMeta(x = tc_obj, meta = tmp)
    return(tc_obj)
}


#' Aggregate trackedCells Objects
#'
#' Aggregate two or more trackedCells-class objects together.
#' Input objects must carry information of cell
#' tracks (otherwise an error will be raised). All tracks form
#' the different experiments/images are returned in a
#' large data.frame. A new unique ID is assigned to specifically
#' identify each cell track from each image/experiment.
#'
#' @param x a \code{trackedCells}-class object where cells have
#' already
#' been tracked
#' @param ... one or more trackedCells-class object(s) where cells
#' have already been tracked
#' @param meta_id_field string, can take one of the following values,
#' c("tiff_file", "experiment", "condition", "replicate"). Indicates
#' the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return An aggregate data.frame including all cells that were
#' tracked over two or more images/experiments.
#' The data.frame includes the following columns: "new.ID",
#' "frame.ID",
#' "X", "Y", "cell.ID", "tiff_name",
#' "experiment", "condition", "replicate". The "new.ID" uniquely
#' identifies a cell in a given image/experiment.
#'
#' @details each trackedCells-class object passed to this function
#' requires a unique identifier (such as a unique
#' tiff_file name). Any of the metadata columns can be used as
#' unique ID for an image/experiment. The function
#' will raise an error if non-unique identifiers are found across
#' the input objects.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' # Please, see the package vignette
#' # for an example of how to use this function.
#' # A pseudo-code example is shown below
#' # Let x0, x1, x2, ... be trackedCells-class objects
#' # with a non-empty tracks slot.
#' x0 <- get(data(TrackCellsDataset))
#' x0 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "CTRL")
#' x1 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DMSO")
#' x2 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DRUG")
#' y <- aggregateTrackedCells(x0, x1, x2, meta_id_field = "condition")
#' utils::head(y, 50)
#'
#'
#' @importFrom graphics par image
#'
#' @export
aggregateTrackedCells <-
    function(
        x, ...,
        meta_id_field = c("tiff_file", "experiment", "condition", "replicate")
) {
    # Inner fxs
    check_trobj <- function(xx) {
        RT <- FALSE
        if (!is.null(xx)) {
            if ("trackedCells" %in% class(xx)) {
                if (nrow(getCellTracks(xx)) > 0) {
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

    meta_id_field <- match.arg(
        arg = meta_id_field,
        choices = c("tiff_file", "experiment", "condition", "replicate"),
        several.ok = FALSE
    )
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
    all_ids <- do.call(c,lapply(big.list, function(xx) {
        getCellTrackMeta(xx)[[meta_id_field]] })
    )
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
    my_tracks[,"new.ID"] <-
        (my.mult * my_tracks[,"new.ID"]) + my_tracks[, "cell.ID"]

    # Adjust as per S request
    keep.colz <- c(
        'new.ID', 'X', 'Y', 'frame.ID', 'cell.ID', 'tiff_file',
        'experiment', 'condition', 'replicate')
    out <- my_tracks[, keep.colz]
    rownames(out) <- NULL

    return(out)
}


#' Filter an Aggregated Table of Cell Tracks
#'
#' Filter an Aggregated Table (data.frame) of cell tracks
#' (from multiple images/experiments) and
#' retain cell tracks from images/experiments of interest
#'
#' @param x data.frame, is an aggregated Table of Cell Tracks.
#' Must include the following columns:
#' "new.ID", "frame.ID", "X", "Y", "cell.ID", "tiff_name",
#' "experiment", "condition", "replicate"
#' @param id_list character vector, indicates the IDs (such as
#' tiff_filenames) to be retained in the
#' output data.frame
#' @param meta_id_field string, can take one of the following values,
#' c("tiff_file", "experiment", "condition", "replicate"). Indicates
#' the meta-data column used as unique ID for the image/experiment.
#' Can be abbreviated. Defaults to "tiff_file".
#'
#' @return data.frame, a filtered aggregated Table of Cell Tracks
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' A <- data.frame(new.ID = seq(1,10,by=1), frame.ID = seq(10,1,by=(-1)),
#'                 X = sample(seq(1,100,by=1), size = 10),
#'                 Y = sample(seq(1,100,by=1), size = 10),
#'                 cell.ID = c(rep(1, 5), rep(2, 5)),
#'                 tiff_file= c(rep("ii", 3), rep("jj", 5), rep('kk', 2)))
#' FilterTrackedCells(A, id_list = c("jj", "kk"), "tiff_file")
#'
#' @export
FilterTrackedCells <- function(x, id_list, meta_id_field) {

    meta_id_field <- match.arg(
        arg = meta_id_field,
        choices = c("tiff_file", "experiment", "condition", "replicate"),
        several.ok = FALSE
    )

    REQd <- c("new.ID", "frame.ID", "X", "Y", "cell.ID", meta_id_field)
    CHK1 <- sum(REQd %in% colnames(x)) == length(REQd)

    stopifnot(CHK1)

    xx <- x[x[, meta_id_field] %in% id_list, ]
    return(xx)
}



#' @title Data preprocessing for random migration (RM)
#'
#' @description This function allows preprocessing of the trajectory
#' data from random
#' migration (RM) experiments.
#'
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed
#' between successive frames in the time-lapse stack. Default is
#' 10 min.
#' @param PixelSize A numeric value of the physical size of a pixel.
#' Default is 1.24.
#' @param FrameN A numeric value of the number of frames.
#' Default is NULL
#' @param ExpName string, name of the experiment. Can be NULL
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @return An CellMig class object with preprocessed data.
#'
#' @examples
#' TrajectoryDataset <- get(data(TrajectoryDataset))
#' rmDF=TrajectoryDataset[seq(1,40,by=1),]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD, FrameN=30)
#'
#' @export
rmPreProcessing = function(
    object, PixelSize=1.24, TimeInterval=10, FrameN=NULL, ExpName = NULL) {
    msg <- NULL
    if ( ! is.data.frame(getCellMigSlot(object, "trajdata"))){
        msg <- c(msg, "input data must be data.frame")
    }
    if (!is.numeric(PixelSize)) {
        stop("PixelSize has to be a positive number")
    } else if ( PixelSize<= 0 ) {
        stop( "PixelSize has to be a positive number" )
    }
    if (!is.numeric(TimeInterval)) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    df <- getCellMigSlot(object, "trajdata")
    object <- setCellMigSlot(object, "adjDS", df)
    df<-df[,c(1,2,3)]    # Removing the unnecessary columns
    spl<-split(df,df[,1])
    cat("This dataset contains:", length(spl), "cell(s) in total\n")
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    tbd<-c()
    LENspl<-length(spl)
    for (i in seq(1, LENspl, by=1)){
        if (length(spl[[i]][,1])<4){ # Removing cells with < 4 tracked frames
            tbd<-c(tbd,i)
        }
    }
    spl[tbd]<-NULL
    df<-do.call(rbind.data.frame, spl)
    L<-length(df[,1])
    df[,seq(4,26,by=1)]<-rep(0,L)
    df[,27]<-rep(NA,L)    # to be used for migration type
    colnames(df)<-c(
        "ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P",
        "Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY",
        "New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F",
        "Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type"
    )
    ID_split <- split(df, df$ID)    #Splitting the data frame based on the ID
    cat("This dataset contains:",length(ID_split),
        "cell(s) with more than three steps in their tracks\n")
    for(j in seq(1,length(ID_split),by=1)){    # Having the ID =group order
        ID_split[[j]][1]=j
    }
    # adjusting x and y (starting from 0 & being multiplied by H)
    for(j in seq_len(length(ID_split))){
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        res <- t(vapply(seq(from = 2, to = MM, by = 1L), function(i){
            ID_split[[j]][i,4]=PixelSize*(
                (ID_split[[j]][i,2])-( ID_split[[j]][1,2])
            )
            ID_split[[j]][i,5]=PixelSize*(
                ( ID_split[[j]][1,3])-(ID_split[[j]][i,3])
            )
            return(as.numeric(ID_split[[j]][i,c(4,5)]))
        }, FUN.VALUE = numeric(2)))
        ID_split[[j]][seq(2,MM,by=1),c(4,5)] <- as.data.frame(res)
        ID_split[[j]][,c(4,5)] <- lapply(ID_split[[j]][,c(4,5)], as.numeric)
    }
    for(j in seq(1,length(ID_split), by=1)){    # removing the old x and y
        ID_split[[j]]=ID_split[[j]][-2] # removing x column
        ID_split[[j]]=ID_split[[j]][-2] # removing y[2] because x[2] is gone.
    }
    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    ID_split<-fixID1(ID_split,TimeInterval=TimeInterval)
    # creating values for cumsum and cumulative directionality ratio
    ID_split<-fixID2(ID_split)
    # creating values for    rel.ang.P    (step to the previous)
    ID_split<-fixID3(ID_split)
    # creating values for    cosine.P    based on rel.ang.P
    ID_split<-fixID4(ID_split)
    # Computing persistence time     (based on rel.ang.P)
    ID_split<-fixID5(ID_split,TimeInterval=TimeInterval)
    # Computing Acceleration
    ID_split<-fixID6(ID_split,TimeInterval=TimeInterval)
    size<-c()
    for(j in seq(1,length(ID_split),by=1)){
        size[j]<-length(ID_split[[j]][,1])
    }
    S<-summary(size)
    names(S)<-NULL
    if (is.null(FrameN)){
        IncompleteTracks<-subset(size,S<S[6])
        for(j in seq(1,length(ID_split),by=1)){
            ID_split[[j]]<-ID_split[[j]][seq(1,S[1],by=1),]
            ID_split[[j]][S[1],seq(4,25,by=1)]=0
            ID_split[[j]][S[1],10]=TimeInterval
            ID_split[[j]][1,10]=0
        }
        cat("The minimum number of steps: ",S[1],"\n")
        cat("The maximum number of steps: ",S[6],"\n")
        cat("Number of cells with a total number of steps less than ",S[6],
            "steps",":",length(IncompleteTracks),"\n")
        cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
        PreprocessedData <- ID_split
        object <- setCellMigSlot(object, "preprocessedDS", PreprocessedData)
    }else{
        if ( ! is.numeric(FrameN)) {
            stop( "FrameN has to be a positive number" )
        } else if ( FrameN<= 0 ) {
            stop( "FrameN has to be a positive number" )
        }
        if (FrameN>S[6]) {
            stop( paste0("No cells have ",FrameN, " steps in their tracks"))
        }
        okTracks<-c()
        for(j in seq(1,length(ID_split),by=1)){
            if (length(ID_split[[j]][,1])>=FrameN){
                okTracks<-c(okTracks,j)
            }
        }
        ID_split=ID_split[okTracks]
        for(j in seq(1,length(ID_split),by=1)){
            ID_split[[j]]<-ID_split[[j]][seq(1,FrameN,by=1),]
            ID_split[[j]][seq(1,FrameN,by=1),1]<-j
            ID_split[[j]][FrameN,seq(4,25,by=1)]=0
            ID_split[[j]][FrameN,10]=TimeInterval
            ID_split[[j]][1,10]=0
        }
        cat("The desired number of steps: ",FrameN,"\n")
        cat("The maximum number of steps: ",S[6],"\n")
        cat("Only: ", length(ID_split), " cells were selected","\n")
        cat("All the tracks of the selected cells are adjusted to have only ",
            FrameN," steps","\n")
        PreprocessedData<-ID_split
        object <- setCellMigSlot(object, "preprocessedDS", PreprocessedData)
    }
    return(object)
}


#' @title Data preprocessing for wound scratch assay (WSA).
#'
#' @description This function allows filtering of cells and
#' preprocessing
#' of the trajectory data from wound scratch assay (WSA) experiments.
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param PixelSize A numeric value of the physical size of a pixel.
#' @param FrameN A numeric value of the number of frames. Default
#' is NULL
#' @param imageH A numeric value of the image height.
#' @param woundH A numeric value of the image height.
#' @param upperE A numeric value of the upper edge of the wound.
#' @param lowerE A numeric value of the lower edge of the wound.
#' @param mar A numeric value of the margin to be used to narrow
#' the clearing zone inside the zone.
#' @param clearW A logical vector that allows removing the cells
#' within the wound. Default is TRUE.
#' @param ExpName string, name of the experiment. Can be NULL
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @return An CellMig class object with filtered, annotated
#' and preprocessed data.
#'
#' @examples
#' WSADataset <- get(data(WSADataset))
#' wasDF=WSADataset[seq(1,30,by=1),]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=20)
#'
#' @export
wsaPreProcessing = function(
    object, PixelSize=1.24, TimeInterval=10, FrameN=NULL,
    imageH=1500, woundH=600, upperE=400, lowerE=1000,
    mar=75, clearW=TRUE, ExpName = NULL
) {
    msg <- NULL
    if ( ! is.data.frame(getCellMigSlot(object, "trajdata"))){
        msg <- c(msg, "input data must be data.frame")
    }
    if ( ! is.numeric(PixelSize)) {
        stop( "PixelSize has to be a positive number" )
    } else if ( PixelSize<= 0 ) {
        stop( "PixelSize has to be a positive number" )
    }
    if ( ! is.numeric(imageH)) {
        stop( "imageH has to be a positive number" )
    } else if ( imageH<= 0 ) {
        stop( "PixelSize has to be a positive number" )
    }
    if ( ! is.numeric(woundH)) {
        stop( "woundH has to be a positive number" )
    } else if ( woundH<= 0 ) {
        stop( "woundH has to be a positive number" )
    }
    if ( ! is.numeric(upperE) ) {
        stop( "upperE has to be a positive number" )
    } else if ( upperE<= 0 ) {
        stop( "upperE has to be a positive number" )
    }
    if ( ! is.numeric(lowerE) ) {
        stop( "lowerE has to be a positive number" )
    } else if ( lowerE<= 0 ) {
        stop( "PixelSize has to be a positive number" )
    }
    if ( ! is.numeric(mar)) {
        stop( "mar has to be a positive number" )
    } else if ( mar<= 0 ) {
        stop( "mar has to be a positive number" )
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    if (clearW == TRUE){
        dff<-getCellMigSlot(object, "trajdata")
        dff<-dff[,c(1,2,3,4)] # Removing the unnecessary columns
        splitFORu<-split(dff,dff[,1])
        for (i in seq(1,length(splitFORu),by=1)){
            if ((
                splitFORu[[i]][1,3]>= (upperE + mar) &
                splitFORu[[i]][1,3]<= (lowerE -mar)
                ) &
                (splitFORu[[i]][1,4]<=20 )
            ) {     # to remove cells within the wound
                splitFORu[[i]]<-NA
            }
        }
        deletedCells<-splitFORu[is.na(splitFORu)]    # To get the deleted cells
        cat(paste0(
            length(deletedCells),
            " cells were inside the wound and they have been removed"), "\n")
        Ftable<-splitFORu[!is.na(splitFORu)]
        Ftable<-do.call(rbind.data.frame, Ftable) # convert list to data frame
        rownames(Ftable)<-NULL
        object <- setCellMigSlot(object, "adjDS", Ftable)
    }else{
        TMP0 <- getCellMigSlot(object, "trajdata")
        object <- setCellMigSlot(object, "adjDS", TMP0)
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    ####### to set the cells orientation
    hh<- upperE + ((lowerE-upperE)/2)
    CellOr<-c()
    finaltable<- getCellMigSlot(object, "adjDS")
    finaltable<-split(finaltable,finaltable[,1])
    for (i in seq(1,length(finaltable),by=1)){
        if (finaltable[[i]][1,3]<=hh){
            CellOr[i]=0
        }else{
            CellOr[i]=1
        }
    }
    object <- setCellMigSlot(object, "cellpos", CellOr)
    df<-getCellMigSlot(object, "adjDS")
    df<-df[,c(1,2,3)] # Removing the unnecessary columns
    spl<-split(df,df[,1])
    cat("This dataset contains:", length(spl), "cell(s) in total\n")
    tbd<-c()
    LENspl<-length(spl)
    for (i in seq(1,LENspl,by=1)){
        if (length(spl[[i]][,1])<4){ # Removing cells with < 4 tracked frames
            tbd<-c(tbd,i)
        }
    }
    spl[tbd]<-NULL
    df<-do.call(rbind.data.frame, spl)
    L<-length(df[,1])
    df[,seq(4,26,by=1)]<-rep(0,L)
    df[,27]<-rep(NA,L) # to be used for migration type
    colnames(df)<-c(
        "ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P",
        "Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY",
        "New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F",
        "Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type"
    )
    ID_split <- split(df, df$ID) #Splitting the data frame based on the ID
    cat("This dataset contains: ",length(ID_split),
        "Cells with more than three steps in their tracks","\n")
    for(j in seq(1,length(ID_split),by=1)){ # Having the ID =group order
        ID_split[[j]][1]=j
    }
    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    ID_split<-fixID1(ID_split,TimeInterval=TimeInterval)
    # creating values for cumsum and cumulative directionality ratio
    ID_split<-fixID2(ID_split)
    # creating values for    rel.ang.P    (step to the previous)
    ID_split<-fixID3(ID_split)
    # creating values for    cosine.P    based on rel.ang.P
    ID_split<-fixID4(ID_split)
    # Computing persistence time     (based on rel.ang.P)
    ID_split<-fixID5(ID_split,TimeInterval=TimeInterval)
    # Computing Acceleration
    ID_split<-fixID6(ID_split,TimeInterval=TimeInterval)

    size<-c()
    for(j in seq(1,length(ID_split),by=1)){
        size[j]<-length(ID_split[[j]][,1])
    }
    S<-summary(size)
    names(S)<-NULL
    if (is.null(FrameN)){
        IncompleteTracks<-subset(size,S<S[6])
        for(j in seq(1,length(ID_split),by=1)){
            ID_split[[j]]<-ID_split[[j]][seq(1,S[1],by=1),]
            ID_split[[j]][S[1],seq(4,25,by=1)]=0
            ID_split[[j]][S[1],10]=TimeInterval
            ID_split[[j]][1,10]=0
        }
        cat("The minimum number of steps: ",S[1],"\n")
        cat("The maximum number of steps: ",S[6],"\n")
        cat("Number of cells with a total number of steps less than ",
            S[6],"steps",":",length(IncompleteTracks),"\n"
        )
        cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
        PreprocessedData<-ID_split
        object<-setCellMigSlot(object, "preprocessedDS", PreprocessedData)
    }else{
        if ( ! is.numeric(FrameN)) {
            stop( "FrameN has to be a positive number" )
        } else if ( FrameN<= 0 ) {
            stop( "FrameN has to be a positive number" )
        }
        if (FrameN>S[6]) {
            stop( paste0("No cells have ",FrameN, " steps in their tracks"))
        }
        okTracks<-c()
        for(j in seq(1,length(ID_split),by=1)){
            if (length(ID_split[[j]][,1])>=FrameN){
                okTracks<-c(okTracks,j)
            }
        }
        ID_split=ID_split[okTracks]
        for(j in seq(1,length(ID_split),by=1)){
            ID_split[[j]]<-ID_split[[j]][seq(1,FrameN,by=1),]
            ID_split[[j]][seq(1,FrameN,by=1),1]<-j
            ID_split[[j]][FrameN,seq(4,25,by=1)]=0
            ID_split[[j]][FrameN,10]=TimeInterval
            ID_split[[j]][1,10]=0
        }
        cat("The desired number of steps: ",FrameN,"\n")
        cat("The maximum number of steps: ",S[6],"\n")
        cat("Only: ", length(ID_split), " cells were selected","\n")
        cat("All the tracks of the selected cells are adjusted to have only ",
            FrameN," steps","\n")
        PreprocessedData<-ID_split
        object <- setCellMigSlot(object, "preprocessedDS", PreprocessedData)
    }
    return(object)
}


#' @title A 2D rose-plot
#'
#' @description Plotting the trajectory data of all cells.
#'
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param Type has to be one of the following: c("p", "l", "b", "o")
#' "p": Points;
#' "l": Lines;
#' "b": Both;
#' "o": Both "overplotted".
#' @param FixedField logical(1) Allows generating a plot with
#' fixed field 800um x 800um. Default is TRUE.
#' @param export if `TRUE` (default), exports plot to JPG file
#' @return A 2D rose-plot showing the tracks of all cells.
#' @param ExpName string, name of the experiment. Can be NULL
#' @details    The visualization shows centered trajectories where
#' the starting point of each track is located at the origin
#' of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' plotAllTracks(object=rmTD, Type="l", FixedField=TRUE,export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points lines
#' @export
plotAllTracks= function(
    object, Type="l", FixedField=TRUE, export=FALSE, ExpName = NULL) {
    if ( ! ( Type %in% c("p","l","b","o") ) ) {
        stop("Type has to be one of the following: p, l, b, o")
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
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
        graphics::plot(
            Object[[1]][seq(1,Step,by=1),2], Object[[1]][seq(1,Step,by=1),3],
            type=Type, xlab="X (um)", ylab="Y (um)",
            col=color[1], las=1, xlim=c(-400,400),cex.lab=0.7,
            ylim=c(-400,400), main=ExpName
        )
        for(n in seq(2,Len,by=1)){
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
        if (export) {
            grDevices::jpeg(
                paste0(ExpName,"_All_tracks_plot.jpg"),
                width = 4, height = 4, units = 'in', res = 300
            )
        }
        graphics::plot(
            Object[[1]][seq(1,Step,by=1),2], Object[[1]][seq(1,Step,by=1),3],
            type=Type, xlab="X (um)", ylab="Y (um)", col=color[1],
            las=1, xlim=c(-400,400), ylim=c(-400,400), main=ExpName
        )
        for(n in seq(2,Len,by=1)){
            graphics::points(
                Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n]
            )
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
        for(j in seq(1,Len,by=1)){
            minX=min(Object[[j]][seq(1,Step,by=1),2])
            minY=min(Object[[j]][seq(1,Step,by=1),3])
            maxX=max(Object[[j]][seq(1,Step,by=1),2])
            maxY=max(Object[[j]][seq(1,Step,by=1),3])
            MinX[j]<-c(minX)
            MaxX[j]<-c(maxX)
            MinY[j]<-c(minY)
            MaxY[j]<-c(maxY)
        }
        RangeX=c(MinX,MaxX)
        RangeY=c(MinY,MaxY)
        graphics::plot(
            Object[[1]][seq(1,Step,by=1),2], Object[[1]][seq(1,Step,by=1),3],
            type=Type, xlab="X (um)", ylab="Y (um)",
            col=color[1], las=1, xlim=range(RangeX),cex.lab=0.7,
            ylim=range(RangeY), main=ExpName)
        for(n in seq(2,Len,by=1)){
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
        if (export) {
            grDevices::jpeg(
                paste0(ExpName,"_All_tracks_plot.jpg"),
                width = 4, height = 4, units = 'in', res = 300
            )
        }
        graphics::plot(
            Object[[1]][seq(1,Step,by=1),2], Object[[1]][seq(1,Step,by=1),3],
            type=Type,xlab="X (um)", ylab="Y (um)", col=color[1],
            las=1, xlim=range(RangeX), ylim=range(RangeY), main=ExpName
        )
        for(n in seq(2,Len,by=1)){
            graphics::points(
                Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n]
            )
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
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param Type has to be one of the following: c("p", "l", "b", "o")
#' @param export if `TRUE` (default), exports plot to JPG file
#' "p": Points;
#' "l": Lines;
#' "b": Both;
#' "o": Both "overplotted".
#' @param FixedField logical(1) Allows generating a plot with
#' fixed field 800um x 800um. Default is TRUE.
#' @param celNum A numeric value showing the desired number of cells
#' to be plotted.
#' @return A 2D rose-plot showing the tracks of sample cells
#' selected randomly based on the desired number of cells selected
#' by the user.
#' @param ExpName string, name of the experiment. Can be NULL
#' @details    The visualization shows centered trajectories where
#' the starting point of each track is located at the
#' origin of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' preProcCellMig <- get(data(preProcCellMig))
#' plotSampleTracks(preProcCellMig, Type="l", FixedField=TRUE,
#'                  celNum=5, export=FALSE, ExpName = NULL)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points lines
#' @export
plotSampleTracks= function(
    object, Type="l", celNum=35,
    FixedField=TRUE,export=FALSE,
    ExpName = NULL
) {
    if ( ! ( Type %in% c("p","l","b","o") ) ) {
        stop("Type has to be one of the following: p, l, b, o")
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    Len<-length(Object)
    if ( ! is.numeric(celNum) ) {
        stop( "celNum has to be a positive number" )
    } else if ( celNum > Len ) {
        stop( "The cellNum should be less than the total number of cells" )
    }
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
    OBJ<-seq(1,Len,by=1)
    cells=sample(OBJ,celNum)
    cells=sort(cells)
    cat(paste0("The plot contains the following cells: "),"\n")
    cat(cells,"\n")
    if ( FixedField == TRUE){
        graphics::plot(
            Object[[cells[1]]][seq(1,Step,by=1),2],
            Object[[cells[1]]][seq(1,Step,by=1),3],
            type=Type, xlab="X (um)", ylab="Y (um)",
            col=color[1], las=1, xlim=c(-400,400),cex.lab=0.7,
            ylim=c(-400,400), main=ExpName
        )
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
        if (export) {
            grDevices::jpeg(
                paste0(ExpName,"_Sample_tracks_plot.jpg"),
                width = 4, height = 4, units = 'in', res = 300
            )
        }
        graphics::plot(
            Object[[cells[1]]][seq(1,Step,by=1),2],
            Object[[cells[1]]][seq(1,Step,by=1),3],
            type=Type,
            xlab="X (um)", ylab="Y (um)", col=color[1],
            las=1, xlim=c(-400,400), ylim=c(-400,400), main=ExpName
        )
        CELLS<-cells[-1]
        for(n in CELLS){
            graphics::points(
                Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n]
            )
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
            minX=min(Object[[j]][seq(1,Step,by=1),2])
            minY=min(Object[[j]][seq(1,Step,by=1),3])
            maxX=max(Object[[j]][seq(1,Step,by=1),2])
            maxY=max(Object[[j]][seq(1,Step,by=1),3])
            MinX<-c(MinX,minX)
            MaxX<-c(MaxX,maxX)
            MinY<-c(MinY,minY)
            MaxY<-c(MaxY,maxY)
        }
        RangeX=c(MinX,MaxX)
        RangeY=c(MinY,MaxY)
        graphics::plot(
            Object[[cells[1]]][seq(1,Step,by=1),2],
            Object[[cells[1]]][seq(1,Step,by=1),3],
            type=Type, xlab="X (um)", ylab="Y (um)",
            col=color[1], las=1, xlim=range(RangeX),cex.lab=0.7,
            ylim=range(RangeY), main=ExpName
        )
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

        if (export) {
            grDevices::jpeg(
                paste0(ExpName,"_Sample_tracks_plot.jpg"),
                width = 4, height = 4, units = 'in', res = 300
            )
        }
        graphics::plot(
            Object[[cells[1]]][seq(1,Step,by=1),2],
            Object[[cells[1]]][seq(1,Step,by=1),3],
            type=Type,
            xlab="X (um)", ylab="Y (um)", col=color[1],
            las=1, xlim=range(RangeX), ylim=range(RangeY), main=ExpName
        )
        CELLS<-cells[-1]
        for(n in CELLS){
            graphics::points(
                Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n]
            )
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
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#' @param interactive logical, shall the 3D plot be generated in
#' a interactive fashion
#'
#' @return A 3D rose-plot showing the tracks of all cells.
#' @details The 3D visualization shows centered trajectories where
#' the starting point of each track is located at the origin
#' of the coordinate system (X=0,Y=0).
#'
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' if (Sys.info()[["sysname"]] != "Darwin") {
#'   # interactive shall be set to TRUE (default)
#'   rmTD <- get(data(preProcCellMig))
#'   plot3DAllTracks(rmTD, VS=3, size=2, interactive = FALSE)
#' }
#'
#' @note This function requires the \code{rgl} package to be installed on your
#' system.
#'
#' @importFrom grDevices rainbow
#'
#' @export
plot3DAllTracks= function(object, VS=3, size=2, interactive = TRUE) {
    if (!requireNamespace("rgl", quietly=TRUE)) {
        stop(
            "Could not load the rgl package. Please install ",
            "the rgl package in order to use plot3DAllTracks()."
        )
    }
    if(!interactive) {
        message("Modality not supported")
        return(NULL)
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if ( ! is.numeric(VS) ) {
        stop( "VS has to be a positive number" )
    } else if ( VS<= 0 ) {
        stop( "VS has to be a positive number" )
    }
    if ( ! is.numeric(size) ) {
        stop( "size has to be a positive number" )
    } else if ( size<= 0 ) {
        stop( "size has to be a positive number" )
    }
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
    for(i in seq(1,Len,by=1)){
        firstNum=((i*Step)-Step)+1
        RowNum=seq(firstNum,(i*Step),by=1)
        plotTable[RowNum,1]= Object[[i]][seq(1,Step,by=1),2]
        plotTable[RowNum,2]=Object[[i]][seq(1,Step,by=1),3]
        plotTable[RowNum,3]=i*VS
        col= c(rep(color[i],Step))
        coll=c(coll,col)
    }
    rgl::plot3d(
        plotTable, col=coll, type="p", size=size, axes=FALSE,
        xlab=" ", ylab=" ",zlab=" "
    )
}


#' @title A 3D rose-plot
#' @description Plotting the trajectory data of particular cells in
#' 3D.
#'
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#' @param cells A numeric vector containing the cell's numbers to be
#' plotted.
#' @param interactive logical, shall a 3D plot built in an interactive way.
#'
#' @return A 3D rose-plot showing the tracks of particular cells.
#' @details    The 3D visualization shows centered trajectories
#' where the starting point of each track is located at the origin
#' of the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' if (Sys.info()[["sysname"]] != "Darwin") {
#'   # interactive shall be set to TRUE (default)
#'   rmTD <- get(data(preProcCellMig))
#'   plot3DTracks(rmTD, VS=3, size=2, cells=seq(1,5,by=1), interactive = FALSE)
#' }
#'
#' @note This function requires the \code{rgl} package to be installed on your
#' system.
#'
#' @importFrom grDevices rainbow
#' @importFrom stats complete.cases
#'
#' @export
plot3DTracks= function(object, VS=3, size=2, cells, interactive=TRUE) {
    if (!requireNamespace("rgl", quietly=TRUE)) {
        stop(
            "Could not load the rgl package. Please install ",
            "the rgl package in order to use plot3DTracks()."
        )
    }
    Object<- getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if (!interactive) {
        message("Modality not supported")
        return(NULL)
    }

    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if ( ! is.numeric(VS) ) {
        stop( "VS has to be a positive number" )
    } else if ( VS<= 0 ) {
        stop( "VS has to be a positive number" )
    }
    if ( ! is.numeric(size) ) {
        stop( "size has to be a positive number" )
    } else if ( size<= 0 ) {
        stop( "size has to be a positive number" )
    }
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
        RowNum=seq(firstNum,(i*Step),by=1)
        plotTable[RowNum,1]= Object[[i]][seq(1,Step,by=1),2]
        plotTable[RowNum,2]=Object[[i]][seq(1,Step,by=1),3]
        plotTable[RowNum,3]=i*VS
        col= c(rep(color[i],Step))
        coll=c(coll,col)
        NewplotTable<-plotTable[stats::complete.cases(plotTable),]

    }
    rgl::plot3d(
        NewplotTable, col=coll, type="p", size=size, axes=FALSE,
        xlab=" ", ylab=" ",zlab=" "
    )
}


#' @title A graphical display of the track of each cell.
#' @description Plotting the trajectory data of each cell.
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param Type has to be one of the following: [p, l, b, o]
#' "p": Points
#' "l": Lines
#' "b": Both
#' "o": Both "overplotted"
#' @param FixedField logical(1) Allows generating individual plots
#' with fixed field. Default is TRUE.
#' @param export if `TRUE` (default), exports plot to JPG file
#' @param ExpName string, name of the experiment. Can be NULL
#' @return 2D rose-plots of the cells' track Separately.
#' @details    The visualization shows centered trajectories where the
#' starting point of each track is located at the origin of
#' the coordinate system (X=0,Y=0).
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' PlotTracksSeparately(rmTD,Type="b", FixedField=FALSE, export = FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom graphics plot points title
#'
#'@export
PlotTracksSeparately= function(
    object,Type="l", FixedField=TRUE, export=FALSE,ExpName = NULL
){
    if ( ! ( Type %in% c("p","l","b","o") ) ) {
        stop("Type has to be one of the following: p, l, b, o")
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    Len<-length(Object)
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
    new.fld <-paste0(ExpName,"_Tracks")
    if (export) {
        cat(
            paste0(
                Len," plots will be generated in a folder called:",
                ExpName,"_Tracks","\n"
            )
        )

        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }

        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
    }
    if ( FixedField == TRUE){
        MinX<-c()
        MaxX<-c()
        MinY<-c()
        MaxY<-c()
        for(j in seq(1,Len,by=1)){
            minX=min(Object[[j]][seq(1,Step,by=1),2])
            minY=min(Object[[j]][seq(1,Step,by=1),3])
            maxX=max(Object[[j]][seq(1,Step,by=1),2])
            maxY=max(Object[[j]][seq(1,Step,by=1),3])
            MinX[j]<-c(minX)
            MaxX[j]<-c(maxX)
            MinY[j]<-c(minY)
            MaxY[j]<-c(maxY)
        }
        RangeX=c(MinX,MaxX)
        RangeY=c(MinY,MaxY)
        for(n in seq(1,Len,by=1)){
            if (export){
                plot_name <-    paste0(ExpName,"_Track_Plot_",n,".jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4, height = 4, units = 'in',
                    res = 300
                )
                graphics::plot(
                    Object[[n]][seq(1,Step,by=1),2],
                    Object[[n]][seq(1,Step,by=1),3],
                    type=Type, xlab="x (um)", ylab="y (um)",
                    col=color[n], las=1, xlim=range(RangeX), ylim=range(RangeY)
                )
                x=c(min(RangeX)-100,max(RangeX)+100)
                y=c(0,0)
                graphics::lines(x, y, type='l', col="black")
                x=c(0,0)
                y=c(min(RangeY)-100,max(RangeY)+100)
                graphics::lines(x, y, type='l', col="black")
                end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
                graphics::points(end,pch=16,col=color[n], cex = 1)
                graphics::title(
                    main=paste0("Cell Number    ", n),col.main="black"
                )
                grDevices::dev.off()
            }
        }
    }else{
        for(n in seq(1,Len,by=1)){
            RangeX= Object[[n]][seq(1,Step,by=1),2]
            RangeY= Object[[n]][seq(1,Step,by=1),3]
            if (export){
                plot_name <- paste0(ExpName,"_Track_Plot_",n,".jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4, height = 4, units = 'in',
                    res = 300
                )
                graphics::plot(
                    Object[[n]][seq(1,Step,by=1),2],
                    Object[[n]][seq(1,Step,by=1),3],
                    type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,
                    xlim=range(RangeX),ylim=range(RangeY)
                )
                x=c(min(RangeX)-100,max(RangeX)+100)
                y=c(0,0)
                graphics::lines(x, y, type='l', col="black")
                x=c(0,0)
                y=c(min(RangeY)-100,max(RangeY)+100)
                graphics::lines(x, y, type='l', col="black")
                end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
                graphics::points(end,pch=16,col=color[n], cex = 1)
                graphics::title(
                    main=paste0("Cell Number    ", n),col.main="black"
                )
                if (export) grDevices::dev.off()
            }
        }
    }
}



#' @title Persistence and Speed
#' @description The PerAndSpeed() generates data and plots for
#' persistence and speed.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param PtSplot A logical vector that allows generating individual
#' plots of persistence time vs speed per cell. Default is TRUE.
#' @param AllPtSplot    A logical vector that allows generating a plot
#' of persistence time vs speed for all cells. Default is TRUE.
#' @param ApSplot A logical vector that allows generating individual
#' plots of angular persistence vs speed per cell. Default is TRUE.
#' @param AllApSplot A logical vector that allows generating a plot
#' of angular persistence vs speed of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#' @param ExpName string, indicates the name of the experiment.
#' Can be NULL
#'
#' @return An CellMig class object with a data frame and plots.
#' The data frame is stored in the PerAanSpeedtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' rmTD <- PerAndSpeed(rmTD,TimeInterval=10, export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off chull
#' @importFrom stats cor.test lm median sd
#' @importFrom graphics plot title abline legend
#' @importFrom sp Polygon
#' @importFrom matrixStats rowMedians
#' @importFrom vioplot vioplot
#' @importFrom utils write.csv
#' @importFrom methods slot
#'
#'
#'@export

PerAndSpeed= function(
    object, TimeInterval=10, PtSplot=TRUE, AllPtSplot=TRUE,
    ApSplot=TRUE, AllApSplot=TRUE, export=FALSE, ExpName=NULL) {
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    if (export) {
        new.fld <-paste0(ExpName,"-PerResults")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }

        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
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
    # creating a table to store the persistence results
    PerResultsTable<-fixPER1(Object,TimeInterval=TimeInterval)

    # creating a table to store the mean velocity with correspondence
    # with the persistence time
    PerResultsTable<-fixPER2(
        Object,PerResultsTable=PerResultsTable, PtSplot=PtSplot,
        AllPtSplot=AllPtSplot, export=export, color=color,
        TimeInterval=TimeInterval,ExpName = ExpName, new.fld = new.fld)
    # calculating the Mean.Square.velocity for each cell
    PerResultsTable<-fixPER3(
        Object,PerResultsTable=PerResultsTable, ApSplot=ApSplot,
        AllApSplot=AllApSplot, export=export,color=color,
        TimeInterval=TimeInterval, ExpName = ExpName, new.fld = new.fld)

    rownames(PerResultsTable)<-c(
        "Cell Number","Mean Persist Time (min)",
        "Mean persist Deviating Time (min)",
        "Persistence Ratio", "Maximum Persistence period (min)",
        "Persistence Time vs Speed (SCC)","RMSS (um per h)",
        "Maximum Speed (um per h)","Minimum Speed (um per h)",
        "Mean Angular Persistence (cosine)",
        "Instantaneous Speed vs Angular Persistence (SCC)","Covered Area (um2)",
        "Segmented Covered Area (um2)","Empty Area (um2)",
        "Number of complete rotations",
        "Number of canceled rotations","Mean Persistence Angle (cosine)",
        "Mean Deviating Angle (cosine)","Median Speed","Speed QBCV",
        "Mean Speed (um per h)","Speed standard deviation (um per h)"
    )
    RM1<-round(
        matrixStats::rowMedians(as.matrix(PerResultsTable),na.rm = TRUE),
        digits=3)
    tmp.range <- c(2,3,4,5,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22)
    PerResultsTable[tmp.range, (length(Object)+1)] <- RM1[tmp.range]
    # Needs fixing. are you starting from 0??
    tmp.range <- seq(1,ncol(PerResultsTable),by=1)-1
    RMSS<-as.numeric(PerResultsTable[7,tmp.range])
    if (export) {
        plot_name <- paste0(ExpName,"_RMSS_profile_of_all_cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    if (export) {
        plot_name <- paste0(ExpName,"_RMSS_profile_of_all_cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    cells<-seq(1, (length(PerResultsTable[1,])-1), by = 1)
    MS<-max(RMSS)
    graphics::plot(
        cells,RMSS,pch=16,type="o",ylab = 'RMSS(um/h)',xlab = 'Cells',las=1,
        ylim = c(0, MS)
    )
    graphics::title(main="RMSS Profile of all cells",cex.main = 1)
    graphics::abline(h=stats::median(RMSS[which(!is.na(RMSS))]),col="red")
    graphics::abline(h=mean(RMSS[which(!is.na(RMSS))]),col="blue")
    graphics::legend(
        1, y=200, legend=c("Mean RMSSs","Median RMSSs"),
        col=c("blue","red"),
        lty=1, cex=0.8
    )
    if (export) grDevices::dev.off()
    if (export) {
        plot_name <-    paste0(ExpName,"_RMSS_ViolinPlot_of_all_cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    graphics::plot(
        1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '',
        ylab = 'RMSS(um/h)', xaxt = 'n',las=1
    )
    graphics::title("RMSS of all cells",cex.main = 1)
    vioplot::vioplot(RMSS, at = 1, add = TRUE, col = "gray")
    if (export) grDevices::dev.off()
    tmp.rangee <- seq(1,length(PerResultsTable[1,]),by=1)-1
    SPEED<-as.numeric(PerResultsTable[19,tmp.rangee])
    if (export) {
        plot_name <-    paste0(ExpName,"_Speed_profile_of_all_cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    cells<-seq(1,(ncol(PerResultsTable)-1),by=1)
    MS<-max(SPEED)
    graphics::plot(
        cells,SPEED,pch=16,type="o",ylab = 'Speed(um/h)',xlab = 'Cells',
        las=1,ylim = c(0, MS)
    )
    graphics::title(main="Speed Profile of all cells",cex.main = 0.7)
    graphics::abline(h=stats::median(SPEED[which(!is.na(SPEED))]),col="red")
    graphics::abline(h=mean(SPEED[which(!is.na(SPEED))]),col="blue")
    graphics::legend(
        1, y=200, legend=c("Mean Speed","Median Speed"),
        col=c("blue","red"),lty=1, cex=0.8
    )
    if (export) grDevices::dev.off()
    if (export) {
        plot_name <-    paste0(ExpName,"_Speed_ViolinPlot_of_all_cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    graphics::plot(
        1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '',
        ylab = 'Speed(um/h)', xaxt = 'n',las=1
    )
    graphics::title("Speed of all cells",cex.main = 1)
    vioplot::vioplot(SPEED, at = 1, add = TRUE, col = "gray")
    if (export) grDevices::dev.off()
    if (export) {
        plot_name <- paste0(
            ExpName,"_Instantaneous_Speed_VS_Persistence Ratio_all_cells.jpg"
        )
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    SPEED<-as.numeric(
        PerResultsTable[21,seq(1,length(PerResultsTable[1,]),by=1)-1])
    PerR<- as.numeric(
        PerResultsTable[4,seq(1,length(PerResultsTable[1,]),by=1)-1])
    MS<-max(SPEED)
    graphics::plot(
        SPEED,PerR,pch=16,type="p",xlab="Average Instantaneous Speed (um/h)",
        ylab=" Persistence Ratio",col="black",las=1,xlim=c(0,MS),ylim=c(0,1)
    )
    reg<-stats::lm(PerR~SPEED)
    graphics::abline(reg,untf=FALSE,col="red")
    # testing the correlation
    s<-cor.test( ~ SPEED+ PerR, method = "spearman",exact=FALSE)
    ss<-unlist(s[4])
    SCC<-round(ss, digits = 3)
    graphics::title(
        main=paste0("All Cells Average Speed vs Persistence Ratio "),
        cex.main =0.7,
        sub=paste0("Spearman's rank correlation coefficient = ",SCC),
        col.sub="red"
    )
    if (export) grDevices::dev.off()
    PerResultsTable[1,(length(Object)+1)]<-"All Cells"
    object <- setCellMigSlot(object, "PerAanSpeedtable", PerResultsTable)
    if (export) {
        utils::write.csv(
            PerResultsTable,
            file = paste0(ExpName,"-PerResultsTable.csv")
        )
        cat(
            "Results are saved as: ", paste0(ExpName,"-PerResultsTable.csv" ),
            "in your directory [use getwd()]\n"
        )
    }
    return(object)
}



#' @title Directionality Table
#' @description Directionality Ratio is the displacement divided by
#' the total length of the total path distance, where displacement is
#' the straight line length between the start point and the endpoint
#' of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param export if `TRUE` (default), exports function output to
#' CSV file
#' @param ExpName string
#'
#' @return An CellMig class object with a data frame stored in the
#' DRtable slot. It contains nine rows: "Cell Number",
#' "Directionality Ratio","Mean Cumulative Directionality Ratio",
#' "Stable Directionality Ratio", "Number of returns","Min CumDR",
#' "Location of Min CumDR, Steps with less CumDR than DR",
#' "Directional Persistence"
#'
#'
#' @details    Directionality Ratio and Directional persistence
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' rmTD <- DiRatio(rmTD, export=FALSE)
#'
#'
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
DiRatio = function(object,TimeInterval=10, export=FALSE, ExpName=NULL) {
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    Step<-length(Object[[1]][,1])
    DRResultsTable<-data.frame()
    DIR.RATIO<-c()
    # calculating  cumsum of distance for each cell
    for(j in seq(1,length(Object),by=1)){
        MM<-Step
        MM2<-MM-1
        # finding the cordinates of the final point in the track.
        end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])

        start<-cbind(0,0)
        # calculating final distance
        final.dis=SpatialTools::dist2(start, end)
        # calculating the total distance
        total.dis=sum(Object[[j]][seq(1,MM2,by=1),6])
        # calc the Dir.Ratio
        Dir.Ratio=round(final.dis/total.dis ,digits = 3)
        mean.Dir.Ratio<-round(
            mean(Object[[j]][seq(1,MM2,by=1),13],na.rm = TRUE) ,digits = 3
        )
        StableDR<- (1-(mean.Dir.Ratio-Dir.Ratio))* mean.Dir.Ratio
        StableDR<-round(StableDR,digits = 3)
        DRResultsTable[1,j]<- j
        DRResultsTable[2,j]<-Dir.Ratio
        DRResultsTable[3,j]<-mean.Dir.Ratio
        DRResultsTable[4,j]<- StableDR
    }
    # Adding min CumDR    and number of angles greater than .75
    for(j in seq(1,length(Object),by=1)){
        MM<-Step
        MM2<-MM-1
        p1<-Object[[j]][seq(1,MM2,by=1),9]
        returns<-subset(p1,p1<(-0.87)) # greater than 150 degrees
        DRResultsTable[5,j]<-length(returns)
        p2<-Object[[j]][seq(1,MM2,by=1),13]
        w<-which.min(p2)
        DRResultsTable[6,j]<-round(p2[w], digits=3)
        DRResultsTable[7,j]<-w
        DR<-as.numeric(DRResultsTable[2,j])
        # number of the steps that have a cumdr less than the final dr
        lessThanMINcumdr<-subset(p2,p2<DR)
        DRResultsTable[8,j]<-length(lessThanMINcumdr)
        Ptime<-Object[[j]][seq(1,MM2,by=1),10]
        # computing the number of persistence steps to be used in computing
        # the persistence ratio
        PerTimLen<-Ptime
        PerTimLen[PerTimLen==0]<-NA
        PerTimLen<-PerTimLen[!is.na(PerTimLen)]
        PerTimLen<-length(PerTimLen)
        PerRatio<-round(PerTimLen/MM2, digits=2)
        DRResultsTable[9,j]<- PerRatio + DRResultsTable[4,j]
    }
    rownames(DRResultsTable)<-c(
        "Cell Number","Directionality Ratio",
        "Mean Cumulative Directionality Ratio","Stable Directionality Ratio",
        "Number of returns","Min CumDR",
        paste0("Location of Min CumDR (out of ",Step-1,")"),
        "Steps with less CumDR than DR","Directional Persistence"
    )
    RM1<-round(
        matrixStats::rowMedians(as.matrix(DRResultsTable),na.rm = TRUE),
        digits=3
    )
    DRResultsTable[,(length(Object)+1)]<-RM1
    DRResultsTable[1,(length(Object)+1)]<-"All Cells"
    object <- setCellMigSlot(object, "DRtable", DRResultsTable)
    if (export) {
        utils::write.csv(
            DRResultsTable,
            file = paste0(ExpName,"-DRResultsTable.csv")
        )
        cat("Results are saved as: ",
            paste0(ExpName,"-DRResultsTable.csv" ),
            "in your directory [use getwd()]","\n")
    }
    return(object)
}



#' @title Directionality Ratio plots
#' @description Directionality Ratio is the displacement divided
#' by the total length of the total path distance, where displacement
#' is the straightline length between the start point and the endpoint
#' of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param export if `TRUE` (default), exports plot to JPG file
#' @param ExpName string, name of the experiment. Can be NULL
#' @return Directionality Ratio plots
#'
#' @details    Directionality Ratio
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' DiRatioPlot(object=rmTD, export=FALSE)
#'
#'
#' @importFrom grDevices rainbow jpeg dev.off rgb
#' @importFrom graphics plot axis title lines polygon
#' @importFrom matrixStats rowMedians rowSds
#'
#'
#' @export
DiRatioPlot = function(object,TimeInterval=10, export=FALSE,ExpName = NULL) {
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    Step<-length(Object[[1]][,1])
    color <-c()
    ExpName<-getCellMigSlot(object, "meta")$expName
    Len<-length(Object)
    if (Len> 1023){
        colnum= Len-1023
        color1 <-grDevices::rainbow(1023)
        colo2 <-grDevices::rainbow(colnum)
        color=c(color1 ,colo2)
    }else{
        color <-grDevices::rainbow(Len)
    }
    if (export) {
        new.fld <-paste0(ExpName,"-DR_Results")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }
        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
    }
    DIR.RATIO.AllCells<-data.frame()
    for(j in seq(1,length(Object),by=1)){
        MM<-Step
        MM2<-MM-1
        Time<-c(seq(1,MM2,by=1))
        t<-60/TimeInterval
        xM<-MM2*TimeInterval
        xMM<-round(xM/60)
        xMMM<-seq(0,xMM,by=1)
        xMMM1<-(c(seq(0,xMM,by=1) )*(60/TimeInterval))
        Xaxis<-seq(1,MM2,by=1)
        if (export){
            plot_name <-    paste0(ExpName,"-D.R.plot-",j,".jpg")
            file_path <- file.path(new.fld, plot_name)
            grDevices::jpeg(
                filename = file_path,width = 4, height = 4, units = 'in',
                res = 300
            )
        }
        p<-graphics::plot(
            Time,Object[[j]][seq(1,MM2,by=1),13], type="l",col=color[j],
            xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1
        )
        graphics::axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
        graphics::title(main=paste0("Cell Number    ", j),col.main=color[j])
        grDevices::dev.off()
        DIR.RATIO.AllCells[seq(1,MM2,by=1),j]<-Object[[j]][seq(1,MM2,by=1),13]
    }
    mycolblue <- grDevices::rgb(
        0, 0, 255, maxColorValue = 255, alpha = 100, names = "blue"
    ) #transparent color
    mean.DIR.RATIO.AllCells<-matrixStats::rowMedians(
        as.matrix(DIR.RATIO.AllCells),na.rm = TRUE
    )
    SD.DIR.RATIO.AllCells<-matrixStats::rowSds(
        as.matrix(DIR.RATIO.AllCells),na.rm = TRUE
    )        # is a function in matrixStats
    meanSDp<-mean.DIR.RATIO.AllCells+SD.DIR.RATIO.AllCells
    meanSDp=ifelse(meanSDp>=1,1,meanSDp)
    meanSDn<-mean.DIR.RATIO.AllCells-SD.DIR.RATIO.AllCells
    meanSDpn<-c(meanSDp,meanSDn)
    if (export) {
        plot_name <-    paste0(ExpName,"directionality ratio for all cells.jpg")
        file_path <- file.path(new.fld, plot_name)
        grDevices::jpeg(
            filename = file_path,width = 4, height = 4, units = 'in', res = 300
        )
    }
    p<-graphics::plot(
        Time,mean.DIR.RATIO.AllCells, type="l",col="black",
        xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1
    )
    graphics::axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
    graphics::lines(Time,meanSDp, col="black")
    graphics::lines(Time,meanSDn, col="black")
    graphics::polygon(
        c(Time, rev(Time)), c(meanSDp, rev(meanSDn)),col = mycolblue ,
        border = NA
    )
    graphics::title(main="Directionality Ratio - All Cells",col.main="black")
    if (export) {
        grDevices::dev.off()
        cat("Plots are saved in a folder in your directory [use getwd()]","\n")
    }
}



#' @title Mean Square Displacement
#' @description The MSD function automatically computes the mean
#' square displacements across several sequential time intervals.
#' MSD parameters are used to assess the area explored by cells over time.
#'
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags
#' for the slope fitting. Default is 0.25, which represents 25 percent
#' of the steps.
#' @param ffLAG A numeric value to be used to get the number of lags
#' for the    Furth formula fitting. Default is 0.25, which represents
#' 25 percent of the steps.
#' @param SlopePlot A logical vector that allows generating individual
#' plots showing the slope of the mean square displacement of the
#' movement of individual cells. Default is TRUE.
#' @param AllSlopesPlot A logical vector that allows generating a plot
#' showing the slope of the mean square displacement of the movement of
#' all cells. Default is TRUE.
#' @param FurthPlot A logical vector that allows generating individual
#' plots fitting the Furth formula using generalized regression by the
#' Nelder–Mead method simplex method per cell. Default is TRUE.
#' @param AllFurthPlot A logical vector that allows generating a plot
#' fitting the Furth formula using generalized regression by the
#' Nelder–Mead method simplex method for all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#' @param ExpName string, anem of the Experiment. Can be NULL
#'
#' @return An CellMig class object with a data frame and plots.
#' The data frame is stored in the MSDtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF <- TrajectoryDataset[seq(1,220,by=1), ]
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
MSD <- function(
    object, TimeInterval=10, sLAG=0.25, ffLAG=0.25, SlopePlot=TRUE,
    AllSlopesPlot=TRUE, FurthPlot=TRUE, AllFurthPlot=TRUE,
    export=FALSE, ExpName=NULL) {
    if (!is.numeric(TimeInterval) | TimeInterval<= 0) {
        stop( "TimeInterval has to be a positive number")
    }
    if (!is.numeric(sLAG) | sLAG<= 0) {
        stop( "sLAG has to be a positive number")
    }
    if (!is.numeric(ffLAG) | ffLAG<= 0) {
        stop( "ffLAG has to be a positive number")
    }
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    Object <- getCellMigSlot(object, "preprocessedDS")
    ExpName <- getCellMigSlot(object, "meta")$expName
    if (export) {
        new.fld <-paste0(ExpName,"-MSDResults")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }
        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
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
    MSDResultsTable<-fixMSD(
        Object,Step=Step,SlopePlot=SlopePlot,AllSlopesPlot=AllSlopesPlot,
        FurthPlot=FurthPlot,AllFurthPlot=AllFurthPlot,sLAG=sLAG,ffLAG=ffLAG,
        color=color,export=export, ExpName = ExpName, new.fld = new.fld)
    rownames(MSDResultsTable)<-c(
        "Cell Number","MSD (lag=1)", "MSD slope", "N-M best fit (Furth) [D]",
        "N-M best fit (Furth) [P]","The significance of fitting D",
        "The significance of fitting P"
    )
    object <- setCellMigSlot(object, "MSDtable", MSDResultsTable)
    if (export) {
        utils::write.csv(
            MSDResultsTable,
            file = paste0(ExpName,"-MSDResultsTable.csv"))
        message("Results were saved in your directory [use getwd()]","\n")
    }
    return(object)
}





#' @title Direction AutoCorrelation
#'
#' @description The DiAutoCor function automatically compute the
#' angular persistence across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags
#' for the slope fitting. Default is 0.25, which represents 25
#' percent of the steps.
#' @param sPLOT A logical vector that allows generating individual
#' plots showing the angular persistence across several sequantial
#' time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot
#' showing the angular persistence across several sequantial time
#' intervals of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV
#' file
#' @param ExpName, string, name of the experiment. Can be NULL
#'
#' @return An CellMig class Object with a data frame and plots.
#' The data frame, which contains six rows: "Cell Number",
#' "Angular Persistence", "Intercept of DA quadratic model",
#' "Mean Direction AutoCorrelation (all lags)",
#' "Stable Direction AutoCorrelation through the track" and
#' "Difference between Mean DA and Intercept DA".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#'
#'
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[seq(1,220,by=1),]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=55)
#' rmTD <- DiAutoCor(rmTD, TimeInterval=10, sLAG=0.25, sPLOT=FALSE,
#'                   aPLOT=FALSE, export=FALSE)
#'
#' @importFrom grDevices rainbow
#' @importFrom stats lm predict median
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics plot lines title abline legend
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
DiAutoCor= function(
    object, TimeInterval=10, sLAG=0.25, sPLOT=TRUE, aPLOT=TRUE,
    export=FALSE, ExpName=NULL) {
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    if ( ! is.numeric(sLAG) ) {
        stop( "sLAG has to be a positive number" )
    } else if ( sLAG<= 0 ) {
        stop( "sLAG has to be a positive number" )
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    if ( ! is.list(Object) ){
        stop(
            "Input data must be a list. ",
            "Please run the PreProcessing step first either ",
            "rmPreProcessing() or wsaPreProcessing()")
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    if (export) {
        new.fld <-paste0(ExpName,"-DIAutoCorResults")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }
        if(!dir.exists(new.fld)) { dir.create(new.fld) }
    }
    Len<-length(Object)
    Step<-length(Object[[1]][,1])
    color <-c()
    if (Len> 1023){
        colnum= Len-1023
        color1 <-grDevices::rainbow(1023)
        colo2 <-grDevices::rainbow(colnum)
        color=c(color1 ,colo2)
    } else {
        color <- grDevices::rainbow(Len)
    }
    DA.ResultsTable<-fixDA(
        Object, Step=Step, sLAG=sLAG, sPLOT=sPLOT, aPLOT=aPLOT, color=color,
        export=export, ExpName = ExpName, new.fld = new.fld)
    rownames(DA.ResultsTable)<-c(
        "Cell Number","Angular Persistence","Intercept of DA quadratic model",
        "Mean Direction AutoCorrelation (all lags)",
        "Stable Direction AutoCorrelation through the track",
        "Difference between Mean DA and Intercept DA"
    )
    object <- setCellMigSlot(object, "DACtable", DA.ResultsTable)
    if (export) {
        utils::write.csv(
            DA.ResultsTable, file = paste0(ExpName,"-DA.ResultsTable.csv"))
        cat("Results are saved in your directory [use getwd()]","\n")
    }
    return(object)
}




#' @title Velocity AutoCorrelation
#'
#' @description The VeAutoCor function automatically compute the
#' changes in both speed and direction across several sequantial time
#' intervals.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags
#' for the slope fitting. Default is 0.25, which represents 25
#' percent of the steps.
#' @param sPLOT A logical vector that allows generating individual
#' plots showing the velocity across several sequantial time
#' intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot
#' showing the velocity across several sequantial time intervals
#' of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to
#' CSV file
#' @param ExpName string, name of the experiment. Can be NULL
#'
#' @return Plots and a data frame, which contains six rows:
#' "Cell Number", "Velocity AutoCorrelation (lag=1)",
#' "2nd normalized Velocity AutoCorrelation",
#' "Intercept of VA quadratic model",
#' "Mean Velocity AutoCorrelation (all lags)", "Mean |Acceleration|"
#' and "Average Speed".

#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[1:300,]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=55)
#' rmTD <- VeAutoCor(rmTD, TimeInterval=10, sLAG=0.25, sPLOT=FALSE,
#'                   aPLOT=FALSE, export=FALSE)
#'
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom stats lm predict median
#' @importFrom graphics plot lines title abline
#' @importFrom utils write.csv
#'
#' @export
VeAutoCor = function(
    object, TimeInterval=10, sLAG=0.25, sPLOT=TRUE, aPLOT=TRUE,
    export=FALSE, ExpName=NULL) {
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    if ( ! is.numeric(sLAG) ) {
        stop( "sLAG has to be a positive number" )
    } else if ( sLAG<= 0 ) {
        stop( "sLAG has to be a positive number" )
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(msg, paste(
            "Input data must be a list. Please run the PreProcessing step ",
            "first, either using rmPreProcessing() or wsaPreProcessing()")
        )
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    if (export) {
        new.fld <-paste0(ExpName,"-VeAutoCorResults")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }
        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
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
    # creating a table that has all the VAC to be able to compute the mean
    VAC.table<-data.frame()
    # to save the fist VAC non-normalized value for all cells
    VAC.first.value<-c()
    # to save the second VAC non-normalized value for all cells
    VAC.second.value<-c()
    for(j in seq_along(Object)){
        meanVAC<-c()
        LAG<-round(Step*sLAG) #taking only the first 12.5%
        for(lag in seq_len(LAG)){
            # starting from 2 to exclude the first cosine which is always 1.
            res <- t(vapply(seq_len(Step - 1), function(i){
                Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2] # newdx
                Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3] # newdy
                return(as.numeric(Object[[j]][i,c(14,15)]))
            }, FUN.VALUE = numeric(2)))
            Object[[j]][seq(1,(Step -1),by=1),c(14,15)] <- res
            Object[[j]][,c(14,15)] <- lapply(Object[[j]][,c(14,15)], as.numeric)
            Object[[j]][,14][is.na(Object[[j]][,14])] <- 0 # replace NA with 0
            Object[[j]][,15][is.na(Object[[j]][,15])] <- 0 # replace NA with 0
            # starting from 2 to exclude the first cosine which is always 1.
            res1 <- vapply(seq_len(Step - lag), function(i){
                (
                    (Object[[j]][i,14]* Object[[j]][i+lag,14]) +
                        (Object[[j]][i,15]* Object[[j]][i+lag,15])
                ) / ((lag*TimeInterval)^2)
            }, FUN.VALUE = numeric(1))
            Object[[j]][seq(1,(Step - lag),by=1),23] <- res1
            meanVAC[lag]<-mean(Object[[j]][seq(1,(Step - lag),by=1),23])
        }
        VAC.first.value[j]<-meanVAC[1]
        NORMmeanVAC<-meanVAC/meanVAC[1]
        VAC.table[seq(1,LAG,by=1),j]<-meanVAC
        assign(paste0("VAC.Cell.",j),meanVAC)
        VAC.second.value[j]<-NORMmeanVAC[2]
        VA.ResultsTable[1,j]<-j
        VA.ResultsTable[2,j]<-round(meanVAC[1],digits=3) # VA (lag =1)
        VA.ResultsTable[3,j]<-round(NORMmeanVAC[2],digits=3) # VA (lag =2)
        lags<-c(seq(1,length(VAC.table[,1]),by=1))
        lags2<- lags^2
        quadratic.model<-c()
        quadratic.m <-stats::lm(VAC.table[,j]~ lags + lags2)
        c<-quadratic.m
        cc<-unlist(c)
        quadratic.model[j]<-cc[1]
        ccc<-unlist(cc[1])
        VA.ResultsTable[4,j]<-round(ccc,digits=3)
        VA.ResultsTable[5,j]<-round(mean(meanVAC),digits=3) # mean VA (all lags)
        timevalues <- seq(1, length(lags), 1)
        predictedcounts <- stats::predict(
            quadratic.m,list(Time=timevalues, Time2=timevalues^2)
        )
        if (sPLOT == TRUE){
            Xaxis<-seq(1,LAG,by=1)
            Yaxis<-meanVAC
            if (export) {
                plot_name <-    paste0(
                    ExpName,"Velocity Autocorrelation.plot.Cell",j,".jpg"
                )
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4, height = 4, units = 'in',
                    res = 300)
            }
            graphics::plot(
                Xaxis,Yaxis, type="o",ylim=range(meanVAC),
                xlim=c(0,lag),col=color[j],
                xlab="Lag",ylab="Velocity    Autocorrelation",
                pch=19,las=1,cex=1.2
            )
            graphics::lines(timevalues, predictedcounts, col = "black", lwd = 3)
            graphics::title(
                main=paste0("Cell Number    ", j, "     VA quadratic model"),
                col.main="darkgreen",
                cex.main=0.8,
                sub=paste0(
                    " Intercept of VA quadratic model = ",round(ccc, digits=3)
                ),col.sub="red"
            )
            if (export) grDevices::dev.off()
        }
        Object[[j]][,c(15,16)]=0
    }
    RM1<-matrixStats::rowMedians(as.matrix(VAC.table),na.rm = TRUE)
    VA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
    VA.ResultsTable[2,(length(Object)+1)]<-round(
        median(VAC.first.value),digits=3
    )
    VA.ResultsTable[3,(length(Object)+1)]<-round(
        median(VAC.second.value),digits=3
    )
    lags<-seq(1,length(VAC.table[,1]),by=1)
    lags2<- lags^2
    quadratic.model<-c()
    quadratic.m <-lm(RM1~ lags + lags2)
    c<-quadratic.m
    cc<-unlist(c)
    quadratic.model[j]<-cc[1]
    VA.ResultsTable[4,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
    VA.ResultsTable[5,(length(Object)+1)]<-round(
        stats::median(as.numeric(
            VA.ResultsTable[5,seq(1,length(Object),by=1)])), digits=3)
    ccc<-unlist(cc[1])
    timevalues <- seq(1, length(lags), 1)
    predictedcounts <- predict(
        quadratic.m,list(Time=timevalues, Time2=timevalues^2)
    )
    if ( aPLOT == TRUE){
        Xaxis<-seq(1,LAG,by=1)
        Yaxis<-RM1
        if (export) {
            plot_name <- paste0(
                ExpName,"-Velocity.Autocorrelation.All.Cells.jpg"
            )
            file_path <- file.path(new.fld, plot_name)
            grDevices::jpeg(
                filename = file_path,width = 4,
                height = 4, units = 'in', res = 300)
        }
        graphics::plot(
            Xaxis,Yaxis, type="o",ylim=range(RM1),xlim=c(0,lag),col="black",
            xlab="Lag",ylab="Velocity    Autocorrelation",pch=19,las=1,cex=1.2
        )
        graphics::lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
        MVA<-round(
            stats::median(as.numeric(
                VA.ResultsTable[5,seq(1,length(Object),by=1)])),
            digits=2)
        graphics::abline(h=MVA,col="blue",lwd = 2)
        graphics::title(
            main=paste0("All Cells - VA quadratic model"), col.main="darkgreen",
            cex.main=0.8, sub=paste0(
                " Intercept of VA quadratic model = ",
                round(ccc, digits=3)),col.sub="red")
        if (export) grDevices::dev.off()
    }
    rownames(VA.ResultsTable)<-c(
        "Cell Number","Velocity AutoCorrelation (lag=1)",
        "2nd normalized Velocity AutoCorrelation",
        "Intercept of VA quadratic model",
        "Mean Velocity AutoCorrelation (all lags)"
    )
    object <- setCellMigSlot(object, "VACtable", VA.ResultsTable)
    if (export) {
        utils::write.csv(
            VA.ResultsTable, file = paste0(ExpName,"-VA.ResultsTable.csv")
        )
        cat("Results are saved in your directory [use getwd()]\n")
    }
    return(object)
}


#' @title Forward Migration
#'
#' @description The ForwardMigration function automatically generates
#' data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param sfptPLOT A logical vector that allows generating individual
#' plots of persistence time vs speed per cell. Default is TRUE.
#' @param afptPLOT    A logical vector that allows generating a plot of
#' persistence time vs speed for all cells. Default is TRUE.
#' @param sfpPLOT A logical vector that allows generating individual
#' plots of angular persistence vs speed per cell. Default is TRUE.
#' @param afpPLOT A logical vector that allows generating a plot of
#' angular persistence vs speed of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV
#' file
#' @param ExpName string, name of the experiment. Can be NULL
#'
#' @return    An CellMig class Object with a data frame and plots.
#' The data frame is stored in the ForMigtable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wsaDF <- WSADataset[seq(1,500,by=1),]
#' wsaTD <- CellMig(wsaDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-ForwardMigration(wsaTD, TimeInterval=10, sfptPLOT=FALSE,
#'                          afptPLOT= FALSE,sfpPLOT= FALSE,
#'                          afpPLOT= FALSE, export=FALSE)
#'
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom stats cor.test lm
#' @importFrom graphics plot title abline
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
ForwardMigration <- function(
    object, TimeInterval = 10, sfptPLOT = TRUE, afptPLOT = TRUE,
    sfpPLOT = TRUE, afpPLOT = TRUE, export = FALSE, ExpName = NULL){
    if (!is.numeric(TimeInterval)) {
        stop("TimeInterval has to be a positive number")
    } else if (TimeInterval <= 0) {
        stop("TimeInterval has to be a positive number")
    }
    Object <- getCellMigSlot(object, "preprocessedDS")
    UPorDO <- getCellMigSlot(object, "cellpos")
    msg        <- NULL
    if (!is.list(Object)) {
        msg <- c(
            msg, "Input data must be a list.",
            "Please run the PreProcessing step first either ",
            "rmPreProcessing() or wsaPreProcessing()"
        )
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    Len     <- length(Object)
    Step    <- length(Object[[1]][, 1])
    color <- c()
    if (Len > 1023) {
        colnum <- Len-1023
        color1 <- grDevices::rainbow(1023)
        colo2    <- grDevices::rainbow(colnum)
        color    <- c(color1, colo2)
    } else {
        color <- grDevices::rainbow(Len)
    }
    if (export) {
        new.fld <-paste0(ExpName, "-ForwardMigrationResults")
        if (dir.exists(new.fld)) {
            unlink(new.fld, recursive = TRUE, force = TRUE)
        }
        if(!dir.exists(new.fld)) {
            dir.create(new.fld)
        }
    }
    # Defining if the cell is ubove (1) the wound or below (0) it
    for (j in seq(1,length(Object),by=1)) {
        Object[[j]][, 25] <- UPorDO[j]
    }
    # creating values for    rel.ang.F    (step to the original)
    Object<-fixFM1(Object,Step=Step)
    # creating values for cosine based on rel.ang.F
    Object<-fixFM2(Object,Step=Step)
    cosine.FP<-fixFM3(Object,Step=Step)
    ## Forward Pesrsistence time, FP deviating time and FP ratio #
    Object<-fixFM4(Object,TimeInterval=TimeInterval,Step=Step)
    # Forward migration results
    FMResultsTable<-fixFM5(Object,TimeInterval=TimeInterval,Step=Step)

    # creating a table to store the mean velocity with
    # correspondence with the FP time
    FMResultsTable<-fixFM6(
        Object,FMResultsTable=FMResultsTable, Step=Step,sfptPLOT=sfptPLOT,
        afptPLOT=afptPLOT,export=export,color=color,sfpPLOT=sfpPLOT,
        TimeInterval=TimeInterval, ExpName = ExpName, new.fld = new.fld)
    RM<-round(
        matrixStats::rowMedians(as.matrix(
            cosine.FP[seq(1,(Step-1),by=1),]),na.rm = TRUE),
        digits=3
    )
    Speed<-data.frame()
    # calculating the Mean.Square.velocity for each cell
    for (j in seq(1,length(Object),by=1)){
        MM<-Step
        MM2<-MM-1
        Speed[seq(1,MM2,by=1),j]<-
            round(sqrt(Object[[j]][seq(1,MM2,by=1),11]),digits = 3)
    }
    RowmeanSpeed<-round(
        matrixStats::rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3
    )
    #testing the correlation
    s<-stats::cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    FMResultsTable[9,(length(Object)+1)]<-VEvsCOSP
    if ( afpPLOT == TRUE){
        if (export){
            plot_name <-    paste0(ExpName,"_All_Cells_FP_VS_Speed.jpg")
            file_path <- file.path(new.fld, plot_name)
            grDevices::jpeg(
                filename = file_path,width = 4, height = 4, units = 'in',
                res = 300
            )
        }
        graphics::plot(
            RowmeanSpeed*60,RM,pch=16,type="p",
            ylab="Forward Persistence Time (min)",
            xlab=" Instantaneous Speed (um/h)",col="black",las=1
        )
        reg<-stats::lm(RM~RowmeanSpeed)
        graphics::abline(reg,untf=FALSE,col="red")
        graphics::title(
            main=paste0("All Cells Speed vs Forward Persistence "),
            cex.main = 0.8,
            sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),
            col.sub="red"
        )
        if (export) grDevices::dev.off()
    }
    RM1<-round(
        matrixStats::rowMedians(as.matrix(FMResultsTable),na.rm = TRUE),
        digits=3
    )
    FMResultsTable[
        c(2,3,4,5,7,8), (length(Object)+1)]<-RM1[c(2,3,4,5,7,8)]
    FMResultsTable[1,(length(Object)+1)]<-"All Cells"
    rownames(FMResultsTable)<-c(
        "Cell Number","Mean Forward Persist Time (min)",
        "Mean Forward Persist Deviating Time (min)",
        "Forward Persistence Ratio",
        "Maximum Forward Persistence period",
        "Forward Persistence Time vs Speed (SCC)",
        "RMSS (um/h)","Mean Forward Angular Persistence (mean cos.F)",
        "Instantaneous Speed vs Forward Persistence (SCC)"
    )
    FMResultsTable<-FMResultsTable[-7,]
    object <- setCellMigSlot(object, "ForMigtable", FMResultsTable)
    if (export) {
        utils::write.csv(
            FMResultsTable,
            file = paste0(ExpName,"-FMResultsTable.csv")
        )
        cat("Results are saved as: ",
            paste0(ExpName,"-FMResultsTable.csv" ),
            "in your directory [use getwd()]\n")
    }
    return(object)
}







#' @title Forward Migration Index
#'
#' @description The FMI function automatically generates data for
#' the forward migration index
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param export if `TRUE` (default), exports function output to
#' CSV file
#' @param ExpName string, name of the experiment. Can be NULL
#'
#' @return    An CellMig class Object with a data frame. The data
#' frame is stored in the FMItable slot.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF=WSADataset[seq(1,300,by=1),]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-FMI(wsaTD,TimeInterval=10, export=FALSE)
#'
#' @importFrom grDevices rainbow jpeg dev.off
#' @importFrom SpatialTools dist2
#' @importFrom matrixStats rowMedians
#' @importFrom utils write.csv
#'
#' @export
FMI= function(object, TimeInterval=10, export=FALSE, ExpName=NULL){
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    if ( ! is.numeric(TimeInterval) ) {
        stop( "TimeInterval has to be a positive number" )
    } else if ( TimeInterval<= 0 ) {
        stop( "TimeInterval has to be a positive number" )
    }
    Object<-getCellMigSlot(object, "preprocessedDS")
    UPorDO<-getCellMigSlot(object, "cellpos")
    msg <- NULL
    if ( ! is.list(Object) ){
        msg <- c(
            msg,
            paste(
                "Input data must be a list.",
                "Please run the PreProcessing step first",
                "either rmPreProcessing() or wsaPreProcessing()"
            )
        )
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
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
    # Defining if the cell is ubove (1) the wound or below (0) it
    for(j in seq(1,length(Object),by=1)){
        Object[[j]][,25]<-UPorDO[j]
    }
    FMIResultsTable<-data.frame()
    # calculating the cumsum of distance for each cell
    for(j in seq_along(Object)){
        MM<-Step
        MM1<-MM-1
        res <- vapply(seq_len(MM), function(i){
            if (
                (Object[[j]][1,25]==0) &&
                (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) &&
                (Object[[j]][i,5]<0)
            ){     # upper cell going up or lower cell going down
                Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
            }
            if (
                (Object[[j]][1,25]==0) &&
                (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) &&
                (Object[[j]][i,5]>0)
            ){
                Object[[j]][i,20]<-cos(Object[[j]][i,19])
            }
            return(Object[[j]][i,20])
        }, FUN.VALUE = numeric(1))
        Object[[j]][seq(1,MM,by=1), 20] <- as.data.frame(res)
        # finding the cordinates of the final point in the track.
        end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])
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
        # finding the cordinates of the final point in the mid track.
        mid<-cbind(Object[[j]][MMM,2],Object[[j]][MMM,3])
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
        returns<-subset(p1,p1<(-0.87))            # greater than 150 degrees
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
    # creating values for    rel.ang.F    (step to the original)
    for(j in seq(1,length(Object),by=1)){
        MM<-Step
        if(
            (Object[[j]][1,25]==0) &&
            (Object[[j]][MM,3]>0) || (Object[[j]][1,25]==1) &&
            (Object[[j]][MM,3]<0)
        ){
            FMIResultsTable[6,j]<- (-1) * FMIResultsTable[6,j]
        } else {
            FMIResultsTable[6,j]<- FMIResultsTable[6,j]
        }
    }
    RM1<-round(
        matrixStats::rowMedians(as.matrix(FMIResultsTable),na.rm = TRUE),
        digits=3
    )
    FMIResultsTable[,(length(Object)+1)]<-RM1
    FMIResultsTable[6,(length(Object)+1)]<-round(
        FMIResultsTable[6,(length(Object)+1)]
    )
    FMIResultsTable[1,(length(Object)+1)]<-"All Cells"
    rownames(FMIResultsTable)<-c(
        "Cell Number","FMI","FMIy", "MTFMI","MTFMIy","Deepness (um)",
        "Number of backwards"
    )
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
    object <- setCellMigSlot(object, "FMItable", FMIResultsTable)
    return(object)
}



#' @title Final Results
#'
#' @description The FinRes function automatically generates a
#' data frame that contains all the results.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param ParCor A logical vector that allows generating a
#' correlation table. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#' to CSV file
#' @param ExpName string, name of the experiment. Can be NULL
#' @return    A data frame that contains all the results.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF <- WSADataset[seq(1,300,by=1), ]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-FMI(wsaTD,TimeInterval=10)
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,)
#' wsaTD <-FinRes(wsaTD,ParCor=FALSE, export=FALSE)
#'
#' @importFrom utils write.csv
#' @importFrom Hmisc rcorr
#'
#' @export
FinRes <- function(
    object, ParCor=TRUE, export=FALSE, ExpName=NULL)
{
    # ==========================================================================
    # Writing results
    # ==========================================================================
    if(!is.null(ExpName)) {
        ExpName <- as.character(ExpName)[1]
        ExpName <- fixExpName(ExpName)
        object <- setExpName(object, ExpName)
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    if (length(getCellMigSlot(object, "results")) > 0)
        object <- setCellMigSlot(object, "results", data.frame() ) # reset
    juxtaposeResults <- function(slt, obj=object) {
        # rbinds slt to object@results, removing the first row of slt
        new_results <- slot(obj, slt)
        tot_results <- slot(obj, "results")
        if (length(new_results) > 0) {
            new_results_1st_row <- new_results[1, ]
            new_results_no_1st_row <- new_results[-1, ]
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
    for (i in partial_result_slots)
        object <- setCellMigSlot(object, "results", juxtaposeResults(i) )
    # ==========================================================================
    # Exporting results
    # ==========================================================================
    if (export) {
        utils::write.csv(
            getCellMigSlot(object, "results"),
            file = paste0(ExpName,"-Final_Results.csv")
        )
        cat(
            "The table of the final results is saved as: ",
            paste0(ExpName, "-Final_Results.csv"),
            " in your directory [use getwd()]","\n"
        )
    }
    # ==========================================================================
    # Calculating correlation table
    # ==========================================================================
    if (ParCor) {
        R <- getCellMigSlot(object, "results")
        R[, ] <- lapply(R[, ], function(x) as.numeric(gsub(",", ".", x)))
        Parameters.Correlation <- Hmisc::rcorr(t(R), type="spearman")
        object <- setCellMigSlot(object, "parCor", Parameters.Correlation$r)
        if (export) {
            utils::write.csv(
                Parameters.Correlation$r,
                file = paste0(ExpName,"-Parameters.Correlation.csv")
            )
            cat(
                "Parameters Correlation table is saved as: ",
                paste0(ExpName, "-Parameters.Correlation.csv"),
                "in your directory [use getwd()]","\n"
            )
        }
    }
    message("\nThese are the parameters in your final results:")
    print(rownames(getCellMigSlot(object, "results")))
    # ==========================================================================
    # Returning whole object
    # ==========================================================================
    return(object)
}


#' @title PCA
#'
#' @description The CellMigPCA function automatically generates
#' Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @param parameters A numeric vector contains the parameters to
#' be included in the Principal Component Analysis. These numbers
#' can be obtained from the outcome of the FinRes() function.
#'
#' @return    PCA Graph of cells and PCA Graph of variables.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF=WSADataset[seq(1,300,by=1),]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-FMI(wsaTD,TimeInterval=10)
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10)
#' wsaTD <-FinRes(wsaTD,ParCor=FALSE)
#' PCAplot<-CellMigPCA(wsaTD,parameters=c(1,4))
#'
#' @importFrom FactoMineR PCA
#'
#' @export
CellMigPCA = function(object, parameters=c(1,2,3)){
    if (!is.list(object) & !is(object, "CellMig")) {
        stop(
            "Input data must be a list. ",
            "Please run the PreProcessing step first, ",
            "either rmPreProcessing() or wsaPreProcessing()"
        )
    }
    if ( length(getCellMigSlot(object, "results")[,1])<1 ){
        stop(
            "There are no results stored. ",
            "Please run trajectory analysis first"
        )
    }
    if ( length(parameters)<2){
        stop("At least two parameters are required to run the PCA")
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    df1<- getCellMigSlot(object, "results")
    #### excluding the last column since it is the avarage of all the cells
    df1=df1[,-length(getCellMigSlot(object, "results")[1,])]
    tt<-t(df1)
    tt1=tt[,parameters]
    res <- FactoMineR::PCA(tt1)
}


#' @title PCA Clusters
#'
#' @description The CellMigPCAclust function automatically generates
#' clusters based on the Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param parameters A numeric vector contains the parameters to be
#' included in the Principal Component Analysis. These numbers can
#' be obtained from the outcome of the FinRes() function.
#' @param export if `TRUE` (default), exports function output to
#' CSV file
#' @param interactive logical, shall the PCA analysis be generated in
#' a interactive fashion
#'
#' @return    PCA Graph of cells and PCA Graph of variables.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' ## The analysis only supports the interactive method!
#' ## If interactive=FALSE, the function will return NULL
#' data(WSADataset)
#' wasDF <- WSADataset[seq(1, 300, by=1), ]
#' wsaTD <- CellMig(wasDF)
#' CellMigPCAclust(wsaTD, parameters=c(1,9), interactive=FALSE)
#' ##
#' ## A real world example is shown below (uncomment)
#' # data(WSADataset)
#' # wasDF <- WSADataset[seq(1,300,by=1),]
#' # wsaTD <- CellMig(wasDF)
#' # wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' # wsaTD <-FMI(wsaTD,TimeInterval=10)
#' # wsaTD <-ForwardMigration(wsaTD,TimeInterval=10)
#' # wsaTD <-FinRes(wsaTD,ParCor=FALSE)
#' # PCAclust <- CellMigPCAclust(wsaTD,parameters=c(1,9))
#'
#' @importFrom FactoMineR PCA
#'
#' @export
CellMigPCAclust = function(
    object, parameters=c(1,2,3), export=FALSE, interactive=TRUE){
    if (!is.list(object) & !is(object, "CellMig")) {
        stop(
            "Input data must be a list. ",
            "Please run the PreProcessing step first, ",
            "either rmPreProcessing() or wsaPreProcessing()"
        )
    }
    if (!interactive){
        return(NULL)
    }
    if ( length(getCellMigSlot(object, "results")[,1])<1 ){
        stop(
            "There are no results stored. Please run trajectory analysis first"
        )
    }
    if ( length(parameters)<2){
        stop("At least two parameters are required to run the PCA")
    }
    ExpName<-getCellMigSlot(object, "meta")$expName
    df1<- getCellMigSlot(object, "results")
    #### excluding the last column since it is the avarage of all the cells
    df1=df1[,-length(getCellMigSlot(object, "results")[1,])]
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


#' @title Aggregating the outcome of several experiments or conditions.
#' @description Aggregate two or more CellMig-class objects together.
#' Input objects must carry information of trajectory analyses
#' (otherwise an error will be raised).
#' All trajectory results form the different experiments/conditions
#' are returned in two data frames.
#' @param x \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param ... one or more CellMig-class object(s) where cells'
#' trajectories have already been analyzed.
#' @param export if `TRUE` (default), exports function output to CSV
#' file
#'
#' @return two data frames:
#' The first data frame shows the average of each parameter per
#' experiment/condition.
#' The second data frame shows the parameters of individual cells
#' of all experiments/conditions.
#' @details    The visualization shows centered trajectories where the
#' starting point of each track is located at the origin of the
#' coordinate system (X=0,Y=0).
#'
#' @author Damiano Fantini and
#' Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#'
#' @examples
#' data(WSADataset)
#' wasDF1 <- WSADataset[seq(1,300,by=1), ]
#' wsaTD1 <- CellMig(wasDF1)
#' wsaTD1 <- wsaPreProcessing(wsaTD1,FrameN=55)
#' wsaTD1 <-FMI(wsaTD1,TimeInterval=10)
#' wsaTD1 <-FinRes(wsaTD1,ParCor=FALSE, export=FALSE)
#' wasDF2 <- WSADataset[seq(500,700,by=1), ]
#' wsaTD2 <- CellMig(wasDF2)
#' wsaTD2 <- wsaPreProcessing(wsaTD2,FrameN=55)
#' wsaTD2 <-FMI(wsaTD2,TimeInterval=10)
#' wsaTD2 <-FinRes(wsaTD2,ParCor=FALSE, export=FALSE)
#' AGG<-aggregateFR(wsaTD1 ,wsaTD2 ,export=FALSE)
#'
#'
#' @importFrom utils write.csv
#'
#' @export

aggregateFR<-function(x, ..., export=FALSE){
    # Make sure that the objects belong to the correct class
    class.nm <- "CellMig"
    stopifnot(class.nm %in% class(x))
    # Read all extra args
    obj.list <- list(...)
    keep <- do.call("c", lapply(obj.list, function(z) {
        class.nm %in% class(z)}))
    stopifnot(sum(keep) > 0)
    if (sum(keep) < length(keep)) {
        message(
            paste0(
                "A number of n=",
                (length(keep) - sum(keep)),
                " objects were excluded from the analysis! ",
                "Wrong class!"
            )
        )
    }
    # Re-compose list
    full.list <- list()
    full.list[[1]] <- getCellMigSlot(x, "results")
    for (i in seq(1,length(keep),by=1)) {
        if(keep[i]) {
            full.list[[length(full.list) + 1]] <-
                getCellMigSlot(obj.list[[i]], "results")
        }
    }
    # how many labs do we need?
    if (length(full.list) > length(LETTERS)) {
        my.labels <- as.character(seq(1, length(full.list), by = 1))
        max.char <- (nchar(length(full.list)) + 1)
        my.prfx <-    max.char - nchar(my.labels)
        my.prfx <-    do.call("c", lapply(my.prfx, function(jj) {
            paste(rep("0", times = jj), collapse = "")}))
        my.labels <- paste0(my.prfx, my.labels)
        my.labels <- paste0("Obj", my.labels, "_")
    } else {
        my.labels <- paste0(LETTERS, "_")
    }
    # Initialize collectors,
    # we need 3
    all.msg <- list()
    all.glob <- list()
    all.indiv <- list()
    for (i in seq(1,length(full.list),by=1)) {
        obj0 <- full.list[[i]]
        if (! is.list(obj0) ){
            all.msg[[length(all.msg) + 1]] <- i
        }
        # replace colNames
        colnames(obj0) <- paste0(my.labels[i], colnames(obj0))
        colnames(obj0) <- gsub("[[:space:]]", "_", colnames(obj0))
        # Aggr Cells Results
        obALL <- data.frame(obj0[, ncol(obj0)], stringsAsFactors = FALSE)
        rownames(obALL) <- rownames(obj0)
        colnames(obALL) <- sub("_$", "", paste0("Group_", my.labels[i]))
        # Individual Cells Results
        obIN <- data.frame(obj0[, -c(ncol(obj0))], stringsAsFactors = FALSE)
        rownames(obIN) <- rownames(obj0)
        colnames(obIN) <- colnames(obj0)[-c(ncol(obj0))]
        # store in tmp list
        all.glob[[length(all.glob) + 1]] <- obALL
        all.indiv[[length(all.indiv) + 1]] <- obIN
    }
    # additional checks
    stopifnot(
        length(all.glob) > 1 &&
            sum(do.call("c", lapply(all.glob, is.data.frame))) > 1 &&
            sum(do.call("c", lapply(all.glob, nrow)) > 1) > 1
    )
    stopifnot(
        length(all.indiv) > 1 &&
            sum(do.call("c", lapply(all.indiv, is.data.frame))) > 1 &&
            sum(do.call("c", lapply(all.indiv, nrow)) > 1) > 1
    )
    # rearrange rows
    all.glob.rows <- table(do.call("c", lapply(all.glob, rownames)))
    all.glob.rows <- names(
        all.glob.rows[all.glob.rows == max(all.glob.rows, na.rm = TRUE)]
    )
    all.indiv.rows <- table(do.call("c", lapply(all.indiv, rownames)))
    all.indiv.rows <- names(
        all.indiv.rows[all.indiv.rows == max(all.indiv.rows, na.rm = TRUE)]
    )
    tmp.all.glob <- list()
    for (xx in all.glob) {
        if (sum(all.glob.rows %in% rownames(xx)) == length(all.glob.rows)) {
            tmp.a <- data.frame(xx[all.glob.rows, ])
            rownames(tmp.a) <- all.glob.rows
            colnames(tmp.a) <- colnames(xx)
            tmp.all.glob[[length(tmp.all.glob) + 1]] <- tmp.a
        }
    }
    tmp.indiv.glob <- list()
    for (xx in all.indiv) {
        if (sum(all.indiv.rows %in% rownames(xx)) == length(all.indiv.rows)) {
            tmp.a <- data.frame(xx[all.indiv.rows, ])
            rownames(tmp.a) <- all.indiv.rows
            colnames(tmp.a) <- colnames(xx)
            tmp.indiv.glob[[length(tmp.indiv.glob) + 1]] <- tmp.a
        }
    }
    # aggregate results
    AllDF <- do.call(cbind, tmp.all.glob)
    InDF <- do.call(cbind, tmp.indiv.glob)
    # Message (msg) aggregation... is this even used?
    if (length(all.msg) > 0) {
        all.msg <- do.call("c", all.msg)
        msg <- paste0(
            "The following objects had uncompatible format: ",
            paste(all.msg, collapse = ", "), ". Result slots ",
            "must be lists/data.frames.")
    } else {
        msg <- "All result slots were included!"
    }
    # Done... now, export if needed.
    if (export) {
        utils::write.csv(
            AllDF,file = "AllDF.csv")
        utils::write.csv(
            InDF,file = "InDF.csv")
        cat("Results are saved in your directory [use getwd()]\n")
    }
    my_list <- list(AllDF,InDF)
    return(my_list)
}




#' @title PCA Clusters of different conditions
#'
#' @description The CellMigPCAclust function automatically generates
#' clusters based on the Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended
#' to all exported tracks and statistics data.
#' @param parameters A numeric vector contains the parameters to be
#' included in the Principal Component Analysis. These numbers can be
#' obtained from the outcome of the FinRes() function.
#' @param export if `TRUE` (default), exports function output to CSV
#' file
#' @param interactive logical, shall the PCA analysis be generated in
#' a interactive fashion
#'
#' @return    PCA Graph of cells and PCA Graph of variables.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' ## The analysis only supports the interactive method!
#' ## If interactive=FALSE, the function will return NULL
#' data(WSADataset)
#' wasDF1 <- WSADataset[seq(1,300,by=1), ]
#' wsaTD1 <- CellMig(wasDF1)
#' wsaTD1 <- wsaPreProcessing(wsaTD1,FrameN=55)
#' wsaTD1 <-FMI(wsaTD1,TimeInterval=10)
#' wsaTD1 <-FinRes(wsaTD1,ParCor=FALSE, export=FALSE)
#' wasDF2 <- WSADataset[seq(500,700,by=1), ]
#' wsaTD2 <- CellMig(wasDF2)
#' wsaTD2 <- wsaPreProcessing(wsaTD2,FrameN=55)
#' wsaTD2 <-FMI(wsaTD2, TimeInterval=10)
#' wsaTD2 <-FinRes(wsaTD2, ParCor=FALSE, export=FALSE)
#' AGG <- aggregateFR(wsaTD1, wsaTD2, export=FALSE)
#' CellMigPCAclustALL(AGG,ExpName="Aggregated_Conditions",
#'                    parameters=c(1,6), export=FALSE, interactive=FALSE)
#' # The previous line returns NULL
#' # In an interactive session, try running the following command (uncomment!)
#' # CellMigPCAclustALL(AGG,ExpName="Aggregated_Conditions",
#' #                    parameters=c(1,6), export=FALSE)
#'
#'
#' @importFrom FactoMineR PCA
#' @importFrom utils write.csv
#' @export
CellMigPCAclustALL <- function(
    object, ExpName="PCA_Clusters", parameters=c(1,2,3),
    export=FALSE, interactive=TRUE
){
    obj <- object[[2]]
    msg <- NULL
    if ( ! is.data.frame(obj)){
        msg <- c(msg, "input data must be data.frame")
    }
    if ( length(parameters)<2){
        stop("At least two parameters are required to run the PCA")
    }
    if(!interactive) { return (NULL) }
    df1<- obj
    #### excluding the last column since it is the avarage of all the cells
    df1=df1[,-length(obj[1,])]
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
            file = paste0(ExpName,"-ClustersALL.csv")
        )
        cat(
            "The table of clusters of the aggregated data frame is saved as: ",
            paste0(ExpName,"-ClustersALL.csv"),
            " in your directory [use getwd()]\n"
        )
    }

}


#' @title Getting the Directionality Table
#' @description Directionality Ratio is the displacement divided by
#' the
#' total length of the total path distance, where displacement is the
#' straight line length between the start point and the endpoint of
#' the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @return A data frame. It contains nine rows: "Cell Number",
#' "Directionality Ratio","Mean Cumulative Directionality Ratio",
#' "Stable Directionality Ratio", "Number of returns","Min CumDR",
#' "Location of Min CumDR, Steps with less CumDR than DR",
#' "Directional Persistence".
#'
#' @details    Directionality Ratio and Directional persistence
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' rmTD <- DiRatio(rmTD, export=FALSE)
#' head(getDiRatio(rmTD))
#'
#' @export
getDiRatio <- function(object)
{
    if (length(getCellMigSlot(object, "DRtable")) < 1)
        stop( "Please run the DiRatio() first" )
    DiRatio<- getCellMigSlot(object, "DRtable")
    return(DiRatio)
}




#' @title Getting the table of Persistence and Speed.
#' @description The PerAndSpeed() generates data and plots for
#' persistence and speed.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @return A data frame of Persistence and Speed.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @examples
#' rmTD <- get(data(preProcCellMig))
#' rmTD <- PerAndSpeed(rmTD,TimeInterval=10, export=FALSE)
#' head(getPerAndSpeed(rmTD))
#'
#' @export
getPerAndSpeed <- function(object)
{
    if (length(getCellMigSlot(object, "PerAanSpeedtable")) < 1) {
        stop( "Please run the PerAndSpeed() first" )
    }
    PerAanSpeedtable <- getCellMigSlot(object, "PerAanSpeedtable")
    return(PerAanSpeedtable)
}

#' @title Getting the Mean Square Displacement
#' @description The MSD function automatically computes the
#' mean square displacements across several sequential time intervals.
#' MSD parameters are used to assess the area explored by cells over
#' time.
#'
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#'
#' @return A data frame of MSD values.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF <- TrajectoryDataset[seq(1,600,by=1), ]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=100)
#' rmTD <- MSD(rmTD, sLAG=0.25, ffLAG=0.25, export=FALSE)
#' head(getMSDtable(rmTD))
#'
#' @export
getMSDtable <- function(object)
{
    if (length(getCellMigSlot(object, "MSDtable")) < 1 )
        stop( "Please run the MSD() first" )
    MSDtable <- getCellMigSlot(object, "MSDtable")
    return(MSDtable)
}

#' @title Getting the Direction AutoCorrelation
#'
#' @description The DiAutoCor function automatically compute the
#' angular persistence across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @return A data frame which contains six rows: "Cell Number",
#' "Angular Persistence", "Intercept of DA quadratic model",
#' "Mean Direction AutoCorrelation (all lags)",
#' "Stable Direction AutoCorrelation through the track"
#' and "Difference between Mean DA and Intercept DA".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#'
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[seq(1,300,by=1),]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=55)
#' rmTD <- DiAutoCor(rmTD, TimeInterval=10, sLAG=0.25, sPLOT=FALSE,
#'                   aPLOT=FALSE, export=FALSE)
#' head(getDACtable(rmTD))
#'
#' @export
getDACtable <- function(object)
{
    if (length(getCellMigSlot(object, "DACtable")) < 1)
        stop( "Please run the DiAutoCor() first" )
    DACtable<- getCellMigSlot(object, "DACtable")
    return(DACtable)
}


#' @title Getting the Velocity AutoCorrelation
#'
#' @description The VeAutoCor function automatically compute the
#' changes in both speed and direction across several sequantial
#' time intervals.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @return A data frame, which contains six rows: "Cell Number",
#' "Velocity AutoCorrelation (lag=1)",
#' "2nd normalized Velocity AutoCorrelation",
#' "Intercept of VA quadratic model",
#' "Mean Velocity AutoCorrelation (all lags)", "Mean |Acceleration|"
#' and "Average Speed".

#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(TrajectoryDataset)
#' rmDF=TrajectoryDataset[seq(1,300,by=1),]
#' rmTD <- CellMig(rmDF)
#' rmTD <- rmPreProcessing(rmTD,FrameN=55)
#' rmTD <- VeAutoCor(rmTD, TimeInterval=10, sLAG=0.25, sPLOT=FALSE,
#'                   aPLOT=FALSE, export=FALSE)
#' head(getVACtable(rmTD))
#'
#' @export
getVACtable<- function(object)
{
    if (length(getCellMigSlot(object, "VACtable")) < 1)
        stop( "Please run the VeAutoCor() first" )
    VACtable<- getCellMigSlot(object, "VACtable")
    return(VACtable)
}


#' @title Getting the Forward Migration
#'
#' @description The ForwardMigration function automatically generates
#' data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @return    A data frame inclusing values of the forward migration
#' analysis.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wsaDF <- WSADataset[seq(1,300,by=1),]
#' wsaTD <- CellMig(wsaDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-ForwardMigration(wsaTD, TimeInterval=10, sfptPLOT=FALSE,
#'                          afptPLOT= FALSE, sfpPLOT= FALSE,
#'                          afpPLOT= FALSE, export=FALSE)
#' head(getForMigtable(wsaTD))
#'
#' @export
getForMigtable<- function(object)
{
    if (length(getCellMigSlot(object, "ForMigtable")) < 1) {
        stop( "Please run the ForwardMigration() first" )
    }
    ForMigtable <- getCellMigSlot(object, "ForMigtable")
    return(ForMigtable)
}





#' @title Getting the Forward Migration Index
#'
#' @description The FMI function automatically generates data
#' for the forward migration index
#' @param object \code{CellMig} class object, which is a list of
#' data frames resulted from the PreProcessing.
#' @return    A data frame for the FMI.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF=WSADataset[seq(1,300,by=1),]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-FMI(wsaTD,TimeInterval=10, export=FALSE)
#' head(getFMItable(wsaTD))
#'
#' @export
getFMItable<- function(object)
{
    if (length(getCellMigSlot(object, "FMItable")) < 1)
        stop( "Please run the FMI() first" )
    FMItable<- getCellMigSlot(object, "FMItable")
    return(FMItable)
}



#' @title Final Results
#'
#' @description The FinRes function automatically generates a data
#' frame that contains all the results.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @return    A data frame that contains all the results.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF <- WSADataset[seq(1,300,by=1), ]
#' wsaTD <- CellMig(wasDF)
#' wsaTD <- wsaPreProcessing(wsaTD,FrameN=55)
#' wsaTD <-FMI(wsaTD,TimeInterval=10)
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,)
#' wsaTD <-FinRes(wsaTD,ParCor=FALSE, export=FALSE)
#' head(getResults(wsaTD))
#'
#' @export
getResults<- function(object)
{
    if (length(getCellMigSlot(object, "results")) < 1)
        stop( "Please run the FinRes() first" )
    results<- getCellMigSlot(object, "results")
    return(results)
}


#' @title Get Available Aggregate Cell Metrics
#'
#' @description Retrieve a list of metrics computed for an aggregated
#' result object.
#' This getter function takes the output of aggregateFR() as input.
#' @param object list of length 2, returned by the aggregateFR()
#' function
#' @return    character vector listing all available metrics
#'
#' @author Damiano Fantini and
#' Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' data(WSADataset)
#' wasDF1 <- WSADataset[seq(1,300,by=1), ]
#' wsaTD1 <- CellMig(wasDF1)
#' wsaTD1 <- wsaPreProcessing(wsaTD1,FrameN=65)
#' wsaTD1 <- FMI(wsaTD1,TimeInterval=10)
#' wsaTD1 <- FinRes(wsaTD1,ParCor=FALSE, export=FALSE)
#' wasDF2 <- WSADataset[seq(1001,1300,by=1), ]
#' wsaTD2 <- CellMig(wasDF2)
#' wsaTD2 <- wsaPreProcessing(wsaTD2,FrameN=65)
#' wsaTD2 <-FMI(wsaTD2,TimeInterval=10)
#' wsaTD2 <-FinRes(wsaTD2,ParCor=FALSE, export=FALSE)
#' AGG <- aggregateFR(wsaTD1 ,wsaTD2 ,export=FALSE)
#' getAvailableAggrMetrics(AGG)
#'
#' @export
getAvailableAggrMetrics <- function(object)
{
    if (!is.list(object) || length(object) < 2) {
        stop ( "Wrong input! Is this the output of aggregateFR()?" )
    }
    results <- rownames(object[[2]])
    return(results)
}

