# To manually check these tests, run BiocGenerics:::testPackage("cellmigRation")

test_OptimizeParams_getOptimizedParams <- function() {
    x <- cellmigRation::TrackCellsDataset
    x <- OptimizeParams(
        tc_obj = x, lnoise_range = c(5,7,10),
        diameter_range = c(12,14,18),
        threshold_range = c(4,7), dryrun = TRUE
    )
    RUnit::checkEquals(
        getOptimizedParameters(x),
        list(lnoise=5, diameter=14, threshold=7)
    )
}

test_CellTracker_getTracks <- function() {
    x <- cellmigRation::TrackCellsDataset
    x <- CellTracker(x, dryrun=TRUE)
    RUnit::checkEquals(
        getTracks(x)[1, ],
        data.frame("cell.ID"=1, "X"=28.9505, "Y"=30.41373, "frame.ID"=1),
        tolerance=1e-4
    )
}

test_ComputeTracksStats_getCellStats <- function() {
    x <- cellmigRation::TrackCellsDataset
    x <- ComputeTracksStats(
        x, time_between_frames = 10,
        resolution_pixel_per_micron = 20
    )
    RUnit::checkEquals(
        getCellsStats(x)[1, ],
        data.frame(
            "Cell_Number"=1, "Speed"=8.878605, "Distance"=532.7163,
            "Displacement"=505.4009, "Persistence"=0.9487242,
            "Degrees"=6.110741, "YFMI"=-0.1627922, "XFMI"=0.9346531,
            "Y_displacement"=-86.72204, "X_displacement"=497.905, "Frames"=6
        ),
        tolerance=1e-4
    )
}

test_aggregateTrackedCells <- function() {
    x0 <- cellmigRation::TrackCellsDataset
    x0 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "CTRL")
    x1 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DMSO")
    x2 <- setCellsMeta(x0, experiment = "my_exp_01", condition = "DRUG")
    y <- aggregateTrackedCells(x0, x1, x2, meta_id_field = "condition")
    RUnit::checkEquals(
        y[1, ],
        data.frame(
            "new.ID"=1001, "X"=28.9505, "Y"=30.41373, "frame.ID"=1, "cell.ID"=1,
            "tiff_file"="cellmigRationDemo.tif", "experiment"="my_exp_01",
            "condition"="CTRL", replicate=NA
        ),
        tolerance=1e-4
    )
}

test_DiRatio <- function() {
    rmTD <- cellmigRation::preProcCellMig
    rmTD <- DiRatio(rmTD, export=FALSE)
    RUnit::checkTrue(is(rmTD, "CellMig"))
}
