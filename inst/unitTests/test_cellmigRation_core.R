test_module_01 <- function() {
    # NextOdd
    RUnit::checkEquals(cellmigRation:::NextOdd(6), 7)

    # circshift
    RUnit::checkEquals(cellmigRation:::circshift(c(1, 2, 3), -1), c(2, 3, 1))

    # AddDimension
    A <- cellmigRation:::AddDimension(x = cbind(c(1,2,3)), y = 1)
    B <- cbind(c(1, 1, 1), c(1, 2, 3))
    RUnit::checkEquals(A, B)

    # MakeHypercube
    A <- cellmigRation:::MakeHypercube(c(1,2), 2)
    B <- cbind(c(1, 2, 1, 2), c(1, 1, 2, 2))
    RUnit::checkEquals(sum(A), sum(B))


    # Add more...
    #checkTrue(is.na(divideBy(4, 0)))
    #checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}

test_module_02 <- function() {
    df <- get(data(TrajectoryDataset))
    rmTD <- CellMig(df)
    RUnit::checkTrue(is(rmTD, "CellMig"))

}

