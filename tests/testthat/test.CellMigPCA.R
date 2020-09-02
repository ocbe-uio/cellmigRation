context("Loading, pre-processing, DiRatio, FinRes")
# ==============================================================================
# Importing and converting data
# ==============================================================================
df <- TrajectoryDataset
rmTD <- CellMig(df)
test_that("Data conversion", expect_output(print(rmTD), "An S4 CellMig object"))
# ==============================================================================
# Pre-processing
# ==============================================================================
# rmTD <- rmPreProcessing(rmTD, FrameN=69) # FIXME: too long! Leaking output!
# test_that("Pre-processing", expect_output(print(rmTD), "An S4 CellMig object"))
# ------------------------------------------------------------------------------
# Other operations
# ------------------------------------------------------------------------------
# rmTD_2 <- DiRatio(rmTD, export=FALSE)
# test_that("DiRatio", expect_output(print(rmTD), "An S4 CellMig object"))
# rmTD_3 <- FinRes(rmTD_2, ParCor=TRUE, export=FALSE)
# test_that("FinRes", {
#     expect_output(
#         object = str(rmTD_3, max.level = 1),
#         expected = "Formal class 'CellMig'"
#     )
# })