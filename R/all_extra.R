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
#' computing cell tracks and stats
#' 
#' @details
#' This stack of Fluorescent Cell images was obtained by 
#' processing TIFF files via cellmigRation. 
#' Imaging experiments were performed as described by Ghannoum S 
#' et al (paper in preparation). Briefly, triple negative breast 
#' cancer BT549 cells were cultured in RPMI supplemented with 10% FCS 
#' and 1% penicillin/streptomycin. Cells were transduced with 
#' NucLight green lentivirus (Essen BioScience), and then sorted by  
#' fluorescence-activated cell sorting (FACS). GFP-positive cells 
#' were seeded at a 1:3 ratio with untransduced BT549 cells in 
#' 96-well image-lock plates (EssenBio) at a density of 1000 
#' total cells per well.
#' Cells were scanned at ten-minute intervals over 24h using an 
#' Incucyte S3 Live-Cell microscope (EssenBio) at 10x magnification 
#' and a Basler Ace 1920-155um camera with CMOS sensor. 
#' A small selection of cells and image stacks is included in
#' this dataset.
#' 
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
#'
#' @details
#' BT549 cell trajectories were computed using cellmigRation.  
#' Imaging experiments were performed as described by Ghannoum S 
#' et al (paper in preparation). Briefly, triple negative breast 
#' cancer BT549 cells were cultured in RPMI supplemented with 10% FCS 
#' and 1% penicillin/streptomycin. Cells were transduced with 
#' NucLight green lentivirus (Essen BioScience), and then sorted by  
#' fluorescence-activated cell sorting (FACS). GFP-positive cells 
#' were seeded at a 1:3 ratio with untransduced BT549 cells in 
#' 96-well image-lock plates (EssenBio) at a density of 1000 
#' total cells per well. Once cells reached the desired density,
#' a thin wound was introduced by scratching the cell monolayer. 
#' Next, cells were scanned at ten-minute 
#' intervals over 24h using an 
#' Incucyte S3 Live-Cell microscope (EssenBio) at 10x magnification 
#' and a Basler Ace 1920-155um camera with CMOS sensor. TIFF images 
#' were imported and processed using the cellmigRation library.
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
#' @details
#' BT549 cell trajectories were computed using cellmigRation.  
#' Imaging experiments were performed as described by Ghannoum S 
#' et al (paper in preparation). Briefly, triple negative breast 
#' cancer BT549 cells were cultured in RPMI supplemented with 10% FCS 
#' and 1% penicillin/streptomycin. Cells were transduced with 
#' NucLight green lentivirus (Essen BioScience), and then sorted by  
#' fluorescence-activated cell sorting (FACS). GFP-positive cells 
#' were seeded at a 1:3 ratio with untransduced BT549 cells in 
#' 96-well image-lock plates (EssenBio) at a density of 1000 
#' total cells per well. Once cells reached the desired density, 
#' they were scanned at ten-minute intervals over 24h using an 
#' Incucyte S3 Live-Cell microscope (EssenBio) at 10x magnification 
#' and a Basler Ace 1920-155um camera with CMOS sensor. 
#' TIFF images 
#' were imported and processed using the cellmigRation library.
#' 
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
#' 
#'
#' @details
#' BT549 cell trajectories were computed using cellmigRation.  
#' Imaging experiments were performed as described by Ghannoum S 
#' et al (paper in preparation). Briefly, triple negative breast 
#' cancer BT549 cells were cultured in RPMI supplemented with 10% FCS 
#' and 1% penicillin/streptomycin. Cells were transduced with 
#' NucLight green lentivirus (Essen BioScience), and then sorted by  
#' fluorescence-activated cell sorting (FACS). GFP-positive cells 
#' were seeded at a 1:3 ratio with untransduced BT549 cells in 
#' 96-well image-lock plates (EssenBio) at a density of 1000 
#' total cells per well.
#' Cells were scanned at ten-minute intervals over 24h using an 
#' Incucyte S3 Live-Cell microscope (EssenBio) at 10x magnification 
#' and a Basler Ace 1920-155um camera with CMOS sensor. 
#' Cells were treated with 100 uM Rac1 Inhibitor (1177865-17-6, 
#' Calbiochem) or left untreated (controls).  TIFF images 
#' were imported and processed using the cellmigRation library.
#' 
#' @usage data(ThreeConditions)
#' @format a list including 3 elements
#'
#' @examples
#' data(ThreeConditions)
#'
#' @keywords internal
"ThreeConditions"

#' @title Trajectories of 11 cells
#' @description Intermediates and Results from Cell Tracking Analyses,
#' used as a representative example of a S4 CellMig object
#' 
#' @details
#' BT549 cell trajectories were computed by using cellmigRation.  
#' Imaging experiments were performed as described by Ghannoum S 
#' et al (paper in preparation). Briefly, triple negative breast 
#' cancer BT549 cells were cultured in RPMI supplemented with 10% FCS 
#' and 1% penicillin/streptomycin. Cells were transduced with 
#' NucLight green lentivirus (Essen BioScience), and then sorted by  
#' fluorescence-activated cell sorting (FACS). GFP-positive cells 
#' were seeded at a 1:3 ratio with untransduced BT549 cells in 
#' 96-well image-lock plates (EssenBio) at a density of 1000 
#' total cells per well. Once cells reached the desired density, 
#' they were scanned at ten-minute intervals over 24h using an 
#' Incucyte S3 Live-Cell microscope (EssenBio) at 10x magnification 
#' and a Basler Ace 1920-155um camera with CMOS sensor. 
#' TIFF images 
#' were imported and processed using the cellmigRation library.
#' 
#' @usage data(preProcCellMig)
#' @format a list including 21 elements
#'
#' @examples
#' data(preProcCellMig)
#'
#' @keywords internal
"preProcCellMig"

