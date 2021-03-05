#' @title Handle non-NULL ExpName
#' @description The fixExpName helps adjusting the name of the experiment
#' in case it is not NULL.
#' @param x string, name of the experiment.
#'
#' @return A string referring to the adjusted experiment name.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixExpName("Hello World")
#'
#' @keywords internal
fixExpName <- function(x) {
    x <- gsub("[[:space:]]", "_", x)
    xx <- strsplit(x, split = "")[[1]]
    idx <- as.numeric(gregexpr("[[:punct:]]", x)[[1]])
    for(i in idx) {
        if (!xx[i] %in% c("_", "-", ".")) {
            xx[i] <- "."
        }
    }
    y <- paste(xx, collapse = "")
    return(y)
}




#' @title Pre-processing First Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixID1(ID_split = list(data.frame()), TimeInterval = 1)
#'
#'@keywords internal
fixID1 <- function(ID_split,TimeInterval) {

    tryCatch({
        for(j in seq_len(length(ID_split))){
            M<- ID_split[[j]][1]
            MM<-length(M[,1])
            res <- t(vapply(seq_len(MM), function(i){
                # creating values for dx
                ID_split[[j]][i,4]=(ID_split[[j]][i+1,2])-
                    ( ID_split[[j]][i,2])
                ID_split[[j]][,4][is.na(ID_split[[j]][,4])] <- 0
                # creating values for dy
                ID_split[[j]][i,5]= ( ID_split[[j]][i+1,3])-
                    (ID_split[[j]][i,3])
                # to remove NA and replace it with 0
                ID_split[[j]][,5][is.na(ID_split[[j]][,5])] <- 0
                # creating values for dis
                ID_split[[j]][i,6]=sqrt(
                    (ID_split[[j]][i,4])^2 + (ID_split[[j]][i,5])^2)
                ID_split[[j]][i,7]=acos((ID_split[[j]][i,4])/
                                            (ID_split[[j]][i,6]))
                # to remove NA and replace it with 0
                ID_split[[j]][,7][is.na(ID_split[[j]][,7])] <- 0
                # creating values for Square Speed
                ID_split[[j]][i,11]=((ID_split[[j]][i,6])/TimeInterval)^2
                return(as.numeric(ID_split[[j]][i,c(4,5,6,7,11)]))
            }, FUN.VALUE = numeric(5)))
            ID_split[[j]][seq(1,MM,by=1),c(4,5,6,7,11)] <- res
            ID_split[[j]][,c(4,5,6,7,11)] <- lapply(
                ID_split[[j]][,c(4,5,6,7,11)], as.numeric
            )
        }
    }, error = function(e) {
        message("An error may have occurred!");
    })

    return(ID_split)
}


#' @title Pre-processing Second Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'@keywords internal
fixID2<-function(ID_split) {
    for(j in seq_along(ID_split)){
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        res1 <- t(vapply(seq_len(MM), function(i){
            # creating values for cumsum
            ID_split[[j]][,12]=cumsum(ID_split[[j]][,6])
            # creating values for cumulative directionality ratio
            ID_split[[j]][i,13]= sqrt(
                ((ID_split[[j]][i+1,2])^2)+((ID_split[[j]][i+1,3])^2)
            ) / (ID_split[[j]][i,12])
            return(as.numeric(ID_split[[j]][i,c(12,13)]))
        }, FUN.VALUE = numeric(2)))
        ID_split[[j]][seq(1,MM,by=1),c(12,13)] <- res1
        ID_split[[j]][,c(12,13)] <- lapply(ID_split[[j]][,c(12,13)], as.numeric)
    }
    return(ID_split)
}




#' @title Pre-processing Third Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'@keywords internal
fixID3<-function(ID_split) {
    for(j in seq_len(length(ID_split))){
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        MM1<-MM-1
        res <- vapply(seq_len(MM1), function(i){
            if(
                (ID_split[[j]][i+1,5]<0) &&
                (ID_split[[j]][i,5]>=0)||(ID_split[[j]][i+1,5]>=0) &&
                (ID_split[[j]][i,5]<0)
            ){
                ID_split[[j]][i,8]= abs(ID_split[[j]][i+1,7])+
                    abs(ID_split[[j]][i,7])
            }
            if(
                (ID_split[[j]][i+1,5]<0) &&
                (ID_split[[j]][i,5]<0)||(ID_split[[j]][i+1,5]>=0) &&
                (ID_split[[j]][i,5]>=0)
            ){
                ID_split[[j]][i,8]=ID_split[[j]][i+1,7]-ID_split[[j]][i,7]
            }
            ID_split[[j]][i,8]<-ifelse(
                (ID_split[[j]][i,8])<= (-pi),
                2*pi+(ID_split[[j]][i,8]),
                (ID_split[[j]][i,8])
            )        # adjusting the rel.ang
            ID_split[[j]][i,8]<-ifelse(
                (ID_split[[j]][i,8])> pi,
                (ID_split[[j]][i,8])-2*pi,
                (ID_split[[j]][i,8])
            )
            return(ID_split[[j]][i, 8])
        }, FUN.VALUE = numeric(1))
        ID_split[[j]][seq(1,MM1,by=1), 8] <- res
    }
    return(ID_split)
}




#' @title Pre-processing Fourth Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'@keywords internal

fixID4<-function(ID_split) {
    cosine.P<-data.frame()
    for(j in seq_along(ID_split)){
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        res <- vapply(seq_len(MM), function(i){
            ID_split[[j]][i,9]<-cos(ID_split[[j]][i,8])
            return(ID_split[[j]][i,9])
        }, FUN.VALUE = numeric(1))
        ID_split[[j]][seq(1,MM,by=1), 9] <- res
        cosine.P[seq(1,MM,by=1),j]<-ID_split[[j]][,9]
    }
    return(ID_split)
}




#' @title Pre-processing Fifth Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'@keywords internal

fixID5<-function(ID_split,TimeInterval) {
    for(j in seq_along(ID_split)){
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        res <- vapply(seq_len(MM), function(i){
            if(abs(ID_split[[j]][i,8])<=1.5707963268){
                ID_split[[j]][i,10]= TimeInterval
            }
            if(abs(ID_split[[j]][i,8])>1.5707963267){
                ID_split[[j]][i,10]= 0
            }
            return(ID_split[[j]][i,10])
        }, FUN.VALUE = numeric(1))
        ID_split[[j]][seq(1,MM,by=1), 10] <- res
    }
    return(ID_split)
}





#' @title Pre-processing Sixst Part
#' @description This function prepare the data in each data frame
#' as a part of the pre-processing.
#' @param ID_split A list of data frames.
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A list of data frames.
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}

#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'@keywords internal

fixID6<-function(ID_split,TimeInterval) {
    for(j in seq_along(ID_split)){ # Computing Acceleration
        M<- ID_split[[j]][1]
        MM<-length(M[,1])
        MM2<-MM-2
        res <- vapply(seq_len(MM2), function(i){
            ID_split[[j]][i,24]= (
                sqrt(ID_split[[j]][i+1,11])-sqrt(ID_split[[j]][i,11])
            )/TimeInterval
            return(ID_split[[j]][i,24])
        }, FUN.VALUE = numeric(1))
        ID_split[[j]][seq(1,MM2,by=1), 24] <- res
    }
    return(ID_split)
}


#' @title Persistence and Speed First Part
#' @description This function is a part of the PerAndSpeed(), which
#' generates data and plots for persistence and speed.
#' @param x \code{CellMig} class object, which is a list of data
#' @param TimeInterval A numeric value of the time elapsed between
#'
#' @return A data frame named "PerResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixPER1(1, 1)
#'
#'@keywords internal
fixPER1<-function(Object,TimeInterval) {
    # creating a table to store the persistence results
    # creating values (NA and persistence time )for    persistence
    # (step to the original)
    PerResultsTable<-data.frame()
    tryCatch({
        for(j in seq(1,length(Object),by=1)){
            MM<-length(Object[[j]][,1])
            MM2<-MM-1
            Ptime<-Object[[j]][seq(1,MM2,by=1),10]
            MeanPerTime<-round(mean(Ptime),digits=2) # mean persistence time
            PerTimLen<-Ptime
            PerTimLen[PerTimLen==0]<-NA
            PerTimLen<-PerTimLen[!is.na(PerTimLen)]
            PerTimLen<-length(PerTimLen)
            # computing persistence ratio
            PerRatio<-round(PerTimLen/MM2, digits=2)
            # computing the direction deviating time
            DD.Ptime<-Object[[j]][seq(1,MM,by=1),10]
            # replacing the 0 with 5 and the 5 with 0
            DD.Ptime[DD.Ptime==0]<-1
            DD.Ptime[DD.Ptime==TimeInterval]<-0
            DD.Ptime[DD.Ptime==1]<-TimeInterval
            MeanDD.PerTime<-round(mean(DD.Ptime),digits=2)
            PerResultsTable[1,j]<-j
            PerResultsTable[2,j]<-MeanPerTime
            PerResultsTable[3,j]<-MeanDD.PerTime
            PerResultsTable[4,j]<-PerRatio
        }
    }, error = function(e) NULL)
    return(PerResultsTable)
}




#' @title Persistence and Speed Second Part
#' @description This function is a part of the PerAndSpeed(), which
#' generates data and plots for persistence and speed.
#' @param x \code{CellMig} class object, which is a list of data
#' @param TimeInterval A numeric value of the time elapsed between
#' @param ExpName String, name of the experiment
#' @param new.fld path to the folder where to save files
#'
#' @return A data frame named "PerResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom stats cor.test lm
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics title plot abline
#'
#' @examples
#' cellmigRation:::fixPER2(1, 1, 1, 1, 1, 1, 1, 1, 1)
#'
#'@keywords internal
fixPER2<-function(
    Object,PerResultsTable,PtSplot,AllPtSplot,export,
    color,TimeInterval,ExpName,new.fld) {
    VelPerTable<-data.frame()

    tryCatch({
        for(j in seq_along(Object)){
            MM=length(Object[[j]][,1])
            Ptime<-Object[[j]][seq(1,MM,by=1),10]
            Ptime0<-c(0,Ptime) #adding a "0" value in the beginning
            Ptime0[Ptime0==TimeInterval]<-NA #replacing "TimeInterval" with NAs
            Ptime0TF<- !is.na(Ptime0)
            MM1<-MM+1
            rowN<-seq(1,MM1,by=1)
            Tval <- rowN[Ptime0TF]

            fillIdx <- cumsum(Ptime0TF)
            s<-Tval[fillIdx] #replacing the NAs with the previous true value
            justTrueVal<-Tval
            s[Ptime0TF]=0 #replacing the true values with 0
            finalID<-s[-1] # removing the first added value

            Vel<-Object[[j]][seq(1,MM,by=1),11]
            Per<-Object[[j]][seq(1,MM,by=1),10]
            Per[Per==0]<-NA
            tab<-data.frame(Vel,Per,finalID)
            tab<-tab[-nrow(tab),]
            finalIDD<-finalID[-MM] # exclude the last row since it has no veloc.
            tabb<-split(tab, finalIDD) #to avoid taking the last row

            meanVEL<-c()
            Per.time<-c()
            res1 <- vapply(seq_along(tabb), function(i){
                meanVEL= round(sqrt(mean(tabb[[i]][,1])),digits=2)
                return(meanVEL)
            }, FUN.VALUE = numeric(1))
            meanVEL= res1
            meanVEL=meanVEL[-1]
            res2 <- vapply(seq_along(tabb), function(i){
                Per.time=sum(tabb[[i]][,2])
                return(Per.time)
            }, FUN.VALUE = numeric(1))
            Per.time= res2
            Per.time=Per.time[-1]
            w<-which.max(Per.time)
            PerResultsTable[5,j]<-Per.time[w]
            Ptime<-Object[[j]][seq(1,MM,by=1),10] # computing the "meanVel.for0"
            # adding a "1" value in the beginning    this value should not be 0
            # because we will not catch 0 persistence if it is in the beginning.
            Ptime00<-c(1,Ptime)
            Ptime00[Ptime00==0]<-NA # replacing the "0" values with NAs
            Ptime00TF<- !is.na(Ptime00)
            MM1<-MM+1
            rowN<-c(seq(1,MM1,by=1))
            Tval0 <- rowN[Ptime00TF]
            TTT<-c()
            res3 <- vapply(seq_len(length (Tval0)-1), function(i){
                TTT<- (Tval0[i+1]- Tval0[i])-1
                return(TTT)
            }, FUN.VALUE = numeric(1))
            TTT=res3
            TTT[TTT==0]<-NA
            F.ID.for.0.per<-TTT[!is.na(TTT)]
            Final.ID.for.0.per<-c()
            for (i in seq(1,length(F.ID.for.0.per),by=1)){
                Final.ID.for.0.per<-c(
                    Final.ID.for.0.per, rep(i,(F.ID.for.0.per[i])))
            }
            Vel<-Object[[j]][seq(1,MM,by=1),11]
            Per<-Object[[j]][seq(1,MM,by=1),10]
            Per[Per==0]<-NA
            tab<-data.frame()
            tab<-data.frame(Vel,Per,finalID)
            tab<-tab[-nrow(tab),]
            finalIDD<-finalID[-MM]
            tabb<-split(tab, finalIDD)         #to avoid taking the last row
            t<-tabb[[1]]
            tabb1<-cbind(t,Final.ID.for.0.per)
            tabbb<-split(tabb1, Final.ID.for.0.per)
            meanVel.for0<-c()
            res4 <- vapply(seq_along(tabbb), function(i){
                meanVel.for0<-round(sqrt(mean(tabbb[[i]][,1])),digits=2)
                meanVel.for0
            }, FUN.VALUE = numeric(1))
            meanVel.for0=res4
            zerooPrep<-rep(0,length(meanVel.for0))
            Per.time1<-c(zerooPrep,Per.time)
            meanVEL1<-c(meanVel.for0,meanVEL)
            colNumP<-j+j-1
            colNumV<-j+j
            rowNum<-c(seq(1,length(meanVEL1),by=1))
            VelPerTable[rowNum,colNumP]<-Per.time1
            VelPerTable[rowNum,colNumV]<-meanVEL1
            PT<-VelPerTable[,j+j]
            c<-  suppressWarnings(
                stats::cor.test(
                    ~ VelPerTable[,j+j-1]+ PT, method = "spearman",
                    exact=FALSE
                )
            ) #testing the correlation
            cc<-unlist(c[4])
            ccPV<-round(cc, digits = 3)
            PerResultsTable[6,j]<-ccPV # Speed vs persistence time Spearman cor
            if ( PtSplot == TRUE){
                if (export){
                    plot_name <-    paste0(
                        ExpName," Persist Time vs Speed",j,".jpg"
                    )
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path, width=4, height=4, units='in',
                        res = 300
                    )
                }
                graphics::plot(
                    VelPerTable[,j+j]*60,VelPerTable[,j+j-1],
                    pch=16,type="p",ylab="Persistence Time (min)",
                    xlab=" Mean Speed during persistence time    (um/h)",
                    col=color[j],las=1
                )
                reg<-stats::lm(PT~VelPerTable[,j+j-1])
                graphics::title(
                    main=paste0(
                        "Cell Number    ", j,"     Speed vs Persistence Time"
                    ), cex.main =0.7,
                    sub=paste0(
                        "Spearman's rank correlation coefficient = ", ccPV),
                    col.sub="red")
                if (export) grDevices::dev.off()
            }
        }
        ## All cells (Persistence times    vs Speed)
        allper<-VelPerTable[,1]
        allvel<-VelPerTable[,2]
        for(j in seq(1,length(Object),by=1)){
            allper<-c(allper,VelPerTable[,j+j-1])
            allvel<-c(allvel,VelPerTable[,j+j])
        }
        allper<-allper[!is.na(allper)]
        allvel<-allvel[!is.na(allvel)]
        allvel<-(allvel/TimeInterval)*60
        all.per.vel.table<-data.frame(allper,allvel)
        reg<-stats::lm(allper~allvel)
        c<-    suppressWarnings(
            #testing the correlation
            stats::cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)
        )
        cc<-unlist(c[4])
        ccP<-round(cc, digits = 3)
        # Speed vs    persistence time    Spearman correlation
        PerResultsTable[6,(length(Object)+1)]<-ccP
        if ( AllPtSplot == TRUE){
            if (export) {
                plot_name <- paste0(
                    ExpName,"_Persist_Time_vs_Speed-All_Cells.jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4, height = 4, units = 'in',
                    res = 300
                )
            }
            graphics::plot(
                allvel,allper,type="p",pch=16,ylab="Persist_Time (min)",
                xlab=" Mean Speed during persistence time (um/h)",col="black",
                las=1
            )
            graphics::abline(reg,untf=FALSE,col="red")
            graphics::title(
                "Speed vs Persist Time (All cells)",
                cex.main = 0.7,
                sub=paste0("Spearman's rank correlation coefficient = ",ccP),
                col.sub="red"
            )
            if (export) grDevices::dev.off()
        }
    }, error = function(e) NULL)
    return(PerResultsTable)
}




#' @title Persistence and Speed Third Part
#' @description This function is a part of the PerAndSpeed(), which
#' generates data and plots for persistence and speed.
#'
#' @param Object \code{CellMig} class object, which is a list of data.
#' @param PerResultsTable A data frame.
#' @param ApSplot A logical vector that allows generating individual
#' plots of angular persistence vs speed per cell. Default is TRUE.
#' @param AllApSplot A logical vector that allows generating a plot
#' of angular persistence vs speed of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output.
#' @param color A vector of colors that will be used for the plots
#' @param TimeInterval A numeric value of the time elapsed between
#' @param ExpName String, name of the experiment
#' @param new.fld path to the folder where to save files
#'
#' @return A data frame named "PerResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom stats median cor.test lm sd
#' @importFrom grDevices chull jpeg dev.off
#' @importFrom methods slot
#' @importFrom graphics plot abline title
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' cellmigRation:::fixPER3(1,1,1,1,1,1,1,1,1)
#'
#'@keywords internal
fixPER3<-function(
    Object, PerResultsTable, ApSplot, AllApSplot, export,
    color, TimeInterval, ExpName, new.fld) {

    err0 <- tryCatch({

        for(j in seq(1,length(Object),by=1)){
            MM<-length(Object[[j]][,1])
            MM2<-MM-1
            Root.Median.Square.Speed<-round(
                sqrt(stats::median(Object[[j]][seq(1,MM2,by=1),11])),digits = 3
            )*60
            PerResultsTable[7,j]<-Root.Median.Square.Speed
            wma<-which.max(sqrt(Object[[j]][,11]))
            wmi<-which.min(sqrt(Object[[j]][seq(1,MM2,by=1),11]))
            PerResultsTable[8,j]<-round(sqrt(Object[[j]][wma,11]),digits=3)* 60
            PerResultsTable[9,j]<-round(sqrt(Object[[j]][wmi,11]),digits=3)* 60
            mean.cosineP<-round(
                mean(Object[[j]][seq(1,(MM2-1),by=1),9],na.rm = TRUE),
                digits = 3)
            PerResultsTable[10,j]<-mean.cosineP
            s<- suppressWarnings(
                stats::cor.test(
                    ~ sqrt(Object[[j]][seq(1,MM2,by=1),11]) +
                        Object[[j]][seq(1,MM2,by=1),9],
                    method = "spearman",exact=FALSE)) #testing the correlation
            ss<-unlist(s[4])
            VEvsCOSP<-round(ss, digits = 3)
            PerResultsTable[11,j]<-VEvsCOSP
            data<-Object[[j]][seq(1,MM2,by=1), c(2,3)]
            data1<-Object[[j]][seq(1,round(MM2/4),by=1),c(2,3)]
            data2<-Object[[j]][seq(round(MM2/4),round(MM2/2),by=1),c(2,3)]
            data3<-Object[[j]][seq(round(MM2/2),round(MM2*3/4),by=1),c(2,3)]
            data4<-Object[[j]][seq(round(MM2*3/4),MM2,by=1),c(2,3)]
            ch    <- grDevices::chull(data)
            ch1 <- grDevices::chull(data1)
            ch2 <- grDevices::chull(data2)
            ch3 <- grDevices::chull(data3)
            ch4 <- grDevices::chull(data4)
            coords    <- data[c(ch, ch[1]), ]    # closed polygon
            coords1 <- data1[c(ch1, ch1[1]), ]    # closed polygon
            coords2 <- data2[c(ch2, ch2[1]), ]    # closed polygon
            coords3 <- data3[c(ch3, ch3[1]), ]    # closed polygon
            coords4 <- data4[c(ch4, ch4[1]), ]    # closed polygon
            p <- suppressWarnings(sp::Polygon(coords))
            p1 <- suppressWarnings(sp::Polygon(coords1))
            p2 <- suppressWarnings(sp::Polygon(coords2))
            p3 <- suppressWarnings(sp::Polygon(coords3))
            p4 <- suppressWarnings(sp::Polygon(coords4))
            pAr <- methods::slot(p, "area")
            p1Ar <- methods::slot(p1, "area")
            p2Ar <- methods::slot(p2, "area")
            p3Ar <- methods::slot(p3, "area")
            p4Ar <- methods::slot(p4, "area")
            EmptyArea=abs(pAr -(p1Ar + p2Ar + p3Ar + p4Ar))
            SegmentedCA<-p1Ar + p2Ar + p3Ar + p4Ar
            PerResultsTable[12,j]=round(pAr,digits=3)
            PerResultsTable[13,j]=round(SegmentedCA,digits=3)
            PerResultsTable[14,j]=round(EmptyArea,digits=3)
            PerResultsTable[15,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3)
            PerResultsTable[16,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3)
            abs(round(sum(Object[[j]][,8])/6.28,digits=3))
            if ( ApSplot == TRUE){
                if (export) {
                    plot_name <- paste0(
                        ExpName,"_Angular_Persistence_vs_Speed",j,".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path,width = 4, height = 4,
                        units = 'in', res = 300)
                }
                Speed=sqrt(Object[[j]][seq(1,MM2,by=1),11])*60
                graphics::plot(
                    Speed, Object[[j]][seq(1,MM2,by=1),9], pch=16, type="p",
                    ylab="Angular Persistence (cosine)",
                    xlab=" Instantaneous Speed (um/h)", col=color[j], las=1)
                reg <- stats::lm(Object[[j]][seq(1,MM2,by=1),9]~Speed)
                graphics::abline(reg,untf=FALSE,col="black")
                graphics::title(
                    main=paste0(
                        "Cell Number    ", j,
                        " Instantaneous Speeds vs Angular Persistence "
                    ),
                    cex.main = 0.7,
                    sub=paste0(
                        "Spearman's rank correlation coefficient = ",VEvsCOSP
                    ),col.sub="red")
                if (export) grDevices::dev.off()
            }
            PoCos<- subset(
                Object[[j]][seq(1,(MM2-1),by=1),9],
                Object[[j]][seq(1,(MM2-1),by=1),9]>0)
            NeCos<- subset(
                Object[[j]][seq(1,(MM2-1),by=1),9],
                Object[[j]][seq(1,(MM2-1),by=1),9]<=0)
            PerResultsTable[17,j]<-round(
                stats::median(PoCos,na.rm = TRUE),digits = 3
            )
            PerResultsTable[18,j]<-round(
                stats::median(NeCos,na.rm = TRUE),digits = 3
            )
            tmp.Rng <- seq(1,length(Object[[j]][,1]),by=1)-1
            AvSp<-(Object[[j]][tmp.Rng,6]/TimeInterval)*60
            PerResultsTable[19,j]<-round(
                stats::median(AvSp,na.rm = TRUE),digits = 3
            )
            s=summary(AvSp)
            PerResultsTable[20,j]<-round((s[5]-s[2])/s[3],digits=4)
            PerResultsTable[21,j]<-round(mean(AvSp),digits = 3)
            PerResultsTable[22,j]<-round(stats::sd(AvSp),digits = 3)
        }
        #(all cells)    persistence vs inst.speed
        MM<-length(Object[[1]][,1])
        MM2<-MM-1
        cosine.P<-data.frame()
        # creating values for    cosine.P    based on rel.ang.P
        for(j in seq_along(Object)){
            M<- Object[[j]][1]
            MM<-length(M[,1])
            res <- vapply(seq_len(MM), function(i){
                cos(Object[[j]][i,8])
            }, FUN.VALUE = numeric(1))
            Object[[j]][seq(1,MM,by=1), 9] <- res
            cosine.P[seq(1,MM,by=1),j]<-Object[[j]][,9]
        }
        RM<-round(
            matrixStats::rowMedians(as.matrix(cosine.P[seq(1,MM2,by=1),]),
                                    na.rm = TRUE), digits=3)
        Speed<-data.frame()
        # calculating the Mean.Square.velocity for each cell
        for (j in seq(1,length(Object),by=1)){
            Speed[seq(1,MM2,by=1), j] <-
                round(sqrt(Object[[j]][seq(1,MM2,by=1),11]), digits = 3)
        }
        RowmeanSpeed<-round(
            matrixStats::rowMedians(as.matrix(Speed),na.rm = TRUE), digits=3)
        #testing the correlation
        s<-stats::cor.test(~RM + RowmeanSpeed, method = "spearman",exact=FALSE)
        ss<-unlist(s[4])
        VEvsCOSP<-round(ss, digits = 3)
        PerResultsTable[11,(length(Object)+1)]<-VEvsCOSP
        nuFlName <- "All_Cells_Average_Angular_Persistence_vs_Average_Speed.jpg"
        if ( AllApSplot == TRUE){
            if (export) {
                plot_name <- paste(ExpName, nuFlName)
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4, height = 4, units = 'in',
                    res = 300
                )
            }
            MS<-max(RowmeanSpeed)*60
            graphics::plot(
                RowmeanSpeed*60,RM,pch=16,type="p",
                ylab="Average Angular Persistence (cosine)",
                xlab=" Average Instantaneous Speed (um/h)",col="black",las=1,
                xlim=c(0,MS)
            )
            NewSpeed=RowmeanSpeed*60
            reg<-stats::lm(RM~NewSpeed)
            graphics::abline(reg,untf=FALSE,col="red")
            mySubLab <- "Spearman's rank correlation coefficient = "
            graphics::title(
                main="All Cells Instantaneous Speed vs Angular Persistence",
                sub=paste0(mySubLab, VEvsCOSP), cex.main = 0.7, col.sub="red")
            if (export) try(grDevices::dev.off(), silent = TRUE)
        }
        # return 0 to err0 if no errors occurred... otherwise return 1.
        0
    }, error = function(e) {
        1
    })

    if (err0 > 0) {
        message("The fixPER3() function is an internal function...")
        message("Are you running it as a standalone function?")
        PerResultsTable <- NULL
    }
    return(PerResultsTable)
}


#' @title Mean Square Displacement
#' @description This function is a part of the MSD function, which
#' computes the mean square displacements across several sequential
#' time intervals.
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param Step A numeric value of the number of trajectory steps.
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
#' @param color A vector of colors that will be used for the plots
#' @param ExpName String, name of the experiment
#' @param new.fld path to the folder where to save files
#'
#' @return A data frame named "MSDResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#'
#' @importFrom stats lm coef
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics par plot title abline lines
#' @importFrom matrixStats rowMedians rowSds
#' @importFrom FME modCost modFit
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' cellmigRation:::fixMSD(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#'
#'@keywords internal
fixMSD<-function(
    Object, Step, SlopePlot, AllSlopesPlot, FurthPlot,
    AllFurthPlot, sLAG, ffLAG, color, export, ExpName, new.fld) {

    MSDResultsTable <- data.frame()

    tryCatch({
        # table that has all MSDs to compute the mean and sd
        MSD.table<-data.frame()
        for(j in seq_along(Object)){
            meanSD<-c()
            LAG<-round(Step*sLAG) # number of lags is based on sLAG
            for(lag in seq(1,LAG,by=1)){
                res <- vapply(seq_len(Step), function(i){
                    ((Object[[j]][i+lag,2]- Object[[j]][i,2])^2) +
                        ((Object[[j]][i+lag,3]- Object[[j]][i,3])^2)
                }, FUN.VALUE = numeric(1))
                Object[[j]][seq(1,Step,by=1), 22] <- res
                Object[[j]][,22][is.na(Object[[j]][,22])] <- 0
                meanSD[lag]<-mean(Object[[j]][seq(1,(Step-LAG),by=1),22])
                Object[[j]][,22]=0
            }
            MSDResultsTable[1,j]<-j
            MSDResultsTable[2,j]<-round(meanSD[1],digits=3)
            MSD.table[seq(1,LAG,by=1),j]<-meanSD
            Yaxis<-assign(paste0("MSD.Cell.",j),meanSD)
            Xaxis<-seq(1,LAG,by=1)
            # plotting the regression line based on sLAG
            NewrowMeans<-meanSD[seq(1,(round(LAG*sLAG)),by=1)]
            NewXaxis<-Xaxis[seq(1,(round(LAG*sLAG)),by=1)]
            reg<-stats::lm(NewrowMeans~ NewXaxis)
            reg1<-round(
                stats::coef(stats::lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],
                digits=2
            )
            MSDResultsTable[3,j]<-reg1
            if (SlopePlot) {
                if (export) {
                    plot_name <- paste0(ExpName, "-MSD.plot.Cell", j, ".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path,
                        width = 4, height = 4, units = 'in', res = 300
                    )
                }
                xn <- expression(paste("MSD (um"^2, ")"))
                graphics::par(mar=c(5.1,5.9,4.1,1.1), mgp=c(3.5,0.5,0), las=0)
                graphics::plot(
                    Xaxis,Yaxis, type="p",col=color[j],xlab="Lag",ylab=xn,
                    pch=19, las=1,log="xy",cex=1.2)
                graphics::title(
                    main=paste0("Cell Number    ", j," -    MSD Slope = ",reg1),
                    col.main="black",cex.main=0.8)
                try({graphics::abline(reg,untf=TRUE,col="black")}, silent=TRUE)
                if (export) grDevices::dev.off()
            }
        }
        if (AllSlopesPlot) {
            RM1<-matrixStats::rowMedians(as.matrix(MSD.table),na.rm = TRUE)
            RSD<-matrixStats::rowSds(as.matrix(MSD.table),na.rm = TRUE)
            xn <- expression(paste("MSD (um"^2, ")"))
            LAG<-round(Step*sLAG)
            Xaxis<-seq(1,LAG,by=1)
            if (export) {
                plot_name <- paste0(ExpName, "-MSD.plot All Cells.jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,
                    width = 4, height = 4, units = 'in', res = 300)
            }
            xn <- expression(paste("MSD (um"^2, ")"))
            graphics::par(mar=c(5.1, 5.5, 4.1, 0.9), mgp=c(3.5, .5, 0), las=0)
            graphics::plot(
                Xaxis,RM1, type="p",col="black",xlab="Lag",ylab=xn,pch=19,
                las=1,cex=1.2,log="xy"
            )
            # best fit based on sLAG *sLAG
            NewrowMeans<-RM1[seq(1,(round(LAG*sLAG)+1),by=1)]
            NewXaxis<-Xaxis[seq(1,(round(LAG*sLAG)+1),by=1)]
            reg<-stats::lm(NewrowMeans~ NewXaxis)
            reg1<-round(
                stats::coef(
                    stats::lm(log10(NewrowMeans)~ log10(NewXaxis))
                )[2],digits=2
            )
            graphics::abline(reg,untf=TRUE,col="red")
            graphics::title(
                main=paste0("All Cells -    MSD Slope = ",reg1),
                col.main="black", cex.main=0.8)
            if (export) grDevices::dev.off()
        }
        # Fitting the Furth formula using generalized regression by
        # the Nelder–Mead method simplex method
        for (j in seq(1,length(MSD.table[1,]),by=1)){
            LAG<-round(Step*ffLAG)
            y=MSD.table[seq(1,LAG,by=1),j]
            t <- seq(1,length(y))
            Data<-c(t,y)
            Data <- matrix(ncol =    2, byrow = FALSE, data =Data)
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
            Fit<-FME::modFit(
                p = c(D = 1, P = 1), lower = c(0,0),upper=c(100,100),
                f = Cost,method = "Nelder-Mead"
            )
            Summary<-summary(Fit)
            MSDResultsTable[4,j]<-round(Fit$par[1],digits=3)
            MSDResultsTable[5,j]<-round(Fit$par[2],digits=3)
            MSDResultsTable[6,j]<-round(Summary$par[1,4],digits=3)
            MSDResultsTable[7,j]<-round(Summary$par[2,4],digits=3)
            if (FurthPlot){
                if (export) {
                    plot_name <- paste0(
                        ExpName, "-MSD.N-M.bestfit.Cell", j, ".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path, width=4, height = 4,
                        units = 'in', res = 300)
                }
                graphics::plot(
                    Data, pch = 16,col=color[j], cex = 1.2,
                    xlab = "Lags", ylab = "MSD")
                x<-seq(0,LAG,1)
                Model <- function(p, x) {
                    return(
                        data.frame(x=x, y=p[1]*4*(x- p[2]*(1-(exp(-x/p[2])))))
                    )
                }
                graphics::lines(Model(Fit$par, x),col="black")
                graphics::title(
                    main=paste(
                        "Cell Number   ", j, "       Nelder-Mead best-fit"),
                    col.main="black", cex.main=0.8,
                    sub=paste0(
                        "D = ",round(Fit$par[1],digits=3),
                        "                    P = ", round(Fit$par[2],digits=3)),
                    col.sub="red")
                if (export) grDevices::dev.off()
            }
        }
        nuCI <- (ncol(MSDResultsTable)+1)
        RM1<-matrixStats::rowMedians(as.matrix(MSD.table),na.rm = TRUE)
        MSDResultsTable[1,nuCI]<-"All Cells"
        MSDResultsTable[2,nuCI]<-round(RM1[1],digits=3)
        RM<-c(0,RM1)
        Xaxis<-seq(1,round(Step*ffLAG),by=1)
        #best fit based on ffLAG *ffLAG
        NewrowMeans<-RM1[seq(1,(round(LAG*ffLAG)+1),by=1)]
        NewXaxis<-Xaxis[seq(1,(round(LAG*ffLAG)+1),by=1)]
        reg<-stats::lm(NewrowMeans~ NewXaxis)
        reg1<-round(stats::coef(stats::lm(
            log10(NewrowMeans) ~ log10(NewXaxis)))[2],digits=2)
        MSDResultsTable[3,nuCI]<-reg1
        y=RM1[seq(1,round(Step*ffLAG),by=1)]
        t <- c(seq(1,length(y),by=1))
        Data<-c(t,y)
        Data <- matrix(ncol =    2, byrow = FALSE, data =Data)
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
        Fit<-FME::modFit(
            p = c(D = 1, P = 1), lower = c(0, 0),upper=c(100,100),f = Cost,
            method = "Nelder-Mead"
        )
        MSDResultsTable[4,nuCI]<-round(Fit$par[1],digits=3)
        MSDResultsTable[5,nuCI]<-round(Fit$par[2],digits=3)
        MSDResultsTable[6,nuCI]<-round(Summary$par[1,4],digits=3)
        MSDResultsTable[7,nuCI]<-round(Summary$par[2,4],digits=3)
        if (AllFurthPlot) {
            if (export) {
                plot_name <- paste0(ExpName, "-MSD.N-M.bestfit.All.Cells.jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,
                    width = 4, height = 4, units = 'in', res = 300
                )
            }
            graphics::plot(
                Data, pch = 16,col="black", cex = 1.2,
                xlab = "Lags", ylab = "MSD")
            x<-seq(0,LAG,1)
            Model <- function(p, x) {
                return(data.frame(x=x, y=p[1]*4*(x- p[2]*(1-(exp(-x/p[2]))))))
            }
            graphics::lines(Model(Fit$par, x),col="red")
            graphics::title(
                main=paste0("All Cells Nelder-Mead best-fit"),col.main="black",
                cex.main=0.8,
                sub=paste0(
                    "D = ",round(Fit$par[1],digits=3),
                    "                    P = ", round(Fit$par[2],digits=3)
                ),
                col.sub="red"
            )
            if (export) grDevices::dev.off()
        }
    }, error = function(e) {NULL})

    return(MSDResultsTable)
}



#' @title Direction AutoCorrelation
#' @description This function is a part of the DiAutoCor function, which
#' computes the angular persistence across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list
#' of data frames resulted from the PreProcessing.
#' @param Step A numeric value of the number of trajectory steps.
#' @param sLAG A numeric value to be used to get the number of lags
#' for the slope fitting. Default is 0.25, which represents 25
#' percent of the steps.
#' @param sPLOT A logical vector that allows generating individual
#' plots showing the angular persistence across several sequantial
#' time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot
#' showing the angular persistence across several sequantial time
#' intervals of all cells. Default is TRUE.
#' @param export if `TRUE` (default), exports function output
#' @param color A vector of colors that will be used for the plots
#' @param ExpName String, name of the experiment
#' @param new.fld path to the folder where to save files
#'
#' @return A data frame named "DA.ResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom stats lm predict median
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics plot lines title abline legend
#' @importFrom matrixStats rowMedians
#'
#' @examples
#' cellmigRation:::fixDA(1, 1, 1, 1, 1, 1, 1, 1, 1)
#'
#'@keywords internal
fixDA<-function(Object,Step,sLAG,sPLOT,aPLOT,color,export,ExpName,new.fld) {
    DA.ResultsTable<-data.frame()
    DA.table<-data.frame()

    tryCatch({
        for(j in seq_along(Object)){
            cos.diff<-c()
            LAG<-round(Step*sLAG) #taking only the first 12.5%
            for(lag in seq_len(LAG)){
                # starting from 2 to exclude the first cosine which is always 1.
                res <- t(vapply(seq(1,Step,by=1), function(i){
                    Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2]
                    Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3]
                    return(as.numeric(Object[[j]][i,c(14,15)]))
                }, FUN.VALUE = numeric(2)))
                Object[[j]][seq(1,Step,by=1),c(14,15)] <- res
                Object[[j]][,c(14,15)] <- lapply(
                    Object[[j]][,c(14,15)], as.numeric)
                Object[[j]][,14][is.na(Object[[j]][,14])] <- 0 # replace NAs
                Object[[j]][,15][is.na(Object[[j]][,15])] <- 0 # idem
                res1 <- vapply(seq_len(Step), function(i){
                    acos(
                        (Object[[j]][i+lag,2]-Object[[j]][i,2]) /
                            sqrt(
                                (Object[[j]][i+lag,2]-Object[[j]][i,2])^2 +
                                    (Object[[j]][i+lag,3]-Object[[j]][i,3])^2
                            )
                    ) # to find the abs angle
                }, FUN.VALUE = numeric(1))
                Object[[j]][seq(1,(Step),by=1),16] <- res1
                # to remove NA and replace it with 0
                Object[[j]][,16][is.na(Object[[j]][,16])] <- 0
                res2 <- vapply(seq_len(Step - lag), function(i){
                    if(
                        (Object[[j]][i+1,15]<0) &&
                        (Object[[j]][i,15]>=0)||(Object[[j]][i+1,15]>=0) &&
                        (Object[[j]][i,15]<0)
                    ){
                        abs(Object[[j]][i+1,16])+abs(Object[[j]][i,16])
                    } else if(
                        (Object[[j]][i+1,15]<0) &&
                        (Object[[j]][i,15]<0)||(Object[[j]][i+1,15]>=0) &&
                        (Object[[j]][i,15]>=0)
                    ){
                        Object[[j]][i+1,16]-Object[[j]][i,16]
                    }
                }, FUN.VALUE = numeric(1))
                Object[[j]][seq(1,(Step - lag),by=1),17] <- res2
                res3 <- t(vapply(seq_len(Step - lag), function(i){
                    Object[[j]][i,17]<-ifelse(
                        (Object[[j]][i,17])<= (-pi),
                        2*pi+(Object[[j]][i,17]),
                        (Object[[j]][i,17])
                    )        # adjusting the ang.diff
                    Object[[j]][i,17]<-ifelse(
                        (Object[[j]][i,17])>= pi,
                        (Object[[j]][i,17])-2*pi,
                        (Object[[j]][i,17])
                    )
                    Object[[j]][i,18]<-cos(Object[[j]][i,17])
                    return(as.numeric(Object[[j]][i,c(17,18)]))
                }, FUN.VALUE = numeric(2)))
                Object[[j]][seq_len(Step - lag),c(17,18)] <- res3
                Object[[j]][,c(17,18)]<-lapply(
                    Object[[j]][,c(17,18)], as.numeric)
                for(i in seq(1,LAG,by=1)){  # computing the cosine mean
                    tmpSQ <- seq(1,(Step-lag-1),by=1)
                    cos.diff[lag]<-mean(Object[[j]][tmpSQ,18], na.rm=TRUE)
                }
            }
            DA.table[seq(1,LAG,by=1) ,j]<-cos.diff
            assign(paste0("DA.Cell.",j),cos.diff)
            DA.ResultsTable[1,j]<-j
            DA.ResultsTable[2,j]<-round(cos.diff[1],digits=3)
            lags<-c(seq(1,length(DA.table[,1]),by=1))
            lags2<- lags^2
            quadratic.model<-c()
            quadratic.m <-stats::lm(DA.table[,j]~ lags + lags2)
            c<-quadratic.m
            cc<-unlist(c)
            quadratic.model[j]<-cc[1]
            DA.ResultsTable[3,j]<-round(unlist(cc[1]),digits=3)
            ccc<-unlist(cc[1])
            DA.ResultsTable[4,j]<-round(mean(DA.table[seq(1,LAG,by=1),j]),
                                        digits=3)
            DA.ResultsTable[5,j]<-round(
                (1-(mean(DA.table[seq(1,LAG,by=1),j])-unlist(cc[1]))
                ) * mean(DA.table[seq(1,LAG,by=1),j]), digits=3
            )
            DA.ResultsTable[6,j]<-round(
                mean(DA.table[seq(1,LAG,by=1),j])-unlist(cc[1]),digits=3
            )
            timevalues <- seq(1, length(lags), 1)
            predictedcounts <- stats::predict(
                quadratic.m,list(Time=timevalues, Time2=timevalues^2)
            )
            if ( sPLOT == TRUE){
                Xaxis<-seq(1,LAG,by=1)
                Yaxis<-cos.diff
                if (export) {
                    plot_name <-    paste0(
                        ExpName,"Direction Autocorrelation.plot.Cell",j,".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path,width = 4, height = 4,
                        units = 'in', res = 300)
                }
                graphics::plot(
                    Xaxis,Yaxis, type="o",ylim=c(-1,1),xlim=c(0,lag),
                    col=color[j], xlab="Lag",ylab="Cosine",pch=19,las=1,cex=1.2)
                xx<-c(0,1)
                yy<-c(1,cos.diff[1])
                graphics::lines(xx,yy, type='l',col=color[j])
                graphics::lines(
                    timevalues, predictedcounts, col = "black", lwd = 3
                )
                graphics::title(
                    main=paste0("Cell Number    ",j, "     DA quadratic model"),
                    col.main="darkgreen",cex.main=0.8,
                    sub=paste0(
                        " Intercept of DA quadratic model = ",
                        round(ccc, digits=3)),col.sub="red")
                if (export) grDevices::dev.off()
            }
            Object[[j]][,seq(15,18,by=1)]=0
        }
        RM1<-matrixStats::rowMedians(as.matrix(DA.table),na.rm = TRUE)
        DA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
        DA.ResultsTable[2,(length(Object)+1)]<-round(RM1[1],digits=3)
        lags<-seq(1,length(DA.table[,1]),by=1)
        lags2<- lags^2
        quadratic.model<-c()
        quadratic.m <-stats::lm(RM1~ lags + lags2)
        c<-quadratic.m
        cc<-unlist(c)
        quadratic.model[j]<-cc[1]
        DA.ResultsTable[3,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
        DA.ResultsTable[4,(length(Object)+1)]<-round(
            stats::median(as.numeric(
                DA.ResultsTable[4,seq(1,length(Object),by=1)])),digits=3)
        DA.ResultsTable[5,(length(Object)+1)]<-round(
            (1 - (stats::median(as.numeric(
                DA.ResultsTable[4,seq(1,length(Object),by=1)])) -
                    unlist(cc[1]))
            ) * stats::median(as.numeric(
                DA.ResultsTable[4,seq(1,length(Object),by=1)])),
            digits=3)
        DA.ResultsTable[6,(length(Object)+1)]<-round(
            median(as.numeric(DA.ResultsTable[4,seq(1,length(Object),by=1)])) -
                unlist(cc[1]),digits=3)
        ccc<-unlist(cc[1])
        timevalues <- seq(1, length(lags), 1)
        predictedcounts <- stats::predict(
            quadratic.m,list(Time=timevalues, Time2=timevalues^2)
        )
        if ( aPLOT == TRUE){
            Xaxis<-seq(1,LAG,by=1)
            Yaxis<-RM1
            if (export) {
                plot_name <- paste0(
                    ExpName,"-Direction Autocorrelation All Cells.jpg"
                )
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path,width = 4,
                    height = 4, units = 'in',res = 300)
            }
            graphics::plot(
                Xaxis, Yaxis, type="o",ylim=c(-1,1), xlim=c(0,lag), col="black",
                xlab="Lag", ylab="Cosine", pch=19, las=1, cex=1.2)
            xx<-c(0,1)
            yy<-c(1,RM1[1])
            graphics::lines(xx,yy, type='l',col="black")
            graphics::lines(
                timevalues, predictedcounts, col = "darkgreen", lwd = 3
            )
            MDA<-round(
                stats::median(as.numeric(
                    DA.ResultsTable[4,seq(1,length(Object),by=1)])),
                digits=2)
            graphics::abline(h=MDA,col="blue",lwd = 2)
            graphics::title(
                main=paste0("All Cells - DA quadratic model"),
                col.main="darkgreen",cex.main=0.8,
                sub=paste(
                    " Intercept of DA quadratic model =",round(ccc, digits=3)
                ),col.sub="red"
            )
            graphics::legend(
                1, y=-0.82, legend=c(
                    "Mean Direction AutoCorrelation", "Quadratic model"),
                col=c("blue","darkgreen"),lty=1, cex=0.8)
            if (export) grDevices::dev.off()
        }
    }, error = function(e) {NULL})
    return(DA.ResultsTable)
}







#' @title Forward Migration First Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param Step A numeric value of the number of trajectory steps.
#' @return An CellMig class Object.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixFM1(1, 1)
#'
#' @keywords internal
fixFM1<-function(Object,Step) {

    # creating values for    rel.ang.F    (step to the original)
    tryCatch({
        for(j in seq_along(Object)){
            MM    <- Step
            MM1 <- MM - 1
            res <- vapply(seq_len(MM1), function(i) {
                if(
                    (Object[[j]][1,25]==0) &&
                    (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) &&
                    (Object[[j]][i,5]<0)
                ){
                    Object[[j]][i,19]= 1.5707963268 - abs(Object[[j]][i,7])
                }
                if(
                    (Object[[j]][1,25]==0) &&
                    (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) &&
                    (Object[[j]][i,5]>0)
                ){
                    Object[[j]][i,19]= abs(Object[[j]][i,7])+1.5707963268
                }
                Object[[j]][i,19]<-ifelse(
                    (Object[[j]][i,19])<= (-pi),
                    2*pi+(Object[[j]][i,19]),
                    (Object[[j]][i,19])
                )        # adjusting the rel.ang
                Object[[j]][i,19]<-ifelse(
                    (Object[[j]][i,19])>= pi,
                    (Object[[j]][i,19])-2*pi,
                    (Object[[j]][i,19])
                )
                return(as.numeric(Object[[j]][i, 19]))
            }, FUN.VALUE = numeric(1))
            Object[[j]][seq(1,MM1,by=1) , 19] <- res
        }
    }, error = function(e) NULL)
    return(Object)
}



#' @title Forward Migration Second Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param Step A numeric value of the number of trajectory steps.
#' @return An CellMig class Object.

#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixFM2(1, 1)
#'
#' @keywords internal
fixFM2<-function(Object,Step) {

    tryCatch({
        for(j in seq_along(Object)){
            MM<-Step
            res <- vapply(seq_len(MM), function(i){
                if(
                    (Object[[j]][1,25]==0) &&
                    (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) &&
                    (Object[[j]][i,5]<0)
                ){     # upper cell going up or lower cell going down
                    Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
                }
                if(
                    (Object[[j]][1,25]==0) &&
                    (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) &&
                    (Object[[j]][i,5]>0)
                ){
                    Object[[j]][i,20]<-cos(Object[[j]][i,19])
                }
                return(as.numeric(Object[[j]][i,20]))
            }, FUN.VALUE = numeric(1))
            Object[[j]][seq(1,MM,by=1) , 20] <- as.data.frame(res)
        }
    }, error = function(e) NULL)
    return(Object)
}




#' @title Forward Migration Third Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param Step A numeric value of the number of trajectory steps.
#' @return A data frame named"cosine.FP".

#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixFM3(1, 1)
#'
#' @keywords internal
fixFM3<-function(Object,Step) {
    cosine.FP<-data.frame()
    tryCatch({
        for(j in seq_along(Object)){
            MM<-Step
            cosine.FP[seq(1,MM,by=1),j]<-Object[[j]][,20]
        }
    }, error = function(e) NULL)
    return(cosine.FP)
}


#' @title Forward Migration Fourth Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param Step A numeric value of the number of trajectory steps.
#' @return An CellMig class Object.
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixFM4(1, 1, 1)
#'
#' @keywords internal
fixFM4<-function(Object,TimeInterval,Step) {

    tryCatch({
        ## Forward Pesrsistence time, FP deviating time and FP ratio #
        for(j in seq(1,length(Object),by=1)){
            MM<-Step
            MM1<-MM-1
            res <- vapply(seq_len(MM1), function(i){
                if(abs(Object[[j]][i,19])<1.5707963268){
                    Object[[j]][i,21]= TimeInterval
                }
                if(abs(Object[[j]][i,19])>=1.5707963267){
                    Object[[j]][i,21]= 0
                }
                return(as.numeric(Object[[j]][i,21]))
            }, FUN.VALUE = numeric(1))
            Object[[j]][seq(1,MM1,by=1),21] <- res
            Object[[j]][MM1,21] <- TimeInterval
        }
    }, error = function(e) NULL)
    return(Object)
}



#' @title Forward Migration Fifth Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param Step A numeric value of the number of trajectory steps.
#' @return A data frame named "FMResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @examples
#' cellmigRation:::fixFM5(1, 1, 1)
#'
#' @keywords internal
fixFM5<-function(Object,TimeInterval,Step) {
    FMResultsTable<-data.frame()

    tryCatch({
        # creating values (NA and forward persistence time )for
        # forward persistence        (step to the forward movement)
        for(j in seq(1,length(Object), by=1)){
            MM<-Step
            F.P.time<-Object[[j]][seq(1,MM,by=1),21]
            # computing mean    Forward Pesrsistence time
            MeanF.P.time<-round(mean(F.P.time),digits=2)
            # comp. the number of FP steps to be used in computing the FP ratio
            F.P.time.Len<-F.P.time
            F.P.time.Len[F.P.time.Len==0]<-NA
            F.P.time.Len<-F.P.time.Len[!is.na(F.P.time.Len)]
            F.P.time.Len<-length(F.P.time.Len)
            # computing FP ratio
            F.P.Ratio<- round(F.P.time.Len/MM, digits=2)
            # computing the FP deviating time
            DD.F.P.time<-Object[[j]][seq(1,MM,by=1),21]
            DD.F.P.time[DD.F.P.time==0]<-1
            DD.F.P.time[DD.F.P.time==TimeInterval]<-0
            DD.F.P.time[DD.F.P.time==1]<-TimeInterval
            MeanDD.F.P.time<-round(mean(DD.F.P.time),digits=2)
            FMResultsTable[1,j]<-j
            FMResultsTable[2,j]<-MeanF.P.time
            FMResultsTable[3,j]<-MeanDD.F.P.time
            FMResultsTable[4,j]<-F.P.Ratio
        }
    }, error = function(e) NULL)
    return(FMResultsTable)
}



#' @title Forward Migration Sixst Part
#'
#' @description This function is a part of the ForwardMigration function,
#' which generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data
#' frames resulted from the PreProcessing.
#' @param FMResultsTable A data frame resulted from the fixFM6().
#' @param TimeInterval A numeric value of the time elapsed between
#' successive frames in the time-lapse stack.
#' @param Step A numeric value of the number of trajectory steps.

#' @param sfptPLOT A logical vector that allows generating individual
#' plots of persistence time vs speed per cell. Default is TRUE.
#' @param afptPLOT  A logical vector that allows generating a plot of
#' persistence time vs speed for all cells. Default is TRUE.
#' @param sfpPLOT A logical vector that allows generating individual
#' plots of angular persistence vs speed per cell. Default is TRUE.
#' @param export if `TRUE` (default), exports function output to CSV
#' file
#' @param color A vector of colors that will be used for the plots
#' @param ExpName String, name of the experiment
#' @param new.fld path to the folder where to save files
#'
#'
#' @return A data frame named "FMResultsTable".
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#' @references
#' \url{https://www.data-pulse.com/dev_site/cellmigration/}
#'
#' @importFrom stats cor.test lm
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics plot title abline
#'
#' @examples
#' cellmigRation:::fixFM6(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#'
#' @keywords internal
fixFM6<-function(
    Object, FMResultsTable, Step, sfptPLOT, afptPLOT, export,
    color, TimeInterval, sfpPLOT, ExpName, new.fld) {

    outY <- NULL
    VelFPTable<-data.frame()
    tryCatch({
        for(j in seq(1,length(Object),by=1)){
            MM=Step #computing the "finalID" to be used for the splitting
            FPtime<-Object[[j]][seq(1,MM,by=1),21]
            FPtime0<-c(0,FPtime) #adding a "0" value in the beginning
            FPtime0[FPtime0==TimeInterval]<-NA #replacing the "5" with NAs
            FPtime0TF<- !is.na(FPtime0)
            MM1<-Step+1
            rowN<-c(seq(1,MM1,by=1))
            Tval <- rowN[FPtime0TF]
            fillIdx <- cumsum(FPtime0TF)
            s<-Tval[fillIdx]   #replacing the NAs with the previous true value
            justTrueVal<-Tval
            s[FPtime0TF]=0    #replacing the true values with 0
            finalID<-s[-1]    # removing the first added value
            Vel<-Object[[j]][seq(1,MM,by=1),11]
            FP<-Object[[j]][seq(1,MM,by=1),21]
            FP[FP==0]<-NA
            tab<-data.frame(Vel,FP,finalID)
            tab<-tab[-nrow(tab),]
            finalIDD<-finalID[-MM] # exclude the last row since it has no veloc.
            tabb<-split(tab, finalIDD) # to avoid taking the last row
            meanVEL<-c()
            FP.time<-c()
            meanVEL <- vapply(seq_along(tabb), function(i){
                round(sqrt(mean(tabb[[i]][,1])),digits=2)
            }, FUN.VALUE = numeric(1))
            meanVEL=meanVEL[-1]
            FP.time <- vapply(seq_along(tabb), function(i){
                sum(tabb[[i]][,2])
            }, FUN.VALUE = numeric(1))
            FP.time=FP.time[-1]
            w<-which.max(FP.time)
            FMResultsTable[5,j]<-FP.time[w]
            FPtime<-Object[[j]][seq(1,MM,by=1) ,21] # computing "meanVel.for0"
            #adding a "1" value in the beginning this value should not be 0
            # because we will not catch 0 persistence if it is in the beginning
            FPtime00<-c(1,FPtime)
            FPtime00[FPtime00==0]<-NA #replacing the "0" values with NAs
            FPtime00TF<- !is.na(FPtime00)
            MM1<-MM+1
            rowN<-seq(1,MM1,by=1)
            Tval0 <- rowN[FPtime00TF]
            TTT <- vapply(seq_len(length (Tval0)-1), function(i){
                TTT<- (Tval0[i+1]- Tval0[i])-1
                return(TTT)
            }, FUN.VALUE = numeric(1))
            TTT[TTT==0]<-NA
            F.ID.for.0.FP<-TTT[!is.na(TTT)]
            Final.ID.for.0.FP<-c()
            for (i in seq(1,length(F.ID.for.0.FP),by=1)){
                Final.ID.for.0.FP <-
                    c(Final.ID.for.0.FP,rep(i,(F.ID.for.0.FP[i])))
            }
            Vel<-Object[[j]][seq(1,MM,by=1) ,11]
            FP<-Object[[j]][seq(1,MM,by=1),21]
            FP[FP==0]<-NA
            tab<-data.frame()
            tab<-data.frame(Vel,FP,finalID)
            tab<-tab[-nrow(tab),]
            finalIDD<-finalID[-MM]
            tabb<-split(tab, finalIDD)         #to avoid taking the last row
            t<-tabb[[1]]
            tabb1<-cbind(t,Final.ID.for.0.FP)
            tabbb<-split(tabb1, Final.ID.for.0.FP)
            meanVel.for0<-c()
            meanVel.for0 <- vapply(seq_along(tabbb), function(i){
                round(sqrt(mean(tabbb[[i]][,1])),digits=2)
            }, FUN.VALUE = numeric(1))
            zerooFP<-rep(0,length(meanVel.for0))
            FP.time1<-c(zerooFP,FP.time)
            meanVEL1<-c(meanVel.for0,meanVEL)
            colNumP<-j+j-1
            colNumV<-j+j
            rowNum<-seq(1,length(meanVEL1),by=1)
            VelFPTable[rowNum,colNumP]<-FP.time1
            VelFPTable[rowNum,colNumV]<-meanVEL1
            c<-stats::cor.test(
                ~ VelFPTable[,j+j-1]+ VelFPTable[,j+j], method = "spearman",
                exact=FALSE
            )   #testing the correlation
            cc<-unlist(c[4])
            ccPV<-round(cc, digits = 3)
            FMResultsTable[6,j]<-ccPV # Speed vs persistence time Spearman cor
            if ( sfptPLOT == TRUE){
                if (export){
                    plot_name <- paste0(ExpName,"_FP_Time _VS_Speed_",j,".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path,width = 4,
                        height = 4, units = 'in',res = 300)
                }
                graphics::plot(
                    VelFPTable[,j+j]*60,VelFPTable[,j+j-1],pch=16,type="p",
                    ylab="Forward Persistence Time (min)",
                    xlab=" Mean Speed during FP time (um/h)",col=color[j],las=1
                )
                reg<-stats::lm(VelFPTable[,j+j]~VelFPTable[,j+j-1])
                graphics::title(
                    main=paste0(
                        "Cell Number    ", j,
                        "     Speed vs Forward Persistence Time"),
                    cex.main = 0.8, col.sub="red",
                    sub=paste0(
                        "Spearman's rank correlation coefficient = ", ccPV))
                if (export) grDevices::dev.off()
            }
        }
        ## All cells (FP times    vs Speed)
        allper<-VelFPTable[,1]
        allvel<-VelFPTable[,2]
        for(j in seq(1,length(Object),by=1)){
            allper<-c(allper,VelFPTable[,j+j-1])
            allvel<-c(allvel,VelFPTable[,j+j])
        }
        allper<-allper[!is.na(allper)]
        allvel<-allvel[!is.na(allvel)]
        all.per.vel.table<-data.frame(allper,allvel)
        reg<-stats::lm(allper~allvel)
        #testing the correlation
        c<-stats::cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)
        cc<-unlist(c[4])
        ccP<-round(cc, digits = 3)
        # Speed vs    persistence time    Spearman correlation
        FMResultsTable[6,(length(Object)+1)]<-ccP
        if ( afptPLOT == TRUE){
            if (export) {
                plot_name <- paste0(ExpName,"_FP_Time_VS_Speed_All_Cells.jpg")
                file_path <- file.path(new.fld, plot_name)
                grDevices::jpeg(
                    filename = file_path, width = 4,
                    height = 4, units = 'in', res = 300)
            }
            graphics::plot(
                allvel*60,allper,type="p",pch=16,ylab="FP Time (min)",
                xlab=" Mean Speed during FP time (um/h)",col="black",las=1
            )
            graphics::abline(reg,untf=FALSE,col="red")
            graphics::title(
                "Speed vs FP Time (All cells)", cex.main = 1, col.sub="red",
                sub=paste0("Spearman's rank correlation coefficient = ",ccP))
            if (export) grDevices::dev.off()
        }
        # calculating the Mean.Square.speed for each cell
        for(j in seq(1,length(Object),by=1)){
            MM<-Step
            MM2<-MM-1
            Root.Mean.Square.Speed<-round(
                sqrt(mean(Object[[j]][seq(1,MM2,by=1),11])),digits = 3
            )
            FMResultsTable[7,j]<-Root.Mean.Square.Speed*60
            mean.cosineFP<-round(
                mean(Object[[j]][seq(1,MM2,by=1),20],na.rm = TRUE),digits = 3
            )
            FMResultsTable[8,j]<-mean.cosineFP
            s<-stats::cor.test(
                ~ sqrt(Object[[j]][seq(1,MM2,by=1),11]) +
                    Object[[j]][seq(1,MM2,by=1),20],
                method = "spearman",exact=FALSE) #testing the correlation
            ss<-unlist(s[4])
            VEvsCOSP<-round(ss, digits = 3)
            FMResultsTable[9,j]<-VEvsCOSP
            if ( sfpPLOT == TRUE){
                if (export){
                    plot_name <-    paste0(ExpName,"_FP_VS_Speed_",j,".jpg")
                    file_path <- file.path(new.fld, plot_name)
                    grDevices::jpeg(
                        filename = file_path,width = 4,
                        height = 4, units = 'in',res = 300)
                }
                graphics::plot(
                    sqrt(Object[[j]][seq(1,MM2,by=1),11])*60,
                    Object[[j]][seq(1,MM2,by=1),20],
                    pch=16,type="p",ylab="Forward Persistence Time (min)",
                    xlab=" Instantaneous Speed (um/h)",col=color[j],las=1
                )
                reg<-stats::lm(
                    Object[[j]][seq(1,MM2,by=1),20] ~
                        sqrt(Object[[j]][seq(1,MM2,by=1),11]))
                graphics::abline(reg,untf=FALSE,col="red")
                graphics::title(
                    main=paste0(
                        "Cell Number ", j, " Speed vs Forward Persistence"),
                    cex.main = 0.8, sub=paste0(
                        "spearman's rank correlation coefficient = ",VEvsCOSP
                    ),col.sub="red")
                if (export) grDevices::dev.off()
            }
        }
        outY <- FMResultsTable
    }, error = function(e) {NULL})
    return(outY)
}









