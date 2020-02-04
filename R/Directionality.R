#' Method DR
#' @title Directionality Table
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straightline length between the start point and the endpoint of the migration trajectory,
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#'
#' @return A data frame with nine rows: "Cell Number","Directionality Ratio","Mean Cumulative Directionality Ratio","Stable Directionality Ratio", "Number of returns","Min CumDR","Location of Min CumDR, Steps with less CumDR than DR","Directional Persistence".
#'
#' @details  Directionality Ratio and Directional persistence
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' df<-Trajectory_dataset
#' prepro<-PreProcessing(df,PixelSize=1.24, TimeInterval=10)
#' DR=DiRatio(prepro, ExpName="Test")
#' }
#'
DiRatio = function(object,TimeInterval=10,ExpName="ExpName") {
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  d=getwd()
  dir.create(paste0(ExpName,"-DRResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-DRResults")))
  Step<-length(object[[1]][,1])
  DRResultsTable<-data.frame()
  DIR.RATIO<-c()
  for(j in 1:length(object)){                             # calculating the cumsum of distance for each cell
    MM<-Step
    MM2<-MM-1
    end<-cbind(object[[j]][MM,2],object[[j]][MM,3])       # finding the cordinates of the final point in the track.
    start<-cbind(0,0)
    final.dis=dist2(start, end)                               # calculating the final distance
    total.dis=sum(object[[j]][1:MM2,6])                     # calculating the total distance
    Dir.Ratio=round(final.dis/total.dis ,digits = 3)          # calculating the Dir.Ratio
    mean.Dir.Ratio<-round(mean(object[[j]][1:MM2,13],na.rm = TRUE) ,digits = 3)
    StableDR<- (1-(mean.Dir.Ratio-Dir.Ratio))* mean.Dir.Ratio
    StableDR<-round(StableDR,digits = 3)

    DRResultsTable[1,j]<- j
    DRResultsTable[2,j]<-Dir.Ratio
    DRResultsTable[3,j]<-mean.Dir.Ratio
    DRResultsTable[4,j]<- StableDR

  }


  for(j in 1:length(object)){                                   #### Adding min CumDR  and number of angles greater than .75
    MM<-Step
    MM2<-MM-1
    p1<-object[[j]][1:MM2,9]
    returns<-subset(p1,p1<(-0.87))                            # greater than 150 degrees
    DRResultsTable[5,j]<-length(returns)

    p2<-object[[j]][1:MM2,13]
    w<-which.min(p2)
    DRResultsTable[6,j]<-round(p2[w], digits=3)
    DRResultsTable[7,j]<-w
    DR<-as.numeric(DRResultsTable[2,j])
    lessThanMINcumdr<-subset(p2,p2<DR)                        # number of the steps that have a cumdr less than the final dr
    DRResultsTable[8,j]<-length(lessThanMINcumdr)

    Ptime<-object[[j]][1:MM2,10]
    PerTimLen<-Ptime                                     # computing the number of persistence steps to be used in computing the persistence ratio
    PerTimLen[PerTimLen==0]<-NA
    PerTimLen<-PerTimLen[!is.na(PerTimLen)]
    PerTimLen<-length(PerTimLen)
    PerRatio<-round(PerTimLen/MM2, digits=2)
    DRResultsTable[9,j]<- PerRatio + DRResultsTable[4,j]

  }


  rownames(DRResultsTable)<-c("Cell Number","Directionality Ratio","Mean Cumulative Directionality Ratio","Stable Directionality Ratio",
                              "Number of returns","Min CumDR",paste0("Location of Min CumDR (out of ",Step-1,")"),"Steps with less CumDR than DR","Directional Persistence")

  RM1<-round(rowMedians(as.matrix(DRResultsTable),na.rm = TRUE),digits=3)
  DRResultsTable[,(length(object)+1)]<-RM1
  DRResultsTable[1,(length(object)+1)]<-"All Cells"
  setwd(d)
  write.csv(DRResultsTable, file = paste0(ExpName,"-DRResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-DRResultsTable.xlsx" ),"in your directory [use getwd()]","\n")
  return(DRResultsTable)
}
#'
#'
#' Method DR
#' @title Directionality Ratio plots
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straightline length between the start point and the endpoint of the migration trajectory,
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#'
#' @return Directionality Ratio plots
#'
#' @details  Directionality Ratio
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' df<-Trajectory_dataset
#' prepro<-PreProcessing(df,PixelSize=1.24, TimeInterval=10)
#' DR=DiRatio(prepro, ExpName="Test")
#' DiRatio.Plot(prepro, ExpName="Test")
#' }
#'
DiRatio.Plot = function(object,TimeInterval=10,ExpName=ExpName) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  Step<-length(object[[1]][,1])
  color <-c()
  Len<-length(object)
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-rainbow(Len)
  }
  d=getwd()

  if (file.exists(paste0(d,"/",paste0(ExpName,"-DRResults")))){
    setwd(paste0(d,"/",paste0(ExpName,"-DRResults")))

  } else {
    stop("Please run the DiRatio() first")
  }

  DIR.RATIO.AllCells<-data.frame()
  for(j in 1:length(object)){
    MM<-Step
    MM2<-MM-1
    Time<-c(1:MM2)
    t<-60/TimeInterval
    xM<-MM2*TimeInterval
    xMM<-round(xM/60)
    xMMM<-c(0:xMM)
    xMMM1<-(c(0:xMM)*(60/TimeInterval))
    Xaxis<-c(1:MM2)
    jpeg(paste0(ExpName,"-D.R.plot-",j,".jpg"))
    p<-plot(Time,object[[j]][1:MM2,13], type="l",col=color[j],xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1)
    axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
    title(main=paste0("Cell Number  ", j),col.main=color[j])
    dev.off()
    DIR.RATIO.AllCells[1:MM2,j]<-object[[j]][1:MM2,13]
  }

  mycolblue <- rgb(0, 0, 255, maxColorValue = 255, alpha = 100, names = "blue")    #transparent color
  mean.DIR.RATIO.AllCells<-rowMedians(as.matrix(DIR.RATIO.AllCells),na.rm = T)
  SD.DIR.RATIO.AllCells<-rowSds(as.matrix(DIR.RATIO.AllCells),na.rm = T)    # is a function in matrixStats
  meanSDp<-mean.DIR.RATIO.AllCells+SD.DIR.RATIO.AllCells
  meanSDp=ifelse(meanSDp>=1,1,meanSDp)
  meanSDn<-mean.DIR.RATIO.AllCells-SD.DIR.RATIO.AllCells
  meanSDpn<-c(meanSDp,meanSDn)

  jpeg(paste0(ExpName,"directionality ratio for all cells.jpg"))
  p<-plot(Time,mean.DIR.RATIO.AllCells, type="l",col="black",xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1)
  axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
  lines(Time,meanSDp, col="black")
  lines(Time,meanSDn, col="black")
  polygon(c(Time, rev(Time)), c(meanSDp, rev(meanSDn)),col = mycolblue , border = NA)
  title(main="Directionality Ratio - All Cells",col.main="black")
  dev.off()
  setwd(d)
  cat("Plots are saved in a folder in your directory [use getwd()]","\n")
}





