#' Method DiAutoCor
#' @title Velocity AutoCorrelation
#'
#' @description The VeAutoCor function automatically compute the changes in both speed and direction across several sequantial time intervals.
#' @param object A trajectory data frame organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param sPLOT A logical vector that allows generating individual plots showing the velocity across several sequantial time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot showing the velocity across several sequantial time intervals of all cells. Default is TRUE.
#' @return Plots and a data frame, which contains six rows: "Cell Number", .
#'
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
#' vac<-VeAutoCor(prepro,TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE)
#' }

VeAutoCor= function(object, TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-VeAutoCorResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-VeAutoCorResults")))

  Len<-length(object)
  Step<-length(object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- rainbow(Len)
  }

  VA.ResultsTable<-data.frame()
  VAC.table<-data.frame()     # creating a table that has all the VAC to be able to compute the mean
  VAC.first.value<-c()        # to save the fist VAC non-normalized value for all cells
  VAC.second.value<-c()        # to save the second VAC non-normalized value for all cells

  for(j in 1:length(object)){
    meanVAC<-c()
    LAG<-round(Step*sLAG)              #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:(Step - 1), function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        object[[j]][i,14]= object[[j]][i+lag,2]-object[[j]][i,2]    # newdx
        object[[j]][i,15]= object[[j]][i+lag,3]-object[[j]][i,3]    # newdy
        return(object[[j]][i,14:15])
      }))
      object[[j]][1:(Step -1),14:15] <- as.data.frame(res)
      object[[j]][,14:15] <- lapply(object[[j]][,14:15], as.numeric)
      object[[j]][,14][is.na(object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      object[[j]][,15][is.na(object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:(Step - lag), function(i){                                               # starting from 2 to exclude the first cosine which is always 1.
        object[[j]][i,23]=((object[[j]][i,14]* object[[j]][i+lag,14])+ (object[[j]][i,15]* object[[j]][i+lag,15]))/ ((lag*T)^2)
        return(object[[j]][i,23])
      })
      object[[j]][1:(Step - lag),23] <- res1
      meanVAC[lag]<-mean(object[[j]][1:(Step - lag),23])

    }
    VAC.first.value[j]<-meanVAC[1]

    NORMmeanVAC<-meanVAC/meanVAC[1]
    VAC.table[1:LAG,j]<-meanVAC
    assign(paste0("VAC.Cell.",j),meanVAC)
    VAC.second.value[j]<-NORMmeanVAC[2]

    VA.ResultsTable[1,j]<-j
    VA.ResultsTable[2,j]<-round(meanVAC[1],digits=3)    # VA (lag =1)
    VA.ResultsTable[3,j]<-round(NORMmeanVAC[2],digits=3)    # VA (lag =2)

    lags<-c(1:length(VAC.table[,1]))
    lags2<- lags^2
    quadratic.model<-c()
    quadratic.m <-lm(VAC.table[,j]~ lags + lags2)
    c<-quadratic.m
    cc<-unlist(c)
    quadratic.model[j]<-cc[1]

    ccc<-unlist(cc[1])
    VA.ResultsTable[4,j]<-round(ccc,digits=3)
    VA.ResultsTable[5,j]<-round(mean(meanVAC),digits=3)      # mean VA (all lags)
    timevalues <- seq(1, length(lags), 1)
    predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

    if (sPLOT == TRUE || sPLOT == T){
      Xaxis<-c(1:LAG)
      Yaxis<-meanVAC
      jpeg(paste0(ExpName,"Velocity Autocorrelation.plot.Cell",j,".jpg"))
      plot(Xaxis,Yaxis, type="o",ylim=range(meanVAC),xlim=c(0,lag),col=color[j],xlab="Lag",ylab="Velocity  Autocorrelation",pch=19,las=1,cex=2)
      lines(timevalues, predictedcounts, col = "black", lwd = 3)
      title(main=paste0("Cell Number  ", j, "   VA quadratic model"),col.main="darkgreen",
            sub=paste0(" Intercept of VA quadratic model = ",round(ccc, digits=3)),col.sub="red")
      dev.off()
    }

    object[[j]][,15:16]=0
  }

  RM1<-rowMedians(as.matrix(VAC.table),na.rm = TRUE)
  VA.ResultsTable[1,(length(object)+1)]<-"All Cells"
  VA.ResultsTable[2,(length(object)+1)]<-round(median(VAC.first.value),digits=3)
  VA.ResultsTable[3,(length(object)+1)]<-round(median(VAC.second.value),digits=3)

  lags<-c(1:length(VAC.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  VA.ResultsTable[4,(length(object)+1)]<-round(unlist(cc[1]),digits=3)
  VA.ResultsTable[5,(length(object)+1)]<-round(median(as.numeric(VA.ResultsTable[5,1:length(object)])),digits=3)

  ccc<-unlist(cc[1])
  timevalues <- seq(1, length(lags), 1)
  predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

  if ( aPLOT == TRUE || aPLOT == T){
    Xaxis<-c(1:LAG)
    Yaxis<-RM1
    jpeg(paste0(ExpName,"-Velocity Autocorrelation All Cells.jpg"))
    plot(Xaxis,Yaxis, type="o",ylim=range(RM1),xlim=c(0,lag),col="black",xlab="Lag",ylab="Velocity  Autocorrelation",pch=19,las=1,cex=2)
    lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
    title(main=paste0("All Cells - VA quadratic model"),col.main="darkgreen",
          sub=paste0(" Intercept of VA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    dev.off()
  }

  rownames(VA.ResultsTable)<-c("Cell Number","Velocity AutoCorrelation (lag=1)","2nd normalized Velocity AutoCorrelation","Intercept of VA quadratic model","Mean Velocity AutoCorrelation (all lags)")

  setwd(d)
  write.csv(VA.ResultsTable, file = paste0(ExpName,"-VA.ResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(VA.ResultsTable)

}

