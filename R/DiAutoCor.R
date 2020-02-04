#' Method DiAutoCor
#' @title Direction AutoCorrelation
#'
#' @description The DiAutoCor function automatically compute the angular persistence across several sequantial time intervals.
#' @param object A trajectory data frame organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param sPLOT A logical vector that allows generating individual plots showing the angular persistence across several sequantial time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot showing the angular persistence across several sequantial time intervals of all cells. Default is TRUE.
#' @return Plots and a data frame, which contains six rows: "Cell Number", "Angular Persistence", "Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)", "Stable Direction AutoCorrelation through the track" and "Difference between Mean DA and Intercept DA".
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
#' dac<-DiAutoCor(prepro,TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE)
#' }

DiAutoCor= function(object, TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-DIAutoCorResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-DIAutoCorResults")))

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

  DA.ResultsTable<-data.frame()
  DA.table<-data.frame()

  for(j in 1:length(object)){
    cos.diff<-c()
    LAG<-round(Step*sLAG)     #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:Step, function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        object[[j]][i,14]= object[[j]][i+lag,2]-object[[j]][i,2]    # newdx
        object[[j]][i,15]= object[[j]][i+lag,3]-object[[j]][i,3]    # newdy
        return(object[[j]][i,14:15])
      }))
      object[[j]][1:Step,14:15] <- as.data.frame(res)


      object[[j]][,14:15] <- lapply(object[[j]][,14:15], as.numeric)
      object[[j]][,14][is.na(object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      object[[j]][,15][is.na(object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:Step, function(i){
        object[[j]][i,16]=acos((object[[j]][i+lag,2]-object[[j]][i,2])/sqrt((object[[j]][i+lag,2]-object[[j]][i,2])^2 +(object[[j]][i+lag,3]-object[[j]][i,3])^2)) # to find the abs angle
        return(object[[j]][i,16])
      })
      object[[j]][1:(Step),16] <- res1
      object[[j]][,16][is.na(object[[j]][,16])] <- 0                                # to remove NA and replace it with 0


      res2 <- sapply(1:(Step - lag), function(i){
        if((object[[j]][i+1,15]<0) && (object[[j]][i,15]>=0)||(object[[j]][i+1,15]>=0) && (object[[j]][i,15]<0)){
          object[[j]][i,17]= abs(object[[j]][i+1,16])+abs(object[[j]][i,16])
        }
        if((object[[j]][i+1,15]<0) && (object[[j]][i,15]<0)||(object[[j]][i+1,15]>=0) && (object[[j]][i,15]>=0) ){
          object[[j]][i,17]=object[[j]][i+1,16]-object[[j]][i,16]
        }

        return(object[[j]][i,17])
      })
      object[[j]][1:(Step - lag),17] <- res2

      res3 <- t(sapply(1:(Step - lag), function(i){
        object[[j]][i,17]<-ifelse((object[[j]][i,17])<= (-pi), 2*pi+(object[[j]][i,17]),(object[[j]][i,17]))    # adjusting the ang.diff
        object[[j]][i,17]<-ifelse((object[[j]][i,17])>= pi,(object[[j]][i,17])-2*pi,(object[[j]][i,17]))
        object[[j]][i,18]<-cos(object[[j]][i,17])
        return(object[[j]][i,17:18])
      }))
      object[[j]][1:(Step - lag),17:18] <- as.data.frame(res3)
      object[[j]][,17:18] <- lapply(object[[j]][,17:18], as.numeric)

      for(i in 1:LAG){
        cos.diff[lag]<-mean(object[[j]][1:(Step-lag)-1,18], na.rm=TRUE)  # computing the cosine mean
      }
    }
    DA.table[1:LAG,j]<-cos.diff
    assign(paste0("DA.Cell.",j),cos.diff)
    DA.ResultsTable[1,j]<-j
    DA.ResultsTable[2,j]<-round(cos.diff[1],digits=3)
    lags<-c(1:length(DA.table[,1]))
    lags2<- lags^2
    quadratic.model<-c()
    quadratic.m <-lm(DA.table[,j]~ lags + lags2)
    c<-quadratic.m
    cc<-unlist(c)
    quadratic.model[j]<-cc[1]
    DA.ResultsTable[3,j]<-round(unlist(cc[1]),digits=3)
    ccc<-unlist(cc[1])
    DA.ResultsTable[4,j]<-round(mean(DA.table[1:LAG,j]),digits=3)
    DA.ResultsTable[5,j]<-round((1-(mean(DA.table[1:LAG,j])-unlist(cc[1])))* mean(DA.table[1:LAG,j]),digits=3)
    DA.ResultsTable[6,j]<-round(mean(DA.table[1:LAG,j])-unlist(cc[1]),digits=3)

    timevalues <- seq(1, length(lags), 1)
    predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

    if ( sPLOT == TRUE || sPLOT == T){
      Xaxis<-c(1:LAG)
      Yaxis<-cos.diff
      jpeg(paste0(ExpName,"Direction Autocorrelation.plot.Cell",j,".jpg"))
      plot(Xaxis,Yaxis, type="o",ylim=c(-1,1),xlim=c(0,lag),col=color[j],xlab="Lag",ylab="Cosine",pch=19,las=1,cex=2)
      xx<-c(0,1)
      yy<-c(1,cos.diff[1])
      lines(xx,yy, type='l',col=color[j])
      lines(timevalues, predictedcounts, col = "black", lwd = 3)
      title(main=paste0("Cell Number  ", j, "   DA quadratic model"),col.main="darkgreen",
            sub=paste0(" Intercept of DA quadratic model = ",round(ccc, digits=3)),col.sub="red")
      dev.off()
    }

    object[[j]][,15:18]=0
  }
  RM1<-rowMedians(as.matrix(DA.table),na.rm = TRUE)
  DA.ResultsTable[1,(length(object)+1)]<-"All Cells"
  DA.ResultsTable[2,(length(object)+1)]<-round(RM1[1],digits=3)
  lags<-c(1:length(DA.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  DA.ResultsTable[3,(length(object)+1)]<-round(unlist(cc[1]),digits=3)
  DA.ResultsTable[4,(length(object)+1)]<-round(median(as.numeric(DA.ResultsTable[4,1:length(object)])),digits=3)
  DA.ResultsTable[5,(length(object)+1)]<-round((1-(median(as.numeric(DA.ResultsTable[4,1:length(object)]))-unlist(cc[1])))* median(as.numeric(DA.ResultsTable[4,1:length(object)])),digits=3)
  DA.ResultsTable[6,(length(object)+1)]<-round(median(as.numeric(DA.ResultsTable[4,1:length(object)]))-unlist(cc[1]),digits=3)

  ccc<-unlist(cc[1])
  timevalues <- seq(1, length(lags), 1)
  predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

  if ( aPLOT == TRUE || aPLOT == T){
    Xaxis<-c(1:LAG)
    Yaxis<-RM1

    jpeg(paste0(ExpName,"-Direction Autocorrelation All Cells.jpg"))
    plot(Xaxis,Yaxis, type="o",ylim=c(-1,1),xlim=c(0,lag),col="black",xlab="Lag",ylab="Cosine",pch=19,las=1,cex=2)
    xx<-c(0,1)
    yy<-c(1,RM1[1])
    lines(xx,yy, type='l',col="black")
    lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
    MDA<-round(median(as.numeric(DA.ResultsTable[4,1:length(object)])),digits=2)
    print(MDA)
    abline(h=MDA,col="blue",lwd = 2)
    title(main=paste0("All Cells - DA quadratic model"),col.main="darkgreen",
    sub=paste0(" Intercept of DA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    legend(1, y=-0.82, legend=c("Mean Direction AutoCorrelation","Quadratic model"), col=c("blue","darkgreen"),lty=1, cex=0.8)
    dev.off()
  }

  rownames(DA.ResultsTable)<-c("Cell Number","Angular Persistence","Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)","Stable Direction AutoCorrelation through the track",
                               "Difference between Mean DA and Intercept DA" )
  setwd(d)
  write.csv(DA.ResultsTable, file = paste0(ExpName,"-DA.ResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(DA.ResultsTable)
}

