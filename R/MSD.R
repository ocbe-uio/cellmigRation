#' Method msd
#' @title Mean Square Displacement
#' @description The MSD function automatically compute the mean square displacements across several sequantial time intervals. MSD parameters are used to assess the area explored by cells over time.
#'
#' @param object A trajectory data frame organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param ffLAG A numeric value to be used to get the number of lags for the  Furth formula fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param SlopePlot A logical vector that allows generating individual plots showing the slope of the mean square displacement of the movement of individual cells. Default is TRUE.
#' @param AllSlopesPlot A logical vector that allows generating a plot showing the slope of the mean square displacement of the movement of all cells. Default is TRUE.
#' @param FurthPlot A logical vector that allows generating individual plots fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method per cell. Default is TRUE.
#' @param AllFurthPlot A logical vector that allows generating a plot fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method for all cells. Default is TRUE.
#'
#' @return A data frame and plots
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' data(Trajectory_dataset)
#' df<-Trajectory_dataset
#' prepro<-PreProcessing(df,PixelSize=1.24, TimeInterval=10)
#' msd<-MSD(prepro,sLAG=0.25, ffLAG=0.25)

MSD = function(object,TimeInterval=10,ExpName="ExpName",sLAG=0.25, ffLAG=0.25, SlopePlot=TRUE,AllSlopesPlot=TRUE,FurthPlot=TRUE,AllFurthPlot=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  if ( ! is.numeric(ffLAG) ) stop( "ffLAG has to be a positive number" ) else if ( ffLAG<= 0 ) stop( "ffLAG has to be a positive number" )
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-MSDResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-MSDResults")))
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

  MSDResultsTable<-data.frame()
  MSD.table<-data.frame()                                                          # creating a table that has all the MSDs to be able to compute the mean and sd
  for(j in 1:length(object)){
    meanSD<-c()
    LAG<-round(Step*sLAG)                                                     # number of lags is based on sLAG
    for(lag in 1:LAG){
      res <- sapply(1:Step, function(i){
        object[[j]][i,22]=(((object[[j]][i+lag,2]- object[[j]][i,2])^2)+((object[[j]][i+lag,3]- object[[j]][i,3])^2))
        return(object[[j]][i,22])
      })
      object[[j]][1:Step, 22] <- as.data.frame(res)
      object[[j]][,22][is.na(object[[j]][,22])] <- 0
      meanSD[lag]<-mean(object[[j]][1:(Step-LAG),22])
      object[[j]][,22]=0
    }
    MSDResultsTable[1,j]<-j
    MSDResultsTable[2,j]<-round(meanSD[1],digits=3)
    MSD.table[1:LAG,j]<-meanSD
    Yaxis<-assign(paste0("MSD.Cell.",j),meanSD)
    Xaxis<-c(1:LAG)
    NewrowMeans<-meanSD[1:(round(LAG*sLAG))]                                    # plotting the regression line based on sLAG
    NewXaxis<-Xaxis[1:(round(LAG*sLAG))]
    reg<-lm(NewrowMeans~ NewXaxis)
    reg1<-round(coef(lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
    MSDResultsTable[3,j]<-reg1
    if ( SlopePlot == TRUE || SlopePlot == T){
      jpeg(paste0(ExpName,"-MSD.plot.Cell",j,".jpg"))
      xn <- expression(paste("MSD (um"^2, ")"))
      par(mar=c(5.1, 5.9, 4.1, 1.1), mgp=c(3.5, 0.5, 0), las=0)
      plot(Xaxis,Yaxis, type="p",col=color[j],xlab="Lag",ylab=xn,pch=19,las=1,log="xy",cex=2)
      title(main=paste0("Cell Number  ", j," -  MSD Slope = ",reg1),col.main="black")
      abline(reg,untf=T,col="black")
      dev.off()
    }

  }
  if ( AllSlopesPlot == TRUE || AllSlopesPlot == T){
    RM1<-rowMedians(as.matrix(MSD.table),na.rm = TRUE)
    RSD<-rowSds(as.matrix(MSD.table),na.rm = T)
    xn <- expression(paste("MSD (um"^2, ")"))
    LAG<-round(Step*sLAG)
    Xaxis<-c(1:LAG)
    jpeg(paste0(ExpName,"-MSD.plot All Cells.jpg"))
    xn <- expression(paste("MSD (um"^2, ")"))
    par(mar=c(5.1, 5.5, 4.1, 0.9), mgp=c(3.5, .5, 0), las=0)
    plot(Xaxis,RM1, type="p",col="black",xlab="Lag",ylab=xn,pch=19,las=1,cex=1.5,log="xy")
    NewrowMeans<-RM1[1:(round(LAG*sLAG)+1)]                                  # best fit based on sLAG *sLAG
    NewXaxis<-Xaxis[1:(round(LAG*sLAG)+1)]
    reg<-lm(NewrowMeans~ NewXaxis)
    reg1<-round(coef(lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
    abline(reg,untf=T,col="red")
    title(main=paste0("All Cells -  MSD Slope = ",reg1),col.main="black")
    dev.off()
  }
  for (j in 1: length(MSD.table[1,])){                      # Fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method
    LAG<-round(Step*ffLAG)
    y=MSD.table[1:LAG,j]
    t <- c(1:length(y))
    Data<-c(t,y)
    Data <- matrix(ncol = 2, byrow = F, data =Data)
    colnames(Data) <- c("time","y")
    parms <- c(D=1,P=1)
    parms["D"] <- 1
    parms["P"] <- 1
    Cost <- function(Par) {
      D <- Par[1]
      P <- Par[2]
      out <- cbind(time = t, y = (D*4*(t-P*(1-(exp(-t/P))))))
      return(modCost(obs = Data, model = out))
    }
    Fit<-modFit(p = c(D = 1, P = 1), lower = c(0,0),upper=c(100,100),f = Cost,method = "Nelder-Mead")
    Summary<-summary(Fit)
    MSDResultsTable[4,j]<-round(Fit$par[1],digits=3)
    MSDResultsTable[5,j]<-round(Fit$par[2],digits=3)
    MSDResultsTable[6,j]<-round(Summary$par[1,4],digits=3)
    MSDResultsTable[7,j]<-round(Summary$par[2,4],digits=3)

    if ( FurthPlot == TRUE || FurthPlot == T){
      jpeg(paste0(ExpName,"-MSD N-M bestfit Cell",j,".jpg"))
      plot(Data, pch = 16,col=color[j], cex = 1.5, xlab = "Lags", ylab = "MSD")
      x<-seq(0,LAG,1)
      Model <- function(p, x) return(data.frame(x = x, y = p[1]*4*(x- p[2]*(1-(exp(-x/p[2]))))))
      lines(Model(Fit$par, x),col="black")
      title(main=paste0("Cell Number  ", j,"    Nelder-Mead best-fit"),col.main="black",
            sub=paste0("D = ",round(Fit$par[1],digits=3),
                       "          P = ", round(Fit$par[2],digits=3)),col.sub="red")
      dev.off()
    }


  }

  RM1<-rowMedians(as.matrix(MSD.table),na.rm = TRUE)
  MSDResultsTable[1,(length(object)+1)]<-"All Cells"
  MSDResultsTable[2,(length(object)+1)]<-round(RM1[1],digits=3)

  RM<-c(0,RM1)
  Xaxis<-c(1:round(Step*ffLAG))
  NewrowMeans<-RM1[1:(round(LAG*ffLAG)+1)]                                  # best fit based on ffLAG *ffLAG
  NewXaxis<-Xaxis[1:(round(LAG*ffLAG)+1)]

  reg<-lm(NewrowMeans~ NewXaxis)
  reg1<-round(coef(lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
  MSDResultsTable[3,(length(object)+1)]<-reg1
  y=RM1[1:round(Step*ffLAG)]
  t <- c(1:length(y))
  Data<-c(t,y)
  Data <- matrix(ncol = 2, byrow = F, data =Data)
  colnames(Data) <- c("time","y")
  parms <- c(D=10,P=1)
  parms["D"] <- 1
  parms["P"] <- 1
  Cost <- function(Par) {
    D <- Par[1]
    P <- Par[2]
    out <- cbind(time = t, y = (D*4*(t-P*(1-(exp(-t/P))))))
    return(modCost(obs = Data, model = out))
  }
  Fit<-modFit(p = c(D = 1, P = 1), lower = c(0, 0),upper=c(100,100),f = Cost,method = "Nelder-Mead")
  MSDResultsTable[4,(length(object)+1)]<-round(Fit$par[1],digits=3)
  MSDResultsTable[5,(length(object)+1)]<-round(Fit$par[2],digits=3)
  MSDResultsTable[6,(length(object)+1)]<-round(Summary$par[1,4],digits=3)
  MSDResultsTable[7,(length(object)+1)]<-round(Summary$par[2,4],digits=3)

  if ( AllFurthPlot == TRUE || AllFurthPlot == T){
    jpeg(paste0(ExpName,"-MSD N-M bestfit All Cells.jpg"))
    plot(Data, pch = 16,col="black", cex = 1.5, xlab = "Lags", ylab = "MSD")
    x<-seq(0,LAG,1)
    Model <- function(p, x) return(data.frame(x = x, y = p[1]*4*(x- p[2]*(1-(exp(-x/p[2]))))))
    lines(Model(Fit$par, x),col="red")
    title(main=paste0("All Cells Nelder-Mead best-fit"),col.main="black",
          sub=paste0("D = ",round(Fit$par[1],digits=3),
                     "          P = ", round(Fit$par[2],digits=3)),col.sub="red")
    dev.off()
  }

  rownames(MSDResultsTable)<-c("Cell Number","MSD (lag=1)", "MSD slope", "N-M best fit (Furth) [D]","N-M best fit (Furth) [P]","The significance of fitting D","The significance of fitting P")
  write.csv(MSDResultsTable, file = paste0(ExpName,"-MSDResultsTable.csv"))

  setwd(d)
  write.csv(MSDResultsTable, file = paste0(ExpName,"-MSDResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(MSDResultsTable)
}
