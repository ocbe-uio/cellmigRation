#' Method Persistence and Speed
#' @title Persistence and Speed
#' @description The PerAndSpeed() generates data and plots for persistence and speed.
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param PtSplot A logical vector that allows generating individual plots of persistence time vs speed per cell. Default is TRUE.
#' @param AllPtSplot  A logical vector that allows generating a plot of persistence time vs speed for all cells. Default is TRUE.
#' @param ApSplot A logical vector that allows generating individual plots of angular persistence vs speed per cell. Default is TRUE.
#' @param AllApSplot A logical vector that allows generating a plot of angular persistence vs speed of all cells. Default is TRUE.
#'
#' @return A data frame and plots
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
#' PerSp<-PerAndSpeed(prepro,TimeInterval=10,ExpName="ExpName")
#' }
#'
PerAndSpeed= function(object,TimeInterval=10,ExpName="ExpName",PtSplot=TRUE,AllPtSplot=TRUE,ApSplot=TRUE,AllApSplot=TRUE) {
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  d=getwd()
  dir.create(paste0(ExpName,"-PerResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-PerResults")))

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


  PerResultsTable<-data.frame()                             # creating a table to store the persistence results
  for(j in 1:length(object)){                             # creating values (NA and persistence time )for  persistence    (step to the original)
    MM<-length(object[[j]][,1])
    MM2<-MM-1
    Ptime<-object[[j]][1:MM2,10]
    MeanPerTime<-round(mean(Ptime),digits=2)             # computing mean persistence time
    PerTimLen<-Ptime                                     # computing the number of persistence steps to be used in computing the persistence ratio
    PerTimLen[PerTimLen==0]<-NA
    PerTimLen<-PerTimLen[!is.na(PerTimLen)]
    PerTimLen<-length(PerTimLen)
    PerRatio<-round(PerTimLen/MM2, digits=2)              # computing persistence ratio


    DD.Ptime<-object[[j]][1:MM,10]                     # computing the direction deviating time
    DD.Ptime[DD.Ptime==0]<-1                             # replacing the 0 with 5 and the 5 with 0
    DD.Ptime[DD.Ptime==TimeInterval]<-0
    DD.Ptime[DD.Ptime==1]<-TimeInterval
    MeanDD.PerTime<-round(mean(DD.Ptime),digits=2)

    PerResultsTable[1,j]<-j
    PerResultsTable[2,j]<-MeanPerTime
    PerResultsTable[3,j]<-MeanDD.PerTime
    PerResultsTable[4,j]<-PerRatio
  }

  VelPerTable<-data.frame()               # creating a table to store the mean velocity with correspondence with the persistence time
  for(j in 1:length(object)){
    MM=length(object[[j]][,1])
    Ptime<-object[[j]][1:MM,10]
    Ptime0<-c(0,Ptime)                      #adding a "0" value in the beginning
    Ptime0[Ptime0==TimeInterval]<-NA        #replacing the "TimeInterval" values with NAs
    Ptime0TF<- !is.na(Ptime0)
    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval <- rowN[Ptime0TF]

    fillIdx <- cumsum(Ptime0TF)
    s<-Tval[fillIdx]                     #replacing the NAs with the previous true value
    justTrueVal<-Tval
    s[Ptime0TF]=0                        #replacing the true values with 0
    finalID<-s[-1]                       # removing the first added value

    Vel<-object[[j]][1:MM,11]
    Per<-object[[j]][1:MM,10]
    Per[Per==0]<-NA
    tab<-data.frame(Vel,Per,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]         # to exclude the last row since it has no velocity
    tabb<-split(tab, finalIDD)     #to avoid taking the last row

    meanVEL<-c()
    Per.time<-c()
    res1 <- sapply(1:length(tabb), function(i){
      meanVEL= round(sqrt(mean(tabb[[i]][,1])),digits=2)
      return(meanVEL)
    })
    meanVEL= res1
    meanVEL=meanVEL[-1]

    res2 <- sapply(1:length(tabb), function(i){
      Per.time=sum(tabb[[i]][,2])
      return(Per.time)
    })
    Per.time= res2
    Per.time=Per.time[-1]
    w<-which.max(Per.time)
    PerResultsTable[5,j]<-Per.time[w]


    Ptime<-object[[j]][1:MM,10]         # computing the "meanVel.for0"
    Ptime00<-c(1,Ptime)                   #adding a "1" value in the beginning  this value should not be 0 because we will not catch 0 persistence if it is in the beginning.
    Ptime00[Ptime00==0]<-NA               #replacing the "0" values with NAs
    Ptime00TF<- !is.na(Ptime00)
    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval0 <- rowN[Ptime00TF]
    #length(Tval0)

    TTT<-c()
    res3 <- sapply(1:(length (Tval0)-1), function(i){
      TTT<- (Tval0[i+1]- Tval0[i])-1
      return(TTT)
    })
    TTT=res3
    TTT[TTT==0]<-NA
    F.ID.for.0.per<-TTT[!is.na(TTT)]
    Final.ID.for.0.per<-c()
    for (i in 1:length(F.ID.for.0.per)){
      Final.ID.for.0.per<-c(Final.ID.for.0.per,rep(i,(F.ID.for.0.per[i])))
    }

    Vel<-object[[j]][1:MM,11]
    Per<-object[[j]][1:MM,10]
    Per[Per==0]<-NA
    tab<-data.frame()
    tab<-data.frame(Vel,Per,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]
    tabb<-split(tab, finalIDD)     #to avoid taking the last row
    t<-tabb[[1]]

    tabb1<-cbind(t,Final.ID.for.0.per)
    tabbb<-split(tabb1, Final.ID.for.0.per)
    meanVel.for0<-c()

    res4 <- sapply(1:length(tabbb), function(i){
      meanVel.for0<-round(sqrt(mean(tabbb[[i]][,1])),digits=2)
      return(meanVel.for0)
    })

    meanVel.for0=res4
    zerooPrep<-rep(0,length(meanVel.for0))
    Per.time1<-c(zerooPrep,Per.time)
    meanVEL1<-c(meanVel.for0,meanVEL)
    colNumP<-j+j-1
    colNumV<-j+j
    rowNum<-c(1:length(meanVEL1))
    VelPerTable[rowNum,colNumP]<-Per.time1
    VelPerTable[rowNum,colNumV]<-meanVEL1
    PT<-VelPerTable[,j+j]
    c<-cor.test( ~ VelPerTable[,j+j-1]+ PT, method = "spearman",exact=FALSE)                 #testing the correlation
    cc<-unlist(c[4])
    ccPV<-round(cc, digits = 3)
    PerResultsTable[6,j]<-ccPV               # Speed vs  persistence time  Spearman correlation


    if ( PtSplot == TRUE || PtSplot == T){
      jpeg(paste0(ExpName," Persist Time vs Speed",j,".jpg"))
      plot(VelPerTable[,j+j],VelPerTable[,j+j-1],pch=16,type="p",ylab="Persistence Time (min)",xlab=" Mean Speed during persistence time  (um/min)",col=color[j],las=1)
      reg<-lm(PT~VelPerTable[,j+j-1])
      #abline(reg,untf=F,col="black")
      title(main=paste0("Cell Number  ", j,"   Speed vs Persistence Time"),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",ccPV),col.sub="red")
      dev.off()

    }

  }

  ## All cells (Persistence times  vs Speed)
  allper<-VelPerTable[,1]
  allvel<-VelPerTable[,2]
  for(j in 1:length(object)){
    allper<-c(allper,VelPerTable[,j+j-1])
    allvel<-c(allvel,VelPerTable[,j+j])
  }
  allper<-allper[!is.na(allper)]
  allvel<-allvel[!is.na(allvel)]
  allvel<-(allvel/TimeInterval)*60
  all.per.vel.table<-data.frame(allper,allvel)
  reg<-lm(allper~allvel)
  c<-cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)                 #testing the correlation
  cc<-unlist(c[4])
  ccP<-round(cc, digits = 3)
  PerResultsTable[6,(length(object)+1)]<-ccP                                          # Speed vs  persistence time  Spearman correlation

  if ( AllPtSplot == TRUE || AllPtSplot == T){
    jpeg(paste0(ExpName,"_Persist_Time_vs_Speed-All_Cells.jpg"))
    plot(allvel,allper,type="p",pch=16,ylab="Persist_Time (min)",xlab=" Mean Speed during persistence time (um/h)",col="black",las=1)
    abline(reg,untf=F,col="red")
    title("Speed vs Persist Time (All cells)",cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccP),col.sub="red")
    dev.off()

  }

  for(j in 1:length(object)){    # calculating the Mean.Square.velocity for each cell
    MM<-length(object[[j]][,1])
    MM2<-MM-1
    Root.Median.Square.Speed<-round(sqrt(median(object[[j]][1:MM2,11])),digits = 3)*60
    PerResultsTable[7,j]<-Root.Median.Square.Speed

    wma<-which.max(sqrt(object[[j]][,11]))
    wmi<-which.min(sqrt(object[[j]][1:MM2,11]))

    PerResultsTable[8,j]<-round(sqrt(object[[j]][wma,11]),digits=3)* 60
    PerResultsTable[9,j]<-round(sqrt(object[[j]][wmi,11]),digits=3)* 60

    mean.cosineP<-round(mean(object[[j]][1:(MM2-1),9],na.rm = TRUE),digits = 3)
    PerResultsTable[10,j]<-mean.cosineP
    s<-cor.test( ~ sqrt(object[[j]][1:MM2,11])+ object[[j]][1:MM2,9], method = "spearman",exact=FALSE)                 #testing the correlation
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    PerResultsTable[11,j]<-VEvsCOSP

    data<-object[[j]][1:MM2,2:3]
    data1<-object[[j]][1:round(MM2/4),2:3]
    data2<-object[[j]][round(MM2/4):round(MM2/2),2:3]
    data3<-object[[j]][round(MM2/2):round(MM2*3/4),2:3]
    data4<-object[[j]][round(MM2*3/4):MM2,2:3]

    ch  <- chull(data)
    ch1 <- chull(data1)
    ch2 <- chull(data2)
    ch3 <- chull(data3)
    ch4 <- chull(data4)

    coords  <- data[c(ch, ch[1]), ]  # closed polygon
    coords1 <- data1[c(ch1, ch1[1]), ]  # closed polygon
    coords2 <- data2[c(ch2, ch2[1]), ]  # closed polygon
    coords3 <- data3[c(ch3, ch3[1]), ]  # closed polygon
    coords4 <- data4[c(ch4, ch4[1]), ]  # closed polygon

    p = Polygon(coords)
    p1 = Polygon(coords1)
    p2 = Polygon(coords2)
    p3 = Polygon(coords3)
    p4 = Polygon(coords4)

    EmptyArea=abs(p@area -(p1@area + p2@area + p3@area + p4@area))
    SegmentedCA<-p1@area + p2@area + p3@area + p4@area

    PerResultsTable[12,j]=round(p@area,digits=3)
    PerResultsTable[13,j]=round(SegmentedCA,digits=3)
    PerResultsTable[14,j]=round(EmptyArea,digits=3)


    PerResultsTable[15,j]=round(sum(abs(object[[j]][,8]))/6.28,digits=3)
    PerResultsTable[16,j]=round(sum(abs(object[[j]][,8]))/6.28,digits=3) - abs(round(sum(object[[j]][,8])/6.28,digits=3))


    if ( ApSplot == TRUE || ApSplot == T){
      jpeg(paste0(ExpName,"_Angular_Persistence_vs_Speed",j,".jpg"))
      Speed=(sqrt(object[[j]][1:MM2,11])/TimeInterval)*60
      plot(Speed,object[[j]][1:MM2,9],pch=16,type="p",ylab="Angular Persistence (cosine)",xlab=" Instantaneous Speed (um/h)",col=color[j],las=1, xlim=c(0,3))
      reg<-lm(object[[j]][1:MM2,9]~Speed)
      abline(reg,untf=F,col="black")
      title(main=paste0("Cell Number  ", j," Instantaneous Speeds vs Angular Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
      dev.off()
    }

    PoCos<- subset(object[[j]][1:(MM2-1),9],object[[j]][1:(MM2-1),9]>0)
    NeCos<- subset(object[[j]][1:(MM2-1),9],object[[j]][1:(MM2-1),9]<=0)

    PerResultsTable[17,j]<-round(median(PoCos,na.rm = TRUE),digits = 3)
    PerResultsTable[18,j]<-round(median(NeCos,na.rm = TRUE),digits = 3)

    AvSp<-(object[[j]][1:length(object[[j]][,1])-1,6]/TimeInterval)*60

    PerResultsTable[19,j]<-round(median(AvSp,na.rm = TRUE),digits = 3)
    s=summary(AvSp)
    PerResultsTable[20,j]<-round((s[5]-s[2])/s[3],digits=4)
    PerResultsTable[21,j]<-round(mean(AvSp),digits = 3)
    PerResultsTable[22,j]<-round(sd(AvSp),digits = 3)
  }

  #(all cells)  persistence vs inst.speed
  MM<-length(object[[1]][,1])
  MM2<-MM-1
  cosine.P<-data.frame()
  for(j in 1:length(object)){              # creating values for  cosine.P  based on rel.ang.P
    M<- object[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      object[[j]][i,9]<-cos(object[[j]][i,8])
      return(object[[j]][i,9])
    })
    object[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-object[[j]][,9]
  }

  RM<-round(rowMedians(as.matrix(cosine.P[1:MM2,]),na.rm = TRUE),digits=3)
  Speed<-data.frame()
  for (j in 1:length(object)){    # calculating the Mean.Square.velocity for each cell
    Speed[1:MM2,j]<-round(sqrt(object[[j]][1:MM2,11]),digits = 3)
  }
  RowmeanSpeed<-round(rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  VEvsCOSP<-round(ss, digits = 3)
  PerResultsTable[11,(length(object)+1)]<-VEvsCOSP

  if ( AllApSplot == TRUE || AllApSplot == T){
    jpeg(paste0(ExpName," All_Cells_Average_Angular_Persistence_vs_Average_Speed.jpg"))
    MS<-max(RowmeanSpeed)*60
    plot(RowmeanSpeed*60,RM,pch=16,type="p",ylab="Average Instantaneous Angular Persistence (cosine)",xlab=" Average Instantaneous Speed (um/h)",col="black",las=1,xlim=c(0,MS))
    NewSpeed=RowmeanSpeed*60
    reg<-lm(RM~NewSpeed)
    abline(reg,untf=F,col="red")
    title(main=paste0("All Cells Instantaneous Speed vs Angular Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
    dev.off()
  }

  rownames(PerResultsTable)<-c("Cell Number","Mean Persist Time (min)","Mean persist Deviating Time (min)","Persistence Ratio",
                               "Maximum Persistence period (min)","Persistence Time vs Speed (SCC)","RMSS (um per h)","Maximum Speed (um per h)","Minimum Speed (um per h)",
                               "Mean Angular Persistence (cosine)","Instantaneous Speed vs Angular Persistence (SCC)","Covered Area (um2)","Segmented Covered Area (um2)","Empty Area (um2)","Number of complete rotations",
                               "Number of canceled rotations","Mean Persistence Angle (cosine)","Mean Deviating Angle (cosine)","Median Speed","Speed QBCV",
                               "Mean Speed (um per h)","Speed standard deviation (um per h)")



  RM1<-round(rowMedians(as.matrix(PerResultsTable),na.rm = TRUE),digits=3)
  PerResultsTable[c(2:5,7:10,12:22),(length(object)+1)]<-RM1[c(2:5,7:10,12:22)]

  RMSS<-as.numeric(PerResultsTable[7,1:length(PerResultsTable[1,])-1])
  jpeg(paste0(ExpName,"_RMSS_profile_of_all_cells.jpg"))
  cells<-c(1:(length(PerResultsTable[1,])-1))
  MS<-max(RMSS)
  plot(cells,RMSS,pch=16,type="o",ylab = 'RMSS(um/h)',xlab = 'Cells',las=1,ylim = c(0, MS))
  title(main="RMSS Profile of all cells",cex.main = 1)
  abline(h=median(RMSS[which(!is.na(RMSS))]),col="red")
  abline(h=mean(RMSS[which(!is.na(RMSS))]),col="blue")
  legend(1, y=200, legend=c("Mean RMSSs","Median RMSSs"), col=c("blue","red"),lty=1, cex=0.8)
  dev.off()

  jpeg(paste0(ExpName,"_RMSS_PiolinPlot_of_all_cells.jpg"))
  plot(1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '', ylab = 'RMSS(um/h)', xaxt = 'n',las=1)
  title("RMSS of all cells",cex.main = 1)
  vioplot(RMSS, at = 1, add = T, col = "gray")
  dev.off()


  SPEED<-as.numeric(PerResultsTable[19,1:length(PerResultsTable[1,])-1])
  jpeg(paste0(ExpName,"_Speed_profile_of_all_cells.jpg"))
  cells<-c(1:(length(PerResultsTable[1,])-1))
  MS<-max(SPEED)
  plot(cells,SPEED,pch=16,type="o",ylab = 'Speed(um/h)',xlab = 'Cells',las=1,ylim = c(0, MS))
  title(main="Speed Profile of all cells",cex.main = 1)
  abline(h=median(SPEED[which(!is.na(SPEED))]),col="red")
  abline(h=mean(SPEED[which(!is.na(SPEED))]),col="blue")
  legend(1, y=200, legend=c("Mean Speed","Median Speed"), col=c("blue","red"),lty=1, cex=0.8)
  dev.off()

  jpeg(paste0(ExpName,"_Speed_PiolinPlot_of_all_cells.jpg"))
  plot(1, 1, xlim = c(0, 2), ylim = c(0, MS), type = 'n', xlab = '', ylab = 'Speed(um/h)', xaxt = 'n',las=1)
  title("Speed of all cells",cex.main = 1)
  vioplot(SPEED, at = 1, add = T, col = "gray")
  dev.off()


  SPEED<-as.numeric(PerResultsTable[21,1:length(PerResultsTable[1,])-1])
  PerR<- as.numeric(PerResultsTable[4,1:length(PerResultsTable[1,])-1])
  MS<-max(SPEED)
  plot(SPEED,PerR,pch=16,type="p",xlab="Average Instantaneous Speed (um/h)",ylab=" Persistence Ratio",col="black",las=1,xlim=c(0,MS),ylim=c(0,1))
  reg<-lm(PerR~SPEED)
  abline(reg,untf=F,col="red")
  RowmeanSpeed<-round(rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-cor.test( ~ SPEED+ PerR, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  SCC<-round(ss, digits = 3)
  title(main=paste0("All Cells Average Speed vs Persistence Ratio "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",SCC),col.sub="red")

  PerResultsTable[1,(length(object)+1)]<-"All Cells"

  setwd(d)
  write.csv(PerResultsTable, file = paste0(ExpName,"-PerResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-PerResultsTable.csv" ),"in your directory [use getwd()]","\n")
  return(PerResultsTable)
}




