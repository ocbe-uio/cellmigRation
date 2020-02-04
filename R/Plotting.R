#' Method Plotting
#' @title A 2D rose-plot
#' @description Plotting the trajectory data of all cells.
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param Type has to be one of the following: [p, l, b, o]
#' “p”: Points
#' “l”: Lines
#' “b”: Both
#' “o”: Both “overplotted”
#'
#' @return A 2D rose-plot showing the tracks of all cells.
#' @details  The visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).
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
#' plotAllTracks(prepro, ExpName="Test",Type="b")
#' }
#'
plotAllTracks= function(object,ExpName="ExpName",Type="l") {
  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  Len<-length(object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
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

  MinX<-c()
  MaxX<-c()
  MinY<-c()
  MaxY<-c()
  for(j in 1:Len){
    minX=min(object[[j]][1:Step,2])
    minY=min(object[[j]][1:Step,3])
    maxX=max(object[[j]][1:Step,2])
    maxY=max(object[[j]][1:Step,3])
    MinX[j]<-c(minX)
    MaxX[j]<-c(maxX)
    MinY[j]<-c(minY)
    MaxY[j]<-c(maxY)
  }
  RangeX=c(MinX,MaxX)
  RangeY=c(MinY,MaxY)
  plot(object[[1]][1:Step,2],object[[1]][1:Step,3],type=Type,xlab="X (um)",ylab="Y (um)",col=color[1],las=1,xlim=range(RangeX),ylim=range(RangeY),main=ExpName)
  for(n in 2:Len){
    points(object[[n]][,2],object[[n]][,3], type=Type,col=color[n])
    end<-cbind(object[[n]][Step,2],object[[n]][Step,3])
    points(end,pch=16,col=color[n], cex = 1)
  }
  x=c(min(RangeX)-100,max(RangeX)+100)
  y=c(0,0)
  lines(x, y, type='l', col="black")
  x=c(0,0)
  y=c(min(RangeY)-100,max(RangeY)+100)
  lines(x, y, type='l', col="black")

  jpeg(paste0(ExpName,"_All_tracks_plot.jpg"))
  plot(object[[1]][1:Step,2],object[[1]][1:Step,3],type=Type,xlab="X (um)",ylab="Y (um)",col=color[1],las=1,xlim=range(RangeX),ylim=range(RangeY),main=ExpName)
  for(n in 2:Len){
    points(object[[n]][,2],object[[n]][,3], type=Type,col=color[n])
    end<-cbind(object[[n]][Step,2],object[[n]][Step,3])
    points(end,pch=16,col=color[n], cex = 1)
  }
  x=c(min(RangeX)-100,max(RangeX)+100)
  y=c(0,0)
  lines(x, y, type='l', col="black")
  x=c(0,0)
  y=c(min(RangeY)-100,max(RangeY)+100)
  lines(x, y, type='l', col="black")
  dev.off()
  cat("The plot is saved in your directory [use getwd()]","\n")
  }
#'
#'
#'
#'
#'
#' @title A 3D rose-plot of all cells
#' @description Plotting the trajectory data of all cells in 3D.
#'
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#'
#'
#' @return A 3D rose-plot showing the tracks of all cells.
#' @details  The 3D visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).

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
#' plot3DAllTracks(prepro, VS=3, size=2)
#' }
#'
plot3DAllTracks= function(object,VS=3,size=2) {
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )
  plotTable<-data.frame()
  Len<-length(object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
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
  col=c()
  coll=c()
  for(i in 1:Len){
    firstNum=((i*Step)-Step)+1
    RowNum=c(firstNum:(i*Step))
    plotTable[RowNum,1]= object[[i]][1:Step,2]
    plotTable[RowNum,2]=object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
  }
  plot3d(plotTable, col=coll, type="p", size=size, axes=F,xlab=" ", ylab=" ",zlab=" ")
}

#'
#'
#'
#'
#'
#'
#' @title A 3D rose-plot
#' @description Plotting the trajectory data of particular cells in 3D.
#'
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param VS A numeric value of the vertical separator between cells.
#' @param size A numeric value of the point's size.
#' @param cells A numeric vector containing the cell's numbers to be plotted.

#'
#' @return A 3D rose-plot showing the tracks of particular cells.
#' @details  The 3D visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).

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
#' plot3DTracks(prepro, VS=3, size=2,cells=c(1,50,150,250,350))
#' }
#'
plot3DTracks= function(object,VS=3,size=2,cells) {
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )
  plotTable<-data.frame()
  Len<-length(object)
  cat(paste0("The plot contains ",length(cells), " Cells"),"\n")
  Step<-length(object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-rainbow(Len)
  }
  cells=sort(cells)
  col=c()
  coll=c()
  for(i in cells){
    firstNum=((i*Step)-Step)+1
    RowNum=c(firstNum:(i*Step))
    plotTable[RowNum,1]= object[[i]][1:Step,2]
    plotTable[RowNum,2]=object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
    NewplotTable<-plotTable[complete.cases(plotTable),]

  }
  plot3d(NewplotTable, col=coll, type="p", size=size, axes=F,xlab=" ", ylab=" ",zlab=" ")
}
#'
#'
#'
#'
#' @title A graphical display of the track of each cell.
#' @description Plotting the trajectory data of each cell.
#' @param object A list of data frames resulted from runnning the function "PreProcessing()".
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param Type has to be one of the following: [p, l, b, o]
#' “p”: Points
#' “l”: Lines
#' “b”: Both
#' “o”: Both “overplotted”
#' @param FixedField logical(1) Allows generating individual plots with fixed field. Default is TRUE.
#'
#' @return 2D rose-plots of the cells' track Separately.
#' @details  The visualization shows centered trajectories where the starting point of each track is located at the origin of the coordinate system (X=0,Y=0).

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
#' PlotTracksSeparately(prepro, ExpName="Test",Type="b",FixedField=FALSE)
#' }
#'
PlotTracksSeparately= function(object,ExpName="ExpName",Type="l",FixedField=TRUE) {
  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  Len<-length(object)
  cat(paste0(Len," plots will be generated in a folder called:",ExpName,"_Tracks","\n"))
  Step<-length(object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <-rainbow(Len)
  }
  dir.create(paste0(ExpName,"_Tracks"))
  d=getwd()
  setwd(paste0(d,"/",paste0(ExpName,"_Tracks")))

  if ( FixedField == TRUE || FixedField == T){
    MinX<-c()
    MaxX<-c()
    MinY<-c()
    MaxY<-c()
    for(j in 1:Len){
      minX=min(object[[j]][1:Step,2])
      minY=min(object[[j]][1:Step,3])
      maxX=max(object[[j]][1:Step,2])
      maxY=max(object[[j]][1:Step,3])
      MinX[j]<-c(minX)
      MaxX[j]<-c(maxX)
      MinY[j]<-c(minY)
      MaxY[j]<-c(maxY)
    }
    RangeX=c(MinX,MaxX)
    RangeY=c(MinY,MaxY)
    for(n in 1:Len){
      jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"))
      plot(object[[n]][1:Step,2],object[[n]][1:Step,3],type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,xlim=range(RangeX),ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      lines(x, y, type='l', col="black")
      end<-cbind(object[[n]][Step,2],object[[n]][Step,3])
      points(end,pch=16,col=color[n], cex = 1)
      dev.off()
    }
  }else{
    for(n in 1:Len){
      RangeX= object[[n]][1:Step,2]
      RangeY= object[[n]][1:Step,3]
      jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"))
      plot(object[[n]][1:Step,2],object[[n]][1:Step,3],type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,xlim=range(RangeX),ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      lines(x, y, type='l', col="black")
      end<-cbind(object[[n]][Step,2],object[[n]][Step,3])
      points(end,pch=16,col=color[n], cex = 1)
      dev.off()
    }

  }
  setwd(d)
}
#'
#'

