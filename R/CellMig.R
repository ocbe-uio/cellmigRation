#' @title The CellMig Class
#'
#' @description The CellMig class is the central object storing all information for both random migration (RM) and wound scratch assay (WSA).
#' It comprises 14 slots
#'
#' @slot trajdata The raw trajectory data matrix organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @slot adjDS A data frame of the trajectory data passed from the WSAprep function.
#' @slot cellpos A binary vector showing on which side of the wound cells are located. "0" reffers to a cell located above the wound whereas "1" reffers to a cell located below the wound.
#' @slot parE A numeric vector contains estimations for the imageH, woundH, upperE and lowerE.
#' @slot preprocessedDS list object of data frames, each data frame shows the trajectories of a single cell.
#' @slot DRtable A data frame of the results of running the DiRatio() function.
#' @slot MSDtable A data frame of the results of running the MSD() function.
#' @slot PerAanSpeedtable A data frame of the results of running the PerAndSpeed() function.
#' @slot DACtable A data frame of the results of running the DiAutoCor() function.
#' @slot VACtable A data frame of the results of running the VeAutoCor() function.
#' @slot ForMigtable A data frame of the results of running the ForwardMigration() function.
#' @slot FMItable A data frame of the results of running the FMI() function.
#' @slot results A data frame of all the results.
#' @slot parCor A data frame for Parameters Correlation.
#'
#' @name CellMig
#' @rdname CellMig
#' @aliases CellMig-class
#' @exportClass CellMig
#' @export
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
CellMig <- setClass("CellMig", slots = c(trajdata= "data.frame",adjDS= "data.frame",cellpos = "vector", parE= "vector", preprocessedDS= "list", DRtable= "data.frame", MSDtable= "data.frame", PerAanSpeedtable= "data.frame", DACtable= "data.frame", VACtable= "data.frame", ForMigtable= "data.frame", FMItable= "data.frame",results="data.frame", parCor="matrix"))

#' validity function for CellMig
#'
#' @param object An CellMig object, which is a trajectory data frame organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @name CellMig
#' @export
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
setValidity("CellMig",
            function(object) {
              msg <- NULL
              if ( nrow(object@trajdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@trajdata) < 4 ){
                msg <- c(msg, "input data must have four columns")
              }
              if (is.null(msg)) TRUE
              else msg
            }
)


setMethod("initialize",
          signature = "CellMig",
          definition = function(.Object, trajdata){
            .Object@trajdata <- trajdata
            validObject(.Object)
            return(.Object)
          }
)



#' @title Data preprocessing for random migration (RM)

#'
#' @description This function allows preprocessing of the trajectory data from random migration (RM) experiments.
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param PixelSize A numeric value of the physical size of a pixel.
#'
#'
#' @return An CellMig class object with preprocessed data.
#' @examples
#' data(Trajectory_dataset)
#' rmTD <- CellMig(Trajectory_dataset)
#' \dontrun{
#' rmTD <- rmPreProcessing(rmTD)
#' }
#' @export
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
rmPreProcessing = function(object,PixelSize=1.24,TimeInterval=10) {
  msg <- NULL
  if ( ! is.data.frame(object@trajdata)){
    msg <- c(msg, "input data must be data.frame")
  }
  if ( ! is.numeric(PixelSize)) stop( "PixelSize has to be a positive number" ) else if ( PixelSize<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  object@adjDS <- object@trajdata
  df<-object@adjDS
  df<-df[,1:3]                                        # Removing the unnecessary columns
  L<-length(df[,1])
  df[,4:26]<-rep(0,L)
  df[,27]<-rep(NA,L)                                  # to be used for migration type
  colnames(df)<-c("ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P","Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY","New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F","Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type")
  ID_split <- split(df, df$ID)                        #Splitting the data frame based on the ID
  cat("This dataset contains: ",length(ID_split),"Cells","\n")

  for(j in 1:length(ID_split)){                        # Having the ID =group order
    ID_split[[j]][1]=j
  }

  for(j in 1:length(ID_split)){                        # adjusting x and y (starting from 0 & being multiplied by H)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(2:MM, function(i){
      ID_split[[j]][i,4]=PixelSize*((ID_split[[j]][i,2])-( ID_split[[j]][1,2]))      # x2-x1
      ID_split[[j]][i,5]=PixelSize*(( ID_split[[j]][1,3])-(ID_split[[j]][i,3]))      # y2-y1
      return(ID_split[[j]][i,4:5])
    }))
    ID_split[[j]][2:MM,4:5] <- as.data.frame(res)
    ID_split[[j]][,4:5] <- lapply(ID_split[[j]][,4:5], as.numeric)
  }



  for(j in 1:length(ID_split)){                    # removing the old x and y
    ID_split[[j]]=ID_split[[j]][-2]            # removing x column
    ID_split[[j]]=ID_split[[j]][-2]            # removing the y column [-2] is used because x column is gone.
  }



  for(j in 1:length(ID_split)){                    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(1:MM, function(i){
      ID_split[[j]][i,4]=(ID_split[[j]][i+1,2])-( ID_split[[j]][i,2])                                         # creating values for dx
      ID_split[[j]][,4][is.na(ID_split[[j]][,4])] <- 0
      ID_split[[j]][i,5]= ( ID_split[[j]][i+1,3])-(ID_split[[j]][i,3])                                        # creating values for dy
      ID_split[[j]][,5][is.na(ID_split[[j]][,5])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,6]=sqrt((ID_split[[j]][i,4])^2 + (ID_split[[j]][i,5])^2)                                # creating values for dis
      ID_split[[j]][i,7]= acos((ID_split[[j]][i,4])/(ID_split[[j]][i,6]))
      ID_split[[j]][,7][is.na(ID_split[[j]][,7])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,11]=((ID_split[[j]][i,6])/TimeInterval)^2                                               # creating values for Square Speed
      return(ID_split[[j]][i,c(4:7, 11)])
    }))

    ID_split[[j]][1:MM,c(4:7, 11)] <- as.data.frame(res)
    ID_split[[j]][,c(4:7, 11)] <- lapply(ID_split[[j]][,c(4:7, 11)], as.numeric)
  }

  for(j in 1:length(ID_split)){
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res1 <- t(sapply(1:MM, function(i){
      ID_split[[j]][,12]=cumsum(ID_split[[j]][,6])                                                            # creating values for cumsum
      ID_split[[j]][i,13]= sqrt(((ID_split[[j]][i+1,2])^2)+((ID_split[[j]][i+1,3])^2))/(ID_split[[j]][i,12])  # creating values for cumulative directionality ratio
      return(ID_split[[j]][i,12:13])
    }))
    ID_split[[j]][1:MM,12:13] <- as.data.frame(res1)
    ID_split[[j]][,12:13] <- lapply(ID_split[[j]][,12:13], as.numeric)
  }

  for(j in 1:length(ID_split)){              # creating values for  rel.ang.P  (step to the previous)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM1<-MM-1
    res <- sapply(1:MM1, function(i){

      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]>=0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]<0) ){
        ID_split[[j]][i,8]= abs(ID_split[[j]][i+1,7])+abs(ID_split[[j]][i,7])
      }
      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]<0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]>=0) ){
        ID_split[[j]][i,8]=ID_split[[j]][i+1,7]-ID_split[[j]][i,7]
      }
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])<= (-pi), 2*pi+(ID_split[[j]][i,8]),(ID_split[[j]][i,8]))    # adjusting the rel.ang
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])> pi,(ID_split[[j]][i,8])-2*pi,(ID_split[[j]][i,8]))
      return(ID_split[[j]][i, 8])
    })
    ID_split[[j]][1:MM1, 8] <- as.data.frame(res)
  }

  cosine.P<-data.frame()
  for(j in 1:length(ID_split)){              # creating values for  cosine.P  based on rel.ang.P
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      ID_split[[j]][i,9]<-cos(ID_split[[j]][i,8])
      return(ID_split[[j]][i,9])
    })
    ID_split[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-ID_split[[j]][,9]
  }


  for(j in 1:length(ID_split)){              # Computing persistence time   (based on rel.ang.P)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      if(abs(ID_split[[j]][i,8])<=1.5707963268){
        ID_split[[j]][i,10]= TimeInterval
      }
      if(abs(ID_split[[j]][i,8])>1.5707963267){
        ID_split[[j]][i,10]= 0
      }
      return(ID_split[[j]][i,10])
    })
    ID_split[[j]][1:MM, 10] <- as.data.frame(res)
  }


  for(j in 1:length(ID_split)){              # Computing Acceleration
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM2<-MM-2
    res <- sapply(1:MM2, function(i){
      ID_split[[j]][i,24]= (sqrt(ID_split[[j]][i+1,11])-sqrt(ID_split[[j]][i,11]))/T
      return(ID_split[[j]][i,24])
    })
    ID_split[[j]][1:MM2, 24] <- as.data.frame(res)
  }

  size<-c()
  for(j in 1:length(ID_split)){
    size[j]<-length(ID_split[[j]][,1])
  }
  S<-summary(size)
  names(S)<-NULL
  IncompleteTracks<-subset(size,S<S[6])
  for(j in 1:length(ID_split)){
    ID_split[[j]]<-ID_split[[j]][1:S[1],]
    ID_split[[j]][S[1],4:25]=0
    ID_split[[j]][S[1],10]=TimeInterval
    ID_split[[j]][1,10]=0
  }

  cat("The minimum number of steps: ",S[1],"\n")
  cat("The maximum number of steps: ",S[6],"\n")
  cat("Number of cells with a total number of steps less than ",S[6],"steps",":",length(IncompleteTracks),"\n")
  cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
  PreprocessedData<-ID_split
  object@preprocessedDS<-PreprocessedData
  return(object)
}



#' @title Data preprocessing for wound scratch assay (WSA).

#'
#' @description This function allows filtering of cells and preprocessing of the trajectory data from wound scratch assay (WSA) experiments.
#' @param object \code{CellMig} class object.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param PixelSize A numeric value of the physical size of a pixel.
#' @param imageH A numeric value of the image hight.
#' @param woundH A numeric value of the image hight.
#' @param upperE A numeric value of the upper edge of the wound.
#' @param lowerE A numeric value of the lower edge of the wound.
#' @param mar A numeric value of the margin to be used to narrow the clearing zone inside the zone.
#' @param clearW A logical vector that allows removing the cells within the wound. Default is TRUE.
#'
#'
#' @return An CellMig class object with filtered, annotated and preprocessed data.
#' @examples
#'
#' data(WSAdataset)
#' wsaTD <- CellMig(WSAdataset)
#' wsaTD <- wsaPreProcessing(wsaTD)
#' @export
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
wsaPreProcessing = function(object,PixelSize=1.24,TimeInterval=10,imageH=1500,woundH=600,upperE=400,lowerE=1000,mar=75,clearW=TRUE) {
  msg <- NULL
  if ( ! is.data.frame(object@trajdata)){
    msg <- c(msg, "input data must be data.frame")
  }
  if ( ! is.numeric(PixelSize)) stop( "PixelSize has to be a positive number" ) else if ( PixelSize<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(imageH)) stop( "imageH has to be a positive number" ) else if ( imageH<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(woundH)) stop( "woundH has to be a positive number" ) else if ( woundH<= 0 ) stop( "woundH has to be a positive number" )
  if ( ! is.numeric(upperE) ) stop( "upperE has to be a positive number" ) else if ( upperE<= 0 ) stop( "upperE has to be a positive number" )
  if ( ! is.numeric(lowerE) ) stop( "lowerE has to be a positive number" ) else if ( lowerE<= 0 ) stop( "PixelSize has to be a positive number" )
  if ( ! is.numeric(mar)) stop( "mar has to be a positive number" ) else if ( mar<= 0 ) stop( "mar has to be a positive number" )
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  if (clearW == TRUE || clearW == T){
    splitFORu<-split(object@trajdata,object@trajdata[,1])
    for (i in 1:length(splitFORu)){
      if ((splitFORu[[i]][1,3]>= (upperE + mar) & splitFORu[[i]][1,3]<= (lowerE -mar)) & (splitFORu[[i]][1,4]<=20 )){   # to remove cells within the wound
        splitFORu[[i]]<-NA
      }
    }
    deletedCells<-splitFORu[is.na(splitFORu)]   ########## To get the deleted cells
    cat(paste0(length(deletedCells)," cells were inside the wound and they have been removed"),"\n")
    Ftable<-splitFORu[!is.na(splitFORu)]
    Ftable<-do.call(rbind.data.frame, Ftable)    ########## to convert the list to a data frame
    rownames(Ftable)<-NULL
    object@adjDS <- Ftable
  }else{
    object@adjDS <- object@trajdata
  }

  ####### to set the cells orientation

  hh<- upperE + ((lowerE-upperE)/2)
  CellOr<-c()
  finaltable<- object@adjDS
  finaltable<-split(finaltable,finaltable[,1])
  for (i in 1:length(finaltable)){
    if (finaltable[[i]][1,3]<=hh){
      CellOr[i]=0
    }else{
      CellOr[i]=1
    }
  }
  object@cellpos<-CellOr

  df<-object@adjDS
  df<-df[,1:3]                                        # Removing the unnecessary columns
  L<-length(df[,1])
  df[,4:26]<-rep(0,L)
  df[,27]<-rep(NA,L)                                  # to be used for migration type
  colnames(df)<-c("ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P","Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY","New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F","Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type")
  ID_split <- split(df, df$ID)                        #Splitting the data frame based on the ID
  cat("This dataset contains: ",length(ID_split),"Cells","\n")

  for(j in 1:length(ID_split)){                        # Having the ID =group order
    ID_split[[j]][1]=j
  }

  for(j in 1:length(ID_split)){                        # adjusting x and y (starting from 0 & being multiplied by H)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(2:MM, function(i){
      ID_split[[j]][i,4]=PixelSize*((ID_split[[j]][i,2])-( ID_split[[j]][1,2]))      # x2-x1
      ID_split[[j]][i,5]=PixelSize*(( ID_split[[j]][1,3])-(ID_split[[j]][i,3]))      # y2-y1
      return(ID_split[[j]][i,4:5])
    }))
    ID_split[[j]][2:MM,4:5] <- as.data.frame(res)
    ID_split[[j]][,4:5] <- lapply(ID_split[[j]][,4:5], as.numeric)
  }



  for(j in 1:length(ID_split)){                    # removing the old x and y
    ID_split[[j]]=ID_split[[j]][-2]            # removing x column
    ID_split[[j]]=ID_split[[j]][-2]            # removing the y column [-2] is used because x column is gone.
  }



  for(j in 1:length(ID_split)){                    # creating values for dx, dy, dis, abs.ang,cumsum, Dir.R
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- t(sapply(1:MM, function(i){
      ID_split[[j]][i,4]=(ID_split[[j]][i+1,2])-( ID_split[[j]][i,2])                                         # creating values for dx
      ID_split[[j]][,4][is.na(ID_split[[j]][,4])] <- 0
      ID_split[[j]][i,5]= ( ID_split[[j]][i+1,3])-(ID_split[[j]][i,3])                                        # creating values for dy
      ID_split[[j]][,5][is.na(ID_split[[j]][,5])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,6]=sqrt((ID_split[[j]][i,4])^2 + (ID_split[[j]][i,5])^2)                                # creating values for dis
      ID_split[[j]][i,7]= acos((ID_split[[j]][i,4])/(ID_split[[j]][i,6]))
      ID_split[[j]][,7][is.na(ID_split[[j]][,7])] <- 0                                                        # to remove NA and replace it with 0
      ID_split[[j]][i,11]=((ID_split[[j]][i,6])/TimeInterval)^2                                               # creating values for Square Speed
      return(ID_split[[j]][i,c(4:7, 11)])
    }))

    ID_split[[j]][1:MM,c(4:7, 11)] <- as.data.frame(res)
    ID_split[[j]][,c(4:7, 11)] <- lapply(ID_split[[j]][,c(4:7, 11)], as.numeric)
  }

  for(j in 1:length(ID_split)){
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res1 <- t(sapply(1:MM, function(i){
      ID_split[[j]][,12]=cumsum(ID_split[[j]][,6])                                                            # creating values for cumsum
      ID_split[[j]][i,13]= sqrt(((ID_split[[j]][i+1,2])^2)+((ID_split[[j]][i+1,3])^2))/(ID_split[[j]][i,12])  # creating values for cumulative directionality ratio
      return(ID_split[[j]][i,12:13])
    }))
    ID_split[[j]][1:MM,12:13] <- as.data.frame(res1)
    ID_split[[j]][,12:13] <- lapply(ID_split[[j]][,12:13], as.numeric)
  }

  for(j in 1:length(ID_split)){              # creating values for  rel.ang.P  (step to the previous)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM1<-MM-1
    res <- sapply(1:MM1, function(i){

      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]>=0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]<0) ){
        ID_split[[j]][i,8]= abs(ID_split[[j]][i+1,7])+abs(ID_split[[j]][i,7])
      }
      if((ID_split[[j]][i+1,5]<0) && (ID_split[[j]][i,5]<0)||(ID_split[[j]][i+1,5]>=0) && (ID_split[[j]][i,5]>=0) ){
        ID_split[[j]][i,8]=ID_split[[j]][i+1,7]-ID_split[[j]][i,7]
      }
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])<= (-pi), 2*pi+(ID_split[[j]][i,8]),(ID_split[[j]][i,8]))    # adjusting the rel.ang
      ID_split[[j]][i,8]<-ifelse((ID_split[[j]][i,8])> pi,(ID_split[[j]][i,8])-2*pi,(ID_split[[j]][i,8]))
      return(ID_split[[j]][i, 8])
    })
    ID_split[[j]][1:MM1, 8] <- as.data.frame(res)
  }

  cosine.P<-data.frame()
  for(j in 1:length(ID_split)){              # creating values for  cosine.P  based on rel.ang.P
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      ID_split[[j]][i,9]<-cos(ID_split[[j]][i,8])
      return(ID_split[[j]][i,9])
    })
    ID_split[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-ID_split[[j]][,9]
  }


  for(j in 1:length(ID_split)){              # Computing persistence time   (based on rel.ang.P)
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      if(abs(ID_split[[j]][i,8])<=1.5707963268){
        ID_split[[j]][i,10]= TimeInterval
      }
      if(abs(ID_split[[j]][i,8])>1.5707963267){
        ID_split[[j]][i,10]= 0
      }
      return(ID_split[[j]][i,10])
    })
    ID_split[[j]][1:MM, 10] <- as.data.frame(res)
  }


  for(j in 1:length(ID_split)){              # Computing Acceleration
    M<- ID_split[[j]][1]
    MM<-length(M[,1])
    MM2<-MM-2
    res <- sapply(1:MM2, function(i){
      ID_split[[j]][i,24]= (sqrt(ID_split[[j]][i+1,11])-sqrt(ID_split[[j]][i,11]))/T
      return(ID_split[[j]][i,24])
    })
    ID_split[[j]][1:MM2, 24] <- as.data.frame(res)
  }

  size<-c()
  for(j in 1:length(ID_split)){
    size[j]<-length(ID_split[[j]][,1])
  }
  S<-summary(size)
  names(S)<-NULL
  IncompleteTracks<-subset(size,S<S[6])
  for(j in 1:length(ID_split)){
    ID_split[[j]]<-ID_split[[j]][1:S[1],]
    ID_split[[j]][S[1],4:25]=0
    ID_split[[j]][S[1],10]=TimeInterval
    ID_split[[j]][1,10]=0
  }

  cat("The minimum number of steps: ",S[1],"\n")
  cat("The maximum number of steps: ",S[6],"\n")
  cat("Number of cells with a total number of steps less than ",S[6],"steps",":",length(IncompleteTracks),"\n")
  cat("All the tracks are adjusted to have only ",S[1]," steps","\n")
  PreprocessedData<-ID_split
  object@preprocessedDS<-PreprocessedData
  return(object)
}




#' Method Plotting
#' @title A 2D rose-plot
#' @description Plotting the trajectory data of all cells.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' data(Trajectory_dataset)
#' rmTD <- CellMig(Trajectory_dataset)
#' \dontrun{
#' rmTD <- rmPreProcessing(rmTD)
#' plotAllTracks(rmTD, ExpName="Test",Type="b")
#' }
#'
#'
plotAllTracks= function(object,ExpName="ExpName",Type="l") {
  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }

  Len<-length(Object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
  Step<-length(Object[[1]][,1])
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
    minX=min(Object[[j]][1:Step,2])
    minY=min(Object[[j]][1:Step,3])
    maxX=max(Object[[j]][1:Step,2])
    maxY=max(Object[[j]][1:Step,3])
    MinX[j]<-c(minX)
    MaxX[j]<-c(maxX)
    MinY[j]<-c(minY)
    MaxY[j]<-c(maxY)
  }
  RangeX=c(MinX,MaxX)
  RangeY=c(MinY,MaxY)
  plot(Object[[1]][1:Step,2],Object[[1]][1:Step,3],type=Type,xlab="X (um)",ylab="Y (um)",col=color[1],las=1,xlim=range(RangeX),ylim=range(RangeY),main=ExpName)
  for(n in 2:Len){
    points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
    points(end,pch=16,col=color[n], cex = 1)
  }
  x=c(min(RangeX)-100,max(RangeX)+100)
  y=c(0,0)
  lines(x, y, type='l', col="black")
  x=c(0,0)
  y=c(min(RangeY)-100,max(RangeY)+100)
  lines(x, y, type='l', col="black")

  jpeg(paste0(ExpName,"_All_tracks_plot.jpg"))
  plot(Object[[1]][1:Step,2],Object[[1]][1:Step,3],type=Type,xlab="X (um)",ylab="Y (um)",col=color[1],las=1,xlim=range(RangeX),ylim=range(RangeY),main=ExpName)
  for(n in 2:Len){
    points(Object[[n]][,2],Object[[n]][,3], type=Type,col=color[n])
    end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
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
#' @title A 3D rose-plot of all cells
#' @description Plotting the trajectory data of all cells in 3D.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' plot3DAllTracks(rmTD, VS=3, size=2)
#' }
#'
plot3DAllTracks= function(object,VS=3,size=2) {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }

  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )
  plotTable<-data.frame()
  Len<-length(Object)
  cat(paste0("The plot contains ",Len, " Cells"),"\n")
  Step<-length(Object[[1]][,1])
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
    plotTable[RowNum,1]= Object[[i]][1:Step,2]
    plotTable[RowNum,2]=Object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
  }
  plot3d(plotTable, col=coll, type="p", size=size, axes=F,xlab=" ", ylab=" ",zlab=" ")
}

#'
#'
#' @title A 3D rose-plot
#' @description Plotting the trajectory data of particular cells in 3D.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' plot3DTracks(rmTD, VS=3, size=2,cells=c(1,50,150,250,350))
#' }
#'
plot3DTracks= function(object,VS=3,size=2,cells) {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(VS) ) stop( "VS has to be a positive number" ) else if ( VS<= 0 ) stop( "VS has to be a positive number" )
  if ( ! is.numeric(size) ) stop( "size has to be a positive number" ) else if ( size<= 0 ) stop( "size has to be a positive number" )
  plotTable<-data.frame()
  Len<-length(Object)
  cat(paste0("The plot contains ",length(cells), " Cells"),"\n")
  Step<-length(Object[[1]][,1])
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
    plotTable[RowNum,1]= Object[[i]][1:Step,2]
    plotTable[RowNum,2]=Object[[i]][1:Step,3]
    plotTable[RowNum,3]=i*VS
    col= c(rep(color[i],Step))
    coll=c(coll,col)
    NewplotTable<-plotTable[complete.cases(plotTable),]

  }
  plot3d(NewplotTable, col=coll, type="p", size=size, axes=F,xlab=" ", ylab=" ",zlab=" ")
}
#'
#'
#' @title A graphical display of the track of each cell.
#' @description Plotting the trajectory data of each cell.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' PlotTracksSeparately(rmTD, ExpName="Test",Type="b",FixedField=FALSE)
#' }
PlotTracksSeparately= function(object,ExpName="ExpName",Type="l",FixedField=TRUE) {
  if ( ! ( Type %in% c("p","l","b","o") ) ) stop("Type has to be one of the following: p, l, b, o")
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  Len<-length(Object)
  cat(paste0(Len," plots will be generated in a folder called:",ExpName,"_Tracks","\n"))
  Step<-length(Object[[1]][,1])
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
      minX=min(Object[[j]][1:Step,2])
      minY=min(Object[[j]][1:Step,3])
      maxX=max(Object[[j]][1:Step,2])
      maxY=max(Object[[j]][1:Step,3])
      MinX[j]<-c(minX)
      MaxX[j]<-c(maxX)
      MinY[j]<-c(minY)
      MaxY[j]<-c(maxY)
    }
    RangeX=c(MinX,MaxX)
    RangeY=c(MinY,MaxY)
    for(n in 1:Len){
      jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"))
      plot(Object[[n]][1:Step,2],Object[[n]][1:Step,3],type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,xlim=range(RangeX),ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      lines(x, y, type='l', col="black")
      end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
      points(end,pch=16,col=color[n], cex = 1)
      title(main=paste0("Cell Number  ", n),col.main="black")
      dev.off()
    }
  }else{
    for(n in 1:Len){
      RangeX= Object[[n]][1:Step,2]
      RangeY= Object[[n]][1:Step,3]
      jpeg(paste0(ExpName,"_Track_Plot_",n,".jpg"))
      plot(Object[[n]][1:Step,2],Object[[n]][1:Step,3],type=Type,xlab="x (um)",ylab="y (um)",col=color[n],las=1,xlim=range(RangeX),ylim=range(RangeY))
      x=c(min(RangeX)-100,max(RangeX)+100)
      y=c(0,0)
      lines(x, y, type='l', col="black")
      x=c(0,0)
      y=c(min(RangeY)-100,max(RangeY)+100)
      lines(x, y, type='l', col="black")
      end<-cbind(Object[[n]][Step,2],Object[[n]][Step,3])
      points(end,pch=16,col=color[n], cex = 1)
      dev.off()
    }

  }
  setwd(d)
}
#'
#'
#'
#'
#'
#' Method Persistence and Speed
#' @title Persistence and Speed
#' @description The PerAndSpeed() generates data and plots for persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param PtSplot A logical vector that allows generating individual plots of persistence time vs speed per cell. Default is TRUE.
#' @param AllPtSplot  A logical vector that allows generating a plot of persistence time vs speed for all cells. Default is TRUE.
#' @param ApSplot A logical vector that allows generating individual plots of angular persistence vs speed per cell. Default is TRUE.
#' @param AllApSplot A logical vector that allows generating a plot of angular persistence vs speed of all cells. Default is TRUE.
#'
#' @return An CellMig class object with a data frame and plots. The data frame is stored in the PerAanSpeedtable slot.
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD <- PerAndSpeed(rmTD,TimeInterval=10,ExpName="ExpName")
#' }
#'
PerAndSpeed= function(object,TimeInterval=10,ExpName="ExpName",PtSplot=TRUE,AllPtSplot=TRUE,ApSplot=TRUE,AllApSplot=TRUE) {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )

  d=getwd()
  dir.create(paste0(ExpName,"-PerResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-PerResults")))

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
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
  for(j in 1:length(Object)){                               # creating values (NA and persistence time )for  persistence    (step to the original)
    MM<-length(Object[[j]][,1])
    MM2<-MM-1
    Ptime<-Object[[j]][1:MM2,10]
    MeanPerTime<-round(mean(Ptime),digits=2)             # computing mean persistence time
    PerTimLen<-Ptime                                     # computing the number of persistence steps to be used in computing the persistence ratio
    PerTimLen[PerTimLen==0]<-NA
    PerTimLen<-PerTimLen[!is.na(PerTimLen)]
    PerTimLen<-length(PerTimLen)
    PerRatio<-round(PerTimLen/MM2, digits=2)              # computing persistence ratio
    DD.Ptime<-Object[[j]][1:MM,10]                        # computing the direction deviating time
    DD.Ptime[DD.Ptime==0]<-1                              # replacing the 0 with 5 and the 5 with 0
    DD.Ptime[DD.Ptime==TimeInterval]<-0
    DD.Ptime[DD.Ptime==1]<-TimeInterval
    MeanDD.PerTime<-round(mean(DD.Ptime),digits=2)

    PerResultsTable[1,j]<-j
    PerResultsTable[2,j]<-MeanPerTime
    PerResultsTable[3,j]<-MeanDD.PerTime
    PerResultsTable[4,j]<-PerRatio
  }

  VelPerTable<-data.frame()               # creating a table to store the mean velocity with correspondence with the persistence time
  for(j in 1:length(Object)){
    MM=length(Object[[j]][,1])
    Ptime<-Object[[j]][1:MM,10]
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

    Vel<-Object[[j]][1:MM,11]
    Per<-Object[[j]][1:MM,10]
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


    Ptime<-Object[[j]][1:MM,10]           # computing the "meanVel.for0"
    Ptime00<-c(1,Ptime)                   # adding a "1" value in the beginning  this value should not be 0 because we will not catch 0 persistence if it is in the beginning.
    Ptime00[Ptime00==0]<-NA               # replacing the "0" values with NAs
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

    Vel<-Object[[j]][1:MM,11]
    Per<-Object[[j]][1:MM,10]
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
      title(main=paste0("Cell Number  ", j,"   Speed vs Persistence Time"),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccPV),col.sub="red")
      dev.off()

    }

  }

  ## All cells (Persistence times  vs Speed)
  allper<-VelPerTable[,1]
  allvel<-VelPerTable[,2]
  for(j in 1:length(Object)){
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
  PerResultsTable[6,(length(Object)+1)]<-ccP                                          # Speed vs  persistence time  Spearman correlation

  if ( AllPtSplot == TRUE || AllPtSplot == T){
    jpeg(paste0(ExpName,"_Persist_Time_vs_Speed-All_Cells.jpg"))
    plot(allvel,allper,type="p",pch=16,ylab="Persist_Time (min)",xlab=" Mean Speed during persistence time (um/h)",col="black",las=1)
    abline(reg,untf=F,col="red")
    title("Speed vs Persist Time (All cells)",cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccP),col.sub="red")
    dev.off()

  }

  for(j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    MM<-length(Object[[j]][,1])
    MM2<-MM-1
    Root.Median.Square.Speed<-round(sqrt(median(Object[[j]][1:MM2,11])),digits = 3)*60
    PerResultsTable[7,j]<-Root.Median.Square.Speed

    wma<-which.max(sqrt(Object[[j]][,11]))
    wmi<-which.min(sqrt(Object[[j]][1:MM2,11]))

    PerResultsTable[8,j]<-round(sqrt(Object[[j]][wma,11]),digits=3)* 60
    PerResultsTable[9,j]<-round(sqrt(Object[[j]][wmi,11]),digits=3)* 60

    mean.cosineP<-round(mean(Object[[j]][1:(MM2-1),9],na.rm = TRUE),digits = 3)
    PerResultsTable[10,j]<-mean.cosineP
    s<-cor.test( ~ sqrt(Object[[j]][1:MM2,11])+ Object[[j]][1:MM2,9], method = "spearman",exact=FALSE)                 #testing the correlation
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    PerResultsTable[11,j]<-VEvsCOSP

    data<-Object[[j]][1:MM2,2:3]
    data1<-Object[[j]][1:round(MM2/4),2:3]
    data2<-Object[[j]][round(MM2/4):round(MM2/2),2:3]
    data3<-Object[[j]][round(MM2/2):round(MM2*3/4),2:3]
    data4<-Object[[j]][round(MM2*3/4):MM2,2:3]

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


    PerResultsTable[15,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3)
    PerResultsTable[16,j]=round(sum(abs(Object[[j]][,8]))/6.28,digits=3) - abs(round(sum(Object[[j]][,8])/6.28,digits=3))


    if ( ApSplot == TRUE || ApSplot == T){
      jpeg(paste0(ExpName,"_Angular_Persistence_vs_Speed",j,".jpg"))
      Speed=(sqrt(Object[[j]][1:MM2,11])/TimeInterval)*60
      plot(Speed,Object[[j]][1:MM2,9],pch=16,type="p",ylab="Angular Persistence (cosine)",xlab=" Instantaneous Speed (um/h)",col=color[j],las=1, xlim=c(0,3))
      reg<-lm(Object[[j]][1:MM2,9]~Speed)
      abline(reg,untf=F,col="black")
      title(main=paste0("Cell Number  ", j," Instantaneous Speeds vs Angular Persistence "),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
      dev.off()
    }

    PoCos<- subset(Object[[j]][1:(MM2-1),9],Object[[j]][1:(MM2-1),9]>0)
    NeCos<- subset(Object[[j]][1:(MM2-1),9],Object[[j]][1:(MM2-1),9]<=0)

    PerResultsTable[17,j]<-round(median(PoCos,na.rm = TRUE),digits = 3)
    PerResultsTable[18,j]<-round(median(NeCos,na.rm = TRUE),digits = 3)

    AvSp<-(Object[[j]][1:length(Object[[j]][,1])-1,6]/TimeInterval)*60

    PerResultsTable[19,j]<-round(median(AvSp,na.rm = TRUE),digits = 3)
    s=summary(AvSp)
    PerResultsTable[20,j]<-round((s[5]-s[2])/s[3],digits=4)
    PerResultsTable[21,j]<-round(mean(AvSp),digits = 3)
    PerResultsTable[22,j]<-round(sd(AvSp),digits = 3)
  }

  #(all cells)  persistence vs inst.speed
  MM<-length(Object[[1]][,1])
  MM2<-MM-1
  cosine.P<-data.frame()
  for(j in 1:length(Object)){              # creating values for  cosine.P  based on rel.ang.P
    M<- Object[[j]][1]
    MM<-length(M[,1])
    res <- sapply(1:MM, function(i){
      Object[[j]][i,9]<-cos(Object[[j]][i,8])
      return(Object[[j]][i,9])
    })
    Object[[j]][1:MM, 9] <- as.data.frame(res)
    cosine.P[1:MM,j]<-Object[[j]][,9]
  }

  RM<-round(rowMedians(as.matrix(cosine.P[1:MM2,]),na.rm = TRUE),digits=3)
  Speed<-data.frame()
  for (j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    Speed[1:MM2,j]<-round(sqrt(Object[[j]][1:MM2,11]),digits = 3)
  }
  RowmeanSpeed<-round(rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  VEvsCOSP<-round(ss, digits = 3)
  PerResultsTable[11,(length(Object)+1)]<-VEvsCOSP

  if ( AllApSplot == TRUE || AllApSplot == T){
    jpeg(paste0(ExpName," All_Cells_Average_Angular_Persistence_vs_Average_Speed.jpg"))
    MS<-max(RowmeanSpeed)*60
    plot(RowmeanSpeed*60,RM,pch=16,type="p",ylab="Average Instantaneous Angular Persistence (cosine)",xlab=" Average Instantaneous Speed (um/h)",col="black",las=1,xlim=c(0,MS))
    NewSpeed=RowmeanSpeed*60
    reg<-lm(RM~NewSpeed)
    abline(reg,untf=F,col="red")
    title(main=paste0("All Cells Instantaneous Speed vs Angular Persistence "),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
    dev.off()
  }

  rownames(PerResultsTable)<-c("Cell Number","Mean Persist Time (min)","Mean persist Deviating Time (min)","Persistence Ratio",
                               "Maximum Persistence period (min)","Persistence Time vs Speed (SCC)","RMSS (um per h)","Maximum Speed (um per h)","Minimum Speed (um per h)",
                               "Mean Angular Persistence (cosine)","Instantaneous Speed vs Angular Persistence (SCC)","Covered Area (um2)","Segmented Covered Area (um2)","Empty Area (um2)","Number of complete rotations",
                               "Number of canceled rotations","Mean Persistence Angle (cosine)","Mean Deviating Angle (cosine)","Median Speed","Speed QBCV",
                               "Mean Speed (um per h)","Speed standard deviation (um per h)")



  RM1<-round(rowMedians(as.matrix(PerResultsTable),na.rm = TRUE),digits=3)
  PerResultsTable[c(2:5,7:10,12:22),(length(Object)+1)]<-RM1[c(2:5,7:10,12:22)]

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
  s<-cor.test( ~ SPEED+ PerR, method = "spearman",exact=FALSE)                 # testing the correlation
  ss<-unlist(s[4])
  SCC<-round(ss, digits = 3)
  title(main=paste0("All Cells Average Speed vs Persistence Ratio "),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",SCC),col.sub="red")

  PerResultsTable[1,(length(Object)+1)]<-"All Cells"
  object@PerAanSpeedtable <-PerResultsTable
  setwd(d)
  write.csv(PerResultsTable, file = paste0(ExpName,"-PerResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-PerResultsTable.csv" ),"in your directory [use getwd()]","\n")
  return(object)
}
#'
#'
#' Method DR
#' @title Directionality Table
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straightline length between the start point and the endpoint of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#'
#' @return An CellMig class object with a data frame stored in the DRtable slot. It contains nine rows: "Cell Number","Directionality Ratio","Mean Cumulative Directionality Ratio","Stable Directionality Ratio", "Number of returns","Min CumDR","Location of Min CumDR, Steps with less CumDR than DR","Directional Persistence".
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD <- DiRatio(rmTD, ExpName="Test")
#' }
DiRatio = function(object,TimeInterval=10,ExpName="ExpName") {
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  d=getwd()
  dir.create(paste0(ExpName,"-DRResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-DRResults")))
  Step<-length(Object[[1]][,1])
  DRResultsTable<-data.frame()
  DIR.RATIO<-c()
  for(j in 1:length(Object)){                             # calculating the cumsum of distance for each cell
    MM<-Step
    MM2<-MM-1
    end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])       # finding the cordinates of the final point in the track.
    start<-cbind(0,0)
    final.dis=dist2(start, end)                               # calculating the final distance
    total.dis=sum(Object[[j]][1:MM2,6])                     # calculating the total distance
    Dir.Ratio=round(final.dis/total.dis ,digits = 3)          # calculating the Dir.Ratio
    mean.Dir.Ratio<-round(mean(Object[[j]][1:MM2,13],na.rm = TRUE) ,digits = 3)
    StableDR<- (1-(mean.Dir.Ratio-Dir.Ratio))* mean.Dir.Ratio
    StableDR<-round(StableDR,digits = 3)

    DRResultsTable[1,j]<- j
    DRResultsTable[2,j]<-Dir.Ratio
    DRResultsTable[3,j]<-mean.Dir.Ratio
    DRResultsTable[4,j]<- StableDR

  }


  for(j in 1:length(Object)){                                   #### Adding min CumDR  and number of angles greater than .75
    MM<-Step
    MM2<-MM-1
    p1<-Object[[j]][1:MM2,9]
    returns<-subset(p1,p1<(-0.87))                            # greater than 150 degrees
    DRResultsTable[5,j]<-length(returns)

    p2<-Object[[j]][1:MM2,13]
    w<-which.min(p2)
    DRResultsTable[6,j]<-round(p2[w], digits=3)
    DRResultsTable[7,j]<-w
    DR<-as.numeric(DRResultsTable[2,j])
    lessThanMINcumdr<-subset(p2,p2<DR)                        # number of the steps that have a cumdr less than the final dr
    DRResultsTable[8,j]<-length(lessThanMINcumdr)

    Ptime<-Object[[j]][1:MM2,10]
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
  DRResultsTable[,(length(Object)+1)]<-RM1
  DRResultsTable[1,(length(Object)+1)]<-"All Cells"
  object@DRtable<-DRResultsTable
  setwd(d)
  write.csv(DRResultsTable, file = paste0(ExpName,"-DRResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-DRResultsTable.xlsx" ),"in your directory [use getwd()]","\n")
  return(object)
}
#'
#'
#' Method DR
#' @title Directionality Ratio plots
#' @description Directionality Ratio is the displacement divided by the total length of the total path distance, where displacement is the straightline length between the start point and the endpoint of the migration trajectory,
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD <-DiRatio(rmTD, ExpName="Test")
#' DiRatio.Plot(rmTD, ExpName="Test")
#' }
DiRatio.Plot = function(object,TimeInterval=10,ExpName=ExpName) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run PreProcessing()")
  }
  Step<-length(Object[[1]][,1])
  color <-c()
  Len<-length(Object)
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
  for(j in 1:length(Object)){
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
    p<-plot(Time,Object[[j]][1:MM2,13], type="l",col=color[j],xlab="Time (hours)",xaxt="n",ylab="Directionality Ratio",lwd=2,las=1)
    axis(1, at=xMMM1, cex.axis=0.8,labels=xMMM)
    title(main=paste0("Cell Number  ", j),col.main=color[j])
    dev.off()
    DIR.RATIO.AllCells[1:MM2,j]<-Object[[j]][1:MM2,13]
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
#'
#'
#'
#'
#'
#'
#'
#'
#' Method msd
#' @title Mean Square Displacement
#' @description The MSD function automatically compute the mean square displacements across several sequantial time intervals. MSD parameters are used to assess the area explored by cells over time.
#'
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param ffLAG A numeric value to be used to get the number of lags for the  Furth formula fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param SlopePlot A logical vector that allows generating individual plots showing the slope of the mean square displacement of the movement of individual cells. Default is TRUE.
#' @param AllSlopesPlot A logical vector that allows generating a plot showing the slope of the mean square displacement of the movement of all cells. Default is TRUE.
#' @param FurthPlot A logical vector that allows generating individual plots fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method per cell. Default is TRUE.
#' @param AllFurthPlot A logical vector that allows generating a plot fitting the Furth formula using generalized regression by the Nelder–Mead method simplex method for all cells. Default is TRUE.
#'
#' @return An CellMig class object with a data frame and plots. The data frame is stored in the MSDtable slot.
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD<-MSD(rmTD,sLAG=0.25, ffLAG=0.25)
#' }

MSD = function(object,TimeInterval=10,ExpName="ExpName",sLAG=0.25, ffLAG=0.25, SlopePlot=TRUE,AllSlopesPlot=TRUE,FurthPlot=TRUE,AllFurthPlot=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  if ( ! is.numeric(ffLAG) ) stop( "ffLAG has to be a positive number" ) else if ( ffLAG<= 0 ) stop( "ffLAG has to be a positive number" )
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-MSDResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-MSDResults")))
  Len<-length(Object)
  Step<-length(Object[[1]][,1])
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
  for(j in 1:length(Object)){
    meanSD<-c()
    LAG<-round(Step*sLAG)                                                     # number of lags is based on sLAG
    for(lag in 1:LAG){
      res <- sapply(1:Step, function(i){
        Object[[j]][i,22]=(((Object[[j]][i+lag,2]- Object[[j]][i,2])^2)+((Object[[j]][i+lag,3]- Object[[j]][i,3])^2))
        return(Object[[j]][i,22])
      })
      Object[[j]][1:Step, 22] <- as.data.frame(res)
      Object[[j]][,22][is.na(Object[[j]][,22])] <- 0
      meanSD[lag]<-mean(Object[[j]][1:(Step-LAG),22])
      Object[[j]][,22]=0
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
    Data <- matrix(ncol =  2, byrow = F, data =Data)
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
  MSDResultsTable[1,(length(Object)+1)]<-"All Cells"
  MSDResultsTable[2,(length(Object)+1)]<-round(RM1[1],digits=3)

  RM<-c(0,RM1)
  Xaxis<-c(1:round(Step*ffLAG))
  NewrowMeans<-RM1[1:(round(LAG*ffLAG)+1)]                                  # best fit based on ffLAG *ffLAG
  NewXaxis<-Xaxis[1:(round(LAG*ffLAG)+1)]

  reg<-lm(NewrowMeans~ NewXaxis)
  reg1<-round(coef(lm(log10(NewrowMeans)~ log10(NewXaxis)))[2],digits=2)
  MSDResultsTable[3,(length(Object)+1)]<-reg1
  y=RM1[1:round(Step*ffLAG)]
  t <- c(1:length(y))
  Data<-c(t,y)
  Data <- matrix(ncol =  2, byrow = F, data =Data)
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
  MSDResultsTable[4,(length(Object)+1)]<-round(Fit$par[1],digits=3)
  MSDResultsTable[5,(length(Object)+1)]<-round(Fit$par[2],digits=3)
  MSDResultsTable[6,(length(Object)+1)]<-round(Summary$par[1,4],digits=3)
  MSDResultsTable[7,(length(Object)+1)]<-round(Summary$par[2,4],digits=3)

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
  object@MSDtable<-MSDResultsTable
  setwd(d)
  write.csv(MSDResultsTable, file = paste0(ExpName,"-MSDResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(object)
}
#'
#'
#'
#' Method DiAutoCor
#' @title Direction AutoCorrelation
#'
#' @description The DiAutoCor function automatically compute the angular persistence across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param sLAG A numeric value to be used to get the number of lags for the slope fitting. Default is 0.25, which represents 25 percent of the steps.
#' @param sPLOT A logical vector that allows generating individual plots showing the angular persistence across several sequantial time intervals. Default is TRUE.
#' @param aPLOT A logical vector that allows generating a plot showing the angular persistence across several sequantial time intervals of all cells. Default is TRUE.
#' @return An CellMig class Object with a data frame and plots. The data frame, which contains six rows: "Cell Number", "Angular Persistence", "Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)", "Stable Direction AutoCorrelation through the track" and "Difference between Mean DA and Intercept DA".
#'
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD <- DiAutoCor(rmTD,TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE)
#' }

DiAutoCor= function(object, TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  Object<-object@preprocessedDS

  if ( ! is.list(Object) ){
    stop("Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-DIAutoCorResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-DIAutoCorResults")))

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
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

  for(j in 1:length(Object)){
    cos.diff<-c()
    LAG<-round(Step*sLAG)     #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:Step, function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2]    # newdx
        Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3]    # newdy
        return(Object[[j]][i,14:15])
      }))
      Object[[j]][1:Step,14:15] <- as.data.frame(res)


      Object[[j]][,14:15] <- lapply(Object[[j]][,14:15], as.numeric)
      Object[[j]][,14][is.na(Object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      Object[[j]][,15][is.na(Object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:Step, function(i){
        Object[[j]][i,16]=acos((Object[[j]][i+lag,2]-Object[[j]][i,2])/sqrt((Object[[j]][i+lag,2]-Object[[j]][i,2])^2 +(Object[[j]][i+lag,3]-Object[[j]][i,3])^2)) # to find the abs angle
        return(Object[[j]][i,16])
      })
      Object[[j]][1:(Step),16] <- res1
      Object[[j]][,16][is.na(Object[[j]][,16])] <- 0                                # to remove NA and replace it with 0


      res2 <- sapply(1:(Step - lag), function(i){
        if((Object[[j]][i+1,15]<0) && (Object[[j]][i,15]>=0)||(Object[[j]][i+1,15]>=0) && (Object[[j]][i,15]<0)){
          Object[[j]][i,17]= abs(Object[[j]][i+1,16])+abs(Object[[j]][i,16])
        }
        if((Object[[j]][i+1,15]<0) && (Object[[j]][i,15]<0)||(Object[[j]][i+1,15]>=0) && (Object[[j]][i,15]>=0) ){
          Object[[j]][i,17]=Object[[j]][i+1,16]-Object[[j]][i,16]
        }

        return(Object[[j]][i,17])
      })
      Object[[j]][1:(Step - lag),17] <- res2

      res3 <- t(sapply(1:(Step - lag), function(i){
        Object[[j]][i,17]<-ifelse((Object[[j]][i,17])<= (-pi), 2*pi+(Object[[j]][i,17]),(Object[[j]][i,17]))    # adjusting the ang.diff
        Object[[j]][i,17]<-ifelse((Object[[j]][i,17])>= pi,(Object[[j]][i,17])-2*pi,(Object[[j]][i,17]))
        Object[[j]][i,18]<-cos(Object[[j]][i,17])
        return(Object[[j]][i,17:18])
      }))
      Object[[j]][1:(Step - lag),17:18] <- as.data.frame(res3)
      Object[[j]][,17:18] <- lapply(Object[[j]][,17:18], as.numeric)

      for(i in 1:LAG){
        cos.diff[lag]<-mean(Object[[j]][1:(Step-lag)-1,18], na.rm=TRUE)  # computing the cosine mean
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

    Object[[j]][,15:18]=0
  }
  RM1<-rowMedians(as.matrix(DA.table),na.rm = TRUE)
  DA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
  DA.ResultsTable[2,(length(Object)+1)]<-round(RM1[1],digits=3)
  lags<-c(1:length(DA.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  DA.ResultsTable[3,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
  DA.ResultsTable[4,(length(Object)+1)]<-round(median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=3)
  DA.ResultsTable[5,(length(Object)+1)]<-round((1-(median(as.numeric(DA.ResultsTable[4,1:length(Object)]))-unlist(cc[1])))* median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=3)
  DA.ResultsTable[6,(length(Object)+1)]<-round(median(as.numeric(DA.ResultsTable[4,1:length(Object)]))-unlist(cc[1]),digits=3)

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
    MDA<-round(median(as.numeric(DA.ResultsTable[4,1:length(Object)])),digits=2)
    abline(h=MDA,col="blue",lwd = 2)
    title(main=paste0("All Cells - DA quadratic model"),col.main="darkgreen",
          sub=paste0(" Intercept of DA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    legend(1, y=-0.82, legend=c("Mean Direction AutoCorrelation","Quadratic model"), col=c("blue","darkgreen"),lty=1, cex=0.8)
    dev.off()
  }

  rownames(DA.ResultsTable)<-c("Cell Number","Angular Persistence","Intercept of DA quadratic model","Mean Direction AutoCorrelation (all lags)","Stable Direction AutoCorrelation through the track",
                               "Difference between Mean DA and Intercept DA" )
  object@DACtable<-DA.ResultsTable
  setwd(d)
  write.csv(DA.ResultsTable, file = paste0(ExpName,"-DA.ResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(object)
}
#'
#'
#'
#' Method VelAutoCor
#' @title Velocity AutoCorrelation
#'
#' @description The VeAutoCor function automatically compute the changes in both speed and direction across several sequantial time intervals.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
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
#' rmTD <- CellMig(Trajectory_dataset)
#' rmTD <- rmPreProcessing(rmTD)
#' rmTD <-VeAutoCor(rmTD,TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE)
#' }

VeAutoCor= function(object, TimeInterval=10,ExpName="ExpName",sLAG=0.25,sPLOT=TRUE,aPLOT=TRUE) {
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  if ( ! is.numeric(sLAG) ) stop( "sLAG has to be a positive number" ) else if ( sLAG<= 0 ) stop( "sLAG has to be a positive number" )
  Object<-object@preprocessedDS
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  dir.create(paste0(ExpName,"-VeAutoCorResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-VeAutoCorResults")))

  Len<-length(Object)
  Step<-length(Object[[1]][,1])
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

  for(j in 1:length(Object)){
    meanVAC<-c()
    LAG<-round(Step*sLAG)              #taking only the first 12.5%
    for(lag in 1:LAG){
      res <- t(sapply(1:(Step - 1), function(i){                      # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,14]= Object[[j]][i+lag,2]-Object[[j]][i,2]    # newdx
        Object[[j]][i,15]= Object[[j]][i+lag,3]-Object[[j]][i,3]    # newdy
        return(Object[[j]][i,14:15])
      }))
      Object[[j]][1:(Step -1),14:15] <- as.data.frame(res)
      Object[[j]][,14:15] <- lapply(Object[[j]][,14:15], as.numeric)
      Object[[j]][,14][is.na(Object[[j]][,14])] <- 0                                # to remove NA and replace it with 0
      Object[[j]][,15][is.na(Object[[j]][,15])] <- 0                                # to remove NA and replace it with 0

      res1 <- sapply(1:(Step - lag), function(i){                                               # starting from 2 to exclude the first cosine which is always 1.
        Object[[j]][i,23]=((Object[[j]][i,14]* Object[[j]][i+lag,14])+ (Object[[j]][i,15]* Object[[j]][i+lag,15]))/ ((lag*T)^2)
        return(Object[[j]][i,23])
      })
      Object[[j]][1:(Step - lag),23] <- res1
      meanVAC[lag]<-mean(Object[[j]][1:(Step - lag),23])

    }
    VAC.first.value[j]<-meanVAC[1]

    NORMmeanVAC<-meanVAC/meanVAC[1]
    VAC.table[1:LAG,j]<-meanVAC
    assign(paste0("VAC.Cell.",j),meanVAC)
    VAC.second.value[j]<-NORMmeanVAC[2]

    VA.ResultsTable[1,j]<-j
    VA.ResultsTable[2,j]<-round(meanVAC[1],digits=3)        # VA (lag =1)
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

    Object[[j]][,15:16]=0
  }

  RM1<-rowMedians(as.matrix(VAC.table),na.rm = TRUE)
  VA.ResultsTable[1,(length(Object)+1)]<-"All Cells"
  VA.ResultsTable[2,(length(Object)+1)]<-round(median(VAC.first.value),digits=3)
  VA.ResultsTable[3,(length(Object)+1)]<-round(median(VAC.second.value),digits=3)

  lags<-c(1:length(VAC.table[,1]))
  lags2<- lags^2
  quadratic.model<-c()
  quadratic.m <-lm(RM1~ lags + lags2)
  c<-quadratic.m
  cc<-unlist(c)
  quadratic.model[j]<-cc[1]
  VA.ResultsTable[4,(length(Object)+1)]<-round(unlist(cc[1]),digits=3)
  VA.ResultsTable[5,(length(Object)+1)]<-round(median(as.numeric(VA.ResultsTable[5,1:length(Object)])),digits=3)

  ccc<-unlist(cc[1])
  timevalues <- seq(1, length(lags), 1)
  predictedcounts <- predict(quadratic.m,list(Time=timevalues, Time2=timevalues^2))

  if ( aPLOT == TRUE || aPLOT == T){
    Xaxis<-c(1:LAG)
    Yaxis<-RM1
    jpeg(paste0(ExpName,"-Velocity Autocorrelation All Cells.jpg"))
    plot(Xaxis,Yaxis, type="o",ylim=range(RM1),xlim=c(0,lag),col="black",xlab="Lag",ylab="Velocity  Autocorrelation",pch=19,las=1,cex=2)
    lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
    MVA<-round(median(as.numeric(VA.ResultsTable[5,1:length(Object)])),digits=2)
    abline(h=MVA,col="blue",lwd = 2)
    title(main=paste0("All Cells - VA quadratic model"),col.main="darkgreen",
          sub=paste0(" Intercept of VA quadratic model = ",round(ccc, digits=3)),col.sub="red")
    dev.off()
  }

  rownames(VA.ResultsTable)<-c("Cell Number","Velocity AutoCorrelation (lag=1)","2nd normalized Velocity AutoCorrelation","Intercept of VA quadratic model","Mean Velocity AutoCorrelation (all lags)")
  object@VACtable<-VA.ResultsTable
  setwd(d)
  write.csv(VA.ResultsTable, file = paste0(ExpName,"-VA.ResultsTable.csv"))
  cat("Results are saved in your directory [use getwd()]","\n")
  return(object)

}




#'
#'
#'
#'
#' Method Forward Migration
#' @title Forward Migration
#'
#' @description The ForwardMigration function automatically generates data and plots for forward persistence and speed.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @param sfptPLOT A logical vector that allows generating individual plots of persistence time vs speed per cell. Default is TRUE.
#' @param afptPLOT  A logical vector that allows generating a plot of persistence time vs speed for all cells. Default is TRUE.
#' @param sfpPLOT A logical vector that allows generating individual plots of angular persistence vs speed per cell. Default is TRUE.
#' @param afpPLOT A logical vector that allows generating a plot of angular persistence vs speed of all cells. Default is TRUE.
#' @return  An CellMig class Object with a data frame and plots. The data frame is stored in the ForMigtable slot.
#'
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' wsaTD <- CellMig(WSAdataset)
#' wsaTD <- wsaPreProcessing(wsaTD)
#' wsaTD <-ForwardMigration(wsaTD,TimeInterval=10,ExpName="ExpName")
#' }

ForwardMigration= function(object, TimeInterval=10,ExpName="ExpName",sfptPLOT =TRUE,afptPLOT =TRUE,sfpPLOT =TRUE,afpPLOT =TRUE){
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  Object<-object@preprocessedDS
  UPorDO<-object@cellpos
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- rainbow(Len)
  }
  dir.create(paste0(ExpName,"-ForwardMigrationResults"))
  setwd(paste0(d,"/",paste0(ExpName,"-ForwardMigrationResults")))

  #if ( length(UPorDO) != length(Object)) stop("UPorDO needs to be a vector with same number of elements as number of cells")

  for(j in 1:length(Object)){                             # Defining if the cell is ubove (1) the wound or below (0) it
    Object[[j]][,25]<-UPorDO[j]
  }


  for(j in 1:length(Object)){                            # creating values for  rel.ang.F  (step to the original)
    MM<-Step
    MM1<-MM-1
    res <- sapply(1:MM1, function(i){

      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){
        Object[[j]][i,19]= 1.5707963268 - abs(Object[[j]][i,7])
      }
      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
        Object[[j]][i,19]= abs(Object[[j]][i,7])+1.5707963268

      }
      Object[[j]][i,19]<-ifelse((Object[[j]][i,19])<= (-pi), 2*pi+(Object[[j]][i,19]),(Object[[j]][i,19]))    # adjusting the rel.ang
      Object[[j]][i,19]<-ifelse((Object[[j]][i,19])>= pi,(Object[[j]][i,19])-2*pi,(Object[[j]][i,19]))
      return(Object[[j]][i, 19])
    })
    Object[[j]][1:MM1, 19] <- as.data.frame(res)
  }


  cosine.FP<-data.frame()
  for(j in 1:length(Object)){              # creating values for  cosine based on rel.ang.F
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM, function(i){
      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){   # upper cell going up or lower cell going down
        Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
      }

      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
        Object[[j]][i,20]<-cos(Object[[j]][i,19])
      }
      return(Object[[j]][i,20])
    })

    Object[[j]][1:MM, 20] <- as.data.frame(res)
    cosine.FP[1:MM,j]<-Object[[j]][,20]
  }

  for(j in 1:length(Object)){             ## Forward Pesrsistence time, FP deviating time and FP ratio #
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM1, function(i){
      if(abs(Object[[j]][i,19])<1.5707963268){
        Object[[j]][i,21]= TimeInterval
      }
      if(abs(Object[[j]][i,19])>=1.5707963267){
        Object[[j]][i,21]= 0
      }
      return(Object[[j]][i,21])
    })
    Object[[j]][1:MM1,21] <- as.data.frame(res)
    Object[[j]][MM1,21] <- TimeInterval

  }


  FMResultsTable<-data.frame()                                  # creating a table to store the forward migration results
  for(j in 1:length(Object)){                                   # creating values (NA and forward persistence time )for  forward persistence    (step to the forward movement)
    MM<-Step
    F.P.time<-Object[[j]][1:MM,21]
    MeanF.P.time<-round(mean(F.P.time),digits=2)             # computing mean  Forward Pesrsistence time
    F.P.time.Len<-F.P.time                                   # computing the number of FP steps to be used in computing the FP ratio
    F.P.time.Len[F.P.time.Len==0]<-NA
    F.P.time.Len<-F.P.time.Len[!is.na(F.P.time.Len)]
    F.P.time.Len<-length(F.P.time.Len)
    F.P.Ratio<- round(F.P.time.Len/MM, digits=2)               # computing FP ratio

    DD.F.P.time<-Object[[j]][1:MM,21]                        # computing the FP deviating time
    DD.F.P.time[DD.F.P.time==0]<-1
    DD.F.P.time[DD.F.P.time==TimeInterval]<-0
    DD.F.P.time[DD.F.P.time==1]<-TimeInterval
    MeanDD.F.P.time<-round(mean(DD.F.P.time),digits=2)

    FMResultsTable[1,j]<-j
    FMResultsTable[2,j]<-MeanF.P.time
    FMResultsTable[3,j]<-MeanDD.F.P.time
    FMResultsTable[4,j]<-F.P.Ratio
  }


  VelFPTable<-data.frame()               # creating a table to store the mean velocity with correspondence with the FP time
  for(j in 1:length(Object)){
    MM=Step      				    ###########computing the "finalID" to be used for the splitting
    FPtime<-Object[[j]][1:MM,21]
    FPtime0<-c(0,FPtime)                   #adding a "0" value in the beginning
    FPtime0[FPtime0==TimeInterval]<-NA                #replacing the "5" values with NAs
    FPtime0TF<- !is.na(FPtime0)
    MM1<-Step+1
    rowN<-c(1:MM1)
    Tval <- rowN[FPtime0TF]


    fillIdx <- cumsum(FPtime0TF)
    s<-Tval[fillIdx]                     #replacing the NAs with the previous true value
    justTrueVal<-Tval
    s[FPtime0TF]=0                       #replacing the true values with 0
    finalID<-s[-1]                       # removing the first added value

    Vel<-Object[[j]][1:MM,11]
    FP<-Object[[j]][1:MM,21]
    FP[FP==0]<-NA
    tab<-data.frame(Vel,FP,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]               # to exclude the last row since it has no velocity
    tabb<-split(tab, finalIDD)           # to avoid taking the last row

    meanVEL<-c()
    FP.time<-c()

    res1 <- sapply(1:length(tabb), function(i){
      meanVEL= round(sqrt(mean(tabb[[i]][,1])),digits=2)
      return(meanVEL)
    })
    meanVEL= res1
    meanVEL=meanVEL[-1]

    res2 <- sapply(1:length(tabb), function(i){
      FP.time=sum(tabb[[i]][,2])
      return(FP.time)
    })
    FP.time= res2
    FP.time=FP.time[-1]
    w<-which.max(FP.time)
    FMResultsTable[5,j]<-FP.time[w]



    FPtime<-Object[[j]][1:MM,21]        # computing the "meanVel.for0"
    FPtime00<-c(1,FPtime)                   #adding a "1" value in the beginning  this value should not be 0 because we will not catch 0 persistence if it is in the beginning.
    FPtime00[FPtime00==0]<-NA                #replacing the "0" values with NAs
    FPtime00TF<- !is.na(FPtime00)

    MM1<-MM+1
    rowN<-c(1:MM1)
    Tval0 <- rowN[FPtime00TF]
    TTT<-c()
    res3 <- sapply(1:(length (Tval0)-1), function(i){
      TTT<- (Tval0[i+1]- Tval0[i])-1
      return(TTT)
    })

    TTT=res3
    TTT[TTT==0]<-NA
    F.ID.for.0.FP<-TTT[!is.na(TTT)]
    Final.ID.for.0.FP<-c()
    for (i in 1:length(F.ID.for.0.FP)){
      Final.ID.for.0.FP<-c(Final.ID.for.0.FP,rep(i,(F.ID.for.0.FP[i])))
    }

    Vel<-Object[[j]][1:MM,11]
    FP<-Object[[j]][1:MM,21]
    FP[FP==0]<-NA

    tab<-data.frame()
    tab<-data.frame(Vel,FP,finalID)
    tab<-tab[-nrow(tab),]
    finalIDD<-finalID[-MM]
    tabb<-split(tab, finalIDD)     #to avoid taking the last row
    t<-tabb[[1]]
    tabb1<-cbind(t,Final.ID.for.0.FP)
    tabbb<-split(tabb1, Final.ID.for.0.FP)
    meanVel.for0<-c()

    res4 <- sapply(1:length(tabbb), function(i){
      meanVel.for0<-round(sqrt(mean(tabbb[[i]][,1])),digits=2)
      return(meanVel.for0)
    })
    meanVel.for0=res4

    zerooFP<-rep(0,length(meanVel.for0))
    FP.time1<-c(zerooFP,FP.time)
    meanVEL1<-c(meanVel.for0,meanVEL)
    colNumP<-j+j-1
    colNumV<-j+j
    rowNum<-c(1:length(meanVEL1))
    VelFPTable[rowNum,colNumP]<-FP.time1
    VelFPTable[rowNum,colNumV]<-meanVEL1
    c<-cor.test( ~ VelFPTable[,j+j-1]+ VelFPTable[,j+j], method = "spearman",exact=FALSE)             #testing the correlation
    cc<-unlist(c[4])
    ccPV<-round(cc, digits = 3)
    FMResultsTable[6,j]<-ccPV               # Speed vs  persistence time  Spearman correlation



    if ( sfptPLOT == TRUE || sfptPLOT == T){
      jpeg(paste0(ExpName," FP Time vs Speed",j,".jpg"))
      plot(VelFPTable[,j+j],VelFPTable[,j+j-1],pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Mean Speed during FP time (um/min)",col=color[j],las=1)
      reg<-lm(VelFPTable[,j+j]~VelFPTable[,j+j-1])
      #abline(reg,untf=F,col="red")
      title(main=paste0("Cell Number  ", j,"   Speed vs Forward Persistence Time"),cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccPV),col.sub="red")
      dev.off()

    }



  }

  ## All cells (FP times  vs Speed)
  allper<-VelFPTable[,1]
  allvel<-VelFPTable[,2]
  for(j in 1:length(Object)){
    allper<-c(allper,VelFPTable[,j+j-1])
    allvel<-c(allvel,VelFPTable[,j+j])
  }
  allper<-allper[!is.na(allper)]
  allvel<-allvel[!is.na(allvel)]
  all.per.vel.table<-data.frame(allper,allvel)
  reg<-lm(allper~allvel)
  c<-cor.test( ~ allper+ allvel, method = "spearman",exact=FALSE)                 #testing the correlation
  cc<-unlist(c[4])
  ccP<-round(cc, digits = 3)
  FMResultsTable[6,(length(Object)+1)]<-ccP                                          # Speed vs  persistence time  Spearman correlation

  if ( afptPLOT == TRUE || afptPLOT == T){
    jpeg(paste0(ExpName," FP Time vs Speed - All Cells.jpg"))
    plot(allvel,allper,type="p",pch=16,ylab="FP Time (min)",xlab=" Mean Speed during FP time (um/min)",col="black",las=1)
    abline(reg,untf=F,col="red")
    title("Speed vs FP Time (All cells)",cex.main = 1,sub=paste0("Spearman's rank correlation coefficient = ",ccP),col.sub="red")
    dev.off()

  }

  for(j in 1:length(Object)){    # calculating the Mean.Square.speed for each cell
    MM<-Step
    MM2<-MM-1
    Root.Mean.Square.Speed<-round(sqrt(mean(Object[[j]][1:MM2,11])),digits = 3)
    FMResultsTable[7,j]<-Root.Mean.Square.Speed

    mean.cosineFP<-round(mean(Object[[j]][1:MM2,20],na.rm = TRUE),digits = 3)
    FMResultsTable[8,j]<-mean.cosineFP
    s<-cor.test( ~ sqrt(Object[[j]][1:MM2,11])+ Object[[j]][1:MM2,20], method = "spearman",exact=FALSE)                 #testing the correlation
    ss<-unlist(s[4])
    VEvsCOSP<-round(ss, digits = 3)
    FMResultsTable[9,j]<-VEvsCOSP

    if ( sfpPLOT == TRUE || sfpPLOT == T){
      jpeg(paste0(ExpName," FP vs Speed",j,".jpg"))
      plot(sqrt(Object[[j]][1:MM2,11]),Object[[j]][1:MM2,20],pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Instantaneous Speed (um/min)",col="black",las=1)
      reg<-lm(Object[[j]][1:MM2,20]~sqrt(Object[[j]][1:MM2,11]))
      abline(reg,untf=F,col="red")
      title(main=paste0("Cell Number  ", j,"   Speed vs Forward Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
      dev.off()

    }

  }
  RM<-round(rowMedians(as.matrix(cosine.FP[1:(Step-1),]),na.rm = TRUE),digits=3)
  Speed<-data.frame()
  for (j in 1:length(Object)){    # calculating the Mean.Square.velocity for each cell
    MM<-Step
    MM2<-MM-1
    Speed[1:MM2,j]<-round(sqrt(Object[[j]][1:MM2,11]),digits = 3)
  }
  RowmeanSpeed<-round(rowMedians(as.matrix(Speed),na.rm = TRUE),digits=3)
  s<-cor.test( ~ RM+ RowmeanSpeed, method = "spearman",exact=FALSE)                 #testing the correlation
  ss<-unlist(s[4])
  VEvsCOSP<-round(ss, digits = 3)
  FMResultsTable[9,(length(Object)+1)]<-VEvsCOSP

  if ( afpPLOT == TRUE || afpPLOT == T){
    jpeg(paste0(ExpName," All Cells FP vs Speed.jpg"))
    plot(RowmeanSpeed,RM,pch=16,type="p",ylab="Forward Persistence Time (min)",xlab=" Instantaneous Speed (um/min)",col="black",las=1)
    reg<-lm(RM~RowmeanSpeed)
    abline(reg,untf=F,col="red")
    title(main=paste0("All Cells Speed vs Forward Persistence "),cex.main = 1,sub=paste0("spearman's rank correlation coefficient = ",VEvsCOSP),col.sub="red")
    dev.off()
  }
  RM1<-round(rowMedians(as.matrix(FMResultsTable),na.rm = TRUE),digits=3)
  FMResultsTable[c(2:5,7:8),(length(Object)+1)]<-RM1[c(2:5,7:8)]
  FMResultsTable[1,(length(Object)+1)]<-"All Cells"
  rownames(FMResultsTable)<-c("Cell Number","Mean Forward Persist Time (min)","Mean Forward Persist Deviating Time (min)","Forward Persistence Ratio",
                              "Maximum Forward Persistence period","Forward Persistence Time vs Speed (SCC)","RMSS (um per min)","Mean Forward Angular Persistence (mean cos.F)","Instantaneous Speed vs Forward Persistence (SCC)")
  FMResultsTable<-FMResultsTable[-7,]
  object@ForMigtable=FMResultsTable


  setwd(d)
  write.csv(FMResultsTable, file = paste0(ExpName,"-FMResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-FMResultsTable.csv" ),"in your directory [use getwd()]","\n")
  return(object)
}
#'
#'
#'
#'
#'
#'
#'
#' Method Forward Migration Index
#' @title Forward Migration Index
#'
#' @description The FMI function automatically generates data for the forward migration index
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @return  An CellMig class Object with a data frame. The data frame is stored in the FMItable slot.
#'
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' wsaTD <- CellMig(WSAdataset)
#' wsaTD <- wsaPreProcessing(wsaTD)
#' wsaTD <-FMI(wsaTD,TimeInterval=10,ExpName="ExpName")
#' }

FMI= function(object, TimeInterval=10,ExpName="ExpName"){
  if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
  Object<-object@preprocessedDS
  UPorDO<-object@cellpos
  msg <- NULL
  if ( ! is.list(Object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  d=getwd()
  Len<-length(Object)
  Step<-length(Object[[1]][,1])
  color <-c()
  if (Len> 1023){
    colnum= Len-1023
    color1 <-rainbow(1023)
    colo2 <-rainbow(colnum)
    color=c(color1 ,colo2)
  }else{
    color <- rainbow(Len)
  }
  for(j in 1:length(Object)){                             # Defining if the cell is ubove (1) the wound or below (0) it
    Object[[j]][,25]<-UPorDO[j]
  }

  FMIResultsTable<-data.frame()
  for(j in 1:length(Object)){                             # calculating the cumsum of distance for each cell
    MM<-Step
    MM1<-MM-1

    res <- sapply(1:MM, function(i){
      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]>0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]<0)){   # upper cell going up or lower cell going down
        Object[[j]][i,20]<-(-1*abs(cos(Object[[j]][i,19])))
      }

      if((Object[[j]][1,25]==0) && (Object[[j]][i,5]<0) || (Object[[j]][1,25]==1) && (Object[[j]][i,5]>0)){
        Object[[j]][i,20]<-cos(Object[[j]][i,19])
      }
      return(Object[[j]][i,20])
    })

    Object[[j]][1:MM, 20] <- as.data.frame(res)
    end<-cbind(Object[[j]][MM,2],Object[[j]][MM,3])       # finding the cordinates of the final point in the track.
    start<-cbind(0,0)
    final.dis=dist2(start, end)
    alpha<-acos(Object[[j]][MM,2]/final.dis)
    bita<- 1.5707963268 - alpha
    FMI<- cos(bita)
    if(Object[[j]][1,25]==0 & Object[[j]][MM,3]>0){
      FMI<-(-FMI)
    }

    if(Object[[j]][1,25]==1 & Object[[j]][MM,3]<0){
      FMI<-(-FMI)
    }

    y=abs(Object[[j]][MM,3])
    cumDis=Object[[j]][MM1,12]
    FMIy=round(y/cumDis,digits=3)
    MMM<-round(MM/2)
    mid<-cbind(Object[[j]][MMM,2],Object[[j]][MMM,3])       # finding the cordinates of the final point in the mid track.
    start<-cbind(0,0)
    mid.dis=dist2(start, mid)
    alphaM<-acos(Object[[j]][MMM,2]/mid.dis)
    bitaM<- 1.5707963268 - alphaM
    MTFMI<- cos(bitaM)
    if(Object[[j]][1,25]==0 & Object[[j]][MMM,3]>0){
      MTFMI<-(-MTFMI)
    }

    if(Object[[j]][1,25]==1 & Object[[j]][MMM,3]<0){
      MTFMI<-(-MTFMI)
    }

    p1<-Object[[j]][,20]
    returns<-subset(p1,p1<(-0.87))      # greater than 150 degrees

    yy=abs(Object[[j]][MMM,3])
    cumMDis=Object[[j]][round(MM1/2),12]
    MTFMIy=round(yy/cumMDis,digits=3)

    FMIResultsTable[1,j]<-j
    FMIResultsTable[2,j]<-round(FMI,digits=3)
    FMIResultsTable[3,j]<-round(FMIy,digits=3)
    FMIResultsTable[4,j]<-round(MTFMI,digits=3)
    FMIResultsTable[5,j]<-round(MTFMIy,digits=3)
    FMIResultsTable[6,j]<-round(abs(Object[[j]][MM,3]),digits=2)
    FMIResultsTable[7,j]<-round(length(returns))
  }
  for(j in 1:length(Object)){                            # creating values for  rel.ang.F  (step to the original)
    MM<-Step
    if((Object[[j]][1,25]==0) && (Object[[j]][MM,3]>0) || (Object[[j]][1,25]==1) && (Object[[j]][MM,3]<0)){
      FMIResultsTable[6,j]<- (-1) * FMIResultsTable[6,j]
    }
    else {
      FMIResultsTable[6,j]<- FMIResultsTable[6,j]
    }

  }
  RM1<-round(rowMedians(as.matrix(FMIResultsTable),na.rm = TRUE),digits=3)
  FMIResultsTable[,(length(Object)+1)]<-RM1
  FMIResultsTable[6,(length(Object)+1)]<-round(FMIResultsTable[6,(length(Object)+1)])
  FMIResultsTable[1,(length(Object)+1)]<-"All Cells"

  rownames(FMIResultsTable)<-c("Cell Number","FMI","FMIy", "MTFMI","MTFMIy","Deepness (um)","Number of backwards")

  write.csv(FMIResultsTable, file = paste0(ExpName,"-FMIResultsTable.csv"))
  cat("Results are saved as: ",paste0(ExpName,"-FMIResultsTable.csv" ),"in your directory [use getwd()]","\n")
  object@FMItable=FMIResultsTable

  return(object)

}
#'
#' Method Final Results
#' @title Final Results
#'
#' @description The FinRes function automatically generates a data frame that contains all the results.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ParCor A logical vector that allows generating a correlation table. Default is TRUE.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data
#' @return  A data frame that contains all the results.
#'
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' wsaTD <- CellMig(WSAdataset)
#' wsaTD <- wsaPreProcessing(wsaTD)
#' wsaTD <-FinRes(wsaTD,ExpName="ExpName",ParCor=FALSE)
#' }

FinRes= function(object,ExpName="ExpName",ParCor=TRUE){
  msg <- NULL
  if ( ! is.list(object) ){
    msg <- c(msg, "Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  cells<-c()
  if (length(object@DRtable)>0){
    cells<-object@DRtable[1,]
    DR<-object@DRtable[-1,]
    object@results <- rbind(object@results,DR)
  }

  if (length(object@MSDtable)>0){
    cells<-object@MSDtable[1,]
    MSD<-object@MSDtable[-1,]
    object@results <- rbind(object@results,MSD)
  }

  if (length(object@PerAanSpeedtable)>0){
    cells<-object@PerAanSpeedtable[1,]
    per<-object@PerAanSpeedtable[-1,]
    object@results <- rbind(object@results,per)
  }


  if (length(object@DACtable)>0){
    cells<-object@DACtable[1,]
    DAC<-object@DACtable[-1,]
    object@results <- rbind(object@results,DAC)
  }

  if (length(object@VACtable)>0){
    cells<-object@VACtable[1,]
    VAC<-object@VACtable[-1,]
    object@results <- rbind(object@results,VAC)
  }

  if (length(object@ForMigtable)>0){
    cells<-object@ForMigtable[1,]
    FM<-object@ForMigtable[-1,]
    object@results <- rbind(object@results,FM)
  }


  if (length(object@FMItable)>0){
    cells<-object@FMItable[1,]
    FMI<-object@FMItable[-1,]
    object@results <- rbind(object@results,FMI)
  }

  Results<-object@results
  colnames(Results)<- cells
  colnames(object@results)<- cells
  write.csv(Results, file = paste0(ExpName,"-Final_Results.csv"))
  if ( ParCor == TRUE || ParCor == T){
    R=Results
    R[,]=lapply(R[,], as.numeric)
    Parameters.Correlation<-rcorr(t(R), type="spearman")
    object@parCor<-Parameters.Correlation$r
    write.csv(Parameters.Correlation$r, file = paste0(ExpName,"-Parameters.Correlation.csv"))
    cat("Parameters Correlation table is saved as: ",paste0(ExpName,"-Parameters.Correlation.csv"),"in your directory [use getwd()]","\n")
  }

  cat("The table of the final results is saved as: ",paste0(ExpName,"-Final_Results.csv")," in your directory [use getwd()]","\n")
  cat("\n", "These are the parameters in your final results:","\n")
  print(rownames(Results))

  return(object)
}

#' Method PCA
#' @title PCA
#'
#' @description The CellMigPCA function automatically generates Principal Component Analysis.
#' @param object \code{CellMig} class object, which is a list of data frames resulted from the PreProcessing.
#' @param ExpName A character string. The ExpName will be appended to all exported tracks and statistics data.
#' @param parameters A numeric vector contains the parameters to be included in the Principal Component Analysis. These numbers can be obtained from the outcome of the FinRes() function.

#' @return  PCA Graph of cells and PCA Graph of variables.
#'
#' @export
#'
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' wsaTD <- CellMig(WSAdataset)
#' wsaTD <- wsaPreProcessing(wsaTD)
#' wsaTD <-FinRes(wsaTD,ExpName="ExpName",ParCor=FALSE)
#' }
#'
CellMigPCA= function(object,ExpName="ExpName",parameters=c(1,2,3)){

  if ( ! is.list(object) ){
    stop("Input data must be a list. Please run the PreProcessing step first either rmPreProcessing() or wsaPreProcessing()")
  }
  if ( length(object@results[,1])<1 ){
    stop("There are no results stored. Please run trajectory analysis first")
  }

  if ( length(parameters)<2){
    stop("At least two parameters are required to run the PCA")
  }
  df1<- object@results
  df1=df1[,-length(object@results[1,])]    #### excluding the last column since it is the avarage of all the cells
  tt<-t(df1)
  tt1=tt[,parameters]
  res <- PCA(tt1)
}
