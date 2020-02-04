#' Method preprocessing
#'
#' @description The PreProcessing function automatically adjusts the X and Y coordinates based on the physical size of the pixel “PixelSize” and let all cells start from a single point (X=0,Y=0). Furthermore, it computes several descriptive parameters, including distances, absolute angles,relative angles, persistence time, squared speed, cumulative distance, Directional ratio and acceleration.
#'
#' @param object A trajectory data frame organized into four columns: cell ID, X coordinates, Y coordinates and Track number, which is the track's path order.
#' @param PixelSize A numeric value of the physical size of a pixel.
#' @param TimeInterval A numeric value of the time elapsed between successive frames in the time-lapse stack.
#'
#' @return A list of data frames, each data frame shows the trajectories of a single cell.
#'
#' @export
#' @author Salim Ghannoum \email{salim.ghannoum@@medisin.uio.no}
#'
#'
#' @examples
#' \dontrun{
#' data(Trajectory_dataset)
#' df<-Trajectory_dataset
#' preproc<-PreProcessing(df,PixelSize=1.24, TimeInterval=10)
#' }
PreProcessing = function(object,PixelSize,TimeInterval) {
            msg <- NULL
            if ( ! is.data.frame(object) ){
              msg <- c(msg, "input data must be data.frame")
            }
            if ( ! is.numeric(PixelSize) ) stop( "PixelSize has to be a positive number" ) else if ( PixelSize<= 0 ) stop( "PixelSize has to be a positive number" )
            if ( ! is.numeric(TimeInterval) ) stop( "TimeInterval has to be a positive number" ) else if ( TimeInterval<= 0 ) stop( "TimeInterval has to be a positive number" )
            df<-object
            df<-df[,1:3]                                        # Removing the unnecessary columns
            L<-length(df[,1])
            df[,4:26]<-rep(0,L)
            df[,27]<-rep(NA,L)                                  # to be used for migration type
            colnames(df)<-c("ID","x","y","X","Y","dx","dy","dis","abs.ang","rel.ang.P","Cos.P","Persist.Time","Square Speed","cumDis","Dir.R","NewDX","NewDY","New.Abs.ang","Ang.Diff","New.Cos.diff","rel.ang.F","Cos.F","Forward.Persist.Time","MSD(lag)","VAC(lag)","Acceleration","M-type")
            ID_split <- split(df, df$ID)          #Splitting the data frame based on the ID
            cat("This dataset contains: ",length(ID_split),"Cells","\n")

            for(j in 1:length(ID_split)){         # Having the ID =group order
              ID_split[[j]][1]=j
            }

            for(j in 1:length(ID_split)){             # adjusting x and y (starting from 0 & being multiplied by H)
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
            return(PreprocessedData)
          }





