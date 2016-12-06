#' Writes the contents of a slot in TSRdata to a BED file
#' @param expName an object of class tssExp with TSR data loaded
#' @param Name the name of the BED file to be written in your working directory
#' @export

setGeneric(
           name="tsrToBed",
           def=function(expName, Name, writeBed) {
               standardGeneric("tsrToBed")
    }
    )

setMethod("tsrToBed",
          signature(expName="tssExp", Name="character", writeBed="logical"),

          function(expName, Name="tsrName", writeBed=TRUE) {
              if (length(expName@tsrData) == 0) {
                  stop("Slot @tsrData is empty.\n\n Please process TSR data before running this command.")
              }

              expName@tsrData -> my.TSR

              message("\nConverting TSRs to BED format...\n")

              my.TSR$plus -> tsr.list.p
              cbind(as.character(tsr.list.p[[1]])) -> tsr.df.p
                  for (j in 1:length(tsr.list.p)) {
                      as.character(tsr.list.p[[j]]) -> this.tsr
                      rbind(tsr.df.p, this.tsr) -> tsr.df.p
                  }

              my.TSR$minus -> tsr.list.m
              cbind(as.character(tsr.list.m[[1]])) -> tsr.df.m
                  for (j in 1:length(tsr.list.m)) {
                      as.character(tsr.list.m[[j]]) -> this.tsr
                      rbind(tsr.df.m, this.tsr) -> tsr.df.m
                  }

              tsr.df <- rbind(tsr.df.p, tsr.df.m)

              tsr.split <- unlist(strsplit(tsr.df, split=":"))

              tsr.mat <- matrix(unlist(tsr.split), ncol=3, byrow=TRUE)

              tsr.coord <- matrix(unlist(strsplit(tsr.mat[,2], split="-")), ncol=2, byrow=TRUE)
              new.df <- cbind(tsr.mat[,1], tsr.coord, c("."), c("."), tsr.mat[,3])

              colnames(new.df) <- c("chr", "start", "end", "ID", "score", "strand")

              df.out <- as.data.frame(new.df)

              if (writeBed==TRUE) {
                  fileName <- paste(Name, "bed", sep=".")
                  write.table(df.out, file=fileName, sep="\t", row.names=FALSE, quote=FALSE)
              }

              return(df.out)

          }
          )
