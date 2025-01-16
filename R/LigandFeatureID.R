#' This function is used for initial searching
#'
#' @param mydata
#' @param Control
#'
#' @return
#' @export
#'
#' @examples

LigandFeatureID<-function(mydata,Fold_cutoff){


  header <- names(mydata)
  mydata$ID<-1:nrow(mydata)

  pattern <- "^[^_]+_[^_]+_\\d+\\.mzXML$"

  header_list <- header[grep(pattern, header, invert = FALSE)]

  Anno <- data.frame(
    File = header_list,
    Protein = sapply(strsplit(header_list, "_"), function(x) x[1]),
    Ion = sapply(strsplit(header_list, "_"), function(x) x[2]),
    Replicate = sapply(strsplit(header_list, "_"), function(x) gsub("\\.mzXML", "", x[3]))
  )

  Anno$Replicate <- as.numeric(Anno$Replicate)


  ###### 5.Determine numbers of protein target

  prname <- unique(Anno$Protein)

  prnum <- length(prname)

  ###### 6.Statistics

  mydata$FOLD<-rep(0,nrow(mydata))
  mydata$PVALUE<-rep(1,nrow(mydata))
  mydata <- as_tibble(mydata)
  index_save<-NULL

  for (i in 1:prnum) {

    ###### Select protein target
    prni <- prname[i]
    column_list <- grep(prni, header, invert = FALSE)

    for (j in 1:nrow(mydata)){
      test<-mydata[j,column_list]
      ctrl<-mydata[j,-c(1,2,column_list,ncol(mydata),ncol(mydata)-1,ncol(mydata)-2)]

      #t test
      fold<-(sum(test)*length(ctrl))/(sum(ctrl)*length(test))
      if (fold<Fold_cutoff){
        next
      }

      if (sd(test) > 0){
        ttest <- t.test(test, ctrl)
        pvalue <- ttest$p.value
        if (pvalue>0.05){
          next
        }
        mydata$FOLD<-fold
        mydata$PVALUE<-pvalue
        index_save<-c(index_save,j)
      }
    }
  }

  ##output the figures

  return(mydata[index_save,])


}


