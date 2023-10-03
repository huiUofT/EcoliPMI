#' This function is used for initial searching
#'
#' @param mydata
#' @param Control 
#'
#' @return
#' @export
#'
#' @examples

LigandFeatureID<-function(mydata,polarity,ctr){


  header <- names(mydata)
  mydata$ID<-1:nrow(mydata)
  
  pattern <- "CAR"
  
  header_list <- header[grep(pattern, header, invert = FALSE)]
  
  Anno <- data.frame(
    File = header_list,
    Protein = sapply(strsplit(header_list, "_"), function(x) x[1]),
    Sample = sapply(strsplit(header_list, "_"), function(x) x[2]),
    Replicate = sapply(strsplit(header_list, "_"), function(x) gsub("\\.mzXML", "", x[3]))
  )
  
  Anno$Replicate <- as.numeric(Anno$Replicate)
  
  
  ###### 5.Determine numbers of protein target
  
  prname <- unique(Anno$Protein)
  
  prname <- prname[prname != "NC"]
  
  prnum <- length(prname)
  
  ###### 6.Statistics
  
  mydata <- as_tibble(mydata)
  
  datafinal <- mydata[, 0]
  peaks<-NULL
  
  for (i in 1:prnum) {
    
    ###### Select protein target
    
    print(i)
    
    prni <- prname[i]
    
    Annopri <- subset(Anno, Protein == prni)
    
    sampname <- unique(Annopri$Sample)
    
    sampname <- sampname[sampname != ctr]
    
    sampnum <- length(sampname)
    
    for (j in 1: sampnum){
      
      ###### Select environmental sample type
      
      print(j)
      
      data <- mydata
      
      new_column_fc <- rep(1, nrow(data))
      
      new_column_p <- rep(1, nrow(data))
      
      data <- cbind(data,
                    new_column_fc, 
                    new_column_p)
      
      sampnj <- sampname[j]
      
      test_anno <- Anno[which(Anno$Protein == prni), ]
      
      test_anno <- test_anno[which(test_anno$Sample == sampnj), ]
      
      
      control <- Anno[which(Anno$Protein == prni), ]
      
      control <- control[which(control$Sample == ctr), ]
      
      FileID_test <- test_anno$File
      
      FileID_control <- control$File
      
      testdata <- data[, colnames(data)%in%FileID_test]
      
      data.crtl <- data[, colnames(data)%in%FileID_control]
      
      
      ###### Calculate fc and pvalue
      
      for (g in 1:nrow(data)){
        
        print(g)
        
        test <- testdata[g,]
        ctrl <- data.crtl[g,]
        
        fold <- (sum(test)*length(ctrl))/(sum(ctrl)*length(test))
        
        data$new_column_fc[g] <- fold
        
        if (sd(test) > 0){
          ttest <- t.test(test, ctrl)
          pvalue <- ttest$p.value
          data$new_column_p[g] <- pvalue
        } else {data$new_column_p[g] <- 1}
        
      }
      
      
      ###### Find significant peaks
      
      sig <- data %>% dplyr::filter(new_column_fc > 100 & new_column_p < 0.01)
        
 
     
      
      ###### Calculate average intensity for each sig peak   
      
      sig$avg <- rep(1, nrow(sig))
      
      data_avg <- sig[, colnames(sig)%in%FileID_test]
      
      for (k in 1:nrow(sig)){
        
        print(prni)
        
        print(sampnj)
        
        print(k)
        
        sam <- data_avg[k,]
        
        avgsam <- sum(sam)/ncol(sam)
        
        sig$avg[k] <- avgsam
        
      }
      
      ###### Prepare data for volcano plot
      
      dataset <- data
      
      dataset$avg <- rep(100, nrow(dataset))
      
      dataset$change <- rep("stable", nrow(dataset))
      
      sig$change <- rep("up", nrow(sig))
      
      ###### Remove sig peaks since the avg has been calculated somewhere else
      
      dataset <- dataset[-which(dataset$new_column_p < 0.01 & 
                                  dataset$new_column_fc > 100),]
      
      dataset <- bind_rows(dataset, sig)
      
      ###### Plotting   
      
      cut_off_pvalue <- -log10(0.01)
      
      cut_off_logFC <- log10(100)
      
      temp_p <- arrange(sig, desc(avg))
      
      tempn <- paste(prni, sampnj, polarity, sep = "_")
      
      p <- ggplot(
        dataset, aes(x = log10(new_column_fc), y = -log10(new_column_p), colour = change, size = avg)) +
        scale_color_manual(values=c("#231f20", "#be1e2d")) +
        geom_point(alpha = 0.7, shape = 16) +
        geom_hline(yintercept = cut_off_pvalue, col="#be1e2d", linetype = 2) +
        geom_vline(xintercept = cut_off_logFC, col="#be1e2d",linetype = 2) +
        theme_bw() +
        theme(panel.grid.major=element_line(colour=NA)) +
        theme(panel.grid.minor=element_line(colour=NA)) +
        geom_text_repel(
          data=temp_p[1:1,],
          aes(label=round(mz,digits=4)), #use mzmed for xcms peak picking
          size=4,
          color="#be1e2d",
          segment.color="black",
          show.legend = FALSE) +
        theme(legend.position = "none") +
        xlab("Log(Fold change)") +
        ylab("-Log P value") +
        theme(axis.text.x = element_text(size = rel(1.5), color = "black")) +
        theme(axis.text.y = element_text(size = rel(1.5), color = "black")) +
        theme(axis.title.x = element_text(size = rel(1.2), color = "black")) +
        theme(axis.title.y = element_text(size = rel(1.2), color = "black")) +
        ggtitle(tempn)
      
      p
      
      tempn <- paste(prni, sampnj, polarity, sep = "_")
      
      tempn <- paste(tempn, ".tiff", sep = "")
      
      ggsave(p, filename = tempn, width = 6, height = 5)
      
      ###### Save data 
      
      postfix_fc <- "fc"
      
      postfix_p <- "pvalue"
      
      colname1 <- paste(prni, sampnj, postfix_fc)
      
      colname2 <- paste(prni, sampnj, postfix_p)
      
      colnames(data)[colnames(data) == "new_column_fc"] <- colname1
      
      colnames(data)[colnames(data) == "new_column_p"] <- colname2
      
      colnames(sig)[colnames(sig) == "new_column_fc"] <- colname1
      
      colnames(sig)[colnames(sig) == "new_column_p"] <- colname2
      
      tempn <- paste(prni, sampnj, polarity, sep = "_")
      
      tempn <- paste(tempn, ".csv", sep = "")
      
      write.table(sig, file = tempn, sep=',', row.names = FALSE)
      peaks<-c(peaks,sig$ID)
      
      datafinal <- cbind(datafinal, data)
      
      datafinal <- datafinal[, !duplicated(colnames(datafinal))]
      
      sig <- NULL
      
    }
    
  }
  
  tempn <- paste(prname, collapse = "_")
  
  tempn <- paste(tempn, polarity, sep = "_")
  
  tempn <- paste(tempn, ".csv", sep = "")
  
  write.table(datafinal, file = tempn, sep=',', row.names = FALSE)
  return(mydata[unique(peaks),])
}


