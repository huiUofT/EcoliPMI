#' This function is used for initial searching
#'
#' @param mydata
#' @param Control
#'
#' @return
#' @export
#'
#' @examples

LigandFeatureID<-function(mydata,polarity,control){


  header <- names(mydata)
  mydata$ID<-1:nrow(mydata)

  pattern <- "^[^_]+_[^_]+_\\d+\\.mzXML$"

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

  prname <- prname[prname != "WT"]

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

    sampname <- sampname[sampname != control]

    sampnum <- length(sampname)

    for (j in 1: sampnum){

      ###### Select environmental sample type

      print(j)

      data <- mydata

      new_column_fc1 <- rep(1, nrow(data))

      new_column_fc2 <- rep(1, nrow(data))

      new_column_p1 <- rep(1, nrow(data))

      new_column_p2 <- rep(1, nrow(data))

      data <- cbind(data,
                    new_column_fc1,
                    new_column_p1,
                    new_column_fc2,
                    new_column_p2)

      sampnj <- sampname[j]

      test_anno <- Anno[which(Anno$Protein == prni), ]

      test_anno <- test_anno[which(test_anno$Sample == sampnj), ]

      control1 <- Anno[which(Anno$Protein == "WT"), ]

      control1 <- control1[which(control1$Sample == sampnj), ]

      control2 <- Anno[which(Anno$Protein == prni), ]

      control2 <- control2[which(control2$Sample == control), ]

      FileID_test <- test_anno$File

      FileID_control1 <- control1$File

      FileID_control2 <- control2$File

      testdata <- data[, colnames(data)%in%FileID_test]

      data.crtl1 <- data[, colnames(data)%in%FileID_control1]

      data.crtl2 <- data[, colnames(data)%in%FileID_control2]

      ###### Calculate fc and pvalue

      for (g in 1:nrow(data)){

        print(g)

        test <- testdata[g,]
        ctrl1 <- data.crtl1[g,]
        ctrl2 <- data.crtl2[g,]

        fold1 <- (sum(test)*length(ctrl1))/(sum(ctrl1)*length(test))
        fold2 <- (sum(test)*length(ctrl2))/(sum(ctrl2)*length(test))

        data$new_column_fc1[g] <- fold1
        data$new_column_fc2[g] <- fold2

        if (sd(test) > 0){
          ttest1 <- t.test(test, ctrl1)
          pvalue1 <- ttest1$p.value
          data$new_column_p1[g] <- pvalue1
        } else {data$new_column_p1[g] <- 1}

        if (sd(test) > 0){
          ttest2 <- t.test(test, ctrl2)
          pvalue2 <- ttest2$p.value
          data$new_column_p2[g] <- pvalue2
        } else {data$new_column_p2[g] <- 1}

      }


      ###### Find significant peaks

      sig <- data %>% dplyr::filter(new_column_fc1 > 5 & new_column_p1 < 0.05) %>%

        dplyr::filter(new_column_fc2 > 5 &  new_column_p2 < 0.05)

      if (nrow(sig)==0){next}#no significant peak

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

      ###### Remove peaks induced due to overexpression

      dataset <- dataset[-which(dataset$new_column_p2 < 0.05 &
                                  dataset$new_column_fc2 < 0.2),]

      dataset <- dataset[-which(dataset$new_column_p2 < 0.05 &
                                  dataset$new_column_fc2 > 5),]

      ###### Remove sig peaks since the avg has been calculated somewhere else

      dataset <- dataset[-which(dataset$new_column_p1 < 0.05 &
                                  dataset$new_column_fc1 > 5),]

      dataset <- bind_rows(dataset, sig)

      ###### Plotting

      cut_off_pvalue <- -log10(0.05)

      cut_off_logFC <- log10(5)

      temp_p <- arrange(sig, desc(avg))

      tempn <- paste(prni, sampnj, polarity, sep = "_")

      p <- ggplot(
        dataset, aes(x = log10(new_column_fc1), y = -log10(new_column_p1), colour = change, size = avg)) +
        scale_color_manual(values=c("#231f20", "#be1e2d")) +
        geom_point(alpha = 0.7, shape = 16) +
        geom_hline(yintercept = cut_off_pvalue, col="#be1e2d", linetype = 2) +
        geom_vline(xintercept = cut_off_logFC, col="#be1e2d",linetype = 2) +
        theme_bw() +
        theme(panel.grid.major=element_line(colour=NA)) +
        theme(panel.grid.minor=element_line(colour=NA)) +
        geom_text_repel(
          data=temp_p[1:5,],
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

      postfix_fc1 <- "fc1"

      postfix_p1 <- "pvalue1"

      postfix_fc2 <- "fc2"

      postfix_p2 <- "pvalue2"

      colname1 <- paste(prni, sampnj, postfix_fc1)

      colname2 <- paste(prni, sampnj, postfix_p1)

      colname3 <- paste(prni, sampnj, postfix_fc2)

      colname4 <- paste(prni, sampnj, postfix_p2)

      colnames(data)[colnames(data) == "new_column_fc1"] <- colname1

      colnames(data)[colnames(data) == "new_column_p1"] <- colname2

      colnames(data)[colnames(data) == "new_column_fc2"] <- colname3

      colnames(data)[colnames(data) == "new_column_p2"] <- colname4

      colnames(sig)[colnames(sig) == "new_column_fc1"] <- colname1

      colnames(sig)[colnames(sig) == "new_column_p1"] <- colname2

      colnames(sig)[colnames(sig) == "new_column_fc2"] <- colname3

      colnames(sig)[colnames(sig) == "new_column_p2"] <- colname4

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


