---
output: raw mass features, putative ligands, and volcano plot
---
The current package is used to identify unknown ligands of proteins from E. coli proteins.
by Hui Peng (UofT), September 25, 2023. 
This example uses the raw data to show the identification of unknown ligands using custom algorithm. Only positive ion mode is analyzed.

Step 1: loading packages
```{r, message=FALSE, warning=FALSE}
  rm(list=ls())
  library(xcms)
  library(EcoliPMI)
  library(devtools)
  library(xlsx)
  library(ChemmineR)
  library(rcdk)
  library(MassSpecWavelet)
  library(Rcpp)
  library(RcppArmadillo)
  library(isopat)
  library(readxl)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(ggrepel)
  library(openxlsx)
  #'set up the path
  path<-getwd()
  
  #------------
  #define the polarity for database searching
  #positive: 1; negative: -1
  #we herein use positive ion mode as an example
  polarity<-1
  setwd(path)
```

Step 2: calibrating mass spectra using lock MS. This step is to further calibrate the mass spectra to <1.5 ppm. If the instrument was externally calibrated, this step could be skipped.
```{r, include= FALSE}
  #MassCalibration(polarity)
```

Step 3: extracting peak features from untargeted metabolomics data
```{r, message=FALSE, warning=FALSE}
  #-------------------------
  #exracting peaks using the Peakextract function
  #10^3 is the intensity cutoff, any peaks below this cutoff were removed as they   are close to instrument detection limits
  #---------------------------
  peaks.raw<-Peakextract(10^3,2.5)
  
  #---------------------------------------
  #This step aims to find out the differentiated metabolites by comparing against 
  #'data has to be organized like NR_XX_X, WT_XX_X, NR_DMSO_X
  Control<-'MeOH'
  peaks<-LigandFeatureID(peaks.raw,'neg',Control)#This is the function extracting differential peaks compared to controL
```
Step 4: identifying primary isotopic peaks, and remove redundant isotopic peaks
```{r, message=FALSE, warning=FALSE}
  #2 ppm is the mass tolerance to find isotopic peaks
  #'0.80 is the correlation coefficient to extract isotopic peaks
  #'3 is the intensity ratio cutoff between primary and other isotopic peaks
  peaks.iso<-FindIsotope(peaks,2,0.80,2)
```
Step 5: loading database, and do initial searching with MS1
```{r, message=FALSE, warning=FALSE}
  setwd(path)
  #----------------
  #import MS database, HMDB and EcoCyc are both recommended, we herein use ECOCYC   database as an example
  Database<-read_excel("Ecocyc.xlsx")
  
  #--------------------
  #this is a simple version of mass calibration using a known compound, we assume   that the mass shift is consitent across the whole mass range which is typically   true in positive ion mode
  peaks.iso$mz<-peaks.iso$mz*(1+0*10^(-6))
  
  #-------------------------------------------
  #database searching, 2 is the ppm mass tolerance for database searching
  #depending on the mass calibration and MS stability, the mass tolerance number    should be revisited accordingly
  Library<-InitialSearch(peaks.iso,2,polarity,Database)
  write.table(Library,file='Metabolites.csv',sep=',',row.names = FALSE)
```
Step 6: extracting fragments and isotopic peaks
```{r, message=FALSE, warning=FALSE}
  setwd(path)
  
  #---------------------------------------
  #This is the function to automatically extract MS2 spectra from DIA window for    each compound   feature. Whenever MS2 spectra is availalbe, it is always preferred to increase the confidence of compound identification
  #10^(-5) is the mass tolerance for feature matching
  #15 is the DIA window width, which should be changed according to MS settings
  FragList<-GetFrag(Library,10^(-5),15)
  setwd(path)
  
  #---------------------------------------
  #isotopic pattern is another useful metrices for compound identification
  #this function is used automatically extract isotopic peaks for each compound     feature
  #0.8 is the correlation coefficient cutoff for isotopic peak matching
  Isotope.Data<-IsotopeFind(Library,10^(-5),0.8)
```
Step 7: extracting adducts. We herein extract the adducts information for each compound, as these adducts information is also useful to increase identificaiton confidence. For instance, if [M+H]+ and [M+NH4]+ are both detected, it largely increases our confidence that M is the molecular ion rather than in source fragments.
```{r, message=FALSE, warning=FALSE}
  setwd(path)
  if (polarity==1){
  Adducts<-c('[M+H]+','[M+NH4]+')
  MW.adducts<-c(1.007825,18.03437)
  Adducts.Find<-AdductsFind(Library,MW.adducts,2*10^(-6),Adducts)
}

if (polarity==-1){
  Adducts<-c('[M-H]-','M+FA-H')
  MW.adducts<-c(-1.007825,44.99765)
  Adducts.Find<-AdductsFind(Library,MW.adducts,2*10^(-6),Adducts)
}
```
Step 8: calculating scores, by combining MS1, MS2, isotopic peak, and adduct information
```{r, message=FALSE, warning=FALSE}
  setwd(path)
  cutoff<-5000
  
  #------------------------------------
  #this is the function for scoring
  #'ms2true' means that MS2 fragment is available for score calculation
  #'ms2false' means that MS2 spectra are not available
  mylib.Target<-DatabaseSearching(Library,cutoff,polarity,Database,Isotope.Data,FragList,iso_list,'ms2true')
```
Step 9: weighted scores
```{r, message=FALSE, warning=FALSE}
  setwd(path)
  #' the path to save results
  path<-getwd()
  #' the path to save raw data
  path.data<-paste0(path,"/data")
  setwd(path.data)
  rawfiles<-list.files()
  xraw<-xcmsRaw(rawfiles[1],includeMSn = 1)
  precursor<-Preclist(xraw)
  
  #--------------------------
  #we can assign different scores to each chemometric information
  weightK<-c(1,50,1,1,1,1)#weight for scores,the second for ms2 score
  setwd(path)
  mylib.score<-Finalscore(mylib.Target,weightK,precursor)
  
  #--------------------------
  #this function aims to predict the retention times of each detected feature
  Target.rt<-ScoreRT(mylib.score,Database)
  output<-Output(Target.rt,-1000)
  
  write.table(output,file='BirA_results.csv',sep=',',row.names = FALSE)
 ``` 
