#' This function is used to extract peak features
#'
#' @param Intensitycut 
#' @param ppm 
#' @return peak features
#' @export
#' @import xcms
#' @description This function takes raw data and then extract all peaks to dataframe

Peakextract <- function(Intensitycut,ppm) {
#' the path to save results
  path<-getwd()
  
#' the path to save raw data
  path.data<-paste0(path,"/data")
  
#' ----------------------------------------------
#' detect peaks from mass spec raw files
#' ---------------------------------------------
  setwd(path.data)
  msfiles<-list.files()
  xset.raw<-xcmsSet(msfiles,method='centWave',ppm=ppm,peakwidth=c(10,30),snthresh=10,nSlaves=1)
  xtest<-xcmsRaw(msfiles[1])
  mzrange<-xtest@mzrange

#' delete peak features with extreme m/z values
  xset<-xset.raw
  index<-which(xset.raw@peaks[,1]<mzrange[1]+1)
  if (length(index)>0){
    xset@peaks<-xset.raw@peaks[-index,]}
  index<-which(xset@peaks[,1]>mzrange[2]-1)
  if (length(index)>0){
    xset@peaks<-xset@peaks[-index,]}
  peaklist<-xset@peaks                                                                    
  len<-length(xset@peaks[,1])
  
#' group peaks across samples
  xset1<-group(xset,bw=60,minsamp=1,minfrac=1/length(msfiles),mzwid=0.001)
  
#' filling missing values
  xset2<-Fillpeak(xset1,10,20,msfiles)

#' peak ID for each group
  test<-unlist(xset2@groupidx)
  len<-length(xset2@groupidx)
  len2<-length(msfiles)
  
#' creating peak matrix
  Allpeak<-array(rep(0,len*(len2+2)),dim=c(len,(len2+2)))
  for (i in 1:len){
    print(c('PeakID...',i,'of...',len))
    temp<-unlist(xset2@groupidx[i])
    len3<-length(temp)
    for (j in 1:len3){
      index1<-xset2@peaks[temp[j],11]
#' m/z vlaues      
      Allpeak[i,1]<-xset2@peaks[temp[j],1]
      
#' retention time in minutes
      Allpeak[i,2]<-xset2@peaks[temp[j],4]/60

#‘ peak intensity
      Allpeak[i,index1+2]<-max(xset2@peaks[temp[j],9],Allpeak[i,index1+2])
    }}
  
#' replacing 0 values
  Allpeak[which(Allpeak==0)]<-100

#' extracting peak intensities across all samples
  colnames(Allpeak)<-c('mz','rt',msfiles)
  Allpeak<-data.frame(Allpeak)
  Allpeak$SampleID<-rep(1,nrow(Allpeak))
  
#' extracting peak intensities
  index.save<-NULL
  for (i in 1:nrow(Allpeak)){
    index<-which.max(Allpeak[i,3:ncol(Allpeak)])
    if (max(Allpeak[i,3:ncol(Allpeak)])>Intensitycut){
      index.save<-c(index.save,i)
      }
    Allpeak$SampleID[i]<-index[1]
  }
  if(length(index.save)>0){
    Allpeak<-Allpeak[index.save,]
  }
  
  
#'reevaluate ms peaks
  Library.new<-mzSmooth(Allpeak,msfiles,90,10)
  index<-which(is.na(Library.new$mz))
  if (length(index)>0){
    Library.new<-Library.new[-index,]
    }  

#' setup the working folder
  setwd(path)
  
  return(Library.new)
}
#‘ devtools::document()



#' This function is used to smooth mz and retention time to get more accurate mz and RT
#'
#' @param Library 
#' @param msfiles 
#' @param rtwin 
#' @param ppmwin 
#'
#' @return
mzSmooth<-function(Library,msfiles,rtwin,ppmwin){
  ppmwin<-ppmwin*10^(-6)
  Library.new<-Library
  for (i in 1:length(msfiles)){
    index<-which(Library$SampleID==i)
    if (length(index)==0){next}
    xrawdata<-xcmsRaw(msfiles[i])
    for (j in 1:length(index)){
      temp.index<-index[j]
      mz<-Library$mz[temp.index]
      rt<-Library$rt[temp.index]*60
      
      #'smooth the exact mass and retention time
      mz.min<-mz*(1-ppmwin)
      mz.min<-max(mz.min, xrawdata@mzrange[1])
      mz.max<-mz*(1+ppmwin)
      mz.max<-min(mz.max,xrawdata@mzrange[2])
      
      #'range of analysis is the entire chromatograph
      rt.min<-max(min(xrawdata@scantime),rt-rtwin)
      
      #' the first and last 3 minutes are skipped to avoid sections with no analytes
      rt.max<-min(max(xrawdata@scantime),rt+rtwin) 
      peaks<-rawEIC(xrawdata,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max))
      scan.max<-which.max(peaks$intensity)
      scan.max<-peaks$scan[scan.max[1]]
      
      #'delete the data if the compound was detected in the last scan
      if (max(peaks$intensity)<1000){
        Library.new$mz[temp.index]<-NA
        next}
      if (scan.max>=(length(xrawdata@scanindex)-2)){
        Library.new$mz[temp.index]<-NA
        next}
      scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
      correctindex<-(scanNum[1]+1):scanNum[2]
      mzlist<-xrawdata@env$mz[correctindex]
      intensity<-xrawdata@env$intensity[correctindex]
      index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
      
      #'the mz of the maxiaml peak
      mz.scan1<-mzlist[index.mz[1]]
      
      #'the mz for scaning point -1 to peak top
      if (scan.max==1){
        mz.scan2<-mz.scan1}else{
          scanNum<-c(xrawdata@scanindex[scan.max-1],xrawdata@scanindex[scan.max])
          correctindex<-(scanNum[1]+1):scanNum[2]
          mzlist<-xrawdata@env$mz[correctindex]
          intensity<-xrawdata@env$intensity[correctindex]
          index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
          
          #' the mz of the maxiaml peak
          mz.scan2<-mzlist[index.mz[1]]
          }
      
      #'the mz for scaning point 1 to peak top
      scanNum<-c(xrawdata@scanindex[scan.max+1],xrawdata@scanindex[scan.max+2])
      correctindex<-(scanNum[1]+1):scanNum[2]
      mzlist<-xrawdata@env$mz[correctindex]
      intensity<-xrawdata@env$intensity[correctindex]
      index.mz<-which(abs(mzlist-mz)<ppmwin*mz)
      
      #'the mz of the maxiaml peak
      mz.scan3<-mzlist[index.mz[1]]
      
      #######replace the mz and rt with newly smoothed results    
      Library.new$mz[temp.index]<-sum(mz.scan1,mz.scan2,mz.scan3)/3
      Library.new$rt[temp.index]<-xrawdata@scantime[scan.max]/60
    }
    }
  return(Library.new)
}

#' Fill missing values across samples
#'
#' @param xset 
#' @param ppm 
#' @param RT 
#'
#' @return a complete peak feature matrix
#' 
Fillpeak<-function(xset,ppm,btw,msfiles){
  ppm<-ppm/10^6
  xset.input<-xset
  
  #'peak feature matrix  
  idsave<-matrix(rep(0,length(xset.input@groupidx)*length(msfiles)*4),nrow=length(xset.input@groupidx)*length(msfiles),ncol=4)
  
  #' find peaks across groups  
  k<-1
  for (i in 1:length(xset.input@groupidx)){
    index<-unlist(xset.input@groupidx[i])
    sampleid<-xset.input@peaks[index,11]
    #' the average mz and Rt retention time across samples    
    mz.value<-mean(xset.input@peaks[index,1])
    rt.value<-mean(xset.input@peaks[index,4])
    
    #' find the peaks across samples
    for (j in 1:length(unlist(phenoData(xset.input)))){
      index<-which(sampleid==j)
      if (length(index)<1){
        idsave[k,]<-c(i,j,mz.value,rt.value)##save groupidx,sampleid
        k<-k+1
      }
    }
  }
  
  #' filling missing values  
  index<-which(idsave[,1]==0)
  idsave<-idsave[-index,]
  
  #' no value to fill  
  if(nrow(idsave)==0){
    return(xset.input)
  }
  
  #' fill peaks  
  newpeak<-matrix(rep(0,nrow(idsave)*5),nrow=nrow(idsave),ncol=5)
  kk<-1
  msfile<-filepaths(xset.input)
  minmz<-min(xset.input@peaks[,1])
  maxmz<-max(xset.input@peaks[,1])
  minrt<-min(xset.input@peaks[,4])
  maxrt<-max(xset.input@peaks[,4])
  if (length(idsave)>0){
    for (k in 1:length(unlist(phenoData(xset.input)))){
      print(c('filling peakID...',k,'of...',length(unlist(phenoData(xset.input)))))
      index<-which(idsave[,2]==k)
      if (length(index)==0){next}
      xraw<-xcmsRaw(msfiles[k])
      for (n in 1:length(index)){
        mz.value<-idsave[index[n],3]
        mzmin<-max(minmz,mz.value-mz.value*ppm)
        mzmax<-min(maxmz,mz.value+mz.value*ppm)
        rt.value<-idsave[index[n],4]
        rtmin<-max(minrt,rt.value-btw)
        rtmax<-min(maxrt,rt.value+btw)
        
        #' extract peak intensities        
        peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        intensity<-max(peak$intensity)
        newpeak[kk,]<-c(mz.value,rt.value,intensity,idsave[index[n],1],k)##mz, rt, intensity,groupidx,sampleid
        kk<-kk+1
      }
    }
  }
  
  #' organize all peak features into a new peak matrix
  peakid<-length(xset.input@peaks[,1])
  for (i in 1:nrow(newpeak)){
    groupid<-newpeak[i,4]
    tempvalue<-unlist(xset.input@groupidx[groupid])
    peakid<-peakid+1
    xset.input@groupidx[groupid]<-list(c(tempvalue,peakid))
  }
  peak.combine<-matrix(0,ncol=11,nrow=nrow(newpeak))
  peak.combine[,1]<-newpeak[,1]
  peak.combine[,4]<-newpeak[,2]
  peak.combine[,9]<-newpeak[,3]
  peak.combine[,11]<-newpeak[,5]
  xset.input@peaks<-rbind(xset.input@peaks,peak.combine)
  return(xset.input)
}