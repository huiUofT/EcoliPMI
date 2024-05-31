#' This function is used to find the primary isotopic peaks
#'
#' @param peaklist 
#' @param mz_tol 
#' @param corrcut 
#' @param intensitycut 
#'
#' @return
#' @export
#'
#' @examples
FindIsotope<-function(peaklist,mz_tol,corrcut,intensitycut){
  #' the path to save results
  path<-getwd()
  
  #' the path to save raw data
  path.data<-paste0(path,"/data")
  

  setwd(path.data)
  msfiles<-list.files()
  mz_tol<-mz_tol*10^(-6)
  index.save<-NULL
  IDsave<-0
  peaklist$PrimaryIso<-peaklist$mz
  
  #'search each mass spec files to find isotopes
  for (i in 1:length(msfiles)){
    print(c('MS file...',i))
    xrawdata<-xcmsRaw(msfiles[i])
    index<-which(peaklist$SampleID==i)
    
    #'no peak in this sample
    if (length(index)<1){next}
    
    #'search isotopes
    for (j in 1:length(index)){
      temp<-index[j]
      mz<-peaklist$mz[temp]
      rt<-peaklist$rt[temp]*60
    
      #'defining searching ranges
      mzrange<-xrawdata@mzrange
      minmz<-mzrange[1]
      maxmz<-mzrange[2]
      mzmin<-max(minmz,mz-mz*mz_tol)
      mzmax<-min(maxmz,mz+mz*mz_tol)
      rtrange<-xrawdata@scantime
      rtmin<-max(min(rtrange),rt-30)
      rtmax<-min(max(rtrange),rt+30)
      
    #'searching the scan number for peak top  
    if (rtmax<0){next}
    peaks<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
    scan.max<-which.max(peaks$intensity)
    scan.max<-peaks$scan[scan.max[1]]
    scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
    correctindex<-(scanNum[1]+1):scanNum[2]
    mzlist<-xrawdata@env$mz[correctindex]
    
    #'mz ranges to find isotopes
    mz.iso<-which(abs(mzlist-mz)<13)
    mz.iso<-mzlist[mz.iso]
    if (length(mz.iso)<1){next}
    
    
    #' searching isotopes using exact masses
    save.iso<-NULL
    for (mm in 1:length(mz.iso)){
      match.iso<-FindHomo(mz.iso[mm],mz)
      
      #'isotopes found
      if(match.iso==1){save.iso<-c(save.iso,mm)}
    }
    if (length(save.iso)<1){next}
    mz.iso<-mz.iso[save.iso]
    
    #'only keep those mz smaller as we are interested in primary isotopic peaks
    index.min<-which(mz.iso<(mz-1))
    if (length(index.min)<1){
      next
    }
    mz.iso<-mz.iso[index.min]
    
    #'using correlations to further find isotopes
    native.peak<-peaks$intensity
    kk<-0
    save.iso<-NULL
    for (k in 1:length(mz.iso)){
      mz0<-mz.iso[k]
      mzmin<-max(minmz,mz0-mz0*mz_tol)
      mzmax<-min(maxmz,mz0+mz0*mz_tol)
      isotope.peak<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      isotope.peak<-isotope.peak$intensity
      
      #'if the intensity of the primary peak is too low, peaks are excluded
      if (max(isotope.peak)<(intensitycut*max(native.peak))){
        next
      }
      
      #'calculating correlation scores
      if (sd(isotope.peak)==0||sd(native.peak)==0){next}
        corr<-cor(native.peak,isotope.peak)
        
        #'isotopes were supported by correlations
        if (corr>corrcut){
          save.iso<-c(save.iso,k)
        }
    }
    
    #find out the primary isotope
    if (length(save.iso)>0){
      mz.iso<-mz.iso[save.iso]
      
      #'the smallest isotopic peak is the exact mass
      peaklist$PrimaryIso[temp]<-min(mz.iso)
    }
    }
  }
  
#' deleting redundant isotopic peaks
  index.save<-NULL
  for (i in 1:(nrow(peaklist)-1)){
    mz<-peaklist$PrimaryIso[i]
    rt<-peaklist$rt[i]
    for (j in (i+1):nrow(peaklist)){
      mz0<-peaklist$PrimaryIso[j]
      rt0<-peaklist$rt[j]
      ppm<-abs(mz-mz0)*1000000/mz
      RTdiff<-abs(rt-rt0)
      
      if (ppm<2.5&&RTdiff<0.5){
        index.save<-c(index.save,j)
      }
    }
  }
  
#'deleteing redundant isotopic peaks  
  if (length(index.save)>0){
    peaklist<-peaklist[-index.save,]
  }
  
  return(peaklist)
  }




#' This function is used to find the primary isotopic peaks
#'
#' @param mz1 
#' @param mz2 
#'
#' @return
#'
#' @examples
FindHomo<-function(mz1,mz2){
  delta.mz<-abs(mz1-mz2)
  for (x in 0:10){
    for (y in 0:10){
      for (z in 0:2){

#' considering C, Cl and Br isotopes        
        delta<-x*1.99795+y*1.99705+1.00335*z
        
#' <0.003,H and carbon
        if (abs(delta.mz-delta)<0.002){return(1)}      
      }
    }
  }
  return(0)
}