#' This function is used to extract isotopic peaks and save them to a List
#'
#' @param Library.new 
#' @param mz_tol 
#'
#' @return
#' @export
IsotopeFind<-function(Library.new,mz_tol,corrcut){####delete the secondary isotopic peaks
  #' the path to save results
  path<-getwd()
  
  #' the path to save raw data
  path.data<-paste0(path,"/data")
  
  
  setwd(path.data)
  msfiles<-list.files()
  mz_tol<-mz_tol
  index.save<-NULL
  IDsave<-0
  peaklist<-Library.new
  peaklist$PrimaryIso<-peaklist$mz
  List.Isotope<-NULL
  
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
      dbID<-index[j]
      
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
      if (scan.max>=length(xrawdata@scanindex)){##the last scanning point
        next
      }
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
      
      native.peak<-peaks$intensity
      kk<-0
      for (k in 1:length(mz.iso)){
        mz0<-mz.iso[k]
        mzmin<-max(minmz,mz0-mz0*mz_tol)
        mzmax<-min(maxmz,mz0+mz0*mz_tol)
        isotope.peak<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        isotope.peak<-isotope.peak$intensity
        if (sd(isotope.peak)==0||sd(native.peak)==0){
          next
          }
        corr<-cor(native.peak,isotope.peak)
        if (corr>corrcut){
          List.Isotope$mz<-c(List.Isotope$mz,mz)
          List.Isotope$sampleID<-c(List.Isotope$sampleID,i)
          List.Isotope$Isotope<-c(List.Isotope$Isotope,mz0)
          List.Isotope$rt<-c(List.Isotope$rt,rt/60)
          List.Isotope$cor<-c(List.Isotope$cor,corr)
          List.Isotope$intensity.iso<-c(List.Isotope$intensity.iso,max(isotope.peak))
          List.Isotope$intensity<-c(List.Isotope$intensity,max(native.peak))
          List.Isotope$ID<-c(List.Isotope$ID,dbID)
        }
        }
      }
  }
  
  setwd(path)
  return(List.Isotope)
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