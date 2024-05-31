#----------------------------------------------------
#this function is used to find adducts for detected features
#----------------------------------------------------
#' @param Library.new
#' @param mw.adducts
#' @param mz_tol
#' @param Adducts
#'
#' @return
#' @export
AdductsFind<-function(Library.new,mw.adducts,mz_tol,Adducts){
  #' the path to save results
  path<-getwd()

  #' the path to save raw data
  path.data<-paste0(path,"/data")

  setwd(path.data)
  msfiles<-list.files()
  index.save<-NULL
  IDsave<-0
  List.Isotope<-NULL
  for (i in 1:length(msfiles)){
    print(c('MS file...',i))
    xrawdata<-xcmsRaw(msfiles[i])
    index<-which(Library.new$SampleID==i)
    if (length(index)<1){next}

    #finding coeluting compound from MS1
    for (j in 1:length(index)){
      temp<-index[j]
      mz<-Library.new$mz[temp]
      rt<-Library.new$rt[temp]*60
      mzrange<-xrawdata@mzrange
      minmz<-mzrange[1]
      maxmz<-mzrange[2]
      mzmin<-max(minmz,mz-mz*mz_tol)
      mzmax<-min(maxmz,mz+mz*mz_tol)
      rtrange<-xrawdata@scantime
      rtmin<-max(min(rtrange),rt-30)
      rtmax<-min(max(rtrange),rt+30)
      if (rtmax<0){next}
      peaks<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
      scan.max<-which.max(peaks$intensity)
      scan.max<-peaks$scan[scan.max[1]]
      if (scan.max>=length(xrawdata@scanindex)){##the last scanning point
        next
      }
      scanNum<-c(xrawdata@scanindex[scan.max],xrawdata@scanindex[scan.max+1])
      correctindex<-(scanNum[1]+1):scanNum[2]
      mz.iso<-xrawdata@env$mz[correctindex]
      if (length(mz.iso)<1){next}


      #'Finding adducts
      save.iso<-data.frame()
      mzadducts<-NULL
      adductssave<-NULL

      #'calculating mass difference between adducts
      for (mm in 1:length(mz.iso)){
        for (kk in 2:length(mw.adducts)){
          mzdiff<-mw.adducts[kk]-mw.adducts[2]
          if (abs(mz-mz.iso[mm]-mzdiff)<mz*mz_tol){
            mzadducts<-c(mzadducts,mz.iso[mm])
            adductssave<-c(adductssave,kk)
          }
        }
        }
      if (length(mzadducts)<1){next}
      save.iso<-data.frame(mz=mzadducts,adducts=Adducts[adductssave])
      mz.iso<-save.iso$mz

      #'corrlations to further confirm adducts
      native.peak<-peaks$intensity
      kk<-0
      for (k in 1:length(mz.iso)){
        mz0<-mz.iso[k]
        mzmin<-max(minmz,mz0-mz0*mz_tol)
        mzmax<-min(maxmz,mz0+mz0*mz_tol)
        isotope.peak<-rawEIC(xrawdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
        isotope.peak<-isotope.peak$intensity
        if (sd(isotope.peak)==0||sd(native.peak)==0){next}###0 values
        corr<-cor(native.peak,isotope.peak)
        if (corr>0.80){
          List.Isotope$mzadducts<-c(List.Isotope$mzadducts,mz)
          List.Isotope$sampleID<-c(List.Isotope$sampleID,i)
          List.Isotope$mzMH<-c(List.Isotope$mzMH,mz0)
          List.Isotope$rt<-c(List.Isotope$rt,rt/60)
          List.Isotope$cor<-c(List.Isotope$cor,corr)
          List.Isotope$intensityMH<-c(List.Isotope$intensityMH,max(isotope.peak))
          List.Isotope$intensityadduct<-c(List.Isotope$intensityadduct,max(native.peak))
          List.Isotope$ID<-c(List.Isotope$ID,temp)
          List.Isotope$Adducts<-c(List.Isotope$Adducts,as.character(save.iso$adducts[k]))
        }
      }
    }
  }

  setwd(path)
  return(List.Isotope)
  }
