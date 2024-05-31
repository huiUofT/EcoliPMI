#-------------------------------------------
#this function is used to extract MS sepctra from DIA windows
#-------------------------------------------
#'
#' @param Library
#' @param path
#' @param msfiles
#' @param mz_tol
#' @export
GetFrag<-function(Library,mz_tol,DIAisowin){
  #' the path to save results
  path<-getwd()

  #' the path to save raw data
  path.data<-paste0(path,"/data")
  setwd(path.data)
  msfiles<-list.files()

  fragments<-NULL
  for (k in 1:length(msfiles)){
    print(c('getfragment...',k))
    xraw<-xcmsRaw(msfiles[k],includeMSn=TRUE)

    index<-which(Library$SampleID==k)
    if (length(index)==0){next}

    #'extracting precursor DIA windows
    precursor<-preclist(xraw)
    len<-length(precursor)
    for (j in 1:length(index)){
      mz<-Library$mz[index[j]]
      DIAwin<-which(abs(mz-precursor)<=(0.5*DIAisowin))
      if (length(DIAwin)<1){next}
      DIAwin<-precursor[DIAwin[1]]

      #' get fragments from each window based correlation >0.85
      fragments<-getfrag(xraw,fragments,index[j],mz_tol,Library,DIAwin,30)
    }
  }
  setwd(path)
  return(fragments)
  }


#' -------------------
#' Finding DIA windows
#'--------------------
#' @param xmsn
#' @return
preclist<-function (xmsn){
  x<-xmsn

  #'extracting DIA windows
  precmz<-xmsn@msnPrecursorMz
  len<-length(precmz)
  precur<-precmz[1]
  for (i in 2:len){
    if (length(which(precur==precmz[i]))==0){
      precur<-c(precur, precmz[i])
    }
    }
  return(precur)
}


#' -------------------------------
#' This function is used to extract fragmetns from each DIA window
#' -------------------------------
#' @param rawdata
#' @param prec_list
#' @param index2
#' @param mz_tol
#' @param Library
#' @param DIAmzwin
#' @param rtwindow
#'
#' @return
#'
getfrag<-function(rawdata,prec_list,index2,mz_tol,Library,DIAmzwin,rtwindow){
  Library$rt<-Library$rt*60
  index<-index2
  precurmz<-Library$mz[index2]

  #' read rawdata for each DIA window
  DIAdata<-ms2copy(rawdata,DIAmzwin)
  mzrange<-DIAdata@mzrange
  minmz<-mzrange[1]
  maxmz<-mzrange[2]
  mzmin<-max(minmz,precurmz-precurmz*mz_tol)
  mzmax<-min(maxmz,precurmz+precurmz*mz_tol)
  rtrange<-DIAdata@scantime
  rtmin<-max(min(rtrange),Library$rt[index]-30)
  rtmax<-min(max(rtrange),Library$rt[index]+30)

  #'extracting precursor ions and peaks
  peaks<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))

  #'finding the scan number of the peak top
  scan.max<-which.max(peaks$intensity)
  scan.max<-peaks$scan[scan.max[1]]
  if (scan.max>=length(DIAdata@scanindex)){##the last scanning point
    return(prec_list)
  }
  print(scan.max)
  scanNum<-c(DIAdata@scanindex[scan.max],DIAdata@scanindex[scan.max+1])
  correctindex<-(scanNum[1]+1):scanNum[2]

  #'finding co-eluting ions
  mz.frag<-DIAdata@env$mz[correctindex]
  if (length(mz.frag)<1){return(prec_list)}
  index<-which(mz.frag<precurmz-10)
  if (length(index)<1){return(prec_list)}
  mz.frag<-mz.frag[index]

  #' using correlations to find fragments
  native.peak<-peaks$intensity
  kk<-0
  for (k in 1:length(mz.frag)){
    mz0<-mz.frag[k]
    mzmin<-max(minmz,mz0-mz0*mz_tol)
    mzmax<-min(maxmz,mz0+mz0*mz_tol)

    #'fragment peaks
    frag.peak<-rawEIC(DIAdata,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
    frag.peak<-frag.peak$intensity
    if (sd(native.peak)==0||sd(frag.peak)==0){next}
    if (max(frag.peak)<2000){next}

    #'correlation calculations
    corr<-cor(native.peak,frag.peak)

    #'save the fragments to list
    if (corr>0.8){
      list.frag<-list(precursor=precurmz,mz=mz.frag[k],rt=Library$rt[index2],intensity=max(frag.peak),score=corr,sampleid=Library$SampleID[index2],libraryid=index2)
      if (length(prec_list)<1){
        prec_list<-list.frag}
      else{
        prec_list$precursor<-c(prec_list$precursor,list.frag$precursor)
        prec_list$mz<-c(prec_list$mz,list.frag$mz)
        prec_list$rt<-c(prec_list$rt,list.frag$rt)
        prec_list$intensity<-c(prec_list$intensity,list.frag$intensity)
        prec_list$score<-c(prec_list$score,list.frag$score)
        prec_list$sampleid<-c(prec_list$sampleid,Library$SampleID[index2])
        prec_list$libraryid<-c(prec_list$libraryid,index2)
      }
    }
    }
  return (prec_list)
  }

#' -------------------------------------------
#' function to copy MS2 data to MS1 matrix, but without precursor information
#'--------------------------------------------
#' @param xmsn
#' @param precursor
#'
#' @return
#'
ms2copy <-function(xmsn,precursor) {
  x<-new("xcmsRaw")
  x@env<-new.env(parent = .GlobalEnv)
  index<-which(xmsn@msnPrecursorMz==precursor)
  x@tic <- xmsn@msnAcquisitionNum[index]
  x@scantime <- xmsn@msnRt[index]
  x@acquisitionNum <- xmsn@msnAcquisitionNum[index]
  x@polarity<-xmsn@polarity[1:length(index)]
  len2<-length(xmsn@msnPrecursorMz)
  index_total<-0
  index3<-0
  for (j in 1:length(index)){
    if (index[j]==len2){
      index2<-(xmsn@msnScanindex[index[j]]+1):length(xmsn@env$msnMz)
    }
    else{
      index2<-(xmsn@msnScanindex[index[j]]+1):xmsn@msnScanindex[index[j]+1]
    }
    index3<-c(index3,index3[length(index3)]+length(index2))
    index_total<-c(index_total,index2)
  }
  index_total<-index_total[-1]
  index3<-index3[-length(index3)]
  x@env$mz <- xmsn@env$msnMz[index_total]
  x@env$intensity <- xmsn@env$msnIntensity[index_total]
  x@mzrange<-c(min(x@env$mz),max(x@env$mz))
  x@scanindex <-as.integer(index3)
  return(x)
}
