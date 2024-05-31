#' This function is used to find precurosr ion window from DIA
#'
#' @param xmsn 
#'
#' @return
#' @export
#'
#' @examples
Preclist<-function (xmsn){
  x<-xmsn
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