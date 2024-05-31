#-----------------------------------
#output chemicals with scores greater than cutoff
#-----------------------------------
#' This function is used to write hits greater than score cutoff
#'
#' @param Target 
#' @param cutoff 
#'
#' @return
#' @export
#'
#' @examples
Output<-function(Target,cutoff){
  index.del<-NULL
  Target$Final<-rep(0,nrow(Target))
  for (i in 1:nrow(Target)){
    score<-Target$allscore[i]
    temp<-strsplit(as.character(score),';')
    temp<-as.numeric(unlist(temp))
    index<-which(temp>cutoff)
    if (length(index)==0){
      index.del<-c(index.del,i)
      next
    }
    
    #'Only output the maximal formula
    index<-which.max(temp)
    formula<-strsplit(as.character(Target$formula.pred[i]),';')
    formula<-unlist(formula)
    Target$Final[i]<-paste(c(formula[index]),collapse = ';')
  }
  if (length(index.del)>0){
  output<-Target[-index.del,]}
  else{
    output<-Target
  }
  return(output)
}