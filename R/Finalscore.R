#------------------------------------------------
#calculate final score
#------------------------------------------------
#' This function is used to sum all scores
#'
#' @param mylib 
#' @param weightK 
#' @param precursor 
#'
#' @return
#' @export
#'
#' @examples
Finalscore<-function(mylib,weightK,precursor){
  mylib$allscore<-rep(0,nrow(mylib))
  for (i in 1:nrow(mylib)){
    formula<-mylib$formula[i]
    if (formula[1]=='0'){next}
    formula<-unlist(strsplit(formula,';'))
    
    MS1score<-unlist(strsplit(mylib$ms1score[i],';'))
    
    #‘adduct score
    if (mylib$adductscore[i]==0){
      adductscore<-rep(0,length(formula))
    }else{
      adductscore<-unlist(strsplit(mylib$adductscore[i],';'))
    }
    
    #’character score
    chascore<-0
    if (length(mylib$chracaterscore)>0){
    if (mylib$chracaterscore[i]==0){
      chascore<-rep(0,length(formula))
    }else{
      chascore<-unlist(strsplit(mylib$chracaterscore[i],';'))
    }}
    
    #‘character score
    neutralscore<-0
    if (length(mylib$neutralscore)>0){
    if (mylib$neutralscore[i]==0){
      neutralscore<-rep(0,length(formula))
    }else{
      neutralscore<-unlist(strsplit(mylib$neutralscore[i],';'))
      neutralscore<-neutralscore[-1]
    }}
    
    #’ion mode score
    if (mylib$ionmodescore[i]==0){
      ionmodescore<-rep(0,length(formula))
    }else{
      ionmodescore<-unlist(strsplit(mylib$ionmodescore[i],';'))
    }
    
    #'MS2score
    if (length(mylib$MS2score)>0){
    if (mylib$MS2score[i]==0){
      MS2score<-rep(0,length(formula))
    }else{
      MS2score<-unlist(strsplit(mylib$MS2score[i],';'))
    }}
    
    #'isotopic peaks score
    if (mylib$isoscore[i]==0){
      isoscore<-rep(0,length(formula))
    }else{
      isoscore<-unlist(strsplit(mylib$isoscore[i],';'))
    }
    
    MS1score<-as.numeric(MS1score)
    MS2score<-as.numeric(MS2score)
    ionmodescore<-as.numeric(ionmodescore)
    neutralscore<-as.numeric(neutralscore)
    chascore<-as.numeric(chascore)
    adductscore<-as.numeric(adductscore)
    
    temp.score<-MS1score*weightK[1]+MS2score*weightK[2]+ionmodescore*weightK[3]+
      neutralscore*weightK[4]+chascore*weightK[5]+adductscore*weightK[6]
    
    #‘ion mode score is -1, put all score to 0
    index.ion<-which(ionmodescore<0)
    if (length(index.ion)>0){
      temp.score[index.ion]<-0
    }
    
    #’isotope score is <0.8, put all score to 0
    isoscore<-as.numeric(isoscore)
    index.iso<-which(isoscore<(log(0.5)))
    if (length(index.iso)>0){
      temp.score[index.iso]<-0
    }
    
    mylib$allscore[i]<-paste(c(temp.score),collapse=';')
  }
  
  index<-which(mylib$mz>min(precursor)-2.5)
  index2<-which(mylib$mz[index]<max(precursor)+2.5)
  
  #’only output the results with MS2 fragments
  #mylib<-mylib[index[index2],]
  return(mylib)
}
