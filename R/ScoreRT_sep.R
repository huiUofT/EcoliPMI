library(readxl)
mydb<-read_excel("Metabolites_STD_KEGG.xlsx")
mydb$Score<-rep(0,nrow(mydb))
for (i in 1:nrow(mydb)){
  smiles.list<-mydb$SMILES[i]
  smiles.list<-strsplit(smiles.list,';')
  smiles.list<-unlist(smiles.list)

  temp.score<-NULL

  for (j in 1:length(smiles.list)){
    ##scoring MS
    ms.list<-mydb$mserror[i]
    ms.list<-strsplit(ms.list,';')
    ms.list<-unlist(ms.list)
    score.ms<--(as.numeric(ms.list[j]))^2/1^2#SD of MS error

    #score RT
    logP.list<-mydb$LOGP[i]
    logP.list<-strsplit(logP.list,';')
    logP.list<-unlist(logP.list)
    logP.list<-as.numeric(logP.list)
    RT.list<-0.473*logP.list[j]+2.57
    score.rt<--(RT.list-mydb$rt[i])^2/0.5^2

    allscore<-score.ms+score.rt
    if (length(temp.score)==0){
      temp.score<-allscore
    }else{
    temp.score<-paste(temp.score,allscore,sep= ';')}
  }
  mydb$Score[i]<-temp.score
}

####ranking db
mydb$ranking<-rep(1,nrow(mydb))
for (i in 1:nrow(mydb)){
  score.list<-mydb$Score[i]
  score.list<-strsplit(score.list,';')
  score.list<-unlist(score.list)
  score.list<-as.numeric(score.list)
  test<-rank(-score.list,ties.method='first')##ranking
  mydb$ranking[i]<-paste(test,collapse = ';')
}

write.table(mydb,file='Metabolites_KEGG_FDR.csv',sep=',',row.names = FALSE)

