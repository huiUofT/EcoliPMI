#------------------------------------------------------------------
#This function is used for final database searching wiht different scores calculated
#------------------------------------------------------------------
#' @param cutoff
#' @param polarity
#' @param Library.new
#' @param Isotope.Data
#' @param Cal.Frag
#' @param iso_list
#' @param Database
#'
#' @import isopat
#' @return
#' @export
DatabaseSearching<-function(Library.new, cutoff,polarity,Database,Isotope.Data,Cal.Frag,iso_list,ms2input){
  element<-c("C","H","N","O","P","S","35Cl","37Cl","79Br","81Br","F")
  number_input<-t(rbind(c(1,0,0,0,0,0,0,0,0,0,0),c(50,80,10,10,5,5,5,5,5,5,30)))
  rownames(number_input)<-element
  isotope<-1.9978#'the value of isotope peaks
  inputvalue<-"yes"##for some compounds we don't need isotope information
  peakwidth<-10#peak width, seconds
  samplenum<-21##the number of sample
  cutratio<-5###cutoff for ratio difference
  timewin<-12###12 sec for rt cutoff
  mw.element<-iso.element(element,iso_list,1)

  #'calculating isotopic peak patterns
  IsotopeDB<-NULL

  #' the path to save results
  path<-getwd()

  #' the path to save raw data
  path.data<-paste0(path,"/data")

  #'save all information to a new data matrix
  setwd(path.data)

  for (i in 1:nrow(Database)){
    print(c('isotope...',i))
    formula<-Database$formula[i]

    #'must contain carbon element
    carbon.index<-grep('C',formula)
    if (length(carbon.index)==0){next}

    formula.split<-form.split(formula)
    formula.cal<-form.comb(formula.split)
    if (max(formula.split)==0) {next}
    iso.pattern<-isopattern_cal(iso_list,formula.cal,1e-4)
    iso.pattern[,2]<-iso.pattern[,2]/max(iso.pattern[,2])

    #'only keeping higher abundant isotopic peaks
    index<-which(iso.pattern[,2]>0.005)
    if (length(index)>1){
      IsotopeDB$mz<-c(IsotopeDB$mz,iso.pattern[index,1])
      IsotopeDB$ID<-c(IsotopeDB$ID,rep(i,length(index)))
    }
    }

  msfiles<-list.files()
  mylib<-Library.new
  mylib$isoscore<-rep(0,nrow(Library.new))
  mylib$ms1score<-rep(0,nrow(Library.new))

  #'saving information
  for (i in 1:length(msfiles)){
    print(c('msfiles...',i))
    index<-which(Library.new$SampleID==i)
    if (length(index)==0){next}
    xrawdata<-xcmsRaw(msfiles[i])
    for (j in 1:length(index)){
      formula<-NULL

      #'determining mass errors
      temp<-unlist(strsplit(Library.new$formula[index[j]],';'))
      for (k in 1:length(temp)){
      formula<-formcal(polarity,number_input,Library.new,index[j],k,iso_list,xrawdata,3,3,Isotope.Data,Database,Cal.Frag,Adducts,IsotopeDB,0)
      if (length(formula)<2){next}
      max.form<-nrow(formula)
      mylib$isoscore[index[j]]<-paste(formula[,6],collapse=';')
      mylib$ms1score[index[j]]<-paste(formula[,8],collapse=';')
      }
    }
  }

  ##-------------------------------------------------
  #neutral loss,Anal. Chem. 2014, 10724
  #-------------------------------------------------
  library(ChemmineOB)
  if (polarity==1){
    Neutral<-data.frame(frag=c('NH3','H2O','CH2O','CH3OH','H2S','HCl','NO2','HCOOH','H3PO4','CO2','CO'),
                        mz=c(17.0265,18.0106,30.0106,32.0262,33.9877,35.9767,45.9929,46.0055,97.97689,43.9898,27.9949))

  }
  if (polarity==-1){
    Neutral<-data.frame(frag=c('CH2O','H2O','NO2','HCOOH','CO2','CO','HBr','Br','HCl'),
                        mz=c(30.0106,18.0106,45.9929,46.0055,43.9898,27.9949,79.9262,78.9183,35.9767))
  }
  mylib.neutral<-Score.neutral(mylib,Cal.Frag,Neutral)
  if (ms2input=='ms2false'){
    mylib.neutral<-mylib}

  #---------------------------------------------
  #MS2 spectra
  #------------------------------------------------
  if (polarity==1){
    Character.ion<-data.frame(frag=c('NH3'),
                              mz=c(17.0265))

  }
  if (polarity==-1){
    Character.ion<-data.frame(frag=c('Br','Br','SO3','PO3
                                     '),
                              mz=c(78.91885,80.91685,79.95735,78.95905))
  }
  mylib.charac<-Score.charac(mylib.neutral,Cal.Frag,Neutral)
  if (ms2input=='ms2false'){
    mylib.charac<-mylib}

  #------------------------------------------------
  #ion mode
  #------------------------------------------------
  mylib.ionmode<-mylib
  mylib.ionmode$ionmodescore<-rep(0,nrow(mylib.ionmode))
  for (i in 1:nrow(mylib.charac)){
    print(c('ionmode',i))
    smiles<-mylib.ionmode$SMILES[i]##SMILES
    if (smiles==0){next}##no library compound
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    adduct<-mylib.ionmode$Adducts[i]##adducts
    if (adduct==0){next}##no library compound
    adduct<-strsplit(adduct,';')
    adduct<-adduct[[1]]
    score.save<-rep(0,length(smiles))
    for (j in 1:length(smiles)){
      score.save[j]<-pattern.ionmode(smiles[j],polarity,'ESI',adduct[j])/2##score
    }
    mylib.ionmode$ionmodescore[i]<-paste(score.save,collapse=';')
  }
  #-------------------------------
  #score adducts
  #------------------------------
  myfinallib<-Score.Adducts(mylib.ionmode,Adducts.Find)##if the adduct is matched to the formula, plus 1

  #-------------------------------------
  #MS2 score
  #-----------------------------------------
  if (ms2input=='ms2false'){
     myfinallib$MS2score<-rep(0,nrow(myfinallib))
  }else{
  myfinallib<-Score.ms2(myfinallib,Cal.Frag,Isotope.Data,Database,polarity)
  }##if the adduct is matched to the formula, plus 1

  return(myfinallib)
}

#' -------------------
#' function to group isotope value for selected element
#' -------------------
#' @param element
#' @param iso_list
#' @param input
#'
#' @return

iso.element<-function(element,iso_list,input){
  index<-NULL
  for (i in 1:length(element)){
    indextemp<-which(iso_list[,input]==element[i])
    if (length(indextemp)<1){
      indextemp<-which(iso_list[,2]==element[i])}
    indextemp<-min(indextemp)
    index<-c(index,indextemp)
    }
  return(iso_list[index,])
}

#‘------------------------------------------------
#’Neutral Loss score
#‘------------------------------------------------
#' @param Library
#' @param fragment
#' @param Neutral
#'
#' @return
Score.neutral<-function(Library,fragment,Neutral){
  #’Determine which compounds show neutral loss
  Neutralloss<-LossCal(Library,Neutral,fragment)

  #‘No neutral loss
  if (length(Neutralloss)==0){
    return(NULL)
    }
  Library$neutralscore<-rep(0,nrow(Library))
  for (i in 1:length(Neutralloss$mz)){
    smiles<-Library$SMILES[Neutralloss$libid[i]]

    #’No library compound
    if (smiles==0){
      next
      }
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    score.save<-NULL

    #‘ The formula for neutral loss
    Loss.pattern<-Neutralloss$fragment[i]
    for (k in 1:length(smiles)){

    #'0.5 for each neutral loss
      score.temp<-pattern.neutral(smiles[k],Loss.pattern)/2
      score.save<-paste(score.save,score.temp,sep=';')
    }
    Library$neutralscore[Neutralloss$libid[i]]<-score.save
    }
  return(Library)
}

#'-------------------------------------
#'neutral loss ID
#'------------------------------------
#'
#' @param Neutral
#' @param Fragment
LossCal<-function(Library,Neutral,Fragment){
  loss.save<-NULL
  for (i in 1:nrow(Library)){
    index<-which(Fragment$libraryid==i)

    #'No fragment
    if (length(index)==0){
      next
      }
    mz.prec<-Fragment$precursor[index[1]]
    mz.frag<-Fragment$mz[index]
    loss.mz<-mz.prec-mz.frag
    if (length(mz.frag)>1){
      for (j in 1:(length(mz.frag)-1)){
        for (m in (j+1):length(mz.frag)){
          loss.mz<-c(loss.mz,abs(mz.frag[m]-mz.frag[j]))
        }
      }
      }
    for (k in 1:length(loss.mz)){
      #'Identifying loss
      index2<-which(abs(loss.mz[k]-Neutral$mz)<0.002)
      if (length(index2)>0){
        loss.save$mz<-c(loss.save$mz,Neutral$mz[index2])
        loss.save$libid<-c(loss.save$libid,i)
        loss.save$fragment<-c(loss.save$fragment,as.character(Neutral$frag[index2]))
      }
    }
    }
  return(loss.save)
}

#'------------------------------------------------
#'Neutral Loss, and characteristic pattern seaching against SMARTS
#'------------------------------------------------
#'
#' @param smiles
#' @param fragment
pattern.neutral<-function(smiles,fragment){
  sdf<-smiles2sdf(smiles)
  sdf<-as(sdf,'SDFset')
  if (fragment=='CO2'){
    Match.loss= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)##carboxylic acid or conjugate salt

    #' matched
    if (max(Match.loss)>0){return(1)
      }
  }
  if (fragment=='H2O'){
    #' Carboxylic acid or conjugate salt
    Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)

    #' Ester or aldehyde
    Match.loss2= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)
    Match.loss<-c(Match.loss1,Match.loss2)

    #' Matched
    if (max(Match.loss)>0){
      return(1)
      }
  }
  if (fragment=='NH3'){
    #'primary amine
    Match.loss= smartsSearchOB(sdf,"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='CH2O'){
    #'Aldehyde
    Match.loss= smartsSearchOB(sdf,"[CX3H1](=O)[#6]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='CH3OH'){
    #'ester or aldehyde
    Match.loss= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='H2S'){
    #'THIOL
    Match.loss= smartsSearchOB(sdf,"[#16X2H]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='HCl'){
    #'Cl
    Match.loss= smartsSearchOB(sdf,"[#6][Cl]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='NO2'){
    #'NITRO GROUP
    Match.loss= smartsSearchOB(sdf,"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='HCOOH'||fragment=='PO3'){
    #'carboxylic acid or conjugate salt
    Match.loss= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)
    #'matched
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='H3PO4'){
    #'Phosphate
    Match.loss= smartsSearchOB(sdf,"[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='CO'){
    #'Carboxilic acid, ester,anhydride and ketone
    Match.loss= smartsSearchOB(sdf,"[CX3](=[OX1])C",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='HBr'||fragment=='Br'){
    #'Cl
    Match.loss= smartsSearchOB(sdf,"[#6][Br]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  if (fragment=='SO3'){
    #'Sulfonic acid and conjugate
    Match.loss= smartsSearchOB(sdf,"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",uniqueMatches=FALSE)
    if (max(Match.loss)>0){return(1)}
  }
  return (0)
}

#------------------------------------------------
#Characteristic ions
#------------------------------------------------
#'
#'
#' @param Library
#' @param fragment
#' @param Character.ion
#'
#' @return
Score.charac<-function(Library,fragment,Character.ion){
  index.match<-NULL
  ions.match<-NULL

  #'determine if there is characteristic ions
  for (i in 1:length(fragment$mz)){
    index<-which(abs(Character.ion$mz-fragment$mz[i])<0.002)
    if (length(index)>0){
      index.match<-c(index.match,i)
      ions.match<-c(ions.match,as.character(Character.ion$frag[index[1]]))}
  }
  Library$chracaterscore<-rep(0,nrow(Library))
  if (length(index.match)==0){return(Library)}
  for (i in 1:length(index.match)){
    index.lib<-fragment$libraryid[index.match[i]]

    #'smiles of compounds
    smiles<-Library$SMILES[index.lib]
    if (smiles==0){next}
    smiles<-strsplit(smiles,';')
    smiles<-smiles[[1]]
    score.save<-NULL

    #'Formulas of fragments
    Loss.pattern<-ions.match[i]
    for (k in 1:length(smiles)){

      #'0.5 for each characteristic fragment
      score.temp<-pattern.neutral(smiles[k],Loss.pattern)/2
      if(is.null(score.save)){
        score.save<-score.temp
        next
      }
      score.save<-paste(score.save,score.temp,sep=';')
    }
    Library$chracaterscore[index.lib]<-score.save}
  return(Library)
}

#------------------------------------------------
#ion mode pattern match
#------------------------------------------------
#'
#'
#' @param smiles
#' @param Ionmode
#' @param Ionsource
#' @param Adducts
#'
#' @return
pattern.ionmode<-function(smiles,Ionmode,Ionsource,Adducts){
  sdf<-smiles2sdf(smiles)
  sdf<-as(sdf,'SDFset')
  if (Ionmode==-1){
    if (Adducts=='[M-H2O-H]-'){

      #'carboxylic acid or conjugate salt
      Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)

      #'phenol and hydroxyl, acidic group including sulfonic and phosphoric acids
      Match.loss2= smartsSearchOB(sdf,"[OX2H]",uniqueMatches=FALSE)
      Match.loss<-c(Match.loss1,Match.loss2)
      if (max(Match.loss)==0){return(-1)
        }
      return (0)
    }
    if (Adducts=='[M+Cl]-'){
      #'chloride,M+Cl
      Match.loss<-smartsSearchOB(sdf,"[#6][Cl]",uniqueMatches=FALSE)
      #'matched
      if (max(Match.loss)==0){
        return(-1)
        }
      return (0)
    }
    if (Adducts=='[M-H]-'){
      #'carboxylic acid or conjugate salt
      Match.loss1= smartsSearchOB(sdf,"[CX3](=O)[OX1H0-,OX2H1]",uniqueMatches=FALSE)

      #'aldehyde
      Match.loss2= smartsSearchOB(sdf,"[CX3H1](=O)[#6]",uniqueMatches=FALSE)

      #'phenol and hydroxyl, acidic group including sulfonic and phosphoric acids
      Match.loss3= smartsSearchOB(sdf,"[OX2H]",uniqueMatches=FALSE)

      #'amine
      Match.loss4= smartsSearchOB(sdf,"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",uniqueMatches=FALSE)

      #'thiol
      Match.loss5= smartsSearchOB(sdf,"[#16X2H]",uniqueMatches=FALSE)

      #'hydroxy acidic
      Match.loss6= smartsSearchOB(sdf,"[$([OH]-*=[!#6])]",uniqueMatches=FALSE)

      #'phosphate acid
      Match.loss7= smartsSearchOB(sdf,"[$([P]=O)][OX2H1]",uniqueMatches=FALSE)
      Match.loss<-c(Match.loss1,Match.loss2,Match.loss3,Match.loss4,Match.loss5,Match.loss6,Match.loss7)
      if (max(Match.loss)==0){return(-1)}##matched
      return (0)
    }
    if (Ionsource=='APCI'){
      #'nitro group
      Match.loss1<-smartsSearchOB(sdf,"[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",uniqueMatches=FALSE)

      #'nitro group
      Match.loss3<-smartsSearchOB(sdf,"[NX4](=O)",uniqueMatches=FALSE)
      Match.loss2<-smartsSearchOB(sdf,"[#6][Br]",uniqueMatches=FALSE)
      if (Adducts=='[M]-'){
        Match.loss<-c(Match.loss1,Match.loss3)
        if (max(Match.loss)==0){return(-1)}
        return (0)}
      if (Adducts=='[M-Br+O]-'){
        Match.loss<-Match.loss2
        if (max(Match.loss)==0){return(-1)}
        return (0)
        }
    }
  }
  if (Ionmode==1){
    #'amine chemicals, but not amide
    Match.loss1= smartsSearchOB(sdf,"[NX3;H2,H1,H0;!$(NC=O)]",uniqueMatches=FALSE)

    #'ester and aldehyde
    Match.loss2= smartsSearchOB(sdf,"[#6][CX3](=O)[OX2H0][#6]",uniqueMatches=FALSE)

    #'ketone
    Match.loss3= smartsSearchOB(sdf,"[#6][CX3](=O)[#6]",uniqueMatches=FALSE)

    #'hydroxyl in alcohol
    Match.loss4= smartsSearchOB(sdf,"[#6][OX2H]",uniqueMatches=FALSE)

    #'phosphate ester
    Match.loss5=smartsSearchOB(sdf,"[$([P]=O)]([OX2][#6])([OX2][#6])([OX2][#6])",uniqueMatches=FALSE)
    Match.loss<-c(Match.loss1,Match.loss2,Match.loss3,Match.loss4,Match.loss5)
    if (max(Match.loss)==0){return(-1)}
  }
  return (0)
}


#------------------------------------------------------------------
#score the adducts
#------------------------------------------------------------------
#'
#'
#' @param mylib
#' @param Adducts.Find
#'
#' @return
Score.Adducts<-function(mylib,Adducts.Find){
  mylib$adductscore<-rep(0,nrow(mylib))
  for (i in 1:nrow(mylib)){
    index<-which(Adducts.Find$ID==i)
    if(length(index)==0){next}
    if (mylib$Adducts[i]==0){next}
    Adducts.temp<-strsplit(mylib$Adducts[i],';')
    adduct<-Adducts.temp[[1]]
    score.save<-rep(0,length(adduct))
    for (k in 1:length(adduct)){
      if(adduct[k]==Adducts.Find$Adducts[index[1]]){
        score.save[k]<-1
      }

      if (adduct[k]=='[M-Br+O]-'||adduct[k]=='[M]-'){
        score.save[k]<-0.5
      }
    }
    mylib$adductscore[i]<-paste(score.save,collapse = ';')
  }
  return(mylib)
}


#------------------------------------------------------------------
#score the MS fragmentation
#------------------------------------------------------------------
#'
#' @param mylib
#' @param Cal.Frag
#' @param Isotope.Data
#' @param mydb
#' @param polarity
#'
#' @return

Score.ms2<-function(mylib,Cal.Frag,Isotope.Data,mydb,polarity){#using dot product
  #'if fragment is not used
  mylib$MS2score<-rep(0,nrow(mylib))
  for (i in 1:nrow(mylib)){
    if (mylib$formula[i]==0){next}
    fragment.id<-which(Cal.Frag$libraryid==i)
    if (length(fragment.id)==0){next}
    mz<-mylib$mz[i]
    intensity<-10000
    Isotope.id<-which(Isotope.Data$ID==i)
    if (length(Isotope.id)>0){
      mz.isotope<-Isotope.Data$Isotope[Isotope.id]
      intensity.isotope<-Isotope.Data$intensity.iso[Isotope.id]
      intensity<-Isotope.Data$intensity[Isotope.id[1]]
    }
    ms1peak<-data.frame(ms=c(mz,mz.isotope),intensity=c(intensity,intensity.isotope))
    ms2peak<-ms1peak
    if (length(fragment.id)>0){
      mz.frag<-Cal.Frag$mz[fragment.id]
      intensity.frag<-Cal.Frag$intensity[fragment.id]
      ms2Act<-data.frame(ms=mz.frag,intensity=intensity.frag)
    }
    ms2Act$intensity<-ms2Act$intensity*100/max(ms2Act$intensity)##normalization

    ID.temp<-unlist(strsplit(mylib$SMILES[i],';'))
    scoresave<-NULL
    for (kk in 1:length(ID.temp)){#preparing reference ms2 spectra from database
      ID<-which(mydb$SMILES==ID.temp[kk])
      if (length(ID)<1){next}
      ID<-ID[1]
      if(polarity==1){
        Fragment_silico<-unlist(strsplit(mydb$Fragpos[ID],';'))
        ms2lib<-data.frame(ms=rep(0,length(Fragment_silico)),intensity=rep(0,length(Fragment_silico)))
        for (k in 1:nrow(ms2lib)){
          mzintensity<-unlist(strsplit(Fragment_silico[k],','))
          ms2lib[k,1]<-as.numeric(mzintensity[1])
          ms2lib[k,2]<-as.numeric(mzintensity[2])
        }
      }
      if(polarity==-1){#preparing reference ms2 spectra from database
        Fragment_silico<-unlist(strsplit(mydb$Fragneg[ID],';'))
        ms2lib<-data.frame(ms=rep(0,length(Fragment_silico)),intensity=rep(0,length(Fragment_silico)))
        for (k in 1:nrow(ms2lib)){
          mzintensity<-unlist(strsplit(Fragment_silico[k],','))
          ms2lib[k,1]<-as.numeric(mzintensity[1])
          ms2lib[k,2]<-as.numeric(mzintensity[2])
        }
      }

      ms2mz<-c(ms2Act$ms,ms2lib$ms)
      ms2Actintensity<-rep(0,length(ms2mz))
      ms2Actintensity[1:nrow(ms2Act)]<-ms2Act$intensity
      ms2libintensity<-rep(0,length(ms2mz))
      ms2libintensity[(nrow(ms2Act)+1):length(ms2mz)]<-ms2lib$intensity
      ppmms2<-15*10^(-6)
      index.match<-NULL
      for (mm in 1:nrow(ms2Act)){
        index<-which(abs(ms2lib$ms-ms2mz[mm])<(ms2mz[mm]*ppmms2))
        if (length(index)>0){
          index.match<-c(index.match,index)
          index2<-which.max(ms2lib$intensity[index])
          if (length(index2)>0){
          ms2libintensity[mm]<-ms2lib$intensity[index[index2]]
          }
        }
      }
      index.match<-index.match+nrow(ms2Act)
      ms2Actlib<-data.frame(mz=ms2mz,Act=ms2Actintensity,Ref=ms2libintensity)

      if (length(index.match)>0&&length(index.match)<nrow(ms2lib)){
        ms2Actlib<-ms2Actlib[-index.match,]
      }

      w_act <- 1/(1+(ms2Actlib[,2]/(sum(ms2Actlib[,2])-0.5)))*ms2Actlib[,2]
      w_lib <- 1/(1+(ms2Actlib[,3]/(sum(ms2Actlib[,3])-0.5)))*ms2Actlib[,3]
      ms2DotProduct <- (w_act %*% w_lib)^2/((w_act %*% w_act) * (w_lib %*% w_lib)) # dot product


      #for ( j in 1:length(Fragment_silico)){##using a simple matching to calculate ms2 score
      #  for (k in 1:nrow(ms2peak)){
      #    mserror<-(ms2peak$ms[k]-as.numeric(Fragment_silico[j]))/ms2peak$ms[k]
      #    if (abs(mserror)<ppm.ms2*10^(-6)){
      #
      #      #'0.5 for each fragment
      #      ms2score<-ms2score+0.5}
      #  }
      #}
      scoresave<-c(scoresave,ms2DotProduct)
    }
    mylib$MS2score[i]<-paste(scoresave,collapse = ';')
  }
  return(mylib)
}


#-----------------------
#split formula to matrix
#-----------------------
#' Title
#'
#' @param formula
#'
#' @return

form.split<-function(formula){
  formula.paste<-rep(0,11)
  element<-c('C','H','N','O','P','S','Cl','37Cl','Br','81Br','F')
  for (i in 1:length(element)){
    if (length(grep(element[i],formula))<1){next}##no element
    temp<-strsplit(formula,element[i])
    temp<-temp[[1]]
    if (length(temp)==1){
      formula.paste[i]<-1##the end of the formula
      next}
    temp<-temp[2]
    temp<-strsplit(temp,'')
    temp<-temp[[1]]
    number.ele<-grep('[0-9]',temp)
    if (length(number.ele)<1){##all 1 atom
      formula.paste[i]<-1
      next
    }
    if (number.ele[1]>1){##1
      formula.paste[i]<-1
      next
    }
    temp.save<-number.ele[1]
    if (length(number.ele)==1){
      formula.paste[i]<-temp[temp.save]
      next}
    for (k in 2:length(number.ele)){
      if (number.ele[k]>(number.ele[k-1]+1)){break}
      temp.save<-c(temp.save,number.ele[k])}
    formula.paste[i]<-paste(temp[temp.save],collapse='')
    }
  return(as.numeric(formula.paste))
}

###########function to read formulae
#' Title
#'
#' @param formulainput
#'
#' @return

form.comb<-function(formulainput){
  formulainput[7]<-sum(formulainput[7:8])
  formulainput[8]<-sum(formulainput[9:10])
  formulainput[9]<-formulainput[11]
  charinput<-c("C","H","N","O","P","S","Cl","Br","F")
  formulainput<-as.character(formulainput)
  form.output<-NULL
  for (i in 1:length(charinput)){
    if (formulainput[i]!=0){
      form.output<-paste(form.output,charinput[i],formulainput[i],sep="")
    }
    }
  return(form.output)
}


#' Title
#'
#' @param iso_list
#' @param compound
#' @param limit
#'
#' @return
#' @export
isopattern_cal<-function(iso_list, compound, limit){
  isos <- iso_list
  isos <- isos[isos[, 4] != 0, ]
  getit <- seq(1, length(isos[, 1]), 1)
  getthat <- c()
  element <- c()
  number <- c()
  ende <- nchar(compound)
  i <- c(1)
  while (i <= ende) {
    warn <- TRUE
    if (substr(compound, i, i) == c("[")) {
      k <- i
      while (any(substr(compound, i, i) == c("]")) != TRUE) {
        i <- c(i + 1)}
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound, i, i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) != TRUE) {
      k <- i
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- c(i - 1)
      element <- c(element, substr(compound, k, m))
      getthat <- c(getthat, getit[isos[, 1] == substr(compound,
                                                      k, m)])
      warn <- FALSE}
    if (any(substr(compound, i, i) == c("0", "1", "2", "3",
                                        "4", "5", "6", "7", "8", "9")) == TRUE) {
      k <- i
      while (any(substr(compound, i, i) == c("0", "1",
                                             "2", "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
        i <- c(i + 1)}
      m <- c(i - 1)
      i <- i - 1
      number <- c(number, as.numeric(substr(compound, k,
                                            m)))
      warn <- FALSE}
    if (warn == TRUE) {
      stop("Calculation interrupted: compound chemical formula incorrect!")
    }
    i <- i + 1}
  for (i in 1:length(element)) {
    if (any(isos[, 1] == element[i]) == FALSE) {
      stop("Calculation interrupted: one element not found in iso_list")}}
  isos <- t(isos[getthat, ])
  monomass <- as.double(0)
  monoabund <- as.double()
  for (i in 1:length(element)) {
    getit <- isos[, isos[1, ] == element[i]]
    if (length(dim(getit)) > 0) {
      monomass <- c(monomass + (as.double(getit[3, as.numeric(getit[4,
                                                                    ]) == max(as.numeric(getit[4, ]))]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4, as.numeric(getit[4,
                                                                     ]) == max(as.numeric(getit[4, ]))])^as.numeric(number[i])))}
    else {
      monomass <- c(monomass + (as.double(getit[3]) * as.numeric(number[i])))
      monoabund <- c(monoabund, (as.double(getit[4])^as.numeric(number[i])))}}
  if (monomass == 0) {
    stop("Calculation interrupted: monoisotopic mass could not be calculated!")}
  if (length(monoabund) == 0) {
    stop("Calculation interrupted: monoisotopic abundance incorrect!")}
  if (length(monoabund) != length(number)) {
    stop("Calculation interrupted: not all elements found in iso_list")}
  getit <- seq(1, length(isos[1, ]), 1)
  road <- matrix(nrow = length(element), ncol = 10, 0)
  rownames(road) = element
  colnames(road) = c("monoiso", rep("nonmono", 9))
  for (i in 1:length(element)) {
    getthat <- getit[isos[1, ] == element[i]]
    road[i, 1] = getthat[as.numeric(isos[4, isos[1, ] ==
                                           element[i]]) == max(as.numeric(isos[4, isos[1, ] ==
                                                                                 element[i]]))]
    getthat <- getthat[as.numeric(isos[4, isos[1, ] == element[i]]) !=
                         max(as.numeric(isos[4, isos[1, ] == element[i]]))]
    if (length(getthat) > 0) {
      road[i, c(2:(1 + length(getthat)))] = getthat}}
  leng1 <- length(isos[1, ])
  peaks <- matrix(ncol = (4 + leng1), nrow = 5e+03, as.double(0))
  leng2 <- length(peaks[1, ])
  leng3 <- length(peaks[, 1])
  peaks[1, 1] = monomass
  peaks[1, 2] = prod(monoabund)
  peaks[1, 3] = 0
  peaks[1, 4] = 0
  peaks[1, 4 + road[, 1]] = number
  colnames(peaks) = c("mass", "abundance", "generation", "stop?",
                      isos[2, ])
  if (length(monoabund) == 1 && monoabund == 1) {
    peaks[1, 4] = 1}
  k <- c(1)
  m <- c(2)
  b_start <- c(2)
  b_end <- c(1)
  g <- c(1)
  minlimit <- limit
  isonumber <- c()
  getit <- rep(1, length(road[1, ]))
  for (i in 1:length(road[, i])) {
    isonumber <- c(isonumber, length(getit[road[i, ] != 0]))}
  for (i in 1:length(road[, 1])) {
    if (isonumber[i] > 1) {
      for (j in 2:isonumber[i]) {
        peaks[m, c(5:(4 + leng1))] = peaks[k, c(5:(4 +
                                                     leng1))]
        peaks[m, 3] = g
        peaks[m, 1] = c(peaks[k, 1] - (as.numeric(isos[3,
                                                       road[i, 1]])) + (as.numeric(isos[3, road[i,
                                                                                                j]])))
        peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                      road[i, j]])/(peaks[m, 4 + road[i, j]] + 1) *
                          peaks[m, 4 + road[i, 1]]/as.numeric(isos[4,
                                                                   road[i, 1]]))
        peaks[m, 4 + road[i, 1]] = c(peaks[m, 4 + road[i,
                                                       1]] - 1)
        peaks[m, 4 + road[i, j]] = c(peaks[m, 4 + road[i,
                                                       j]] + 1)
        if (peaks[m, 2] < minlimit) {
          peaks[m, 4] = 1}
        m <- m + 1}}}
  g <- c(g + 1)
  k <- c(k + 1)
  b_end <- c(m - 1)
  while (any(peaks[peaks[, 3] == (g - 1), 4] == 0) && m < leng3 &&
         b_end >= b_start) {
    while (k <= b_end) {
      if (peaks[k, 2] >= limit) {
        for (i in 1:length(road[, 1])) {
          if (isonumber[i] > 1) {
            if (peaks[k, 4 + road[i, 1]] != 0 && m <
                leng3) {
              for (j in 2:isonumber[i]) {
                peaks[m, c(5:(4 + leng1))] = peaks[k,
                                                   c(5:(4 + leng1))]
                peaks[m, 3] = g
                peaks[m, 2] = c(peaks[k, 2] * as.numeric(isos[4,
                                                              road[i, j]])/(peaks[m, 4 + road[i,
                                                                                              j]] + 1) * peaks[m, 4 + road[i, 1]]/as.numeric(isos[4,
                                                                                                                                                  road[i, 1]]))
                peaks[m, 1] = round(peaks[k, 1] - (as.numeric(isos[3,
                                                                   road[i, 1]])) + (as.numeric(isos[3,
                                                                                                    road[i, j]])), digits = 9)
                peaks[m, 4 + road[i, 1]] = c(peaks[m,
                                                   4 + road[i, 1]] - 1)
                peaks[m, 4 + road[i, j]] = c(peaks[m,
                                                   4 + road[i, j]] + 1)
                if (peaks[m, 2] < minlimit) {
                  peaks[m, 4] = 1}
                m <- m + 1}}}}
        k <- c(k + 1)}
      else {
        k <- c(k + 1)}}
    b_end <- c(m - 1)
    b_start <- c(k)
    if (b_end > b_start) {
      getit <- seq(b_start, b_end, 1)
      getthat <- c()
      back <- getit[order(peaks[getit, 1], decreasing = FALSE)]
      peaksort <- peaks[back, ]
      for (i in 1:(length(peaksort[, 1]) - 1)) {
        if (peaksort[i, 1] == peaksort[i + 1, 1]) {
          if (round(peaksort[i, 2], digits = 3) == round(peaksort[i +
                                                                  1, 2], digits = 3)) {
            if (all(peaksort[i, c(5:leng2)] == peaksort[i +
                                                        1, c(5:leng2)])) {
              getthat <- c(getthat, back[i + 1])}}}}
      leng4 <- length(getthat)
      if (leng4 > 0) {
        peaks <- peaks[-getthat, ]
        m <- c(m - leng4)
        b_end <- c(b_end - leng4)
        leng3 <- c(leng3 - leng4)
        rm(peaksort)}}
    g <- c(g + 1)}
  if (m >= leng3) {
    warning("Storage maximum for number of peaks reached!")}
  rm(road, isonumber, k, b_start, b_end, leng1, leng2, leng3)
  peaks <- peaks[peaks[, 1] != 0, ]
  if (m > 2) {
    peaks <- peaks[peaks[, 2] != 0, ]
    peaks <- peaks[order(peaks[, 1], decreasing = FALSE),]
    }
  return(peaks)
}


#' Title
#'
#' @param number_input
#' @param Peakinfo
#' @param index.input
#' @param iso_list
#' @param rawdata
#' @param polarity
#' @param ppm
#' @param ppm.ms2
#' @param Isotope.Data
#' @param Database
#' @param Fragment
#' @param adducts
#' @param IsotopeDB
#' @param mwoffset
#'
#' @return
#' @export
#'
# The summary data frame ds is used to plot larger red points on top
formcal<-function(polarity,number_input,Peakinfo,index.input,indexID,iso_list,rawdata,ppm,ppm.ms2,Isotope.Data,Database,Fragment,adducts,IsotopeDB,mwoffset){##number input is the limit of element composition number, isolist is the mw of element
  #########restrict the element number
  mz<-Peakinfo$mz[index.input]
  Adducts.pred<-unlist(strsplit(Peakinfo$Adducts[index.input],';'))
  Adducts.pred<-Adducts.pred[indexID]
  mserror<-unlist(strsplit(Peakinfo$mserror[index.input],';'))
  mserror<-as.numeric(mserror[indexID])
  SmiID<-unlist(strsplit(Peakinfo$SMILES[index.input],';'))
  SmiID<-SmiID[indexID]

  ###################constrain adducts###################
  formula.index<-which(Database$SMILES==SmiID)
  formula.index<-formula.index[1]
  if (length(formula.index)==0){return(NULL)}##no possible formulas
  formula.pred<-as.character(Database$formula[formula.index])
  SMILES.pred<-as.character(Database$SMILES[formula.index])
  DatabaseID<-formula.index

  ###############delete those formulas with rare elements#############
  test.element<-'CHONPSClBrIF0123456789'
  index<-NULL
  for (i in 1:length(formula.pred)){
    element.test<-strsplit(formula.pred[i],'')
    element.test<-unlist(element.test)
    match.element<-sapply(element.test,grep,test.element)
    if (length(match.element)==length(unlist(match.element))){##no rare elements
      index<-c(index,i)
    }}
  if(is.null(index)){return(NULL)}
  formula.pred<-formula.pred[index]
  Adducts.pred<-Adducts.pred[index]
  SMILES.pred<-SMILES.pred[index]
  mserror<-mserror[index]


  ####################isotope distribution##############################
  formula.split<-sapply(formula.pred,form.split)
  RMSE.iso1<-0
  if (length(ncol(formula.split))>0){
    form.offset<-sapply(Adducts.pred,DefineAdducts)
    formula.split<-formula.split+form.offset[1:11,]
    index.save<-NULL
    for (kk in 1:ncol(formula.split)){#Check if the number of elment is lesser than -1
      if (min(formula.split[,kk])<0){
        index.save<-c(index.save,kk)
      }
    }
    if (length(index.save)>0){
      formula.split<-formula.split[,-index.save]
    }
    formula.split<-matrix(formula.split,nrow=11)
    if (length(formula.split)==0){return(NULL)}
    RMSE.iso1<-apply(formula.split,2,deiso.formula,iso_list,index.input,Peakinfo,Isotope.Data,polarity)
  }
  if(RMSE.iso1[1]==0){return(NULL)}###error in isotope calculation because of rare element
  RMSE.iso<-exp((-0.5)*RMSE.iso1)

  ###############calculate the ms1 error
  MS1.score<-rep(0,length(mserror))
  index.small<-which(abs(mserror)<1)
  if (length(index.small)>0){MS1.score[index.small]<-1}##when smaller than 1ppm, no difference
  index.big<-which(abs(mserror)>1)
  if (length(index.big)>0){
    MS1.score[index.big]<-exp((-0.5)*(mserror[index.big]/0.5)^2)}##0.5ppm for delta


  ###############calculate the subset of formulas and ms2 error##
  fragments<-which(Fragment$libraryid==index.input)
  index.fragments<-which(Fragment$intensity[fragments]>2000)##only abundant fragments are used
  ms2score<-rep(0,length(formula.pred))
  for (j in 1:length(formula.pred)){
    if (length(index.fragments)<1){break}##no fragments
    fragments<-fragments[index.fragments]
    number.limit<-form.split(formula.pred[j])##subset of the precusor formula
    form.offset<-sapply(Adducts.pred,DefineAdducts)
    number.limit<-number.limit+form.offset[1:11,]
    number.limit[8]<-number.limit[7]
    number.limit[10]<-number.limit[9]
    element.DAT<-c("C","H","N","O","P","S","35Cl","37Cl","79Br","81Br","F")
    mw.element<-iso.element(element.DAT,iso_list,1)##calculate mw for each element
    number.low<-rep(0,11)
    numberset<-cbind(number.low,number.limit,element.DAT)
    number_list<-c(1,3,4,3,2,3,2,2,2,2,2,2)#the ratio of element to carbon
    ###################fragments##################
    score.temp<-0

    for (k in 1:length(fragments)){
      mz.frag<-Fragment$mz[fragments[k]]
      mz.pred<-mz.frag+0.0005485799*polarity###the mz used for formula prediction, should consider electrons
      formula.frag<-formpred(mz.pred,ppm.ms2,numberset,mw.element[,3],number_list,0)##primary calculation
      if (max(formula.frag)==0){score.temp<-score.temp-1}else{
        score.temp<-score.temp+1
      }
      }
    ms2score[j]<-score.temp/length(fragments)
    }

  #################total score####################################
  MS1.score<-log(MS1.score+0.01)
  RMSE.iso<-log(RMSE.iso+0.01)
  all.score<-MS1.score+RMSE.iso+ms2score
  if (length(all.score)<1){
    return(NULL)}
  rank1<-order(-all.score)
  index<-rank1[1:length(rank1)]
  if (length(index)>1){
    formula.final<-cbind(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])}else{
      formula.final<-c(formula.pred[index],Adducts.pred[index],SMILES.pred[index],DatabaseID[index],MS1.score[index],RMSE.iso[index],ms2score[index],all.score[index])}
  formula.final<-matrix(formula.final,ncol=8)
  return(formula.final)
  }###just output the best one


deiso.formula<-function(vector1,iso_list,index.input,Peakinfo,Isotope.Data,polarity){
  mz_tol<-5*10^(-6)
  formula.cal<-form.comb(vector1[1:11])
  if (max(vector1)==0) {return(10)}
  iso.pattern<-isopattern_cal(iso_list,formula.cal,1e-4)
  iso.pattern<-iso.pattern[,1:2]
  iso.pattern[,2]<-iso.pattern[,2]/(max(iso.pattern[,2]))
  iso.pattern<-iso.pattern[which(iso.pattern[,2]>0.01),]
  iso.pattern<-combine.mz(iso.pattern,0.002)###combine those peaks could not be distinguished by mass spec
  iso.pattern[,1]<-iso.pattern[,1]-0.0005485799*polarity##plus the decoy adducts


  index<-which(Isotope.Data$mz==Peakinfo$mz[index.input])
  if(length(index)==0)(return(10))

  iso.intensity<-Isotope.Data$intensity.iso[index]###isotopic peaks from real spectra
  iso.mz<-Isotope.Data$Isotope[index]
  mz<-Peakinfo$mz[index.input]
  index<-which(abs(iso.pattern[,1]-mz)<mz_tol*mz)
  if (length(index)==0){return(10)}
  normal.ratio<-iso.pattern[index[1],2]###the ratio for mz
  index<-which(abs(iso.mz-mz)<mz_tol*mz)
  if (length(index)==0){return(10)}
  iso.intensity<-iso.intensity*normal.ratio/iso.intensity[index[1]]


  match.1<-data.frame(iso.pattern)###match expt spectra to predicted spectra
  match.1$match<-rep(0,nrow(iso.pattern))
  for (i in 1:nrow(iso.pattern)){
    index<-which(abs(match.1$mass[i]-iso.mz)<mz_tol*match.1$mass[i])
    if (length(index)>0){
      match.1$match[i]<-iso.intensity[index[1]]
    }
  }

  match.2<-data.frame(mz=iso.mz,intensity=iso.intensity)###match expt spectra to predicted spectra
  match.2$match<-rep(0,nrow(match.2))
  for (i in 1:nrow(match.2)){
    index<-which(abs(match.2$mz[i]-match.1$mass)<mz_tol*match.2$mz[i])
    if (length(index)>0){
      match.2$match[i]<-match.1$abundance[index[1]]
    }
  }

  RMSE.1<-sum((match.1[,2]-match.1[,3])^2)
  RMSE.2<-sum((match.2[,2]-match.2[,3])^2)
  return(RMSE.1+RMSE.2)
}


##############function to combine mz using mz_tol###########
combine.mz<-function(mz.list,mz.tol){
  mz.new<-NULL
  for (i in 1:nrow(mz.list)){
    if (i==1){
      mz.new<-mz.list[1,]
    }else{
      if (length(mz.new)==ncol(mz.list)){
        index<-which(abs(mz.new-mz.list[i,1])<mz.tol)
      }else{
        index<-which(abs(mz.new[,1]-mz.list[i,1])<mz.tol)}
      if (length(index)<1){
        mz.new<-rbind(mz.new,mz.list[i,])
      }else{
        if (length(mz.new)==ncol(mz.list)){
          mz.value<-sum(mz.new[1],mz.list[i,1])/2
          intensity<-sum(mz.new[2],mz.list[i,2])
          mz.new<-c(mz.value,intensity)
        }
        else{
          mz.value<-sum(mz.new[index,1],mz.list[i,1])/2
          intensity<-sum(mz.new[index,2],mz.list[i,2])
          mz.new[index[1],]<-c(mz.value,intensity)}
      }}}
  return(mz.new)}


#####################calculate the score for a given isotope distribution##############
#####################for actural library, adjusted by the predicted mz at first, and then normalized by summed value from the library##########
######for the first adjustment, if use maxmal peak or summed peak for adjustment, there may be big effects from noise########
######for the second adjustment, if use maximal peak for adjustment, there will be bias towards low isotopes###########
isotope.score<-function(iso.pattern,iso.info){##iso.pattern and iso.info is the library and small number of peaks,
  index1<-ceiling(length(iso.info)/4)
  index<-which(abs(iso.pattern[,1]-iso.info[index1,1])<0.006)#index is the ID of mz for alignment
  factor<-iso.pattern[index,2]/sum(iso.pattern[,2])
  iso.pattern[,2]<-iso.pattern[,2]/sum(iso.pattern[,2])###this step is to use sum value for adjustment to avoid any bias
  iso.info[,2]<-iso.info[,2]*factor/(iso.info[index1,2])###normalized based on known mz
  RMSE<-0
  index.act<-rep(0,length(iso.info)/2)
  for (i in 1:(length(iso.pattern)/2)){##calculate the RMSE from library
    index.temp<-which(abs(iso.info[,1]-iso.pattern[1])<0.006)
    if (length(index.temp)>0){
      RMSE.temp<-(iso.pattern[2]-iso.info[index.temp,2])^2
      index.act[index.temp]<-i}else{##record the actual peaks have been calcualted
        RMSE.temp<-(iso.pattern[2])^2}
    RMSE<-RMSE+RMSE.temp}
  for (i in 1:(length(iso.info)/2)){##calculate the RMSE from unaligned actual peaks
    if (index.act[i]==0){##not aligned to library
      RMSE<-RMSE+(iso.info[2])^2}}
  return(RMSE)}


####################calculate probability of distribution of a given formulae based on prior information
density.form<-function(vector1){##vector1 is the number of element, C, H, O, N, P, S, Cl, 37Cl, Br, 81Br
  Hnumber<-vector1[2]+sum(vector1[7:11])
  Hnumber<-Hnumber/vector1[1]
  vector1<-vector1/vector1[1]
  prob<-exp((-0.5)*(Hnumber-1.6)^2/1.6^2)+exp((-0.5)*(vector[5]-0)^2/0.15^2)+exp((-0.5)*(vector[5]-0.1)^2/0.4^2)
}

#------------------------------
#define adducts
#------------------------------
DefineAdducts<-function(adducts){
  element<-NULL
  if (adducts=='[M-H]-'){
    element<-c(0,-1,0,0,0,0,0,0,0,0,0,-1.007825)
  }
  if (adducts=='[M-H-H2O]-'){
    element<-c(0,-3,0,-1,0,0,0,0,0,0,0,-19.01839)
  }
  if (adducts=='[M+Cl]-'){
    element<-c(0,0,0,0,0,0,1,0,0,0,0,34.96885)
  }
  if (adducts=='[M+CH2O2-H]-'){
    element<-c(1,1,0,2,0,0,0,0,0,0,0,44.99765)
  }
  if (adducts=='[M]-'){
    element<-c(0,0,0,0,0,0,0,0,0,0,0,0)
  }
  if (adducts=='[M-Br+O]-'){
    element<-c(0,0,0,1,0,0,0,0,-1,0,0,-62.92342)
  }
  if (adducts=='[M+H]+'){
    element<-c(0,1,0,0,0,0,0,0,0,0,0,1.007825)
  }
  if (adducts=='[M+H-H2O]+'){
    element<-c(0,-1,0,-1,0,0,0,0,0,0,0,-17.00274)
  }
  if (adducts=='[M+CH3OH+H]+'){
    element<-c(1,5,0,1,0,0,0,0,0,0,0,33.034040)
  }
  if (adducts=='[M+NH4]+'){
    element<-c(0,4,1,0,0,0,0,0,0,0,0,18.034374)
  }
  return(element)
}

formpred<-function(mz,ppm,numberset,mz_list,number_list,mwoffset){
  numberset<-matrix(as.numeric(numberset[,1:2]),ncol=2,nrow=nrow(numberset))
  numberset[2,1]<-max(numberset[2,1],floor(numberset[1,1]*0.5-sum(numberset[7:11,2])))##the minimum number of bromine
  formulae<-rep(0,nrow(numberset))
  for.pred<-itercal(numberset,mz_list,ppm,mz,mwoffset)##mwoffset is for prediction of rare element for fragment
  if (length(for.pred)<2)
    for.pred<-NULL
  return(for.pred)
}
