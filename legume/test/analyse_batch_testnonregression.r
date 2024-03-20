



#determiner le path du fichier actuel et le recuper 
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#marche pas hors de rstudio/ligne de commande? (https://stackoverflow.com/questions/47044068/get-the-path-of-current-script)

#dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\test"



source(paste(dir, "fonctions_analyses.r",sep="\\"))
source(paste(dir, "fonctions_mef.r",sep="\\"))


#dir0 <- paste(dir, "control_v1.0",sep="\\")
dir0 <- paste(dir, "previouscheck",sep="\\")
dirlast <-  paste(dir, "lastcheck",sep="\\")
#dirlast <-  paste("C:\\outputs", "lastcheck",sep="\\")
#dirlast <- paste(dir, "test2",sep="\\")


#recupere dico de reference (a refaire si ajoute vaiables...)
setwd(dir0)
dtotoref <- read.csv("dicolast.csv")#("dico_controlv0.csv")

setwd(dirlast)#(dir0)#
ls_files <- list.files(dirlast)#(dir0)#


#recupere la liste des toto file names du dossier de travail
ls_toto <- ls_files[grepl('toto', ls_files)]
ls_paramSD <- ls_files[grepl('paramSD', ls_files)]



#creation du dataFrame dtoto et recup des info fichier


#11 col (avec sd)
cols_ <- strsplit(ls_toto, '_')
test_long <- as.numeric(lapply(cols_, length)) #pour separer selon nb de champs (avec sd)

dtoto <- as.data.frame(t(as.data.frame(cols_[test_long==11])))#as.data.frame(t(as.data.frame(strsplit(ls_toto, '_'))))#
row.names(dtoto) <- 1: length(dtoto[,1])
dtoto <- dtoto[,c(2,3,4,5,6,7,8,10)]
names(dtoto) <- c('usm','lsystem','mix','damier','scenario','Mng', 'seed','sd')
dtoto$name <- ls_toto[test_long==11]
dtoto$seed <- as.numeric(as.character(dtoto$seed))#substr(as.character(dtoto$seed), 1, 1)
dtoto$scenario <- substr(as.character(dtoto$scenario), 9, nchar(as.character(dtoto$scenario)))

#dtoto <- rbind(temp, dtoto) #merge des 2
dtoto$keysc <- paste(dtoto$scenario, dtoto$mix, dtoto$Mng, dtoto$sd)# ajout d'une cle unique par scenario
#dtoto$damier <- as.numeric(substr(as.character(dtoto$damier), 7, 7))




#split de dtoto et stockage dans une liste de scenatios
sp_dtoto <- split(dtoto, dtoto$keysc)

#??info densite->lu dans toto?
#reprendre les calculs pa usm




#didcols <- as.data.frame(residcols[seq(1,21,3), ])
#didcols$damier <- unique(as.character(dtoto$damier))
#write.csv(didcols, "didcols.csv", row.names=F)

#fichier d'id colones esp1 pour damier 8 (pour les cas ou bug/oubli dans les noms de colonnes)
#didcols <- read.csv("C:/devel/l-egume/legume/multisim/didcols.csv") 




for (key in names(sp_dtoto))#key <- names(sp_dtoto)[1]
{
  ls_toto_paquet <- sp_dtoto[[key]]$name
  
  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  ltoto <- read_ltoto(ls_toto_paquet)
  #version locale du paquet de doto
  dtoto <- sp_dtoto[[key]]
  
  #recup du nom des esp
  mix <- strsplit(ls_toto_paquet[1], '_')[[1]][4] #suppose paquet fait par traitement
  esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
  esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'
  
  #visu des rendement moyen m2 / a un DOY 
  surfsolref <- NULL
  nbplt <- NULL
  nbplt1 <- NULL
  nbplt2 <- NULL
  
  #DOYScoupe <- c(165,199,231,271,334)#Avignon
  DOYScoupe <- c(187,229,282,334)#Lusignan an1
  #DOYScoupe <- c(463,504,557,609)#Lusignan an2
  DOYdeb <- 60#
  idDOYScoupe <- DOYScoupe - DOYdeb
  Ytot <- NULL
  Ycoupe <- NULL
  
  YEsp1 <- NULL
  YEsp2 <- NULL
  
  QNfix <- NULL
  QNupttot <- NULL
  QNuptleg <- NULL
  
  PARi1 <- NULL
  PARi2 <- NULL
  Surf1 <- NULL
  Surf2 <- NULL
  LRac1 <- NULL
  LRac2 <- NULL
  MRac1 <- NULL
  MRac2 <- NULL
  
  for (i in 1:length(ls_toto_paquet))#(ls_toto))
  {
    name <- ls_toto_paquet[i]
    damier <- strsplit(name, '_')[[1]][5]
    dat <- ltoto[[name]]
    s <- dat[dat$V1=='pattern',3]#m2
    surfsolref <- cbind(surfsolref, as.numeric(as.character(s)))
    nb <- length(dat)-2
    nbplt <- cbind(nbplt, nb)
    
    #Y Totaux
    MSaerien <- as.matrix(dat[dat$V1=='MSaerien' & dat$steps %in% DOYScoupe,3:(3+nb-1)], ncol=nb)
    ProdIaer <- rowSums(MSaerien) / s
    Ycoupe <- rbind(Ycoupe, ProdIaer)
    Ytot <- cbind(Ytot, sum(ProdIaer))#cumul des 5 coupes
  
    
    #N totaux et fixation
    Qfix <- as.matrix(dat[dat$V1=='Qfix',3:(3+nb-1)], ncol=nb)
    Qfix <- as.numeric(rowSums(Qfix) / s)
    Nuptake_sol_tot <- as.matrix(dat[dat$V1=='Nuptake_sol',3:(3+nb-1)], ncol=nb)
    Nuptake_sol_tot <- as.numeric(rowSums(Nuptake_sol_tot) / s)
    QNfix <-cbind(QNfix, sum(Qfix))
    QNupttot <- cbind(QNupttot, sum(Nuptake_sol_tot))
  
    #YEsp1
    #esp <- 'Fix2'#'Fix3'#'Fix1'#'Fix' #pourquoi c'est ce nom au lieu de Fix???
    #esp2 <- 'nonFixSimTest'#'nonFix1'#'nonFix0' #pourquoi c'est ce nom au lieu de Fix???
    
    nomcol <- names(ltoto[[name]])
    #if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    #{
    #  idcols <- as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
    #} else
    #{
      idcols <- grepl(esp, nomcol) & !grepl(esp2, nomcol)#contient esp1 et pas esp2
    #}
    
    dat1 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
    nb1 <- length(dat1)-2
    nbplt1 <- cbind(nbplt1, nb1)
    if (nb1>0)
    {
      MS1 <- as.matrix(dat1[dat1$V1=='MSaerien' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)
      ProdIaer1 <- rowSums(MS1) / s
      Nuptake_sol_leg <- as.matrix(dat1[dat1$V1=='Nuptake_sol',3:(3+nb1-1)], ncol=nb1)
      Nuptake_sol_leg <- as.numeric(rowSums(Nuptake_sol_leg) / s)
      jPARi1 <- rowSums(as.matrix(dat1[dat1$V1=='PARiPlante' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jSurf1 <- rowSums(as.matrix(dat1[dat1$V1=='SurfPlante' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jLRac1 <- rowSums(as.matrix(dat1[dat1$V1=='RLTot' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jMRac1 <- rowSums(as.matrix(dat1[dat1$V1=='MS_rac_fine' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jMPiv1 <- rowSums(as.matrix(dat1[dat1$V1=='MS_pivot' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      
    } else
    {
      ProdIaer1 <- 0 #pas de plante de l'esp1
      Nuptake_sol_leg <- 0
      jPARi1 <- 0
      jSurf1 <- 0
      jLRac1 <- 0
      jMRac1 <- 0
      jMPiv1 <- 0
    }
    YEsp1 <- cbind(YEsp1, sum(ProdIaer1))#cumul des 5 coupes
    QNuptleg <- cbind(QNuptleg, sum(Nuptake_sol_leg))
    PARi1 <- cbind(PARi1, sum(jPARi1))
    Surf1 <- cbind(Surf1, sum(jSurf1))
    LRac1 <- cbind(LRac1, max(jLRac1))
    MRac1 <- cbind(MRac1, max(jMRac1)+max(jMPiv1))
    
    #YEsp2
    #if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    #{
    #  idcols <- !as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
    #  idcols[1:2] <- FALSE #remet a faux les deux premieres colonnes
    #} else
    #{
      idcols <- grepl(esp2, nomcol)#contient esp2
    #}
    
    dat2 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
    nb2 <- length(dat2)-2
    nbplt2 <- cbind(nbplt2, nb2)
    if (nb2>0)
    {
      MS2 <- as.matrix(dat2[dat2$V1=='MSaerien' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)
      ProdIaer2 <- rowSums(MS2) / s
      jPARi2 <- rowSums(as.matrix(dat2[dat2$V1=='PARiPlante' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jSurf2 <- rowSums(as.matrix(dat2[dat2$V1=='SurfPlante' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jLRac2 <- rowSums(as.matrix(dat2[dat2$V1=='RLTot' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jMRac2 <- rowSums(as.matrix(dat2[dat2$V1=='MS_rac_fine' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jMPiv2 <- rowSums(as.matrix(dat2[dat2$V1=='MS_pivot' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      
    }
    else
    {
      ProdIaer2 <- 0 #pas de plante de l'esp2
      jPARi2 <- 0
      jSurf2 <- 0
      jLRac2 <- 0
      jMRac2 <- 0
      jMPiv2 <- 0
    }
    YEsp2 <- cbind(YEsp2, sum(ProdIaer2))#cumul des 5 coupes
    #YEsp2 <- cbind(YEsp2, sum(ProdIaer2))#cumul des 5 coupes
    PARi2 <- cbind(PARi2, sum(jPARi2))
    Surf2 <- cbind(Surf2, sum(jSurf2))
    LRac2 <- cbind(LRac2, max(jLRac2))
    MRac2 <- cbind(MRac2, max(jMRac2)+max(jMPiv2))
    
  }
  
  dtoto$surfsolref <- as.numeric(surfsolref)
  dtoto$nbplt <- as.numeric(nbplt)
  dtoto$nbplt1 <- as.numeric(nbplt1)
  dtoto$nbplt2 <- as.numeric(nbplt2)
  dtoto$Ytot <- as.numeric(Ytot)
  dtoto$densite <- dtoto$nbplt/dtoto$surfsolref
  dtoto$densite1 <- dtoto$nbplt1/dtoto$surfsolref
  dtoto$YEsp1 <- as.numeric(YEsp1)
  dtoto$densite2 <- dtoto$nbplt2/dtoto$surfsolref
  dtoto$YEsp2 <- as.numeric(YEsp2)
  dtoto$Semprop1 <- dtoto$densite1/dtoto$densite
  dtoto$Yprop1 <- dtoto$YEsp1 / (dtoto$YEsp1 +dtoto$YEsp2)
  dtoto$Yprop2 <- dtoto$YEsp2 / (dtoto$YEsp1 +dtoto$YEsp2)
  dtoto$QNfix <- as.numeric(QNfix)
  dtoto$QNupttot <- as.numeric(QNupttot)
  dtoto$QNuptleg <- as.numeric(QNuptleg)
  dtoto$QNtot <- dtoto$QNfix + dtoto$QNupttot

  #new var
  dtoto$Pari1 <- as.numeric(PARi1)
  dtoto$Pari2 <- as.numeric(PARi2)
  dtoto$Surf1 <- as.numeric(Surf1)
  dtoto$Surf2 <- as.numeric(Surf2)
  dtoto$PhiSurf1 <- as.numeric(PARi1) / (as.numeric(Surf1) + 10e-12)#Phi Surf
  dtoto$PhiSurf2 <- as.numeric(PARi2) / (as.numeric(Surf2) + 10e-12)
  dtoto$PhiMass1 <- as.numeric(PARi1) / (as.numeric(YEsp1) + 10e-12)#Phi Mass
  dtoto$PhiMass2 <- as.numeric(PARi2) / (as.numeric(YEsp2) + 10e-12)
  dtoto$LRac1 <- as.numeric(LRac1)
  dtoto$LRac2 <- as.numeric(LRac2)
  dtoto$MRac1 <- as.numeric(MRac1)
  dtoto$MRac2 <- as.numeric(MRac2)
  dtoto$UptNLen1 <- (as.numeric(QNupttot) - as.numeric(QNuptleg)) / (as.numeric(LRac1) + 10e-12)#Uptake par Len
  dtoto$UptNLen2 <- as.numeric(QNuptleg) / (as.numeric(LRac2) + 10e-12)
  dtoto$UptNMass1 <- (as.numeric(QNupttot) - as.numeric(QNuptleg)) / (as.numeric(MRac1) + 10e-12)#Uptake par Mass root
  dtoto$UptNMass2 <- as.numeric(QNuptleg) / (as.numeric(MRac2) + 10e-12)
  
  
  
  
  #remise du dtoto locl dans sp_dtoto
  sp_dtoto[[key]] <- dtoto
  
}

#reagrege dtoto
dtoto <- do.call("rbind", sp_dtoto)
row.names(dtoto) <- 1:length(dtoto$usm)

#prevoir calcul autres variables! (racines, stressN, eau..., biomasse pa coupe)

#write.csv(dtoto, "dico_controlv0.csv", row.names=F)
write.csv(dtoto, "dicolast.csv", row.names=F)

sp_dtoto <- split(dtoto, dtoto$keysc)
sp_dtotoref <- split(dtotoref, dtotoref$keysc)





#
#nomvar <- "QNtot" #"Ytot"# 
#keysc <- names(sp_dtoto)[4]

plot_nonregression <- function(dtoto, dtotoref, nomvar)
{
  #split
  sp_dtoto <- split(dtoto, dtoto$keysc)
  sp_dtotoref <- split(dtotoref, dtotoref$keysc)
  
  #x,y globaux
  x <- dtoto[,nomvar]#sp_dtoto[[keysc]][,nomvar]
  y <- dtotoref[,nomvar]#sp_dtotoref[[keysc]][,nomvar]
  
  plot(x, y, xlab=nomvar, ylab=paste(nomvar,"ref0"), main=paste("non_reg",nomvar), xlim=c(0,1.5*max(x)), ylim=c(0,1.5*max(x)))
  #droite 1:1
  abline(0,1)
  
  res_rmse <- signif(rmse(x,y),3) #RMSE
  res_rrmse <- signif(rrmseCoucheney(x,y),3)#rRMSE
  res_rmses <- rmsesCoucheney(x,y)#RMSEs
  res_rmseu <- rmseuCoucheney(x,y)#RMSEu
  res_prmses <- signif(pRMSEs(rmse(x,y), res_rmses),3)#pRMSEs
  res_prmseu <- signif(pRMSEu(rmse(x,y), res_rmseu),3)#pRMSEu
  res_EF <- signif(efficiencyCoucheney(x,y),3)#Efficiency
  res_r2 <- signif(summary(lm(y~x))$r.squared,3)#r2
  #faire une fonction qui renvoie une liste
  
  text(1.5*max(x) , 0.4*max(x), adj=1, cex=0.8, paste("rmse: ", res_rmse))
  text(1.5*max(x) , 0.3*max(x), adj=1,  cex=0.8, paste("pRMSEu: ", res_prmseu))
  text(1.5*max(x) , 0.2*max(x),  adj=1, cex=0.8, paste("EF: ", res_EF))
  text(1.5*max(x) , 0.1*max(x), adj=1,  cex=0.8, paste("R2: ", res_r2))
  
  # droite de regression globale
  
  
  #par couleur les differents traitements
  col=1
  for (keysc in names(sp_dtoto))
  {
    x <-sp_dtoto[[keysc]][,nomvar]
    y <- sp_dtotoref[[keysc]][,nomvar]
    points(x, y, col=col)
    col=col+1
  }
}





#ceation d'un rapport pour non regression

nomrap <- paste( 'rapport_nonregression',Sys.Date(),basename(dir0),basename(dirlast),'.pdf', sep="_")
pdf(paste(dir,nomrap, sep='\\'), onefile=T)


#Constant Final Yield Law (Fig 6.1a Louarn & Faverjon 2018)
cexval <- 1.5
x <- sp_dtoto[["-1--1 Fix2-Fix2 Lusignan30IrrN -"]]
plot(x$densite2, x$Ytot, xlab='Seeding density (plant.m-2)', ylab='Aboveground production (g.m-2)', ylim=c(0,2500),cex=cexval,cex.lab=cexval, main="Constant Final Yield - Fix2")
x <- sp_dtoto[["-1--1 Fix2-Fix2 Lusignan30 -"]]
points(x$densite2, x$Ytot, pch=16,cex=cexval)
x <- sp_dtoto[["-1--1 Fix2-Fix2 Lusignan30Irr -"]]
points(x$densite, x$Ytot, col='grey', pch=16,cex=cexval)

#avec mineralisation des residu active
x <- sp_dtoto[["1-1 Fix2-Fix2 Lusignan30IrrN -"]]
plot(x$densite2, x$Ytot, xlab='Seeding density (plant.m-2)', ylab='Aboveground production (g.m-2)', ylim=c(0,2500),cex=cexval,cex.lab=cexval, main="Constant Final Yield - Fix2 - Resid+")

#pour TB
x <- sp_dtoto[["-1--1 giga-giga Lusignan30IrrN -"]]
plot(x$densite2, x$Ytot, xlab='Seeding density (plant.m-2)', ylab='Aboveground production (g.m-2)', ylim=c(0,2200),cex=cexval,cex.lab=cexval, main="Constant Final Yield - giga")



#Neutral binary mixture avec fixation
tabmoy <- Build_AverageScTable(dtoto, keysc="55-1 Fix2-nonFixSimTest Lusignan30IrrN2 -")
YtotvsProp(tabmoy, nom="55-1 Fix2-nonFixSimTest Lusignan30IrrN2 -")


#Neutral mixture sans fixation with intraspecific variability
tabmoy <- Build_AverageScTable(dtoto, keysc="55-3 Fix2-nonFixSimTest Lusignan30IrrN2 SD2-2")
YtotvsProp(tabmoy, nom="55-3 Fix2-nonFixSimTest Lusignan30IrrN2 SD2-2")

#lecture fichier toto 50/50 (dat) 
key <- "55-3 Fix2-nonFixSimTest Lusignan30IrrN2 SD2-2"
damier <- "damier4" #50/50
ls_toto_paquet <- sp_dtoto[[key]]$name
ls_toto_paquetOK <- ls_toto_paquet[grepl(damier, ls_toto_paquet)]
ltoto <- read_ltoto(ls_toto_paquetOK)
IDfichier <- 1
nomfichier <- names(ltoto)[IDfichier]
dat <- ltoto[[nomfichier]]
num_usm <- strsplit(nomfichier, '_')[[1]][2]
scenar <- strsplit(nomfichier, '_')[[1]][6]
graine <- strsplit(nomfichier, '_')[[1]][8]
secenarSD <- strsplit(nomfichier, '_')[[1]][10]
esps <- strsplit(nomfichier, '_')[[1]][4]

#lecture du fichier paramSD de l'USM dans tabSD
nomSD <- ls_paramSD[grepl(num_usm, ls_paramSD)]
param_name <- "Len"
tabSD <- read.table(nomSD, header=T, sep=';')
MStot <- dat[dat$V1=='MStot',3:(3+nb-1)] #ajout de MStot
tabSD$MStotfin <- as.numeric(MStot[275,])

#split de tabSD par espece et ajout des deciles
sp_tabSD <- split(tabSD, tabSD$name)
lscol10 <- colorRampPalette(c("blue", "red"))( 11 ) #palette de couleur des deciles
sp <- unique(as.character(tabSD$name))[1]#"Fix2"#"nonFixSimTest"#
valparams <- sp_tabSD[[sp]][,c(param_name)]
sp_tabSD[[sp]]$decile <- Which_decile(valparams)
sp <- unique(as.character(tabSD$name))[2]#"nonFixSimTest"#
valparams <- sp_tabSD[[sp]][,c(param_name)]
sp_tabSD[[sp]]$decile <- Which_decile(valparams)


#plot Area par espece
for (sp in names(sp_tabSD))
{
  res <- Build_EvolProportions(MStot, sp_tabSD, sp)
  don <- res[,2:12] #t et les deciles
  titre <- paste(sp, num_usm, scenar, secenarSD, damier, graine)#esps,
  My_AreaPlot(don, titre=titre, xlab="t", ylab=paste("decile ", param_name), lscol=rev(lscol10))
}



#non egression par variable
ls_var <- c("Ytot", "QNtot", "Yprop2")

for (nomvar in ls_var)
{
  plot_nonregression(dtoto, dtotoref, nomvar)
}

dev.off()
#fin rapport


#pRMSEu pas bon! ou erreur d'arrondi??
#faire graph global puis couleur


#plot(dobssim$obs, dobssim$sim, xlim=c(0, 1.5*max(c(dobssim$obs,dobssim$sim))), ylim=c(0, 1.5*max(c(dobssim$obs,dobssim$sim))), xlab='', ylab='', main=main,...)



#regression / r2 / rRMSE...








#faire fonction qui fait le graph de comparaison pour n'importe quelle vaiable et pour tous les keysc
  
#faire une fonction rapport avec tous les graphs...
#voir sortie simulee: densite/Ytot











