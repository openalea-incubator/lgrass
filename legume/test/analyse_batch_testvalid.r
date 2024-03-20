library(ggplot2)
library(ggpubr) #to arrange ggplots


#determiner le path du fichier actuel et le recuper 
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#marche pas hors de rstudio/ligne de commande? (https://stackoverflow.com/questions/47044068/get-the-path-of-current-script)


#dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\test"


library("readxl")
source(paste(dir, "fonctions_analyses.r",sep="\\"))
source(paste(dir, "fonctions_mef.r",sep="\\"))


#le path avec les derniere simulation
dirlast <-  paste(dir, "lastvalidBis",sep="\\")
#dirlast <-  paste(dir, "lastvalid",sep="\\")
#dirlast <-  paste(dir, "test_champ",sep="\\") #pour visu dossier sorties champ
#dirlast <-  "C://inputs//inputs mayssa//new tests"
setwd(dirlast)#(dir0)#


#le path avec les fichier obs
pathobs <- paste(dir, "obs", sep="\\")



#liste les fichiers de simulation
ls_files <- list.files(dirlast)#(dir0)#


#unzip the files
ls_zip <- ls_files[grepl('.zip', ls_files)]
for (file_ in ls_zip)
{unzip(file_, exdir=dirlast)}
ls_files <- list.files(dirlast)#reliste les fichier apres dezippage


#recupere la liste des toto file names du dossier de travail
ls_toto <- ls_files[grepl('toto', ls_files) & !grepl('.zip', ls_files)]
ls_paramSD <- ls_files[grepl('paramSD', ls_files)]







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





#rapport des graphs dynamiques (rapport produit dans dossier du script R dans l-egume)
nomrap <- paste( 'rapport_eval-Dyn',Sys.Date(),basename(dirlast),'.pdf', sep="_")
pdf(paste(dir,nomrap, sep='\\'), onefile=T)


for (key in names(sp_dtoto))#
{
  #analyse par usm
  #key <- names(sp_dtoto)[1]#dileg luz


  ls_toto_paquet <- sp_dtoto[[key]]$name
  
  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  ltoto <- read_ltoto(ls_toto_paquet)
  #version locale du paquet de doto
  dtoto <- sp_dtoto[[key]]
  
  #recup du nom des esp et traitement
  mix <- strsplit(ls_toto_paquet[1], '_')[[1]][4] #suppose paquet fait par traitement
  esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
  esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'
  meteo <- strsplit(ls_toto_paquet[1], '_')[[1]][9]
  damier <- strsplit(ls_toto_paquet[1], '_')[[1]][5]
  
  #recup des obs correspondant
  namexl <- paste0(meteo, "_obs.xls")#"morpholeg14_obs.xls"
  trait <- if (esp == esp2 & damier=="homogeneous0") "ISO" else "HD-M2" #sera a adapter selon les melanges ou a renommer "timbale-krasno"
  trait <- if (esp == esp2 & damier=="homogeneous0" & meteo == "DigitLuz10") "LD" else trait
  if (meteo == "DivLeg15" | meteo == "LusignanDivLeg" | meteo == "LusignanAsso16")
  {
    trait <- "HD"
  }
  #trait <- if (esp == esp2 & damier=="homogeneous0" & meteo=="morpholegRGR15") "LD"
  onglet <- paste0(meteo, "_",trait,"_",esp)#"morpholeg14_ISO_timbale" #marche pour isole; a revoir pour autres
  #onglet <- "F_E1D1_B_R1" #"force!!
  obs <- read_excel(paste(pathobs,namexl,sep="\\"), sheet = onglet, col_names = TRUE, na = "")
  
  #calcul de la moyenne des simuls
  #pour l'espece 1
  simmoy <- build_simmoy(ltoto, lsusm=names(ltoto), esp)
  
    
  #fait les graph dynamiques
  #recup de surfsolref
  dat <- ltoto[[1]]
  surfsolref <- dat[dat$V1=='pattern',3] #m2
  #pour l'espece 1
  dynamic_graphs(simmoy, name=onglet, obs=obs, surfsolref=surfsolref) 
  
  #complement pour l'espece 2 si besoin
  if (esp != esp2 & damier!="homogeneous0")
  { 
    onglet2 <- paste0(meteo, "_",trait,"_",esp2)#"morpholeg14_ISO_timbale" #marche pour isole; a revoir pour autres
    obs2 <- read_excel(paste(pathobs,namexl,sep="\\"), sheet = onglet2, col_names = TRUE, na = "")
    simmoy2 <- build_simmoy(ltoto, lsusm=names(ltoto), esp2) 
    dynamic_graphs(simmoy2, name=onglet2, obs=obs2, surfsolref=surfsolref)
  }

}

dev.off()

#! NA dans fichier obs genere des bug de format + noms d'onglets a repredre...






#### test nouvelle boucle avec ggplot2 (gg_plotsim, gg_addplotobs, gg_plotObsSim)


#LAI <- moysimval(ltoto, lsusm=names(ltoto), var='SurfPlante')/ surfsolref
#LAIsd <- moysimval(ltoto, lsusm=names(ltoto), var='SurfPlante',optSD=T)/ surfsolref



#exemple analyse pour 1 usm pour usm 'key'
key <- names(sp_dtoto)[3]#dileg luz


ls_toto_paquet <- sp_dtoto[[key]]$name

#recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
ltoto <- read_ltoto(ls_toto_paquet)
#version locale du paquet de doto
dtoto <- sp_dtoto[[key]]

#recup du nom des esp et traitement
mix <- strsplit(ls_toto_paquet[1], '_')[[1]][4] #suppose paquet fait par traitement
esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'
meteo <- strsplit(ls_toto_paquet[1], '_')[[1]][9]
damier <- strsplit(ls_toto_paquet[1], '_')[[1]][5]

#recup des obs correspondant
namexl <- paste0(meteo, "_obs.xls")#"morpholeg14_obs.xls"
trait <- if (esp == esp2 & damier=="homogeneous0") "ISO" else "HD-M2" #sera a adapter selon les melanges ou a renommer "timbale-krasno"
trait <- if (esp == esp2 & damier=="homogeneous0" & meteo == "DigitLuz10") "LD" else trait
if (meteo == "DivLeg15" | meteo == "LusignanDivLeg" | meteo == "LusignanAsso16")
{
  trait <- "HD"
}
#trait <- if (esp == esp2 & damier=="homogeneous0" & meteo=="morpholegRGR15") "LD"
onglet <- paste0(meteo, "_",trait,"_",esp)#"morpholeg14_ISO_timbale" #marche pour isole; a revoir pour autres
#onglet <- "F_E1D1_B_R1" #"force!!
obs <- read_excel(paste(pathobs,namexl,sep="\\"), sheet = onglet, col_names = TRUE, na = "")



#calcul de la moyenne des simuls
# pour l'espece 1
simmoy <- build_simmoy(ltoto, lsusm=names(ltoto), esp)
simsd <- build_simmoy(ltoto, lsusm=names(ltoto), esp, optSD=T)
name <- onglet



#remet obs a dimension des sim (necessaire avec ggplot)
#garde seulement les obs dans le range des sim
obsOK <- obs[obs$DOY %in% simmoy$STEPS,]
obsMerge <- simmoy[,c("STEPS","TT")]
names(obsMerge)[1] <- "DOY"
obsMerge <- merge(obsMerge, obsOK, by="DOY",all=T)



#correspondance des noms de variable sim/obs
corresp <- data.frame(sim=c("NBI","LAI","MSA"), obs=c("NBI-quart","surf_tot","MSaerien"))



#cree une liste ls_plt avec les figures/plot simule dans simmoy
ls_plt <- vector("list", length=length(names(simmoy)))
names(ls_plt)  <- names(simmoy)

for (var_ in names(simmoy))
{
  ls_plt[[var_]] <- gg_plotsim(var_, simmoy, simsd, name)
}
#ls_plt[["NBI"]]




#ajout des obs pour les variables dispo dans la liste des plot simule (liste dans corresp)
for (var_ in names(ls_plt))
{
  if (var_ %in% corresp$sim)
  {
    ls_plt[[var_]] <- gg_addplotobs(ls_plt[[var_]], var_, obsMerge, corresp)
  }
}
#ls_plt[["NBI"]]



#arrange gg plots en pannels
ggarrange(ls_plt[["NBI"]], ls_plt[["NBphyto"]], ls_plt[["LAI"]] + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 1, nrow =3 )


ggarrange(ls_plt[["RDepth"]], ls_plt[["Hmax"]], ls_plt[["MSA"]] + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 1, nrow =3 )

ggarrange(ls_plt[["FTSW"]], ls_plt[["NNI"]], ls_plt[["R_DemandC_Root"]] + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 1, nrow =3 )

plot(simmoy$NNI, ylim=c(0,2), type='l')
ls_plt[["NNI"]]
ls_plt[["Npc_aer"]]
ls_plt[["MSArec"]]
ls_plt[["MSAnonrec"]]
ls_plt[["MSA"]]
ls_plt[["MSpiv"]]
ls_plt[["MSracfine"]]


#ajout d'une liste de variables optionelles a simmoy

#avec ponder surfsolref
for (var_ in c("dMSenNonRec", "dMSenPiv", "dMSenFeuil", "dMSenTige", "dMSenRoot","perteN_Piv","perteN_NonRec","NaerienNonRec", "MSsenaerien", "perteN_aerien", "dMSmortGel","dNmortGel", "dMSmortPlant_aer"))#"graineC", "graineN", "CreservPiv", "NreservPiv"
{
  simmoy[,var_] <- moysimval(ltoto, names(ltoto), esp, var=var_)/surfsolref
  simsd[,var_] <- moysimval(ltoto, lsusm=names(ltoto), esp, var=var_, optSD=T)/surfsolref
  ls_plt[[var_]] <- gg_plotsim(var_, simmoy, simsd, name)
}
ls_plt[["dMSenNonRec"]]
ls_plt[["dMSenPiv"]]
ls_plt[["dMSenFeuil"]]
ls_plt[["dMSenTige"]]
ls_plt[["dMSenRoot"]]
ls_plt[["MSsenaerien"]]
ls_plt[["perteN_Piv"]]
ls_plt[["perteN_NonRec"]]
ls_plt[["perteN_aerien"]]
ls_plt[["NaerienNonRec"]]
ls_plt[["dMSmortGel"]]
ls_plt[["dNmortGel"]]
ls_plt[["MSArec"]]
ls_plt[["MSA"]]
ls_plt[["dMSmortPlant_aer"]]
ls_plt[["MSpiv"]]


#avec ponder nbplt
var_ <- "Npc_piv"
for (var_ in c("Npc_piv", "Npc_aerNonRec"))
{
  simmoy[,var_] <- moysimval(ltoto, names(ltoto), esp, var=var_)
  simsd[,var_] <- moysimval(ltoto, lsusm=names(ltoto), esp, var=var_, optSD=T)
  ls_plt[[var_]] <- gg_plotsim(var_, simmoy, simsd, name)
}

ls_plt[["Npc_piv"]]
ls_plt[["Npc_aerNonRec"]]
ls_plt[["Npc_aer"]]
#pas bon car manque /nbplt!!! + pb nbp_piv

ls_plt[["R_DemandC_Root"]]
ls_plt[["RDepth"]]


var_ <- "dMSenTige"#"dMSenFeuil"#"dMSenRoot"#"dMSenPiv"#"dMSenNonRec"
x <- cumsum(simmoy[,var_])
plot(x, ylab=var_, type="l", col=2)


#test bilan
var_ <- "dMSenNonRec"
x <- cumsum(simmoy[,var_])
x2 <- simmoy[,"MSAnonrec"]
plot(x, ylab=var_, type="l", col=2, ylim=c(0,40))
points(x2, col=3,type="l")
points(x+x2, col=1, type="l")

var_ <- "dMSenPiv"
x <- cumsum(simmoy[,var_])
x2 <- simmoy[,"MSpiv"]
plot(x, ylab=var_, type="l", col=2, ylim=c(0,120))
points(x2, col=3,type="l")
points(x+x2, col=1, type="l")
plot(x2)

#test les obs_sim


#liste des plot obs_sim par variable presente dans corresp
name <- onglet

plt_obssim <- vector("list", length=length(corresp$sim))
names(plt_obssim ) <- as.character(corresp$sim)

for (var_ in as.character(corresp$sim))
{
  nomvarobs <- as.character(corresp[corresp$sim==var_,c("obs")])
  obssim <- na.omit(data.frame(obs=simmoy[,var_], sim=obsMerge[,nomvarobs]))
  plt_obssim[[var_]] <- gg_plotObsSim(obssim, var_, name=onglet)
}
#plt_obssim[["MSA"]]




ggarrange(plt_obssim[["NBI"]], plt_obssim[["LAI"]], plt_obssim[["MSA"]] + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow =2 )



# to check : decalage d'i jour des DOY sim
# normalisation des donnees obs










#if efface fenetre graphique de R!

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/


# to manipulate layers of existing ggplots

#https://cran.r-project.org/web/packages/gginnards/vignettes/user-guide-2.html

# to manipulate multifacet multi plots
#http://zevross.com/blog/2019/04/02/easy-multi-panel-plots-in-r-using-facet_wrap-and-facet_grid-from-ggplot2/



########################## pour tests


#faut lui donner moyenne et ecartypeen amont!
var_ <- "NBI"#"MSA"#"FTSW"#"NNI"#"LAI"#
nomvarobs <- as.character(corresp[corresp$sim==var_,c("obs")])



min <- 0
max <- 1.5*max(simmoy[,var_])

plot_var <- ggplot(data = simmoy, aes(x = STEPS)) +
  geom_line(aes(y = simmoy[,var_]), color="blue")+
  geom_ribbon(aes(ymin=simmoy[,var_]-simsd[,var_],ymax=simmoy[,var_]+simsd[,var_]),fill="blue",alpha=0.2)+
  geom_hline(yintercept=0)+
  ylim(min,max)+
  geom_text(x=1.20*min(simmoy$STEPS), y=0.98*max, label=name)+
  theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))+
  labs(title = "obs",subtitle = "sim",x = "DOY", y = var_)+
  theme(plot.title=element_text(size=10,color = "red"),plot.subtitle = element_text(size=10,color = "blue"))

plot_var



#ajout conditionnel des points d'obs
#if norma fait plante ggplot?? -> syntaxe specifique avec if inclu dans construction graph

plot_var2 <- plot_var + {if(var_ %in% corresp$sim) geom_point(aes(obsMerge$DOY, obsMerge[,nomvarobs]), fill="red",color="red" , size=2)}

plot_var2


#ls_plt[["NBI"]] <- plot_var
#ls_plt[[var_]] <- plot_var2


#plot_var3 <- gg_addplotobs(plot_var, "NBI", obsMerge, corresp)

#ajouter les SD des mesures quand disponible
#decaler les DOYsim de -1??







# faire une fonction gg_plotsim
# gg_addplotobs
# combine_ggplots




#obs-sim

var_ <- "NBI"#"MSA"#
nomvarobs <- as.character(corresp[corresp$sim==var_,c("obs")])
obssim <- na.omit(data.frame(obs=simmoy[,var_], sim=obsMerge[,nomvarobs]))
name <- onglet


min <- 0
max <- 25
reg   <- lm(obs ~ sim, data = obssim)
coeff <- coefficients(reg)
eq    <- paste0("y = ", round(coeff[2],2), "x + ", round(coeff[1],2))
RMSE_ <- round(rmse(obssim$obs,obssim$sim), 2)
rmses_ <- rmsesCoucheney(obssim$obs,obssim$sim)
rmse_ <- rmseuCoucheney(obssim$obs,obssim$sim)
pRMSEs_ <- round(pRMSEs(rmse_, rmses_), 2)
rRMSE <- round(rrmseCoucheney(obssim$obs,obssim$sim), 2)
EF <- round(efficiencyCoucheney(obssim$obs,obssim$sim), 2)

plot_ObsSim <- ggplot(obssim, aes(x = obs, y = sim)) +
  ggtitle(name)+
  geom_abline(intercept = 0, slope = 1, color = "black")+
  geom_point(aes(color = "obs"))+
  geom_smooth(method=lm, se = FALSE, color = "red")+
  ylim(min,max)+
  xlim(min,max)+
  geom_text(x=0.2*max, y=0.95*max, label=eq)+
  geom_text(x=0.2*max, y=0.9*max, label=paste("RMSE: ",RMSE_))+
  geom_text(x=0.2*max, y=0.85*max, label=paste("rRMSE: ",rRMSE))+
  geom_text(x=0.2*max, y=0.8*max, label=paste("pRMSEs: ",pRMSEs_))+
  geom_text(x=0.2*max, y=0.75*max, label=paste("EF: ",EF))+
  labs(x = paste("Obs ", var_), y = paste("Sim ", var_))


plot_ObsSim 






####################################
# A reprendre marche pas pour derniere simuls

#rapport des graphs obs-sim
nomrap <- paste( 'rapport_eval-Obs-Sim',Sys.Date(),basename(dirlast),'.pdf', sep="_")
pdf(paste(dir,nomrap, sep='\\'), onefile=T)



#construction du ls_dobsim d'une espece pour les plante isolee
esp_ <- 'timbale'#'giga'#'formica'#'sevanskij'#'leo'#'canto'#'kayanne'#
#ls_expe <- c('morpholeg14', 'morpholeg15')
ls_expe <- names(sp_dtoto)[grepl(esp_, names(sp_dtoto)) & grepl("homogeneous0", names(sp_dtoto)) & grepl("morpholeg", names(sp_dtoto))]#cle comportant le bon geno
#ls_expe <- names(sp_dtoto)[grepl(esp_, names(sp_dtoto)) & grepl("damier4", names(sp_dtoto))]#cle comportant le bon geno
ls_var <- c('NBI','nb_phyto_tot','surf_tot','Hmax','MSaerien')#,'long_pivot')
ls_varsim <- c('NBI','NBphyto','LAI', 'Hmax','MSA')#, 'RDepth')
#plante car ls_expe pas bon!


ls_dobssim <- build_ls_dobssim(esp_, ls_expe, ls_var, ls_varsim)
#erreur dans la recuperation de Rdepth de 2015 (dates avant le debut??) 
#erreur d'affihage qd tableau vide?


#figure des obs sim de l'espece
layout(matrix(1:6,2,3, byrow=T))
for (var in ls_varsim) #var <- "NBI"#"MSA"#"LAI"#"NBphyto"#
{
  #var <- "NBI"
  mobssim <- merge_dobssim(ls_dobssim[[var]])
  plot_obssim(mobssim, name=paste(esp_, "ISO/LD", var), displayusm=T)
}


#construction du ls_dobsim d'une espece pour les couverts denses
esp_ <- 'timbale'#'giga'#'formica'#'sevanskij'#'leo'#'canto'#'kayanne'#
#ls_expe <- c('morpholeg14', 'morpholeg15')
#ls_expe <- names(sp_dtoto)[grepl(esp_, names(sp_dtoto)) & grepl("homogeneous0", names(sp_dtoto))]#cle comportant le bon geno
ls_expe <- names(sp_dtoto)[grepl(esp_, names(sp_dtoto)) & grepl("damier4", names(sp_dtoto))]#cle comportant le bon geno
ls_var <- c('NBI','nb_phyto_tot','surf_tot','Hmax','MSaerien')#,'long_pivot')
ls_varsim <- c('NBI','NBphyto','LAI', 'Hmax','MSA')#, 'RDepth')


ls_dobssim <- build_ls_dobssim(esp_, ls_expe, ls_var, ls_varsim)
#erreur dans la recuperation de Rdepth de 2015 (dates avant le debut??) 
#erreur d'affihage qd tableau vide?



#figure des obs sm de l'espece
layout(matrix(1:6,2,3, byrow=T))
for (var in ls_varsim) #var <- "NBI"#"MSA"#"LAI"#"NBphyto"#
{
  mobssim <- merge_dobssim(ls_dobssim[[var]])
  plot_obssim(mobssim, name=paste(esp_, "HD", var), displayusm=T)
}


dev.off()




#finir le rapport

#tester avec melange a 2 especes -> OK pour dynamique / bug pas non plus pour obs-sim

#seprer simul plante isole et melanges dans les eval....? dans ls_expe


#fait par lucas
#separer NBI pousse ini et repousse dans eval + pour HD sim=comptage sur les decile superrieur de tige qui poussent

#visiblement pb d'unite dans les figures de NBphyto et LAI de combileg??? -> pb des 7 premieres tiges??




