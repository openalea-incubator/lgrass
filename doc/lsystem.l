#include <lpfgall.h>
#include <stdio.h>
#include <math.h>
#include <minmax.h>


#define GROUNDSIZE 60
#define GROUND SP()Down(90)f(GROUNDSIZE)Left(90)F(GROUNDSIZE)Left(90)F(GROUNDSIZE)Left(90)F(GROUNDSIZE)Left(90)F(GROUNDSIZE)Left(90)F(GROUNDSIZE)EP()


#define SORTIE "Sortie.xls"
#define CHAINE "C:\Chaine.txt"
#define GEOMETRIE "Geometrie.txt"
#define LUMIERE "lumiere.txt"           // Declaration du tableau retour des lumieres visualea
#define LUMIEREINIT "lumiere_init.txt"  // Initialisation tableau retour des lumieres visualea

FILE*fps;
FILE*fpss;
FILE*fpsss;

#define Defoliation 2  /* Pour les groupes */
#define Croissance 1
#define Photosynth 3
#define FOLD 7
#define DECOMPOECRIT 5
#define ECRITURE 6
#define LECTURE 8
#define ALLOCATION 9



#include "para_lgrass.txt"
#include LUMIERE;

double LSTEP=2;
double Biomasse;
double Densitefeuille=1;      /* g/unité de longeur - pour calculer biomasse recoltée */
double sommlineairecoupee;    /* Somme des longueurs de feuilles coupee */
double Torsiongaine=0;        /*Angle torsion des gaines */
float curve=0.3*LSTEP;              /*Angle courbure des limbe */
float INSERTIONLIMBE=40;      /* Angle Insertion limbe (par rapport à la verticale) */
float DIVTALLE=30;

float numerodurun;            /* Numero du run pour quand lance X rune successifs */

double compteurtemps; /*compteurs de temps à partir emergence feuille 4*/

double Taillepointelineraire=5;         /* mm, pointe extreme ou considere relation lineaire */
double STEPS=120+iter*100;
double dt;
double TPS; 
double SeuilRef[100][400][30];      /* Seuil de ref quand emerge, pour comparaison au cas ou defoli */
int indicetalle[50];                /* Numero des talle en fonction de leur ordre d'apparition et pas d egenealogie, pour ca que var globale */
float Seuil[1000][30];
double Tailfingain[1000][4000][30];   /* Premier talle , second feuille dans talle, troisieme plante */
int NBtalles[50];                     /* Nombre totales de talle VISIBLES par plante */
double A=110;                       /* Taille finale premiere feuille */
double divergetalle;                /*  angle de divergence de la talle par rapport a sa talle-mère */                                    
double Lumiere;

double Surfoliairetotale;        /* Calcul de la surface folaire totale de interpretation */
double STOPTALLAGE=110000;         /*Parametre de surface ou fermeture couvert dans lequel tire le random pour stoper tallage*/
double SENSETALLAGE=10000;

/*----------------------- Groupe Defoliation ----------------------------------------------------*/

double HCOUP=80/*+(5*iter)*/;
double DCOUP;
double PREMIERCOUPE=200;     /*(cortar)*/     /* Date de la premiere coupe */
double FREQUCOUPE=8000;              /* Peridodicité coupe en °C.j */
double DUREECOUPE=600000;             /* Duree sur laquelle effecute coupe */

/*-----------------------------------------------------------------------------------------*/
/*------------------------- Groupe Ecriture -----------------------------------------------*/


double tpsecriture=STEPS-10;
double Groupedebut=1;

double realisationcroiss=1;

double Besointotencroiss;















/*----------------------- Groupe Photosynthese - Turbid medium ----------------------------------------------------*/

double Photosynthese;      /* Tps ou lance groupe de decomposition photosynthese */
double Premphoto=2480000000;           /* Date de la premiere coupe */
double Frequphoto=100;              /* Peridodicité coupe en °C.j */
double Dureephoto=1000;             /* Duree sur laquelle effecute coupe */

float DensiteNRJincident=0.003101; /* densité NRJetique du rayonnement, en KJ par mm2 de sol. Reactiver le sinus dans calcul de Io sur par mm2 de surface orientée vers rayonnement */
const int nbstrates=20;         /* Nb de strates */
const int hstrate=10;          /* Hauteur de chaque strate */
float Largstrate;      /*Arrette des strates - Carré - largeur et longeur */                                
float Longstrate;

#define CAGE Down(90)f(Largstrate)Left(90)F(Longstrate)Left(90)F(Largstrate)Left(90)F(Longstrate)Left(90)F(Largstrate)Left(90)F(Longstrate)

const int nbdirsun=1;          /* Nombre de direction du soleil */
float Elevation[nbdirsun];     /* Elevation du soleil pour calcul interception */
float Azimuthsoleil;           /* Azimth du soleil */

float SURFTOTPARSTRATE[nbstrates];       /* Surf. foliaire totale ds chaque strate */ 
float CLASSINCLINAISON[9][nbstrates+1];  /* Classe d'inclinaison : distribution en 10 classes des inclinaisons  - La strate n+1 correspond au total*/
float Intercepte[nbstrates];
float Bornescage[4];
float Transmis[nbstrates];
float SURFPROJPARSTRATE[nbstrates];





float lumiere[lum_dim_plt][lum_dim_feuil][lum_dim_elt][lum_dim_talle];


/*----------------------- ----------- ----------------------------------------------------*/




struct Feuil
{
  double age;
  float Agecroiss;
  float Taillefeuille;
  int   id_rang;
  int   id_talle;
  int   id_plante;
  float Ymax;
  float Difftps;
  float Reservoir;
  float Taillefinalel;
  float Taillefinaleg;
  float Taillelimbe;
  float Taillegaine;
  double Phase;
  double rapportK;                  /* = rapport gaine / limbe */
  double coupe;                   /* ce qui a été coupé par defoliation*/
  double angleinsert;             /* angle limbe*/
  double geno;
  double Assimibrut;
  double Besoinencroiss;

};

struct Internode            /* Structure de données des Entrenoeuds */
{
 float val;
 int id_rang;
 int id_talle;
 int id_plante;
 float angletal;
 float Assimdescend;
 float Demandescend;
 float Assimmonte;
 float Demandmonte;
};


struct Feuildefoli   /* Structure de données des micro modules C pour defoliation */
{
  double age;
  float Agecroiss;
  float Taillefeuille;
  int   id_rang;
  int   id_talle;
  int   id_plante;
  float Ymax;
  float Difftps;
  float Reservoir;
  float Taillefinalel;
  float Taillefinaleg;
  float Taillelimbe;
  double Taillegaine;
  double Phase;
  double indice;
  double rapportK;
  double status;
  double angleinsert;
  double dejacoupe;              /* ce qui a deja ete coupe par defoli lors evenement precedents */
  double geno;
  double Assimibrut;
  double Besoinencroiss;


};


struct FeuilTrans      /* Structure de données du micro module B qui doit les retransferer à feuil apres defoli */
{
  double age;
  float Agecroiss;
  float Taillefeuille;
  int   id_rang;
  int   id_talle;
  int   id_plante;
  float Ymax;
  float Difftps;
  float Reservoir;
  float Taillefinalel;
  float Taillefinaleg;
  float Taillelimbe;
  double Taillegaine;
  double Phase;
  double longueur;
  double Cutstatus;
  double rapportK;
  double coupe;
  double angleinsert;
  double dejacoupe;
  double geno;
  double Assimibrut;
  double Besoinencroiss;


};



struct FeuilTransPhoto      /* Structure de données du micro module B qui doit les retransferer à feuil apres defoli */
{
  double age;
  float Agecroiss;
  float Taillefeuille;
  int   id_rang;
  int   id_talle;
  int   id_plante;
  float Ymax;
  float Difftps;
  float Reservoir;
  float Taillefinalel;
  float Taillefinaleg;
  float Taillelimbe;
  double Taillegaine;
  double Phase;
  double longueur;
  double Cutstatus;
  double rapportK;
  double coupe;
  double angleinsert;
  double dejacoupe;
  double geno;
  double Assimibrut;
  double Besoinencroiss;
}; 
   
struct Feuilphoto   /* Structure de données des micro modules Cp pour photosynth */
{
  double age;
  float Agecroiss;
  float Taillefeuille;
  int   id_rang;
  int   id_talle;
  int   id_plante;
  float Ymax;
  float Difftps;
  float Reservoir;
  float Taillefinalel;
  float Taillefinaleg;
  float Taillelimbe;
  double Taillegaine;
  double Phase;
  double indice;
  double rapportK;
  double status;
  double angleinsert;
  double dejacoupe;              /* ce qui a deja ete coupe par defoli lors evenement precedents */
  double geno;
  double Assimibrut;
  float SStratPropr[nbstrates];          /* Longueur des feuilles repartie dans chacune des strates */
  float SprojPropr[nbstrates];            /* Surf. projetée de morceau de feuille dans chaque strates */
  double Besoinencroiss;

};






/************************************************ FONCTIONS *****************************************************************/
/****************************************************************************************************************************/


/* Fonc de Croissance */

float FoncCroiss(float t, float Ymax, float te, float tm) 
{ 
 float tend=te; 
 float tmax=tm; 
 float y;

    if (t>tend) 
    {
    y=Ymax;
    }
    else
    {
     y=Ymax*(1+((tend-t)/(tend-tmax)))*(pow((t/tend),(tend/(tend-tmax))));
    }
return y;
}

 
/* Solveur d'equation qui permet inverser fonction beta pour avoir T à partir d ela taille */


float TpourY(float Y, float Ym, float Te, float Tm)                              
{
float T; 
float H;
for (int i=0; i<800; i=i+1)
  {
  T=FoncCroiss(i, Ym, Te, Tm); 
   if (T<Y) 
   H=i;
  }
return H;
}

/* Ci dessous : fonction de determination du rang de la feuille precedente la plus haute */


int rangduplushaut(double X[1000][4000][30], int dim2, int dim3, int rangpropre)
{
double grand=0;
int rangrand=7;
double y;
 
if (rangpropre>4)
{
 for (int i=1; i<4; i++)
  {
   y=X[dim2][rangpropre-i][dim3];      
   if (y>grand) {grand=y;rangrand=i;} 
  }
  return rangrand;
}
else return 1;
}




/****************************************** STARTERS ET INITIALISATION *******************************************************/

ignore:RollR Down SetColor MoveTo Left;

Start:
{ 
	for (int i=0;i<lum_dim_plt; i++){for(int j=0; j<lum_dim_feuil; j++){for(int k=0; k<lum_dim_elt; k++){for(int l=0; l<lum_dim_talle; l++){lumiere[i][j][k][l]=0;}}}}
	#include LUMIEREINIT; 







  sran(4509);
  
  
  DCOUP=PREMIERCOUPE;
  Photosynthese=Premphoto;
  for (int i=0; i<1000; i++){for(int j=0; j<30; j++){Seuil[i][j]=5.5;}};
  for (int i=0;i<1000; i++){for(int j=0; j<4000; j++){for(int k=0; k<30; k++){Tailfingain[i][j][k]=5.5;}}};
  for (int i=0; i<50; i++){indicetalle[i]=1;};
  for (int i=0;i<100; i++){for(int j=0; j<400; j++){for(int k=0; k<30; k++){SeuilRef[i][j][k]=5.5;}}};
  for (int i=0; i<50; i++){NBtalles[i]=0;};
  for (int i=0; i<nbstrates; i++){SURFTOTPARSTRATE[i]=0;};
  for (int i=0; i<9; i++){for(int j=0; j<nbstrates+1; j++){CLASSINCLINAISON[i][j]=0;}};
  for (int i=1; i<5; i++){Bornescage[i]=0;};
  for (int i=0; i<nbstrates; i++){Transmis[i]=1000;};

  



  divergetalle=DIVTALLE;
  TPS=0;
  Biomasse=0; 
  sommlineairecoupee=0;
  Lumiere=44;
  UseGroup(Groupedebut);
  dt=1;
  
 compteurtemps=-80;
 
 fps=fopen(SORTIE, "w");
 fprintf(fps, "TPS \t Surf \n"); 
 
 
  numerodurun=numerodurun+1;


}

StartEach:
{ 



Besointotencroiss=0;



if(TPS==tpsecriture)
{
 fpss=fopen(CHAINE, "w");
   fpsss=fopen(GEOMETRIE, "w");
   fprintf(fpsss,"Gain/Limb \t Hauteur \t Diamgain \t id-plant \t id-talle \t id-feuil \t id-indice \t XPos \t YPos \t ZPos \t XHead \t YHead \t ZHead \t XUp \t YUp \t ZUp \n");
//   fprintf(fpsss,"Gain/Limb \t Hauteur \t Diamgain \t id-plant \t id-talle \t id-feuil \t XPos \t YPos \t ZPos \t XHead \t YHead \t ZHead \t XUp \t YUp \t ZUp \n");

}
if(TPS==tpsecriture+2*dt)
{
 fclose(fpss);
 fclose(fpsss);
}


/*
for (int i=0; i<30; i++)
   {
    if (TPS==21*i)  
    {
      fprintf(fps,"%i \t", NBtalles[1]);
    }
   }

if (TPS==660)  
   {
     fprintf(fps,"%i \n", NBtalles[1]);
   }
*/ 
 
   /*Printf("Transmis Strat bas: %f\n",Transmis[0]);
   Printf("numero du run: %f\n",numerodurun);*/

  if (TPS>DCOUP+10 && TPS<(PREMIERCOUPE+DUREECOUPE)) {DCOUP=DCOUP+FREQUCOUPE;}
  if (TPS>Photosynthese+10 && TPS<(Premphoto+Dureephoto)) {Photosynthese=Photosynthese+Frequphoto;}

  if (TPS>0 && TPS<50) { Lumiere=43;}
  if (TPS>=50 && TPS<600) { Lumiere=44;}
  if (TPS>600  && TPS<700) { Lumiere=43;}
  if (TPS>=700) { Lumiere=44;}

  TPS=TPS+dt;
  Printf("TPS: %f\n",TPS);
  if (compteurtemps>=-60) {compteurtemps=compteurtemps+dt;}
 


 if (TPS==Photosynthese)
  {

  for (int i=0; i<nbstrates; i++){SURFTOTPARSTRATE[i]=0;};
  for (int i=0; i<nbstrates; i++){SURFPROJPARSTRATE[i]=0;};
  for (int i=0; i<9; i++){for(int j=0; j<nbstrates+1; j++){CLASSINCLINAISON[i][j]=0;}};
  for (int i=0; i<nbdirsun; i++){Elevation[i]=/*((i*10)+(i+1)*10)/2*/56.72;};
  for (int i=1; i<5; i++){Bornescage[i]=0;};
  for (int i=0; i<nbstrates; i++){Transmis[i]=1000;};

  Azimuthsoleil=-57.3;


   switch(CurrentGroup())
   {
    case Croissance: 
    UseGroup(Photosynth);
    Forward();
    dt=0.125;
    break;
   }
  }

  if (TPS==Photosynthese+1)
  {
   switch(CurrentGroup())
   {
    case Photosynth:
    
    UseGroup(FOLD);
    Backward();
    dt=0.125;
    break;
   }
  }

if (TPS==Photosynthese+2)
  {
   switch(CurrentGroup())
   {
    case FOLD:
    
    UseGroup(ALLOCATION);
    Forward();
    dt=0.25;
    break;
   }
  }
  if (TPS==Photosynthese+3)
  {
   switch(CurrentGroup())
   {
    case ALLOCATION:
    UseGroup(Croissance);
    Forward();
    dt=1;
    break;
   }
  }

  if (TPS==DCOUP)
  {
   switch(CurrentGroup())
   {
    case Croissance: 
    UseGroup(Defoliation);
    Forward();
    dt=0.125;
    break;
   }
  }
  if (TPS==DCOUP+1)
  {
   switch(CurrentGroup())
   {
    case Defoliation:
    UseGroup(Croissance);
    Forward();
    dt=1;
    break;
   }
  }

   if (TPS==tpsecriture)
  {
   switch(CurrentGroup())
   {
    case Croissance:
    
    UseGroup(DECOMPOECRIT);
    Forward();
    dt=0.5;
    break;
   }
  }
  if (TPS==tpsecriture+1)
  {
   switch(CurrentGroup())
   {
    case DECOMPOECRIT:
    
    UseGroup(ECRITURE);
    Forward();
    dt=1;
    break;
   }
  }







if (TPS==5)
  {
   switch(CurrentGroup())
   {
    case LECTURE:
    
    UseGroup(Croissance);
    Backward();
    dt=1;
    break;
   }
  }

if (TPS==109)
  {
   switch(CurrentGroup())
   {
    case FOLD:
    UseGroup(ALLOCATION);
    Forward();
    dt=1;
    break;
   }
  }
if (TPS==159)
  {
   switch(CurrentGroup())
   {
    case ALLOCATION:
    
    UseGroup(Croissance);
    Forward();
    dt=1;
    break;
   }
  }




}

EndEach:
{



Besointotencroiss=0;


for (int i=0; i<100; i++)
   {
    if (TPS==20*i)  
    {
      fprintf(fps,"%f\t %f\n",TPS, Surfoliairetotale);
    }
   }







if (TPS==Photosynthese+0.75)
   {
         
     
       for(int i=0; i<nbstrates; i++) 
        {
         Printf("Surface proj entre %i et %i : %f\n\n", (i*hstrate), ((i+1)*hstrate),SURFPROJPARSTRATE[i]);
        } 
	
       /*


        
          for(int i=0; i<9; i++) 
        {
         Printf("Longeur entre %i et %i : %f\n", i*10, (i+1)*10 ,CLASSINCLINAISON[i][3]);
        } */
   }
}

End:
{ 
fclose(fps);
}

/****************************************************************************************************************************/

derivation length: STEPS;

/***************************************************************************************************************************/


module Feuille(Feuil);
module Entrenoeud(Internode);
module Apex(int, int,int,int);  /* rang, talle, plante,geno*/
module Talle(int,int,int);      /* Num de talle et de plante,geno */          
module ApexTal(int,int);        /* Second : geno */
module Gaine(double, double, int, int); /* Longueur et angle tal transmis à sig, rang et talle*/
module Centralerect(float, float, float, float, float); /* Gestion erectness des talle:  age, hauteur etape avant, status pour savoir si dois faire up ou down, angle, geno  */
module Neon(float);             /* Neon et Defoli pour defoliation */
module Affichage(float);        /* Gestion affichage commentaires etc */
module B(FeuilTrans);           /* B et C pour phase de defoliation */
module C(Feuildefoli);
module Marqueur(float);         /* module utilisé pour marqué le bout de la feuille apres defoli */
module Graine(int, float, float, int); /*age, num talle, plante, geno pour heterogeneite germination population */
module Cp(Feuilphoto);          /* Module comme le C mais pour Photosynth */
module Pool(float, float);             /* Pool de Carbone */
module Racine(float);           /* Racines */
module Reserve(float);          /* Reserve en C/equivalent-biomasse */







module Planteur(int, int);      /* Module "planteur" : genere une population avec X plantes en ligne et Y en colonnes */
module Soleil(float, float);     /* Pour représentation de position soleil dans voute celeste */


/*---------pour ecriture ---------*/

module Cdecomp(double,float, float, int, int, int,float, float, float, float, float, float,double, double,double,double,double,double,double,double, double,double); /* Module C de decomposition attention : sans les deux dernier parametre qui sont tableau */
module Bdecomp(double,float,float, int,int,int,float, float,float,float, float, float,double,double,double,double,double,double,double,double,double,double,double); /* Module B de decomposition */
module Entrendecomp(float,int,int,int,float,float,float, float, float);



/***************************************************************************************************************************/


Axiom: 
/*
SB() Soleil(0,0) EB()

SB() Neon(0) EB()
SB()MoveTo(0, 0, 0)SetColor(21)SetWidth(2)GROUND EB()
SB()MoveTo(GROUNDSIZE+0, 0, 0)SetColor(21)SetWidth(2)GROUND EB()*/
Planteur(2,2)

/*SB() Affichage(0) EB()

Graine(0,0,1,1)*/





;

/************************************************ CONTEXT SENSITIVE RULES *************************************************/
/**************************************************************************************************************************/


/*

Gaine(longueur, angletal, rang, talle) > SB() GetPos(X, Y, Z):
{
  if (Y<0.8)
  {
   produce 
   Gaine(longueur,0, rang, talle);
  } 
if (Y>0.8)
  {
   produce 
   Gaine(longueur,divergetalle, rang, talle);
  } 
}
*/

/*

Entrenoeud(X) > SB() Gaine(longueur, angletalle, rang, talle)  :
{   
    X.angletal==angletalle;				
    produce Entrenoeud(X);
}
*/

/*

Centralerect(age, rang, Rot, angle,geno) > GetHead(X, Y, Z) :
{ 

if (Lumiere==43)
{ 
  if (Y>0.2)
  {
  produce Centralerect(age+dt, rang, Rot, angle+2,geno);
  }

  if (Y<0.1)
  {
  produce Centralerect(age+dt, rang, Rot, angle-2,geno);
  }

  else
  {
  produce Centralerect(age+dt, rang, Rot, angle,geno);
  }
}
if (Lumiere==44)
{ 
  if (Y<0.9 && Y>0.1)
  {
  produce Centralerect(age+dt, rang, Rot, angle-2,geno);
  }
  if (Y<0.1)
  {
  produce Centralerect(age+dt, rang, Rot, (angle-2)*(-1),geno);
  }
  else
  {
  produce Centralerect(age+dt, rang, Rot, angle,geno);
  }
}
}
  
*/





C(Z) > SB() GetPos(X, Y, Za):

{
  if (Y>HCOUP)
  {
   produce C(Z) SB() Marqueur(0) EB()Cut();
  }
  else
  {
   produce C(Z);
  }
}

C(Z) > SB() Marqueur(onsenfout):
{ 
  FeuilTrans D={Z.age, Z.Agecroiss, Z.Taillefeuille,Z.id_rang, Z.id_talle,Z.id_plante, Z.Ymax, Z.Difftps, Z.Reservoir, Z.Taillefinalel, Z.Taillefinaleg, Z.Taillelimbe, Z.Taillegaine, Z.Phase, Z.Taillefeuille-(Z.dejacoupe+(Z.indice*LSTEP)),1,Z.rapportK,(Z.indice*LSTEP),Z.angleinsert, Z.dejacoupe, Z.geno, Z.Assimibrut, Z.Besoinencroiss};
  produce B(D);
}









/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------    GROUPE CROISSANCE --------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/


group Croissance : 



Entrenoeud(X) > SB()Feuille(Y) :
{ 
  if (abs(Seuil[Y.id_talle][Y.id_plante]-Y.Taillefeuille)<2 )
   {  
    X.val=1;				
    produce Entrenoeud(X);
   }
}

Entrenoeud(X1) > Entrenoeud(X2)  :
{ 
  if (X2.val==1)
   {   	
      X1.val=2;	
      produce Entrenoeud(X1);
   }
}

Entrenoeud(X) < Apex(id_rang, id_talle, id_plante, id_geno)   :
{ 
  if (X.val==1)
   {   
      float Premiecroiss;
       if(id_geno==1)
        {
         Premiecroiss=A;
        }
       if(id_geno==2)
        {
         Premiecroiss=A;
        } 			    
      Feuil Y={0,Premiecroiss,0,id_rang,id_talle,id_plante,Premiecroiss,0,0,0.8*Premiecroiss,0.2*Premiecroiss,0,0,0,0.2,0,0,id_geno,0,0}; 
      Internode Creat={0,id_rang,id_talle, id_plante, X.angletal,0,0,0,0};
      produce RollR(180)Entrenoeud(Creat)SB()/*MoveTo(id_rang*40, 0, id_talle*40)*/Feuille(Y)EB()SB()ApexTal(0,id_geno)EB()Apex(id_rang+1, id_talle,id_plante,id_geno);
   }
}

/*                  RollR(180) comme ca tout le phytomere tourne a chaque nouvelle feuille   */

Entrenoeud(X) < ApexTal(rang,id_geno) : /* Declenchement des bourgeons axillaires */
{
  if (X.val==2 && id_geno<6)
    { 
      float Sitefillingproba=ran(STOPTALLAGE);                     
      
      if (Sitefillingproba>Surfoliairetotale)       /* Pour la gestoin du site filling : plus la surface augmente et plus la probabilité de faire un talle diminue */
      { 

         
         if (Surfoliairetotale<=SENSETALLAGE)
         { 
           indicetalle[X.id_plante]=indicetalle[X.id_plante]+1; 
           if (indicetalle[X.id_plante]>=2)
            {  
            Seuil[indicetalle[X.id_plante]-1][X.id_plante]=Tailfingain[X.id_talle][X.id_rang-1][X.id_plante];
            produce Right(0)Talle(indicetalle[X.id_plante]-1,X.id_plante,id_geno);
            }
         }

         if (Surfoliairetotale>SENSETALLAGE)
         { 
           float deuxiemede=ran(5);
           if (deuxiemede>=2)
            {
              indicetalle[X.id_plante]=indicetalle[X.id_plante]+1; 
              if (indicetalle[X.id_plante]>=2)
              {  
               Seuil[indicetalle[X.id_plante]-1][X.id_plante]=Tailfingain[X.id_talle][X.id_rang-1][X.id_plante];
               produce Right(0)Talle(indicetalle[X.id_plante]-1,X.id_plante,id_geno);
              }
            }
           else
            {
             produce /*F(0)*/ ApexTal(rang,6); /* Bourgeon avorté/dormant.Pour que si jamais la proba. ne declenchais pas une talle il ne reessaye pas à chaque step. */
            }
          }


       }
       
      else
      {
       produce /*F(0)*/ ApexTal(rang,6) ; /* Bourgeon avorté.Pour que si jamais la proba. ne declenche pas une talle il ne reessaye pas à chaque step. */
      }
    }
}


Entrenoeud(X) > SB() Gaine(longueurg, angletalg, rangg, talleg)Feuille(Y)EB()SB()Right(angleencours):
{
 X.angletal=angleencours;
 produce Entrenoeud(X);
}

Gaine(longueur, angletal, rang, talle) > Feuille(Y)EB()Entrenoeud(X):
{    
  produce Gaine(Y.Taillegaine, X.angletal, rang, talle);
}
 
Left(t)>Gaine(longueur, angletal, rang, talle):
{
 if(angletal>0)
 {
 produce Left(angletal);
 }
}


Right(angletal)>Centralerect(ageC, rang, Rot, angle,geno):
{
 if (ageC>=120 && angletal<=divergetalle/*geno*/)
 {
  produce Right(angletal+(dt/3));
 }
 else
 {
  produce Right(angletal);
 }
}

/* Considere le Right des talles comme un module. Comme ca fait varier sa valeur */
/* ICi fait varier en fonction de l'AGE de la TALLE */
/* Ci dessous transfert de l'info de l'angle de la talle pour ecartement simultané de sa feuille engainante */







/*--------------------------------------------CROISSANCE PRODUCTION --------------------------------------------------------*/


production:






Racine(temps):{produce Racine(temps+dt);}







/* Je remplace mon module B qui se souvient de la taille des feuilles restantes par des vraies feuilles qui reprennent leur croissance ou elle
en etaient */

B(D):
{ 

 if (TPS==DCOUP+dt)
 {
  sommlineairecoupee=sommlineairecoupee+(D.coupe);  
  Biomasse=sommlineairecoupee*Densitefeuille;
 }

 /* Calcule du cumul des tissus coupé, fait le calcul dans B car avec les X.coupe de feuille je serais baisé car les feuilles 
qui sont en cutstatus 2 sont pas remplacée - Sinon pas de soucis avec TPS+dt, il le fait à Dcoup+1 car on est dans groupe croissance
donc pas acces à dt==0.125. Mais sans probleme, c'est bien comme ca) */
 
 
  Feuil P={D.age,D.Agecroiss,D.Taillefeuille,D.id_rang,D.id_talle,D.id_plante,D.Ymax,D.Difftps,D.Reservoir,D.Taillefinalel,D.Taillefinaleg,D.Taillelimbe,D.Taillegaine, D.Phase,D.rapportK, D.dejacoupe+D.coupe, D.angleinsert, D.geno, D.Assimibrut, D.Besoinencroiss};
 
  if (D.Cutstatus==0)
  { 
   if (D.Phase==1)
    {
    produce Gaine(D.Taillegaine, 0, D.id_rang, D.id_talle)SB()GetPos(0, 0, 0)EB()Feuille(P);
    }
   if (D.Phase==0)
    {
    produce Feuille(P);
    }
  }
 
  if (D.Cutstatus==1)
  { 
    if (D.Phase==1)
    {
     if (D.Taillelimbe<0)
      {
       
       P.Taillegaine=D.longueur; /*sinon P.taillegaine a en memeoire la longeur avant coupe et ca fera un etape ou gaine aura longueur incorrcte */
       produce Gaine(D.longueur, divergetalle, D.id_rang, D.id_talle)SB()GetPos(0, 0, 0)EB()Feuille(P); /* pour le cas precis ou gaine coupée dans sa croissance */
      } 
    produce Gaine(D.Taillegaine, 0, D.id_rang, D.id_talle)SB()GetPos(0, 0, 0)EB()Feuille(P);
    }
   if (D.Phase==0)
    {
    produce Feuille(P);
    }
  }
  if (D.Cutstatus==2)
  { 
   Feuil hh={0,0,0,D.id_rang,D.id_talle,D.id_plante,0,D.Difftps,0,0,0,0,0, D.Phase,D.rapportK, D.dejacoupe+D.coupe, D.angleinsert, D.geno, D.Assimibrut, D.Besoinencroiss};
   produce Gaine(D.Taillegaine, 0, D.id_rang, D.id_talle)Feuille(hh);
  } 
}



C(Z): {produce ;}
Cp(Z): {produce ;}

Graine(age, Id_talle, Id_plante, Id_geno):
{
 if (age==1)
   {
    produce Racine(0) Pool(0,0) SB()Reserve(0)EB() Talle(Id_talle, Id_plante, Id_geno);
   }
 else
   {
    produce Graine(age+dt, Id_talle, Id_plante, Id_geno);
   }
}     


Feuille(X):
{  

   
   double LEDGAINE; 
   
  if (X.Taillefeuille<X.Ymax+10)
  {
    X.age=X.age+dt;
  } 
   
  if(X.angleinsert<=INSERTIONLIMBE && X.Taillegaine>Tailfingain[X.id_talle][X.id_rang-1][X.id_plante]+1)
  {
    X.angleinsert=/*X.angleinsert+1*/ INSERTIONLIMBE;/* pour avoir limbe qui tombe pas brusquement */
  }

     

 /*if ((X.id_rang==2 || X.id_rang==6)  && X.age==100) 
  {
      Seuil[X.id_talle][X.id_plante]=Seuil[X.id_talle][X.id_plante]*1.3; 
  } */
 /* Coeff 1.16 pour faire en sorte que la feuille 2 qui a meme seuil donne pas meme taille*/
  


  if ((X.Taillegaine)>Seuil[X.id_talle][X.id_plante]) 
  {
      Seuil[X.id_talle][X.id_plante]=X.Taillegaine;
  } 
  if (X.Taillegaine>Tailfingain[X.id_talle][X.id_rang][X.id_plante]) 
  {
     Tailfingain[X.id_talle][X.id_rang][X.id_plante]=X.Taillegaine;
  } 
  

  if (X.Phase==0)
  {

    if(X.geno==1)   /* Pour FL */
    {
     

      X.Reservoir=1*((pow(2.718,(X.age+102)/31.5)-30)-(pow(2.718,(X.age+19.82)/20.99)));
    







     X.rapportK=0.1+(X.age/600);   /* Pour que rapport gaine/limbe varie fonction temps passé dans gaine */
     X.Agecroiss=X.age;
    }












    if(X.geno==2)  /* Pour FC  */
    {
     //X.Reservoir=0.55*((615*pow(2.718,X.age*0.0053))-736.9-(0.1207*pow(2.718,X.age*0.062)));
     
     X.Reservoir=0.75*((pow(2.718,(X.age+102)/31.5)-30)-(pow(2.718,(X.age+19.82)/20.99)));
     X.rapportK=0.1+(X.age/600);   /* Pour que rapport gaine/limbe varie fonction temps passé dans gaine */
     X.Agecroiss=X.age;
    }
 }



  if (X.Phase==1)
  {
    X.Ymax=X.Reservoir;
    X.Taillefinalel=(1-X.rapportK)*X.Reservoir;
    X.Taillefinaleg=X.Reservoir*X.rapportK;
    X.Agecroiss=X.age+X.Difftps;
    
  }
  

if(X.geno==1)   /* Pour FC */
    {
     LEDGAINE=((X.Taillefinaleg*5)*2)-X.Ymax;  
    }
if(X.geno==2)  /* Pour FL  */
    {
     LEDGAINE=((X.Taillefinaleg*5)*2)-X.Ymax;  
    }

  
   double Accroiss;

  if ((X.Taillefeuille>Seuil[X.id_talle][X.id_plante]-10) && (X.Taillefeuille<Seuil[X.id_talle][X.id_plante]) && X.Phase==0)
  {
    X.Difftps=TpourY(Seuil[X.id_talle][X.id_plante], X.Reservoir,1.35*(((X.Reservoir/0.2)*(0.2+((X.age/(60*10))-0.1))*2)-X.Reservoir),0.76*1.35*(((X.Reservoir/0.2)*(0.2+((X.age/(60*10))-0.1))*2)-X.Reservoir))-X.age;
    
  }

   
 if (X.Taillegaine<(Tailfingain[X.id_talle][X.id_rang-1][X.id_plante]+13))
  { 
   

   float inversea=TpourY(X.Taillefeuille, X.Ymax, 1.35*LEDGAINE,0.76*1.35*LEDGAINE);
  
    if(X.Phase==0)
    {
    Accroiss=FoncCroiss(X.Agecroiss, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE)-FoncCroiss(X.Agecroiss-1, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE);
    }
    if(X.Phase==1)
    {
     Accroiss=realisationcroiss*(FoncCroiss(inversea+dt, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE)-FoncCroiss(inversea, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE));
    }
    X.Taillefeuille=X.Taillefeuille+Accroiss;
    //if(X.id_rang==1){Printf("Accroiss: %f  Lfeuil: %f\n", Accroiss, X.Reservoir);}              
 }

 

if (X.Taillegaine>0.98*X.Taillefinaleg && X.Taillegaine<=Tailfingain[X.id_talle][X.id_rang-1][X.id_plante]-5 && X.Phase==1)
  {
   Accroiss=1;
   X.Taillefeuille=X.Taillefeuille+1;
  }

  if (X.Taillegaine>=(Tailfingain[X.id_talle][X.id_rang-1][X.id_plante]+13) && X.Phase==1)
  {
   Accroiss=0;
   X.Taillefeuille=X.Taillefeuille;
  }
  


if(TPS==tpsecriture-1)
{
 float inverse=TpourY(X.Taillefeuille,X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE);
 X.Besoinencroiss=Densitefeuille*Accroiss;//(FoncCroiss(inverse+dt, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE)-FoncCroiss(inverse, X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE));
 
}
/*if(TPS==STEPS-1)
{
  Printf("Besoin en croiss masse : %f\n", X.Besoinencroiss);
}*/





float inverse=TpourY(X.Taillefeuille,X.Ymax, 1.35*LEDGAINE, 0.76*1.35*LEDGAINE);
Besointotencroiss=Besointotencroiss+(Densitefeuille*Accroiss);


   
  
 






  double R=0;

  if ((X.Taillelimbe+X.coupe)<(0.75*X.Taillefinalel))
  {
    R=0.973;
  }
  if (((X.Taillelimbe+X.coupe)>(0.75*X.Taillefinalel)) && ((X.Taillelimbe+X.coupe)<X.Taillefinalel))
  {
    R=0.965;
  }
  
  X.Taillelimbe=X.Taillelimbe+Accroiss*R;

  if (X.Taillelimbe>=0)
  {
   X.Taillegaine=X.Taillefeuille-(X.Taillelimbe+X.coupe);
   }
  if (X.Taillelimbe<0)
  {
   X.Taillegaine=X.Taillefeuille-X.coupe;
  }
  




 
    if (((X.Taillelimbe+X.Taillegaine)>Seuil[X.id_talle][X.id_plante]) && X.Phase==0)
   {  
    /* Mis X.tgaine+X.tlimbe, comme si si coupé en croissance c'est taille reelle qui detecte */ 
     
     SeuilRef[X.id_talle][X.id_rang][X.id_plante]=Seuil[X.id_talle][X.id_plante];
     Seuil[X.id_talle][X.id_plante]=Seuil[X.id_talle][X.id_plante]*1.15;
     X.Phase=1;
     
     if (X.id_rang==4)
     {
      compteurtemps=-51;
     }

     if (X.id_rang==1)
     {
      NBtalles[X.id_plante]=NBtalles[X.id_plante]+1;
     }
     produce Left(0)Gaine(X.Taillegaine, 0, X.id_rang, X.id_talle)SB()GetPos(0, 0, 0)EB()Feuille(X);
   }














   double Agesenesc=400+(X.Ymax/20)*30; 
   if (X.age>Agesenesc+100)
   {                                                       
    produce Cut();                                                /* Pour que la feuille soit vraiment coupée supprimee de chaine, en ENTIER */
   }
   else 
   {
    produce Feuille(X);
   }
}

Centralerect(age, rang, Rot, angle,geno):
{
  produce Centralerect(age+dt, rang, Rot, angle,geno);
}




/*--------------------------------------------CROISSANCE DECOMPOSITION -----------------------------------------------------*/



decomposition:

Planteur(Nbligne, Nbcolonne):   /* Planteur : si population aléatoire de genotype faire tirage ran(2) avec 1 et 2 qui sont geno. */
{ int numplant;
  numplant=1;

  for (int i=0; i<=Nbligne; i++)
  {
   for (int j=0; j<=Nbcolonne; j++)
    {
     double numgeno=1/*int(ran(2))*/;
     if(numgeno==0){numgeno=1;};
     nproduce MoveTo(-50+i*50,0,-50+j*50) Graine(0,0,numplant,numgeno);
     numplant=numplant+1;
    }
   }
   produce;
}


Talle(rang,plante,geno):
{  
 float Premiecroiss;
 if(geno==1)
 {
 Premiecroiss=A;
 }
 if(geno==2)
 {
 Premiecroiss=A;
 }
 Feuil X={0,Premiecroiss,0,1, rang,plante,Premiecroiss,0,0,Premiecroiss*0.8,Premiecroiss*0.2,0,0,0,0.2,0,0,geno,0,0};
    if (rang>=0)
     {
        double Rotation=20+ran(300);
        Internode Creat={0,1,rang,plante,0,0,0,0,0};
        produce 

        /*MoveTo(-GROUNDSIZE/2*rang, 0,GROUNDSIZE/2)*/RollR(Rotation)Centralerect(0,rang,Rotation,1,geno)GetHead(0, 0, 0)
        
        Entrenoeud(Creat)
        SB()Feuille(X)EB()
        SB() ApexTal(0,geno) EB()
        Apex(2, rang,plante,geno);  
       
      }
}

/*--------------------------------------------CROISSANCE  INTERPRETATION -----------------------------------------------------*/


interpretation:



float couleur=2;


Feuille(X):
{
double surfacefeuille=0;

if (X.id_talle==0 && X.id_rang==1){Surfoliairetotale=0;}
/* Comme ca remet à zero le calcul de surfacetotale à chaque interpretation quand arrive 
a premiere feuille. Malgré tout il continue à accumuler toutes les feuilles derriere !!! */



if (X.Phase==0)
   { 
    double numsteps=(X.Taillefeuille-X.coupe)/LSTEP;
      for (int i=0; i<numsteps; i++)
       {
         double large;           /* Largeur à chaque i */
         double L0;
         double taille=X.Taillefeuille;
         L0=X.Taillefinalel-taille;                            /* Point de depart de parcours de la feuille */
         
         if (i<((taille-Taillepointelineraire)/LSTEP))
         {
          large=(-1.38*pow(((L0+(i*LSTEP))/X.Taillefinalel), 2)+0.66*((L0+(i*LSTEP))/X.Taillefinalel)+0.903)*(((1.9645*log(double(X.Ymax)))-5.701)/2);  
          /*Derniere expression remplace largeur max, comme ca largeur max fontion de longeur ! */
         }
         else
         {
         double largini=(-1.38*pow(((L0+((i-1)*LSTEP))/X.Taillefinalel), 2)+0.66*((L0+((i-1)*LSTEP))/X.Taillefinalel)+0.903)
         *(((1.9645*log(double(X.Ymax)))-5.701)/2);
         large=largini-((i-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire));
         }
      
         if(i<(10/LSTEP))
         {
          double rank=X.id_rang;
          nproduce  StartGC() CurrentContour(Limbo) SetWidth(large)SetColor(couleur)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
         }
         if(i>(10/LSTEP))
         {        
         nproduce StartGC() CurrentContour(Limbo)SetWidth(large)SetColor(couleur)F(LSTEP) ;
         nproduce EndGC();
         }
       }
     produce;
    }

if (X.Phase==1) 
   {
    int numsteps=ceil((X.Taillefeuille-X.coupe)/LSTEP);
    int numstepsg=ceil((X.Taillegaine)/LSTEP);
    for (int i=0; i<numsteps; i++)
    {
        
       /*if (i==1)
       {
        nproduce SetWidth(3)SetColor(couleur) F(LSTEP);
       }
       if (i<numstepsg)
       {
        nproduce SetWidth(3)SetColor(couleur) F(LSTEP); 
       }
       if (i==numstepsg-1 && X.Taillegaine>Tailfingain[X.id_talle][X.id_rang-1][X.id_plante]+1)
       {
        nproduce StartGC() CurrentContour(Gaino)SetWidth(3)SetColor(couleur)F(LSTEP); 
       }*/

       if (i>=numstepsg-1) 
       {

       double large;           /* Largeur à chaque i */
       double L0;
       double taille=X.Taillelimbe+X.coupe;

       L0=X.Taillefinalel-taille;  /* Point de depart de parcours de la feuille */ 
                         
          if ((i-numstepsg)<((taille-Taillepointelineraire)/LSTEP))
          {
           large=(-1.38*pow(((L0+((i-numstepsg)*LSTEP))/X.Taillefinalel), 2)+0.66*((L0+((i-numstepsg)*LSTEP))/X.Taillefinalel)+0.903)*(((1.9645*log(double(X.Ymax)))-5.701)/2);
          }
          else
          {
           double largini=(-1.38*pow(((L0+(((i-numstepsg)-1)*LSTEP))/X.Taillefinalel), 2)+0.66*((L0+(((i-numstepsg)-1)*LSTEP))/X.Taillefinalel)+0.903)*(((1.9645*log(double(X.Ymax)))-5.701)/2);
           large=largini-(((i-numstepsg)-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
          }
          if(i<numstepsg-1) {large=0;}

          surfacefeuille=surfacefeuille+large*2*LSTEP;
          Surfoliairetotale=Surfoliairetotale+large*2*LSTEP; 
          /* Facteur deux dans calcul surface car le contour fait deux de large, donc sa largeur réelle est sous estimée expres  dans le calcul de "large" destiné au setwidth */
          


         double Agesenesc=400+(X.Ymax/20)*30; 

         if(i==numstepsg) 
           
         {
           int precedentplusgrand=rangduplushaut(Tailfingain, X.id_talle, X.id_plante, X.id_rang);
           if(X.Taillegaine>Tailfingain[X.id_talle][X.id_rang-precedentplusgrand][X.id_plante])
           {
           nproduce StartGC() CurrentContour(Limbo)Left(X.angleinsert)SetColor(couleur)F(LSTEP)EndGC();
           }
         }
         
         if (X.age<=Agesenesc)
         {
           if ((i>numsteps-2) && X.coupe>2)          /* Pour avoir cicatrice au bout feuille coupées */
            {
              nproduce StartGC() CurrentContour(Limbo) SetWidth(large)Left(curve)SetColor(couleur)F(LSTEP)EndGC();
          
            }
            
           if(X.id_rang==1 && X.id_talle>0)
            {  
               if (i>(SeuilRef[X.id_talle][X.id_rang][X.id_plante]-2)/LSTEP)
                {
                  nproduce StartGC() CurrentContour(Limbo) SetWidth(large)RollR(/*(X.geno-1)*3*/0)Left(curve)SetColor(couleur)F(LSTEP)EndGC();
                }  
               if (i<(SeuilRef[X.id_talle][X.id_rang][X.id_plante]-2)/LSTEP)
                {
                  double rank=X.id_rang;
               
                 if(i<(10/LSTEP))
                  {
                   nproduce  StartGC() CurrentContour(FoldLimbo) SetWidth(large)SetColor(couleur)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                  }
                  else
                  {
	          nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(couleur)F(LSTEP)EndGC();
                  }
                }
             }      
            else
            {

              int precedentplusgrand=rangduplushaut(Tailfingain, X.id_talle, X.id_plante, X.id_rang);

              if (i>(Tailfingain[X.id_talle][X.id_rang-precedentplusgrand][X.id_plante]/LSTEP))
              {
               nproduce StartGC() CurrentContour(Limbo) SetWidth(large)RollR(/*(X.geno-1)*3*/0)Left(curve)SetColor(couleur)F(LSTEP)EndGC();
              }  
              if (i<(Tailfingain[X.id_talle][X.id_rang-precedentplusgrand][X.id_plante]/LSTEP))
              {
               double rank=X.id_rang;
               
                 if(i<(10/LSTEP))
                  {
                   nproduce  StartGC() CurrentContour(FoldLimbo) SetWidth(large)SetColor(couleur)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                  }
                  else
                  {
	          nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC();
                  }
              }                        /* Pour avoir feuille courbe et "s'ouvre" que quand en dehors de la gaine precendete la plus grande */
            }
         }
         else 
         if (X.age>Agesenesc && X.age<=Agesenesc+50)
         {
          nproduce StartGC() CurrentContour(Limbo) SetWidth(large)Left(curve)SetColor(couleur)F(LSTEP)EndGC();
          /* Pour le brunissement avant mort, cut dans interpet sert a rien, mis dans prod. */ 
         }

         
         


     } 
  }
static char rang[80];
sprintf(rang,"T: %1.0d F:%1.0d Long: %1.0f ",X.id_talle, X.id_rang, X.Taillefeuille);
static char surf[80];
sprintf(surf,"Surface %6.2f", surfacefeuille); 
if (X.id_talle>=0)
{
 nproduce /*Label(rang) *//*f(10)Label(surf)*/;   /* Module affichage sortie ecran maquette de i) indice de chaque feuille ii) Surface feuille */
}
produce ;              
  }
       }


Gaine(longueur, angletal, rang, talle):
{ 
 double rank=rang;
 double large=(2/3.1415)*(1.9645*log(double(longueur*5))-5.701);
 
    int numstepsg=ceil(longueur/LSTEP);
    for (int i=0; i<numstepsg; i++)
    {
         if(i<(10/LSTEP))
    	 {
           if(i==1)
            {
             nproduce SetWidth(large) CurrentContour(Gaino)StartGC()SetColor(3)Down(Torsiongaine*pow(-1,rank))/*Left(angletal)*/F(LSTEP)EndGC();
            }
     	 nproduce SetWidth(large) CurrentContour(Gaino)StartGC()SetColor(3)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
     	 }
         else 
         {
         nproduce SetWidth(large) CurrentContour(Gaino)StartGC()SetColor(3)F(LSTEP)EndGC();    
         }
 } 
produce ;
}




Centralerect(age, rang,Rot, angle,geno):
{  
  
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/ Down(angle) SetWidth(2)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/ RollR(Rot-90)SetColor(14)SetWidth(1)F(0.1);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.1);
    }
  }
produce;
}




Affichage(para):
{ 











  static char Biomasser[80];
  sprintf(Biomasser,"Biomasse Recoltee = %7.0f", Biomasse);
  static char surface[80];
  sprintf(surface,"Surface foliaire tot. = %7.0f", Surfoliairetotale);
  static char Tallage1[80];
  sprintf(Tallage1,"NB Talles = %7.0i", NBtalles[1]);
  static char Tallage2[80];
  sprintf(Tallage2,"NB Talles = %7.0i", NBtalles[2]);


 if (Surfoliairetotale>STOPTALLAGE)
  {
    nproduce SetWidth(20) SetColor(17);
  } 
  Surfoliairetotale=0;
  nproduce MoveTo(-GROUNDSIZE, 7, GROUNDSIZE/2)F(0)/*Label(Biomasser)*/;
  nproduce MoveTo(-GROUNDSIZE, 10, GROUNDSIZE/2)F(0)/*Label(surface)*/;   /*surf foliaire totale */
  nproduce MoveTo(-GROUNDSIZE*0.5, 10, GROUNDSIZE/2)F(0)Label(Tallage1);
  nproduce MoveTo(GROUNDSIZE*0.5+30, 10, GROUNDSIZE)F(0)/*Label(Tallage2)*/;



  nproduce MoveTo(-GROUNDSIZE, 80, GROUNDSIZE*0.25)Right(90)SetColor(Lumiere)SetWidth(7) /*F(GROUNDSIZE)*/;
  nproduce MoveTo(-GROUNDSIZE, 80, GROUNDSIZE*0.50)SetColor(Lumiere)SetWidth(7) /*F(GROUNDSIZE)*/;
  nproduce MoveTo(-GROUNDSIZE, 80, GROUNDSIZE*0.75)SetColor(Lumiere)SetWidth(7) /*F(GROUNDSIZE)*/;
  
  produce;
}

/* module d'affichage pour les 3 neons et pour affichage de surface etc...*/


/*---------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   






/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------    GROUPE DEFOLIATION  -----------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group Defoliation :



/*-----------------------------------------DEFOLIATION PRODUCTION --------------------------------------------------------*/

production:

SB() GetPos(X, Y, Z) EB(): 
{ 
 if (/*Y>HCOUP && (*/TPS==DCOUP+(2*dt)/*)*/) 
  produce;
}

Centralerect(age, rang,Rot, angle,geno):
{
produce Centralerect(age, rang,Rot, angle,geno);
}


Gaine(longueur, angletal, rang, talle)>Feuille(X):
{
  if (X.id_plante>=0)
  {
  produce ; 
  }                       /*Pour se debarasser des vieilles gaines */
}

B(D):
{ 
 
if (TPS==DCOUP+dt*5) 
    { 
    if (D.coupe>1)
      {
       D.Cutstatus=1;
       produce B(D);
      }
    }
if (TPS==DCOUP+dt*6) 
    {
     if (D.Taillegaine>D.longueur && D.Taillefeuille>(D.Ymax*0.9)) /* La condition *0.9 c'est pour les cas ou une gaine en croissance est sectionnée, il faut alors qu'elle finisse sa croissance, donc elle a besoin de la feuille qu'elle porte */
      {
       D.Cutstatus=2;
       D.Taillegaine=D.longueur;
       Seuil[D.id_talle][D.id_plante]=D.Taillegaine; 
       Tailfingain[D.id_talle][D.id_rang][D.id_plante]=D.Taillegaine; 
                                                          /* Pour avoir feuilles suivantes moins longues !! */
      }
      else
      {
       /*D.Taillefeuille=D.longueur;*/
       D.Taillegaine=D.Taillegaine;
       D.Taillelimbe=D.Taillelimbe-D.coupe;
       
      }
      produce B(D);
    }
}


/*-------------------------------------------DEFOLIATION  DECOMPOSITION -----------------------------------------------------*/

decomposition:

Feuille(X):
{
 if (X.id_plante==1 || X.id_plante==2 || X.id_plante==4)
  {
  int numsteps=ceil((X.Taillefeuille-X.coupe)/LSTEP);
  for (int i=0; i<numsteps; i++)
  {
   Feuildefoli Z={X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle, X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, ((X.Taillefeuille-X.coupe)/LSTEP)-i,X.rapportK, 0, X.angleinsert, X.coupe, X.geno, X.Assimibrut, X.Besoinencroiss};
   nproduce  C(Z) SB() GetPos(0, 0, 0) EB();
  }   
  Feuildefoli Z2={X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle,X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, 0, X.rapportK, 0, X.angleinsert, X.coupe, X.geno, X.Assimibrut, X.Besoinencroiss};   
  FeuilTrans D={X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle,X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, X.Taillefeuille,0,X.rapportK,0, X.angleinsert, X.coupe, X.geno, X.Assimibrut, X.Besoinencroiss};
  nproduce C(Z2) B(D);
  produce ;
 }
}




/*------------------------------------- DEFOLIATION  INTERPRETATION -----------------------------------------------------*/


interpretation:

C(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{
   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil((Z.Taillegaine)/LSTEP);
    
 



  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;
  
  
 
    taille=Z.Taillelimbe;
    L0=Z.Taillefinalel-taille;   
   


   /*----------------- Calcul largeur ---------------------------------------*/


  if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
    else
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }

  /*---------------------------------------------------------------------------------------*/    

   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-(10/LSTEP))
      {
       double rank=Z.id_rang;
       produce  StartGC() CurrentContour(Limbo) SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
      }
      else
      {
       produce StartGC()CurrentContour(Limbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
    
       if (Z.indice>(numsteps-numstepsg))  /* GAINE - Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-(10/LSTEP))
           {
           produce SetWidth((2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701))SetColor(2)CurrentContour(Gaino)StartGC()Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
           }
           else
           { 
	   produce SetWidth((2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701))SetColor(2)CurrentContour(Gaino)StartGC()F(LSTEP)EndGC();
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       {
        produce SetWidth((2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701))Left(Z.angleinsert)SetColor(1)F(LSTEP);   /* LIGULE */
       }

       if (Z.indice<(numsteps-numstepsg))     /* LIMBE */
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-(10/LSTEP))
                  {
                   produce  StartGC() CurrentContour(FoldLimbo) SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                  }
                  else
                  {
	          produce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               produce StartGC() CurrentContour(Limbo) Left(curve)SetWidth(large)SetColor(1)F(LSTEP)EndGC();  
              }                 
          
         }
       }

   }

}
         

B(D):
{
 produce SetColor(17) Sphere(1);
}

Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/Down(angle) SetWidth(3)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/RollR(Rot-90)SetColor(14)SetWidth(2)F(0.3);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.3);
    }
  }
  produce;
}

Marqueur(onsenfout):
{
 produce SetColor(7) Sphere(1);
}


Affichage(para):
{  
   produce MoveTo(0/*+GROUNDSIZE*/, HCOUP, 0) SetColor(17)SetWidth(0.3)GROUND;
}
/*---------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   






/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------    GROUPE PHOTOSYNTH -----------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group Photosynth :



/*------------------------------------------------PHOTO PRODUCTION --------------------------------------------------------*/

production:

Soleil(Elev, Azim):
{


  /* --------------Calcul longueur et largeur de la cage------------------- */

      if(TPS==Photosynthese)  
      {
       Bornescage[1]=-100;
       Bornescage[2]=100;
       Bornescage[3]=-100;
       Bornescage[4]=100;
      }


      Longstrate=Bornescage[1]-Bornescage[2];
      Largstrate=Bornescage[3]-Bornescage[4];


   /*-------------------------Interception réalisée par strate -------------------------- */

   if (TPS==Photosynthese+3*dt)
   {
    
    
    float Io=DensiteNRJincident*abs(Largstrate*Longstrate)/**sin(Elevation[0]*(3.14/180))*/;     /* Reactive le sinus si densité NRJ pas par unité de surface de SOl mais par unite de surf face à rayonnement *//* Densité NRJetique incident */
    
    for (int i=0; i<nbstrates; i++){Transmis[i]=0;};
    for (int i=0; i<nbstrates; i++){Intercepte[i]=0;};

     /* Calcul de l'interception */


     for (int i=nbstrates-1; i>=0; i--)  /* Extinction et Interception par Strates */
     {

      float Proportiontrous=1-(SURFPROJPARSTRATE[i]/abs(Largstrate*Longstrate));                                                       /*hstrate*(SURFTOTPARSTRATE[i]/(Largstrate*Longstrate*hstrate));*/ /* Chemin parcouru dans layer de Densité Optique */
                                                                          /* Devrait ajuster chemin en fonction de direction soleil, à faire */ /*Revient au meme de faire *LAI que de faire *LAD cumulé * longueur vu que m2/m2 est adimentionnel comme m3/m3 */
       if (i==(nbstrates-1))
       {
        Transmis[i]=max(0,Io*Proportiontrous);         /* avant : Io*exp(-coeffKmoy[i]*LAD)  */
        Intercepte[i]=Io-Transmis[i];
       }
      else
       {
        Transmis[i]=max(0,Transmis[i+1]*Proportiontrous);
        Intercepte[i]=Transmis[i+1]-Transmis[i];
       }Printf("Transmis %f\n", Transmis[i]);
        
     }
   }

 produce Soleil(Elevation[0], Azimuthsoleil);
}







Entrenoeud(X):                    /* Flushing du carbone du photosynth precedent */
{
X.Assimdescend=0; 
produce Entrenoeud(X);
}



GetPos(X, Y, Z): 
{ 
 if (/*Y>HCOUP && (*/TPS==DCOUP+(2*dt)/*)*/) 
  produce;
}

Centralerect(age, rang,Rot, angle,geno):
{
produce Centralerect(age, rang,Rot, angle,geno);
}


Gaine(longueur, angletal, rang, talle)>Feuille(X):
{
  if (X.id_plante>=0)
  {
  produce ; 
  }                      
}
/*Pour se debarasser des vieilles gaines */




Cp(Z1) << Cp(Z) > SB() GetPos(X, Y, Za)EB() SB()GetLeft(Ac, Bc, Cc) EB()SB()GetHead(Ah, Bh, Ch) EB():
{
  

  /* Calcul proprietés directionnelles */
  float Azimuth;
  float Inclinaison;
  
   

  if ((Z.indice+(Z.dejacoupe/LSTEP))>=ceil((Z.Taillefeuille-Z.Taillegaine)/LSTEP)) /* Pour GAINES */
  {   
     /* Doit determiner vecteur normal tel que coplanaire directeur et solaire */
    
     /* Vecteur directeur soleil */ 
     float  Xsoleil=cos(Elevation[0]*(3.14/180))*sin(Azimuthsoleil*(3.14/180));
     float  Ysoleil=-sin(Elevation[0]*(3.14/180));
     float  Zsoleil=cos(Elevation[0]*(3.14/180))*-cos(Azimuthsoleil*(3.14/180));
     
     /* Produit vectoriel rayon x tube */   
      float Xvectoriel=Ysoleil*Ch-Zsoleil*Bh;
      float Yvectoriel=Zsoleil*Ah-Xsoleil*Ch;
      float Zvectoriel=Xsoleil*Bh-Ysoleil*Ah;

      float Rapp=((-Ch*Xvectoriel/Ah)+Zvectoriel)/((-Bh*Xvectoriel/Ah)+Yvectoriel);
  
      /* Calcul coordonnées vecteur normal surface - signe sans importance */  
      float Zinconnu=sqrt(1/(1+pow(Rapp, 2)+pow((Bh*Rapp-Ch)/Ah, 2)));
      float Yinconnu=-Zinconnu*Rapp;
      float Xinconnu=Zinconnu*((-Ch/Ah)+((Bh/Ah)*Rapp));
  
      /* Inclinaison de la normale à la surface par rapport à verticale : 0°= surface horizontale */
      Inclinaison=min(acos(Yinconnu/sqrt(Xinconnu*Xinconnu+Yinconnu*Yinconnu+Zinconnu*Zinconnu))*(180/3.14),180-acos(Yinconnu/sqrt(Xinconnu*Xinconnu+Yinconnu*Yinconnu+Zinconnu*Zinconnu))*(180/3.14)); /* Inclinaison de la normale à la surface par rapport à verticale : 0°= surface horizontale */
    
      /* Calcul azimuth : angle par rapport à Z. 0° = vers nous. 90°=a gauche. A 180° pres */

      if(Xinconnu<0 && Zinconnu>=0){Azimuth=atan(abs(Xinconnu/Zinconnu))*(180/3.14);}
      if(Xinconnu<0 && Zinconnu<0){Azimuth=180-atan(abs(Xinconnu/Zinconnu))*(180/3.14);}
      if(Xinconnu>=0 && Zinconnu<0){Azimuth=180+atan(abs(Xinconnu/Zinconnu))*(180/3.14);}
      if(Xinconnu>=0 && Zinconnu>=0){Azimuth=360-atan(abs(Xinconnu/Zinconnu))*(180/3.14);}

     /*Printf("XRay: %f\n YRay:%f\n ZRay:%f\n", Xsoleil, Ysoleil, Zsoleil);
     Printf("XHead: %f\n YHead:%f\n ZHead:%f\n", Ah, Bh, Ch);
     Printf("Xvect: %f\n Yvect:%f\n Zvect:%f\n", Xvectoriel, Yvectoriel, Zvectoriel);
     Printf("Xincon: %f\n Yincon:%f\n Zincon:%f\n", Xinconnu, Yinconnu, Zinconnu);
     Printf("Inclin :%f\n", Inclinaison);
     Printf("Azim : %f\n ", Azimuth);*/
  
  }
  else  /* Pour LIMBES */
  {

    /* Inclinaison de la normale à la surface par rapport à verticale : 0°= surface horizontale */
  
    Inclinaison=min(acos(Bc/sqrt(Ac*Ac+Bc*Bc+Cc*Cc))*(180/3.14),180-acos(Bc/sqrt(Ac*Ac+Bc*Bc+Cc*Cc))*(180/3.14)); 
  
    float Az=-Ac; float Bz=-Bc; float Cz=-Cc; /* Pour que vecteur normal a feuille regarde vers le haut au lieu bas */

    if(Az<0 && Cz>=0){Azimuth=atan(abs(Az/Cz))*(180/3.14);}
    if(Az<0 && Cz<0){Azimuth=180-atan(abs(Az/Cz))*(180/3.14);}
    if(Az>=0 && Cz<0){Azimuth=180+atan(abs(Az/Cz))*(180/3.14);}
    if(Az>=0 && Cz>=0){Azimuth=360-atan(abs(Az/Cz))*(180/3.14);}

   /* Printf("Inclin :%f\n", Inclinaison);
      Printf("Azim : %f\n ", Azimuth);
      Printf("Ac: %f\n Bc: %f\n Cc: %f\n", Ac, Bc, Cc);*/
   }
 

   float Type;             /* Pour dire si gaine ou limbe le C(Z) : 1=gaine, 2=Limbe */
   float large=0;            /* Largeur de chaque element de surface */
   double L0;
   double taille;

   taille=Z.Taillelimbe;
   L0=Z.Taillefinalel-taille;   
    
   if (Z.indice<(ceil((Z.Taillefeuille-Z.Taillegaine)/LSTEP)))
   { 
    Type=2;  /* Limbe*/
   }
   else
   {
    Type=1;
   }


   /*----------------- Calcul largeur pour surface ---------------------------------------*/



   if( Z.indice < (Z.Taillelimbe/LSTEP)+1)
   {
     if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
      {
       large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    // large=0;
       }
      else
      {
       double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
       large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
     
       //large=0;
      }
    }

     if( Z.indice >= (Z.Taillelimbe/LSTEP)+1)
   {
      if (Z.Taillegaine>max(Tailfingain[Z.id_talle][Z.id_rang-1][1]+1,Tailfingain[Z.id_talle][Z.id_rang-2][1]+1) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang+1][1] )
       { 
        large=(2*(((2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701))/4)*3.1415)/2; /* divisé par 4 parce si met setWidth 2 alors fait un tube de 1 de DIAMETRE, donc j'ai le rayon et je fait 2PiR, divisé par deux ensuite car multiplié par 2 dans calcul surf. */
      //large=0;
       
      }
     else
      {    
      large=0;
      }
    }
    
    /* ATTENTION là calcul large de derniere gaine de brin et c'est tout */


    /* Pour pas introduire les morceau de feuilles/gaine dans pseudotiges dans calcul : si plus haut que gaine d'avant alors large=0 */
     
 

    /*------------------- Bornes de la cage d'analyse adaptées à plante ---------------------- */

  /* if (TPS==Photosynthese+dt)
   {  
       if (X>Bornescage[1]){Bornescage[1]=X;}
       if (X<Bornescage[2]){Bornescage[2]=X;}       
       if (Za>Bornescage[3]){Bornescage[3]=Za;}
       if (Za<Bornescage[4]){Bornescage[4]=Za;}  
   }*/  

   /*-------------- Incrementation classe inclinaison et surfaces par strates ---------- */

   if (TPS==Photosynthese+2*dt)
   {
    
     if(ceil(Z1.indice)==2)
       {
       /*Printf("Mon X:%f Mon Y:%f Mon Z:%f\n", X, Y, Za);
         Printf("Larg : %f Long : %f\n", Largstrate, Longstrate);
         Printf("X+: %f X-: %f Z+: %f Z-: %f\n", Bornescage[1],Bornescage[2],Bornescage[3],Bornescage[4]);*/
       }


     for(int i=0; i<9; i++) 
       {
         if ( Inclinaison>(i*10) && Inclinaison<(i+1)*10 )
           {     
             CLASSINCLINAISON[i][nbstrates]=CLASSINCLINAISON[i][nbstrates]+(LSTEP*large*2);  /* Classe inclinaison globale */
           }
       }

      for(int i=0; i<nbstrates; i++) 
       {
         if(Y>i*hstrate && Y<=(i+1)*hstrate || i==(nbstrates-1) && Y>i*hstrate)
         { 
           for(int k=0; k<9; k++) 
             {
                   if ( Inclinaison>(k*10) && Inclinaison<(k+1)*10 )                     
                    {    
                     CLASSINCLINAISON[k][i]=CLASSINCLINAISON[k][i]+(LSTEP*large*2);                /* Classe inclinaison par strate */
                    }
             }
           for(int j=0; j<nbstrates; j++)
            {
                   if (j==i)
                    {
                      Z.SStratPropr[j]=Z1.SStratPropr[j]+(LSTEP*large*2);
                      SURFTOTPARSTRATE[j]=SURFTOTPARSTRATE[j]+(LSTEP*large*2);    /* Fait remonter info vers bout de feuille */
                    }
                   else
                    {
                      Z.SStratPropr[j]=Z1.SStratPropr[j];                         /* Fait remonter info vers bout de feuille */
                    }
              
             }
          }
        } 

        /*------------------------------- Calcul Coeff de Projection ------------------------------ */

         float coeffKpropr=0;       
                                           /* Calcul du coefficient de projection du morceau pour chaque classe d'elevation solaire */
         for (int j=0; j<nbdirsun; j++)     /* Elevation soleil */
          { 
             float Gfunc=0;
             Gfunc=(abs((cos(Inclinaison*3.14/180)*sin(Elevation[j]*3.14/180))+(sin(Inclinaison*3.14/180)*cos(Elevation[j]*3.14/180)*cos((Azimuth-Azimuthsoleil)*3.14/180))))/sin(Elevation[j]*3.14/180);
             coeffKpropr=coeffKpropr+Gfunc;
          }
       
         float coeffprojmoy=coeffKpropr/nbdirsun;   /* Coefficient d'Extinction moyen sur les differentes direction solaires.*/ 
         //Printf("CoeffK %f\n", coeffKpropr);
         

         for(int i=0; i<nbstrates; i++) 
         {
            if(Y>i*hstrate && Y<=(i+1)*hstrate || i==(nbstrates-1) && Y>i*hstrate)
             { 
               
              if( Z.indice < (Z.Taillelimbe/LSTEP)+1)
               {
                SURFPROJPARSTRATE[i]=SURFPROJPARSTRATE[i]+((LSTEP*large*2)*coeffprojmoy);
                Z.SprojPropr[i]=Z1.SprojPropr[i]+((LSTEP*large*2)*coeffprojmoy);
               }
              if( Z.indice >= (Z.Taillelimbe/LSTEP)+1)  /* Prend que moitié surface en projection si gaine !! */
               {
                SURFPROJPARSTRATE[i]=SURFPROJPARSTRATE[i]+((LSTEP*large)*coeffprojmoy);
                Z.SprojPropr[i]=Z1.SprojPropr[i]+((LSTEP*large)*coeffprojmoy);  /* Calcul surf projetée totale de la feuille en question dans strate i */
               }
             }
            else
             { 
               Z.SprojPropr[i]=Z1.SprojPropr[i];
             } 
         }
      
              
     }

   produce Cp(Z);
}





Cp(Z) << B(D):
{
  float Assimimax=40;   /* taux Assimilation à saturation lumineuse */
  float epsilon=1;    /* light use efficiency gCO2 cm-2 h-1/ J m-2 h-1 */


 if (TPS==Photosynthese+dt*4)
  {
   for(int i=0; i<nbstrates; i++)
   {
    float PARintercepte=0;
    

    if(SURFTOTPARSTRATE[i]>0)
      {
       PARintercepte=Intercepte[i]*(Z.SprojPropr[i]/SURFPROJPARSTRATE[i]);    /* Repartition lum pro rata de surf projetée propre */ 
       D.Assimibrut=D.Assimibrut+(Assimimax*(1-exp(-epsilon*PARintercepte/Assimimax))); /* taux Assimilation brut CO2 : gCO2 cm-2 h-1 */
      }
     else
      {
       D.Assimibrut=D.Assimibrut;   /* Pour shunter les cas ou strates vide ou on a donc division par 0 du SURFTOTSTRATE */
      } 
    
    
     /* 
     Printf("Somme projection : %f\n", Z.SprojPropr[i]);   
     
     Printf("SURFPROJ : %f\n", SURFPROJPARSTRATE[i]);
     Printf("INTER : %f\n", Intercepte[i]);
     Printf("PARi : %f\n", PARintercepte);
      */    
   }
    
    
  }
 
  produce B(D);
}






/*------------------------------------------------ PHOTO DECOMPOSITION -----------------------------------------------------*/

decomposition:

Feuille(X):
{ 
  int numsteps=ceil((X.Taillefeuille-X.coupe)/LSTEP);
  for (int i=0; i<numsteps; i++)
  {
   Feuilphoto Z={X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle, X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, ((X.Taillefeuille-X.coupe)/LSTEP)-i,X.rapportK, 0, X.angleinsert, X.coupe, X.geno, /*X.Assimibrut*/0, X.Besoinencroiss};
   for (int i=0; i<nbstrates; i++){Z.SStratPropr[i]=0;Z.SprojPropr[i]=0;}
   nproduce   Cp(Z) SB() GetPos(0, 0, 0) EB()SB()GetLeft(0, 0, 0) EB()SB()GetHead(0, 0, 0) EB();  /* GetLeft==GetUP >>> Pour pas avoir a faire rotation 90 !!!!!!!!!!!*/
  }   
  FeuilTrans D={X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle,X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, X.Taillefeuille,0,X.rapportK,0, X.angleinsert, X.coupe, X.geno, /*X.Assimibrut*/0, X.Besoinencroiss};
  nproduce B(D);
  produce ;
}




/*------------------------------------------------PHOTO  INTERPRETATION -----------------------------------------------------*/


interpretation:


   Cp(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{
   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);
        

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*----------------- Calcul largeur ---------------------------------------------*/


   
    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
     else if (Z.indice+(Z.dejacoupe/LSTEP)<(numsteps-numstepsg))  
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }
     else if (Z.indice+(Z.dejacoupe/LSTEP)>(numsteps-numstepsg)) 
    {
     large=(2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701);
    }




   /* ----------------------------------------------------------------------------*/


   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-(10/LSTEP))
      {
       double rank=Z.id_rang;
       produce  StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
      }
      else
      {
       produce StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
       
     /*if(ceil(Z.indice)==numsteps) { produce RollR(0)Down(90); }*/


       if (Z.indice>(numsteps-numstepsg))  /* Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-(10/LSTEP))
           {
           produce SetWidth(large)SetColor(2)Down(Torsiongaine*pow(-1,rank))F(LSTEP);
           }
           else
           { 
	   produce SetWidth(large)SetColor(2)F(LSTEP);
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       {
        produce SetWidth(large)Left(Z.angleinsert)SetColor(1)F(LSTEP); 
       }

       if (Z.indice<(numsteps-numstepsg)) 
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-(10/LSTEP))
                  {
                   nproduce   StartGC()CurrentContour(FoldLimbo) SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                   
                   produce;
                  }
                  else
                  {
	            nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
                   
                    produce;
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               nproduce  StartGC() CurrentContour(Limbo3)Left(curve)SetWidth(large)SetColor(1)F(LSTEP)EndGC();  
               produce;
              }                 
          
         }
       }

   }

}


         

B(D):
{ 
  static char Assimil[30];
  sprintf(Assimil,"Assimilation : %4.2f", D.Assimibrut); 
  if(D.Assimibrut>2)
  {
   nproduce SB()SetWidth(0)F(8)SetColor(D.Assimibrut)Label(Assimil)EB();
  }
  produce SetColor(D.Assimibrut)Sphere(2);
}

Soleil(Elev, Azim):    /*Représentation du soleil dans la voute et de Azimuth=0° */
{
  if(TPS==Photosynthese+2*dt || TPS==Photosynthese+3*dt || TPS==Photosynthese+4*dt)
  {
  nproduce MoveTo(-GROUNDSIZE/2,0,GROUNDSIZE/2) ;
  nproduce SB()RollToVert Down(90) SetColor(26) F(50) EB();  /* Indique Azimuth=0° */
  nproduce SB()RollToVert() RollR(Azim) Down(90-Elev) SetColor(25) SetWidth(0.5) F(200)  Sphere(10)EB();
  }
  produce;
}



Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/Down(angle) SetWidth(3)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/RollR(Rot-90)SetColor(14)SetWidth(2)F(0.3);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.3);
    }
  }
 
  produce;
}

Marqueur(onsenfout):
{
 produce SetColor(7) Sphere(1);
}

Affichage(para):
{ 

for (int i=0; i<nbstrates; i++)
 {
  nproduce MoveTo(Bornescage[1], i*hstrate,Bornescage[4]) SB() SetColor(33)SetWidth(1) F(0) CAGE EB();
 }
produce;  
}


/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------GROUP FOLD----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group FOLD :



production:

Cp(Z) > B(D):
{
    Z.Assimibrut=D.Assimibrut; 
    Z.Besoinencroiss=D.Besoinencroiss; 
    produce Cp(Z);  
}



Cp(Z1) >> Cp(Z2):      // Descente de la valeur du C assimilé le long de la feuille
{
   Z1.Assimibrut=Z2.Assimibrut;
   Z1.Besoinencroiss=Z2.Besoinencroiss;
   produce Cp(Z1);
}


Entrenoeud(X) >> SB()Cp(Z)EB() Apex(rang, talle,plante,geno):   // Passage dans la "tige"
{
 X.Assimdescend=Z.Assimibrut;
 X.Demandescend=Z.Besoinencroiss;
 produce Entrenoeud(X);
}

Entrenoeud(X1) >> SB()Cp(Z)EB() SB()ApexTal(par,geno)EB() Entrenoeud(X2)  :  
{
 X1.Assimdescend=Z.Assimibrut+X2.Assimdescend;
 X1.Demandescend=Z.Besoinencroiss+X2.Demandescend;
 produce Entrenoeud(X1);
}


Entrenoeud(X1) >> SB()Cp(Z)EB()SB()Right(angle)Centralerect(para,para2, para3,para4,para5)GetHead(Aa,Bb, Cc)Entrenoeud(X2)EB()Entrenoeud(X3):  
{
 X1.Assimdescend=Z.Assimibrut+X2.Assimdescend+X3.Assimdescend;
 X1.Demandescend=Z.Besoinencroiss+X2.Demandescend+X3.Demandescend;
 produce Entrenoeud(X1);
}

Pool(ParaCpool, ParaCdemand) >> Centralerect(para,para2, para3,para4,para5)GetHead(Aa,Bb, Cc)Entrenoeud(X1):
{
 produce Pool(X1.Assimdescend,X1.Demandescend);
}

/*------------------------------------ FOLD INTERPRETATION ---------------------------------------------------------------------------------------*/

interpretation:

 Cp(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{
   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);
        

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*----------------- Calcul largeur ---------------------------------------*/


    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
    else
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }
    


/* ----------------------------------------------------------------------------*/


   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-(10/LSTEP))
      {
       double rank=Z.id_rang;
       produce  StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
      }
      else
      {
       produce StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
    
       if (Z.indice>(numsteps-numstepsg))  /* Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-(10/LSTEP))
           {
           produce SetWidth(large)SetColor(2)Down(Torsiongaine*pow(-1,rank))F(LSTEP);
           }
           else
           { 
	   produce SetWidth(large)SetColor(2)F(LSTEP);
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       {
        produce SetWidth(large)Left(Z.angleinsert)SetColor(1)F(LSTEP); 
       }

       if (Z.indice<(numsteps-numstepsg)) 
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-(10/LSTEP))
                  {
                   nproduce   StartGC()CurrentContour(FoldLimbo) SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                   
                   produce;
                  }
                  else
                  {
	            nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
                   
                    produce;
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               nproduce  StartGC() CurrentContour(Limbo3)Left(curve)SetWidth(large)SetColor(1)F(LSTEP)EndGC();  
               produce;
              }                 
          
         }
       }

   }

}


         

B(D):
{
  produce SetColor(D.Assimibrut)Sphere(2);
}

Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/Down(angle) SetWidth(3)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/RollR(Rot-90)SetColor(14)SetWidth(2)F(0.3);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.3);
    }
  }
 
  produce;
}

Pool(Ccontenu,demande):
{
       static char Assim[30];
       sprintf(Assim,"Ass TOTALE: %6.2f Demand TOTALE: %6.2f", Ccontenu,demande); 
       
       nproduce SB()Down(180)SetColor(Ccontenu/5)SetWidth(4)F(30)Sphere(/*Ccontenu*/10)SetWidth(0)F(15)Label(Assim)EB();    
       produce;
}




Marqueur(onsenfout):
{
 produce SetColor(7) Sphere(1);
}

Affichage(para):
{   
 static char affich[30];
 sprintf(affich,"FOLD"); 
 nproduce MoveTo(-120,40,0) SetColor(255)Label(affich);
 nproduce MoveTo(-100,50,0) SetWidth(3)Down(180)F(10)SetWidth(6)StartGC()SetWidth(0.1)F(10)EndGC();
 produce;
}



/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------GROUP DECOMPOECRITURE----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group DECOMPOECRIT:


decomposition:

Feuille(X):
{ 
  int numsteps=ceil((X.Taillefeuille-X.coupe)/LSTEP);
  for (int i=0; i<numsteps; i++)
  { 
   nproduce   Cdecomp(X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle, X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, ((X.Taillefeuille-X.coupe)/LSTEP)-i,X.rapportK, 0, X.angleinsert, X.coupe, X.geno, /*X.Assimibrut*/0, X.Besoinencroiss) SB() GetPos(0, 0, 0) EB()SB()GetUp(0, 0, 0) EB()SB()GetHead(0, 0, 0) EB();  /* GetLeft==GetUP >>> Pour pas avoir a faire rotation 90 !!!!!!!!!!!*/
  }   
  nproduce Bdecomp(X.age, X.Agecroiss, X.Taillefeuille,X.id_rang, X.id_talle,X.id_plante, X.Ymax, X.Difftps, X.Reservoir, X.Taillefinalel, X.Taillefinaleg, X.Taillelimbe, X.Taillegaine, X.Phase, X.Taillefeuille,0,X.rapportK,0, X.angleinsert, X.coupe, X.geno, /*X.Assimibrut*/0, X.Besoinencroiss);
  produce ;
}


/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------GROUP ECRITURE----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group ECRITURE:


production:


Cp(Z)>SB() GetPos(Xpos, Ypos, Zpos) EB()SB()GetUp(Xleft, Yleft, Zleft) EB()SB()GetHead(Xhead, Yhead, Zhead) EB():
 {



/* Gain/Limb \t Hauteur \t Diamgain \t id-plant \t id-talle \t id-feuil \t XPos \t YPos \t ZPos \t XHead \t YHead \t ZHead \t XLeft \t YLeft \t ZLeft \t\n" */

  int numsteps=ceil((Z.Taillefeuille)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*---------------------------------------------- Calcul largeur ---------------------------------------*/


    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
    else
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }

  /*------------------------------------------ Calcul largeur ---------------------------------------*/



float idGL=0;    /* Identifiant gaine ou limbe : 1=gain, 2=limb */


if (Z.indice+(Z.dejacoupe/LSTEP)<=(numsteps-numstepsg)) {idGL=2;}
else if (Z.indice+(Z.dejacoupe/LSTEP)>(numsteps-numstepsg)) {idGL=1;}
float diamgain=0;
float hauteurg=0;

if(idGL==1) {
hauteurg=/*Z.Taillegaine*/LSTEP;
diamgain=(2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701);}

if(idGL==2) {
hauteurg=LSTEP;
diamgain=large*2;
}  /* Longueur de chaque morceau entre les points */


if(TPS==tpsecriture+2*dt){

fprintf(fpsss, "%f\t %f\t %f\t %i\t %i\t %i\t %2.0f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", idGL,hauteurg,diamgain,Z.id_plante,Z.id_talle,Z.id_rang,Z.indice,Xpos, Ypos, Zpos, Xhead, Yhead,Zhead, Xleft, Yleft, Zleft);
//fprintf(fpsss, "%f\t %f\t %f\t %i\t %i\t %i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", idGL,hauteurg,diamgain,Z.id_plante,Z.id_talle,Z.id_rang,Xpos, Ypos, Zpos, Xhead, Yhead,Zhead, Xleft, Yleft, Zleft);

//fprintf(fpsss, "%f\t %f\t %f\t %i\t %i\t %i\t %i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n", idGL,hauteurg,diamgain,/*Z.id_plante,/*Z.id_talle+50*/444,/*Z.id_rang+10*/33,/*Z.indice+1300*/555,Xpos, Ypos, Zpos, Xhead, Yhead,Zhead, Xleft, Yleft, Zleft);


}
 }

 



SB():
 {if(TPS==tpsecriture+1*dt)
 {fprintf(fpss, "SB()");}
 }

EB():
 {if(TPS==tpsecriture+1*dt){
 fprintf(fpss, "EB()");}
 }
 
RollToVert():
 {if(TPS==tpsecriture+1*dt){
 fprintf(fpss, "RollToVert");}
 }

Left(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Left(%1.0f)", para);}
 }

Right(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Right(%1.0f)", para);}
 }

RollR(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "RollR(%1.0f)", para);}
 }

RollL(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "RollL(%1.0f)", para);}
 }

Up(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Up(%1.0f)", para);}
 }

Down(para):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Up(%1.0f)", para);}
 }
/*
GetPos(X, Y, Z):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "GetPos(%f,%f,%f)", X,Y,Z);}
 }

GetLeft(X, Y, Z):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "GetLeft(%f,%f,%f)", X,Y,Z);}
 }

GetHead(X, Y, Z):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "GetHead(%f,%f,%f)", X,Y,Z);}
 }

*/
Centralerect(age, rang,Rot, angle,geno):
{if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Centralerect(%3.0f,%2.0f,%2.0f,%2.0f,%2.0f)", age,rang,Rot,angle,geno);}
}

Cdecomp(a,b,c,d,e,foo,goo,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v):
 {
if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"Cdecomp(%1.0f,%1.0f,%1.0f,%i,%i,%i, %1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%2.3f)\n",a,b,c,d,e,foo,goo,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v);
  Feuilphoto Z={a,b,c,d,e,foo,goo,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v};
  produce Cp(Z);}
 }


Apex(id_rang, id_talle, id_plante, id_geno) :
{ if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Apex(%i,%i,%i,%i)\n",id_rang,id_talle,id_plante,id_geno);}
}


ApexTal(id_rang,id_geno):
{if(TPS==tpsecriture+1*dt){
 fprintf(fpss, "ApexTal(%i,%i)\n",id_rang,id_geno);}
}


Entrenoeud(X):
 {if(TPS==tpsecriture+1*dt){
  fprintf(fpss, "Entrendecomp(%1.0f,%i,%i,%i,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f)\n",X.val, X.id_rang, X.id_talle, X.id_plante, X.angletal, X.Assimdescend, X.Demandescend, X.Assimmonte, X.Demandmonte);}
 }

Bdecomp(age,Agecroiss,Taillefeuille,id_rang,id_talle,id_plante,Ymax,Difftps,Reservoir,Taillefinalel,Taillefinaleg,Taillelimbe,Taillegaine,Phase,longueur,cutstatus,rapportK,coupe, angleinsert, dejacoupe, geno, Assimibrut, Besoinencroiss):
{ if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"Bdecomp(%1.0f,%1.0f,%1.0f,%i,%i,%i,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.3f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%2.3f)\n",age,Agecroiss,Taillefeuille,id_rang,id_talle,id_plante,Ymax,Difftps,Reservoir,Taillefinalel, Taillefinaleg, Taillelimbe, Taillegaine, Phase,Taillefeuille,cutstatus,rapportK,coupe, angleinsert, dejacoupe, geno, Assimibrut, Besoinencroiss);
  FeuilTrans D={age,Agecroiss,Taillefeuille,id_rang,id_talle,id_plante,Ymax,Difftps,Reservoir,Taillefinalel, Taillefinaleg, Taillelimbe, Taillegaine, Phase,Taillefeuille,0,rapportK,0, angleinsert, coupe, geno, 0, Besoinencroiss};
  produce B(D);}
}


MoveTo(Xm, Ym, Zm):
{
 if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"MoveTo(%f,%f,%f)", Xm,Ym,Zm);
}
}

Pool(ccontenu, demande):
{
 if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"Pool(%f,%f)", ccontenu, demande);
}
}

Racine(ccontenu):
{
 if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"Racine(%f)", ccontenu);
}
}

Reserve(ccontenu):
{
 if(TPS==tpsecriture+1*dt){
  fprintf(fpss,"Reserve(%f)", ccontenu);
}
}







/*

Alban
*/






/*------------------------------------ ECRITURE INTERPRETATION ---------------------------------------------------------------------------------------*/


interpretation:


   Cp(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{
   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);
        

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*----------------- Calcul largeur ---------------------------------------------*/


   
    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
     else if (Z.indice+(Z.dejacoupe/LSTEP)<(numsteps-numstepsg))  
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }
     else if (Z.indice+(Z.dejacoupe/LSTEP)>=(numsteps-numstepsg)) 
    {
     large=(2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701);
    }




   /* ----------------------------------------------------------------------------*/


   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-(10/LSTEP))
      {
       double rank=Z.id_rang;
       produce  StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(37)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
      }
      else
      {
       produce StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(37)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
   
       
     /*if(ceil(Z.indice)==numsteps) { produce RollR(0)Down(90); }*/


       if (Z.indice>(numsteps-numstepsg))  /* Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-(10/LSTEP))
           {
           produce SetWidth(large)SetColor(37)Down(Torsiongaine*pow(-1,rank))F(LSTEP);
           }
           else
           { 
	   produce SetWidth(large)SetColor(37)F(LSTEP);
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       { 
        produce SetWidth(large)Left(Z.angleinsert)SetColor(37)F(LSTEP); 
       }

       if (Z.indice<(numsteps-numstepsg)) 
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-(10/LSTEP))
                  {
                   nproduce   StartGC()CurrentContour(FoldLimbo) SetWidth(large)SetColor(37)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                   
                   produce;
                  }
                  else
                  {
	            nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(37)F(LSTEP)EndGC(); 
                   
                    produce;
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               nproduce  StartGC() CurrentContour(Limbo3)Left(curve)SetWidth(large)SetColor(37)F(LSTEP)EndGC();  
               produce;
              }                 
          
         }
       }

   }

}


         

B(D):
{ 
  static char Assimil[30];
  sprintf(Assimil,"Assimilation : %4.2f", D.Assimibrut); 
  if(D.Assimibrut>2)
  {
   nproduce SB()SetWidth(0)F(8)SetColor(D.Assimibrut)Label(Assimil)EB();
  }
  produce SetColor(D.Assimibrut)Sphere(2);
}





Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/Down(angle) SetWidth(3)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/RollR(Rot-90)SetColor(14)SetWidth(2)F(0.3);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.3);
    }
  }
 
  produce;
}

Marqueur(onsenfout):
{
 produce SetColor(7) Sphere(1);
}



/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------GROUP LECTURE----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   

group LECTURE:


production:
/* Regeneration des tableau  *//* Alban */


Cp(Z)<<Cp(Z1):
{if(TPS==3){
 Z1.Assimibrut=Z1.Assimibrut+Z.Assimibrut;
 produce Cp(Z1);}
}


Cp(Z)<B(D) :
{

if (D.Taillegaine>Seuil[D.id_talle][D.id_plante]) 
  {
      Seuil[D.id_talle][D.id_plante]=D.Taillegaine;
  } 
  if (D.Taillegaine>Tailfingain[D.id_talle][D.id_rang][D.id_plante]) 
  {
     Tailfingain[D.id_talle][D.id_rang][D.id_plante]=D.Taillegaine;
  }
  if (D.id_talle>indicetalle[D.id_plante])
   {
     indicetalle[D.id_plante]=D.id_talle+1;
   }
   SeuilRef[D.id_talle][D.id_rang][D.id_plante]=12;  
   if((D.id_talle+1)>=NBtalles[D.id_plante]){NBtalles[D.id_plante]=D.id_talle+1;}
   D.Assimibrut=Z.Assimibrut;
   
   produce B(D);

}




decomposition:

Cdecomp(age, Agecroiss, Taillefeuille,id_rang, id_talle, id_plante, Ymax, Difftps, Reservoir, Taillefinalel, Taillefinaleg, Taillelimbe, Taillegaine, Phase, indice,rapportK, status, angleinsert, dejacoupe, geno, assimilation,v):
 {

   int numsteps=ceil((Taillefeuille-dejacoupe)/LSTEP);
   int numstepsg=ceil(Taillegaine/LSTEP);

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Taillelimbe;
  L0=Taillefinalel-taille;   

 
  

   /*---------------------------------------------- Calcul largeur ---------------------------------------*/


    if ((indice+(dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Taillelimbe/LSTEP)-(indice+(dejacoupe/LSTEP)))*LSTEP))/Taillefinalel), 2)+0.66*((L0+((Taillelimbe/LSTEP)-(indice+(dejacoupe/LSTEP)))*LSTEP)/Taillefinalel)+0.903)*(((1.9645*log(double(Ymax)))-5.701)/2);
    }
    else
    {
      double largini=(-1.38*pow(((L0+((((Taillelimbe/LSTEP)-(indice+(dejacoupe/LSTEP)))-1)*LSTEP))/Taillefinalel), 2)+0.66*((L0+((((Taillelimbe/LSTEP)-(indice+(dejacoupe/LSTEP)))-1)*LSTEP))/Taillefinalel)+0.903)*(((1.9645*log(double(Ymax)))-5.701)/2);
      large=largini-((((Taillelimbe/LSTEP)-(indice+(dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }

  /*------------------------------------------ Calcul largeur ---------------------------------------*/




  float Assimibrut;
  Assimibrut=LSTEP*2*large*lumiere[int(id_plante)][int(id_rang)][int(indice)][int(id_talle)];  // Doit convertir en Integer
  Feuilphoto Z={age, Agecroiss, Taillefeuille,id_rang, id_talle, id_plante, Ymax, Difftps, Reservoir, Taillefinalel, Taillefinaleg, Taillelimbe, Taillegaine, Phase, indice,rapportK, status, angleinsert, dejacoupe, geno, Assimibrut,v}; // Assimibrut=lumiere recue de Canestra via visualea ; ficheir lumiere
  produce Cp(Z);                                                  
 }




Bdecomp(age,Agecroiss,Taillefeuille,id_rang,id_talle,id_plante,Ymax,Difftps,Reservoir,Taillefinalel,Taillefinaleg,Taillelimbe,Taillegaine,Phase,longueur,cutstatus,rapportK,coupe, angleinsert, dejacoupe, geno, Assimibrut, Besoinencroiss):
{
  FeuilTrans D={age,Agecroiss,Taillefeuille,id_rang,id_talle,id_plante,Ymax,Difftps,Reservoir,Taillefinalel, Taillefinaleg, Taillelimbe, Taillegaine, Phase,longueur,cutstatus,rapportK,coupe, angleinsert, dejacoupe, geno, Assimibrut, Besoinencroiss};
  produce B(D);
}

Entrendecomp(val, id_rang, id_talle, id_plante, angletal, Assimdescend, Demandescend, Assimmonte, Demandmonte):
{
 Internode Creat={val, id_rang, id_talle, id_plante, angletal, Assimdescend, Demandescend, Assimmonte, Demandmonte};
 produce Entrenoeud(Creat);
}











/*------------------------------------ LECTURE INTERPRETATION ---------------------------------------------------------------------------------------*/


interpretation:


   Cp(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{

   static char assm[80];
    sprintf(assm,"Ass: %3.3f ",Z.Assimibrut);
   if(Z.id_rang==4 && Z.indice==10 && Z.id_talle==0) {nproduce Label(assm);}









   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);
   

   
   float color=(Z.Assimibrut*100);
   
   
  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*----------------- Calcul largeur ---------------------------------------------*/


   
    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
     else if (Z.indice+(Z.dejacoupe/LSTEP)<(numsteps-numstepsg))  
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
     }
     else if (Z.indice+(Z.dejacoupe/LSTEP)>=(numsteps-numstepsg)) 
    {
     large=(2/3.1415)*(1.9645*log(double(Z.Taillegaine*5))-5.701); 
    }

  


   /* ----------------------------------------------------------------------------*/


   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-(10/LSTEP))
      {  
       double rank=Z.id_rang;
       produce  StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(color)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
       
      }
      else
      {
       produce StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(color)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
   
       if (Z.indice>(numsteps-numstepsg))  /* Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-(10/LSTEP))
           {
           produce SetWidth(large)SetColor(color)Down(Torsiongaine*pow(-1,rank))F(LSTEP);
           }
           else
           { 
	   produce SetWidth(large)SetColor(color)F(LSTEP);
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       { 
        produce SetWidth(large)Left(Z.angleinsert)SetColor(color)F(LSTEP); 
       }

       if (Z.indice<(numsteps-numstepsg)) 
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-(10/LSTEP))
                  {
                   nproduce   StartGC()CurrentContour(FoldLimbo) SetWidth(large)SetColor(color)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                   
                   produce;
                  }
                  else
                  {
	            nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(color)F(LSTEP)EndGC(); 
                   
                    produce;
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               nproduce  StartGC() CurrentContour(Limbo3)Left(curve)SetWidth(large)SetColor(color)F(LSTEP)EndGC();  
               produce;
              }                 
          
         }
       }

   }

}


         

B(D):
{ 
  static char Assimil[30];
  sprintf(Assimil,"Assimilation : %4.2f Demande Croiss: %2.3f", D.Assimibrut, D.Besoinencroiss); 
  if(D.Assimibrut>0)
  {
   nproduce SB()SetWidth(0)F(8)SetColor(D.Assimibrut)Label(Assimil)EB();
  }
  produce SetColor(D.Assimibrut)Sphere(2);
}



/*

Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce Down(angle) SetWidth(0)SetColor(21)F(0);
    }
    if (i==2)
    {
     nproduce  RollR(Rot-90)SetColor(21)SetWidth(0)F(0);
    }
    else
    {
     nproduce SetColor(21)SetWidth(0)F(0);
    }
  }
 
  produce;
}

Marqueur(onsenfout):
{
  produce SetColor(14) Sphere(1);
}

Affichage(para):
{ 


produce;  
}

*/


/*---------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------GROUP ALLOCATION----------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------*/   
group ALLOCATION :



production:


Pool(ParaCpool,ParaCdemand) < Reserve(Contenu):
{
 if(ParaCpool>ParaCdemand /*&& TPS==Photosynthese+2*/)
  {
   produce Reserve(Contenu+(ParaCpool-ParaCdemand)); 
  }
  else
  {
   produce Reserve(Contenu);
  }
}
/* Mise en reserve d'Excedent */


Pool(ParaCpool,ParaCdemand):
{
 if(ParaCpool>=ParaCdemand)
  {
   produce Pool(ParaCdemand,ParaCdemand);
  }   
 if(ParaCpool<ParaCdemand)
  {
   produce Pool(ParaCpool,ParaCdemand);
  }
}

Pool(ParaCpool,ParaCdemand)Centralerect(para,para2, para3,para4,para5)GetHead(Aa,Bb, Cc) << Entrenoeud(X1):
{                             /* Si excedent alors montée=demande vu que excedent en reserce, sinon montée=dispo */
   X1.Assimmonte=ParaCpool;
   X1.Demandmonte=ParaCdemand;
  produce Entrenoeud(X1);
}


Entrenoeud(X1)Right(angle)Centralerect(para,para2, para3,para4,para5)GetHead(Aa,Bb, Cc) << Entrenoeud(X2):
{
 X2.Assimmonte=X1.Assimmonte;
 X2.Demandmonte=X1.Demandmonte;
 produce Entrenoeud(X2);
}

Entrenoeud(X1) << Entrenoeud(X2):
{
 X2.Assimmonte=X1.Assimmonte;
 X2.Demandmonte=X1.Demandmonte;
 produce Entrenoeud(X2);
}


Entrenoeud(X1) << Cp(Z) :  
{
 Z.Assimibrut=X1.Assimmonte;
 Z.Besoinencroiss=X1.Demandmonte;
 produce Cp(Z);
}


Cp(Z1) << Cp(Z2):     
{
   Z2.Assimibrut=Z1.Assimibrut;
   Z2.Besoinencroiss=Z1.Besoinencroiss;
   produce Cp(Z2);
}


Cp(Z) < B(D):
{
    D.Assimibrut=((D.Besoinencroiss/Z.Besoinencroiss)*Z.Assimibrut)/D.Besoinencroiss; /* Poportion de realisation du potentiel = AssimTot/Besointot*/
    D.Besoinencroiss=Z.Besoinencroiss; 
    produce B(D);  
}























interpretation:

 Cp(Z):    /* Copier-Coller de l'Interpretation de Feuille vu que doit avoir meme architecture */
{
   int numsteps=ceil((Z.Taillefeuille-Z.dejacoupe)/LSTEP);
   int numstepsg=ceil(Z.Taillegaine/LSTEP);
        

  float large;            /* Largeur de chaque element de surface */
  double L0;
  double taille;

  taille=Z.Taillelimbe;
  L0=Z.Taillefinalel-taille;   

 
  

   /*----------------- Calcul largeur ---------------------------------------*/


    if ((Z.indice+(Z.dejacoupe/LSTEP))>(Taillepointelineraire/LSTEP))
    {
     large=(-1.38*pow(((L0+(((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))*LSTEP)/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
    }
    else
    {
      double largini=(-1.38*pow(((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel), 2)+0.66*((L0+((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-1)*LSTEP))/Z.Taillefinalel)+0.903)*(((1.9645*log(double(Z.Ymax)))-5.701)/2);
      large=largini-((((Z.Taillelimbe/LSTEP)-(Z.indice+(Z.dejacoupe/LSTEP)))-((taille-Taillepointelineraire)/LSTEP))*((largini*LSTEP)/Taillepointelineraire)); 
    }
    


/* ----------------------------------------------------------------------------*/


   if (Z.Phase==0)
   { 
      if(Z.indice>=numsteps-10)
      {
       double rank=Z.id_rang;
       produce  StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
      }
      else
      {
       produce StartGC()CurrentContour(Limbo3)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
      }
   }

   if (Z.Phase==1) 
   {
    
       if (Z.indice>(numsteps-numstepsg))  /* Car indice part de bout de feuille */
       {
         double rank=Z.id_rang;
           if(Z.indice>=numsteps-10)
           {
           produce SetWidth(large)SetColor(2)Down(Torsiongaine*pow(-1,rank))F(LSTEP);
           }
           else
           { 
	   produce SetWidth(large)SetColor(2)F(LSTEP);
           } 
       }

       if (ceil(Z.indice)==(numsteps-numstepsg) && Z.Taillegaine>Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]-2)
       {
        produce SetWidth(large)Left(Z.angleinsert)SetColor(1)F(LSTEP); 
       }

       if (Z.indice<(numsteps-numstepsg)) 
       {

         if (Z.age<=8000)
         {    
              if (Z.indice>(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
                 double rank=Z.id_rang;
                 if(Z.indice>=numsteps-10)
                  {
                   nproduce   StartGC()CurrentContour(FoldLimbo) SetWidth(large)SetColor(23)Down(Torsiongaine*pow(-1,rank))F(LSTEP)EndGC();
                   
                   produce;
                  }
                  else
                  {
	            nproduce StartGC()CurrentContour(FoldLimbo)SetWidth(large)SetColor(23)F(LSTEP)EndGC(); 
                   
                    produce;
                  }
              }  
              if (Z.indice<(numsteps-((Tailfingain[Z.id_talle][Z.id_rang-1][Z.id_plante]+1)/LSTEP)))
              {
               nproduce  StartGC() CurrentContour(Limbo3)Left(curve)SetWidth(large)SetColor(1)F(LSTEP)EndGC();  
               produce;
              }                 
          
         }
       }

   }

}


B(D):
{
  static char Remonte[30];
  sprintf(Remonte,"Dispo Tot: %6.2f Dmd Tot: %6.2f", D.Assimibrut,D.Besoinencroiss); 
  produce SetColor(D.Assimibrut)Sphere(2)f(5)Label(Remonte);
}


Affichage(para):
{   
 static char Alloc[30];
 sprintf(Alloc,"ALLOCATION"); 
 nproduce MoveTo(-120,40,0) SetColor(27)Label(Alloc);
 nproduce MoveTo(-100,50,0) SetWidth(3)F(10)SetWidth(6)StartGC()SetWidth(0.1)F(10)EndGC();
 produce;
}


Centralerect(age, rang,Rot, angle,geno):
{ 
  for (int i=0; i<2; i++)
  {
    if (i==1)
    {
     nproduce /*RollToVert()*/Down(angle) SetWidth(3)SetColor(14)F(0);
    }
    if (i==2)
    {
     nproduce  /*RollToVert()*/RollR(Rot-90)SetColor(14)SetWidth(2)F(0.3);
    }
    else
    {
     nproduce SetColor(14)SetWidth(2)F(0.3);
    }
  }
 
  produce;
}

Racine(contenu):
{
produce f(-60)SetColor(15) Sphere(9)f(60);
}

Reserve(Contenu):
{ 
 static char Reserv[30];
 sprintf(Reserv,"En reserve: %4.2f", Contenu); 
 
 produce SetWidth(10)SetColor(33)f(-60)F(20)Label(Reserv);
 
}


Pool(Ccontenu,demande):
{
       static char Assim[30];
       sprintf(Assim,"Ass TOTALE: %6.2f Demand TOTALE: %6.2f", Ccontenu,demande); 
       nproduce SB()Down(180)SetColor(Ccontenu/5)SetWidth(4)F(30)Sphere(/*Ccontenu*/10)SetWidth(0)F(15)Label(Assim)EB();
       produce;
}










endgroup

































