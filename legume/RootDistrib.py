#########
## G. Louarn
## 25/06/2015
##
## calcul des distributions de longueur de racine dans une grille de sol 3D
## systeme racinaire representes comme une liste de cylindres de densite homogene
## calcul approxime aux bounding box de chaque cylindre (A faire : ponderer par distance au centre)
#########


#from numpy import linspace
#from scipy import array,ones,sum, shape, reshape
import numpy as np
from copy import deepcopy

#from rpy import r
#r.matrix.local_mode(0)


def lims_soil(pattern8, dxyz, unit='cm'):
    # pattern en cm ou m;dxyz de l'objet sol en m
    #! suppose  les espacement reguliers
    # lims en unite unit du pattern
    
    if unit=='cm':
        cor_unit=100.
    elif unit=='m':
        cor_unit=1.

    xlims = np.linspace(pattern8[0][0], pattern8[1][0],len(dxyz[0])+1)
    ylims = np.linspace(pattern8[0][1], pattern8[1][1],len(dxyz[1])+1)
    zlims = np.linspace(0, sum(dxyz[2])*cor_unit,len(dxyz[2])+1)
    lims = [xlims,ylims,zlims]
    return lims 
    # a mettre en fonction dans module sol?
    #a reprendre pour que ne soit pas que sur grille reguliere! (ajout de dx_i successifs)

#lims = lims_soil(pattern8, dxyz)

def BBOX(p,r,H):
    "bounding box d'un cylindre avec p (centre face inferieure), rayon r, hauteur H"
    p1 = p - np.array([r,r,0])
    p2 = p + np.array([r,r,H])
    return p1, p2

def VolBBOX(p1,p2):
    "volume d'un box"
    return (p2[0]-p1[0])*(p2[1]-p1[1])*(p2[2]-p1[2])


#pmin,pmax = BBOX(array([0.,0.,0.]),r=1.5,H=3.)


def cor_points(pmin,pmax,lims):
    """revoie points d'extremite  d'un segment   corrige pour rester dans une grille finie (hypothese de couvert infini)"""
    pcormin, pcormax = deepcopy(pmin), deepcopy(pmax)
    #coord X
    if (pmax[0]-pmin[0])>= (lims[0][-1]-lims[0][0]):#segment plus grand que largeur max-> prend les extremite de la grille
        pcormin[0] = lims[0][0]
        pcormax[0] = lims[0][-1]
    elif lims[0][0]<=pmin[0]<lims[0][-1] and lims[0][0]<=pmax[0]<lims[0][-1]: #tous deux compris dans la grille
        pass #pcor deja OK
    elif lims[0][0]>pmin[0]: #pmin inferieur a borne inf mais inferieur a largeur max ->miroir
        delta = lims[0][0]-pmin[0]
        pcormin[0] = lims[0][-1] - delta%(lims[0][-1]-lims[0][0])
    elif pmax[0]>=lims[0][-1] :#pmax sup a borne sup mais inferieur a largeur max -> miroir
        delta = pmax[0]-lims[0][-1]
        pcormax[0] = lims[0][0] + delta%(lims[0][-1]-lims[0][0])

    #coord Y
    if (pmax[1]-pmin[1])>= (lims[1][-1]-lims[1][0]):#segment plus grand que largeur max-> prend extremite du segment
        pcormin[1] = lims[1][0]
        pcormax[1] = lims[1][-1]
    elif lims[1][0]<=pmin[1]<lims[1][-1] and lims[1][0]<=pmax[1]<lims[1][-1]: #tous deux compris dans la grille
        pass #pcor deja OK
    elif lims[1][0]>pmin[1]: #pmin inferieur a borne inf
        delta = lims[1][0]-pmin[1]
        pcormin[1] = lims[1][-1] - delta%(lims[1][-1]-lims[1][0])
    elif pmax[1]>=lims[1][-1] :#pmax sup a borne sup
        delta = pmax[1]-lims[1][-1]
        pcormax[1] = lims[1][0] + delta%(lims[1][-1]-lims[1][0])

    #Coord z : devrait etre gere en amont en limitant taille des cylindre -> pas traite, mais peut gerer en forcant aux bornes sup et inf
    return pcormin, pcormax

#a faire: 3 segments selon x, y, z de proportion par voxel selon x, y ,z (une fois que sait que point dans la grille)
#2 cas selon sin min> max ou pas
#! ponderation de differentes bbox selon leur volume initial, pas selon points corriges!

def frac_voxelsBBox(pmin,pmax,lims):
    """ calcul les fractions de voxel utilise selon les trois axes x,y,z - pmin et pmax definisse une BBOX inclue dans la grille ; lims = limite des voxels selon x,y,z"""
    #pmin et pmax deja corrige par cor_points
    vv = []
    for axe in range(3):#= 2
        if pmin[axe]<=pmax[axe]: #verifie si pmin est plus petit que pmax (ou est le debut de la BBOX?)
            in_ = 0
        else:
            in_ = 1

        v = []
        for i in range(0, len(lims[axe])-1):
            frac = 0
            if lims[axe][i]<=pmin[axe]<lims[axe][i+1] and in_==0: #debut segment
                frac = (lims[axe][i+1] - pmin[axe]) / (lims[axe][i+1] - lims[axe][i])  
                in_=1

            if lims[axe][i]<=pmax[axe]<lims[axe][i+1] and in_==1: #fin segment
                if frac>0.: #point de debut dans le meme voxel
                    frac = 1-(pmax[axe] - lims[axe][i]) / (lims[axe][i+1] - lims[axe][i])
                else:
                    frac = (pmax[axe] - lims[axe][i]) / (lims[axe][i+1] - lims[axe][i])
                in_=0
            elif (lims[axe][i]<=pmin[axe]<lims[axe][i+1])!=True and (lims[axe][i]<=pmax[axe]<lims[axe][i+1])!=True and in_==0:#pas de point ds le voxel et segent pas commence
                frac =0.
            elif (lims[axe][i]<=pmin[axe]<lims[axe][i+1])!=True and (lims[axe][i]<=pmax[axe]<lims[axe][i+1])!=True and in_==1:#pas de point ds le voxel et segent commence
                frac =1.

            v.append(frac)
        vv.append(np.array(v))
    return vv


#pmin,pmax = BBOX(array([0.,0.,0.]),r=3.5,H=3.)
#fracs = frac_voxelsBBox(pmin,pmax,lims)


def fracBBOX(fracs, m1):
    """calcule la matrice des fractions"""
    m = deepcopy(m1)*1.
    for z in range(m.shape[0]):
        m[z,:,:] = m[z,:,:]*fracs[2][z]

    for y in range(m.shape[1]):
        m[:,y,:] = m[:,y,:]*fracs[1][y]

    for x in range(m.shape[2]):
        m[:,:,x] = m[:,:,x]*fracs[0][x]
    return m
    #m1 pourrait etre calcule selons les lims (comme dans updateRootDistrib)
    # a optimiser: renvoer des ous matrice avec voxels extremes plutot que grille complete
    # -> temps de calcul et memoire pour les grosses scenes?



def updateRootDistrib(RLtot, syst_rac, lims, optNorm=0):
    """ Distribution de longueur totale de racine (RLtot, m) dans grille (lims) pour une liste de cylindre decrivant le volume d'une systeme racinaire"""
    # for a single root system
    m_1 = np.ones([len(lims[2])-1, len(lims[1])-1,len(lims[0])-1])#matrice equivalente d'un objet sol
    if syst_rac==[]:#secondaires pas develope
        mL = deepcopy(m_1)*0.
        #ou met-on RLtot?
        #pass 
        #a gerer -> pour renvoyer mL (voir ancienne version) -> une facon de gerer = mettre un rectangle pour le pivot -> jamais vide du coup
    else:
        #calcul de volumes de BBox et des fraction de voxel occupee par chacun
        Vols = []
        ls_fracs = []
        for i in range(len(syst_rac)):
            pmin,pmax = BBOX(np.array([syst_rac[i][0], syst_rac[i][1], syst_rac[i][2]]),r=syst_rac[i][3],H=syst_rac[i][4])
            pminc,pmaxc =cor_points(pmin,pmax,lims)
            fracs = frac_voxelsBBox(pminc,pmaxc,lims)
            ls_fracs.append(fracBBOX(fracs, m_1))
            Vols.append(VolBBOX(pmin,pmax))

        #calcul des Longueur de racine par voxel
        frac_roots = np.array(Vols)/sum(Vols) #fraction proportionnelle au volume => densite constante sauf dans zones d'overlap
        L_roots = frac_roots*RLtot

        mL = deepcopy(m_1)*0.
        for i in range(len(ls_fracs)):
            mL += ls_fracs[i]*L_roots[i]/(sum(ls_fracs[i])+10e-15)

    #si renvoi distrib normalisee
    if optNorm==1:
        mL = mL / max(np.sum(mL),10e-15)#

    return mL

def calc_ls_roots_fromNorm(ls_rootsN, RLtot):
    """ to update RLtot absolute distribution keeping the same relative distributions """
    # for a list of root systems
    new_ls_roots = deepcopy(ls_rootsN)
    for nump in range(len(new_ls_roots)):
        new_ls_roots[nump] = new_ls_roots[nump] * RLtot[nump] #invar['RLTotNet'][nump] * 100

    return new_ls_roots


def propRootDistrib(ls_roots):
    """ proportion of plant root length in each voxel for a list of root system grids"""
    # for a list of root systems
    ls_props = []
    for nump in range(len(ls_roots)):
        mat_ =  ls_roots[nump]/max(sum(ls_roots[nump]),10e-15)
        ls_props.append(mat_)

    return ls_props
    #rq: inclu dans propRootDistrib_upZ

def propRootDistrib_upZ(ls_roots, depth=None, dz_sol=5.):
    """ proportion of plant root length in each voxel up to depth Z for a list of root system grids"""
    # for a list of root systems
    ls_props = []
    for nump in range(len(ls_roots)):
        if depth is None: #pas de profondeur renseignee = prend tout le profil
            mat_ =  ls_roots[nump]/max(sum(ls_roots[nump]), 10e-15)
            ls_props.append(mat_)
        elif depth <= dz_sol: #premiere couche seulement
            mat_ini = ls_roots[nump]
            mat_keep = mat_ini[0,:,:]
            mat_2 = mat_ini*0.
            mat_2[0,:,:] = mat_keep
            mat_ = mat_2 / max(np.sum(mat_2), 10e-15)
            ls_props.append(mat_)
        else:
            mat_ini = ls_roots[nump]
            row_to_keep = min(int(depth / dz_sol) + 1, mat_ini.shape[0])
            mat_keep = mat_ini[0:row_to_keep, :, :]
            mat_2 = mat_ini * 0.
            mat_2[0:row_to_keep, :, :] = mat_keep
            mat_ = mat_2 / max(np.sum(mat_2), 10e-15)
            ls_props.append(mat_)

    return ls_props
    #ou defaut=30 cm?



def VoxWithRoots(ls_roots):
    """ 1 in voxels with plant roots for a list of root system grids """
    ls_ones = []
    for nump in range(len(ls_roots)):
        mat_1 = deepcopy(ls_roots[nump])
        mat_1[mat_1>0] = 1
        ls_ones.append(mat_1)

    return ls_ones

# def VisuRootDistrib2D(RLmap, lims, zlim=None):
#     #carte 2D
#     x = lims[1][0:-1] + (lims[1][1]-lims[1][0])/2.
#     y = -lims[2][0:-1]
#     y.sort()
#     mirorRLmap = deepcopy(RLmap)
#     for i in range(shape(RLmap)[0]):#remet a l'endroit
#         mirorRLmap[-1-i,:] = RLmap[i,:]
#
#     z= r.matrix(reshape(mirorRLmap, shape(RLmap)[0]*shape(RLmap)[1],1), shape(RLmap)[1], shape(RLmap)[0], byrow=True)
#     if zlim!=None:
#         r.filled_contour(x,y,z, zlim=zlim)#echelle des z donnee
#     else:
#         r.filled_contour(x,y,z)#defaut: echelle en relatif
#     #Afaire : exprimer en densite (rediviser par les volumes de voxel
#     #autre couleur map?
#
# #retire car dependence a r et rpy -> passer en rpy2 et dans un autre fichier si veut l'utiliser


def build_ls_roots_mult(RLTot, dic_systrac, lims, optNorm=0):
    """ """
    #si optNorm == 1, renvoie des distributions normalisees
    #build_ls_roots_mult?? ou ca?
    ls_roots=[]
    for nump in range (len(RLTot)):#s'assurer que nump pris dans l'ordre??
        ls_roots.append(updateRootDistrib(RLTot[nump], dic_systrac[nump], lims, optNorm))

    return ls_roots

    

def convd(d):
    """ pour recupere sortie d'un dico profil en liste"""
    ks = list(d.keys())
    ks.sort()
    res=[]
    for k in ks: res.append(d[k])
    return res


#plus utilise pour racines
#def filtre_ratio(ratio, tresh=1.):
#    """ pour filtrer les valeur < a un seuil - array de 0 et 1 """
#    res = ratio>=tresh
#    return res*1.


#plus utilise pour racine
#def updateRootLenprofil(RLtot_, RprospectProfil_, RLProfil_):
#    if sum(RprospectProfil_)==0.:#secondaires pas developpees
#        RLProfil_[0] = RLtot_ #toute la surface d'echange dans le premier horizon
#    else:
#        vols = array(RprospectProfil_)*array(RprospectProfil_)
#        props = vols*1./sum(vols).tolist()
#        for k in range(len(RprospectProfil_)):
#            RLProfil_[k] = RLtot_*props[k]
#
#    return RLProfil_


##Test: 1) produire liste de Root-system avec le modele L-egume et les visualiser dans une grille 2D

##lecture d'un profil racinaire produit par L-system
#import IOtable

#f = file(r'H:\devel\grassland\grassland\L-gume\lsAxes.csv', 'r')#test_rac1.csv
#tab = IOtable.table_csv(f)
#f.close()

#MSrac_fine = 2. #g
#SRL = 200. #m.g-1
#Ltot = MSrac_fine*SRL#m

#pattern8 = [[-50,-50], [50.,50.]]
#Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
#dz_sol = 5. #cm
#ncouches_sol = 20

#dxyz = [[Lsol/10.]*10, [largsol/10.]*10, [dz_sol/100.]*ncouches_sol] #m
#lims = lims_soil(pattern8, dxyz)

#syst_rac = tab
#RL = updateRootDistrib(Ltot, syst_rac, lims)
#RLmap = sum(RL, axis=2) #some selon axe des X
#VisuRootDistrib2D(RLmap, lims,zlim=[0,12])
##VisuRootDistrib2D(RLmap, lims)


#from openalea.plantgl.all import *
#import sys
#path_ = r'H:\devel\grassland\grassland\luzerne'
#path2_ = r'H:\devel\grassland\grassland\L-gume'
#sys.path.insert(0, path_)
#sys.path.insert(0, path2_)
#from Obj3Dutils import *

#def VisuRootDistrib3D(RL, lims, dxyz, MaScene = Scene()):
#    #plot de distribution de densite relative par voxel
#    bb =Box(Vector3(1.,1.,1.))
#
#    #MaScene = Scene()
#    MonViewer = Viewer
#    for z in range(RL.shape[0]):
#        for y in range(RL.shape[1]):
#            for x in range(RL.shape[2]):
#                bb2 = transformation(bb,dxyz[0][0],dxyz[1][0],dxyz[2][0],0, 0, 0,lims[0][x]/100.,lims[1][y]/100.,-lims[2][z]/100.)
#                if RL[z,y,x]>0.:
#                    col_ = couleur (4, 1-RL[z,y,x]/RL.max())
#                    bb2s = Shape(bb2, Material(Color3(col_[0],col_[1],col_[2]), transparency=0.65))
#                    MaScene.add(bb2s)
#
#    MonViewer.display(MaScene)








##Test 2) avec sol 1D simplifier 
#pattern8 = [[-50,-50], [50.,50.]]
#Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
#dz_sol = 5. #cm
#ncouches_sol = 20

#dxyz = [[Lsol/10.]*1, [largsol/10.]*1, [dz_sol/100.]*ncouches_sol] #m
#lims = lims_soil(pattern8, dxyz)

#syst_rac = tab
#RL = updateRootDistrib(Ltot, syst_rac, lims)
#RLmap = sum(RL, axis=2) #some selon axe des X

##calcul de RL OK
##pb avec la visu2D (mais c'est du 1D!)
##VisuRootDistrib2D(RLmap, lims,zlim=[0,12])
#r.plot(RLmap[:,0], lims[2][0:-1])


##A faire
##pour ponderer segment par distance au centre -> integrer une fonction trigo?
#from scipy.integrate import quad, dblquad
#from scipy import sqrt, pi


##fonction sin -> rapide mais pas equivalent a vrai cyclindre (moins dense au centre
#I = quad(cos, -3.14/2, 3.14/2) #aire de la derivee de sin = cos


##double integrale pour aire d'un vrai cercle/cylindre, mais lent et ne fonctionne pas (retombe pas sur 3.14?)
#I = dblquad(lambda x,r: sqrt(r**2-x**2), -1, 1., lambda x: 0, lambda x: 1.)


#def sin_demi(x):
#    #fonction qui ramene en -r et r pour 0,pi ressemble beaucoup a un cercle (pourrait optimiser l'exposant)
#    return (sin(x))**0.5

#I = quad(sin_demi, 0, pi)[0]

##complication: gerer les overlap pour bien appliquer les ponderations lorsque deborde... si gere comme avec bbox-> cree discontinuite
##ponderation a calculer a partir des vv de frac_voxelsBBox
