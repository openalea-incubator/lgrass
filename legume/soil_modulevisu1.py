from openalea.plantgl.all import *
import sys
path_ = r'C:\devel\l-egume\legume'
sys.path.insert(0, path_)
from Obj3Dutils import *


def couleur (echelle, indice):
    """renvoie RGB pour differentes echelles en fonction d'un indice entre 0 et 1"""
    # echelle 1 : Bleu au blanc
    if echelle == 1 :
        R = int(255*indice)
        G = int(255*indice)
        B = 255

    # echelle 2 : Noir au blanc en gradient de bleu
    if echelle == 2:
        if indice <0.2:
            R = 0
            G = 0
            B = int(255*indice/0.2)
        elif indice >= 0.2:
            R = int(255*(indice-0.2)/0.8)
            G = int(255*(indice-0.2)/0.8)
            B = 255

    # echelle 3 : Bleu au rouge
    if echelle == 3:
        if indice <0.4:
            R = int(170*indice/0.4)
            G = int(170*indice/0.4)
            B = 255
        elif indice >= 0.4 and indice < 0.6:
            R = int(170+(255-170)*(indice-0.4)/(0.6-0.4))
            G = 170
            B = int(255-(255-170)*(indice-0.4)/(0.6-0.4))
        elif indice >0.6 :
            R = 255
            G = int(170*(1-indice)/0.4)
            B = int(170*(1-indice)/0.4)

    # echelle 4 : noir au blanc en gradiend de jaune
    if echelle == 4 :
        if indice <0.7:
            R = int(255*indice/0.7)
            G = int(255*indice/0.7)
            B = 0
        elif indice >= 0.7:
            R = 255
            G = 255
            B = int(150*(indice-0.7)/0.3)

    # echelle 3 : Rouge au bleu
    if echelle == 5:
        indice = 1.-indice
        if indice <0.4:
            R = int(170*indice/0.4)
            G = int(170*indice/0.4)
            B = 255
        elif indice >= 0.4 and indice < 0.6:
            R = int(170+(255-170)*(indice-0.4)/(0.6-0.4))
            G = 170
            B = int(255-(255-170)*(indice-0.4)/(0.6-0.4))
        elif indice >0.6 :
            R = 255
            G = int(170*(1-indice)/0.4)
            B = int(170*(1-indice)/0.4)

    return R,G,B


def plot_soil_properties (S, vals, MaScene=Scene(), col_scale=5):#dxyz, m_soil_vox, asw_t pourraient etre remplacee par objet sol
    """ S=soil object; vals = matrice 3D de valeur entre 0 et 1 de propriete de sol a visualiser / e.g S.ftsw_t """
    bx = Box(Vector3(1.,1.,1.))
    for z in range(len(S.dxyz[2])):
        for x in range(len(S.dxyz[0])):
            for y in range(len(S.dxyz[1])):
                dims = [S.dxyz[0][x]/2., S.dxyz[1][y]/2., S.dxyz[2][z]/2.]
                p_ini = S.m_soil_vox[z][x][y]
                col = couleur (col_scale, max(0.,vals[z][x][y]))
                b = transformation(bx, dims[0], dims[1], dims[2], 0,0,0, p_ini[0], p_ini[1], -p_ini[2]-S.dxyz[2][z]/2.)
                MaScene.add(Shape(b, Material(Color3(col[0],col[1],col[2]), transparency=0.65)))

    return MaScene

