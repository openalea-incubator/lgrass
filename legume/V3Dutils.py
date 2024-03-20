### set de fonctions utiles pour manipuler les coordonnees 3D et les vecteurs

#import scipy
import numpy as np

def XyzToPol (coordxy) :
    """ converti les coordonnees carthesiennes d'un point (x,y,z) en coordonnees polaires (r,azi,incli)"""
    x,y,z = coordxy[0], coordxy[1], coordxy[2]
    xy2 = x*x+y*y
    r = np.sqrt(xy2+z*z)
    if r==0 :
        incli =0
    else :
        incli = np.arcsin(z/r)
    if (x==0 and y==0):
        azi = 0
    elif (y>=0) :
        azi = np.arccos(x/np.sqrt(xy2))
    else :
        azi = -np.arccos(x/np.sqrt(xy2))
    return np.array([r,azi,incli])


def PolToXyz (coordpol) :
    """ converti les coordonnees polaires (r,azi,incli) d'un point en coordonnees carthesiennes(x,y,z)"""
    r,azi,incli = coordpol[0], coordpol[1], coordpol[2]
    z = r * np.sin(incli)
    l = np.sqrt (r*r-z*z)
    x = l * np.cos(azi)
    y = l * np.sin (azi)
    return np.array([x, y, z])

def RotateAxis (coordxy, r_azi, r_incli):
    """ calcule les nouvelles coord d'un point apres rotation autour des axes y (r_incli en radians) et z (r_azi en radians)"""
    #incli d'abord = rotation autour de y
    Pol_ini = XyzToPol (np.array([coordxy[0], coordxy[2], -coordxy[1]]))
    xyz_r = PolToXyz (np.array([Pol_ini[0],Pol_ini[1]+r_incli,Pol_ini[2]]))
    #azi ensuite = rotation autour de z
    Pol_sec = XyzToPol (np.array([xyz_r[0],-xyz_r[2],xyz_r[1]]))
    xyz_r2 = PolToXyz (np.array([Pol_sec[0],Pol_sec[1]+r_azi,Pol_sec[2]]))
    return xyz_r2

def Translate (coordxy, t):
    """ calcule les nouvelles coord d'un point apres trnaslation de t """
    if len(coordxy)==len(t):
        return coordxy + t
    else:
        print('vector lengths do not match')


#scipy.dot -> produit scalaire
#scipy.cross -> produit vectoriel
def produit_scalaire (u, v) :
    """produit scalaire de deux cecteurs u et v"""
    p = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    return p


def produit_vectoriel (u, v) :
    """produit vectoriel de deux cecteurs u et v"""
    p = np.array([ u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] ])
    return p

def normalised_v(vec):
    """ mise a 1 de la norme de vec """
    if vec[2] > 0. :
        z = np.sqrt((vec[2]*vec[2])/((vec[2]*vec[2])+(vec[1]*vec[1])+(vec[0]*vec[0])))
        y = z*vec[1]/vec[2]
        x = z*vec[0]/vec[2]
    elif vec[2] == 0. :
        vec[2] = 10e-12
        z = np.sqrt((vec[2]*vec[2])/((vec[2]*vec[2])+(vec[1]*vec[1])+(vec[0]*vec[0])))
        y = z*vec[1]/vec[2]
        x = z*vec[0]/vec[2]
    else :
        z = -np.sqrt((vec[2]*vec[2])/((vec[2]*vec[2])+(vec[1]*vec[1])+(vec[0]*vec[0])))
        y = z*vec[1]/vec[2]
        x = z*vec[0]/vec[2]

    return np.array([x,y,z])


def norme_v(vec):
    """ calcule la norme d'un vecteur """
    return np.sqrt((vec[2]*vec[2])+(vec[1]*vec[1])+(vec[0]*vec[0]))

def distance(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)

def plane_eq(n, p):
    """compute the parameters a b c d (ax + by, +cz +d =0) of a plane defined by its normalised normal n and a point p"""
    a, b, c = n[0], n[1], n[2]
    d = -n[0]*p[0]-n[1]*p[1]-n[2]*p[2]
    return [a,b,c,d]

def intersec_D_plane(plane_par, v, p0):
    """ compute the intersection bewteen a plane and a line defined by p0,v - return -1 if there is no intersection"""
    a,b,c,d = plane_par
    n = np.array([a, b, c])
    ps = produit_scalaire (n, v)
    if ps == 0:
        print("plane and line colinear")
        return -1
    else:
        t = -(a*p0[0]+b*p0[1]+c*p0[2]+d)/ps#-(produit_scalaire (n, p0)+d)/ps
        if t<0:
            return p0+t*v#-1 #depend si considere demi droite a partir de p0 ou droite complete
        else:
            return p0+t*v



## test
#c=scipy.array([1.2,1,0.9])
#XyzToPol (c)
#PolToXyz (XyzToPol (c))
#RotateAxis (c, scipy.pi, scipy.pi)
#Translate (c, scipy.array([1.,1.,1.]))
#intersec_D_plane([0., 0., 1., 0.], array([0.,0.,1.]), array([1.,2.,1.]))

