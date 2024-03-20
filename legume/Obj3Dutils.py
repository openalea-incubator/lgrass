### set de fonctions utiles pour manipuler les geometry plantGL

import openalea.plantgl.all as pgl
from V3Dutils import *
import numpy as np

def mesh(geometry):
    """ renvoie le mesh d'une geometry"""
    tessel = pgl.Tesselator()
    geometry.apply(tessel)
    mesh_ = tessel.triangulation
    return mesh_


def tri(p1, p2, p3):
    """ renvoie le TriangleSet d'un triangle a partir de ses 3 points """
    points= pgl.Point3Array([ pgl.Vector3(p1[0],p1[1],p1[2]), pgl.Vector3(p2[0],p2[1],p2[2]), pgl.Vector3(p3[0],p3[1],p3[2])])
    indices= pgl.Index3Array([ pgl.Index3(0,1,2)])
    return pgl.TriangleSet(points, indices)

def quadform(p1, p2, p3, p4, opt=None):
    """ renvoie le TriangleSet a 4 points ; si opt != None inverse les points pour representation des triangles"""
    points= pgl.Point3Array([ pgl.Vector3(p1[0],p1[1],p1[2]), pgl.Vector3(p2[0],p2[1],p2[2]), pgl.Vector3(p3[0],p3[1],p3[2]), pgl.Vector3(p4[0],p4[1],p4[2])])
    if opt == None:
        indices= pgl.Index3Array([ pgl.Index3(0,1,3), pgl.Index3(1,2,3)])
    else:
        indices= pgl.Index3Array([ pgl.Index3(0,1,2), pgl.Index3(0,2,3)])

    return pgl.TriangleSet(points, indices)


def turtle36():
    ## demi triangles hauts
    ## surface du turtle36 de 1 de rayon 2.9577
    #liste points
    ls_pt = []
    dazi = np.pi/6.
    for elv in [0, np.pi/6., np.pi/3.]:
        r = np.cos(elv)
        z = np.sin(elv)
        for i in range(12):
            teta = i*dazi
            x = r*np.cos(teta)
            y = r*np.sin(teta)
            ls_pt.append(pgl.Vector3(x,y,z))

    ls_pt.append(pgl.Vector3(0,0,1.))

    #liste d'index
    ls_id = []
    #1er ring
    for i in range(11):
        ls_id.append(pgl.Index3(i, i+12, i+13))

    ls_id.append(pgl.Index3(11, 11+12, 12))
    #2e ring
    for i in range(11):
        ls_id.append(pgl.Index3(i+12, i+24, i+25))

    ls_id.append(pgl.Index3(23, 23+12, 24))
    #3e ring
    for i in range(11):
        ls_id.append(pgl.Index3(i+24, i+25, 36))

    ls_id.append(pgl.Index3(35, 24, 36))
    
    # triangleset
    points= pgl.Point3Array(ls_pt)
    indices= pgl.Index3Array(ls_id)
    return  pgl.TriangleSet(points, indices) #36 tiangles en 3 rings


def transformation(obj, sx, sy, sz, rx, ry, rz, tx, ty, tz ): 
    """ Return a scaled, rotated and translated 3D object - Similar to 'transformation' in PovRay """ 
    s_obj = pgl.Scaled (pgl.Vector3(sx,sy,sz), obj)
    r_obj = pgl.EulerRotated (rx, ry, rz, s_obj)
    t_obj = pgl.Translated (pgl.Vector3(tx,ty,tz), r_obj)
    return t_obj


def conv_cyl(p1, p2, r):
    """ calcule des parametres pour le positionnement d'un cylindre a partir des coordonnees 
    des deux points extremes et du rayon """
    vec = XyzToPol (p2-p1)
    return p1, vec[0], r, vec[1], vec[2] #p1, longueur, rayon, azi, incli

def euler_normal(AA, BB, CC):
    """ compute the normal at the plane defined by euler angles AA, BB, CC """
    tri0 = tri(np.array([0.,0.,0.]), np.array([1.,0.,0.]), np.array([0.,1.,0.]))
    v = mesh(transformation(tri0, 1,1,1,AA+ np.pi/2.,CC,-BB,0,0,0))
    v.computeNormalList()
    return v.normalAt(0)

def tri_ortho(p1,p2,p3):
    """ calcul  orthocentre d'un triangle - meme resulat que .faceCenter(id) sur un triangle set"""
    #p4 milieu  p2-p3
    p4 = np.array(p2)+0.5*(np.array(p3)-np.array(p2))
    return np.array(p1)+2.*(p4-np.array(p1))/3.

def compute_ortho_list(ind_ls, pt_ls, epsilon = 0.001):
    """ calcule liste des orthocentre d'un triangle set """
    surf = compute_surface_list(ind_ls, pt_ls)
    ortho_ls = []
    for i in range(len(ind_ls)):
        if surf[i]>epsilon:
            p1, p2, p3 = pt_ls[ind_ls[i][0]], pt_ls[ind_ls[i][1]], pt_ls[ind_ls[i][2]] 
            ortho = tri_ortho(p1,p2,p3)
            ortho_ls.append(pgl.Vector3(ortho[0], ortho[1], ortho[2]))

    return ortho_ls

def compute_normal_list(ind_ls, pt_ls, epsilon = 0.001):
    """ calcule liste des normales d'un triangle set - contourne pb des nb normale different des nb faces avec .computeNormalList()"""
    surf = compute_surface_list(ind_ls, pt_ls)
    n_ls = []
    for i in range(len(ind_ls)):
        if surf[i]>epsilon:
            p1, p2, p3 = np.array(pt_ls[ind_ls[i][0]]), np.array(pt_ls[ind_ls[i][1]]), np.array(pt_ls[ind_ls[i][2]])
            pv = produit_vectoriel ((p2-p1), (p3-p1))
            n_ls.append(pv/norme_v(pv))

    return n_ls

def mesh_points(geometry):
    """ get the mesh points of a geometry """
    g = mesh(geometry)
    return list(map(np.array, g.pointList))

def triangle_area(p1, p2, p3):
    """ compute surface area of a triangle """
    u = p2-p1
    v = p3-p1
    return 0.5*norme_v(produit_vectoriel (u, v))

def compute_surface_list(ind_ls, pt_ls):
    """ calcule liste des surfaces d'un triangle set """
    s_ls = []
    for i in range(len(ind_ls)):
        p1, p2, p3 = np.array(pt_ls[ind_ls[i][0]]), np.array(pt_ls[ind_ls[i][1]]), np.array(pt_ls[ind_ls[i][2]])
        s_ls.append(triangle_area(p1, p2, p3))

    return s_ls


def leg_leaf(Lmax, largmax, alpha=0., gamma=0., unifol=0):
    gamma = gamma * np.pi / 180  # en radians
    lf, la, pe, br, crois = 21. / 21., 6.5 / 21., 10. / 21., 3.6 / 21., 0.1 / 21.  # leaf Trudeau modifie
    leaf = quadform(np.array([0., 0., 0.]), np.array([0.5, 0.5, 0.]), np.array([0., 1., 0.]), np.array([-0.5, 0.5, 0.]),
                    opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    up = transformation(leaf, 1, 1, 1, 0, 0, gamma, 0, pe / lf * Lmax, 0)
    right = transformation(leaf, 1, 1, 1, -np.pi / 180 * 80, 0, gamma, br / 2. * Lmax, crois * Lmax, 0)
    left = transformation(leaf, 1, 1, 1, np.pi / 180 * 80, 0, gamma, -br / 2. * Lmax, crois * Lmax, 0)
    if unifol == 0:
        return pgl.Group([up, right, left])  # groupe les differents geom
    else:
        return up  # 1 foliole pour la premiere feuille


def leg_leaf_lucas(Lmax, largmax, alpha=0., gamma=0., nfol=3, angfol=10., ecfol=6., anginit=45., geom=True, opt_Trud=0):  # angfol : pi/nombre de rangs n?ssaires pour boucler le demi-cercle // ecfol : longueur de rachis entre chaque paire de folioles (mm).
    nr = (nfol - 1) / 2 if nfol % 2 == 1 else nfol / 2  # nombre de paires de folioles lateraux
    angfol = 10  # 180/nr #demi-cercle complet obtenu sur le nombre de rangs de paires de folioles. Attention, ne marche que si l'ecartement ecfol est constant!
    gamma = gamma * np.pi / 180  # en radians
    anginit = anginit * np.pi / 180
    angfol = angfol * np.pi / 180
    ecfol = ((1.88 * nfol ** -0.54) * Lmax) * 10  # relation empirique entre nombre de folioles et ecartement entre deux rang? pour le sainfoin.

    ls_pts = []

    lf, la, pe, br, crois = 21. / 21., 6.5 / 21., 10. / 21., 3.6 / 21., 0.1 / 21.  # leaf Trudeau modifie
    leaf = quadform(np.array([0., 0., 0.]), np.array([0.5, 0.5, 0.]), np.array([0., 1., 0.]), np.array([-0.5, 0.5, 0.]),opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    if opt_Trud ==1: #Trudeau leaf
        alph = 15 #3.14 / 4.  #degre
        leaf = mesh_leaflet(lf, la, alph, 10)
        leaf = transformation(leaf, largmax / la, Lmax / lf,Lmax / lf, 0, 0, 0, 0, 0, 0)

    angup = (angfol * -(nr - 1)) - np.pi - anginit  # angle de placement du foliole central, au bout de la chaine
    up = transformation(leaf, 1, 1, 1, 0, 0, gamma, 0, (pe / lf * Lmax) + (ecfol * (np.cos(angup) + np.cos(anginit))), ecfol * (np.sin(angup) - np.sin(anginit)))
    ls_pts.append(np.array([0, (pe / lf * Lmax) + (ecfol * (np.cos(angup) + np.cos(anginit))), ecfol * (np.sin(angup) - np.sin(anginit))]))
    # right = transformation(leaf, 1,1,1,-3.14/180*80,0,gamma, br/2.*Lmax, crois*Lmax,0)
    # left = transformation(leaf, 1,1,1,3.14/180*80,0,gamma, -br/2.*Lmax, crois*Lmax,0)
    listfol = [up] if nfol % 2 == 1 else []
    for i in range(int(nr)):  # nombre de paires de folioles lateraux
        ang = (angfol * -i) - np.pi - anginit
        ecfolopp = (np.sin(ang) - np.sin(anginit)) * ecfol
        ecfoladj = (np.cos(ang) + np.cos(anginit)) * ecfol
        listfol.append(transformation(leaf, 1, 1, 1, -np.pi / 180 * 80, 0, gamma, br / 2. * Lmax, ecfoladj, ecfolopp))
        listfol.append(transformation(leaf, 1, 1, 1, np.pi / 180 * 80, 0, gamma, -br / 2. * Lmax, ecfoladj, ecfolopp))
        ls_pts.append(np.array([br / 2. * Lmax, ecfoladj, ecfolopp]))
        ls_pts.append(np.array([-br / 2. * Lmax, ecfoladj, ecfolopp]))

    if geom == True:
        return pgl.Group(listfol)  # groupe les differents geom
    else:
        return ls_pts


def geomstip(Lmax, largmax, alpha=0., gamma=0.):
    gamma = gamma * np.pi / 180  # en radians
    stip = quadform(np.array([0., 0., 0.]), np.array([0.5, 0.5, 0.]), np.array([0., 1., 0.]), np.array([-0.5, 0.5, 0.]),
                    opt=2)  # prends pas alpha en compte
    stip = transformation(stip, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    if Lmax >= largmax:
        right = transformation(stip, 1, 1, 1, -3.14 / 180 * alpha, -gamma, 0, 0, 0, 0)
        left = transformation(stip, 1, 1, 1, 3.14 / 180 * alpha, gamma, 0, 0, 0, 0)
    else:
        right = transformation(stip, 1, 1, 1, -3.14 / 180 * alpha, 0, gamma, 0, 0, 0)
        left = transformation(stip, 1, 1, 1, 3.14 / 180 * alpha, 0, gamma, 0, 0, 0)
    return pgl.Group([right, left])
    # par convention choisi longueur dans direction du petiole?
    # faire porter gamma sur la plus grande direction -> marche normalement pour pois
    # reprendre gammaFeuil pour le gamma (pas IncPet petiole)


def leg_grass(Lmax, largmax, gamma=0., angfol=10., nfol=8, anginit=45., geom=True):
    anginit = anginit * 3.14 / 180
    angfol = angfol * 3.14 / 180
    ecfol = Lmax / nfol  # longueur de segment de feuille

    ls_pts = []

    leaf = quadform(np.array([-0.5, 0., 0.]), np.array([-0.5, 1., 0.]), np.array([0.5, 1., 0.]), np.array([0.5, 0., 0.]),
                    opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, ecfol, 1., 0, 0, 0, 0, 0, 0)
    bottom = transformation(leaf, 1, 1, 1, 0, 0, 3.14 / 2 - anginit, 0, 0, 0)
    listfol = [bottom]

    # ajout points 2 premiers points debuts + 1er segments
    ls_pts.append(np.array([0, 0, 0]))
    ls_pts.append(np.array([0, np.sin(anginit) * ecfol, np.cos(anginit) * ecfol]))

    for i in range(1, nfol):  # nombre de segment restants
        ang = anginit - (angfol * -i)  # (angfol * -i) - 3.14 - anginit
        z = ls_pts[-1][2] + np.cos(ang) * ecfol
        y = ls_pts[-1][1] + np.sin(ang) * ecfol
        listfol.append(transformation(leaf, 1, 1, 1, 0, 0, 3.14 / 2 - ang, 0, ls_pts[-1][1], ls_pts[-1][2]))
        ls_pts.append(np.array([0, y, z]))
        #print i, ecfol, angfol, ang, distance(ls_pts[-1], ls_pts[-2])

    if geom == True:
        return pgl.Group(listfol)  # groupe les differents geom
    else:
        return ls_pts[1:]  # retire premier point -> surface attribuee aux extremite de segment

#pourrait reconstruire avec liste des points et liste des angles fournie
#ou pour visu seulement mettre un objet feuille groupe? (pas tous les segments?)


def leg_grass_withoutgeom(Lmax, largmax, gamma=0., angfol=10., nfol=8, anginit=45.):
    #fonction simplifiee de leg_grass qui renvoie que les points (aucun calcul de geometrie)
    anginit = anginit * np.pi / 180
    angfol = angfol * np.pi / 180
    ecfol = Lmax / nfol  # longueur de segment de feuille

    ls_pts = []
    ls_cos = []
    ls_sin = []

    # ajout points 2 premiers points debuts + 1er segments
    ls_pts.append(np.array([0, 0, 0]))
    cosi = np.cos(anginit)
    sini = np.sin(anginit)
    ls_pts.append(np.array([0, sini * ecfol, cosi * ecfol]))
    ls_cos.append(cosi)
    ls_sin.append(sini)

    for i in range(1, nfol):  # nombre de segment restants
        ang = anginit - (angfol * -i)  # (angfol * -i) - 3.14 - anginit
        cosi = np.cos(ang)
        sini = np.sin(ang)
        z = ls_pts[-1][2] + cosi * ecfol
        y = ls_pts[-1][1] + sini * ecfol
        ls_pts.append(np.array([0, y, z]))
        ls_cos.append(cosi)
        ls_sin.append(sini)

        # print i, ecfol, angfol, ang, distance(ls_pts[-1], ls_pts[-2])

    pointes = ls_pts[1:] #ce qui est renvoue par leg_grass # retire premier point -> surface attribuee aux extremite de segment

    return pointes, ls_pts , ls_cos, ls_sin
    #avec ecfol=1. , renvoie directe les sinus et cosinus


def larg_norm_trudeau(L):
    if L<0.996:
        return -12.268*L**4 + 22.958*L**3 -16.929*L**2 +6.2135*L #leaf Trudeau
    else:
        return 0

def larg_fol(Lrel, Lmax, largmax):
    "larg rel depuis base foliole"
    return larg_norm_trudeau(Lrel/Lmax)*largmax


def mesh_leaflet(Lmax, largmax, alpha=0., n=8):
    #liste de pts
    ls_pt = [pgl.Vector3(0.,0.,0.)]
    for i in range(1, n):
        Lrel = float(i)/float(n)
        l = larg_fol(Lrel, Lmax, largmax)
        ls_pt.append(pgl.Vector3(-l/2.*np.cos(alpha), Lrel*Lmax, l/2.*np.sin(alpha)))
        ls_pt.append(pgl.Vector3(0., Lrel*Lmax, 0.))
        ls_pt.append(pgl.Vector3(l/2*np.cos(alpha), Lrel*Lmax, l/2*np.sin(alpha)))

    ls_pt.append(pgl.Vector3(0., Lmax, 0.))

    #liste d'index
    ls_ind = [pgl.Index3(0,1,2), pgl.Index3(0,2,3)]
    for i in range(1, n):
        if i< n-1:
            ls_ind.append(pgl.Index3(i*3-2, (i+1)*3-2, (i+1)*3-1))
            ls_ind.append(pgl.Index3(i*3-1, i*3-2, (i+1)*3-1))
            ls_ind.append(pgl.Index3(i*3, i*3-1, (i+1)*3-1))
            ls_ind.append(pgl.Index3(i*3, (i+1)*3-1, (i+1)*3))
        elif i == n-1:
            ls_ind.append(pgl.Index3(i*3-1, i*3-2, i*3+1))
            ls_ind.append(pgl.Index3(i*3, i*3-1, i*3+1))

    return pgl.TriangleSet(pgl.Point3Array(ls_pt), pgl.Index3Array(ls_ind))