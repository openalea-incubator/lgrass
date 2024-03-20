'''

 RIRI5: A 3D Turbid medium model adapted from RATP model- simplified for a  5 directions sky turtle 
 ******************************
 Authors: G. Louarn / D. Combes - Mai 2014
 
 - Considers diffuse light only (soc / uoc)
 - Considers beam elevations of a minimal hemispheric turtle with 6 faces: 1 vertical direction + 1 elevation ring
 - Considers only beams crossing by voxel diagonals into (x,y) and (y,z) planes : 1 vertical beam + 4 for diagonal beams = 5 directions sampled
 - This imposes a voxel size and a 3D drid respecting : dx = dy = dz/tan(0.4637) = 2*dz
 

'''



from numpy import linspace
from scipy import array, exp, zeros, ones, set_printoptions, pi, radians, cos, sin, tan, arccos, histogram
from copy import deepcopy
#from numpy.random import seed, normal
import numpy as np


def def_na_lims(pattern8, dz_ini, Hmax, opt='3D', unit='cm'):
    """ Calculate the number of voxels per direction and the limits of the voxels in a 3D grid respecting the ratio asumption of RIRI5, given a pattern8, dz_ini and Hmax 
    
    :param pattern8: A list of two [x,y] coordinates setting the upper corner and lower corner of the scene/grid
    :type pattern8: list
    :param dz_ini: targeted dz of the grid, which is then adpated if it does not fit with the ratio asumption 
    :type dz_ini: float
    :param Hmax: targeted Hmax of the grid, which is then adpated if it does not fit with the ratio asumption
    :type Hmax: float
    :param opt: either '3D' for a 3D grid or '1D' for a vertical one dimention grid
    :type opt: string
    :param unit: unit used for input lengths : either 'cm' or 'm'
    :type unit: string
    
    :return:
        * "number of voxels" : [nx, ny, nz]
        * "voxel size" : [dx, dy, dz]
        * "3 lists of voxel limits" : [xlims,ylims,zlims]
        * "grid origin" : [xorigin, yorigin, zorigin]
        * "reference surface of a voxel (m2)" 
    
    """
    if opt=='1D':
        nxa,nya = 1,1
        nza = int(Hmax/dz_ini)+1
        dxa,dya = pattern8[1][0] - pattern8[0][0], pattern8[1][1] - pattern8[0][1]
        dza = dz_ini
    elif opt=='3D':
        nxa = int((pattern8[1][0] - pattern8[0][0])/(2*dz_ini))
        nya = int((pattern8[1][1] - pattern8[0][1])/(2*dz_ini))
        dxa = (pattern8[1][0] - pattern8[0][0])/float(nxa)
        dya = (pattern8[1][1] - pattern8[0][1])/float(nya)
        dza = 0.5*dxa #dz recalcule pour 
        nza = int(Hmax/dza)+1

    xlims = linspace(pattern8[0][0], pattern8[1][0],nxa+1)
    ylims = linspace(pattern8[0][1], pattern8[1][1],nya+1)
    zlims = array([x * dza for x in range(0, nza+1)])

    na = [nxa,nya,nza]
    lims = [xlims,ylims,zlims]
    dxyz = [dxa,dya,dza]
    origin_grid = array([pattern8[0][0], pattern8[0][1],zlims[-1]])
    
    if unit=='cm':
        surf_refVOX = dxa*dya/10000.#m2
    else:
        surf_refVOX = dxa*dya #suppose m
    
    return na, dxyz, lims, origin_grid, surf_refVOX

#na, dxyz, lims, origin_grid = def_na_lims(pattern8, dz_aerien, Hmax)


## dans quel voxel?
def WhichVoxel(p, origin_grid, na, dxyz):
    """ Determine in which voxel of a 3D grid is located a point 'p'
    
    :param p: point [x,y,z] to be located
    :type p: array
    :param origin_grid: grid origin [xorigin, yorigin, zorigin]
    :type origin_grid: array
    :param na: number of voxels per direction : [nx, ny, nz]
    :type na: list
    :param dxyz: voxel size : [dx, dy, dz]
    :type dxyz: list
    
    :return: A list of integer corresponding to the voxel ids [x,y,z] in the 3D grid (for Z, id=0 corresponds to the top of the grid/canopy)
   
    
    """
    # en z, id=0 = haut du couvert 
    
    p1_rel = p - origin_grid
    vox = [int(p1_rel[0]//dxyz[0]), int(p1_rel[1]//dxyz[1]), int(-p1_rel[2]//dxyz[2])] #// division entiere

    #test si dans grille et sinon retouve indice correspondant
    test = [0 <= vox[0] <na[0], 0 <= vox[1] <na[1], 0 <= vox[2] <na[2]]
    if test != [True,True,True]:#si hors grille
        #recup les indices du fautif
        matches = [i for i in range(3) if test[i]==False]
        #faire le calcul du bon id equivalent
        for i in matches:
            if i!=2: #x,y -> couvert infini
                vox[i] = vox[i]%na[i] # marche pour > et <0 et na[i]=1
            else: #z ->  
                if vox[i]<0:#au dessus: met dans la derniere strate
                    vox[i] = 0
                else: #en dessous: met dans la strate au dessus du sol
                    vox[i] = na[i]-1
    return vox


#na
#WhichVoxel(array([-12.,-10.,130.]), origin_grid, na, dxyz)
#WhichVoxel(array([-200.,100.,230.]), origin_grid, na, dxyz)#ok
#WhichVoxel(array([-200.,10.,10.]), origin_grid, na, dxyz)#ok
#WhichVoxel(array([-200.,10.,208.]), origin_grid, na, dxyz)#ok



def get_tripletY(ydeb, x, ny, nz, sens='+'):
    """ List of triplet voxel IDs for beams in Y plane (y,z)
    
    :param ydeb: initial beam position, ID along y
    :type ydeb: integer
    :param x: x voxel ID
    :type x: integer
    :param ny: number of voxels along y in the plane
    :type ny: integer
    :param nz: number of voxels along z in the plane
    :type nz: integer
    :param sens: Direction of the diagonal beams with respect to vertical; either '+' or '-' 
    :type sens: string
    
    :return: A list of 3 lists with voxel indices [zz, yy, xx] of the beams path
    
    """
    #""" pour un plan dans y, rayon allant dans le + ou le - """
    
    zz = list(range(nz))
    xx = [x]*nz
    res = []
    count = ydeb
    if sens=='+':
        for k in range(nz):
            if count > ny-1:
                res.append(0)
                count=1
            else:
                res.append(count)
                count = count+1

    elif sens=='-':
        for k in range(nz):
            if count < 0:
                res.append(ny-1)
                count=ny-2
            else:
                res.append(count)
                count=count-1

    return [zz, res, xx]

#get_tripletY(3, 0, 4, 4, '-')
#get_tripletY(3, 0, 4, 4, '+')
#get_tripletY(0, 0, 4, 4, '-')
#get_tripletY(0, 0, 4, 4, '+')
#get_tripletYbis(0, 0, 4, 10, sens='+')

#utilise un count (initialiser avec y deb) et meme demarche
def get_tripletX(xdeb, y, nx, nz, sens='+'):
    """ List of triplet voxel IDs for beams in X plane (x,z)
    
    :param xdeb: initial beam position, ID along x
    :type xdeb: integer
    :param y: y voxel ID
    :type y: integer
    :param nx: number of voxels along x in the plane
    :type nx: integer
    :param nz: number of voxels along z in the plane
    :type nz: integer
    :param sens: Direction of the diagonal beams with respect to vertical; either '+' or '-' 
    :type sens: string
    
    :return: A list of 3 lists with voxel indices [zz, yy, xx] of the beams path
    
    
    """
    #""" pour un plan dans y, rayon allant dans le + ou le - """
    
    zz = list(range(nz))
    yy = [y]*nz
    res = []
    count = xdeb
    if sens=='+':
        for k in range(nz):
            if count > nx-1:
                res.append(0)
                count=1
            else:
                res.append(count)
                count=count+1
    elif sens=='-':
        for k in range(nz):
            if count < 0:
                res.append(nx-1)
                count=nx-2
            else:
                res.append(count)
                count=count-1

    return [zz, yy, res]

#get_tripletX(3, 0, 6, 4, '+')
#get_tripletX(1, 0, 6, 4, '+')
#get_tripletX(1, 0, 6, 4, '-')
#get_tripletX(0, 2, 4, 10, sens='-')

def get_tripletVert(x,y,nz):
    """ List of triplet voxel IDs for vertical beams
    
    :param x: x voxel ID
    :type x: integer
    :param y: y voxel ID
    :type y: integer
    :param nz: number of voxels along z in the plane
    :type nz: integer
    
    :return: A list of 3 lists with voxel indices [zz, yy, xx] of the beams path
    
    """
    #""" pour triplets verticaux """
    
    xx = [x]*nz
    yy = [y]*nz
    zz = list(range(nz))
    return [zz,yy,xx]



def get_ls_triplets(mat, opt='VXpXmYpYm'):
    """ List of triplet voxel IDs for all beams in X plane, Y plane and vertical direction
    
    :param mat: a 3D grid of format
    :type mat: array
    :param opt: A string to indicate which beam direction to be considered; either 'VXpXmYpYm' (5 directions) or 'V' (vertical beam only)
    :type opt: string
    
    :return: A list of 3 lists with voxel indices [zz, yy, xx] of the beams path
    
    """
    #""" liste de triplets par direction; V= vertical; Xp,Yp,Xm,Ym 4 rayons a 45 degre dans le plan normal aux faces des voxels"""
    
    nz,ny,nx = mat.shape
    ls_tripletYp, ls_tripletYm =[], []
    ls_tripletXp, ls_tripletXm =[], []
    ls_tripletVert = []
    for i in range(nx):
        for j in range(ny):
            ls_tripletYp.append( get_tripletY(j, i, ny, nz, '+') )
            ls_tripletYm.append( get_tripletY(j, i, ny, nz, '-') )
            ls_tripletXp.append( get_tripletX(i, j, nx, nz, '+') )
            ls_tripletXm.append( get_tripletX(i, j, nx, nz, '-') )
            ls_tripletVert.append(get_tripletVert(i,j,nz))
    
    ls_triplet = []
    if opt.find('V')>=0:
        ls_triplet.append(ls_tripletVert)
    if opt.find('Xp')>=0:
        ls_triplet.append(ls_tripletXp)
    if opt.find('Xm')>=0:
        ls_triplet.append(ls_tripletXm)
    if opt.find('Yp')>=0:
        ls_triplet.append(ls_tripletYp)
    if opt.find('Ym')>=0:
        ls_triplet.append(ls_tripletYm)
    
    return ls_triplet

#triplets = get_ls_triplets(ma, opt='V')
#triplets = get_ls_triplets(ma, opt='VXpXmYpYm')


def k_teta_DC(mean_incl, elevations=[9.23,10.81,26.57,31.08,47.41,52.62,69.16,90]):
    """ Computes extinction coefficient from mean elevation angle according to SIRASCA model (p210-211, "Crop stucture and light microclimate" book)
    
        :param mean_incl: mean elevation angle of leaves (degree)
        :type mean_incl: float
        :param elevations: list of beam elevation classes to consider for integration (degres); default=elevation for turtle with 46 directions
        :type elevations: list
        
        :return: A list of ...
   
        
    """
    
    #computes extinction coefficient from mean inclination
    # equations from SIRASCA. Calculations for each elevation angle (p210-211 crop stucture and light microclimate)
    # default elevation for turtle 46 directions 
    #elevations (beam, degre)-> 90 = vertical ; 0=horizontal
    #mean_incl (feuille, degre) -> 90 = vertical ; 0=horizontal
    
    nh=len(elevations)
    xinc = radians(mean_incl)
    lst_xk=[]
    for ih in elevations:
        hh = radians(ih)
        if (hh >=xinc):
            xk = cos(xinc)
        else:
            xmm = -tan(hh)/tan(xinc)
            #print xmm
            xmm = arccos(xmm)
            xk = cos(xinc)*(2*(xmm - tan(xmm))/pi - 1.)

        lst_xk.append(xk)
    return lst_xk
  

#k_teta_DC(45)
#k_teta_DC(25., elevations=[5., 15., 25., 35., 45., 55., 65., 75., 85.])



def k_teta_distf(teta_beam, dist_class_teta_f):
    """ Computes extinction coefficient for a given beam integrating the ditribution of leaf area of a species into leaf elevation classes
    
        :param teta_beam: Beam elevation with respect to horizontal plane (degree)
        :type teta_beam: float
        :param dist_class_teta_f: relative distribution of leaf area in elevation classes ; 9 classes in degrees, with centers [5,15,25,35,45,55,65,75,85]
        :type dist_class_teta_f: list
        
        :return: Calculated extinction coefficient
        
    """
    #""" pour tenir compte d'une distribution d'incli quelconque """
    
    # teta_beam, ang_class_f (en degre, par rapport au plan horizontal)
    #dist_class_teta_f: relative distribution des surface par classe d'incli, tous les 10 degre, centre [5,15,25,35,45,55,65,75,85]
    
    ang_class_f = [5,15,25,35,45,55,65,75,85] 
    res = []
    for i in range(len(ang_class_f)):
        kf = k_teta_DC(float(ang_class_f[i]), [teta_beam])[0]*dist_class_teta_f[i]
        res.append(kf)

    return sum(res)



def disttetaf(mf, sdf, nbs=10000, seed=0):
    """ Distribute leaf angle normal ditribution into discrete regular elevation classes - A python version of 'dist_tetaf.r'
    
        :param mf: mean elevation (degre)
        :type mf: float
        :param sdf: stadard deviation of the normal distribution (degree)
        :type sdf: float
        :param nbs: number of random drawings to generate de districution
        :type nbs: integer
        :param seed: seed number for the random number generator
        :type seed: integer
        
        :return: A list with relative leaf area distribution in 9 elevation angle classes : [5,15,25,35,45,55,65,75,85]
        
    """
    
    #""" python version of dist_tetaf.r"""
    
    #mf = 5.
    #sdf = 10.
    #n = 100#00
    #seed = 0
    
    #np.random.seed(seed)
    x = np.random.normal(loc=mf, scale=sdf, size=nbs) #numpy.random.normal
    x
    x[x<0] = -x[x<0]
    x[x>90] = 90 - (x[x>90] - 90)
    x[x < 0] = -x[x < 0]
    x[x > 90] = 90 - (x[x > 90] - 90)
    x
    res = histogram(x, bins=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    return res[0] / float(nbs)
    #print 't', disttetaf(5., 10.)



def calc_extinc_ray_multi(ls_vlai, I0, ls_k):
    """ Calculate light extinction and absoption per species for a beam going through a series of voxels
    
        :param ls_vlai: List of partial LAI (Leaf Are Index, unit: m2.m-2) per species for the series of voxels crossed by the beam
        :type ls_vlai: array
        :param I0: Incoming light from the beam angular sector (unit: )
        :type I0: float
        :param ls_k: List of extinction coefficients per species (unitless)
        :type ls_k: array
        
        :return:
            * "res_trans: " : array
            * "res_abs_i: " : list of arrays
        
    """
    
    #""" calcul extinction (profil de Iout) et absorption par espece pour un rayon traversant une serie de voxel en serie decrite par profils de LAI (du haut en bas) et k """
    #""" pour partage entre especes dans un voxel homogene  - sinoquet et al. 2000 """
    #""" prend un k par entite/espece, mais peut facilement prendre profils de k """

    I = I0
    res_trans = []
    res_abs_i = []
    for strate in range(len(ls_vlai[0])):
        lai_ki = ls_vlai[:,strate] * ls_k
        sum_lai_ki = sum(lai_ki)
        if sum_lai_ki!=0. :
            LIEi = lai_ki/sum_lai_ki * (1 - exp(-sum_lai_ki))
        else: #si tout les lai sont a zero, renvoir zero
            LIEi = lai_ki

        Iout = I*(1-sum(LIEi))

        res_trans.append (Iout)
        res_abs_i.append(I*LIEi)
        I=Iout

    res_abs_i = array(res_abs_i).T #chaque ligne correspond a une entite/espece dans l'odre des ls_vlai
    return res_trans,  res_abs_i



#m = array([[0.,1.,3.,3.,1.,0.],[0.,1.,3.,3.,1.,0.], [0.,1.,3.,3.,1.,0.], [0.,1.,3.,3.,1.,0.]])
#ma = array([m, m,m, m])/4.
#ma2 = array([m, m,m, m])/8.

#ves = array([ma[:,2,2], ma[:,1,1]]) #liste de lai par especeet voxel sur le chemin du rayon
#ks = array([1., 0.4]) #liste de k par espece
#I0 = 1000.
#res_trans, res_abs_i = calc_extinc_ray_multi(ves, I0, ks) # proche d'un lai de 4 : 1000.*exp(-1.*4.)


def calc_extinc_allray_multi(ls_mlai, ls_triplets_dir ,ls_distf , I0, optsky=None):
    """ Calculate light extinction and absoption per species for a series of beams 
    
        :param ls_mlai: List of 3D grids for LAI distribution (Leaf Are Index, unit: m2.m-2) per species
        :type ls_mlai: list
        :param ls_triplets_dir: List of list of triplet voxel indices [zz, yy, xx] for the different beam directions considered
        :type ls_triplets_dir: list of 
        :param ls_distf: List of dist_class_teta_f (relative distribution of leaf area in elevation classes) by species
        :type ls_distf: list
        :param I0: Incoming light (unit: )
        :type I0: float
        :param optsky: Option for diffuse sky, relative light distribution with elevation angle; either 'soc' or 'uoc'
        :type optsky: string
        
        :return:
            * "res_trans: " : array
            * "res_abs_i: " : list of arrays
        
    """
    #""" calcul extinction (profil de Iout) et absorption par espece pour des liste de voxel, regroupe par direction """
    #""" presupose que ls_triplets tout ou partie d'un turtle 6 => dx=dy=dz / tan(0.4637) = 2*dz !"""
    ##res = deepcopy(ls_triplets_dir) #liste des triplets de voxel par direction, associe a res_trans et res_abs_i

    alfa_turtle6 = 0.4637 #radians / 26.57 degree
    # distribution de I0 entre sources
    n_dir = len(ls_triplets_dir)
    if optsky==None or n_dir<5: #cas par defaut = toutes directions le meme poids / a utiliser notamment qd veut tester une seule dir
        ls_poids = array([1./float(n_dir)]*n_dir)  
    elif n_dir == 5 and optsky=='uoc': #diffus couvert / poids pour un turtle 6 ramene a 5 direction
        ls_poids = array([1./6.]+[5./(6.*4)]*4)
    elif n_dir == 5 and optsky=='soc': #diffus clair
        ls_poids = [1./6.]+[5./(6.*4)]*4
        effet_sin = [sin(pi/2)]+ [sin(alfa_turtle6)]*4
        x = array(ls_poids)*effet_sin
        ls_poids = x/sum(x)

    #print 'ls poids:', ls_poids

    #coeff de ponderation des lai par longueur de parcours dans le voxel
    # -> juste un effet dz (du coup dans mon cas tout dz est parcourru)(effet longueur proprement dit integre dans le k)

    #coeff d'extinction par direction et espece kteta
    ls_k_teta = []
    for i in range(len(ls_distf)):
        ls_k_teta.append([k_teta_distf(90., ls_distf[i]), k_teta_distf(26.57, ls_distf[i])])

    ls_k = []
    for i in range(len(ls_distf)):
        if ls_triplets_dir[0][0][1] == ls_triplets_dir[0][0][2]: #y a 1 ray vertical, place en premier selon get_ls_triplets
            ls_k.append([ls_k_teta[i][0]] + [ls_k_teta[i][1]]*(n_dir-1))
        else:
            ls_k.append([ls_k_teta[i][1]]*n_dir)

    ls_k = array(ls_k)
    #print 'ls k:', ls_k

    #calculs
    res_trans_form = zeros(ls_mlai[0].shape)#sortie pour recuperer les cumuls trans
    res_abs_i_form  = zeros(ls_mlai.shape)#sortie pour recuperer les cumul abs par espece
    for dir in range(len(ls_triplets_dir)):
        for tri in range(len(ls_triplets_dir[dir])):
            ves = ls_mlai[:,ls_triplets_dir[dir][tri][0], ls_triplets_dir[dir][tri][1], ls_triplets_dir[dir][tri][2]]
            res_trans, res_abs_i = calc_extinc_ray_multi(ves, I0*ls_poids[dir], ls_k[:,dir])#[1.,0.6]) #pour korcer des k 
            
            res_trans_form[ls_triplets_dir[dir][tri][0], ls_triplets_dir[dir][tri][1], ls_triplets_dir[dir][tri][2]] += res_trans
            res_abs_i_form[:, ls_triplets_dir[dir][tri][0], ls_triplets_dir[dir][tri][1], ls_triplets_dir[dir][tri][2]] += res_abs_i
    
    return res_trans_form, res_abs_i_form

#sortir distribution de I0/ls_poids et ls_k (a passer en argument) pour une fonction fortan uniquement calcul


def calc_extinc_allray_multi_reduced(ls_mlai, ls_triplets_dir, ls_distf, I0, optsky=None, opt='VXpXmYpYm'):
    """ Calculate light extinction and absoption per species for a series of beams, avoiding calculation in empty voxels above the canopy
    
        :param ls_mlai: List of 3D grids for LAI distribution (Leaf Are Index, unit: m2.m-2) per species
        :type ls_mlai: list
        :param ls_triplets_dir: List of list of triplet voxel indices [zz, yy, xx] for the different beam directions considered
        :type ls_triplets_dir: list of 
        :param ls_distf: List of dist_class_teta_f (relative distribution of leaf area in elevation classes) by species
        :type ls_distf: list
        :param I0: Incoming light (unit: )
        :type I0: float
        :param optsky: Option for diffuse sky, relative light distribution with elevation angle; either 'soc' or 'uoc'
        :type optsky: string
        :param opt: A string to indicate which beam direction to be considered; either 'VXpXmYpYm' (5 directions) or 'V' (vertical beam only)
        :type opt: string
        
        :return:
            * "res_trans: " : array
            * "res_abs_i: " : list of arrays
    
    """

    # combien/quelles lignes a zeros de LAI au dessus
    laicum = np.sum(ls_mlai, axis=0)
    laicumvert = np.sum(laicum, axis=(1, 2))
    nb0 = 0  # nb de couches sans feuilles/LAI
    for i in range(len(laicumvert)):
        if laicumvert[i] == 0.:
            nb0 += 1
        else:
            break

    # ajouter un if sur nb0 ou le faire a chaque fois?
    # redim des m_lai et triplets
    shp = np.shape(ls_mlai)
    reduced_mlai = []
    for plt in range(shp[0]):
        reduced_mlai.append(ls_mlai[plt, (nb0 - 1):shp[1], :, :])

    reduced_triplets = get_ls_triplets(reduced_mlai[0], opt)
    reduced_mlai = array(reduced_mlai)

    # ls_distf = [riri.disttetaf(abs(45.), 0.), riri.disttetaf(abs(45.), 0.)]
    res_trans_form_red, res_abs_i_form_red = calc_extinc_allray_multi(reduced_mlai, reduced_triplets, ls_distf, I0, optsky)

    # expand matrices de sortie de nb0-1 couches sans faire les calculs
    shpnew = np.shape(res_abs_i_form_red)
    if nb0<=0:
        print("shoot canopy too high : m_lai out of voxels!")
    
    # zz = zeros((shpnew[0], (nb0-1), shpnew[2], shpnew[3]))
    # oo = ones(((nb0-1), shpnew[2], shpnew[3]))*res_trans_form_red[0,0,0]
    res_trans_form = np.concatenate((ones(((nb0 - 1), shpnew[2], shpnew[3])) * res_trans_form_red[0, 0, 0], res_trans_form_red), axis=0)
    res_abs_i_form = np.concatenate((zeros((shpnew[0], (nb0 - 1), shpnew[2], shpnew[3])), res_abs_i_form_red), axis=1)

    return res_trans_form, res_abs_i_form


#### R:FR L Faverjon

def schnute(x,a,b,c,d,x1,x2):
    """ General Schnute function
    
    """
    #"""fonction schnute generale appliquee a x"""
    Y=pow((pow(c,b)+(pow(d,b)-pow(c,b)))*((1-exp(-a*(x-x1)))/(1-exp(-a*((x2-x1))))),(1/b))
    return Y

def rfr_calc_relatif(relI0,a=3.09,b=1.59,c=0,d=1.12,x1=0.,x2=2.):
    """ Calculate Red:Far red ratio from relative transmitted PAR - function and parameter values according to Escobar et al (2009) 
        (Agricultural and Forest Meteorology, 149(8), 1244-1253)
    
        :param relI0: A 3D grid for relative transmitted PAR radiation
        :type relI0: array
    
        :return:
            * "res_rfr: " : array
    """
    
    #"""calcul du ratio rouge clair:rouge sombre narrowband a partir du PAR transmis relatif par rapport au rayonnement incident (relI0)
    #les parametres par defaut a=3.09,b=1.59,c=0,d=1.12,x1=0.,x2=2. de la fonction schnute sont ajustes sur un couvert de sorgho, 
    #Escobar et al 2009 - Agricultural and Forest Meteorology, 149(8), 1244-1253  """
    res_rfr = schnute(relI0,a,b,c,d,x1,x2)
    return res_rfr

#applique a la matrice des transmis res_trans_form
#res_rfr = rfr_calc_relatif(res_trans_form)


