
from riri5 import RIRI5 as riri
import numpy as np


############
# exemple 1 : 3 especes homogene horizontalement, differents LAI, meme incli, meme hauteurs
############

######
# definir info plante et grille

Hmax1 = 100. #cm
Hmax2 = 100. #80.
Hmax3 = 100. #60.

lai1 = 1.5 #m2/m2
lai2 = 1.
lai3 = 0.5

incli1 = 45. #degre
incli2 = 45.
incli3 = 45.

#grille matrice homogene 1D
HmaxMat = 110.
nz = 11
base = np.array([[1.]])
mat1 = np.array([base]*nz)
I0=100.

######
# definir liste de distribution de lai

m_lai1 = mat1*0.
m_lai2 = mat1*0.
m_lai3 = mat1*0.
#laisse premiere couche vide
m_lai1[1:11,0,0] = lai1/(nz-1)
m_lai2[1:11,0,0] = lai2/(nz-1)
m_lai3[1:11,0,0] = lai3/(nz-1)


ls_mlai = np.array([m_lai1, m_lai2, m_lai3])


# distributtions incli en classes d'angles
sd = 0. #ecart type distribution
dist1 = riri.disttetaf(abs(incli1), sd)
dist2 = riri.disttetaf(abs(incli2), sd)
dist3 = riri.disttetaf(abs(incli3), sd)

ls_distf = [dist1, dist2, dist3]


# coord voxels passages rayons
opt_direction = 'V' #'VXpXmYpYm'#
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)
#triplets = riri.get_ls_triplets(m_lai1, opt='V')
#triplets = riri.get_ls_triplets(m_lai1, opt='VXpXmYpYm')


#####
# calcul

#calcul avec reduction automatique: plante si rien a reduire!! prevu pour avoir une marge sur les hauteurs
res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0, opt=opt_direction)

# sinon: calcul sans reduction
#res_trans, res_abs_i = riri.calc_extinc_allray_multi(ls_mlai, triplets ,ls_distf , I0=100) #calcul sans reduction

print("exemple 1 : 3 especes homogene horizontalement, differents lai, meme incli, meme hauteurs")
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
print("")






############
# exemple 2 : 3 especes homogene horizontalement, meme LAI, differentes incli, meme hauteurs
############

######
# definir info plante et grille

Hmax1 = 100. #cm
Hmax2 = 100. #80.
Hmax3 = 100. #60.

lai1 = 1. #m2/m2
lai2 = 1.
lai3 = 1.

incli1 = 45. #degre
incli2 = 20.
incli3 = 70.

#grille matrice homogene 1D
HmaxMat = 110.
nz = 11
base = np.array([[1.]])
mat1 = np.array([base]*nz)
I0=100.

######
# definir une liste de distribution de lai

m_lai1 = mat1*0.
m_lai2 = mat1*0.
m_lai3 = mat1*0.
m_lai1[1:11,0,0] = lai1/(nz-1)
m_lai2[1:11,0,0] = lai2/(nz-1)
m_lai3[1:11,0,0] = lai3/(nz-1)


ls_mlai = np.array([m_lai1, m_lai2, m_lai3])


# distributtions incli en classes d'angles
sd = 0. #ecart type distribution
dist1 = riri.disttetaf(abs(incli1), sd)
dist2 = riri.disttetaf(abs(incli2), sd)
dist3 = riri.disttetaf(abs(incli3), sd)

ls_distf = [dist1, dist2, dist3]


# coord voxels passages rayons
opt_direction = 'V' #'VXpXmYpYm'
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)



#calcul sans reduction
#res_trans, res_abs_i = riri.calc_extinc_allray_multi(ls_mlai, triplets ,ls_distf , I0=100) #calcul sans reduction
#calcul avec reduction
res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0, opt=opt_direction)
# plante si rien a reduire!! prevu pour avoir une marge sur les hauteurs


print("exemple2: 3 especes homogene horizontalement, meme LAI, differentes incli, meme hauteurs")
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
print("")

#relance avec autres directions
opt_direction = 'VXpXmYpYm' #'V'
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)


res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0=100, opt=opt_direction)


print("exemple2: 3 especes homogene horizontalement, meme LAI, differentes incli, meme hauteurs")
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
print('!!! dx=dy=2dz pas respecte !!')
print("")







############
# exemple 3 : 3 especes homogene horizontalement, meme LAI, meme incli, differentes hauteurs
############

######
# definir info plante et grille

Hmax1 = 100. #cm
Hmax2 = 80.
Hmax3 = 60.

lai1 = 1. #m2/m2
lai2 = 1.
lai3 = 1.

incli1 = 45. #degre
incli2 = 45.
incli3 = 45.

#grille matrice homogene 1D
HmaxMat = 110.
nz = 11
base = np.array([[1.]])
mat1 = np.array([base]*nz)


######
# definir une liste de distribution de lai

nbcouches1 = int(nz*Hmax1/HmaxMat)
nbcouches2 = int(nz*Hmax2/HmaxMat)
nbcouches3 = int(nz*Hmax3/HmaxMat)

m_lai1 = mat1*0.
m_lai2 = mat1*0.
m_lai3 = mat1*0.

m_lai1[(11-nbcouches1):11,0,0] = lai1/nbcouches1
m_lai2[(11-nbcouches2):11,0,0] = lai2/nbcouches2
m_lai3[(11-nbcouches3):11,0,0] = lai3/nbcouches3


ls_mlai = np.array([m_lai1, m_lai2, m_lai3])


# distributtions incli en classes d'angles
sd = 0. #ecart type distribution
dist1 = riri.disttetaf(abs(incli1), sd)
dist2 = riri.disttetaf(abs(incli2), sd)
dist3 = riri.disttetaf(abs(incli3), sd)

ls_distf = [dist1, dist2, dist3]


# coord voxels passages rayons
opt_direction = 'V' #'VXpXmYpYm'#
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)
#!!! presuppose dx=dy=dz / tan(0.4637) = 2*dz pour traversee de biais!
#ps respecte ici


res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0=100, opt=opt_direction)


print("exemple 3: 3 especes homogene horizontalement, meme LAI, meme incli, differentes hauteurs")
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
print("")




############
# exemple 4 : 3 especes homogene horizontalement, meme LAI, meme incli, differentes hauteurs, grilles 3D
############

######
# definir info plante et grille

Hmax1 = 100. #cm
Hmax2 = 80.
Hmax3 = 60.

lai1 = 1. #m2/m2
lai2 = 1.
lai3 = 1.

incli1 = 45. #degre
incli2 = 45.
incli3 = 45.

#grille matrice homogene 1D
HmaxMat = 110.
cote = 100.
nx, ny = 5, 5
nz = 11
dx, dy = cote/nx, cote/ny
dz = HmaxMat/nz
base = np.array([[1.]*nx]*ny)
mat1 = np.array([base]*nz)
I0=100.

######
# definir une liste de distribution de lai

nbcouches1 = int(nz*Hmax1/HmaxMat)
nbcouches2 = int(nz*Hmax2/HmaxMat)
nbcouches3 = int(nz*Hmax3/HmaxMat)

m_lai1 = mat1*0.
m_lai2 = mat1*0.
m_lai3 = mat1*0.

m_lai1[(11-nbcouches1):11,:,:] = lai1/nbcouches1 #lai1/(nbcouches1*nx*ny)
m_lai2[(11-nbcouches2):11,:,:] = lai2/nbcouches2 #lai2/(nbcouches2*nx*ny)
m_lai3[(11-nbcouches3):11,:,:] = lai3/nbcouches3 #lai3/(nbcouches3*nx*ny)
#! chaque voxel doit avoir une valeur en LAI partiel vertical!, mais pas horizontalement (somme en colonne vertical fait LAI local)

ls_mlai = np.array([m_lai1, m_lai2, m_lai3])


# distributtions incli en classes d'angles
sd = 0. #ecart type distribution
dist1 = riri.disttetaf(abs(incli1), sd)
dist2 = riri.disttetaf(abs(incli2), sd)
dist3 = riri.disttetaf(abs(incli3), sd)

ls_distf = [dist1, dist2, dist3]


# coord voxels passages rayons
opt_direction = 'VXpXmYpYm'#'V' #
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)
#!!! presuppose dx=dy=dz / tan(0.4637) = 2*dz pour traversee de biais!
#ps respecte ici


res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0/(nx*ny), opt=opt_direction)
#valeur de I0 est pour entree en un voxel, adapter a surface ou nb de voxels!!


print("exemple 4: 3 especes homogene horizontalement, meme LAI, meme incli, differentes hauteurs, grille 3D bonne dimension")
#print("ls_mlai", ls_mlai)
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
if dz == 0.5*dx and dz==0.5*dy:
    print("dx=dy=2dz bien respecte")
print("")




## avec autres option rayons vertical

# coord voxels passages rayons
opt_direction = 'V' #'VXpXmYpYm'#
triplets = riri.get_ls_triplets(m_lai1, opt=opt_direction)
#!!! presuppose dx=dy=dz / tan(0.4637) = 2*dz pour traversee de biais!
#ps respecte ici


res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai, triplets ,ls_distf , I0/(nx*ny), opt=opt_direction)
#valeur de I0 est pour entree en un voxel, adapter a surface ou nb de voxels!!


print("exemple 4: 3 especes homogene horizontalement, meme LAI, meme incli, differentes hauteurs, grille 3D bonne dimension")
#print("ls_mlai", ls_mlai)
print('Hmax (cm):', [Hmax1, Hmax2, Hmax3])
print('LAI:', [lai1, lai2, lai3])
print('incli (degre):', [incli1, incli2, incli3])
print('opt:',opt_direction)
print('I0:',I0)
#print('res_trans', res_trans)
#print('res_abs1', res_abs_i[0])
#print('res_abs2', res_abs_i[1])
#print('res_abs3', res_abs_i[2])
print('res_abs par espece:', [np.sum(res_abs_i[0]), np.sum(res_abs_i[1]), np.sum(res_abs_i[2])])
print('transmis sol:', np.sum(res_trans[-1]))
if dz == 0.5*dx and dz==0.5*dy:
    print("dx=dy=2dz bien respecte")
print("")
