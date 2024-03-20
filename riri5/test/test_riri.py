from riri5 import RIRI5 as riri
import numpy as np

#############
#test 1 : calcul  avec tous les voxels sur matrice heterogene horizontalement

m = np.array([[0.,1.,3.,3.,1.,0.],[0.,1.,3.,3.,1.,0.], [0.,1.,3.,3.,1.,0.], [0.,1.,3.,3.,1.,0.]]) #32 repartis en 6*4 heterogene
ma = np.array([m, m,m, m])/4. #lai esp1 = 32/24
ma2 = np.array([m, m,m, m])/8. #lai esp2 = 16/24
ls_mlai = np.array([ma, ma2])
surf_sol = 6*4

ves = np.array([ma[:,2,2], ma[:,1,1]]) #liste de lai par especeet voxel sur le chemin du rayon
ks = np.array([1., 0.4]) #liste de k par espece
I0 = 100.

# calcul pour 1 rayon
#pour 1 rayon
res_trans1, res_abs_i1 = riri.calc_extinc_ray_multi(ves, I0, ks) # proche d'un lai de 4 : 1000.*exp(-1.*4.)


#calcul des listes de triplets selon la matrice
opt_direction = 'VXpXmYpYm'#'V'
triplets = riri.get_ls_triplets(ma, opt=opt_direction)


#distributtions d'angles
incli1, incli2 = 45., 45.
ls_distf = [riri.disttetaf(abs(incli1), 0.), riri.disttetaf(abs(incli2), 0.)]


# calul pour tous les rayons
res_trans, res_abs_i = riri.calc_extinc_allray_multi(ls_mlai, triplets ,ls_distf , I0)


print("exemple 1 : 2 especes heterogene horizontalement grille3D, differents lai, meme incli")
print('mLAI1:', ma)
print('mLAI2:', ma2)
print('LAI:', [np.sum(ma)/surf_sol,  np.sum(ma2)/surf_sol])# en suppsant cellules de 1 m2
print('incli (degre):', [incli1, incli2])
print('opt:',opt_direction)
print('I0:',I0)

print('res_abs par espece:', [np.sum(res_abs_i[0])/surf_sol, np.sum(res_abs_i[1])/surf_sol])
print('transmis sol:', np.sum(res_trans[-1])/surf_sol)
print("!!! dx=dy=2dz doit etre respecte (pas renseigne) !!!")
print("")


##############
#test 2: calcul  avec tous les voxels sur matrice creuse homogene horizontalement: calc_extinc_allray_multi

ma3 = np.zeros([8,4,4])
ma3[6:8,:,:] = 1.#1./(4*4*2)
ma4 = np.zeros([8,4,4])
ma4[5:8,:,:] = 0.5#0.5/(4*4*3)
ls_mlai2 = np.array([ma3, ma4])
surf_sol = 4*4

triplets = riri.get_ls_triplets(ma3, opt='VXpXmYpYm')
ls_distf = [riri.disttetaf(abs(45.), 0.), riri.disttetaf(abs(45.), 0.)]
res_trans, res_abs_i = riri.calc_extinc_allray_multi(ls_mlai2, triplets ,ls_distf , I0)

print("exemple 2 : 2 especes homogene horizontalement grille3D, differents lai, differentes hauteur, meme incli")
print('mLAI1:', ma3)
print('mLAI2:', ma4)
print('LAI:', [np.sum(ma3)/surf_sol,  np.sum(ma4)/surf_sol])# en suppsant cellules de 1 m2
print('incli (degre):', [incli1, incli2])
print('opt:',opt_direction)
print('I0:',I0)

print('res_abs par espece:', [np.sum(res_abs_i[0])/surf_sol, np.sum(res_abs_i[1])/surf_sol])
print('transmis sol:', np.sum(res_trans[-1])/surf_sol)
print("!!! dx=dy=2dz doit etre respecte (pas renseigne) !!!")
print("")

#rq: calcul fait pour un hmax (hauteur voxel) suppose respecte les hypothese de riri5: 


############
#test 3: calcul sur matrice creuse avec avec fonction qui redimensionne (plus rapide!): calc_extinc_allray_multi_reduced

res_trans, res_abs_i = riri.calc_extinc_allray_multi_reduced(ls_mlai2, triplets ,ls_distf , I0, optsky=None, opt=opt_direction)


