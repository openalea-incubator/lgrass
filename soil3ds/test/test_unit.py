
import os
import sys
import soil3ds
path_ = os.path.dirname(os.path.abspath(soil3ds.__file__))  # local absolute path of L-egume
print(('path', path_))
sys.path.insert(0, path_)

from soil3ds.soil_moduleW import * #soil3ds installe comme module

################## TESTS ###############################

def test_uni1():
    """
    
    """
    ###############
    ## DEMO 1
    ## sol 1D (S)
    ## 1 seul systeme racinaire en dynamique
    ###############

    ## meteo
    Et0 = 4.
    precips = [0.]*29+[90.]# precipitations = 20mm tous les 30 jours
    precips = precips*20
    Precip = precips[0]#0.
    Irrig = 0.

    ## sol
    par_sol = {'13':{'KST': '5', 'teta_ad': 0.0078506705623363031, 'soil number': '13', 'teta_fc': 0.35816688969184857, 'WCST': '0.503', 'teta_wp': 0.11250451238171076, 'soil type': '13 - loam', 'teta_sat': 0.503, 'gamma_theo': '0.08489778'}}

    dz_sol = 5. #cm
    ncouches_sol = 40
    soil_type = 13 #"13 - loam"
    ZESX = 0.3#m
    CFES = 1.
    Uval = 3.#U quantite d'eau dans une couche superieure en mm (5 par default)
    b = 0.63

    ##vegetation
    espi_t = [0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006]
    ls_epsi = [espi_t[0]]

    ## soil initialisation
    S = Soil(par_sol, soil_number = [soil_type]*ncouches_sol, dxyz = [[1.], [1.], [dz_sol/100.]*ncouches_sol], vDA=[1.25]*ncouches_sol, ZESX=ZESX, CFES=CFES)

    #RL_profil = table[0] #profil exporte de l-py via table
    RL_profil = RLprof_t(10, ncouches_sol)
    R1 = vert_roots(S.dxyz, RL_profil)#*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
    ls_roots = [R1]#[R1, R2]#[R3]#

    ##bouble journaliere
    ls_transp, evapo_tot, D, state,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig,previous_state=[0.,0.,0.], ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    for i in range(1,len(espi_t)):
        previous_asw = sum3(S.asw_t)
        Precip = 0.#precips[i]
        Irrig = 0.
        ls_epsi = [espi_t[i]]

        RL_profil = RLprof_t(10+i, ncouches_sol)#table[i]
        R1 = vert_roots(S.dxyz, RL_profil)#*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
        ls_roots = [R1]#[R1, R2]#[R3]#

        ls_transp, evapo_tot, D, state,  m_frac_transpi, m_frac_evap, ls_ftsw =  S.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, state, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=0.63, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)

        print(i, ls_ftsw, sum3(S.asw_t),  sum3(m_frac_evap), sum3(m_frac_transpi[0]), state)

    #    ##visu
    #    #bx = Box(Vector3(1.,1.,1.))
    #    Monviewer=Viewer
    #    MaScene=Scene()
    #    MaScene = S.plot_soil_properties (vals=S.ftsw_t, MaScene=Scene(), col_scale=5)#plot_soilWC (S.dxyz, S.m_soil_vox, S.ftsw_t, MaScene)
    #    Monviewer.display(MaScene)



def test_uni2():
    """
    
    """
    #################
    ## DEMO 2
    ## avec sol 3D (S2)
    ## 2 systemes racinaires fixes en competition (homogenes horizontalement) R3/R4
    ## en statique
    #################


    ##meteo
    Et0 = 2.
    #precips = [0.]*29+[90.]# precipitations = 20mm tous les 30 jours
    #precips = precips*20
    precips = [0.]*500#pas de precipitations (dessechement)
    Precip = precips[0]#0.
    Irrig=0.

    ## sol
    par_sol = {'13':{'KST': '5', 'teta_ad': 0.0078506705623363031, 'soil number': '13', 'teta_fc': 0.35816688969184857, 'WCST': '0.503', 'teta_wp': 0.11250451238171076, 'soil type': '13 - loam', 'teta_sat': 0.503, 'gamma_theo': '0.08489778'}}

    dz_sol = 5. #cm
    ncouches_sol = 10
    soil_type = 13 #"13 - loam"
    ZESX = 0.3 #m
    CFES=1.
    b=0.63
    Uval = 3.#sum(S.asw_t[0])#U quantite d'eau dans une couche superieure en mm (5 par default)
    stateEV = [0.,0.,0.]

    ### soil initialisation
    S2 = Soil(par_sol, soil_number = [soil_type]*ncouches_sol, dxyz = [[0.2]*5, [0.2]*5, [0.2]*10], vDA=[1.25]*ncouches_sol, ZESX=ZESX, CFES=CFES)

    ## vegetation
    ls_epsi = [0.4, 0.4]##[0.4]# !! correspondance avec les nb de root systems!
    ## 3D root distributions (fixe)
    R3 = np.zeros(np.shape(S2.m_1))
    R3[0:5, 1:3, 1:3] = np.ones([5,2,2])*0.25
    R4 = np.zeros(np.shape(S2.m_1))
    R4[0:7, 0:4, 2:4] = np.ones([7,4,2])*0.15

    ls_roots = [R3, R4]

    ##bouble journaliere
    ls_transp, evapo_tot, D, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S2.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, stateEV, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S2.update_ftsw()

    for i in range(40):
        Precip = precips[i]
        ls_transp, evapo_tot, D, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S2.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, stateEV, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
        S2.update_ftsw()

    #    ##visu
    #    Monviewer=Viewer
    #    MaScene=Scene()
    #    MaScene = S2.plot_soil_properties (vals=S2.ftsw_t, MaScene=Scene(), col_scale=5)
    #    Monviewer.display(MaScene)

    print((S2.ftsw_t))

print("test 1")
test_uni1()

#print("test 2")
#test_uni2()
