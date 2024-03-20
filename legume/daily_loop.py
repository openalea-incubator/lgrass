#from scipy import *
import time
import IOtable
import IOxls
import ShootMorpho as sh
import RootDistrib as rtd
import RootMorpho as rt
from copy import deepcopy
import numpy as np
import os

try:
    from soil3ds import soil_moduleN as solN #import de la version develop si module soil3ds est installe
except:
    import soil_moduleN3 as solN #soil_moduleN2_bis as solN #! renommer car dans nouvelle version Lpy, mot module est reserve et fait planter!


#daily loop
# decoupe daily_growth_loop initial en 4 fonctions pour donner acces au calcul du sol depuis l'exterieur

def daily_growth_loop(ParamP, invar, outvar, ls_epsi, meteo_j, mng_j, nbplantes, surfsolref, ls_ftswStress, ls_NNIStress, ls_TStress, lsApex, lsApexAll, opt_stressW=1, opt_stressN=1, opt_stressGel=0):
    """ daily potential growth loop (computes epsi, DM production / allocation / Ndemand) """

    epsilon = 10e-10  #
    isTTcut = False
    if mng_j['Coupe'] == 1.:
        isTTcut = True
        #print("cut ???")

    #gel a gerer par plante!
    invar['graineC'], invar['graineN'] = sh.reserves_graine(invar, ParamP)
    isGelDam = [0]*nbplantes #liste par plante comme invar['isGelDam']
    for nump in range(nbplantes):
        ##update gel status par plante
        if (meteo_j['TmoyDay'] > ParamP[nump]['Tgel'] and opt_stressGel == 1) or opt_stressGel == 0 or invar['graineN'][nump]>0.:
            isGelDam[nump] = 0 #laisse a zero
        else:
            isGelDam[nump] = 1
            #print("gel ???")

    # calcul de ls_epsi
    # invar['parap'] = array(list(map(sum, invar['PARaPlante'])))
    # invar['parip'] = array(list(map(sum, invar['PARiPlante'])))
    # # qatot= sum(res_trans[-1][:][:])*3600.*24/1000000. + sum(invar['parip'])#(MJ.day-1) #approximatif! a reprendre avec un vrai bilan radiatif
    # # print sum(res_trans[-1][:][:]), sum(res_trans[-1][:][:])*3600.*24/1000000., sum(res_trans[-1][:][:])*3600.*24/1000000.  +   sum(invar['parip'])
    # # ls_epsi = invar['parip']/qatot.tolist() #a reprendre : approximatif slmt! -> changera un peu avec un vrai bilan radiatif
    # # transmi_sol = 1-sum(ls_epsi)
    # # epsi = 1-transmi_sol #a reprendre pour differencier cible et vois #
    # transmi_sol = sum(res_trans[-1][:][:]) / (meteo_j['I0'] * surfsolref)  # bon
    # epsi = 1. - transmi_sol  # bon
    # ls_epsi = epsi * invar['parip'] / (sum(invar['parip']) + 10e-15)
    #ls_epsi = step_epsi(invar, res_trans, meteo_j, surfsolref)

    #print('graine', graineC, graineN, invar['NBI'], IOxls.get_lsparami(ParamP, 'DurGraine'),invar['TT'])

    # calcul de Biomasse tot
    stressHRUE = np.array(ls_ftswStress['WaterTreshRUE'])
    stressNRUE = np.array(ls_NNIStress['NTreshRUE'])
    if opt_stressW==0:
        stressHRUE = 1.
    if opt_stressN==0:
        stressNRUE = 1.

    stressFIX = 1 - np.array(invar['Ndfa']) * np.array(IOxls.get_lsparami(ParamP, 'NODcost'))  # coeff 0.15 = 15% reduction RUE a 100% fixation -> a passer en paarmetre
    stressTRUE = np.array(ls_TStress['stressTRUE'])#1.#

    invar['RUEpot'] = np.array(IOxls.get_lsparami(ParamP, 'RUE')) * stressTRUE
    invar['RUEactu'] = invar['RUEpot'] * stressHRUE * stressNRUE* stressFIX
    invar['PARaPlanteU'] = np.array(ls_epsi) * 0.95 * meteo_j['I0'] * 3600. * 24 / 1000000. * surfsolref  # facteur 0.95 pour reflectance / PARa used for calculation
    dM = invar['PARaPlanteU'] * invar['RUEactu'] + invar['graineC']
    # dM2 = array(dpar) * array(get_lsparami(ParamP, 'RUE'))

    # allocation
    froot = sh.rootalloc(IOxls.get_lsparami(ParamP, 'alloc_rootB'), IOxls.get_lsparami(ParamP, 'alloc_rootA'), invar['MS_aer_cumul'])  # fraction aux racines
    for nump in range(nbplantes):
        if invar['germination'][nump] < 2:  # tout aux racines avant apparition de la premiere feuille
            froot[nump] = 0.99


    Frac_remob = np.array(IOxls.get_lsparami(ParamP, 'frac_remob'))
    invar['CreservPiv'] = Frac_remob * invar['MS_pivot'] #fonction du compart pivot a t-1
    invar['remob'] = sh.Cremob(np.array(IOxls.dic2vec(nbplantes, invar['DemCp'])), invar['R_DemandC_Shoot'], invar['MS_pivot'], Frac_remob)  # vraiment marginal
    invar['CreservPiv'] -= invar['remob']
    rac_fine = dM * froot * np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine'))  # * rtd.filtre_ratio(invar['R_DemandC_Shoot'])
    pivot = dM * froot * (1 - np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine'))) - invar['remob']
    aer = dM - rac_fine - pivot #+ invar['remob']
    aer[aer==0.] += invar['remob'][aer==0.] # pour cas ou organes en croissance apres coupe, mais sans feuille (demande mini = remob)
    ffeuil = np.array(IOxls.dic2vec(nbplantes, invar['DemCp_lf'])) / (np.array(IOxls.dic2vec(nbplantes, invar['DemCp'])) + epsilon)  # fraction aux feuilles
    feuil = aer * ffeuil
    tige = aer * (1 - ffeuil)
    senaerien = np.array(invar['dMSenFeuil']) + np.array(invar['dMSenTige'])

    MS_aerien_tm1 = invar['MS_aerien'] #recupere MS aerien du t-1
    Npc_aerien_tm1 = invar['Npc_aer'] #Npc aerien du t-1

    invar['Mtot'].append(dM.tolist())
    invar['Mrac_fine'].append(rac_fine.tolist())  # matrice des delta MSrac fine par date
    invar['Mpivot'].append(pivot.tolist())  # matrice des delta MSpivot par date
    invar['Maerien'].append(aer.tolist())  # matrice des delta MSaerien par date
    invar['Mfeuil'].append(feuil.tolist())  # matrice des delta MSfeuil par date
    invar['Mtige'].append(tige.tolist())  # matrice des delta MSfeuil par date
    invar['Msenaerien'].append(senaerien.tolist())
    #invar['MS_pivot'] = list(map(sum, IOtable.t_list(invar['Mpivot'])))  # vecteur des MSpivot cumule au temps t
    invar['MS_pivot'] = np.array(invar['MS_pivot']) + pivot
    invar['MS_aerien'] = list(map(sum, IOtable.t_list(invar['Maerien'])))  # vecteur des MSaerien cumule au temps t
    invar['MS_feuil'] = list(map(sum, IOtable.t_list(invar['Mfeuil'])))  # vecteur des MSfeuil cumule au temps t
    invar['MS_tige'] = list(map(sum, IOtable.t_list(invar['Mtige'])))  # vecteur des MStige cumule au temps t
    invar['MS_senaerien'] = list(map(sum, IOtable.t_list(invar['Msenaerien'])))  # vecteur des MS_senaerien cumule au temps t
    invar['MS_aer_cumul'] += aer
    invar['MS_tot'] = list(map(sum, IOtable.t_list(invar['Mtot'])))
    invar['MS_rac_fine'] = list(map(sum, IOtable.t_list(invar['Mrac_fine'])))  # vecteur des MSracines_fines cumule au temps t
    invar['DiampivMax'] = np.sqrt(invar['MS_pivot'] * np.array(IOxls.get_lsparami(ParamP, 'DPivot2_coeff')))
    # invar['RLTot'] = array(map(sum, IOtable.t_list(invar['Mrac_fine']))) * array(IOxls.get_lsparami(ParamP, 'SRL')) #somme de toutes les racinesfines produites par plante
    invar['NBsh'], invar['NBI'] = sh.calcNB_NI(lsApex, nbplantes, seuilcountTige=0.25, seuilNItige=0.25)
    nbsh_2, nb1_2 = sh.calcNB_NI(lsApexAll, nbplantes, seuilcountTige=0.25, seuilNItige=0.25)  # recalcul sur tous les axes pour eviter bug des arret de tiges
    #nbsh_2, nb1_2 = sh.calcNB_NI(lsApexAll, nbplantes, seuilcountTige=0.,seuilNItige=0.25)  # recalcul sur tous les axes pour eviter bug des arret de tiges


    #print('nbsh',invar['NBsh'], nbsh_2)
    for nump in range(nbplantes):
        if nb1_2[nump] > invar['NBI'][nump]:
            invar['NBI'][nump] = nb1_2[nump]

        if nbsh_2[nump] > invar['NBsh'][nump]: # pour compter aussi les tiges arretes (au moins pour graminees)
            invar['NBsh'][nump] = nbsh_2[nump]


    invar['L_Sp'] = np.array(invar['MS_feuil']) / (np.array(invar['MS_aerien']) - np.array(invar['MS_feuil']) + epsilon)
    #pas tres realiste / a revoir (allometrie?)

    # print("MS AERIEN",invar['MS_aerien'],invar['MS_aer_cumul'])
    # print invar['Mtot']

    #ls_demandeN = array(invar['DemandN_Tot']) * 0.001 + 1e-15  # en kg N.plant-1 #[1e-12]*nbplantes #sera a renseigner -> la, force a zero - devra utiliser invar['DemandN_Tot'] qui est mis a jour + loin #en kg N
    #Npc_aer = array(invar['Naerien']) / (aer + array(invar['MS_aerien'])) * 100.  # Npc avec accroissement de biomasse pour calculer la demande

    #Npc_aer = array(invar['Naerien']) / (array(invar['MS_aerien'])) * 100.  # aer deja dans MS_aerien! -> mis a jour
    #distingue calcul partie coupee pour jour de coupe ou 1er step ou jour de gel
    if MS_aerien_tm1 == []:
        #jour 1
        Npc_aer = np.array(invar['Naerien']) / (np.array(invar['MS_aerien'])) * 100. #aer deja dans MS_aerien! -> mis a jour
        Npc_piv = np.array(invar['Npivot']) / (pivot + np.array(invar['MS_pivot'])) * 100.
        Npc_rac_fine = np.array(invar['Nrac_fine']) / (rac_fine + np.array(invar['MS_rac_fine'])) * 100.
    elif isTTcut:
        # jour de coupe = prends info du jour
        # MS_aerien_tm1 = array(aer)
        Npc_aer = np.array(invar['Naerien']) / (np.array(invar['MS_aerien'])) * 100.  # aer deja dans MS_aerien! -> mis a jour
        Npc_piv = invar['Npc_piv']
        Npc_rac_fine = invar['Npc_rac_fine']
    else:#pas jour de coupe
        Npc_piv = invar['Npc_piv']
        Npc_rac_fine = invar['Npc_rac_fine']
        if sum(isGelDam) == 0: #jour sans gel
            #jour d'avant pour etre synchro Naerien / MSaerein
            #Npc_aer = Npc_aerien_tm1
            Npc_aer = np.array(invar['Naerien']) / (np.array(invar['MS_aerien'])) * 100.
        else: #jour avec gel -> plante a plante
            Npc_aer = np.array(invar['Naerien']) / (np.array(invar['MS_aerien'])) * 100.  # aer deja dans MS_aerien! -> mis a jour
            for nump in range(nbplantes):
                if isGelDam[nump] == 0: #plante pas gelee
                    Npc_aer[nump] = Npc_aerien_tm1[nump]

    #Npc_piv = array(invar['Npivot']) / (pivot + array(invar['MS_pivot'])) * 100.
    #Npc_rac_fine = array(invar['Nrac_fine']) / (rac_fine + array(invar['MS_rac_fine'])) * 100.


    #reserve Piv
    invar['NreservPiv'] = np.array(invar['Npivot']) * (Npc_piv - np.array(IOxls.get_lsparami(ParamP, 'NminPiv'))) / Npc_piv
    invar['NreservPiv'][invar['NreservPiv'] < 0.] = 0.  # verifier que depasse pas zero!!


    ls_demandeN_aer, NcritTot_, MStot_ = solN.demandeNdefaut2(MSp=np.array(invar['MS_aerien'])-np.array(aer), dMSp=aer, Npc=Npc_aer, surfsolref=surfsolref, a=np.array(IOxls.get_lsparami(ParamP, 'ADIL')), b1=np.array(IOxls.get_lsparami(ParamP, 'BDILi')), b2=np.array(IOxls.get_lsparami(ParamP, 'BDIL')))
    #ls_demandeN_aer, NcritTot_, MStot_ = solN.demandeNdefaut2(MSp=array(MS_aerien_tm1), dMSp=aer, Npc=Npc_aer, surfsolref=surfsolref, a=array(IOxls.get_lsparami(ParamP, 'ADIL')), b1=array(IOxls.get_lsparami(ParamP, 'BDILi')), b2=array(IOxls.get_lsparami(ParamP, 'BDIL')))

    ls_demandeN_aer = ls_demandeN_aer * 0.001 #+ 1e-15  # en kg N.plant-1
    ls_demandN_piv = solN.demandeNroot(np.array(invar['MS_pivot']), pivot, Npc_piv, surfsolref, np.array(IOxls.get_lsparami(ParamP, 'NoptPiv'))) * 0.001 + epsilon #+ 1e-15  # en kg N.plant-1
    ls_demandN_rac_fine = solN.demandeNroot(np.array(invar['MS_rac_fine']), rac_fine, Npc_rac_fine, surfsolref, np.array(IOxls.get_lsparami(ParamP, 'NoptFR'))) * 0.001 #+ 1e-15  # en kg N.plant-1

    ls_demandeN_bis = ls_demandeN_aer + ls_demandN_piv + ls_demandN_rac_fine #+ epsilon
    fracNaer = ls_demandeN_aer  / (ls_demandeN_bis + epsilon)
    fracNpiv = ls_demandN_piv / (ls_demandeN_bis + epsilon)
    fracNrac_fine = ls_demandN_rac_fine / (ls_demandeN_bis + epsilon)

    invar['DemandN_TotAer'] = ls_demandeN_aer
    #print('in pot Npc', invar['Naerien'], invar['MS_aerien'], Npc_aer, Npc_aerien_tm1,ls_demandeN_aer, NcritTot_, MStot_)

    # print invar['Maerien']#invar['MS_aerien']
    # print aer
    #test si bourgeons en vie (could add conditions related water stress in another variable)
    for nump in range(nbplantes):
        if sum(invar['SurfPlante'][nump]) == 0. and invar['NBsh'][nump] == 0. and invar['NBB'][nump] == 0. and invar['NBD1'][nump] == 0.:
            invar['aliveB'][nump] += 1. #ajoute 1 jour ou condition morte remplie
        else:
            invar['aliveB'][nump] = 0. #ou remise a zero


    #senescence perenne et flux N
    invar['dMSenNonRec'], invar['dMSenPiv'], invar['perteN_NonRec'], invar['perteN_Piv']  = sh.Turnover_compart_Perenne(invar, ParamP)
    invar['Npc_aerNonRec'] = invar['NaerienNonRec'] / invar['MS_aerienNonRec'] *100.
    if invar['Npc_aer'] == []:#1er step
        invar['perteN_aerien'] = (np.array(invar['dMSenFeuil']) + np.array(invar['dMSenTige'])) * 0.
    else:
        invar['perteN_aerien'] = (np.array(invar['dMSenFeuil']) + np.array(invar['dMSenTige'])) * invar['Npc_aer']/100.

    #print('perteN_aerien', invar['perteN_aerien'], invar['dMSenFeuil'] , invar['dMSenTige'], invar['Npc_aer'], len(invar['Npc_aer']))
    #connecter avec invar ['Msenaerien'] + ajouter une variable de cumul senenscence par coupe

    #print('dMSenNonRec', invar['dMSenNonRec'])
    #print('piv', invar['MS_pivot'], invar['dMSenPiv'])

    # ajout des bilan C plante pour sorties / m2
    outvar['BilanC_PARa'].append(sum(invar['PARaPlanteU']) / surfsolref)
    outvar['BilanC_RUE'].append(sum(dM) / (sum(invar['PARaPlanteU'])+epsilon))
    outvar['BilanCdMStot'].append(sum(dM) / surfsolref)
    outvar['BilanCdMrac_fine'].append(sum(rac_fine) / surfsolref)
    outvar['BilanCdMpivot'].append(sum(pivot) / surfsolref)
    outvar['BilanCdMaer'].append(sum(aer) / surfsolref)
    outvar['BilanCdMSenFeuil'].append(sum(invar['dMSenFeuil']) / surfsolref)
    outvar['BilanCdMSenTige'].append(sum(invar['dMSenTige']) / surfsolref)
    outvar['BilanCdMSenNonRec'].append(sum(invar['dMSenNonRec']) / surfsolref)
    outvar['BilanCdMSenPiv'].append(sum(invar['dMSenPiv']) / surfsolref)

    #test des 3 autres fonctions en interne au sein d'une fonction globale
    temps = [aer, rac_fine, pivot, fracNaer, fracNpiv, fracNrac_fine, MS_aerien_tm1, isTTcut,NcritTot_,epsilon,meteo_j] #variable temporaires pour passer entre fonctions (passer ds invar?)

    return invar, outvar, ls_demandeN_bis, temps


def step_epsi(invar, res_trans, lsFeuilBilanR, meteo_j, surfsolref):
    """ calculate ls_epsi from res_trans and invar"""
    #invar['parap'] = array(list(map(sum, invar['PARaPlante'])))
    #invar['parip'] = array(list(map(sum, invar['PARiPlante'])))
    #MAJ invar['parap'] , invar['parip']
    sh.calc_para_Plt(invar, lsFeuilBilanR)

    # qatot= sum(res_trans[-1][:][:])*3600.*24/1000000. + sum(invar['parip'])#(MJ.day-1) #approximatif! a reprendre avec un vrai bilan radiatif
    # print sum(res_trans[-1][:][:]), sum(res_trans[-1][:][:])*3600.*24/1000000., sum(res_trans[-1][:][:])*3600.*24/1000000.  +   sum(invar['parip'])
    # ls_epsi = invar['parip']/qatot.tolist() #a reprendre : approximatif slmt! -> changera un peu avec un vrai bilan radiatif
    # transmi_sol = 1-sum(ls_epsi)
    # epsi = 1-transmi_sol #a reprendre pour differencier cible et vois #
    transmi_sol = np.sum(res_trans[-1][:][:]) / (meteo_j['I0'] * surfsolref)  # bon
    epsi = 1. - transmi_sol  # bon
    ls_epsi = epsi * invar['parip'] / (sum(invar['parip']) + 10e-15)
    return ls_epsi, invar



def Update_stress_loop(ParamP, invar, invar_sc, temps, DOY, nbplantes, surfsolref, ls_epsi, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, ls_demandeN_bis, ls_ftswStress, ls_TStress, dicOrgans, dicFeuilBilanR, lsApex, start_time, cutNB, deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant, I_I0profilInPlant, NlClasses, NaClasses, NlinClasses, outvar):
    """ Update daily N uptake/fixation from soil WN balance and plant demands / prepares stress variables for next step / write output variables   """

    aer, rac_fine, pivot, fracNaer, fracNpiv, fracNrac_fine, MS_aerien_tm1, isTTcut, NcritTot_, epsilon,meteo_j = temps# temps[0], temps[1],temps[2], temps[3],temps[4], temps[5],temps[6] #unpacks variables temporaires passes entre fonction -> a repasser dans invar!!!

    # water
    invar['transpi'] = ls_transp
    invar['cumtranspi'] += np.array(ls_transp)

    #print('test demandeN', (invar['DemandN_TotAer'] * 1000. + invar['Naerien'])*100. / (invar['MS_aerien'] + aer), invar['MS_aerien'] , aer, invar['DemandN_TotAer'] * 1000. , invar['Naerien'])
    # Uptake N et allocation
    invar['Nuptake_sol'] = np.array(list(map(np.sum, ls_Act_Nuptake_plt))) * 1000 + invar['graineN']  # g N.plant-1 #test ls_demandeN_bis*1000.#
    try:
        NremobC = invar['remob'] * invar['Npc_piv'] / 100.  # remobilise N pivot qui part avec le C
        invar['Naerien'] += invar['Nuptake_sol'] * fracNaer + NremobC # uptake N va dans partie aeriennes au prorata des demandes
        invar['Npivot'] += invar['Nuptake_sol'] * fracNpiv - NremobC
    except:  # 1er step
        NremobC = 0.
        invar['Naerien'] += invar['Nuptake_sol'] * fracNaer + NremobC
        invar['Npivot'] += invar['Nuptake_sol'] * fracNpiv
        #print('rem')

    invar['Nrac_fine'] += invar['Nuptake_sol'] * fracNrac_fine

    # Fixation et allocation
    maxFix = sh.Ndfa_max(invar['TT'], IOxls.get_lsparami(ParamP, 'DurDevFix')) * np.array(IOxls.get_lsparami(ParamP, 'MaxFix')) / 1000. * aer  #invar['MS_aerien']# * invar['dTT']
    stressHFix = np.array(ls_ftswStress['WaterTreshFix']) * maxFix  # effet hydrique
    invar['Qfix'] = sh.ActualFix(ls_demandeN_bis * 1000., invar['Nuptake_sol'], stressHFix)  # g N.plant-1
    invar['Ndfa'] = invar['Qfix'] / (invar['Qfix'] + invar['Nuptake_sol'] + 1e-15)

    delta_besoinN_aerien = invar['DemandN_TotAer'] * 1000. - invar['Qfix'] * fracNaer - invar['Nuptake_sol'] * fracNaer - NremobC  # besoin N are sont ils couverts? g N.plant-1
    delta_besoinN_aerien[delta_besoinN_aerien < 0.] = 0.#max(0., delta_besoinN_aerien) #si negatif (e.g. avec remob!)
    NremobN = np.minimum(delta_besoinN_aerien, invar['NreservPiv'])  # si pas couvert remobilisation N du pivot directement
    NremobN[NremobN < 0.] = 0.  # verifie que pas de negatif
    #print('bilan', invar['DemandN_TotAer'] * 1000., invar['Qfix'] * fracNaer + invar['Nuptake_sol'] * fracNaer + NremobC + delta_besoinN_aerien, invar['Qfix'] * fracNaer , invar['Nuptake_sol'] * fracNaer , NremobC ,NremobN, delta_besoinN_aerien, fracNaer, invar['Nuptake_sol'])

    # print 'Npivot', invar['Npivot'][0:2]
    # print 'NreservPiv', invar['NreservPiv'][0:2]
    # print 'delta_besoinN', delta_besoinN_aerien[0:2]
    # print 'NremobN', NremobN[0:2]

    invar['Naerien'] += invar['Qfix'] * fracNaer + NremobN #- invar['dNmortGel']
    invar['Npivot'] += invar['Qfix'] * fracNpiv - NremobN
    invar['NreservPiv'] -= NremobN
    invar['Nrac_fine'] += invar['Qfix'] * fracNrac_fine  # total : vivantes et mortes

    #correction si gel
    #if sum(invar['isGelDam'])!=0: #certaines plantes gel
    #    for nump in range(nbplantes):
    #        if invar['isGelDam'][nump] == 1:
    #            invar['Naerien'][nump] = invar['dMSmortGel'][nump] * invar['Npc_aer'][nump] / 100.
    #            #invar['Npc_aer'] = invar['dNmortGel']/invar['dMSmortGel']
    #            #invar['dMSmortGel'][nump] = MSA
    #            #invar['dNmortGel'][nump]

    #print(invar['Naerien'], invar['Qfix'] * fracNaer, NremobN, invar['dMSmortGel'])
    # effet feedback N pas fait (priorite) -> necessaire???
    # mise a jour Npc et calcul NNI

    invar['Npc_aer'] = np.array(invar['Naerien']) / (aer + np.array(invar['MS_aerien'])) * 100.  # %
    invar['Npc_piv'] = np.array(invar['Npivot']) / (pivot + np.array(invar['MS_pivot'])) * 100.  # %
    invar['Npc_rac_fine'] = np.array(invar['Nrac_fine']) / (rac_fine + np.array(invar['MS_rac_fine'])) * 100.  # %
    #print('in stres Npc', invar['Naerien'], invar['MS_aerien'], aer, invar['Npc_aer'], delta_besoinN_aerien)
    #print('besoinN', delta_besoinN_aerien, invar['DemandN_TotAer'] * 1000., invar['Qfix'] * fracNaer, invar['Nuptake_sol'] * fracNaer, NremobC)

    # print 'Npc_piv', invar['Npc_piv'][0:2]

    critN_inst = NcritTot_#solN.critN(sum(aer + array(invar['MS_aerien'])) / (surfsolref * 100.))  # azote critique couvert
    critN_inst[critN_inst > 7.] = 7. #seuil de crit pour calcul stress NNI plus bas pour petites plantes
    #print('critN' , critN_inst, invar['Npc_aer'])
    invar['NNI'] = invar['Npc_aer'] / critN_inst


    # update des indices de stress hydrique par plante pour step suivant
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = [], [], [], [], [], [], [], [], [], []  # liste de parametres
    for nump in range(nbplantes):
        p1.append([ParamP[nump]['WaterTreshExpSurfs'], ParamP[nump]['WaterTreshExpSurfd']])
        p2.append([ParamP[nump]['WaterTreshDevIIs'], ParamP[nump]['WaterTreshDevIId']])
        p3.append([ParamP[nump]['WaterTreshDevIs'], ParamP[nump]['WaterTreshDevId']])
        p4.append([ParamP[nump]['WaterTreshFixs'], ParamP[nump]['WaterTreshFixd']])
        p5.append([ParamP[nump]['WaterTreshRUEs'], ParamP[nump]['WaterTreshRUEd']])
        p6.append([ParamP[nump]['NTreshRUEs'], ParamP[nump]['NTreshRUEd']])
        p7.append([ParamP[nump]['NTreshExpSurfs'], ParamP[nump]['NTreshExpSurfd']])
        p8.append([ParamP[nump]['NTreshDevs'], ParamP[nump]['NTreshDevd']])
        p9.append([ParamP[nump]['NTreshDevIIs'], ParamP[nump]['NTreshDevIId']])
        p10.append([ParamP[nump]['TempTreshRUEb'], ParamP[nump]['TempTreshRUEh']])


    ls_ftswStress = {}
    ls_ftswStress['WaterTreshExpSurf'] = list(map(sh.FTSW_resp, ls_ftsw, p1))
    ls_ftswStress['WaterTreshDevII'] = list(map(sh.FTSW_resp, ls_ftsw, p2))
    ls_ftswStress['WaterTreshDevI'] = list(map(sh.FTSW_resp, ls_ftsw, p3))
    ls_ftswStress['WaterTreshFix'] = list(map(sh.FTSW_resp, ls_ftsw, p4))
    ls_ftswStress['WaterTreshRUE'] = list(map(sh.FTSW_resp, ls_ftsw, p5))

    # update des indices de stress N par plante pour step suivant
    ls_NNIStress = {}
    ls_NNIStress['NTreshRUE'] = list(map(sh.NNI_resp, invar['NNI'], p6))
    ls_NNIStress['NTreshExpSurf'] = list(map(sh.NNI_resp, invar['NNI'], p7))
    ls_NNIStress['NTreshDev'] = list(map(sh.NNI_resp, invar['NNI'], p8))
    ls_NNIStress['NTreshDevII'] = list(map(sh.NNI_resp, invar['NNI'], p9))

    #update Temperature stress
    ls_TStress = {}
    ls_TStress['stressTRUE'] = list(map(sh.linear_stress2, [meteo_j['TmoyDay']]*nbplantes, p10))#



    # print invar['TT'], Ndfa_max(invar['TT'], IOxls.get_lsparami(ParamP, 'DurDevFix')), maxFix, stressHFix
    # print invar['TT'], ls_demandeN_bis, invar['Nuptake_sol'], stressHFix
    # print sum(mapN_Rain), sum(mapN_Irrig), sum(mapN_fertNO3), sum(mapN_fertNH4), meteo_j['Tsol']
    # print ls_demandeN_bis, ls_demandeN, Npc_temp, array(map(sum, ls_Act_Nuptake_plt)), invar['Naerien'] #pour convertir en g N
    # print invar['Npc_bis']
    # print ls_demandeN_bis[0], ls_demandeN
    # print solN.critN(sum(aer+array(invar['MS_aerien']))#, invar['Npc_bis']

    # calcul offre/demandeC
    #tab = IOtable.conv_dataframe(IOtable.t_list(lsOrgans))
    # OffCp = calcOffreC (tab, 'plt')#pas utilise??!
    # invar['DemCp'] = calcDemandeC(tab, 'plt')#attention, pour que calcul soit bon, faut le STEPS  suivant mis a jour!-> a faire en StartEach
    # invar['L_Sp'] = sh.calcLeafStemRatio(ParamP, tab, lsAxes)

    # calcul surf par tige/axe
    invar_sc['plt']['Surf'], invar_sc['plt']['SurfVerte'], invar_sc['sh']['Surf'], invar_sc['sh']['SurfVerte'], \
    invar_sc['ax']['Surf'], invar_sc['ax']['SurfVerte'], invar_sc['plt']['PARaF'], invar_sc['sh']['PARaF'], \
    invar_sc['ax']['PARaF'], invar_sc['ax']['MaxPARaF'] = sh.calcSurfLightScales(dicFeuilBilanR, ParamP)

    invar_sc['ax']['AgePiv'] = sh.AgePivScales(dicOrgans, ParamP)


    # calcul de fraction de PARa par pivot
    invar_sc['ax']['fPARaPiv'] = rt.calc_daxfPARaPiv(nbplantes, invar_sc['ax']['AgePiv'], invar_sc['plt']['PARaF'], invar_sc['ax']['PARaF'])

    # calcul demande par pivot
    invar_sc['ax']['DemCRac'], invar_sc['ax']['NRac'] = rt.calc_DemandC_roots(ParamP, invar_sc['ax']['AgePiv'],
                                                                              invar['Udevsol'],
                                                                              invar_sc['ax']['QDCmoyRac'])

    # calcul biomasse, diametres pivots indivs, QDC des racines, increment de longueur et SRL
    daxPiv = rt.distrib_dM_ax(invar_sc['ax']['fPARaPiv'], pivot, Frac_piv_sem=IOxls.get_lsparami(ParamP, 'Frac_piv_sem'),
                              Frac_piv_loc=IOxls.get_lsparami(ParamP,
                                                             'Frac_piv_loc'))  # rt.distrib_dM_ax(invar_sc['ax']['fPARaPiv'], pivot)
    invar_sc['ax']['MaxPiv'] = IOxls.add_dic(daxPiv, invar_sc['ax']['MaxPiv'])
    invar_sc['ax']['DiampivMax'] = rt.calc_DiamPivMax(ParamP, invar_sc['ax']['MaxPiv'])

    invar_sc['ax']['OfrCRac'] = rt.distrib_dM_ax(invar_sc['ax']['fPARaPiv'], rac_fine,
                                                 Frac_piv_sem=IOxls.get_lsparami(ParamP, 'Frac_piv_sem'),
                                                 Frac_piv_loc=IOxls.get_lsparami(ParamP, 'Frac_piv_loc'))
    invar_sc['ax']['QDCRac'] = rt.calc_QDC_roots(invar_sc['ax']['OfrCRac'], invar_sc['ax']['DemCRac'])
    invar_sc['ax']['QDCmoyRac'] = rt.calc_QDCmoy_roots(invar_sc['ax']['QDCRac'], invar_sc['ax']['QDCmoyRac'],
                                                       invar_sc['ax']['AgePiv'], invar['Udevsol'])
    invar_sc['ax']['StressHmoyRac'] = rt.calc_StressHmoy_roots(invar_sc['ax']['StressHRac'],
                                                               invar_sc['ax']['PonderStressHRac'],
                                                               invar_sc['ax']['StressHmoyRac'],
                                                               invar_sc['ax']['AgePiv'], invar[
                                                                   'Udevsol'])  # (dStressH, dPonder, dStressHmoy, dAgePiv, dTT)

    invar_sc['ax']['dlRac'] = rt.calc_dLong_roots(ParamP, invar_sc['ax']['NRac'], invar['Udevsol'],
                                                  invar_sc['ax']['QDCRac'], invar_sc['ax']['StressHRac'],
                                                  invar_sc['ax'][
                                                      'PonderStressHRac'])  # passe STEPS, mais devrait filer les dTT de chaque plante
    invar_sc['ax']['cumlRac'] = IOxls.add_dic(invar_sc['ax']['dlRac'], invar_sc['ax']['cumlRac'])
    invar['RLen1'], invar['RLen2'], invar['RLen3'], invar['RLentot'] = rt.cumul_plante_Lrac(nbplantes,
                                                                                            invar_sc['ax']['cumlRac'])
    dl1, dl2, dl3, dltot = rt.cumul_plante_Lrac(nbplantes,
                                                invar_sc['ax']['dlRac'])  # calcul des delta de longueur par plante
    invar['dRLen2'].append(dl2)  # stocke les dl du jour pour cacalcul senescence de plus tard
    invar['dRLen3'].append(dl3)
    # invar['SRL'] = invar['RLentot']/(invar['MS_rac_fine'][0]+10e-15)
    # print invar_sc['ax']['QDCRac']

    # print 'graine', graineC, dltot, invar['Surfcoty'], invar['Mcoty']#

    dur2 = (np.array(IOxls.get_lsparami(ParamP, 'GDs2')) + np.array(IOxls.get_lsparami(ParamP, 'LDs2'))) / 20.  # en jours a 20 degres!
    dur3 = (np.array(IOxls.get_lsparami(ParamP, 'GDs3')) + np.array(IOxls.get_lsparami(ParamP, 'LDs3'))) / 20.  # en jours a 20 degres!
    invar['dRLenSentot'], invar['dMSenRoot'] = rt.calc_root_senescence(invar['dRLen2'], invar['dRLen3'], dur2, dur3, np.array(invar['SRL']))
    invar['RLentotfromDev'] = np.array(invar['RLentotfromDev']) + dltot - invar['dRLenSentot']
    invar['MS_rac_fineNet'] = np.array(invar['MS_rac_fineNet']) + rac_fine - invar['dMSenRoot']
    invar['SRL'] = invar['RLentotfromDev'] / (invar['MS_rac_fineNet'][0] + 10e-15)

    paramSRLmin = np.array(IOxls.get_lsparami(ParamP, 'SRLmin'))
    invar['RLentotfromRootMass'] = invar['MS_rac_fine'] * paramSRLmin
    invar['RLTotNet'] = deepcopy(invar['RLentotfromDev'])

    #update RLTotNet / RLtoyNet si option
    for nump in range(nbplantes):
        if invar['SRL'][nump] < paramSRLmin[nump]:
            invar['RLTotNet'][nump] = invar['RLentotfromRootMass'][nump]
            invar['SRL'][nump] = paramSRLmin[nump]


    invar['perteN_rac_fine'] = invar['dMSenRoot'] * invar['Npc_rac_fine'] / 100.
    # sortir une variable cumule d'N des rac mortes? -> compement a invar['Nrac_fine'] qui comprend les deux


    # calcul senesc a faire a l'echelle des axes plutot? -> a priori pas necessaire

    invar['R_DemandC_Root'] = rt.calc_QDplante(nbplantes, invar_sc['ax']['QDCRac'], invar_sc['ax']['cumlRac'], invar['RLTotNet'])#invar['RLentot'])
    invar['R_DemandC_Shoot'] = aer / (np.array(IOxls.dic2vec(nbplantes, invar['DemCp'])) + epsilon)#10e-15)
    #print('R_DemandC_Shoot', invar['R_DemandC_Shoot'], aer, invar['DemCp'], array(IOxls.dic2vec(nbplantes, invar['DemCp'])))

    # if '0_0_0' in invar_sc['ax']['NRac'].keys():
    #    print invar_sc['ax']['NRac']['0_0_0']
    #    print invar_sc['ax']['QDCRac']['0_0_0']
    #    print invar_sc['ax']['dlRac']['0_0_0']
    # print invar['RLentot'], invar['MS_rac_fine'], invar['RLentot'][0]/(invar['MS_rac_fine'][0]+0.00000001)

    # calcul demandN -> a depalcer dans le starteach comme pour C?? -> pas utilise actuellement
    if lsApex != []:
        I_I0profilInPlant = sh.cumul_lenIN(lsApex, dicOrgans, I_I0profilInPlant, deltaI_I0, nbI_I0)

    # pas utilise
    for nump in range(nbplantes):
        invar['DemandN_Feuil'][nump] = sum(I_I0profilLfPlant[nump] * NaClasses)
        invar['DemandN_Pet'][nump] = sum(I_I0profilPetPlant[nump] * NlClasses)
        invar['DemandN_Stem'][nump] = sum(I_I0profilInPlant[nump] * NlinClasses)
        # invar['DemandN_Tot'][nump] = invar['DemandN_Feuil'][nump] + invar['DemandN_Pet'][nump] + invar['DemandN_Stem'][nump]

    invar['DemandN_Tot'] = ls_demandeN_bis * 1000.
    # print invar['DemandN_Tot'][0], sum(ls_Act_Nuptake_plt[0]), sum(ls_Act_Nuptake_plt[0])/(invar['DemandN_Tot'][0]+10e-12), sum(S.m_NO3)

    Npc = (np.array(invar['DemandN_Feuil']) + np.array(invar['DemandN_Pet']) + np.array(invar['DemandN_Stem'])) * 100. / np.array(invar['MS_aerien'])

    # sorties
    outvar = increment_dailyOutput(outvar, invar, DOY, nbplantes, start_time, ls_epsi, aer, ls_ftsw, ls_transp, Npc, cutNB)
    # to do: passer ls_epsi, aer, ls_ftsw, ls_transp, Npc dans invar!

    return invar, invar_sc, outvar, I_I0profilInPlant, ls_ftswStress, ls_NNIStress, ls_TStress



def increment_dailyOutput(outvar, invar, DOY, nbplantes, start_time, ls_epsi, aer, ls_ftsw, ls_transp, Npc, cutNB):
    """ add daily invar values into outvar """

    # temps de calcul
    past_time = time.time() - start_time

    # sorties
    outvar['TT'].append(['TT', DOY] + invar['TT'])
    outvar['TTudev'].append(['TTudev', DOY] + invar['TTudev'])
    outvar['dTT'].append(['dTT', DOY] + invar['dTT'])
    outvar['Udev'].append(['Udev', DOY] + invar['Udev'])
    outvar['Udevstress'].append(['Udevstress', DOY] + invar['Udevstress'])
    outvar['TTphyllo'].append(['TTphyllo', DOY] + invar['TTphyllo'].tolist())
    outvar['time'].append(['time', DOY] + [past_time] * nbplantes)
    outvar['cutNB'].append(['cutNB', DOY] + [cutNB] * nbplantes)
    outvar['SurfPlante'].append(['SurfPlante', DOY] + list(map(sum, invar['SurfPlante'])))
    outvar['PARaPlante'].append(['PARaPlante', DOY] + invar['PARaPlanteU'].tolist())  # append(['PARaPlante',DOY]+invar['parap'].tolist())
    outvar['PARiPlante'].append(['PARiPlante', DOY] + invar['parip'].tolist())
    outvar['epsi'].append(['epsi', DOY] + ls_epsi.tolist())
    outvar['dMSaer'].append(['dMSaer', DOY] + aer.tolist())
    outvar['Hplante'].append(['Hplante', DOY] + invar['Hplante'])
    outvar['Dplante'].append(['Dplante', DOY] + invar['Dplante'])
    outvar['RLTot'].append(['RLTot', DOY] + invar['RLentot'])
    outvar['RDepth'].append(['RDepth', DOY] + invar['RDepth'])
    outvar['MS_aerien'].append(['MSaerien', DOY] + invar['MS_aerien'])
    outvar['MS_feuil'].append(['MSfeuil', DOY] + invar['MS_feuil'])
    outvar['MS_tige'].append(['MStige', DOY] + invar['MS_tige'])
    outvar['MS_tot'].append(['MStot', DOY] + invar['MS_tot'])
    outvar['countSh'].append(['countSh', DOY] + invar['countSh'])
    outvar['countShExp'].append(['countShExp', DOY] + invar['countShExp'])
    outvar['demandC'].append(['demandC', DOY] + IOxls.dic2vec(nbplantes, invar['DemCp']))
    outvar['Leaf_Stem'].append(['Leaf_Stem', DOY] + invar['L_Sp'].tolist())
    outvar['NBsh'].append(['NBsh', DOY] + invar['NBsh'])
    outvar['NBI'].append(['NBI', DOY] + invar['NBI'])
    outvar['FTSW'].append(['FTSW', DOY] + ls_ftsw)
    outvar['Etransp'].append(['Etransp', DOY] + ls_transp)
    outvar['DemandN_Feuil'].append(['DemandN_Feuil', DOY] + invar['DemandN_Feuil'])
    outvar['DemandN_Pet'].append(['DemandN_Pet', DOY] + invar['DemandN_Pet'])
    outvar['DemandN_Stem'].append(['DemandN_Stem', DOY] + invar['DemandN_Stem'])
    outvar['DemandN_Tot'].append(['DemandN_Tot', DOY] + invar['DemandN_Tot'].tolist())
    outvar['Npc'].append(['Npc', DOY] + Npc.tolist())
    outvar['NBD1'].append(['NBD1', DOY] + invar['NBD1'])
    outvar['NBB'].append(['NBB', DOY] + invar['NBB'])
    outvar['NBBexp'].append([['NBBexp', DOY] + invar['NBBexp']])
    outvar['R_DemandC_Root'].append(['R_DemandC_Root', DOY] + invar['R_DemandC_Root'])
    outvar['SRL'].append(['SRL', DOY] + invar['SRL'].tolist())
    outvar['DemandN_Tot_Aer'].append(['DemandN_Tot_Aer', DOY] + invar['DemandN_TotAer'].tolist())
    outvar['Naerien'].append(['Naerien', DOY] + invar['Naerien'].tolist())
    outvar['Npc_aer'].append(['Npc_aer', DOY] + invar['Npc_aer'].tolist())  # -> ancien Npc_bis
    outvar['Npc_piv'].append(['Npc_piv', DOY] + invar['Npc_piv'].tolist())
    outvar['Npc_rac_fine'].append(['Npc_rac_fine', DOY] + invar['Npc_rac_fine'].tolist())
    outvar['Nuptake_sol'].append(['Nuptake_sol', DOY] + invar['Nuptake_sol'].tolist())
    outvar['NNI'].append(['NNI', DOY] + invar['NNI'].tolist())
    outvar['Ndfa'].append(['Ndfa', DOY] + invar['Ndfa'].tolist())
    outvar['Qfix'].append(['Qfix', DOY] + invar['Qfix'].tolist())
    outvar['dMSenFeuil'].append(['dMSenFeuil', DOY] + invar['dMSenFeuil'])
    outvar['dMSenTige'].append(['dMSenTige', DOY] + invar['dMSenTige'])
    outvar['MS_pivot'].append(['MS_pivot', DOY] + invar['MS_pivot'])
    outvar['MS_rac_fine'].append(['MS_rac_fine', DOY] + invar['MS_rac_fine'])
    outvar['R_DemandC_Shoot'].append(['R_DemandC_Shoot', DOY] + invar['R_DemandC_Shoot'].tolist())
    outvar['RUE'].append(['RUE', DOY] + invar['RUEactu'].tolist())
    outvar['RUEpot'].append(['RUEpot', DOY] + invar['RUEpot'].tolist())
    outvar['DemCp'].append(['DemCp', DOY] + IOxls.dic2vec(nbplantes, invar['DemCp']))
    outvar['remob'].append(['remob', DOY] + invar['remob'].tolist())
    outvar['dRLenSentot'].append(['dRLenSentot', DOY] + invar['dRLenSentot'].tolist())
    outvar['dMSenRoot'].append(['dMSenRoot', DOY] + invar['dMSenRoot'].tolist())
    outvar['RLTotNet'].append(['RLTotNet', DOY] + np.array(invar['RLTotNet']).tolist())
    outvar['MS_rac_fineNet'].append(['MS_rac_fineNet', DOY] + invar['MS_rac_fineNet'].tolist())
    outvar['perteN_rac_fine'].append(['perteN_rac_fine', DOY] + invar['perteN_rac_fine'].tolist())
    outvar['NBphyto'].append(['NBphyto', DOY] + invar['NBphyto'])
    outvar['aliveB'].append(['aliveB', DOY] + invar['aliveB'])
    outvar['NBapexAct'].append(['NBapexAct', DOY] + invar['NBapexAct'])  # pour correction du nb phyto par rapport au comptage observe
    outvar['transpi'].append(['transpi', DOY] + invar['transpi'])
    outvar['cumtranspi'].append(['cumtranspi', DOY] + invar['cumtranspi'].tolist())
    outvar['dMSmortGel'].append(['dMSmortGel', DOY] + invar['dMSmortGel_aer'])
    outvar['dNmortGel'].append(['dNmortGel', DOY] + invar['dNmortGel_aer'])
    outvar['MS_aerienNonRec'].append(['MSaerienNonRec', DOY] + invar['MS_aerienNonRec'].tolist())
    outvar['MS_aerienRec'].append(['MSaerienRec', DOY] + invar['MS_aerienRec'].tolist())
    outvar['NaerienNonRec'].append(['NaerienNonRec', DOY] + invar['NaerienNonRec'].tolist())
    outvar['graineC'].append(['graineC', DOY] + invar['graineC'].tolist())
    outvar['graineN'].append(['graineN', DOY] + invar['graineN'].tolist())
    outvar['CreservPiv'].append(['CreservPiv', DOY] + invar['CreservPiv'].tolist())
    outvar['NreservPiv'].append(['NreservPiv', DOY] + invar['NreservPiv'].tolist())
    outvar['dMSenPiv'].append(['dMSenPiv', DOY] + invar['dMSenPiv'].tolist())
    outvar['dMSenNonRec'].append(['dMSenNonRec', DOY] + invar['dMSenNonRec'].tolist())
    outvar['perteN_Piv'].append(['perteN_Piv', DOY] + invar['perteN_Piv'].tolist())
    outvar['perteN_NonRec'].append(['perteN_NonRec', DOY] + invar['perteN_NonRec'].tolist())
    outvar['perteN_aerien'].append(['perteN_aerien', DOY] + invar['perteN_aerien'].tolist())
    outvar['Npc_aerNonRec'].append(['Npc_aerNonRec', DOY] + invar['Npc_aerNonRec'].tolist())
    outvar['MS_senaerien'].append(['MSsenaerien', DOY] + invar['MS_senaerien'])
    outvar['dMSmortPlant_aer'].append(['dMSmortPlant_aer', DOY] + invar['dMSmortPlant_aer'].tolist())
    outvar['dMSmortPlant_pivot'].append(['dMSmortPlant_pivot', DOY] + invar['dMSmortPlant_pivot'].tolist())
    outvar['dMSmortPlant_racfine'].append(['dMSmortPlant_racfine', DOY] + invar['dMSmortPlant_racfine'].tolist())
    outvar['dNmortPlant_aer'].append(['dNmortPlant_aer', DOY] + invar['dNmortPlant_aer'].tolist())
    outvar['dNmortPlant_pivot'].append(['dNmortPlant_pivot', DOY] + invar['dNmortPlant_pivot'].tolist())
    outvar['dNmortPlant_racfine'].append(['dNmortPlant_racfine', DOY] + invar['dNmortPlant_racfine'].tolist())
    outvar['RLentotfromRootMass'].append(['RLentotfromRootMass', DOY] + invar['RLentotfromRootMass'].tolist())
    outvar['RLentotfromDev'].append(['RLentotfromDev', DOY] + np.array(invar['RLentotfromDev']).tolist())
    outvar['ConcNmoy'].append(['ConcNmoy', DOY] + invar['ConcNmoy'])


    # !! ces 4 sorties lucas ne sont pas au format attentdu!
    #outvar['phmgPet'].append(['phmgPet', DOY] + list(map(max, invar['phmgPet'])))
    #outvar['phmgEntr'].append(['phmgEntr', DOY] + list(map(max, invar['phmgEntr'])))
    #outvar['phmgPet_m'].append(['phmgPet_m', DOY] + list(map(min, invar['phmgPet_m'])))
    #outvar['phmgEntr_m'].append(['phmgEntr_m', DOY] + list(map(min, invar['phmgEntr_m'])))

    return outvar
    # to do: passer ls_epsi, aer, ls_ftsw, ls_transp, Npc dans invar!


def sol_dicout_endsim(S, outvar, DOYdeb, DOYend, opt_residu=0):
    """ close soil balance and prepare output dictionnary - call in End()"""
    # fermeture des bilans
    S.CloseWbalance(print_=0)
    S.CloseNbalance(print_=0)
    S.CloseCbalance(print_=0)
    # sys.stdout.close() #!! -> fait planter les print en re-run
    # en faire une fonction +  ecriture des disctionnaires en Rdata?

    dicout = {}
    dicout['NRain'] = S.bilanN['cumRain']
    dicout['NIrrig'] = S.bilanN['cumIrrig']
    dicout['fertNO3'] = S.bilanN['cumfertNO3']
    dicout['fertNH4'] = S.bilanN['cumfertNH4']
    dicout['HumusNMin'] = S.bilanN['cumMinN']
    if opt_residu == 1:
        dicout['Res1'] = S.bilanN['cumNRes1']
        dicout['Res2'] = S.bilanN['cumNRes2']
        dicout['Res3'] = S.bilanN['cumNRes3']
        dicout['ResidueMinN'] = S.bilanN['cumNRes1'] + S.bilanN['cumNRes2'] + S.bilanN['cumNRes3']
        for i in range(len(S.bilanN['NminfromNres'])):
            dicout['NminfromNres' + str(i)] = S.bilanN['NminfromNres'][i]  # ajout des sorties Nmin par residu

    dicout['Lix'] = S.bilanN['cumLix']
    # dicout['N2O'] = S.bilanN['cumN2O']
    dicout['UptPlt'] = list(map(np.sum, S.bilanN['cumUptakePlt']))
    dicout['azomes'] = S.bilanN['azomes']
    # serait a formater dans module sol

    # y ajoute les elements du bilan de C plante - mais devrait sortir un fichier separe!
    dicout['PARa'] = outvar['BilanC_PARa']
    dicout['RUE'] = outvar['BilanC_RUE']
    dicout['dMStot'] = outvar['BilanCdMStot']
    dicout['dMSrac'] = outvar['BilanCdMrac_fine']
    dicout['dMSpiv'] = outvar['BilanCdMpivot']
    dicout['dMSaer'] = outvar['BilanCdMaer']
    dicout['dMSenFeuil'] = outvar['BilanCdMSenFeuil']
    dicout['dMSenTige'] = outvar['BilanCdMSenTige']
    dicout['dMSenNonRec'] = outvar['BilanCdMSenNonRec']
    dicout['dMSenPiv'] = outvar['BilanCdMSenPiv']

    # WB
    dicout['cumEV'] = S.bilanW['cumEV']
    dicout['cumTransp'] = S.bilanW['cumTransp']
    dicout['cumD'] = S.bilanW['cumD']

    dicout['DOY'] = np.array(range(DOYdeb, DOYend))

    return dicout


def write_vgl_outf(outf, path_out, ls_outf_names, ls_objw, ls_keyvar_pot, outfvar):
    """ write output csv files - call in End() """

    outvarfile, outBilanNfile, outHRfile, resrootfile, lsorgfile, outMngfile, outsdfile = ls_outf_names
    outvar, dicout, out_HR, res_root, savelsOrgans, mng, res_sd = ls_objw
    ls_fileOUT = []  # liste des fichiers ecrits en sortie et affichee en fin de simul

    # ecriture variables journaliere dtoto
    if outf['outvarfile'] != 0.:
        outvarpath = os.path.join(path_out, outvarfile)  # r'H:\devel\grassland\grassland\L-gume\toto.csv'

        # valide cles active via fichier d'entree
        ls_keyvar = ['colnames']
        for i in range(1, len(ls_keyvar_pot)):
            k = ls_keyvar_pot[i]
            if outfvar[k] != 0.:
                ls_keyvar.append(k)

        # ecrit fichier avec variables selectionnees
        IOtable.write_dicttables(outvarpath, outvar, ls_keyvar)
        ls_fileOUT.append(outvarpath)

    # ecriture bilan sol
    if outf['outBilanNfile'] != 0.:
        IOtable.write_dict(dicout, path_out, outBilanNfile)
        ls_fileOUT.append(os.path.join(path_out, outBilanNfile))

    # ecriture profils sol
    if outf['outHRfile'] != 0. or outf['outFTSWfile'] != 0. or outf['outNO3file'] != 0. or outf['outNH4file'] != 0.:
        outHRpath = os.path.join(path_out, outHRfile)  # r'H:\devel\grassland\grassland\L-gume\outHR.csv'
        f = open(outHRpath, 'w')  # file (outHRpath, 'w')
        IOtable.ecriture_csv(out_HR, f)  # ls_systrac[0]
        f.close()
        ls_fileOUT.append(outHRpath)

    # ecriture en sortie du profil racinaire
    if outf['resrootfile'] != 0.:
        resrootpath = os.path.join(path_out, resrootfile)  # r'H:\devel\grassland\newres.csv'
        f = open(resrootpath, 'w')  # file (resrootpath, 'w')#r'H:\devel\grassland\newres.csv'
        IOtable.ecriture_csv(res_root, f)
        f.close()
        ls_fileOUT.append(resrootpath)

    # ecriture info organes
    if outf['lsorgfile'] != 0.:
        lsorgpath = os.path.join(path_out, lsorgfile)  # r'H:\devel\grassland\grassland\L-gume\lsAxes.csv'
        f = open(lsorgpath, 'w')  # file (lsorgpath, 'w')
        IOtable.ecriture_csv_fromlist(savelsOrgans, f)  # savelsOrgans=liste des lsOrgans pour chaque iteration #lsOrgans #lsAxes #lsApex #lsApexStop #ls_systrac[0]
        f.close()
        ls_fileOUT.append(lsorgpath)

    # ecriture mng
    if outf['outMngfile'] != 0.:
        IOtable.write_dict(mng, path_out, outMngfile)
        ls_fileOUT.append(os.path.join(path_out, outMngfile))

    # ecriture parametres avec opt sd
    if outf['outsdfile'] != 0.:  # valeur de parametres par plante
        IOtable.write_dict(res_sd, path_out, outsdfile)
        ls_fileOUT.append(os.path.join(path_out, outsdfile))

    return ls_fileOUT


def distrib_residue_mat_frominvar(ls_mat_res, S, ls_roots, profres, ParamP, invar, opt_stressGel):
    """ Distribute senescing tissues in ls_mat_res - After plant senescence/per residu type - from invar and ParamP of legume model """

    dz_sol = S.dxyz[2][0]*100. #cm
    #couches2keep = min(int(profres / dz_sol) + 1, len(S.dxyz[2]))

    root_prop0 = rtd.propRootDistrib_upZ(ls_roots, depth=dz_sol, dz_sol=dz_sol) #proportion horizon de surface
    root_propProf = rtd.propRootDistrib_upZ(ls_roots, depth=profres, dz_sol=dz_sol) #proportion horizon de profHum

    for nump in range(len(invar['dMSenRoot'])):
        groupe_resid = int(ParamP[nump]['groupe_resid'])

        ## senescence feuilles
        mat_res = root_prop0[nump] * invar['dMSenFeuil'][nump] #en surface
        ls_mat_res[groupe_resid * 4 + 0] += mat_res #ajout au groupe 1 = feuille

        ## senescence tiges
        mat_res = root_prop0[nump] * invar['dMSenTige'][nump]  # en surface
        ls_mat_res[groupe_resid * 4 + 1] += mat_res  # ajout au groupe 2 = tiges

        # turnover aerien non rec
        mat_res = root_prop0[nump] * invar['dMSenNonRec'][nump] #en surface
        ls_mat_res[groupe_resid * 4 + 0] += mat_res #ajout au groupe 1 = feuille

        ## senescence racines
        mat_res = root_propProf[nump] * invar['dMSenRoot'][nump] #sur profhum
        ls_mat_res[groupe_resid * 4 + 2] += mat_res #ajout au groupe 3 = racine fines

        ## senescence pivot
        mat_res = root_propProf[nump] * invar['dMSenPiv'][nump]  # sur profhum
        ls_mat_res[groupe_resid * 4 + 3] += mat_res  # ajout au groupe 4 = pivot

        #si gel, ajout mort gel
        if sum(invar['isGelDam']) != 0 and opt_stressGel == 1:
            mat_res = root_prop0[nump] * invar['dMSmortGel_aer'][nump]  # en surface
            ls_mat_res[groupe_resid * 4 + 0] += mat_res  # ajout au groupe 1 = feuille

        #si nouvelle plante morte
        if sum(invar['dMSmortPlant_aer']) > 0.:
            #parties aeriennes
            mat_res = root_prop0[nump] * invar['dMSmortPlant_aer'][nump]  # en surface
            ls_mat_res[groupe_resid * 4 + 0] += mat_res  # ajout au groupe 1 = feuille

            # parties racines
            mat_res = root_propProf[nump] * invar['dMSmortPlant_racfine'][nump]  # sur profhum
            ls_mat_res[groupe_resid * 4 + 2] += mat_res  # ajout au groupe 3 = racines

            # parties pivot
            mat_res = root_propProf[nump] * invar['dMSmortPlant_pivot'][nump]  # sur profhum
            ls_mat_res[groupe_resid * 4 + 3] += mat_res  # ajout au groupe 4 = feuille

    return ls_mat_res


def merge_residue_mat(ls_mat_res, vCC, S):
    # ajout dans la matrice des residus au sol

    if sum(list(map(np.sum, ls_mat_res))) > 0.:  # si de nouveaux residus (ou supeieur a un seuil
        for i in range(len(ls_mat_res)):
            mat_res = ls_mat_res[i]
            if np.sum(mat_res) > 0.:
                S.mixResMat(mat_res, i, vCC[i])

    return S


def update_residue_mat(ls_mat_res, vCC, S, ls_roots, profres, ParamP, invar, opt_residu, opt_stressGel):
    """ Distribute senescing tissues in ls_mat_res - After plant senescence/per residu type """
    # ajout dans la matrice des residus

    # dz_sol = S.dxyz[2][0]*100. #cm
    # #couches2keep = min(int(profres / dz_sol) + 1, len(S.dxyz[2]))
    #
    # root_prop0 = rtd.propRootDistrib_upZ(ls_roots, depth=dz_sol, dz_sol=dz_sol) #proportion horizon de surface
    # root_propProf = rtd.propRootDistrib_upZ(ls_roots, depth=profres, dz_sol=dz_sol) #proportion horizon de profHum
    #
    # for nump in range(len(invar['dMSenRoot'])):
    #     groupe_resid = int(ParamP[nump]['groupe_resid'])
    #
    #     ## senescence feuilles
    #     mat_res = root_prop0[nump] * invar['dMSenFeuil'][nump] #en surface
    #     ls_mat_res[groupe_resid * 4 + 0] += mat_res #ajout au groupe 1 = feuille
    #
    #     ## senescence tiges
    #     mat_res = root_prop0[nump] * invar['dMSenTige'][nump]  # en surface
    #     ls_mat_res[groupe_resid * 4 + 1] += mat_res  # ajout au groupe 2 = tiges
    #
    #     # turnover aerien non rec
    #     mat_res = root_prop0[nump] * invar['dMSenNonRec'][nump] #en surface
    #     ls_mat_res[groupe_resid * 4 + 0] += mat_res #ajout au groupe 1 = feuille
    #
    #     ## senescence racines
    #     mat_res = root_propProf[nump] * invar['dMSenRoot'][nump] #sur profhum
    #     ls_mat_res[groupe_resid * 4 + 2] += mat_res #ajout au groupe 3 = racine fines
    #
    #     ## senescence pivot
    #     mat_res = root_propProf[nump] * invar['dMSenPiv'][nump]  # sur profhum
    #     ls_mat_res[groupe_resid * 4 + 3] += mat_res  # ajout au groupe 4 = pivot
    #
    #     #si gel, ajout mort gel
    #     if sum(invar['isGelDam']) != 0 and opt_stressGel == 1:
    #         mat_res = root_prop0[nump] * invar['dMSmortGel_aer'][nump]  # en surface
    #         ls_mat_res[groupe_resid * 4 + 0] += mat_res  # ajout au groupe 1 = feuille
    #
    #     #si nouvelle plante morte
    #     if sum(invar['dMSmortPlant_aer']) > 0.:
    #         #parties aeriennes
    #         mat_res = root_prop0[nump] * invar['dMSmortPlant_aer'][nump]  # en surface
    #         ls_mat_res[groupe_resid * 4 + 0] += mat_res  # ajout au groupe 1 = feuille
    #
    #         # parties racines
    #         mat_res = root_propProf[nump] * invar['dMSmortPlant_racfine'][nump]  # sur profhum
    #         ls_mat_res[groupe_resid * 4 + 2] += mat_res  # ajout au groupe 3 = racines
    #
    #         # parties pivot
    #         mat_res = root_propProf[nump] * invar['dMSmortPlant_pivot'][nump]  # sur profhum
    #         ls_mat_res[groupe_resid * 4 + 3] += mat_res  # ajout au groupe 4 = feuille

    ls_mat_res = distrib_residue_mat_frominvar(ls_mat_res, S, ls_roots, profres, ParamP, invar, opt_stressGel)

    # inclusion dans objet sol
    # if opt_residu == 1:  # option residu activee: mise a jour des cres
    #     if sum(list(map(sum, ls_mat_res))) > 0.:  # si de nouveaux residus (ou supeieur a un seuil
    #         for i in range(len(ls_mat_res)):
    #             mat_res = ls_mat_res[i]
    #             if sum(mat_res) > 0.:
    #                 S.mixResMat(mat_res, i, vCC[i])

    if opt_residu == 1:  # option residu activee: mise a jour des cres
        S = merge_residue_mat(ls_mat_res, vCC, S)

    # calcul senesc a faire a l'echelle des axes plutot? -> a priori pas necessaire
    return [ls_mat_res, S]








#daily loop separated in 4 sub-funtionen 4 fonctions
#disentangle plant and soil steps at first! -> to allow external calls

#invar, outvar, ls_epsi, ls_demandeN_bis, temps = daily_growth_loop(ParamP, invar, outvar, res_trans, meteo_j, nbplantes, surfsolref, ls_ftswStress, ls_NNIStress, lsApex, lsApexAll)
#S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt = step_bilanWN_sol(S, par_SN, lims_sol, surfsolref, stateEV, Uval, b_, meteo_j,  mng_j, ParamP, invar, ls_epsi, ls_systrac, ls_demandeN_bis, opt_residu)
#invar, invar_sc, outvar, I_I0profilInPlant, ls_ftswStress, ls_NNIStress = Update_stress_loop(ParamP, invar, invar_sc, temps, DOY, nbplantes, surfsolref, ls_epsi, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, ls_demandeN_bis, ls_ftswStress, lsOrgans, lsApex, start_time, cutNB, deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant, I_I0profilInPlant, NlClasses, NaClasses, NlinClasses, outvar)
#ls_mat_res, S = update_residue_mat(ls_mat_res, vCC, S, carto, lims_sol, ParamP, invar, opt_residu)


#pourquoi residu semblent pas affecter les sorties et ls_mat_res remis a zero juste pres??? (parametrage de durre de vie?)
#completer les commentaires...
#separer ecriture de outvar dans une fonction specifique??



