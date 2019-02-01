import os
#from openalea.lpy import Lsystem
import lgrass
#
from openalea.lpy import *
from generateScene import *
from openalea.plantgl.all import *

import pandas

li=list()

NSTEP = 20

def test_run():
    lpy_filename = os.path.join(lgrass.__path__[0], "lgrass.lpy")
    lsys = Lsystem(lpy_filename)
    axiom = lsys.axiom

    runL = run(lsys, axiom = axiom , nbstep = NSTEP)

    #lstring = lsys.derive(axiom, NSTEP)

    s_leg = runL[1].sceneInterpretation(runL[0]).deepcopy()
    s_leg.save("mysave.bgeom")
    print type(s_leg)
    #Viewer.display(s_leg)
    #return s_leg

def test_run2():
    lpy_filename = os.path.join(lgrass.__path__[0], "lgrass.lpy")
    lsys = lsystem(lpy_filename)
    print ' Lsystem'
    lstring = lsys.axiom
    #print lsys.TPS



    for i in xrange(NSTEP):
        print lstring
        print 'iter ',i

        for m in lstring:
            if m.name=='Feuille':
                print m.ParamFeuille.age
                li.append([lsys.TPS,m.ParamFeuille.id_plante,m.ParamFeuille.id_talle,m.ParamFeuille.id_rang,m.ParamFeuille.DigMS,m.ParamFeuille.NDFMS,m.ParamFeuille.DigNDF,m.ParamFeuille.biomass])
                m.ParamFeuille.age+=0

        lstring = lsys.derive(lstring,i,1)
        s_leg = lsys.sceneInterpretation(lstring)

    s_leg.save("mysave2.bgeom")
    table = pandas.DataFrame(li,columns=('TPS','Plante', 'Talle', 'Feuille', 'DigMS', 'NDFMS', 'DigNDF', 'Biomasse'))
    table.to_csv('Values.csv', sep=';',index=False)
    #Viewer.display(s_leg)
    #return s_leg


if __name__ == '__main__':
    test_run2()
    #Viewer.display(test_run())