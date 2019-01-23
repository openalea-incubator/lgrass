import os
#from openalea.lpy import Lsystem
import lgrass

from openalea.lpy import *
from generateScene import *
from openalea.plantgl.all import *

NSTEP = 200

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

        lstring = lsys.derive(lstring,i,1)
        s_leg = lsys.sceneInterpretation(lstring)

    s_leg.save("mysave2.bgeom")
    #Viewer.display(s_leg)
    #return s_leg


if __name__ == '__main__':
    test_run2()
    #Viewer.display(test_run())