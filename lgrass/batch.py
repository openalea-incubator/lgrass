import os
#from openalea.lpy import Lsystem
import lgrass

from openalea.lpy import *
from generateScene import run
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
    lsys = Lsystem(lpy_filename)
    axiom = lsys.axiom

    for i in range(NSTEP):

        #lsys.TPS = TPS

        lsys.DebExp=i+1
        print(lsys.DureeExp)
        runL = run(lsys, axiom = axiom , nbstep = 1)

        #lstring = lsys.derive(axiom, NSTEP)

        s_leg = runL[1].sceneInterpretation(runL[0]).deepcopy()

        lsys = runL[1]
        axiom = runL[0]
    s_leg.save("mysave2.bgeom")
    #Viewer.display(s_leg)
    #return s_leg


if __name__ == '__main__':
    test_run2()
    #Viewer.display(test_run())