# -*- coding: latin-1 -*-
"""

    setup
    ~~~~~
    
    Setup script for installation.
    
    See README.md for installing procedure.

    :copyright: TODO.
    :license: Cecill-CeCILL-C.
    
    .. seealso:: Verdenal et al, 2008.
"""

"""
    Information about this versioned file:
        $LastChangedBy: cchambon $
        $LastChangedDate: 2017-09-14 10:38:00 +0200 (jeu., 14 sept. 2017) $
        $LastChangedRevision: 21 $
        $URL: https://subversion.renater.fr/respi-wheat/trunk/setup.py $
        $Id: setup.py 21 2017-09-14 08:38:00Z cchambon $
"""

import ez_setup
import pkg_resources

ez_setup.use_setuptools()

import sys, os
from setuptools import setup, find_packages

import lgrass

if sys.version_info < (2, 7):
    print('ERROR: lgrass requires at least Python 2.7 to run.')
    sys.exit(1)

if sys.version_info >= (3, 0):
    print('WARNING: lgrass has not been tested with Python 3.')

pkg_resources.require('numpy>=1.11.0', 'pandas>=0.18.0',  'OpenAlea.Lpy', 'OpenAlea.PlantGL', 'OpenAlea.Mtg')

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "L-grass",
    version=lgrass.__version__,
    packages = find_packages(),
    include_package_data = True,
    author = "A. Verdenal, D. Combes, V. Migault",
    author_email = "didier.combes@inra.fr, abraham.escobar-gutierrez@inra.fr",
    description = "A model of rye grass morphogenesis",
    long_description = read('README.md'),
    license = "CeCILL-C",
    keywords = "functional-structural plant model, rye-grass, morphogenesis, self regulated architecture, leaf, roots ",
    url = "https://sourcesup.renater.fr/projects/lgrass/",
    download_url = "",
)
