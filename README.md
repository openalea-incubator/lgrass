[![License](https://img.shields.io/badge/license-CeCILL--C-blue )](https://img.shields.io/badge/license-CeCILL--C-blue )

[![Documentation Status](https://readthedocs.org/projects/lgrass/badge/?version=latest)](https://lgrass.readthedocs.io/en/latest/?badge=latest)

# L-Grass
L-Grass is a Functional Structural Plant Model (FSPM) of perennial rye grass morphogenesis and phenology.

L-Grass simulates:
* Above-ground 3D architecture and morphogenesis (leaf extension, growth and tillering) as a self-regulated system. See [Verdenal *et al.* (2008)](https://doi.org/10.1071/FP08050).
* Morphogenesis of the root system in 3D. The model is inspired from ArchiSimple model and accounts for C partioning between shoot and roots. See [Migault thesis (2015)] (https://hal.inrae.fr/tel-02799991)
* Floral transition of individual apex and heading of the spike according to temperature, photoperiod and vegetative morphogenesis. See [Rouet *et al.* (2022)](https://doi.org/10.1093/insilicoplants/diac012).

# Table of Contents
- [Installation](#installation)
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
    + [Users](#users)
    + [Developers](#developers)
- [Usage](#usage)
- [Credits](#credits)
  * [Authors](#authors)
  * [Contributors](#contributors)
  * [Funding](#funding)
- [License](#license)


# Installation

*Lgrass* has been tested on Windows 10 64 bit. *Lgrass* was iniatially implemented in Python 2.7 and is now under Python 3.

## Prerequisites
*Lgrass* has the following dependencies (see documentation in the links provided, instructions for their installation are given in [Installing](#installing)):
* To run the model: 
    * [Python](http://www.python.org) >= 3.7 (old releases in Python 2.7)
    * [NumPy](http://www.numpy.org/)
    * [SciPy](http://www.scipy.org/)
    * [Pandas](http://pandas.pydata.org/)
    * [Openalea.MTG](https://github.com/openalea/mtg)
    * [Openalea.Plantgl](https://github.com/openalea/plantgl)
    * [Openalea.Lpy](https://github.com/openalea/lpy)
    * [Alinea.Caribu](https://github.com/openalea-incubator/caribu) 
    * [ephem](https://pypi.org/project/ephem/)
    * [geopy](https://pypi.org/project/geopy/)
    * [xlrd](https://pypi.org/project/xlrd/)
    
## Installing
*Lgrass* has to be installed in a conda environment containing all dependencies.

* Install Miniconda 2 or 3 for Python 3.7: https://docs.conda.io/en/latest/miniconda.html
* Create a new environment in an Anaconda prompt:
   `conda create -n Lgrass python=3.7 openalea.mtg openalea.plantgl openalea.lpy alinea.caribu ephem geopy xlrd -c conda-forge -c fredboudon`
   
### Users
* Download Lgrass package from https://github.com/openalea-incubator/lgrass/archive/master.zip
* Extract archive
* Open an Anaconda promt in your local copy of *Lgrass*,
* Activate the conda environment: `conda activate Lgrass`
* Run command: `python setup.py install --user` 

### Developers
* Fork Lgrass repository in your own account
* Clone your fork with git : `git clone https://github.com/my_account/lgrass.git`
* Open an Anaconda promt in your cloned copy of *Lgrass*,
* Activate the conda environment: `conda activate Lgrass`
* Run command: `python setup.py develop --user` 

#### Updating submodules

If you want to update *Lgrass*:

    git pull origin/master

# Usage

* For a basic usage of *Lgrass*, go to `lgrass/lgrass` and open `lgrass.lpy` with the Lpy software which was installed in your conda environment. Then click on `Run`
* For external call of *lgrass*, examples are provided in `lgrass/example`

The version of lgrass integrating the flowering stages (Lgrass-F) is available in: 

* `lgrass/example/calibration_GEVES_2021`

## calibration_GEVES_2021
This example deals with the implementation of the reproductive development in the model and its calibration by Rouet et coll. 
This work led to the research article [Rouet *et al.* (2022)](https://doi.org/10.1093/insilicoplants/diac012).
Results were obtained from the tag paper_ISPLANTS_2022. To run the model used for the paper, please download the code archives at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6873725.svg)](https://doi.org/10.5281/zenodo.6873725)

# Credits
## Authors
* Alban VERDENAL - model designing, development and validation
* Vincent MIGAULT - model designing, development and validation
* Simon ROUET - model designing, development and validation - [SimonRouet](https://github.com/SimonRouet)
* Abraham ESCOBAR-GUTIERREZ - model designing and validation, scientific project management
* Didier COMBES - model designing, development and validation, scientific project management - [dicombes](https://github.com/dicombes)

## Contributors
* Romain BARILLOT - [rbarillot](https://github.com/rbarillot)
* Thibault RAQUET - [traquet](https://github.com/traquet)
* Jean-Louis DURAND

## Funding
* [INRAE](https://www.inrae.fr/): salaries of permanent staff 
* [itk](https://www.itk.fr/en/) company and [ANRT](http://www.anrt.asso.fr/fr): funded the [Cifre](http://www.anrt.asso.fr/fr/cifre-7843) PhD thesis of V. Migault
* Région Nouvelle-Aquitaine: funding of PhD thesis of S. Rouet
to be continued...

# License
This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
 
