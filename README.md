=========================
README for L-Grass
=========================

This is L-Grass model, a model of rye grass morphogenesis.

See 
Migault V. 2015. Insertion de la morphogenèse racinaire dans L grass, un modèle structure-fonction de graminées fourragères.
Verdenal A. 2009. De la simulation de la morphogénèse de l’appareil aérien du ray-grass anglais (Lolium perenne L.). Exploration d’un schéma cybernétique inspiré du concept d’auto-organisation et applications.
Verdenal A, Combes D, Escobar-Gutiérrez AJ. 2008. A study of ryegrass architecture as a self-regulated system, using functional-structural plant modelling. Functional Plant Biology 35: 911–924.

## 1. Getting Started

These instructions will get you a copy of *L-grass* up and running on your local 
machine.

### 1.1 Prerequisites

To install and use *L-grass*, you need first to install the dependencies.

*L-grass* has been tested on Windows.
 
#### 1.1.1 Install the dependencies on Windows 10 64 bit

1. Install Python  

    * go to https://www.python.org/downloads/windows/download, 
    * click on "Latest Python 2 Release [...]", 
    * download "Windows x86-64 MSI installer" and install it selecting the following options:
        * install for all users,
        * default destination directory,
        * install all subfeatures, including subfeature "Add python.exe to Path".

2. Install OpenAlea Vplants:  

    * go to http://openalea.gforge.inria.fr/dokuwiki/doku.php?id=download:windows, 
    * download `OpenAlea 1.2` and `Vplants 1.2`,
    * install both installers: 
		  
3. Install NumPy:  

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy, 
    * download `NumPy+MKL` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `NumPy+MKL` for Python 2 64 bit,
        * install `NumPy+MKL` from the downloaded wheel file.  
          For example, if you downloaded file "numpy-1.13.1+mkl-cp27-cp27m-win_amd64.whl", 
          type: `pip install "numpy-1.13.1+mkl-cp27-cp27m-win_amd64.whl"`.

4. Install Pandas  

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#pandas, 
    * download `Pandas` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Pandas` for Python 2 64 bit,
        * install `Pandas` from the downloaded wheel file.  
          For example, if you downloaded file "pandas-0.20.3-cp27-cp27m-win_amd64.whl", 
          type: `pip install "pandas-0.20.3-cp27-cp27m-win_amd64.whl"`.
          
5. Install Sphinx

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#misc, 
    * download `Sphinx` for Python 2,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Sphinx` for Python 2,
        * install `Sphinx` from the downloaded wheel file.  
          For example, if you downloaded file "Sphinx-1.6.3-py2.py3-none-any.whl", 
          type: `pip install "Sphinx-1.6.3-py2.py3-none-any.whl"`.
          
### 1.2 Installing

__Note__: We suppose you already installed the dependencies for your operating system. Otherwise follow these [instructions](prerequisites "Prerequisites").

You can install *L-grass* either in "install" or "develop" mode.

#### 1.2.1 Install *L-grass* in "install" mode

Install *L-grass* in "install" mode if you're not going to develop, edit or debug 
it, i.e. you just want to used it as third party package.

To install *L-grass* in "end-user" mode:

* open a command line interpreter,
* go to your local copy of project *L-grass*,
* run command: `python setup.py install --user`.

#### 1.2.2 Install *L-grass* in "develop" mode

Install *L-grass* in "develop" mode if you want to get *L-grass* installed and then 
be able to frequently edit the code and not have to re-install *L-grass* to have the 
changes to take effect immediately.

To install *L-grass* in "develop" mode:

* open a command line interpreter,
* go to your local copy of project *L-grass*,
* run command: `python setup.py develop --user`.

### 1.3 Running

To run a simulation example, two options:

* 1. open Lpy platform,
	 load lgrass.lpy file from lgrass folder,
	 Use Run or animate button to launch a simulation
  2. Run lgrass from a python script (see main.py in test folder for an example)

See the user guide for a step by step explanation of how to set and run model *L-grass*.

## 2. Reading the docs

To build the user and reference guides:

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the top directory of your local copy of the project,
* run this command: `python setup.py build_sphinx`,
* and direct your browser to file `doc/_build/html/index.html`.

## 3. Testing

The test allows to verify that the model implementation accurately 
represents the developer’s conceptual description of the model and its solution.

The test:

* runs the model on 200 steps,
* concatenates the outputs of the model in pandas dataframes,
* writes the outputs dataframes to CSV files,
* compares actual to expected outputs,
* raises an error if actual and expected outputs are not equal up to a given tolerance.     

To run the test :

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the directory `test` of your local copy of the project,
* and run this command: `python test_lgrass.py`.

## Built With

* [Python](http://www.python.org/), [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/), 
* [Sphinx](http://sphinx-doc.org/): building of the documentation, 

## Contributing

First, send an email to <lgrass-request@groupes.renater.fr> to be added to the project.  

Then,
 
* check for open issues or open a fresh issue to start a discussion around a
  feature idea or a bug: https://sourcesup.renater.fr/tracker/?group_id=3957.
* If you feel uncomfortable or uncertain about an issue or your changes, feel
  free to email <lgrass@groupes.renater.fr>.

## Contact

For any question, send an email to <lgrass-request@groupes.renater.fr>.

## Versioning

We use a Git repository on [SourceSup](https://sourcesup.renater.fr) for 
versioning: https://sourcesup.renater.fr/projects/lgrass/.  
If you need an access to the current development version of the model, please send 
an email to <lgrass-request@groupes.renater.fr>.
For versionning, use a git client and get git clone git+ssh://git@git.renater.fr:2222/lgrass.git. SSH will is required

## Authors

**Alban VERDENAL**, **Didier COMBES** and **Vincent MIGAULT** - see file [AUTHORS](AUTHORS) for details

## License

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details
