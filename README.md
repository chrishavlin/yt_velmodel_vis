# yt_velmodel_vis
yt visualizations of seismic velocity models

# dependencies:
the yt package (https://yt-project.org) is required.

`geopandas` is also required. It has a number of dependencies itself including `pandas` and `shapely`, so easiest
to install with pip or conda.

* For conda: `conda install geopandas`
* For pip: `pip install geopandas`

# setting up yt_velmodel_vis

After you have the above dependencies installed, you can install this package in two steps:
1.  package installation with `pip`
2.  (optional) setting environment variables and initializing the local database


## 1. package installation

Navigate to the top level of the repository and then use `pip` to install.

To install in the default pip location:
```
pip install .
```

To install within a conda environment, activate the environment first. If you are installing system wide, you will likely need `sudo` here.

To install the package in current user's home:
```
pip install -u .
```

Or to install in debug mode (changes to source code will update), while in top directory of repo:

```
pip install -e .
```

If using conda, remember to activate the appropriate environment first.


## 2. (optional) setting environment variables and building local database

The example scripts and some `yt_velmodel_vis` classes require datafiles to be stored in a local filesystem database (cache). The location of this database is controlled by the directory stored in the environment variable `YTVELMODELDIR`.

So after installing the package, set this environment variable to a directory where you'd like to install the local filesystem database. It can be an existing or new directory: new data files will be added as subdirectories within. `yt_velmodel_vis` classes will recursively search this directory for model files.

After setting the environment variable, open a python shell and initialize the database with

```
from yt_velmodel_vis import datamanager as dm
dm.initializeDB()
```

To add an environment variable on a unix system using a bash shell, add the following to the `.bashrc` or `.bash_aliases` file:

```
export YTVELMODELDIR=/path/to/yt_data/velocity_models
```

If you do not set `YTVELMODELDIR`, then models will be installed in `~/.ytvelmodeldata/`.

# Directory Structure
`yt_velmodel_vis` : primary package directory.

`notebooks` : python notebooks with examples

`scripts` : scripts with examples
