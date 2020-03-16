# yt_velmodel_vis
yt visualizations of seismic velocity models

# dependencies:
the yt package (https://yt-project.org) is required.

`geopandas` is also required. It has a number of dependencies itself including `pandas` and `shapely`, so easiest
to install with pip or conda.

* For conda: `conda install geopandas`
* For pip: `pip install geopandas`

# setting up yt_velmodel_vis


## 1. package installation

After you have the above dependencies installed, you can install this package using pip:

#### general system install
installs in default pip location, may require sudo if installing system-wide
```
pip install .
```
to install within a conda environment, activate the environment first.

#### user-specific install
installs the package in current user's home
```
pip install -u .
```

#### in debug mode
to install in debug mode (changes to source code will update), while in top directory of repo:

```
pip install -e .
```
to install within a conda environment, activate the environment first.

## 2. example data and environment variables

To run the example scripts here, you need to build a local filesystem database and fetch data from IRIS. This can be done by running the following script:

```
python scripts/dataSetup.py -top_dir /path/to/your/dir
```

You can then set the following environment variable:

`YTVELMODELDIR`

`yt_velmodel_vis` classes will recursively search this dir for model files. TO add an environment variable on a unix system using a bash shell, add the following to the `.bashrc` or `.bash_aliases` file:
```
export YTVELMODELDIR=/path/to/yt_data/velocity_models
```

# Directory Structure
`yt_velmodel_vis` : primary package directory.

`notebooks` : python notebooks with examples

`scripts` : scripts with examples
