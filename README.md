# yt_velmodel_vis
yt visualizations of seismic velocity models


# install
1. install yt
2. install netcdf4
3. use pip to install this package:

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

# environment variables
There are a few environment variables that can be useful:

`YTVELMODELDIR` : top level velocity model directory. yt_velmodel_vis classes will recursively search this dir for model files. If not using this environment variable, simply provide the full file path when loading velocity models. In bash:
```
export YTVELMODELDIR=/path/to/yt_data/velocity_models
```

# Directory Structure
`yt_velmodel_vis` : primary package directory.

`notebooks` : python notebooks with examples

`scripts` : scripts with examples
