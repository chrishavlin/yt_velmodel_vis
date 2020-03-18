# `yt_velmodel_vis` Modules

This directory contains the modules in the `yt_velmodel_vis` package.

## `seis_model`

module for loading and manipulating IRIS netcdf files. Includes interpolation routines from spherical to cartesian coordinates.

## `shapeplotter`

module for loading shapefiles and building yt line and point sources from shapefile point, line and polygon data.

## `datamanager`

simple filesystem database/cache manager. Includes an initialization class, datamanager.initializeDB() that will initialize the database, including fetching from IRIS. This now gets run via setup.py when a user installs.

## `animator`

collection of functions for building animations (sparsely populated as of now)

## `transfer_function_helper`

a copy of what was in yt's repo. Initially did this to easily see how it worked, can be deleted.


# Data

The directory `yt_velmodel_vis/datafiles` contains the shapefile data included in the repo. These get copied over to the filesystem database by  `datamanager.initializeDB()`. Other required data (from IRIS) is downloaded when the database is initialized.
