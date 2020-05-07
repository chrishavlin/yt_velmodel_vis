from setuptools import setup

setup(name='yt_velmodel_vis',
      version='0.1',
      description='support for yt visualizations of seismic velocity models',
      url='https://github.com/chrishavlin/yt_velmodel_vis',
      author='Chris Havlin',
      author_email='chris.havlin@gmail.com',
      license='MIT',
      packages=['yt_velmodel_vis'],
      package_data={'yt_velmodel_vis':
                      ['datafiles/harvard-glb-volc-shapefile/*',
                       'datafiles/cb_2018_us_state_20m/*']
                   },
      install_requires=['yt','geopandas','imageio','scipy','numpy',
                        'netCDF4','h5py'],
      zip_safe=False)
