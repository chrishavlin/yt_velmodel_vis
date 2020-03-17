from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install

class PostDevelopCommand(develop):
    """Post-installation for development mode: initialize the filesystem db"""
    def run(self):
        develop.run(self)
        from yt_velmodel_vis import datamanager as dm
        dm.initializeDB()


class PostInstallCommand(install):
    """Post-installation for installation mode: initialize the filesystem db"""
    def run(self):
        install.run(self)
        from yt_velmodel_vis import datamanager as dm
        dm.initializeDB()

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
      install_requires=['yt','geopandas'],
      cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
      },
      zip_safe=False)
