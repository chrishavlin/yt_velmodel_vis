'''
loadsNWUS11-S_percent.nc and plots some slices of dVs at different depths,
saves to provided directory.
'''
import yt
from yt_velmodel_vis import seis_model as sm
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

#if __name__=='__main__':
fname='NWUS11-S_percent.nc'
out_dir='./output'
model=sm.netcdf(fname)
model.coordTransform('sphere2cart')
data={}
data['dvs']=model.data.variables['dvs'][:]
data['dvs'][data['dvs']>900.00000]=0
data['dvs'][data['dvs']>=0]=0
data['dvs']=data['dvs']*-1

sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox']
ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
    periodicity=(True,True,True),unit_system="mks")


center=np.array(ds.domain_center)
dm={'left':np.array(ds.domain_left_edge),
    'right':np.array(ds.domain_right_edge)}


dm['left']=dm['left']+200.*1000
dm['right']=dm['right']-200*1000
boxregion=ds.region(ds.domain_center,dm['left'],dm['right'])


for isoval in [2,2.5,3,3.5]:
    surface = ds.surface(boxregion, "dvs", isoval)
    p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.add_collection(p3dc)
    max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
    centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
    bounds = np.zeros([3,2])
    bounds[:,0] = centers[:] - max_extent/2
    bounds[:,1] = centers[:] + max_extent/2
    print(bounds)
    ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])
    
    plt.savefig('output/iso_dvs_'+str(isoval)+'.png')
    plt.close()
