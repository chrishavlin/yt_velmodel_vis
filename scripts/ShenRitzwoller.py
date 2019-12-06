from yt_velmodel_vis import seis_model as sm
import numpy as np
import os
import yt

# if __name__=='__main__':
fname='Shen_Ritzwoller_2016.nc'
out_dir='./output'
model=sm.netcdf(fname)
model.coordTransform('sphere2cart')

# model.generateDepthSlices('vsv',depths=[50,75,100,125],
#     outDir=out_dir,SavePrefix='ShenRitz16',data_min=3.8,data_max=5)


### VOLUME RENDERING
# data={}
# data['vsv']=model.data.variables['vsv'][:]
# data['vsv'][data['vsv']==model.data['vsv'].missing_value]=np.nan
# data['vsv'] = np.transpose(data['vsv'], (1, 2, 0))
# yaw_angle=0
# roll_angle=0
# im=0
# for rot_angle in range(-180,190,10):
#     sc_mult=1.0 # scale multiplier
#     bbox = model.cart['bbox']
#     print(bbox)
#     ds = yt.load_uniform_grid(data,data['vsv'].shape,sc_mult,bbox=bbox,nprocs=1,
#                             periodicity=(False,False,False),unit_system="mks")
#     sc = yt.create_scene(ds,'vsv')
#     sc.camera.rotate(rot_angle*np.pi/180.)
#     sc.camera.yaw(yaw_angle*np.pi/180.)
#     sc.camera.roll(roll_angle*np.pi/180.)
#
#     tf = yt.ColorTransferFunction((1, 5))
#     tf.add_layers(4., w=0.1)
#
#     source = sc.sources['source_00']
#     source.tfh.set_log(False)
#     source.set_transfer_function(tf)
#     if im<10:
#         imstr='0'+str(im)
#     else:
#         imstr=str(im)
#     ang=imstr+'_'+str(rot_angle)+'_'+str(yaw_angle)+'_'+(str(roll_angle))
#     im=im+1
#     sc.save(os.path.join(out_dir,'ShenRitz_vol/ShenRitz16_vol_'+ang+'.png'),sigma_clip=4)

# data={}
# data['vsv']=model.data.variables['vsv'][:]
# data['vsv'][data['vsv']==model.data['vsv'].missing_value]=np.nan
# data['vsv'] = np.transpose(data['vsv'], (1, 2, 0))
# yaw_angle=0
# roll_angle=0
# im=0
# rot_angle=-110;
#
# sc_mult=1.0 # scale multiplier
# bbox = model.cart['bbox']
# print(bbox)
# ds = yt.load_uniform_grid(data,data['vsv'].shape,sc_mult,bbox=bbox,nprocs=1,
#                         periodicity=(False,False,False),unit_system="mks")
# sc = yt.create_scene(ds,'vsv')
# sc.camera.rotate(rot_angle*np.pi/180.)
# sc.camera.yaw(yaw_angle*np.pi/180.)
# sc.camera.roll(roll_angle*np.pi/180.)
#
# tf = yt.ColorTransferFunction((1, 5))
# tf.add_layers(2., w=1)
#
# source = sc.sources['source_00']
# # source.tfh.set_log(False)
# source.set_transfer_function(tf)
# if im<10:
#     imstr='0'+str(im)
# else:
#     imstr=str(im)
#
#
# ang=imstr+'_'+str(rot_angle)+'_'+str(yaw_angle)+'_'+(str(roll_angle))
# im=im+1
# sc.save(os.path.join(out_dir,'ShenRitz_vol/ShenRitz16_vol_'+ang+'.png'))


vsv=model.data.variables['vsv'][:]
Vmap={'vsv':vsv.data.ravel()}

lon=model.data.variables['longitude'][:].data
lat=model.data.variables['latitude'][:].data
depth=model.data.variables['depth'][:].data
[Dg,latg,long]=np.meshgrid(depth,lat,lon)
# 
# for iDepth in range(0,depth.size):
#     d_max=depth[iDepth]
#     d_min=depth[iDepth+1]
#

Dg=Dg.ravel()
latg=latg.ravel()
long=long.ravel()

# build unstructured mesh coords
Coords=[[Dg[i],latg[i],long[i]] for i in range(0,Dg.size)]

# assign elements to coords

# print(yt.load_unstructured_mesh.__doc__)

coordinates = np.array([[0.0, 0.0, 0.5], [0.0, 1.0, 0.5],[0.5, 1, 0.5], [0.5, 0.5, 0.0],[0.5, 0.5, 1.0]])
