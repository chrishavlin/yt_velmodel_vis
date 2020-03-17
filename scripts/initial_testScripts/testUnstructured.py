from yt_velmodel_vis import seis_model as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import yt,os
import numpy as np

fname='NWUS11-S_percent.nc'
model=sm.netcdf(fname)
model.moveToUnstructured(fields=['dvs'],coordfields=['depth'])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x=[]
y=[]
z=[]
cvar=[]
for ele in range(0,200000,500):#[0,10000,200000,50,10000]:
    # print(ele)
    for vert in model.unstructured['cnnct'][ele]:
        x.append(model.unstructured['vertices'][vert,0])
        y.append(model.unstructured['vertices'][vert,1])
        z.append(model.unstructured['vertices'][vert,2])

    # node-depths
    for cval in  model.unstructured['data']['connect1', 'depth'][ele]:
        cvar.append(cval)

# x = np.array(x)
# y = np.array(y)
# z = np.array(z)
# cvar =np.sqrt(x**2+y**2+z**2)

ax.scatter(x,y,z,c=cvar,marker='.')
ax.set_xlim([min(x),0])
ax.set_ylim([min(y),0])
ax.set_zlim([max(z),0])
plt.show()


# https://yt-project.org/doc/examining/loading_data.html#unstructured-grid-data
# model.unstructured['data']['connect1', 'dvs']=np.float64(model.unstructured['data']['connect1', 'dvs'])
# model.unstructured['cnnct']=np.float64(model.unstructured['cnnct'])

ds = yt.load_unstructured_mesh(model.unstructured['cnnct'], model.unstructured['vertices'],node_data=model.unstructured['data'])

sc = yt.create_scene(ds,('connect1', 'depth'))
# # override the default colormap
# ms = sc.get_source()
# ms.cmap = 'Eos A'
#
# # adjust the camera position and orientation
# cam = sc.camera
# cam.focus = ds.arr([0.0, 0.0, 0.0], 'code_length')
# cam_pos = ds.arr([-3.0, 3.0, -3.0], 'code_length')
# north_vector = ds.arr([0.0, -1.0, -1.0], 'dimensionless')
# cam.set_position(cam_pos, north_vector)
#
# # increase the default resolution
# cam.resolution = (800, 800)

# render and save
# sc.save('output/testUnstructured.png')
# tf = yt.ColorTransferFunction((0,5))
# tf.add_gaussian(3.,1,[1.0,0.0,0.0,0.1])
# source = sc.sources['source_00']
# # source.tfh.set_log(False)
# source.set_transfer_function(tf) # apply the transfer function!



ms = sc.get_source()
ms.cmap = 'Eos A'
sc.save(os.path.join('output','testUnstructured.png'))#,sigma_clip=2.)
