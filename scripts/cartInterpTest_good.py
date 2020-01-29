import yt
from yt_velmodel_vis import seis_model as sm
import numpy as np
import os
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

fname='NWUS11-S_percent.nc'
out_dir='./output'
model=sm.netcdf(fname)

dvs_n=os.path.join(out_dir,fname.split('.')[0]+'_dvs_interp.npy')
if os.path.isfile(dvs_n):
    model.coordTransform('sphere2cart')
    print("Loading cartesian-interpolated data from file")
    dvs = np.load(dvs_n)
else:
    model.interp2cartesian(fields=['dvs'],res=[10000,10000,10000], input_units='m',max_dist=100000)
    dvs=model.interp['data']['dvs']
    np.save(dvs_n,dvs)

data={}
data['dvs']=dvs
bnds=[-6,10.]
tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
# # [center location, peak width, (red, green, blue, alpha)]
TF_gausians=[[-2,1,(1.,0.0,0.0,0.8)],
             [2,1,(0.,0.0,1.0,0.8)]]
for gau in TF_gausians:
    tf.add_gaussian(gau[0],gau[1],gau[2])

# plot the tf with a histogram
x = np.linspace(bnds[0],bnds[1],tf.nbins)
y = tf.funcs[3].y
w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                   tf.funcs[3].y]).T
fig = plt.figure(figsize=[6, 3])
ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
d_hist=ax.hist(data['dvs'].ravel(),bins=100,density=True,log=False)
ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
       log=False, color=colors, bottom=[0])
plt.savefig(os.path.join(out_dir,'cartInterpTest_tf.png'))


#
sc_mult=1.0 # scale multiplier
bbox = model.cart['bbox'] # list-like [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
print(bbox)
ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                        periodicity=(True,True,True),unit_system="mks")

sc = yt.create_scene(ds,'dvs')

# Draw the domain boundary
sc.annotate_domain(ds, color=[1, 1, 1, 0.001])

# manually draw fixed radius lines
lat_rnge=[np.min(model.data.variables['latitude']),np.max(model.data.variables['latitude'])]
lon_rnge=[np.min(model.data.variables['longitude']),np.max(model.data.variables['longitude'])]
# fixed radius, longitude vary latitude
Nsegments=100
phi=(90.-np.linspace(lat_rnge[0],lat_rnge[1],Nsegments))*np.pi / 180.
min_dep=50.
max_dep=1200.
print("adding fixed radius lines")
for depth in [min_dep,max_dep]:
    radvals=(6371.0-depth)*1000.*np.ones(phi.shape) # meters
    if depth==min_dep:
        clrval=[1.,1.,1.,0.05]
    else:
        clrval=[0.,1.,0.,0.05]

    for lon_value in lon_rnge:
        # crds['lon'][crds['lon']<0.0]=crds['lon'][crds['lon']<0.0]+360.
        if lon_value < 0:
            lon_value = lon_value + 360.
        lonval=(lon_value)*np.pi/180*np.ones(phi.shape)
            # phi: angle from north in radians (0 = north pole)
            # theta: longitudinal angle in radians
            # radius: radius, in any units
        x,y,z=sm.sphere2cart(phi,lonval,radvals)
        print([x.min(),x.max(),y.min(),y.max(),z.min(),z.max()])
        clrs=[]
        segments=[]
        for i in range(0,len(x)-1):
            segment=[[x[i],y[i],z[i]],
                    [x[i+1],y[i+1],z[i+1]]]
            segments.append(segment)
            clrs.append(clrval)
        pts=np.array(segments)
        lines=yt.visualization.volume_rendering.api.LineSource(pts,np.array(clrs))
        sc.add_source(lines)

# fixed radius, vary latitude
pos=sc.camera.position
# pos[2]=np.mean(bbox[2])
# pos[]
sc.camera.set_position(pos,north_vector=np.array([0.0, 0.0, 1.0]))
# cam.position = pos
# cam.focus = np.array([np.mean(bbox[0]),np.mean(bbox[1]),np.mean(bbox[2])])
# cam.north_vector = np.array([0.0, 0.0, 1.0])

#
source = sc.sources['source_00']
source.tfh.set_log(False)
# source.tfh.grey_opacity = False
source.set_transfer_function(tf) # apply the transfer function!
sc.save(os.path.join(out_dir,'cartInterpTest.png'),sigma_clip=2.)
