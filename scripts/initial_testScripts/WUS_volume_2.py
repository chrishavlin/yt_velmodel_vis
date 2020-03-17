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

if __name__=='__main__' or gogo:
    fname='NWUS11-S_percent.nc'
    out_dir='./output'
    model=sm.netcdf(fname)


    model.coordTransform('sphere2cart')
    data={}
    # hvals={'depth':[20.,40.]} # values to put on boundaries
    # data['dvs']=model.addBoundaryHighlights('dvs',hvals,newFill=0.)
    data['dvs']=np.array(model.data.variables['dvs'][:])
    # data['dvs'][data['dvs']>100.00000]=0.01
    # data['dvs']=data['dvs']*-1
    # data['dvs'][data['dvs']<=0.]=0.01

    # depth = 0
    data['dvs'][0,:,:]=20.
    data['dvs'][1,:,:]=20.
    data['dvs'][2,:,:]=20.

    # depth = max
    # data['dvs'][-2:,:,:]=40.
    # data['dvs'][-1:,:,:]=40.

    #
    # data['dvs'][:,0,:]=20.
    # data['dvs'][:,1,:]=20.
    # data['dvs'][:,2,:]=20.

    # data['dvs'] = np.transpose(data['dvs'], (1, 2, 0))
    print(data['dvs'].max())
    print(data['dvs'].min())

    # build the transfer function
    bnds=[data['dvs'].min(),60.]
    tf = yt.ColorTransferFunction((bnds[0],bnds[1]))
    # tf.add_layers(5, w=0.05)

    # add the boundary value highlights
    # tf.add_gaussian(40.,2.,[0.,1.0,0.0,.2])
    tf.add_gaussian(20.,2.5,[1.0,0.0,0.0,0.1])
    # tf.add_gaussian(40.,5.,[0.,1.0,0.,.5])

    # build list of gaussians to add to tf
    # [center location, peak width, (red, green, blue, alpha)]
    # TF_gausians=[[1.5,.25,(0.,0.0,1.0,0.6)],
    #              [-1.5,.25,(1.,0.5,0.0,0.5)],
    #              [-4,.25,(1.,0.0,0.0,0.8)]]
    TF_gausians=[[-6,5,(1.,0.0,0.0,0.8)],
                 [4,3,(0.,0.0,1.0,0.5)],
                 [5,0.5,(0.,0.0,1.0,0.9)]]
    TF_gausians=[[2,3,(0.,0.0,1.0,0.5)]]
    # TF_gausians=[[-1,.5,(1.,0.0,0.0,0.6)]]
    # TF_gausians=[[1,.5,(0.,0.0,1.0,0.9)]]
    # for gau in TF_gausians:
    #     tf.add_gaussian(gau[0],gau[1],gau[2])
    print("hmm")

    print(tf)

    # plot the tf with a histogram
    x = np.linspace(bnds[0],bnds[1],tf.nbins)
    y = tf.funcs[3].y
    w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
    colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                       tf.funcs[3].y]).T
    print("yup")

    fig = plt.figure(figsize=[6, 3])
    ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
    d_hist=ax.hist(data['dvs'].ravel(),bins=100,density=True,log=False)

    print("hello")

    ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
           log=False, color=colors, bottom=[0])
    print("hi")
    plt.savefig(os.path.join(out_dir,'WUS_volume_2_tf.png'))
    print("fig saved")


    # yt spherical expects R, theta, phi (depth, ~lat,~lon)
    sc_mult=1.0 # scale multiplier
    bbox = model.cart['bbox']
    print(bbox)
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
                            periodicity=(True,True,True),unit_system="mks")

    sc = yt.create_scene(ds,'dvs')
    # sc.annotate_axes()
    # cam = sc.add_camera()


    # Choose a vector representing the viewing direction.
    # L = [0.5, 0.5, 0.5]
    # # Define the center of the camera to be the domain center
    # c = ds.domain_center[0]
    # # Define the width of the image
    # W = 1.5*ds.domain_width[0]
    # # Define the number of pixels to render
    # Npixels = 512
    # cam = sc.add_camera()
    # cam.width = W
    # cam.center = c
    # cam.normal_vector = L
    # cam.north_vector = [0, 0, 1]

    pos=sc.camera.position
    # pos[2]=np.mean(bbox[2])
    # pos[]
    sc.camera.set_position(pos,north_vector=np.array([0.0, 0.0, 1.0]))
    # cam.position = pos
    # cam.focus = np.array([np.mean(bbox[0]),np.mean(bbox[1]),np.mean(bbox[2])])
    # cam.north_vector = np.array([0.0, 0.0, 1.0])

    # sc.camera.zoom(.3)
    # wid=np.array([bbox[0][1]-bbox[0][0],bbox[1][1]-bbox[1][0],bbox[2][1]-bbox[2][0]])
    # wid=wid*0.1
    # wid_q=[ds.quan(wid[i],'m') for i in range(0,3)]
    # sc.camera.set_width(wid_q)#ds.quan(50*1000., 'm'))
    # sc.camera.north_vector = np.array([0.0, 0.0, 1.0])
    # sc.camera.rotate(10*np.pi/180.)
    # sc.camera.roll(90.*np.pi/180.)
    # sc.camera.pitch(np.pi/4)
    # sc.camera.yaw(45.*np.pi/180)
    # cam = sc.add_camera()
    # cam.resolution=1000

    source = sc.sources['source_00']
    source.tfh.set_log(False)
    # source.tfh.grey_opacity = False
    source.set_transfer_function(tf) # apply the transfer function!
    sc.save(os.path.join(out_dir,'WUS_volume_2.png'),sigma_clip=2.)

    # bbox = model.cart['bbox']
    # ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,
    #     periodicity=(True,True,True),unit_system="mks")
    #
    #
    # center=np.array(ds.domain_center)
    # dm={'left':np.array(ds.domain_left_edge),
    #     'right':np.array(ds.domain_right_edge)}
    #
    #
    # dm['left']=dm['left']+200.*1000
    # dm['right']=dm['right']-200*1000
    # boxregion=ds.region(ds.domain_center,dm['left'],dm['right'])
    #
    # isoval=1
    # surface = ds.surface(boxregion, "dvs", isoval)
    # p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
    #
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.add_collection(p3dc)
    # max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
    # centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
    # bounds = np.zeros([3,2])
    # bounds[:,0] = centers[:] - max_extent/2
    # bounds[:,1] = centers[:] + max_extent/2
    # print(bounds)
    # ax.auto_scale_xyz(bounds[0,:], bounds[1,:], bounds[2,:])
    #
    # plt.savefig('output/WUS_isosurf_'+str(isoval)+'.png')
    # plt.close('all')
