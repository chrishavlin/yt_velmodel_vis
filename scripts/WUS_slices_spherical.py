'''
loadsNWUS11-S_percent.nc and plots some slices of dVs at different depths,
saves to provided directory.
'''
import yt
from yt_velmodel_vis import seis_model as sm
import numpy as np
import os

if __name__=='__main__':
    fname='NWUS11-S_percent.nc'
    out_dir='./output'
    model=sm.netcdf(fname)

    data={}
    data['dvs']=model.data.variables['dvs'][:]
    data['dvs'][data['dvs']>900.00000]=np.nan
    # yt spherical expects R, theta, phi (depth, ~lat,~lon)
    sc_mult=1.0 # scale multiplier
    bbox = model.spherical['bbox']
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,periodicity=(False,False,False),geometry='spherical')


    c_phi=np.mean(model.spherical['bounds']['phi'])
    c_theta=np.mean(model.spherical['bounds']['theta'])

    for d_slice in [100.,150.,200.,250.,350.,400.]:
        R = (6371. - d_slice)/6371.
        cnt=[R,c_theta,c_phi]
        print("\nSlicing at "+str(d_slice))
        print(cnt)
        slc=yt.SlicePlot(ds,'r','dvs',center=cnt)
        slc.set_log("dvs", False)
        # slc.set_xlabel("latitude")
        # slc.set_ylabel("longitude")
        slc.set_cmap('dvs', 'gist_heat')
        slc.annotate_title("Depth = "+str(int(d_slice))+' km')
        slc.save(os.path.join(out_dir,'SphericalSlice_'+str(int(d_slice))))
