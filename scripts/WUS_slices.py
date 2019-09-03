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

    bnds={
    'longitude':[-130,-104],
    'latitude':[32,54],
    'depth':[50,1200]
    }

    for bnd in [0,1]:
        if bnds['latitude'][bnd] > 180:
            bnds['latitude'][bnd]=bnds['latitude'][bnd]-360

    data={}
    data['dvs']=model.data.variables['dvs'][:]
    data['dvs'][data['dvs']>900.00000]=np.nan
    data['dvs'] = np.transpose(data['dvs'], (2, 1, 0))

    sc_mult=1.0 # scale multiplier
    bbox = np.array([bnds['longitude'],bnds['latitude'],bnds['depth']])
    ds = yt.load_uniform_grid(data,data['dvs'].shape,sc_mult,bbox=bbox,nprocs=1,periodicity=(False,False,False))


    clon=np.mean(bnds['longitude'])
    clat=np.mean(bnds['latitude'])
    for d_slice in [100.,150.,200.,250.]:
        cnt=[clon,clat,d_slice]
        slc=yt.SlicePlot(ds,'z','dvs',center=cnt)
        slc.set_log("dvs", False)
        slc.set_xlabel("latitude")
        slc.set_ylabel("longitude")
        slc.set_cmap('dvs', 'gist_heat')
        slc.annotate_title("Depth = "+str(int(d_slice))+' km')
        slc.save(os.path.join(out_dir,'Slice_'+str(int(d_slice))))
