from yt_velmodel_vis import shapeplotter as sp
from yt_velmodel_vis import seis_model as sm
import matplotlib.pyplot as plt
import numpy as np


# fixed lat rings
fig=plt.figure()
lons=np.linspace(0,359,100)
for lat in [-80,-20,-40,0,20,40,80]:
    Tr=sp.shapeTrace(np.full(lons.shape,lat),lons)
    Tr.mplot(fig,c='k',ptype='plot')

# fixed lon rings
lats=np.linspace(-90,90,100)
lons=np.linspace(0,359,10)
for lon in lons:
    Tr=sp.shapeTrace(lats,np.full(lats.shape,lon))
    Tr.mplot(fig,c='k',ptype='plot')


# define extents of spherical chunk
lat1=10
lon1=10
lat2=50
lon2=100
rmax=5500
rmin=3000


rads=np.linspace(0,6371,20)
Trs=[]
Trs.append(sp.shapeTrace(np.full(rads.shape,lat1),np.full(rads.shape,lon1),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat1),np.full(rads.shape,lon2),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat2),np.full(rads.shape,lon1),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat2),np.full(rads.shape,lon2),rads))
for Tr in Trs:
    Tr.mplot(fig,c='g',ptype='plot')

lats=np.linspace(lat1,lat2,20)
Trs=[]
Trs.append(sp.shapeTrace(lats,np.full(lats.shape,lon1),rmax))
Trs.append(sp.shapeTrace(lats,np.full(lats.shape,lon2),rmax))
Trs.append(sp.shapeTrace(lats,np.full(lats.shape,lon1),rmin))
Trs.append(sp.shapeTrace(lats,np.full(lats.shape,lon2),rmin))

lons=np.linspace(lon1,lon2,20)
Trs.append(sp.shapeTrace(np.full(lons.shape,lat1),lons,rmax))
Trs.append(sp.shapeTrace(np.full(lons.shape,lat2),lons,rmax))
Trs.append(sp.shapeTrace(np.full(lons.shape,lat1),lons,rmin))
Trs.append(sp.shapeTrace(np.full(lons.shape,lat2),lons,rmin))

rads=np.linspace(rmin,rmax,20)
Trs.append(sp.shapeTrace(np.full(rads.shape,lat1),np.full(rads.shape,lon1),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat1),np.full(rads.shape,lon2),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat2),np.full(rads.shape,lon1),rads))
Trs.append(sp.shapeTrace(np.full(rads.shape,lat2),np.full(rads.shape,lon2),rads))

for Tr in Trs:
    Tr.mplot(fig,c='r',ptype='plot')


plt.show()
