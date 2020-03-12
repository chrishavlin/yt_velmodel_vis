import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from yt_velmodel_vis import datamanager as dm


# shpfile = gpd.datasets.get_path('naturalearth_lowres')

db=dm.filesysDB()
shpfile=db.validateFile('cb_2018_us_state_20m.shp')

df=gpd.read_file(shpfile)

pt_df=df[df.geometry.type=='MultiPolygon'].geometry.tolist()
for multi in pt_df:
     for poly in multi:
         if poly.boundary.type=='LineString':
            lons=np.array(poly.boundary.xy[0])
            lons[lons<0]=lons[lons<0]+360.
            plt.plot(lons,poly.boundary.xy[1])
         else:
            for ln in poly.boundary:
                lons=np.array(ln.xy[0])
                lons[lons<0]=lons[lons<0]+360.
                plt.plot(lons,ln.xy[1])


pt_df=df[df.geometry.type=='Polygon'].geometry.tolist()
for poly in pt_df:
    if poly.boundary.type=='LineString':
        lons=np.array(poly.boundary.xy[0])
        lons[lons<0]=lons[lons<0]+360.
        plt.plot(lons,poly.boundary.xy[1])
    else:
        for ln in poly.boundary:
            lons=np.array(ln.xy[0])
            lons[lons<0]=lons[lons<0]+360.
            plt.plot(lons,ln.xy[1])

plt.show()
