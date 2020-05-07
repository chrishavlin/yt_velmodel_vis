import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.visualization.image_writer import map_to_colors as mtoc


class TFsegment(object):
    def __init__(self,TFob,**kwargs):

        bnds=kwargs.get('bounds',None)
        if bnds is None:
            self.dvbins=TFob.dvbins
            self.dvbins_c=TFob.dvbins_c
        else:
            # find the dvbins for bounds for this part of the TF
            self.findSubBins(TFob.dvbins,TFob.dvbins_c,bnds)
        self.nbins=self.dvbins_c.size
        self.bounds=bnds

        # set color map for this TF addition
        cmap_name=kwargs.get('cmap','Blue-Red_r')
        self.cmap=colorMapHelper(cmap_name,self.dvbins_c)
        return

    def findSubBins(self,dvbins,dvbins_c,bounds):
        self.dvbins_c=dvbins_c[(dvbins_c>=bounds[0])&(dvbins_c<=bounds[1])]
        dbin=self.dvbins_c[1]-self.dvbins_c[0]
        self.dvbins=self.dvbins_c-dbin/2.
        finalval=self.dvbins[-1]+dbin
        self.dvbins=np.append(self.dvbins,finalval)
        return

class dv(object):
    """ transfer function generator for perturbation fields (dvs or dvp)
    """
    def __init__(self,data,**kwargs):
        """
        Parameters
        ----------
        data : ndarray
            the data to generate transfer function
        **kwargs : type
            bounds (list):
                min/max for transfer function bounds, defaults to
                [data.min(), data.max()]

        """

        self.data=data[~np.isnan(data)]
        self.bounds=kwargs.get('bounds',[self.data.min(),self.data.max()])
        self.tf=yt.ColorTransferFunction((self.bounds[0],self.bounds[1]))
        self.dvbins= np.linspace(self.bounds[0],self.bounds[1],self.tf.nbins)
        self.dvbins_c=(self.dvbins[0:-1]+self.dvbins[1:])/2 # bin centers
        self.histData,_=self.calcHist()

        return

    def calcHist(self,bin_edges=None,density=True):
        """calculates histogram data for current field

        Parameters
        ----------
        bin_edges : array-like
            the bin edges to use for the histogram (if None, will use full range
            of bins)
        density : boolean
            returns normalized histogram if True (default is True).

        Returns
        -------
        ndarray
            histogram values

        ndarray
            histogram bin edges

        """

        if bin_edges is None:
            bounds=self.bounds
            bin_edges=self.dvbins

        return np.histogram(self.data,bins=bin_edges,density=density)

    def piecewiseLinear(self,dv_pts=[],**kwargs):
        """ piecewise linear transfer function

        piecewise linear function for transmission coeff., alpha, vs dv:

        at dv <= dv_pts[0],
            alpha = alpha_max
        at dv_pts[0] < dv <= dv_pts[1]
            alpha = alpha_max -(alpha_max-alpha_min) / (dv_pts[1]-dv_pts[0]) * (dv -dv_pts[0])
        at dv_pts[1] < dv <= dv_pts[2]
            alpha = alpha_min
        at dv_pts[2] < dv <= dv_pts[3]
            alpha = alpha_min + (alpha_max-alpha_min) / (dv_pts[3]-dv_pts[2]) * (dv -dv_pts[2])
        at dv > dv_pts[3]
            alpha = alpha_max

        Parameters
        ----------
        dv_pts : list of floats
            length-4 list:

        **kwargs : keyword dict
            alpha_max : float, max alpha value (default 1.0)
            alpha_min : float, min alpha value (default 0.)

        Returns
        -------
        type
            Description of returned object.

        """

        if len(dv_pts)<4:
            raise ValueError("'dv_pts' must be length 4")

        # create a segment for this part of the TF, allow unique colormap
        TFseg=TFsegment(self,**kwargs)
        dv_c=TFseg.dvbins_c
        dvbins=TFseg.dvbins

        alpha_min=kwargs.get('alpha_min',0.)
        alpha_max=kwargs.get('alpha_max',1.)

        dv_span_1=dv_pts[1]-dv_pts[0]
        dv_span_2=dv_pts[3]-dv_pts[2]
        dal=alpha_max - alpha_min

        alpha=np.expand_dims(np.piecewise(dv_c,
            [dv_c <= dv_pts[0],(dv_c>dv_pts[0])&(dv_c<=dv_pts[1]),
             (dv_c>dv_pts[1])&(dv_c<=dv_pts[2]),(dv_c>dv_pts[2])&(dv_c<=dv_pts[3]),
             dv_c > dv_pts[3]],
            [alpha_max,lambda dv_c : alpha_max - dal*(dv_c - dv_pts[0]) / dv_span_1,
            alpha_min, lambda dv_c : alpha_min + dal*(dv_c - dv_pts[2]) / dv_span_2,
            alpha_max]), axis=1)

        RGB=np.append(TFseg.cmap.RGB,alpha,axis=1)
        self.tf=addTFAsSteps(self.tf,RGB,dvbins)

        return

    def addTFsegment(self,alpha,TFseg=None,**kwargs):

        if TFseg is None:
            dv_c=self.dvbins_c
            dvbins=self.dvbins
            cmap=colorMapHelper(kwargs.get('cmap','Blue-Red_r'),dv_c)
        else:
            dv_c=TFseg.dvbins_c
            dvbins=TFseg.dvbins
            cmap=TFseg.cmap

        if len(dv_c)!=len(alpha):
            raise ValueError('length of alpha must equal the length of segment bins')

        if alpha.ndim == 1:
            alpha=np.expand_dims(alpha,1)

        RGB=np.append(cmap.RGB,alpha,axis=1)
        self.tf=addTFAsSteps(self.tf,RGB,dvbins)
        return



    def addHist(self,ax=None,mplot=True,**kwargs):

        if mplot:
            if ax is None:
                f=plt.figure()
                ax=plt.axes()

            bins=self.dvbins # the bins used by TF
            bin_size=bins[1]-bins[0]
            N_bins=int((self.data.max()-self.data.min())/bin_size)
            ax.hist(self.data,bins=N_bins,**kwargs)
            return ax

    def addTFtoPlot(self,ax=None,mplot=True,**kwargs):

        x=self.dvbins
        w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
        tf=self.tf
        colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                           tf.funcs[3].y]).T
        if mplot:
            if ax is None:
                f=plt.figure()
                ax=plt.axes()

            ax.plot(x,tf.funcs[3].y)
            ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
                   log=False, color=colors, bottom=[0])
        return ax

def addTFAsSteps(tf,RGBa,bin_edges):
    """adds a general RGBa array to transfer function as step functions.

    Parameters
    ----------
    tf : yt ColorTransferFunction
        the transfer function, already initialized
    RGBa : ndarray, shape (4,N)
        The RGBa values for each step, values between 0,1
    bin_edges : ndarray, shape (1,N+1)
        The bin edges for each RGBa value

    Returns
    -------
    tf
        modified transfer function

    """

    for i_step in range(0,len(RGBa)):
        bin_start=bin_edges[i_step]
        bin_stop=bin_edges[i_step+1]
        tf.add_step(bin_start,bin_stop,RGBa[i_step])

    return tf

class colorMapHelper(object):
    def __init__(self,cmap_name,data_vals,**kwargs):
        self.cmap=cmap_name

        # store just the RGB values of this colormap
        if cmap_name=='step_RtoB':
            # constant color value above/below center point
            center_point=kwargs.get('step_center',0)
            R=np.zeros(data_vals.shape)
            R[data_vals<center_point]=255

            B=np.zeros(data_vals.shape)
            B[data_vals>=center_point]=255

            G=np.zeros(data_vals.shape)
            self.RGB=np.dstack((R,G,B))[0]
        else:
            RGB=mtoc(np.linspace(0,1,data_vals.size),self.cmap)
            self.RGB=RGB[0,:,0:3]

        self.RGB=self.RGB/255.
        return

    # def eval(self,data_val,alpha_val=255,normalized=True):
    #
    #     RGB = None
    #     if normalized:
    #         denom=255.
    #     else:
    #         denom=1.0
    #
    #     try:
    #         RGB=self.lut[data_val]
    #         return np.append(RGB,alpha_val)/denom
    #     except KeyError:
    #         if data_val < self.data_bounds[0]:
    #             RGB=self.lut[self.data_bounds[0]]
    #             return np.append(RGB,alpha_val)/denom
    #         elif data_val > self.data_bounds[1]:
    #             RGB=self.lut[self.data_bounds[1]]
    #             return np.append(RGB,alpha_val)/denom
    #
    #     if RGB is None:
    #         raise KeyError(("data value "+str(data_val)+" is not defined for ",
    #                         "the present colormap."))
    #     return
    #
    #
    #     a=mtoc(data_vals,self.cmap)
    #     if data_val < 0:
    #         R=1.
    #         G=0.
    #         B=0.
    #     else:
    #         R=0.
    #         G=0.
    #         B=1.
    #     RGBa=(R,G,B,al)
    #
    #     return RGBa
