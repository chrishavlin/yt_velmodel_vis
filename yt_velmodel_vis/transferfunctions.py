import yt
import numpy as np
import matplotlib.pyplot as plt

class dv(object):
    """ transfer function generator for perturbation fields
    """
    def __init__(self,data,**kwargs):
        """Short summary.

        Parameters
        ----------
        data : ndarray
            the data to generate transfer function
        **kwargs : type
            bounds (list):
                min/max for transfer function bounds, defaults to
                [data.min(), data.max()]

        Returns
        -------
        type
            Description of returned object.

        """

        self.data=data[~np.isnan(data)]
        self.bounds=kwargs.get('bounds',[self.data.min(),self.data.max()])
        self.tf=yt.ColorTransferFunction((self.bounds[0],self.bounds[1]))
        self.dvbins= np.linspace(self.bounds[0],self.bounds[1],self.tf.nbins)
        self.cmap=colorMapHelper(self.bounds,self.tf.nbins)

        return

    def piecewiseLinear(self,dv_pts=[],**kwargs):
        """ piecewise linear transfer function

        piecewise linear function for opacity, alpha, vs dv:

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

        alpha_min=kwargs.get('alpha_min',0.)
        alpha_max=kwargs.get('alpha_max',1.)

        dv_c=(self.dvbins[0:-1]+self.dvbins[1:])/2 # bin centers
        alpha=np.full(dv_c.shape,alpha_min)

        alpha[dv_c <= dv_pts[0]] = alpha_max
        alpha[dv_c > dv_pts[3]] = alpha_max

        dv_span=dv_pts[1]-dv_pts[0]
        dal=alpha_max - alpha_min
        alpha[dv_c>dv_pts[0]] = alpha_max - dal*(dv_c[dv_c>dv_pts[0]] - dv_pts[0]) / dv_span


        alpha[dv_c>dv_pts[1]] = alpha_min

        dv_span=dv_pts[3]-dv_pts[2]
        alpha[dv_c>dv_pts[2]] = alpha_min + dal*(dv_c[dv_c>dv_pts[2]] - dv_pts[2]) / dv_span

        alpha[dv_c>dv_pts[3]] = alpha_max

        setattr(self,'alpha',alpha)
        setattr(self,'dvbins_c',dv_c)


        # add as steps (better way to do this?)
        for i_step in range(0,len(alpha)):
            al=alpha[i_step]
            dv_1=self.dvbins[i_step]
            dv_2=self.dvbins[i_step+1]
            dv_c=self.dvbins_c[i_step]
            RGBa=self.cmap.eval(dv_c,al)
            self.tf.add_step(dv_1,dv_2,RGBa)


        return

    def addHist(self,ax=None,mplot=True,**kwargs):

        if mplot:
            if ax is None:
                f=plt.figure()
                ax=plt.axes()
            ax.hist(self.data,bins=self.dvbins,**kwargs)
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

class colorMapHelper(object):
    def __init__(self,data_bounds,n_bins):
        self.bounds=data_bounds
        self.n_bins=n_bins
        return

    def eval(self,data_val,al):

        if data_val < 0:
            R=1.
            G=0.
            B=0.
        else:
            R=0.
            G=0.
            B=1.
        RGBa=(R,G,B,al)

        return RGBa
