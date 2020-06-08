"""
Model viewer.
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

def _plot(ax,parameter,yobj,lab_opts):
    """
    plot driver

    ParVars
    ----------
    ax : matplotlib.pyplot.ax
        Axis on which to plot
    xobj : XObj
        X-object to plot against
    yobj : YObj
        Y-object to plot
    lab_opts : dict
        Dictionary with label options
    """

    # plot
    axs.plot(v_plot,F[j](v_plot),**plt_opts)

    # x-axis attributes (same for all variables)
    #   ...set_xlabel... (see below)
    axs.set_xticks(_part(mu_v,v_hat))
    axs.set_xlim([0.,v_bar])
    axs.set_xticklabels([],**lab_opts)

    # y-axis attributes
    axs.set_ylabel(ylabs[j],**lab_opts)
    axs.set_ylim([0.,yticks[j]])
    axs.set_yticks([0.,yticks[j]])
    axs.set_yticklabels([r'$0$',yticklabels[j]],**lab_opts)

    # clean-up
    axs.grid()

    return ax

class XObj:
    """x-object"""

    def __init__(self,parameter,xlab,xtick,xticklab):

        # plotting attributes
        self.xlab = xlab
        self.xtick = xtick
        self.xticklab = xticklab

    def __str__(self):
        return 'FIX ME'

class YObj:
    """y-object"""

    def __init__(self,func,ylab,ytick,yticklab):
        """
        Parameters
        ----------
        func : callable
            Scalar valued function
        """

        # plotting attributes
        self.ylab = ylab
        self.ytick = ytick
        self.yticklab = yticklab

    def __str__(self):
        return 'FIX ME'

def mplot(self,xobj,yobj,
        fig_opts = {'figh' : 11., 'figw' : 8.5},
        plt_opts = {'linewidth' : 4, 'color' : 'black'},
        lab_opts = {'fontsize' : 20},
        vertical = False
        )
    """
    (m)odel (plot)ter

    Parameters
    ----------
    xobj : list (of XObj)
        X-objects
    yobj : list (of YObj)
        Y-objects
    fig_opts : dict
        Figure height
    plt_opts : dict
        Plot options
    lab_opts : dict
        Label options
    vertical : boolean
        If True, stack plots vertically and only annotate the last x-axis,
        Default is False.
    """

    # plot options
    rc('text',usetex=True)

    if vertical: 

        if len(xobj) > 1:
            raise Exception("More than one XObj supplied for vertical plot")

        fig, axs = plt.subplots(nrows=len(yobjs),ncols=1)
        for ax in axs:
            _plot_core()

        # common xlabel
        axs[-1].set_xlabel(r'$x$',**lab_opts)
        fig.savefig(yobj.name + '.pdf')
        
    else: # one figure for each plot

        if len(xobj) is not len(yobj):
            raise Exception("Number of XObjs not equal to number of YObjs")

        for yobj in yobjs:
            fig, axs = plt.subplots()
            axs = _plot_core(axs)
            fig.set_tight_layout(True)
            fig.savefig(yobj.name + '.pdf')
