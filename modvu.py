"""
Model viewer.
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

def _plot_core(ax,parameter,mod_obj,lab_opts):
    """
    plot driver

    ParVars
    ----------
    ax : matplotlib.pyplot.ax
        Axis on which to plot
    parameter : ParVar
        ParVars against which to plot the model object
    mod_obj : ModObj
        Model object to plot

    lab_opts : dict
        Dictionary with label options
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

class ParVar:
    """(Par)ameter/(Var)iable class"""

    def __init__(self,parameter,xlab,xtick,xticklab):

        # plotting attributes
        self.xlab = xlab
        self.xtick = xtick
        self.xticklab = xticklab

class ModObj:
    """(Mod)el (Obj)ect (e.g. policy, value, etc.)"""

    def __init__(self,func,ylab,ytick,yticklab):
        """
        Parameters
        ----------
        func : callable
            Scalar valued function of the ParVar
        """

        # plotting attributes
        self.ylab = ylab
        self.ytick = ytick
        self.yticklab = yticklab

    def __str__(self):
        return 'FIX ME'

    def plot(self,
            mod_obj     = [],
            fig_opts    = {'figh' : 11., 'figw' : 8.5},
            plt_opts    = {'linewidth' : 4, 'color' : 'black'},
            lab_opts    = {'fontsize' : 20},
            vertical    = False
            )
        """
        Parameters
        ----------
        fig_opts : float
            Figure height
        mod_obs : list (of str)
            List of model object names (ModObj.name). If None, then asks
            user to choose model objects
        vertical : boolean
            If True, stack plots vertically and only annotate the last x-axis,
            Default is False.
        """

        # if unspecified, ask user which model objects to plot
        if not mod_objs:
            for mod_obj in self.model_objects:
                if input('Plot ' + mod_obj.name) + '? (y/N)') is 'y':
                    mod_objs.append(mod_obj.name)

        # plot options
        rc('text',usetex=True)

        if vertical: # one figure with one column with one x-axis for all plots 

            fig, axs = plt.subplots(nrows=len(mod_objs),ncols=1)
            for ax in axs:
                _plot_core()

            # common xlabel
            axs[-1].set_xlabel(r'firm value $v$',**lab_opts)
            fig.savefig(mod_obj.name + '.pdf')
            
        else: # one figure for each plot

            for mod_obj in mod_objs:
                fig, axs = plt.subplots()
                axs = _plot_core(axs)
                fig.set_tight_layout(True)
                fig.savefig(mod_obj.name + '.pdf')
