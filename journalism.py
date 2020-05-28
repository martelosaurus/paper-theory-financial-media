import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

class Parameter:
    """Parameter class"""

    def __init__(self,parameter,xlab,xtick,xticklab):

        # plotting attributes
        self.xlab = xlab
        self.xtick = xtick
        self.xticklab = xticklab

class ModelObject:
    """model object (e.g. policy, value, etc.)"""

    def __init__(self,func,ylab,ytick,yticklab):
        """
        Parameters
        ----------
        func : callable
            Scalar valued function of model parameters
        """

        # plotting attributes
        self.ylab = ylab
        self.ytick = ytick
        self.yticklab = yticklab

class JournalismEquilibrium:
    """Equilibrium class for journalism model"""

    def __init__(self,alpha,kappa,chi,b_bar,c_bar,v_bar):
        """
        Parameters
        ----------
        alpha : float
            The journalist's skill
        kappa : float
            Convexity of the trading cost
        chi : float
            Mass of readers
        b_bar : float
            Highest permissible level of obfuscation
        c_bar : float
            Highest opportunity cost for the journalist
        v_bar : float
            Length of the support of firm values
        """

        # ---------------------------------------------------------------------
        # parameters

        # organic parameters
        self.alpha  = alpha
        self.kappa  = kappa
        self.chi    = chi
        self.b_bar  = b_bar
        self.c_bar  = c_bar
        self.v_bar  = v_bar

        # synthetic parameters
        self.mu_v   = v_bar/2. 
        self.v_hat  = (1.-alpha)*b_bar
        self.K0     = (1.+2.*chi)/(2.*kappa*c_bar*(1.+chi)**2.)
        b_bar_max   = self.mu_v/(3.*(1.-alpha))
        c_bar_min   = 16.*(1.+2.*chi)*self.mu_v**2./(9.*kappa*(1.+chi)**2.)
        
        # ---------------------------------------------------------------------
        # parameter checks

        if alpha > 1. or alpha < 0.:
            raise Exception('alpha > 1. or alpha < 0.')

        if kappa < 0.:
            raise Exception('kappa < 0.')

        if chi < 0.:
            raise Exception('chi < 0.')

        if b_bar > b_bar_max:
            raise Exception('b_bar > b_bar_max')

        if c_bar < c_bar_min:
            raise Exception('c_bar < c_bar_min')

        # ---------------------------------------------------------------------
        # parameters
        self.model_parameters = []

        firm_value  = Parameter(
            r'$0$',
            r'$\mu_{v}-\widehat{v}$',
            r'$\mu_{v}$',
            r'$\mu_{v}+3\widehat{v}$',
            r'$\overline{v}$'
        )

        firm_bias   = Parameter(
            r'$0$',
            r'$\\sup_{\\Delta_{R}}^{-1}(\\overline{c})$',
            r'$\\sup_{\\Delta_{R}}^{-1}(0)$'
        )

        # ---------------------------------------------------------------------
        # objects
        self.model_objects = []

        # ---------------------------------------------------------------------
        # equilibrium reporting probability"""
        @np.vectorize
        def _pi_R(b):
            """equilibrium reporting probability"""
            coef1 = (1.+2.*chi)/(2.*kappa*(1.+chi)**2.)
            Delta_R = coef1*((v-mu_v)**2.-((1.-alpha)*b)**2.
            if Delta_R < 0.:
                return 0.
            elif Delta_R < c_bar:
                return Delta_R/c_bar
            else:
                return 1.

        self.model_objects.append(Object(
            name = 'equilibrium reporting probability',
            func = _pi_R,
            ylab = r'reporting probability $\pi_{R}^{*}$',
            1.,
            r'$1$'
            ))

        # ---------------------------------------------------------------------
        # equilibrium bias
        @np.vectorize
        def _b(v):
            """equilibrium bias"""
            if v < -v_hat:
                return b_hat
            elif v < 0.:
                return (1./(3.*v_hat))*(v-mu_v)*b_hat
            elif v < 3*v_hat:
                return (1./v_hat)*(v-mu_v)*b_hat
            else: 
                return b_hat

        self.model_objects.append(Object(
            name = 'equilibrium bias',
            func = _b,
            ylab = r'bias $\b^{*}$',
            b_bar,
            r'$\overline{\b}$'
            ))

        # ---------------------------------------------------------------------
        # equilibrium reporting probability
        @np.vectorize
        def _pi_R(v):
            """equilibrium reporting probability"""
            if v < -v_hat:
                return K0*((v-mu_v)**2.-v_hat**2.)
            elif v < 0.:
                return (8./9.)*K0*(v-mu_v)**2.
            elif v < 3*v_hat:
                return 0.
            else: 
                return K0*((v-mu_v)**2.-v_hat**2.)

        self.model_objects.append(Object(
            name = 'equilibrium reporting probability',
            func = _pi_R,
            ylab = r'reporting probability $\pi_{R}^{*}$',
            1.,
            r'$1$'
            ))

        # ---------------------------------------------------------------------
        # equilibrium report
        @np.vectorize
        def _sigma_J(v):
            """equilibrium report"""
            return v+(1.-alpha)*_b(v)

        self.model_objects.append(Object(
            name = 'equilibrium report',
            func = _sigma_J,
            ylab = r'report $s_{J}^{*}$',
            v_bar+(1.-alpha)*b_bar,
            r'$\overline{v}+(1-\alpha)\overline{\b}$'
            ))

    def __str__(self):
        return 'FIX ME'

    def plot(self,
            mod_objs = [],
            f_name='demo.pdf',
            fig_opts = {'figh' : 11., 'figw' : 8.5},
            plt_opts = {'linewidth' : 4, 'color' : 'black'},
            lab_opts = {'fontsize' : 20},
            together = False,
            vertical = False
            )
        """
        Parameters
        ----------
        f_name : str
            Figure file name
        fig_opts : float
            Figure height
        mod_obs : list (of str)
            List of model object names (ModelObject.name). If None, then asks
            user to choose model objects
        one_page : boolean
            If True, then plots all of the model objects in one figure. 
            Default is False
        vertical : boolean
            If True, stack plots vertically and only annotate the last x-axis,
            Default is False.

        Notes
        -----
        If vertical = True, then it is assumed that one_page = True. 
        The converse is not necessarily true.
        """
    
        if not mod_objs:
            for mod_obj in self.model_objects:
                if input('Plot ' + mod_obj.name) + '? (y/N)') is 'y':
                    mod_objs.append(mod_obj.name)

        # plot options
        rc('text',usetex=True)

        if vertical:
            fig, axs = plt.subplots(nrows=len(mod_objs),ncols=1)

            # common xlabel
            axs[-1].set_xlabel(r'firm value $v$',**lab_opts)
        elif one_page:
            fig.set_size_inches(**fig_size)
        else:

        for mod_obj in mod_objs:

            # plot
            axs[j].plot(v_plot,F[j](v_plot),**plt_opts)

            # x-axis attributes (same for all variables)
            #   ...set_xlabel... (see below)
            axs[j].set_xticks(_part(mu_v,v_hat))
            axs[j].set_xlim([0.,v_bar])
            axs[j].set_xticklabels([],**lab_opts)

            # y-axis attributes
            axs[j].set_ylabel(ylabs[j],**lab_opts)
            axs[j].set_ylim([0.,yticks[j]])
            axs[j].set_yticks([0.,yticks[j]])
            axs[j].set_yticklabels([r'$0$',yticklabels[j]],**lab_opts)

            # clean-up
            axs[j].grid()

        # clean-up
        fig.set_tight_layout(True)
        fig.savefig(f_name)
        fig.show()
