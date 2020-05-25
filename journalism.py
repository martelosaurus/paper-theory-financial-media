import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

# ------------------------------------------------------------------------------
# parameters
rc('text',usetex=True)
n_plot = 100
v_plot = np.linspace(0.,1.,n_plot)
plt_opts = {'linewidth' : 4, 'color' : 'black'}
lab_opts = {'fontsize' : 20}

# ------------------------------------------------------------------------------
# model functions
def check(inequality,parameters):
    #TODO: fix scope
    """
    Evaluate an inequality and raise an exception if it evaluates to False

    Parameters
    ----------
    inequality : str
        String representation of a python logical statement. 
    parameters : dict
        Dictionary of variables

    Examples
    --------
    >>> x = 1, y = 2
    >>> check('x>y')
    
    will raise the exception.
    """
    try:
        _eval = eval(inequality,parameters)
    except NameError:
        print(inequality + ' contains an undefined parameter.')
    if not _eval:
        raise Exception(inequality + ' does not hold. Fix your parameters.')

class _plot_attr:

    def __init__(self,ylab,ytick,yticklab):

        # plotting attributes
        self.ylab = ylab
        self.ytick = ytick
        self.yticklab = yticklab

# ------------------------------------------------------------------------------
# equilibrium 
class Equilibrium:

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
            Length of the support of firm value 
        """

        # organic
        self.alpha  = alpha
        self.kappa  = kappa
        self.chi    = chi
        self.b_bar  = b_bar
        self.c_bar  = c_bar
        self.v_bar  = v_bar

        # synthetic
        self.mu_v   = v_bar/2. 
        self.v_hat  = (1.-alpha)*b_bar
        self.K0     = (1.+2.*chi)/(2.*kappa*c_bar*(1.+chi)**2.)
        b_bar_max   = self.mu_v/(3.*(1.-alpha))
        c_bar_min   = 16.*(1.+2.*chi)*self.mu_v**2./(9.*kappa*(1.+chi)**2.)
        
        bias_plot = _plot_attr(r'bias $\b^{*}$',b_bar,r'$\overline{\b}$')
        prob_plot = _plot_attr(r'reporting probability $\pi_{R}^{*}$',1.,r'$1$')
        rprt_plot = _plot_attr(
                r'report $s_{J}^{*}$',
                v_bar+(1.-alpha)*b_bar,
                r'$\overline{v}+(1-\alpha)\overline{\b}$')

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

        # check parameter inequalities
        # for inequality in inequalities:
            # check(inequality)

        @np.vectorize
        def _part(v):
            """Partition: returns Boolean list"""
            return np.array([-mu_v,-v_hat,0.,3*v_hat,mu_v])

        # equilibrium bias
        self._b():
                [
                    lambda v : b_hat, 
                    lambda v : (1./(3.*v_hat))*(v-mu_v)*b_hat,
                    lambda v : (1./v_hat)*(v-mu_v)*b_hat,
                    lambda v : b_hat
                ]
            
        # equilibrium reporting probability
        @np.vectorize
        self._pi_R:
             [
                    lambda v : K0*((v-mu_v)**2.-v_hat**2.),
                    lambda v : (8./9.)*K0*(v-mu_v)**2.,
                    lambda v : 0.,
                    lambda v : K0*((v-mu_v)**2.-v_hat**2.)
                ]

        @np.vectorize
        def _sigma_J(alpha):
            """Equilibrium report (conditional on reporting)"""
            return v+(1.-alpha)*_b(v)

    def plot(self,n_plot=100,f_name='figure4.pdf',figh=11.,figw=8.5):
        """
        Plot the equilibrium quantities

        Parameters
        ----------
        n_plot : int
            Number of plotting nodes
        f_name : str
            Figure file name
        figh : float
            Figure height
        figw : float
            Figure width
        """
        
        # parameters
        v_plot = np.linspace(0.,self.v_bar,n_plot)

        # figures
        fig, axs = plt.subplots(nrows=3)
        fig.set_figheight(figh)
        fig.set_figwidth(figw)
        fig.set_tight_layout(True)

        # common axis attributes
        F = [_b, _pi_R, _sigma_J] 
        for j in [0,1,2]:

            # plot
            axs[j].plot(v_plot,F[j](v_plot),**plt_opts)

            # x-axis attributes (same for all variables)
            #   ...set_xlabel... (see below)
            axs[j].set_xticks(_part(mu_v,v_hat))
            axs[j].set_xlim([0.,v_bar])
            axs[j].set_xticklabels([
                    r'$0$',
                    r'$\mu_{v}-\widehat{v}$',
                    r'$\mu_{v}$',
                    r'$\mu_{v}+3\widehat{v}$',
                    r'$\overline{v}$'
                ],**lab_opts)

            # y-axis attributes
            axs[j].set_ylabel(ylabs[j],**lab_opts)
            axs[j].set_ylim([0.,yticks[j]])
            axs[j].set_yticks([0.,yticks[j]])
            axs[j].set_yticklabels([r'$0$',yticklabels[j]],**lab_opts)

            # clean-up
            axs[j].grid()

        # common xlabel
        axs[-1].set_xlabel(r'firm value $v$',**lab_opts)

        # clean-up
        fig.savefig(f_name)
        fig.show()
