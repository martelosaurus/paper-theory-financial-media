# TODO: support multivariate Variables
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from copy import deepcopy

rc('font',size=14)

# ==============================================================================
# MODEL FUNCTIONS
# ==============================================================================
# VARIABLES AND VARIABLE PLOTTING

class Variable:
    """Variable class"""

    def __init__(self,name='',symb='',vran=None,func=None,tick=None):
        """
        Parameters
        ----------
        name : str
                Variable name (e.g. "quantity")
        symb : str
                Variable symbol in TeX (e.g. "$q$")
        vran : (float, float)
                Variable's min and max values resp.
        func : callable
                Function <<<???>>>
        tick : dict<<<???>>>
                x ticks
        """		
        self.name = name
        self.symb = symb
        self.vran = vran
        self.func = func
        self.tick = tick

    def of(self,var):
        """
        Parameters
        ----------
        var : Variable
                Variable of which to evaluate 
        """
        _self = deepcopy(self)
        return _self 

def plot(x,Y,
    _linestyles=['-','--'],
    _linewidths=[1,1],
    _colors=['k','k'],
    xticks=None,
    yticks=None,
    n_plt=100
):
    """
    For each y in Y, plots y against x

    Parameters
    ----------
    x, Y : Variable, (list of) Variables
                    Variables to plot
    xticks, yticks :
            Allow user to override the ticks
    linestyles : (list of) str
            Linestyles (see plt.plot)
    colors : (list of) str
            Colors (see plt.plot)
    linewidths : (list of) float
            Linewidths (see plt.plot)
    n_plt : int
            Number of knots at which to plot the variables

    Notes
    -----
    Default behavior is to cycle through linespecs and colors
    """		

    linestyles = _linestyles.copy()
    linewidths = _linewidths.copy()
    colors = _colors.copy()

    # TODO: collect these properties into a dict

    # create x-vector
    x_plt = np.linspace(x.vran[0],x.vran[1],n_plt)

    # plot 
    for y in Y:
        y_plt = np.linspace(y.vran[0],y.vran[1],n_plt)
        c = colors.pop()
        lw = linewidths.pop()
        ls = linestyles.pop()
        plt.plot(x_plt,y_plt,
            color	=   c,
            linewidth	=   lw,
            linestyle 	=   ls
        )	

    # x-axis
    plt.xlabel(x.name + ' ' + x.symb)
    plt.xticks(**x.tick)
    plt.xlim(x.vran)

    # ylabel (or legend)
    if len(Y) > 1:
        plt.legend([y.name + " " + y.symb for y in Y])
    else:
        plt.ylabel(Y[0].name + " " + Y[0].symb)
        plt.yticks(**Y[0].tick)
        plt.ylim(Y[0].vran)

    # misc
    plt.grid()	

    # save
    fname_x = x.name.replace(' ','_')
    fname_y = '_'.join([y.name.replace(' ','_') for y in Y])
    plt.savefig('Model_Figures/' + fname_x + '_' + fname_y + '.pdf')
    plt.close()

# ==============================================================================
# EQUILIBRIUM

class Equilibrium:
    """Equilibrium class for journalism model"""

    def __init__(self,
        alpha = 1./2.,
        kappa = 2./3.,
        chi   = 1./2.,
        b_bar = 4./3.,
        c_bar = 10.,
        v_bar = 4.,
        varphi = 0.
    ):
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
        varphi : float
                Loss aversion parameter
        """

        # ----------------------------------------------------------------------
        # parameters

        # organic parameters
        self.alpha  = alpha
        self.kappa  = kappa
        self.chi    = chi
        self.b_bar  = b_bar
        self.c_bar  = c_bar
        self.v_bar  = v_bar
        self.varphi = varphi

        # synthetic parameters
        self.mu_d   = v_bar/2. 
        self.v_hat  = (1.-alpha)*b_bar
        self.K0     = (1.+2.*chi)/(2.*kappa*(1.+chi)**2.)
        self.rho	= chi/(1.+chi)
        b_bar_max   = self.mu_d/(3.*(1.-alpha))
        c_bar_min   = 16.*(1.+2.*chi)*self.mu_d**2./(9.*kappa*(1.+chi)**2.)
        
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

        @np.vectorize
        def _beta(delta):
            """<<<EQUILIBRIUM>>> bias"""
            if v < -v_hat:
                return b_hat
            elif v < 0.:
                return (1./(3.*v_hat))*delta*b_hat
            elif v < 3*v_hat:
                return (1./v_hat)*delta*b_hat
            else: 
                return b_hat

        @np.vectorize
        def _sigma_J(delta,beta):
            """journalist's report"""
            return delta+beta

        @np.vectorize
        def _pi_R(delta,beta):
            """reporting probability (with loss aversion)"""
            Delta_R = K*(delta**2.-beta**2.)-varphi*(delta>0)*delta**2.
            if Delta_R < 0.:
                return 0.
            elif Delta_R < c_bar:
                return Delta_R/c_bar
            else:
                return 1.

        @np.vectorize
        def _Lambda(delta,beta):
            """price quality"""
            var = lambda x: (kappa*(1.-rho)*sigma_u)**2.+(rho*x)**2.
            return (-(_pi_R(delta,beta)*var(beta)+
                (1.-_pi_R(delta,beta))*var(delta)))

        @np.vectorize
        def _Omega(delta,beta):
                """drift"""
                return (_pi_R(delta,beta)*(rho*beta)+
                        (1.-_pi_R(delta,beta))*(-rho*delta))

        # ======================================================================
        # VARIABLES

        # ----------------------------------------------------------------------
        # news
        self.delta = Variable(
            name = 'news',
            symb = r'$\delta$',
            vran = (0.,1.),
            func = lambda x: x,
            tick = {
                'ticks' : [-self.mu_d,self.mu_d], 
                'labels' : [r'$-\mu_{d}$',r'$+\mu_{d}$']
            }
        )

        # ----------------------------------------------------------------------
        # bias
        self.beta = Variable(
            name = 'bias',
            symb = r'$b$',
            vran = (0.,1.),
            func = lambda x : x,
            tick = {
                'ticks' : [0.,self.b_bar], 
                'labels' : [r'$0$',r'$\\overline{b}$']
            }
        )

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> bias
        self.beta_eq = deepcopy(self.beta)

        self.beta_eq.name = 'equilibrium bias'
        self.beta_eq.symb = r'$\beta(\delta)$'

        # ----------------------------------------------------------------------
        # reporting probability
        self.pi = Variable(
            name = 'reporting probability',
            symb = r'$\pi_{R}^{*}$',
            vran = (0.,1.),
            func = _pi_R,
            tick = {'ticks' : [0.,1.], 'labels' : [r'$0$',r'$1$']}
        )

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> reporting probability
        self.pi_eq = deepcopy(self.pi)

        self.pi_eq.func = lambda delta: _pi(delta,_beta(delta))
        self.pi_eq.name = 'equilibrium reporting probability'
        self.pi_eq.symb = r'$\pi_{R}^{*}(v)$'

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> manager's report
        self.s_F = Variable(
            name = 'manager''s report',
            symb = r'$s_{F}^{*}$',
            vran = (0.,1.),
            func = lambda delta: _s_F(delta,_beta(delta)),
            tick = {'ticks' : None, 'labels' : None}
        )

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> journalist's report
        self.s_J = Variable(
            name = 'journalist report',
            symb = r'$s_{J}^{*}$',
            vran = (0.,1.),
            func = lambda delta: _s_J(delta,_beta(delta)),
            tick = {'ticks' : None, 'labels' : None}
        )

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> price quality
        self.Lambda = Variable(
            name = 'price quality',
            symb = r'$\Lambda$',
            vran = (0.,1.),
            func = lambda delta: _Lambda(delta,_beta(delta)),
            tick = {'ticks' : None, 'labels' : None} 
        )

        # ----------------------------------------------------------------------
        # <<<EQUILIBRIUM>>> drift
        self.Omega = Variable(
            name = 'drift',
            symb = r'$\Omega$',
            vran = (0.,1.),
            func = lambda delta: _Omega(delta,_beta(delta)),
            tick = {'ticks' : None, 'labels' : None }
        )
