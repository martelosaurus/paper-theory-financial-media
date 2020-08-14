import numpy as np

# =====================================================================
# MODEL FUNCTIONS

@np.vectorize
def _beta(delta,e):
    """EQUILIBRIUM bias"""
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
    """reporting probability"""
    Delta_R = K*(delta**2.-beta**2.)
    if Delta_R < 0.:
        return 0.
    elif Delta_R < c_bar:
        return Delta_R/c_bar
    else:
        return 1.

@np.vectorize
def _price(delta,beta):
    """price"""
    return

@np.vectorize
def _Lambda(delta,beta):
    """price quality"""
    return 

@np.vectorize
def _Omega(delta,beta):
    """price correction"""
    return np.abs(delta-_price(delta,beta))

# ---------------------------------------------------------------------
# PDF of news
def _f(delta,beta):
    """

    """

# ==============================================================================
# VARIABLES AND VARIABLE PLOTTING

class Variable:
    """Variable class"""

    def __init__(self,name,ran,labels,ticks,ticklabels):
        """
        Parameters
        ----------
        ran : (float, float)
                Variable's min and max values resp.
        ticks : (list of) float
                x ticks
        ticklabs : (list of) str
                x tick labels
        """		

def plot(x,Y,linespecs=['-'],colors=['k'],linewidths=[1]):
    """
    For each y in Y, plots y against x

    Parameters
    ----------
    x, Y : Variable, (list of) Variables
            Variables to plot

    Notes
    -----
    Default behavior is to cycle through linespecs and colors
    """		
    if not len(linespecs) is len(colors) is len(linewidths) is len(Y):
        raise Exception('PC load letter')

    for y in Y:
        plt.plot(x._plt,y._plt)	
    plt.xlabel(x.name)

    if len(Y) > 1:
        plt.legend([y.name for y in Y])
    else:
        plt.ylabel(Y[0].name)

    fname_x = x.name.replace(' ','_')
    fname_y = '_'.join([y.name.replace(' ','_') for y in Y])
    plt.savefig(fname_x + '_' + fname_y + '.pdf')

# ==============================================================================
# EQUILIBRIUM

class Equilibrium:
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

        # ----------------------------------------------------------------------
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

        # ----------------------------------------------------------------------
        # equilibrium bias
        self.
        _b(v):
        name = 'equilibrium bias'
        ylab = r'bias $\b^{*}$'
        self._b = np.vectorize(_b)

        # ----------------------------------------------------------------------
        # equilibrium reporting probability
        _pi_R(v):
        name = 'equilibrium reporting probability'
        ylab = r'reporting probability $\pi_{R}^{*}$'
        self._pi_R = np.vectorize(_pi_R)

        # ----------------------------------------------------------------------
        # equilibrium report
        _sigma_J(v):
        name = 'equilibrium report'
        ylab = r'report $s_{J}^{*}$'

        # ----------------------------------------------------------------------
        # equilibrium reporting probability"""
        _pi_R(b):
        name = 'equilibrium reporting probability'
        ylab = r'reporting probability $\pi_{R}^{*}$'
        self._pi_R = np.vectorize(_pi_R)
