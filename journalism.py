import numpy as np
from copy import deepcopy

# ==============================================================================
# MODEL FUNCTIONS
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

def plot(x,Y,
	linespecs=['-'],
	colors=['k'],
	linewidths=[1],
	xticks=None,
	yticks=None
):
    """
    For each y in Y, plots y against x

    Parameters
    ----------
    x, Y : Variable, (list of) Variables
            Variables to plot
	xticks, yticks :
		Allow user to override the ticks

    Notes
    -----
    Default behavior is to cycle through linespecs and colors
    """		
    if not len(linespecs) is len(colors) is len(linewidths) is len(Y):
        raise Exception('PC load letter')

	# plot 
    for y in Y:
        plt.plot(x._plt,y._plt)	

	# x-axis
    plt.xlabel(x.name)
	plt.xticks(**x.tick)
	plt.xlim(x.vran)

    if len(Y) > 1:
        plt.legend([y.name for y in Y])
    else:
        plt.ylabel(Y[0].name)
	
	# misc
	plt.grid()	

	# save
    fname_x = x.name.replace(' ','_')
    fname_y = '_'.join([y.name.replace(' ','_') for y in Y])
    plt.savefig(fname_x + '_' + fname_y + '.pdf')

# ==============================================================================
# EQUILIBRIUM

class Equilibrium:
    """Equilibrium class for journalism model"""

    def __init__(self,alpha,kappa,chi,b_bar,c_bar,v_bar,varphi):
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
        self.mu_v   = v_bar/2. 
        self.v_hat  = (1.-alpha)*b_bar
        self.K0     = (1.+2.*chi)/(2.*kappa*(1.+chi)**2.)
		self.rho 	= chi/(1.+chi)
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
		return -(_pi_R(delta,beta)*var(beta)+(1.-_pi_R(delta,beta))*var(delta))

	@np.vectorize
	def _Omega(delta,beta):
		"""drift"""
		return _pi_R(delta,beta)*(rho*beta)+(1.-_pi_R(delta,beta))*(-rho*delta)

	# ==========================================================================
	# VARIABLES

	# --------------------------------------------------------------------------
	# news
	delta = Variable(
		func = lambda x: x,
		tick = {
			'ticks' : [-mu_d,mu_d], 
			'labels' : ['$-\\mu_{d}$','$+\\mu_{d}$']
		}, 
		name = 'news',
		symb = r'$\delta$'
	)

	# --------------------------------------------------------------------------
	# bias
	b = Variable(
		func = lambda x : x,
		tick = {
			'ticks' : [0.,b_bar], 
			'labels' : ['$0$','$\\overline{b}$']
		}, 
		name = 'bias',
		symb = r'$\b$'
	)

	# --------------------------------------------------------------------------
	# <<<EQUILIBRIUM>>> bias
	b_eq = deepcopy(b)

	b_eq.func = np.vectorize(_b),
	b_eq.name = 'equilibrium bias',
	b_eq.symb = r'$\\beta(\\delta)$'
	
	# --------------------------------------------------------------------------
	# reporting probability
	pi = Variable(
		func = _pi_R,
		tick = {'ticks' :, 'labels' : }, 
		vran = (0.,1.),
		name = 'reporting probability',
		symb = r'$\pi_{R}^{*}$',
	)

	# --------------------------------------------------------------------------
	# <<<EQUILIBRIUM>>> reporting probability
	pi_eq = deepcopy(pi)
	
	pi_eq.func = lambda delta, beta: _pi(delta,_beta(delta))
	pi_eq.name = 'equilibrium reporting probability',
	pi_eq.symb = r'$\pi_{R}^{*}(v)$',
	
	# --------------------------------------------------------------------------
	# journalist's report
	s_J = Variable(
		func = _sigma_J,
		tick = {'ticks' :, 'labels' : }, 
		vran = (0.,1.),
		name = 'journalist report',
		symb = r'$s_{J}^{*}$'
	)

	# --------------------------------------------------------------------------
	# price quality
	Lambda = Variable(
		func = _Lambda,
		tick = {'ticks' :, 'labels' : }, 
		vran = None,
		name = 'price quality',
		symb = r'$\\Lambda$'
	)

	# --------------------------------------------------------------------------
	# drift
	Omega = Variable(
		func = _Omega,
		tick = {'ticks' :, 'labels' : }, 
		vran = None,
		name = 'drift',
		symb = r'$\\Omega$'
	)
