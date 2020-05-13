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
# functions
@np.vectorize
def _part(mu_v,v_hat):
    """Partition: returns Boolean list"""
    return np.array([-mu_v,-v_hat,0.,3*v_hat,mu_v])

@np.vectorize
def _omega(mu_v,omega_hat,v_hat):
    """equilibrium obfuscation"""
    return [
            lambda v : omega_hat, 
            lambda v : (1./(3.*v_hat))*(v-mu_v)*omega_hat,
            lambda v : (1./v_hat)*(v-mu_v)*omega_hat,
            lambda v : omega_hat
        ]
    return np.piecewise()
    
@np.vectorize
def _pi_R(mu_v,K0,v_hat):
    """equilibrium reporting probability"""
    return [
            lambda v : K0*((v-mu_v)**2.-v_hat**2.),
            lambda v : (8./9.)*K0*(v-mu_v)**2.,
            lambda v : 0.,
            lambda v : K0*((v-mu_v)**2.-v_hat**2.)
        ]
    return np.piecewise()

@np.vectorize
def _sigma_J(alpha):
    """Equilibrium report (conditional on reporting)"""
    return v+(1.-alpha)*_omega(v)

# ------------------------------------------------------------------------------
# equilibrium 
class Equilibrium:

    def __init__(self,alpha,kappa,chi,omega_bar,c_bar,v_bar):
        """
        Parameters
        ----------
        alpha : float
        kappa : float
        chi : float
        omega_bar : float
            Highest permissible level of obfuscation
        c_bar : float
            Highest opportunity cost for the journalist
        v_bar : float
            
        """

        # organic
        self.alpha      = alpha
        self.kappa      = kappa
        self.chi        = chi
        self.omega_bar  = omega_bar
        self.c_bar      = c_bar
        self.v_bar      = v_bar

        # synthetic
        self.mu_v       = v_bar/2. 
        self.v_hat      = (1.-alpha)*omega_bar
        self.K0         = (1.+2.*chi)/(2.*kappa*c_bar*(1.+chi)**2.)
        
        ylabs = [
                    r'obfuscation $\omega^{*}$',
                    r'reporting probability $\pi_{R}^{*}$',
                    r'report $s_{J}^{*}$'
                ]
        yticks = [
                    omega_bar,
                    1.,
                    v_bar+(1.-alpha)*omega_bar]
        yticklabels = [
                    r'$\overline{\omega}$',
                    r'$1$',
                    r'$\overline{v}+(1-\alpha)\overline{\omega}$',
                ]

        # parameter checks
        omega_bar_max = self.mu_v/(3.*(1.-alpha))
        if omega_bar > omega_bar_max:
            print('Highest permissible level of obfuscation too high')
            raise AttributeError

        c_bar_min = 16.*(1.+2.*chi)*self.mu_v**2./(9.*kappa*(1.+chi)**2.)
        if c_bar < c_bar_min:
            print('Highest opportunity cost for the journalist too low')
            raise AttributeError

    @np.vectorize
    def omega(var_dict):
        """(wrapper for) equilibrium obfuscation"""
        return _omega(var_dict)
        
    @np.vectorize
    def pi_R(var_dict):
        """(wrapper for) equilibrium reporting probability"""
        return _pi_R(var_dict)
        
    @np.vectorize
    def sigma_J(var_dict):
        """(wrapper for) equilibrium report (conditional on reporting)"""
        return _sigma_J(var_dict)
        
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
        v_plot = np.linspace(0.,v_bar,n_plot)

        # figures
        fig, axs = plt.subplots(nrows=3)
        fig.set_figheight(figh)
        fig.set_figwidth(figw)
        fig.set_tight_layout(True)

        # common axis attributes
        F = [_omega, _pi_R, _sigma_J] 
        #for j in [0,1,2]:
        for ax in axs:

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
