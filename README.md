# paper-theory-financial-media
Numerical solutions for A Theory of Financial Media

```python
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

        firm_value  = ParVar(
            r'$0$',
            r'$\mu_{v}-\widehat{v}$',
            r'$\mu_{v}$',
            r'$\mu_{v}+3\widehat{v}$',
            r'$\overline{v}$'
        )

        firm_bias   = ParVar(
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
```
