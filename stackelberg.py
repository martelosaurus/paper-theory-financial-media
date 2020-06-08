class Stackelberg:
    """
	class for Stackelberg model

	Notes
	-----
	There are two dates, t=1,2, and two firms, j=1,2. On date 1, firm 1 produces 	q1 at marginal cost c>0. On date 2, firm 2 observes q1, and produces q2 at 
	marginal cost c>0. Firms sell their quantities at the market price 
	P(q1,q2) = a-q1-q2, where a>0 is the choke price. 
	"""

    def __init__(self,a,c):
        """
        Parameters
        ----------
        a : float
            Choke price
        c : float
            Marginal cost
        """

        # ---------------------------------------------------------------------
        # parameters

        # organic parameters
        self.a = a
        self.c = c

        # synthetic parameters
        self.mu_v   = v_bar/2. 
        
        # ---------------------------------------------------------------------
        # parameter checks

        if a < 0.:
            raise Exception('a<0')

        if c < 0.:
            raise Exception('c<0')

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

		@np.vectorize
		def _q1_star(a,c):
			"""firm 1's equilibrium quantity"""

		

        # ---------------------------------------------------------------------
        # equilibrium reporting probability
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
