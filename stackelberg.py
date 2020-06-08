class Stackelberg:
    """
	class for Stackelberg model

	Model Description
	-----------------
	There are two dates, t=1,2, and two firms, j=1,2. On date 1, firm 1 produces 	
        q1 at marginal cost c>0. On date 2, firm 2 observes q1, and produces q2 at 
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
        # x objects
        self.xobj = []

        q1 = ParVar(
            r'$0$',
            r'$\mu_{v}-\widehat{v}$',
            r'$\mu_{v}$',
            r'$\mu_{v}+3\widehat{v}$',
            r'$\overline{v}$'
        )

        a = ParVar(
            r'$0$',
            r'$\\sup_{\\Delta_{R}}^{-1}(\\overline{c})$',
            r'$\\sup_{\\Delta_{R}}^{-1}(0)$'
        )

		b = ParVar(
			'hello',
		)

        # ---------------------------------------------------------------------
        # auxiliary functions
		@np.vectorize
		def _P(q1,q2,a):
			"""market price"""
			return a-q1-q2

		@np.vectorize
		def _pi1(q1,q2,a,c):
			"""firm 1's profit"""
			return (P(q1,q2,a)-c)*q1

		@np.vectorize
		def _pi2(q1,q2,a,c):
			"""firm 2's profit"""
			return (P(q1,q2,a)-c)*q2

        # ---------------------------------------------------------------------
        # y objects
        self.yobj = []

		@np.vectorize
		def _q2(q1,a,c):
			"""firm 2's response to firm 1's"""
			return ((a-c)-q1)/2.

		@np.vectorize
		def _q1_star(a,c):
			"""firm 1's equilibrium quantity"""
			return (a-c)/2.

		@np.vectorize
		def _q2_star(a,c):
			"""firm 2's equilibrium quantity"""
			return (a-c)/4.

        self.yobj(Object(
            name = 'equilibrium reporting probability',
            func = _pi_R,
            ylab = r'reporting probability $\pi_{R}^{*}$',
            1.,
            r'$1$'
            ))

