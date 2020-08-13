from journalism import *

e = Equilibrium(
	alpha = 1./2.,
	kappa = 2./3.,
	chi   = 1./2.,
	b_bar = 4./3.,
	c_bar = 10.,
	v_bar = 4.
)

plot(e.v,[e.b])
