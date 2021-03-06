from journalism import *

E1 = Equilibrium()
E2 = Equilibrium(varphi = 1.)

#-------------------------------------------------------------------------------
# policy functions

plot(E1.delta,[E1.pi_of_delta])
plot(E1.beta,[E1.pi_of_beta])

#-------------------------------------------------------------------------------
# manager's bias
plot(E1.delta,[E1.beta])

# manager's report (and straight-shooter)
plot(E1.delta,[E1.delta,E1.s_F])

# journalist's report (and straight-shooter)
plot(E1.delta,[E1.delta,E1.s_J])

#-------------------------------------------------------------------------------
# equilibrium with loss aversion

# journalist's reporting probability
plot(E2.delta,[E1.pi_of_delta,E2.pi_of_delta])

#-------------------------------------------------------------------------------
# empirical implications

# probability of reporting and bias
plot(E1.beta,[E1.pi_of_delta],E1.delta)

# drift
plot(E1.pi_of_delta,[E1.Omega],E1.delta)
plot(E1.beta,[E1.Omega])

# price quality
plot(E1.pi_of_delta,[E1.Lambda],E1.delta)
plot(E1.beta,[E1.Lambda])
