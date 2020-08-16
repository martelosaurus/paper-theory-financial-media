from journalism import *

E1 = Equilibrium()
E2 = Equilibrium(varphi = 1.)

#-------------------------------------------------------------------------------
# policy functions

plot(E1.delta,[E1.pi.of(E1.delta)])
plot(E1.beta,[E1.pi.of(E1.beta)])

#-------------------------------------------------------------------------------
# equilibrium

# journalist's reporting probability 
plot(E1.delta,[E1.pi.of(E1.delta),E1.delta])

# manager's bias
plot(E1.delta,[E1.beta.of(E1.delta)])

# manager's report (and straight-shooter)
plot(E1.delta,[E1.s_F.of(E1.delta),E1.delta])

# journalist's report (and straight-shooter)
plot(E1.delta,[E1.s_J.of(E1.delta),E1.delta])

#-------------------------------------------------------------------------------
# equilibrium with loss aversion

# journalist's reporting probability
plot(E2.delta,[E2.pi.of(E2.delta)])

#-------------------------------------------------------------------------------
# empirical implications

# probability of reporting and bias
plot(E1.beta.of(E1.delta),[E1.pi.of(E1.delta),E1.pi.of(E1.delta)])

# drift
plot(E1.pi.of(E1.delta),[E1.Omega.of(E1.delta)])
plot(E1.beta.of(E1.delta),[E1.Omega.of(E1.delta)])

# price quality
plot(E1.pi.of(E1.delta),[E1.Lambda.of(E1.delta)])
plot(E1.beta.of(E1.delta),[E1.Lambda.of(E1.delta)])
