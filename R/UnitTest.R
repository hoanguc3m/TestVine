library("VineCopula")

# Student copula
family = 2
par = 0.8
u = 0.5; v = 0.5;
par2 = 5
log(BiCopPDF(u, v, family, par, par2 = par2))
BiCopDeriv(u, v, family, par, par2 = par2, deriv = "par2", log = TRUE)

BiCopHfunc2(u,v, family = family, par = par, par2 = par2)
# [1] 0.07324421
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "par")
# [1] -0.2328563
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "u2")
# [1] -0.1983125
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "par2")
# [1] 0.0004019854


# Gaussian
u = 0.2
v = 0.7
par = 0.6
family = 1

BiCopHfunc2(u,v, family = family, par = par)
# [1] 0.07418304
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "par")
# [1] -0.2822272
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "u2")
# [1] -0.3028041


# Student copula
family = 2
par = 0.6
par2 = 5
BiCopHfunc2(u,v, family = family, par = par, par2 = par2)
# [1] 0.07324421
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "par")
# [1] -0.2328563
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "u2")
# [1] -0.1983125
BiCopHfuncDeriv(u, v, family = family, par = par, par2 = par2, deriv = "par2")
# [1] 0.0004019854


# Clayton
par = 0.6
family = 3

BiCopHfunc2(u,v, family = family, par = par)
# [1] 0.1068516
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "par")
# [1] -0.1137902
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "u2")
# [1] -0.1386487

A = (u^(-par) + v^(-par) -1)
h = A^(-1/par-1) * v^(-par-1)
dh_v = (par+1) * A^(-1/par-2) *v^(2 *(-par-1)) - (par+1) * A^(-1/par-1) *v^(-par-2)



# Gumbel
par = 1.6
family = 4

BiCopHfunc2(u,v, family = family, par = par)
# [1] 0.1025011
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "par")
# [1] -0.1337794
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "u2")
# [1] -0.3150611

# Frank
par = 0.6
family = 5

BiCopHfunc2(u,v, family = family, par = par)
# [1] 0.1801582
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "par")
# [1] -0.0339533
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "u2")
# [1] -0.08862073


# Joe
par = 1.6
family = 6

BiCopHfunc2(u,v, family = family, par = par)
# [1] 0.1629389
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "par")
# [1] -0.07728109
BiCopHfuncDeriv(u, v, family = family, par = par, deriv = "u2")
# [1] -0.3067064

