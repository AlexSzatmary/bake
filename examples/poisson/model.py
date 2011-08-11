#!/usr/bin/python

# x_l and x_r are the left and right bounds on the domain
x_l = 0.
x_r = 1.
# u_l and u_r are the left and right boundary conditions: u(x_l) = u_l
u_l = 0.
u_r = 0.
# number of nodes, total; n-2 are actually used on the domain.
n = 11
# h is the grid spacing
h = (x_r - x_l)/(n - 1)

# Forcing function for the right hand side; this is a regular function of x
def f(x):
    return @forcing_function@
