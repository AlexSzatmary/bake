#!/usr/bin/python

import math

# x_l and x_r are the left and right bounds on the domain
x_l = @x_l@
x_r = @x_r@
# u_l and u_r are the left and right boundary conditions: u(x_l) = u_l
u_l = @u_l@
u_r = @u_r@
# number of nodes, total; n-2 are actually used on the domain.
n = @n@
# h is the grid spacing
h = (x_r - x_l)/(n - 1)

# Forcing function for the right hand side; this is a regular function of x
def f(x):
    return @forcing_function@
