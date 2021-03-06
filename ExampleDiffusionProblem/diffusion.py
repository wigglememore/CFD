#-------------------- intro --------------------
#     Program to solve Laplace equation over the domain 0 <= x <= 1,  0 <= y <= 1
#     with boundary conditions phi=x*y

#-------------------- library imports --------------------
import numpy as np

#-------------------- parameters --------------------
nx = 25.0 #size of grid in x
ny = 25.0 #size of grid in y

omega = 1.7 #<1 is under relaxed, >1 is over relaxed
max_iterations = 500

#-------------------- initialisation --------------------
dx = 1.0/nx
dy = 1.0/ny

s_u = np.zeros(nx, ny)
s_p = np.zeros(nx, ny)
phi = np.zeros(nx + 1, ny + 1)

a_e = np.ones(nx, ny)
a_w = np.ones(nx, ny)
a_n = np.ones(nx, ny)
a_s = np.ones(nx, ny)

#-------------------- edge boundary conditions --------------------
# west boundary
a_w[:, 0] = 0.0
s_u[:, 0] = 0.0
s_p[:, 0] = -2.0

#north boundary
a_n[0, :] = 0.0
s_p[0, :] = -2.0
for i in range(0, nx):
    s_u[ny+1, i] = 2.0 * ((dx/2) + dx*)i - 1))

#east boundary
a_e[:, nx-1] = 0.0
s_p[:, nx-1] = -2.0
for j in range(0, ny):
    s_u[j, nx-1] = 2.0 * ((dy/2) + dy*(j - 1))

#south boundary
a_s[ny-1, :] = 0.0
s_u[ny-1, :] = 0.0
s_p[ny-1, :] = -2.0

#-------------------- corner boundary conditions --------------------
#northwest
a_n[0, 0] = 0.0
a_w[0, 0] = 0.0
s_u[0, 0] = 2.0 * (dx/2)
s_p[0, 0] = -4.0

#northeast
a_n[0, nx-1] = 0.0
a_e[0, nx-1] = 0.0
s_u[0, nx-1] = 2.0 * ((dx/2) + dx*(nx-1)) + 2.0 * ((dy/2) + dy*(ny - 1))
sp[0, nx-1] = -4.0

#southeast
a_s[ny-1, nx-1] = 0.0
a_e[ny-1, nx-1] = 0.0
s_u[ny-1, nx-1] = 2.0 * (dy/2)
s_p[ny-1, nx-1] = -4.0

#southwest
a_s[ny-1, 0] = 0.0
a_w[ny-1, 0] = 0.0
s_u[ny-1, 0] = 0.0
s_p[ny-1, 0] = -4.0

a_p = a_e + a_w + a_n + a_s - s_p

#-------------------- iterative solver --------------------
for iterations in range(1, max_iterations):
    for i in range(0, nx-1):
        for j in range(0, ny-1):
            phi[i, j] = phi[i, j] + (omega/a_p[i,j]) + (a_e[i,j])*phi[i+1,j] + a_w[i,j]*phi[i-1,j] + a_n[i,j]*phi[i,j+1] + a_s[i,j]*phi[i,j-1] + s_u[i,j] - a_p[i,j]*phi[i,j])

#-------------------- write csv for paraview --------------------
