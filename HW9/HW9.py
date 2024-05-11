from vpython import *
import numpy as np

miu0 = 4*pi*10**(-7)
R_s = 0.06
R_l = 0.12
d = 0.1

n_ring = 1000
n_radius = 1000


# now calculating the flux for small ring
flux_s = 0
dr = R_s/n_radius
for i in range(n_radius):

    r = i * dr
    pos_r = vec(r,0,d)
    B = vec(0,0,0)

    for j in range(n_ring):

        theta = 2*pi*j/n_ring
        pos_s = R_l*vec(cos(theta),sin(theta),0)
        ds = 2*pi*R_l/n_ring*vec(-sin(theta),cos(theta),0)
        dB = miu0/(4*pi)*cross(ds,pos_r-pos_s)/(mag(pos_r-pos_s))**3
        B += dB
    
    flux_s += 2*pi*r*dr*B.z

print(flux_s)

# now calculating the flux for large ring
flux_l = 0
dr = R_l/n_radius
for i in range(n_radius):

    r = i * dr
    pos_r = vec(r,0,d)
    B = vec(0,0,0)

    for j in range(n_ring):

        theta = 2*pi*j/n_ring
        pos_s = R_s*vec(cos(theta),sin(theta),0)
        ds = 2*pi*R_s/n_ring*vec(-sin(theta),cos(theta),0)
        dB = miu0/(4*pi)*cross(ds,pos_r-pos_s)/(mag(pos_r-pos_s))**3
        B += dB
    
    flux_l += 2*pi*r*dr*B.z

print(flux_l)