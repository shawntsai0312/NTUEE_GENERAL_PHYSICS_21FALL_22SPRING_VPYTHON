from numpy.lib.function_base import add_newdoc_ufunc
from vpython import *
import numpy as np
from histogram import *

N = 200
# He atom are 10 times bigger for easiear collision but not too big for accuracy
m, size = 4E-3/6E23, 31E-12*10
# 2L is the cubic container's original length, width, and height
L = ((24.4E-3/(6E23))*N/(2*pi))**(1/3.0) + size  # V0=pi*L^2*2L
# Boltzmann Constant and initial temperature
k, T = 1.38E-23, 298.0
t, dt = 0, 3E-13
# the initial root mean square velocity
vrms = (2*k*1.5*T/m)**0.5
stage = 0                                   # stage number
atoms = []                                  # list to store atoms
times = 0
delta_P = 0
v_wall = L/(20000*dt)
C = 0
gamma = 5/3

# histogram setting
deltav = 50.                                # slotwidth for v histogram
vdist = graph(x=800, y=0, ymax=N*deltav/1000., width=500,
              height=300, xtitle='v', ytitle='dN', align='left')

# for plot of the curve for the atom speed distribution
theory_low_T1 = gcurve(color=color.cyan)
theory_low_T2 = gcurve(color=color.green)
dv = 10.
for v in arange(0., 4201.+dv, dv):            # theoretical speed distribution
    theory_low_T1.plot(pos=(v, (deltav/dv)*N*4.*pi *
                       ((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation1 = ghistogram(graph=vdist, bins=arange(
    0., 4200., deltav), color=color.red)  # for the simulation speed distribution
observation2 = ghistogram(graph=vdist, bins=arange(
    0., 4200., deltav), color=color.blue)  # for the simulation speed distribution


def keyinput(evt):
    global stage
    s = evt.key
    if s == 'n':
        stage += 1


# initialization
scene = canvas(width=500, height=500,
               background=vector(0.2, 0.2, 0), align='left')
#container = box(length=2*L, height=2*L, width=2 * L, opacity=0.2, color=color.yellow)
container = cylinder(radius=L, pos=vector(-L, 0, 0),
                     axis=vector(2*L, 0, 0), opacity=0.2, color=color.yellow)
scene.bind('keydown', keyinput)
# particle position array and particle velocity array, N particles and 3 for x, y, z
p_a, v_a = np.zeros((N, 3)), np.zeros((N, 3))

for i in range(N):
    # particle is initially random positioned in container
    random_phi = 2*pi*random()
    random_r = (L-size)*random()
    p_a[i] = [container.pos.x+container.length*random(), random_r *
              cos(random_phi), random_r*sin(random_phi)]
    # the last atom is with yellow color and leaves a trail
    if i == N-1:
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]),
                      radius=size, color=color.yellow, make_trail=True, retain=50)
    # other atoms are with random color and leaves no trail
    else:
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius=size, color=vector(
            random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    # particle initially same speed but random direction
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]
    atoms.append(atom)


# the function for handling velocity after collisions between two atoms
def vcollision(a1p, a2p, a1v, a2v):
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime


while True:
    scene.center = container.pos + vector(L,L,0)
    t += dt
    rate(1000)
    times += 1
    T = 0
    for i in range(N):
        T += m*(v_a[i, 0]**2+v_a[i, 1]**2+v_a[i, 2]**2)/(3*N*k)
    if times == 100:
        V = pi*container.radius**2*container.length
        A = 2*pi*container.radius**2+2*pi*container.radius*container.length
        p = (delta_P/100/dt)/A
        print('T=', T, 'p=', p, 'V=', V, 'p*V=', p*V,
              'N*k*T=', N*k*T, 'p*(V**gamma)=', p*(V**gamma))
        delta_P = 0
        times = 0

    p_a += v_a*dt  # calculate new positions for all atoms
    for i in range(N):
        # to display atoms at new positions
        atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2])

    if stage == 0:
        # freeze histogram for stage = 0
        observation1.plot(data=np.sqrt(np.sum(np.square(v_a), -1)))
    if stage == 1:
        container.length -= 2*v_wall*dt
        if container.length <= L:
            stage = 2
    if stage == 2:
        # freeze histogram for stage = 2
        observation2.plot(data=np.sqrt(np.sum(np.square(v_a), -1)))
        if C == 0:
            for v in arange(0, 4201+dv, dv):
                theory_low_T2.plot(pos=(v, (deltav/dv)*N*4.*pi *
                                        ((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
            C = 1
    if stage == 3:
        container.length = 2*L

# find collisions between pairs of atoms, and handle their collisions
    # array for vector from one atom to another atom for all pairs of atoms
    r_array = p_a-p_a[:, np.newaxis]
    # distance array between atoms for all pairs of atoms
    rmag = np.sqrt(np.sum(np.square(r_array), -1))
    # if smaller than 2*size meaning these two atoms might hit each other
    hit = np.less_equal(rmag, 2*size)-np.identity(N)
    # change hit to a list
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist()
    for ij in hitlist:                                                  # i,j encoded as i*Natoms+j
        # atom pair, i-th and j-th atoms, hit each other
        i, j = divmod(ij, N)
        # remove j,i pair from list to avoid handling the collision twice
        hitlist.remove(j*N+i)
        # only handling collision if two atoms are approaching each other
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0:
            v_a[i], v_a[j] = vcollision(
                p_a[i], p_a[j], v_a[i], v_a[j])  # handle collision

    # find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        delta_vy = 0
        delta_vz = 0
        if p_a[i][0] >= container.pos.x + container.length - size and v_a[i][0] > 0:
            # x right
            if stage == 1:
                delta_P += 2*m*(v_a[i][0]+v_wall)
                v_a[i][0] = -v_a[i][0]-2*v_wall
            else:
                delta_P += 2*m*abs(v_a[i][0])
                v_a[i][0] = - v_a[i][0]
                

        if p_a[i][0] <= container.pos.x + size and v_a[i][0] < 0:
            # x left
            delta_P += 2*m*abs(v_a[i][0])
            v_a[i][0] = - v_a[i][0]
            

        if (p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2 >= (container.radius-size)**2 and v_a[i][1]*(p_a[i][1]-container.pos.y) + v_a[i][2]*(p_a[i][2]-container.pos.z) >= 0:
            delta_P += 2*m*abs(v_a[i][1]*(p_a[i][1]-container.pos.y) +
                               v_a[i][2]*(p_a[i][2]-container.pos.z))/sqrt((p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2)
            delta_vy = 2*(v_a[i][1]*(p_a[i][1]-container.pos.y)+v_a[i][2] * (p_a[i][2]-container.pos.z)) / ((p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2)*(p_a[i][1]-container.pos.y)
            delta_vz = 2*(v_a[i][1]*(p_a[i][1]-container.pos.y)+v_a[i][2] * (p_a[i][2]-container.pos.z)) / ((p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2)*(p_a[i][2]-container.pos.z)
            v_a[i][1] -= delta_vy
            v_a[i][2] -= delta_vz