from typing import Container
from numpy.lib.function_base import add_newdoc_ufunc
from vpython import *
import numpy as np

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
v_wall = L/(2000*dt)
C = 0
gamma = 5/3
conpress_ratio = 0.43

energy_t = graph(width = 450, align = 'right')
bullet_Ek_t = gcurve(graph = energy_t, color=color.blue, width=4)
container_Ek_t = gcurve(graph = energy_t, color=color.red, width=4)
helium_U_t = gcurve(graph = energy_t, color=color.green, width=4)
total_E_t = gcurve(graph = energy_t, color=color.black, width=4)

def keyinput(evt):
    global stage
    s = evt.key
    if s == 'n':
        stage += 1


# initialization
scene = canvas(width=500, height=500,
               background=vector(0.2, 0.2, 0), align='left',color=color.black)
container = cylinder(radius=L, pos=vector(-L, 0, 0),
                     axis=vector(2*L, 0, 0), opacity=0.2, color=color.red)
barrel = cylinder(radius=container.radius+0.2*L,pos=container.pos,axis=vector(11*L,0,0),opacity=0.2)
bullet = cylinder(radius=container.radius-0.1*L, pos=container.pos +
              vector(container.length*conpress_ratio, 0, 0), axis=vector(2*L, 0, 0),opacity=0)
scene.bind('keydown', keyinput)

p_a, v_a = np.zeros((N, 3)), np.zeros((N, 3))

bullet.m = (2/4)*(N/6E23)
bullet.v = vector(0,0,0)

container.m = bullet.m*600

for i in range(N):
    # particle is initially random positioned in container
    random_phi = 2*pi*random()
    random_r = (L-size)*random()
    p_a[i] = [container.pos.x+container.length*random(), random_r *
              cos(random_phi), random_r*sin(random_phi)]
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

while bullet.pos.x-container.pos.x <= 10*L:
    scene.center = bullet.pos+bullet.axis/2
    t += dt
    rate(1000)
    times += 1
    T = 0
    for i in range(N):
        T += m*(v_a[i, 0]**2+v_a[i, 1]**2+v_a[i, 2]**2)/(3*N*k)

    bullet_Ek_t.plot(t,1/2*bullet.m*mag(bullet.v)**2)
    container_Ek_t.plot(t,1/2*container.m*mag(bullet.v*(bullet.m/container.m))**2)
    helium_U_t.plot(t,3/2*N*k*T)
    total_E_t.plot(t,1/2*bullet.m*mag(bullet.v)**2+1/2*container.m*mag(bullet.v*(bullet.m/container.m))**2+3/2*N*k*T)

    if times == 1000:
        V = pi*container.radius**2*container.length
        A = 2*pi*container.radius**2+2*pi*container.radius*container.length
        p = (delta_P/1000/dt)/A
        print('bullet speed = ',bullet.v.x)
        delta_P = 0
        times = 0

    p_a += v_a*dt  # calculate new positions for all atoms
    for i in range(N):
        atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2])
    
    bullet.pos += bullet.v * dt

    if stage == 1:
        container.length -= 2*v_wall*dt
        if container.length <= 2*L*conpress_ratio:
            stage = 2
    if stage == 2:
        # freeze histogram for stage = 2
        if C == 0:
            C = 1
    if stage == 3:
        bullet.opacity = 0.2    

# find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:, np.newaxis]
    rmag = np.sqrt(np.sum(np.square(r_array), -1))
    hit = np.less_equal(rmag, 2*size)-np.identity(N)
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist()
    for ij in hitlist:
        i, j = divmod(ij, N)
        hitlist.remove(j*N+i)
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0:
            v_a[i], v_a[j] = vcollision(
                p_a[i], p_a[j], v_a[i], v_a[j])

    # find collisions between the atoms and the walls, and handle their elastic collisions
    for i in range(N):
        delta_vy = 0
        delta_vz = 0
        stage_3_delta_v = 0
        stage_3_bullet_delta_v = 0

        if stage == 3 and p_a[i][0] >= bullet.pos.x - size and v_a[i][0] > 0:
            v_cm = ( m * v_a[i][0] + bullet.m * bullet.v.x ) / (m + bullet.m)
            stage_3_delta_v = 2 * (v_cm-v_a[i][0])
            stage_3_bullet_delta_v = 2 * (v_cm-bullet.v.x)
            v_a[i][0] += stage_3_delta_v
            bullet.v.x += stage_3_bullet_delta_v
            delta_P += 2*m*stage_3_delta_v

        if stage != 3 and p_a[i][0] >= container.pos.x + container.length - size and v_a[i][0] > 0:
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
            delta_vy = 2*(v_a[i][1]*(p_a[i][1]-container.pos.y)+v_a[i][2] * (p_a[i][2]-container.pos.z)) / (
                (p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2)*(p_a[i][1]-container.pos.y)
            delta_vz = 2*(v_a[i][1]*(p_a[i][1]-container.pos.y)+v_a[i][2] * (p_a[i][2]-container.pos.z)) / (
                (p_a[i][1]-container.pos.y)**2 + (p_a[i][2]-container.pos.z)**2)*(p_a[i][2]-container.pos.z)
            v_a[i][1] -= delta_vy
            v_a[i][2] -= delta_vz
