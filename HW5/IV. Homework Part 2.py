from vpython import *
from diatomic import *

N = 20 # 20 molecules
L = ((24.4E-3/(6E23))*N)**(1/3.0)/50 # 2L is the length of the cubic container box, the number is made up
m = 14E-3/6E23 # average mass of O and C
k, T = 1.38E-23, 298.0 # some constants to set up the initial speed
initial_v = (3*k*T/m)**0.5 # some constant

scene = canvas(width = 400, height =400, align = 'left', background = vec(1, 1, 1))
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.4, color = color.yellow )
energies = graph(width = 600, align = 'left', ymin=0)
c_avg_com_K = gcurve(color = color.green)
c_avg_v_P = gcurve(color = color.red)
c_avg_v_K = gcurve(color = color.purple)
c_avg_r_K = gcurve(color = color.blue)

COs=[]

for i in range(N): # initialize the 20 CO molecules
    O_pos = vec(random()-0.5, random()-0.5, random()-0.5)*L # random() yields a random number between 0 and 1
    CO = CO_molecule(pos=O_pos, axis = vector(1.0*d, 0, 0)) # generate one CO molecule
    CO.C.v = vector(initial_v*random(), initial_v*random(), initial_v*random()) # set up the initial velocity of C randomly
    CO.O.v = vector(initial_v*random(), initial_v*random(), initial_v*random()) # set up the initial velocity of O randomly
    COs.append(CO) # store this molecule into list COs
    
times = 0 # number of loops that has been run
dt = 5E-16
t = 0
total_com_K=0
total_v_K=0
total_v_P=0
total_r_K=0
while True:
    rate(3000)
    t+=dt
    times+=1
    for CO in COs:
        CO.time_lapse(dt)
        
    for i in range(N-1): # the first N-1 molecules
        for j in range(i+1,N): # from i+1 to the last molecules, to avoid double checking
            ##TO DO:check and handle the collisions between the atoms of different molecules
            ##C,C
            if mag(COs[i].C.pos-COs[j].C.pos)<=2*size and dot(COs[i].C.pos-COs[j].C.pos,COs[i].C.v-COs[j].C.v)<=0:
                (COs[i].C.v,COs[j].C.v)=collision(COs[i].C,COs[j].C)
            ##C,O
            if mag(COs[i].C.pos-COs[j].O.pos)<=2*size and dot(COs[i].C.pos-COs[j].O.pos,COs[i].C.v-COs[j].O.v)<=0:
                (COs[i].C.v,COs[j].O.v)=collision(COs[i].C,COs[j].O)
            ##O,C
            if mag(COs[i].O.pos-COs[j].C.pos)<=2*size and dot(COs[i].O.pos-COs[j].C.pos,COs[i].O.v-COs[j].C.v)<=0:
                (COs[i].C.v,COs[j].C.v)=collision(COs[i].C,COs[j].C)
            ##O,O
            if mag(COs[i].O.pos-COs[j].O.pos)<=2*size and dot(COs[i].O.pos-COs[j].O.pos,COs[i].O.v-COs[j].O.v)<=0:
                (COs[i].O.v,COs[j].O.v)=collision(COs[i].O,COs[j].O)
            
    
    for CO in COs:
        ##TO DO:check and handle the collision of the atoms of all molecules on all 6 walls
        ##C.x
        if CO.C.pos.x+size>=L or CO.C.pos.x-size<=-L:
            CO.C.v.x*=-1
        ##C.y
        if CO.C.pos.y+size>=L or CO.C.pos.y-size<=-L:
            CO.C.v.y*=-1
        ##C.z
        if CO.C.pos.z+size>=L or CO.C.pos.z-size<=-L:
            CO.C.v.z*=-1
        ##O.x
        if CO.O.pos.x+size>=L or CO.O.pos.x-size<=-L:
            CO.O.v.x*=-1
        ##O.y
        if CO.O.pos.y+size>=L or CO.O.pos.y-size<=-L:
            CO.O.v.y*=-1
        ##O.z
        if CO.O.pos.z+size>=L or CO.O.pos.z-size<=-L:
            CO.O.v.z*=-1
            
    ## sum com_K, v_K, v_P, and r_K for all molecules, respectively, to get total_com_K, total_v_K, total_v_P, total_r_K at the current moment
    for CO in COs:
        total_com_K+=CO.com_K()
        total_v_K+=CO.v_K()
        total_v_P+=CO.v_P()
        total_r_K+=CO.r_K()
    ## calculate avg_com_K to be the time average of total_com_K since the beginning of the simulation, and do the same for others..
    avg_com_K=total_com_K/times
    avg_v_K=total_v_K/times
    avg_v_P=total_v_P/times
    avg_r_K=total_r_K/times
    ## plot avg_com_K, avg_v_K, avg_v_P, and avg_r_K
    c_avg_com_K.plot(t,avg_com_K)
    c_avg_v_K.plot(t,avg_v_K)
    c_avg_v_P.plot(t,avg_v_P)
    c_avg_r_K.plot(t,avg_r_K)
