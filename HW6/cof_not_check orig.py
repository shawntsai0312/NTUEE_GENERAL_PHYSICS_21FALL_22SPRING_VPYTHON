from vpython import *
import numpy as np
from histogram import *
#import math



height=150
g=9.8
bsize = 0.05
theta = pi/2
C_drag = 0.45
air_density = 1.22
mass = 0.2



N = 200
stop=True
m, size = 4E-3/6E23*4, 31E-12*10
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size 
k, T = 1.38E-23, 1000.0
t, dt = 0, 3E-13
dP=0
B_a_test=10
B_mass=1
B_V=vector(0,0,0)
vrms = (2*k*1.5*T/m)**0.5
stage = 0
atoms = []
sceneb= canvas(width=600, height=600, align = 'right', center =vec(0,height/2,0), background=vec(0.2,0.2,0.2))

v_W=0
oscillation = graph(width = 450, align = 'right')
contain_v=vector(0,0,0)
contain_pos=vector(-L,0,0)
contain_m=B_mass*10
cof_Vol=2
lenght_barrel=12*L
scene = canvas(width=800, height=500, background=vector(0.2,0.2,0), align = 'left')
container = cylinder(radius=L,pos=vector(-L,0,0),axis=vector(2*L,0,0), opacity=0.2, color = color.yellow )

container.v=0
# bullet=box(pos=vector(L+L/cof_Vol,0,0),length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.red )
bullet=sphere(pos=vector(2*L/cof_Vol,0,0),radius=L, opacity=0.2, color = color.red )
barrel=cylinder(opacity=0.2,pos=vector(-L,0,0),axis=vector(lenght_barrel+2*L+L/cof_Vol,0,0),radius=1.1*L,color=color.white)
# barrel_big=cylinder(opacity=0.2,pos=vector(-L/cof_Vol,0,0),axis=vector(2*L/cof_Vol,0,0),radius=2*L)
# histogram setting
deltav = 50.
vdist = graph(x=800, y=0, ymax = N*deltav/1000.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left') 
theory_low_T = gcurve(color=color.cyan) # for plot of the curve for the atom speed distribution
dv = 10.
for v in arange(0.,4201.+dv,dv): # theoretical speed distribution
    theory_low_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation = ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.red)# for the simulation speed distribution
observation1= ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.blue)#adiabatic
observation2= ghistogram(graph = vdist, bins=arange(0.,4200.,deltav), color=color.yellow)#free expansion
#initialization

p_a, v_a = np.zeros((N,3)), np.zeros((N,3)) # particle position array and particle velocity array, N particles and 3 for x, y, z 
for i in range(N):
    p_a[i] = [1/3*L*random() , 1/3*L*random() , 1/3*L*random() ] # particle is initially random positioned in container 
    if i== N-1: # the last atom is with yellow color and leaves a trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50) 
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random())) 
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction 
    atoms.append(atom)
def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms 
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime
def keyinput(evt):
    global stage
    s=evt.key
    if s =='n':
        stage=stage+1
    if(stage==2):
        print('ready to fire')
    elif(stage==3):
        print('shot')
scene.bind('keydown',keyinput)
def calculateT():
    global v_a,N,m,T
    Ek=0
    
    Ek+=m*np.sum(np.square(v_a))/2
    T=2*Ek/(3*N*k)
    for v in arange(0.,4201.+dv,dv): # theoretical speed distribution
        theory_low_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))


times=0
while stop:

    container.pos+=contain_v*dt
    barrel.pos=container.pos
    t += dt 
    bullet.pos+=B_V*dt
    times+=1
    rate(10000)
    W = mag(container.axis)/2
    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N): atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
    if stage ==0 :
         observation.plot(data = np.sqrt(np.sum(np.square(v_a),-1))) ## freeze histogram for stage != 1
    ### find collisions between pairs of atoms, and handle their collisions
    r_array = p_a-p_a[:,np.newaxis]
    rmag = np.sqrt(np.sum(np.square(r_array),-1)) 
    hit = np.less_equal(rmag,2*size)-np.identity(N) 
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist() 

    for ij in hitlist:
        i, j = divmod(ij,N)
        hitlist.remove(j*N+i)
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0 :
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j]) # handle collision                                                                      # only handling collision if two atoms are approaching each other                                                                             #find collisions between the atoms and the walls, and handle their elastic collisions
    if stage<=2:
        for i in range(N):
            if p_a[i][0] >= container.pos.x + container.length - size and p_a[i][0]*v_a[i][0] > 0 :
                
                v_a[i][0] = - (abs(v_a[i][0])+2*v_W)*v_a[i][0]/abs(v_a[i][0])
                dP+=2*m*(abs(v_a[i][0]+2*abs(v_W)))

            if p_a[i][0] <= -L+ size and p_a[i][0]*v_a[i][0] > 0 :
                
                v_a[i][0] = - (abs(v_a[i][0])+2*v_W)*v_a[i][0]/abs(v_a[i][0])
                dP+=2*m*(abs(v_a[i][0]+2*abs(v_W)))

            if sqrt(p_a[i][1]**2+p_a[i][2]**2)>=L-size and (p_a[i][1]*v_a[i][1]+p_a[i][2]*v_a[i][2])>0:
                v_a[i][1] += -2*(v_a[i][1]*p_a[i][1]+v_a[i][2]*p_a[i][2])/(p_a[i][1]**2+p_a[i][2]**2)*p_a[i][1]
                v_a[i][2] += -2*(v_a[i][1]*p_a[i][1]+v_a[i][2]*p_a[i][2])/(p_a[i][1]**2+p_a[i][2]**2)*p_a[i][2]
                dP+=2*m*abs(v_a[i][1]*p_a[i][1]+v_a[i][2]*p_a[i][2])/sqrt(p_a[i][1]**2+p_a[i][2]**2)

            """if abs(p_a[i][1]) >= L - size and p_a[i][1]*v_a[i][1] > 0 :
                v_a[i][1] = - v_a[i][1]
                dP+=2*m*(abs(v_a[i][1]))
            if abs(p_a[i][2]) >= L - size and p_a[i][2]*v_a[i][2] > 0 :
                v_a[i][2] = - v_a[i][2]
                dP+=2*m*(abs(v_a[i][2]))"""
    elif stage>2:
        for i in range(N):
            if p_a[i][0] <=  + size+container.pos.x and p_a[i][0]*v_a[i][0] > 0 :
                contain_v.x=2*m/(m+contain_m)*v_a[i][0]+(contain_m-m)/(contain_m+m)*contain_v.x
                VAtom=(m-contain_m)/(contain_m+m)*v_a[i][0]+2*contain_m/(m+contain_m)*contain_v.x
                dP+=2*m*(abs(v_a[i][0]+2*abs(v_W)))
                v_a[i][0]=VAtom
            if p_a[i][0]>=(bullet.pos.x-L-size):
                VB=2*m/(m+B_mass)*v_a[i][0]+(B_mass-m)/(B_mass+m)*B_V.x
                VAtom=(m-B_mass)/(B_mass+m)*v_a[i][0]+2*B_mass/(m+B_mass)*B_V.x
                B_V=vec(VB,0,0)
                v_a[i][0]=VAtom
            if sqrt(p_a[i][1]**2+p_a[i][2]**2)>=L-size and (p_a[i][1]*v_a[i][1]+p_a[i][2]*v_a[i][2])>0:
                v_a[i][1]=-v_a[i][1]
                v_a[i][2]=-v_a[i][2]

    if stage==1:
        v_W = L / (20000.0*dt)
        container.length-=v_W*dt
        if container.length<=2*L/cof_Vol:
            stage=2
    if stage==2:###place the bullet
        v_W=0
        
        observation1.plot(data = np.sqrt(np.sum(np.square(v_a),-1)))
    if stage==3:#shot the bullet
        
        #container.length=2*L
        
        observation2.plot(data = np.sqrt(np.sum(np.square(v_a),-1)))
    if times%1000==0:
        calculateT()
        print('T=',T)
        p=dP/(1000*dt)/(2*(container.length*container.height+container.length*container.width+container.height*container.width))
        V=container.length*container.width*container.height
        print('P=',p,'V=',V,'p*V=',p*V,'NkT=',N*k*T,'p*V^5/3=',p*(V**(5/3)))
        print('speed of bullet=',B_V.x)
        dP=0
        print('Ek total for gas=',N*k*T*3/2,'Ek bullet=',1/2*B_mass*mag(B_V)**2,"sum=",contain_m/2*mag(contain_v)**2+N*k*T*3/2+1/2*B_mass*mag(B_V)**2)
        times=0
    if bullet.pos.x>=13*L:
        stop=False







floor = box(canvas=sceneb, length=1500, height=0.01, width=10, color=color.green)
ball = sphere(canvas=sceneb, radius = bsize, color=color.blue, make_trail = True, retain = 50)
ball.velocity = mag(B_V)
ball.pos = vec(-15, bsize, 0)
ball.v = vec(ball.velocity*cos(theta), ball.velocity*sin(theta), 0)
a1 = arrow(canvas=sceneb, color = color.red, shaftwidth = 0.1)
a1.pos = ball.pos
a1.axis = ball.v

vfunc = gcurve(graph=oscillation, color=color.blue, width=4)

dt = 0.001
i=0
k = 0
t=0
dis=0
ymax=0

while i <= 0:
    rate(1000)
    t+=dt
    vfunc.plot(pos=(t, ball.v.mag))
    dis+=ball.v.mag*dt
    ball.v += vec(0, -g, 0) * dt - C_drag*(ball.v.mag**2)*air_density*pi*(bsize**2)*(1/2)*dt*ball.v.norm()/mass
    ball.pos+=ball.v*dt
    a1.pos = ball.pos
    a1.axis = ball.v
    if ball.pos.y>ymax:
        ymax = ball.pos.y
    if ball.pos.y<bsize:
        i+=1
    sceneb.center = ball.pos
    