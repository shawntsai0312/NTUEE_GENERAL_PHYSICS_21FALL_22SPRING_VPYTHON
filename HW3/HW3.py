from vpython import *
G=6.673E-11
mass = {'earth': 5.97E24, 'moon': 7.36E22, 'sun':1.99E30}
radius = {'earth': 6.371E6*10, 'moon': 1.317E6*10, 'sun':6.95E8*10} #10 times larger for better view
earth_orbit = {'r': 1.495E11, 'v': 2.9783E4}
moon_orbit = {'r': 3.84E8, 'v': 1.022E3}
theta = 5.145*pi/180.0
def G_force(m1,m2,pos_vec1,pos_vec2):
    return -G * m1 * m2 / mag2(pos_vec2-pos_vec1) * norm(pos_vec2-pos_vec1)

scene = canvas(width=800, height=600, background=vector(0.5,0.5,0))

earth = sphere(radius = radius['earth'], m = mass['earth'], texture = {'file':textures.earth})
moon = sphere(radius = radius['moon'], m = mass['moon'])
sun = sphere(pos = vector(0,0,0), radius = radius['sun'], m = mass['sun'] , color = color.orange, emissive = True)

moon.pos = vector(-moon_orbit['r']*earth.m/(earth.m+moon.m)*cos(theta),-moon_orbit['r']*earth.m/(earth.m+moon.m)*sin(theta),0) + vector(earth_orbit['r'], 0, 0) 
earth.pos = vector(moon_orbit['r']*moon.m/(earth.m+moon.m)*cos(theta),moon_orbit['r']*moon.m/(earth.m+moon.m)*sin(theta),0) + vector(earth_orbit['r'], 0, 0) 

moon.v = vector(0, 0, -moon_orbit['v']*earth.m/(earth.m+moon.m)) + vector(0, 0, -earth_orbit['v'])
earth.v = vector(0,0,moon_orbit['v']*moon.m/(earth.m+moon.m)) + vector(0, 0, -earth_orbit['v'])

sev = arrow(color=color.green,shaftwidth=radius['earth']/20)
mev = arrow(color=color.yellow,shaftwidth=radius['earth']/20)
sev.pos = earth.pos
mev.pos = earth.pos
sev.axis = norm(cross(earth.pos-sun.pos,earth.v))*radius['earth']*2
mev.axis = norm(cross(moon.pos-earth.pos,moon.v-earth.v))*radius['earth']*(-2)

a1 = arrow(color=color.red,shaftwidth=radius['earth']/20)
a1.pos = earth.pos
a1.axis = mev.axis

t=0
dt=60*60*6
T=0
times=1
sum_T=0 
while True:
    rate(1000)   
    t+=dt
    moon.a = (G_force(earth.m , moon.m , earth.pos , moon.pos) + G_force(sun.m , moon.m , sun.pos , moon.pos)) / moon.m
    earth.a = (G_force(earth.m , moon.m , moon.pos , earth.pos) + G_force(sun.m , earth.m , sun.pos , earth.pos)) / earth.m
    moon.v += moon.a * dt
    earth.v += earth.a * dt
    moon.pos += moon.v * dt
    earth.pos += earth.v * dt
    sev.pos = earth.pos
    mev.pos = earth.pos
    sev.axis = norm(cross(earth.pos-sun.pos,earth.v))*radius['earth']*2
    mev.axis = norm(cross(moon.pos-earth.pos,moon.v-earth.v))*radius['earth']*(-2)
    scene.center = earth.pos
    k = dot(a1.axis,mev.axis)/mag(a1.axis)/mag(mev.axis)
    phi = acos(k)
    a1.pos = earth.pos
    
    if round(phi*180)==0 and T>=500000000 :
        sum_T+=T
        print("第",times,"次的週期約為",round(T/60/60/24),"天;目前平均約為",round(sum_T/times/60/60/24),"天")
        T=0
        times+=1
    else :
        T+=dt
