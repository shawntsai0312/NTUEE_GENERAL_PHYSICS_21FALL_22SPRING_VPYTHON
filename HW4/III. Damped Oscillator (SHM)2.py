from vpython import *
size, m = 0.02, 0.2 # ball size = 0.02 m, ball mass = 0.2kg
L, k = 0.2, 20 # spring original length = 0.2m, force constant = 20 N/m
fa = 0.1
b = 0.05 * m * sqrt(k/m)
omega_d = sqrt(k/m)
'''
scene = canvas(width=600, height=400, fov = 0.03, align = 'left', center=vec(0.3, 0, 0), background=vec(0.5,0.5,0))
wall_left = box(length=0.005, height=0.3, width=0.3, color=color.blue) # left wall
ball = sphere(radius = size, color=color.red) # ball
spring = helix(radius=0.015, thickness =0.01)
oscillation1 = graph(width = 400, align = 'left', xtitle='t',ytitle='x',background=vec(0.5,0.5,0))
x=gcurve(color=color.red,graph = oscillation1)
'''
oscillation2 = graph(width = 400, align = 'left', xtitle = 't', ytitle = 'average_power', background=vec(0.5,0.5,0))
p = gdots(color=color.cyan, graph = oscillation2)
class obj: pass ####
ball, spring = obj(), obj() ####
spring.pos = vector(0, 0, 0)
ball.pos = vector(L, 0 , 0) # ball initial position
ball.v = vector(0, 0, 0) # ball initial velocity
ball.m = m
T = 2*pi / omega_d
t, dt, n = 0, 0.001, 1
power = 0.0
while True:
    #rate(1000)
    spring.axis = ball.pos - spring.pos
    spring_force = - k * (mag(spring.axis) - L) * norm(spring.axis)
    force = vector(fa*sin(omega_d*t), 0, 0)
    ball.a = (force + spring_force - b * ball.v ) / m
    ball.v += ball.a*dt
    ball.pos += ball.v*dt
    t += dt
    #x.plot(pos=(t,ball.pos.x - L))
    power += dot(force, ball.v)*dt
    if t / T >= n:
        p.plot(pos = (t, power / T))
        n += 1.0
        power = 0
