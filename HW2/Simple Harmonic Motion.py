from vpython import *
g = 9.8
size, m = 0.05, 0.2
L, k = 0.5, 15
scene = canvas(width=500, height=500, center=vec(0, -0.2, 0), background=vec(0.5,0.5,0))
ceiling = box(length=0.8, height=0.005, width=0.8, color=color.blue)
ball = sphere(radius = size, color=color.red)
spring = helix(radius=0.02, thickness =0.01) # default pos = vec(0, 0, 0)
ball.v = vec(0, 0, 0)
ball.pos = vec(0, -L, 0)
dt = 0.001
while True:
 rate(1000)
 spring.axis = ball.pos - spring.pos # new: extended from spring endpoint to ball

 spring_force = - k * (mag(spring.axis) - L) * spring.axis.norm() # to get spring force vector
 ball.a = vector(0, - g, 0) + spring_force / m # ball acceleration = - g in y + spring force /m
 ball.v += ball.a*dt
 ball.pos += ball.v*dt
