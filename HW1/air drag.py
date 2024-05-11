from vpython import *
g=9.8 # g = 9.8 m/s^2
size = 0.25 # ball radius = 0.25 m
height = 15.0 # ball center initial height = 15 m
C_drag = 1.2
scene = canvas(width=600, height=600, center =vec(0,height/2,0), background=vec(0.5,0.5,0))
floor = box(length=30, height=0.01, width=10, color=color.blue)
ball = sphere(radius = size, color=color.red, make_trail = True)
ball.pos = vec(-15, size, 0)
ball.v = vec(16, 16, 0) # ball initial velocity
dt = 0.001 # time step
while ball.pos.y >= size: # until the ball hit the ground
 rate(1000) # run 1000 times per real second
 ball.v += vec(0, -g, 0) * dt - C_drag*ball.v*dt
 ball.pos += ball.v*dt
msg = text(text = 'final speed = ' + str(ball.v.mag), pos = vec(-10, 15, 0))
