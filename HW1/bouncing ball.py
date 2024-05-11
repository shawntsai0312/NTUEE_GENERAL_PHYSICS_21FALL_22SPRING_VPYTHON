from vpython import *
g=9.8 # g = 9.8 m/s^2
size = 0.25 # ball radius = 0.25 m
scene = canvas(center = vec(0,5,0), width=600, background=vec(0.5,0.5,0))
floor = box(length=30, height=0.01, width=4, color=color.blue)
ball = sphere(radius = size, color=color.red, make_trail=True, trail_radius = size/3)

ball.pos = vec( -15.0, 10.0, 0.0) # ball initial position
ball.v = vec(2.0, 0.0 , 0.0) # ball initial velocity
dt = 0.001
while ball.pos.x < 15.0: # simulate until x=15.0m
 rate(1000)
 ball.pos += ball.v*dt
 ball.v.y += - g*dt
 if ball.pos.y <= size and ball.v.y < 0: # new: check if ball hits the ground
  ball.v.y = - ball.v.y # if so, reverse y component of velocity
