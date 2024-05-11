from vpython import *
scene = canvas(width=1000, height=1000, background=vec(0.5,0.5,0)) # open a window
a1 = arrow(color = color.green, shaftwidth = 0.05)
b1 = arrow(color = color.blue, shaftwidth = 0.05)
a1.pos = vec(1, 1, 0)
a1.axis = vec(1, -1, 0)
b1.pos = a1.pos + a1.axis
b1.axis = vec(2, 1, 0)
c1 = arrow(color = color.yellow, shaftwidth=0.05)
c1.pos = a1.pos
c1.axis = a1.axis + b1.axis
