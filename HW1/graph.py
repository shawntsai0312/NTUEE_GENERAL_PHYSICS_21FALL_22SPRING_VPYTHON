from vpython import *
scene1 = canvas(width = 200, align = 'left', background = vec(0, 0.5, 0.5))
scene2 = canvas(width = 300, height = 300, align = 'left', background = vec(0.5, 0.5, 0))
box(canvas = scene1)
sphere(canvas = scene2)
oscillation = graph(width = 450, align = 'right')
funct1 = gcurve(graph = oscillation, color=color.blue, width=4)
funct2 = gvbars(graph = oscillation, delta=0.4, color=color.red)
funct3 = gdots(graph = oscillation, color=color.orange, size=3)
t = 0
while t < 80:
 rate(25)
 t = t+1
 funct1.plot( pos=(t, 5.0+5.0*cos(-0.2*t)*exp(0.015*t)) )
 funct2.plot( pos=(t, 2.0+5.0*cos(-0.1*t)*exp(0.015*t)) )
 funct3.plot( pos=(t, 5.0*cos(-0.03*t)*exp(0.015*t)) )
