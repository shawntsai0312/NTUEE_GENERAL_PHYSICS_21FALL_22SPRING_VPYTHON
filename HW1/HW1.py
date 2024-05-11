from vpython import*
g = 9.8
size = 0.25
height = 0.0
C_drag = 0.9
theta = pi/4
scene = canvas(width = 600,height = 300,algin='left',background=vec(0.5,0.5,0))
floor = box(length = 30,height = 0.01,width = 10,color = color.blue)
ball = sphere(radius = size,color = color.red,make_trail = True)
ball.pos = vec(-15,size,0)
ball.v = vec(20*cos(theta),20*sin(theta),0)
dt = 0.001
count = 0
v1 = arrow(color=color.green,shaftwidth=0.05)
v1.pos = ball.pos
v1.axis = ball.v/2
max_height=0
distance=0
vt_graph = graph(width=450,align='right')
funct1 = gcurve(graph=vt_graph,color=color.yellow,width=4)
t=0

while count<3:
    rate(1/dt)
    ball.v += vec(0,-g,0)*dt-C_drag*ball.v*dt
    ball.pos += ball.v*dt
    v1.pos = ball.pos
    v1.axis = ball.v/2
    distance+=mag(ball.v)*dt
    funct1.plot(t,mag(ball.v))
    t+=dt
    if max_height<=ball.pos.y :
        max_height = ball.pos.y
    if ball.pos.y<=size and ball.v.y<0:
        ball.v.y=-ball.v.y
        count+=1        

displacement = ball.pos-vec(-15,size,0)


msg = text(text='lagest height = '+str(max_height),pos=vec(-10,-4,0))
msg = text(text='displacement = '+str(displacement),pos=vec(-10,-5.5,0))
msg = text(text='distance = '+str(distance),pos=vec(-10,-7,0))
