from vpython import *
dt, k, g ,L = 0.0001, 150000, vector(0,-9.8,0), 2
scene = canvas(width=800, height=600, center=vec(0, -1, 0), background=vec(0.5,0.5,0))

def af_col_v(m1, m2, v1, v2, x1, x2): # function after collision velocity
 v1_prime = v1 + 2*(m2/(m1+m2))*(x1-x2) * dot (v2-v1, x1-x2) / dot (x1-x2, x1-x2)
 v2_prime = v2 + 2*(m1/(m1+m2))*(x2-x1) * dot (v1-v2, x2-x1) / dot (x2-x1, x2-x1)
 return (v1_prime, v2_prime)

graph1 = graph(width=450,align='right')
Kt = gcurve(graph=graph1,color=color.yellow,width=4)
Ut = gcurve(graph=graph1,color=color.green,width=4)
graph2 = graph(width=450,align='right')
aKt = gcurve(graph=graph2,color=color.yellow,width=4)
aUt = gcurve(graph=graph2,color=color.green,width=4)

balls = []
pivots = []
rods = []

N=2

for i in range(0,N):
 balls.append(sphere(pos = vec(0.4*(i-2)-sqrt(2*2-1.95*1.95),-1.95,0), radius = 0.2, color = color.white))
for i in range(N,5):
 balls.append(sphere(pos = vec(0.4*(i-2),-2,0), radius = 0.2, color = color.white))
for i in range(0,5):
 pivots.append(sphere(pos = vec(0.4*(i-2),0,0), radius = 0.05, color = color.white))
 rods.append(cylinder(pos = vec(0.4*(i-2),0,0), axis = balls[i].pos-pivots[i].pos, radius = 0.02))
 balls[i].m = 1
 balls[i].v = vec(0,0,0)

avg_K = 0
avg_U = 0
t = 0
spring_force = [0,0,0,0,0]
while(True):
    rate(5000)
    t+=dt
    total_K = 0
    total_U = 0
    for i in range(0,5):
        rods[i].axis = balls[i].pos-pivots[i].pos
        spring_force[i] = - k * (mag(rods[i].axis) - L) * rods[i].axis.norm()
        balls[i].a = g + spring_force[i]/balls[i].m
        balls[i].v += balls[i].a*dt
        balls[i].pos += balls[i].v*dt
        total_K += balls[i].m*mag(balls[i].v)*mag(balls[i].v)/2
        total_U += balls[i].m*mag(g)*(balls[i].pos.y+2)
    avg_K += total_K
    avg_U += total_U
    Kt.plot(t,total_K)
    Ut.plot(t,total_U)
    aKt.plot(t,avg_K/t)
    aUt.plot(t,avg_U/t)
    for i in range(0,4):
        if (mag(balls[i].pos - balls[i+1].pos) <= 0.4 and dot(balls[i].pos-balls[i+1].pos, balls[i].v-balls[i+1].v) <= 0) :
            (balls[i].v, balls[i+1].v) = af_col_v (balls[i].m, balls[i+1].m, balls[i].v, balls[i+1].v, balls[i].pos, balls[i+1].pos)

    
