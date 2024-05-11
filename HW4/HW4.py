import numpy as np
from vpython import *

A, N = 0.10, 50
m, k, d = 0.1, 10.0, 0.4
dt, avgCount = 0.0003, 10

#scene = canvas(title='Spring Wave', width=800, height=300, background=vec(0.5,0.5,0), center = vec((N-1)*d/2, 0, 0))
g = graph(width=1200, height=600, align='left')
f = gcurve(graph=g, color=color.yellow, width=2)

class Obj:
    pass

balls = [Obj for i in range(N)]
springs = [Obj for i in range(N - 1)]

for i in range(1,N//2):
    Unit_K = 2.0 * pi / (N * d)
    Wavevector = i * Unit_K
    phase = Wavevector * arange(N) * d
    ball_pos, ball_v, spring_len = np.arange(N) * d + A * np.sin(phase), np.zeros(N), np.ones(N) * d

    total_t=0.0
    count=0

    while count < 2 * avgCount:
        total_t += dt
        
        spring_len[:-1] = ball_pos[1:] - ball_pos[:-1]
        ball_v[1:] += (spring_len[1:] - spring_len[:-1]) * k / m * dt
        spring_len[N - 1] = ball_pos[0] - ball_pos[N - 1] + N * d
        ball_v[0] += (spring_len[0] - spring_len[N - 1]) * k / m * dt
        
        if ball_pos[0] * (ball_pos[0] + ball_v[0] * dt) < 0 and total_t > 0.001:
            count += 1
        ball_pos += ball_v * dt
    
    T=total_t/avgCount
    print("當n =",i,"時，角頻率約為",2*pi/T)
    f.plot(Wavevector, 2.0 * pi / T)

    
