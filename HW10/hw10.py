from vpython import*

fd = 120  # 120Hz
T = 1/fd

R = 30
L = 0.2
C = 20*10**(-6)

Q = 0
i = 0
E_12 = 0

t = 0
dt = 1.0/(fd * 5000)  # 5000 simulation points per cycle

scene1 = graph(align='left', xtitle='t',
               ytitle='i (A) blue, v (100V) red,', background=vector(0.2, 0.6, 0.2))
scene2 = graph(align='left', xtitle='t', ytitle='Energy (J)',
               background=vector(0.2, 0.6, 0.2))

i_t = gcurve(color=color.blue, graph=scene1)
v_t = gcurve(color=color.red, graph=scene1)
E_t = gcurve(color=color.red, graph=scene2)

flag_i = False
flag_v = False
flag_E = False

i_prev = -1
v_prev = -1

i_max_t = 0
v_max_t = 0

while t <= 12*T:
    rate(10000)
    v = 36*sin(2*pi*fd*t)
    vR = i*R
    vC = Q/C
    vL = v-vC-vR
    i += vL/L*dt
    Q += i*dt
    E = 1/2*C*vC**2+1/2*L*i**2
    i_t.plot(t, i)
    v_t.plot(t, v/100)
    E_t.plot(t, E)
    t += dt
    E_12 = E
    if t >= 9*T:
        if flag_i == False:
            if i_prev < i:
                i_prev = i
            else:
                flag_i = True
                i_max_t = t
                print(f"I magnitude = {i_prev}")
        if flag_v == False:
            if v_prev < v:
                v_prev = v
            else:
                flag_v = True
                v_max_t = t
                print(f"V magnitude = {v_prev}")
        
if flag_i == True and flag_v == True:
    print(f"phase difference = {(v_max_t-i_max_t)/T*360}")

while t > 12*T and t <= 20*T:
    rate(10000)
    v = 0
    vR = i*R
    vC = Q/C
    vL = v-vC-vR
    i += vL/L*dt
    Q += i*dt
    E = 1/2*C*vC**2+1/2*L*i**2
    i_t.plot(t, i)
    v_t.plot(t, v/100)
    E_t.plot(t, E)
    t += dt
    if E < E_12*0.1 and flag_E == False:
        flag_E=True
        print(t-12*T)
