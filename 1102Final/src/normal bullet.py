from vpython import*

'''----------------------circuit rendering----------------------'''
inductor_length = 3E-2
inductor_radius = 1E-2
wire_length = 3E-2
R = 10
L = 1.2E-4
C = 220E-6
V0 = 75
N = 146

scene = canvas(width=800, height=800,
               background=vector(0.2, 0.2, 0), align='left')

inductor = helix(canvas=scene, pos=vector(-inductor_length/2, 0, 0),
                 axis=vector(inductor_length, 0, 0,), radius=inductor_radius, coils=10)
capacitor = box(pos=vector(0, 0, inductor_radius + wire_length),
                size=vector(1E-3, 8E-3, 8E-3))
wire1 = cylinder(pos=vector(-inductor_length/2, 0, inductor_radius),
                 axis=vector(0, 0, wire_length), radius=inductor_radius/20)
wire2 = cylinder(pos=vector(inductor_length/2, 0, inductor_radius),
                 axis=vector(0, 0, wire_length), radius=inductor_radius/20)
wire3 = cylinder(pos=vector(-inductor_length/2, 0, inductor_radius+wire_length),
                 axis=vector(inductor_length, 0, 0), radius=inductor_radius/20)

'''----------------------bullet rendering----------------------'''
miu0 = 4*pi*10**(-7)
bullet_mass = 3E-4
bullet_length = 2E-3
bullet_radius = 1E-3
bullet_volume = bullet_length * pi * (bullet_radius)**2
bullet_magnetic_susceptibility = 700


bullet = cylinder(pos=vector(-inductor_length/2, 0, 0),
                  axis=vector(bullet_length, 0, 0), radius=bullet_radius)
bullet.v = vector(0, 0, 0)
'''----------------------i(t)----------------------'''
lambda_1 = (-R+sqrt(R**2-4*L/C))/(2*L)
lambda_2 = (-R-sqrt(R**2-4*L/C))/(2*L)

k_1 = V0/(L*(lambda_2-lambda_1))
k_2 = V0/(L*(lambda_1-lambda_2))

'''----------------------magnetic field----------------------'''


def magnetic_field(pos_x, current):
    cos1 = (inductor_length/2 + pos_x) / sqrt((inductor_length/2 + pos_x)**2+inductor_radius**2)
    cos2 = (inductor_length/2 - pos_x) / sqrt((inductor_length/2 - pos_x)**2+inductor_radius**2)
    return miu0*N/inductor_length*current/2*(cos1+cos2)


'''---------------------- v-t graph----------------------'''
graph1 = graph(width=450, align='left')
Vt = gcurve(graph=graph1, color=color.yellow, width=4)

t = 0
dt = 1E-4
dx = 1E-6
final_v = 0

# scene.center = vector(0,0,0)
scene.camera.pos = vector(0, 1E-1, 0)
scene.camera.axis = vector(0, -1, 0)

while bullet.pos.x <= inductor_length/2*3:
    rate(2000)
    i = abs(k_1*exp(lambda_1*t)+k_2*exp(lambda_2*t))
    B = magnetic_field(bullet.pos.x, i)
    grad_B = miu0 * (N/inductor_length) * i / 2 * inductor_radius**2 * (1/((inductor_length/2+bullet.pos.x)**2+inductor_radius**2)**(3/2) +
                                                                        1/((inductor_length/2-bullet.pos.x)**2+inductor_radius**2)**(3/2))
    F = bullet_volume*bullet_magnetic_susceptibility*grad_B*B/miu0
    a = F/bullet_mass
    bullet.v.x += a*dt
    bullet.pos.x += bullet.v.x*dt
    final_v = bullet.v.x
    # print("i=", i)
    # print("v=", bullet.v.x)
    Vt.plot(t, bullet.v.x)
    t += dt
print("capacitor = ", C)
print("final speed = ", final_v)
