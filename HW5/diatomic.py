from vpython import *
size, m_o, m_c, k_bond = 31E-12, 16.0/6E23, 12.0/6E23, 18600.0 # These numbers are all made up
d = 2.5*size
dt = 1E-16
class CO_molecule:
    def __init__(self, pos, axis):
        self.O = sphere(pos = pos, radius = size, color = color.red)
        self.C = sphere(pos = pos+axis, radius = size, color = color.blue)
        self.bond = cylinder(pos = pos, axis = axis, radius = size/2.0, color = color.white)
        self.O.m = m_o
        self.C.m = m_c
        self.O.v = vector(0, 0, 0)
        self.C.v = vector(0, 0, 0)
        self.bond.k = k_bond
        
    def bond_force_on_O(self): # return bond force acted on the O atom
         return self.bond.k*(mag(self.bond.axis)-d)*norm(self.bond.axis)
        
    def time_lapse(self, dt): # by bond's force, calculate a, v and pos of C and O, and bond's pos and axis after dt
        self.C.a = - self.bond_force_on_O() / self.C.m
        self.O.a = self.bond_force_on_O() / self.O.m
        self.C.v += self.C.a * dt
        self.O.v += self.O.a * dt
        self.C.pos += self.C.v * dt
        self.O.pos += self.O.v * dt
        self.bond.axis = self.C.pos - self.O.pos
        self.bond.pos = self.O.pos

    def com(self): # return position of center of mass
        return (self.O.pos*self.O.m+self.C.pos*self.C.m)/(self.O.m+self.C.m)
    
    def com_v(self): # return velocity of center of mass
        return (self.O.v*self.O.m+self.C.v*self.C.m)/(self.O.m+self.C.m)

    def v_P(self): # return potential energy of the bond for the vibration motion
        return 1/2*self.bond.k*(mag(self.bond.axis)-d)*(mag(self.bond.axis)-d)
    
    def v_K(self): # return kinetic energy of the vibration motion
        O_v_K=1/2*self.O.m*(dot(self.O.v-self.com_v(),norm(self.bond.axis)))**2
        C_v_K=1/2*self.C.m*(dot(self.C.v-self.com_v(),norm(self.bond.axis)))**2
        return O_v_K+C_v_K

    def r_K(self): # return kinetic energy of the rotational motion
        O_K=1/2*self.O.m*(mag(self.O.v-self.com_v()))**2
        C_K=1/2*self.C.m*(mag(self.C.v-self.com_v()))**2
        return O_K+C_K-self.v_K()

    def com_K(self): #return kinetic energy of the translational motion of the center of mass
        return 1/2*(self.O.m+self.C.m)*(mag(self.com_v()))**2

def collision(a1, a2):
    v1prime = a1.v - 2 * a2.m/(a1.m+a2.m) *(a1.pos-a2.pos) * dot (a1.v-a2.v, a1.pos-a2.pos) / mag(a1.pos-a2.pos)**2
    v2prime = a2.v - 2 * a1.m/(a1.m+a2.m) *(a2.pos-a1.pos) * dot (a2.v-a1.v, a2.pos-a1.pos) / mag(a2.pos-a1.pos)**2
    return v1prime, v2prime
    
if __name__ == '__main__':
    a = CO_molecule(pos=vector(0, 0, 0), axis = vector(2.6*size, 0, 0))
    a.O.v = vector(1.0, 1.0, 0)
    a.C.v = vector(2.0, -1.0, 0)
    a.time_lapse(dt)
    print(a.bond_force_on_O(), a.com(), a.com_v(), a.v_P(), a.v_K(), a.r_K(), a.com_K())

