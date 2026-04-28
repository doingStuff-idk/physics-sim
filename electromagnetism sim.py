import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def gradient(scalar_field,dx,dy,dz):
    grad_z, grad_y, grad_x = np.gradient(scalar_field, dz, dy, dx)
    return grad_z,grad_y,grad_z

def div(vector_field,dx,dy,dz):
    Fx, Fy, Fz = vector_field
    dFx_dx = np.gradient(Fx,dx,axis=0)
    dFy_dy = np.gradient(Fy, dy, axis=1)
    dFz_dz = np.gradient(Fz, dz, axis=2)
    return dFz_dz+dFy_dy+dFx_dx

def curl(vector_field,dx,dy,dz):
    Fx,Fy,Fz = vector_field
    dFx_dy = np.gradient(Fx,dy,axis=1)
    dFy_dz = np.gradient(Fy, dz, axis=2)
    dFz_dx = np.gradient(Fz, dx, axis=0)
    dFx_dz = np.gradient(Fx,dz,axis=2)
    dFy_dx = np.gradient(Fy, dx, axis=0)
    dFz_dy = np.gradient(Fz, dy, axis=1)

    curl_x = dFz_dy-dFy_dz
    curl_y = dFx_dz-dFz_dx
    curl_z = dFy_dx-dFx_dy
    return curl_x,curl_y,curl_z

class Particle:
    def __init__(self,mass, position, velocity,charge):
        self.mass = mass
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity,dtype=float)
        self.charge = charge
#all in metric
#distances between points in space
dx,dy,dz = 2e-12,2e-12,2e-12
#points in space
nx,ny,nz = 50,50,50
#volume element
dv = dx*dy*dz

eps = 8.85e-12
mu = 0

#make sure locs isnt in the nx,ny,nz cubes so derivatives work nice
steps = 1000
particle1 = Particle(9.1e-31,np.array([25*dx-5.29e-11,25*dy,25*dz]),np.array([0,0,0]),-1.602e-19)
particle2 = Particle(1.672e-27,np.array([25*dx,25*dy,25*dz]),np.array([0,0,0]),1.602e-19)
particles = [particle1,particle2]
locs = [np.zeros((steps,3)),np.zeros((steps,3))]
locs[0][0] = particles[0].position
locs[1][0] = particles[1].position

#boundry is 0
Electric = np.zeros((nx, ny, nz, 3))
Magnetic = np.zeros((nx, ny, nz, 3))
#EVERY TIME NI+1 IS USED FPOR THE DERIVATIVE EQUATE IT TO 0
#ensuring guassian starter conditions so the gauss equations hold true for later
for i in particles:
    x = np.floor(particles[i].position[0]/dx)
    y = np.floor(particles[i].position[1]/dy)
    z = np.floor(particles[i].position[2]/dz)
    if x+1<nx:
        Electric[x+1, y, z, 0] = (particles[i].charge/dv)*(dx/eps)
    if x+1>=nx:
        #leaving it here for later
        Electric[x, y, z,0] = 0
#order of steps
#calculate fields
#apply forces and update locs and vels
correcB = []
correcE = []
for l in range(1,steps):
    for x in nx:
        for y in ny:
            for z in nz:
                #now we will look at the derivatives, which here are forward (deriv at x,y,z is diff between x+1 and x), this sucks
                derEX,derEY,derEZ,derBX,derBY,derBZ = 0,0,0,0,0,0
                if x + 1 >= nx:
                    derEX=0
                    derBX=0
                if x+1<nx:
                    derEX = (E[x + 1, y, z] - E[x, y, z]) / dx
                    derBX = (B[x + 1, y, z] - B[x, y, z]) / dx
                if y + 1 >= ny:
                    derEY=0
                    derBY=0
                if y + 1 < ny:
                    derEY = (E[x, y + 1, z] - E[x, y, z]) / dy
                    derBY = (B[x, y + 1, z] - B[x, y, z]) / dy
                if z + 1 >= nz:
                    derEZ=0
                    derBZ=0
                if z+1<nz:
                    derEZ = (E[x, y, z + 1] - E[x, y, z]) / dz
                    derBZ = (B[x, y, z + 1] - B[x, y, z]) / dz

                if derEX!=0 or derEY!=0 or derEZ!=0:
                    correcB.append([[x,y,z],[derEZ[1]-derEY[2],derEX[2]-derEZ[0],derEY[0]-derEX[1]]])
                    J = 0
                    for i in particles:
                        tempx = np.floor(particles[i].position[0] / dx)
                        tempy = np.floor(particles[i].position[1] / dy)
                        tempz = np.floor(particles[i].position[2] / dz)
                        if ([tempx,tempy,tempz] == [x,y,z]):
                            dens = particles[i].charge/dv
                            J = [particles[i].velocity[0]*dens,particles[i].velocity[1]*dens,particles[i].velocity[2]*dens]
                    if J!=0 or derBX!=0 or derBY!=0 or derBZ!=0:
                        tempV = [(derBY[2]-derBZ[1])/(mu*eps)-j[0]/eps,(derBZ[0]-derBX[2])/(mu*eps)-j[1]/eps,(derBX[1]-derBY[0])/(mu*eps)-j[2]/eps]
                        correcE.append([[x,y,z],tempV])
    #i have found all the corrections needed to be done to the fields. now we need to use them to update our fields
    #NEED LENGTH OF CORRECTION MATRICES, GO THROGUH THEM CHANGE ELEMENT ELEMENT, DONT TRY TO BE SMART AND EFFICIENT, YOU AINT TURIN

