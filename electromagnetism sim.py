import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Particle:
    def __init__(self,mass, position, velocity,charge):
        self.mass = mass
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity,dtype=float)
        self.charge = charge
#all in metric
#distances between points in space
dx,dy,dz = 6e-12,6e-12,6e-12
#points in space
nx,ny,nz = 50,50,50
#volume element
dv = dx*dy*dz

eps = 8.85e-12
mu = 1.256e-6
steps = 10
dt = 1e-17

#make sure locs isnt in the nx,ny,nz cubes so derivatives work nice
particle1 = Particle(9.1e-31,np.array([25*dx-5.29e-11,25*dy,25*dz]),np.array([0,0,0]),-1.602e-19)
particle2 = Particle(1.672e-27,np.array([25*dx,25*dy,25*dz]),np.array([0,0,0]),1.602e-19)
particles = [particle1,particle2]
locs = [np.zeros((steps,3)),np.zeros((steps,3))]
locs[0][0] = particles[0].position
locs[1][0] = particles[1].position

#boundry is 0
E = np.zeros((nx, ny, nz, 3))
B = np.zeros((nx, ny, nz, 3))
#EVERY TIME NI+1 IS USED FPOR THE DERIVATIVE EQUATE IT TO 0
#ensuring guassian starter conditions so the gauss equations hold true for later
for i in particles:
    sx = int(np.floor(i.position[0]/dx))
    sy = int(np.floor(i.position[1]/dy))
    sz = int(np.floor(i.position[2]/dz))
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if x!=sx or y!=sy or z!=sz:
                    bluh =((x-sx)**2+(y-sy)**2+(z-sz)**2)
                    factor = 1/(np.sqrt(((x-sx)**2+(y-sy)**2+(z-sz)**2)))**3
                    E[x,y,z] = E[x,y,z]+np.array([x-sx,y-sy,z-sz])*factor
#order of steps
#calculate fields
#apply forces and update locs and vels
correcB = []
correcE = []
sizecorB = 0
sizecorE = 0
for l in range(1,steps):
    print(l)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                #now we will look at the derivatives, which here are forward (deriv at x,y,z is diff between x+1 and x), this sucks
                derEX,derEY,derEZ,derBX,derBY,derBZ = np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3)
                if x+1<nx:
                    derEX = (E[x + 1, y, z] - E[x, y, z]) / dx
                    derBX = (B[x + 1, y, z] - B[x, y, z]) / dx
                if y + 1 < ny:
                    derEY = (E[x, y + 1, z] - E[x, y, z]) / dy
                    derBY = (B[x, y + 1, z] - B[x, y, z]) / dy
                if z+1<nz:
                    derEZ = (E[x, y, z + 1] - E[x, y, z]) / dz
                    derBZ = (B[x, y, z + 1] - B[x, y, z]) / dz

                if derEX.any() or derEY.any() or derEZ.any():
                    correcB.append([[x,y,z],[derEZ[1]-derEY[2],derEX[2]-derEZ[0],derEY[0]-derEX[1]]])
                    sizecorB = sizecorB + 1
                    j = [0,0,0]
                    for i in particles:
                        tempx = np.floor(i.position[0] / dx)
                        tempy = np.floor(i.position[1] / dy)
                        tempz = np.floor(i.position[2] / dz)
                        if ([tempx,tempy,tempz] == [x,y,z]):
                            dens = i.charge/dv
                            j = [i.velocity[0]*dens,i.velocity[1]*dens,i.velocity[2]*dens]
                    if j!=0 or derBX!=0 or derBY!=0 or derBZ!=0:
                        tempV = [(derBY[2]-derBZ[1])/(mu*eps)-j[0]/eps,(derBZ[0]-derBX[2])/(mu*eps)-j[1]/eps,(derBX[1]-derBY[0])/(mu*eps)-j[2]/eps]
                        correcE.append([[x,y,z],tempV])
                        sizecorE = sizecorE + 1
    #i have found all the corrections needed to be done to the fields. now we need to use them to update our fields
    #NEED LENGTH OF CORRECTION MATRICES, GO THROGUH THEM CHANGE ELEMENT ELEMENT, DONT TRY TO BE SMART AND EFFICIENT, YOU AINT TURIN
    for i in range(sizecorE):
        x,y,z = correcE[i][0][0],correcE[i][0][1],correcE[i][0][2]
        arr = np.array([correcE[i][1][0],correcE[i][1][1],correcE[i][1][2]])
        E[x,y,z] = E[x,y,z] + arr*dt
    for i in range(sizecorB):
        x,y,z = correcB[i][0][0],correcB[i][0][1],correcB[i][0][2]
        arr = np.array([correcB[i][1][0],correcB[i][1][1],correcB[i][1][2]])
        B[x,y,z] = B[x,y,z] + arr*dt
    #DONE WITH UPDATING FIELDS, NOW NEED TO CALCULATE FORCE ACTING ON EACH PARTICLE, I.E RETRIEVE THE FIELD AT THE LOC OF THE PARTILCE. MAYBE BORIS ALGO?
    #HERE WE GO
    for i in particles:
        x = int(np.floor(i.position[0] / dx))
        y = int(np.floor(i.position[1] / dy))
        z = int(np.floor(i.position[2] / dz))
        VxB = [i.velocity[1]*B[x,y,z,2]-i.velocity[2]*B[x,y,z,1],i.velocity[2]*B[x,y,z,0]-i.velocity[0]*B[x,y,z,2],i.velocity[0]*B[x,y,z,1]-i.velocity[1]*B[x,y,z,0]]
        acc=(E[x,y,z]+VxB)*i.charge/i.mass
        i.position = i.position + i.velocity * dt
        i.velocity = i.velocity + acc * dt
        if (i.position[0]>= (nx-1) * dx):
            i.position[0] = 10*dx
            #punishment
            i.velocity[0] = 0
        if (i.position[1] >= (ny-1) * dy):
            i.position[1] = 10*dy
            i.velocity[1] = 0
        if (i.position[2] >= (nz-1) * dz):
            i.position[2] = 10*dz
            i.velocity[2] = 0
        #this position updating part using lorentz seems to work, now only locs updating needed, fk me never mind need to make sure particle is inside the boundries
        #THIS IS USELESS START LOOKING AT ENERGY, PROBABLY NOT CONSERVED
        #oh, it worked, my bad code compiled, huh]
    locs[0][l] = particles[0].position
    locs[1][l] = particles[1].position