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

class Sim:
    def __init__(self,grid_shape,grid_spacing,dt):
        self.nx, self.ny, self.nz = grid_shape
        self.dx,self.dy,self.dz = grid_spacing
        self.dt = dt
        self.eps0 = 8.854187817e-12
        self.mu0 = 4 * np.pi * 1e-7
        self.c = 1 / np.sqrt(self.eps0 * self.mu0)
        stable = self.c*dt/min(self.dx,self.dy,self.dz)
        if(stable>=1):
            print("Ain't no good, this simulation won't be stable. light is too fast")
        #a field is a function R^3->R^3, for a computer R^3->R is a cube with each grid point having a number.
        #we devide the actual field into three such arrays. less confusing. each gives one component of the actual field at a point
        self.Ex = np.zeros(grid_shape)
        self.Ey = np.zeros(grid_shape)
        self.Ez = np.zeros(grid_shape)
        self.Bx = np.zeros(grid_shape)
        self.By = np.zeros(grid_shape)
        self.Bz = np.zeros(grid_shape)
        self.rho = np.zeros(grid_shape)
        self.Jx = np.zeros(grid_shape)
        self.Jy = np.zeros(grid_shape)
        self.Jz = np.zeros(grid_shape)
        self.particles= []
        self.time =0.0
        self.history = []
    def add_particle(self, particle):
        self.particles.append(particle)

    def charge_magic(self):
        self.rho.fill(0)
        self.Jx.fill(0)
        self.Jy.fill(0)
        self.Jz.fill(0)
        volume = self.dx*self.dy*self.dz
        for p in self.particles:
            roundx = int(np.round(p.position[0]/self.dx))
            roundy = int(np.round(p.position[1] / self.dy))
            roundz = int(np.round(p.position[2] / self.dz))
            # really important boundry stuff, what happens when particle skidaddle outside of the box
            # this one laps around, pacman style, aborted pacman, doing reflection
            roundx = np.clip(roundx, 0, self.nx - 1)
            roundy = np.clip(roundy, 0, self.ny - 1)
            roundz = np.clip(roundz, 0, self.nz - 1)
            if 0<=roundx<=self.nx and 0<=roundy<=self.ny and 0<=roundz<=self.nz:
                self.rho[roundx, roundy, roundz] += p.charge / volume
                self.Jx[roundx, roundy, roundz] += p.charge* p.velocity[0] / volume
                self.Jx[roundx, roundy, roundz] += p.charge * p.velocity[1] / volume
                self.Jx[roundx, roundy, roundz] += p.charge * p.velocity[2] / volume

    def update_B(self):
        E = [self.Ex,self.Ey,self.Ez]
        curlE_dirX,curlE_dirY, curlE_dirZ = curl(E,self.dx,self.dy,self.dz)
        self.Bx -= curlE_dirX*self.dt
        self.By -= curlE_dirY*self.dt
        self.Bz -= curlE_dirZ*self.dt
    def update_E(self):
        B = [self.Bx, self.By, self.Bz]
        curlB_dirX, curlB_dirY, curlB_dirZ = curl(B, self.dx, self.dy, self.dz)
        self.Ex += self.dt * (self.c * self.c * curlB_dirX - self.Jx / self.eps0)
        self.Ey += self.dt * (self.c * self.c * curlB_dirY - self.Jy / self.eps0)
        self.Ez += self.dt * (self.c * self.c * curlB_dirZ - self.Jz / self.eps0)

    def field_at_pos(self,pos):
        roundx = int(np.round(pos[0] / self.dx))
        roundy = int(np.round(pos[1] / self.dy))
        roundz = int(np.round(pos[2] / self.dz))

        #really important boundry stuff, what happens when particle skidaddle outside of the box
        #this one laps around, pacman style, pacman is a nogo, doing reflection instead
        roundx = np.clip(roundx, 0, self.nx - 1)
        roundy = np.clip(roundy, 0, self.ny - 1)
        roundz = np.clip(roundz, 0, self.nz - 1)

        E = np.array([self.Ex[roundx,roundy,roundz],self.Ey[roundx,roundy,roundz],self.Ez[roundx,roundy,roundz]])
        B = np.array([self.Bx[roundx,roundy,roundz],self.By[roundx,roundy,roundz],self.Bz[roundx,roundy,roundz]])
        return E,B
    def update_particles(self):
        for p in self.particles:
            E,B = self.field_at_pos(p.position)
            p.position += p.velocity * self.dt
            p.velocity += (p.charge/p.mass)*(E+np.cross(p.velocity,B))*self.dt
            # really important boundry stuff, what happens when particle skidaddle outside of the box
            # this one laps around, pacman style
            #Edit, pacman worked horribly, everything just got out of control, this will do reflect the particles
            if p.position[0] < 0:
                p.position[0] = -p.position[0]
                p.velocity[0] = -p.velocity[0]
            elif p.position[0] > self.nx * self.dx:
                p.position[0] = 2 * (self.nx * self.dx) - p.position[0]
                p.velocity[0] = -p.velocity[0]
            if p.position[1] < 0:
                p.position[1] = -p.position[1]
                p.velocity[1] = -p.velocity[1]
            elif p.position[1] > self.ny * self.dy:
                p.position[1] = 2 * (self.ny * self.dy) - p.position[1]
                p.velocity[1] = -p.velocity[1]
            if p.position[2] < 0:
                p.position[2] = -p.position[2]
                p.velocity[2] = -p.velocity[2]
            elif p.position[2] > self.nz * self.dz:
                p.position[2] = 2 * (self.nz * self.dz) - p.position[2]
                p.velocity[2] = -p.velocity[2]

    def step(self):
        self.charge_magic()
        #i hate computers and i hate discrete things, they are cute in analysis, they are anger inducing here.
        self.update_B()
        self.update_E()
        self.update_particles()
        self.time +=self.dt
        self.history.append({'time': self.time,'particles': [{'pos': p.position.copy(), 'vel': p.velocity.copy()} for p in self.particles]})

    def run_the_thing(self,number_of_steps):
        print(f"Running awesome minkowski magic with {number_of_steps} steps")
        print(f"Grid is: {self.nx}x{self.ny}x{self.nz}")
        print(f"Total time: {self.dt * number_of_steps:.2e} s")
        for i in range(number_of_steps):
            self.step()
            if (i + 1) % 10 == 0:
                print(f"  Step {i + 1}/{number_of_steps}")
        print("DONE")


if __name__ == "__main__":
    # Grid setup (3D!)
    nx, ny, nz = 50, 50, 50
    dx, dy, dz = 2e-9, 2e-9, 2e-9  # 2 nanometer spacing (coarser for speed)

    # Time step (must satisfy Courant condition)
    c = 3e8  # speed of light
    dt = 1 * min(dx, dy, dz) / c  # Factor of 0.3 for safety

    sim = Sim((nx, ny, nz), (dx, dy, dz), dt)

    # Physical constants
    e = 1.602176634e-19  # Elementary charge
    m_e = 9.10938356e-31  # Electron mass
    m_p = 1.6726219e-27  # Proton mass

    # Add particles in 3D configuration
    # Two electrons orbiting
    sim.add_particle(Particle(
        position=[25e-9, 25e-9, 40e-9],
        velocity=[5e5, 5e5, 0],
        charge=-e,
        mass=m_e
    ))

    sim.add_particle(Particle(
        position=[75e-9, 75e-9, 60e-9],
        velocity=[-5e5, -5e5, 0],
        charge=-e,
        mass=m_e
    ))

    # Add a proton in the middle
    sim.add_particle(Particle(
        position=[50e-9, 50e-9, 50e-9],
        velocity=[0, 0, 2e5],
        charge=e,
        mass=m_e
    ))

    # Another electron coming from top
    sim.add_particle(Particle(
        position=[50e-9, 50e-9, 90e-9],
        velocity=[0, 0, -3e5],
        charge=-1*e,
        mass=10*m_e
    ))

    # Run simulation
    num_steps = 160
    sim.run_the_thing(num_steps)

    # 3D Animation
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlim(0, nx * dx * 1e9)
    ax.set_ylim(0, ny * dy * 1e9)
    ax.set_zlim(0, nz * dz * 1e9)
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_zlabel('z (nm)')
    ax.set_title('3D EM Particle Simulation')

    # Color code particles by charge
    colors = ['red' if p.charge > 0 else 'blue' for p in sim.particles]
    sizes = [100 if p.charge > 0 else 50 for p in sim.particles]

    # Initialize with first frame data
    first_positions = np.array([p['pos'] * 1e9 for p in sim.history[0]['particles']])
    scatter = ax.scatter(first_positions[:, 0], first_positions[:, 1], first_positions[:, 2],
                         c=colors, s=sizes, alpha=0.8)
    time_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes)

    # Trails for each particle
    trails = [ax.plot([], [], [], c=colors[i], alpha=0.3, linewidth=1)[0]
              for i in range(len(sim.particles))]


    def animate(frame):
        state = sim.history[frame]
        positions = np.array([p['pos'] * 1e9 for p in state['particles']])

        # Update particle positions
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])

        # Update trails (show last 30 frames)
        trail_length = 30
        for i, trail in enumerate(trails):
            trail_data = np.array([sim.history[f]['particles'][i]['pos'] * 1e9
                                   for f in range(max(0, frame - trail_length), frame + 1)])
            if len(trail_data) > 0:
                trail.set_data(trail_data[:, 0], trail_data[:, 1])
                trail.set_3d_properties(trail_data[:, 2])

        time_text.set_text(f'Step: {frame}/{len(sim.history) - 1}\nTime: {state["time"] * 1e12:.2f} ps')

        # Rotate view slowly
        ax.view_init(elev=20, azim=frame * 0.5)

        return scatter, *trails, time_text


    anim = FuncAnimation(fig, animate, frames=len(sim.history),
                         interval=50, blit=False, repeat=True)
    plt.show()

    print(f"\nCourant number: {c * dt / min(dx, dy, dz):.3f}")
    print(f"Grid: {nx}×{ny}×{nz} = {nx * ny * nz:,} points")
    print(f"Time per step: ~{(dx / c) * 1e15:.2f} femtoseconds")