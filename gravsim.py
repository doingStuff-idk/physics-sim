import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D


class Body:
    def __init__(self,mass, location, velocity):
        self.mass = mass
        self.location = np.array(location, dtype=float)
        self.velocity = np.array(velocity,dtype=float)

#acceleration exerted on self by other
    def acc(self, other, G=6.67430e-11):
        if self is other:
            return np.array([0.0,0.0,0.0])
        direction = other.location-self.location
        norm = np.linalg.norm(direction)
        if (norm!=0):
            force_mag = other.mass*G/(norm**3)
            force = force_mag*direction
            return force
        return np.array([0.0,0.0,0.0])


def load_bodies_from_file(filename):
    bodies = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            first_bracket_start = line.find('[')
            first_bracket_end = line.find(']')
            second_bracket_start = line.find('[', first_bracket_end)
            second_bracket_end = line.find(']', second_bracket_start)

            mass_str = line[:first_bracket_start].strip()
            mass = float(mass_str)

            loc_str = line[first_bracket_start + 1:first_bracket_end]
            location = [float(x) for x in loc_str.split(',')]

            vel_str = line[second_bracket_start + 1:second_bracket_end]
            velocity = [float(x) for x in vel_str.split(',')]
            bodies.append(Body(mass, location, velocity))

    return bodies
#actual physics sim
dt=400
steps = 1000000
objects = load_bodies_from_file('bodies.txt')
n = len(objects)
acceleration = np.zeros((n,3))
Locations = np.zeros((steps,n,3))
Velocities = np.zeros((steps,n,3))

for i in range(n):
    Locations[0,i] = objects[i].location
    Velocities[0, i] = objects[i].velocity

current_step = 1
while (current_step < steps):
    for i in range(n):
        acceleration[i] = np.zeros(3)
        for j in range(n):
            acceleration[i] += objects[i].acc(objects[j])

#now we will calculate the change of position, THEN the change in velocity
        objects[i].location += dt*objects[i].velocity
        objects[i].velocity += dt*acceleration[i]
        Locations[current_step,i] = objects[i].location
        Velocities[current_step,i] = objects[i].velocity
    current_step+=1

#claude animation junk
def animate_simulation(Locations, dt, skip_frames=100):
    """
    Animate the gravity simulation in 3D.

    Parameters:
    - Locations: array of shape (steps, n_bodies, 3)
    - dt: time step used in simulation
    - skip_frames: only show every Nth frame (makes animation faster)
    """
    n_steps, n_bodies, _ = Locations.shape

    # Create figure and 3D axis
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Find bounds for the plot (based on all positions)
    all_x = Locations[:, :, 0].flatten()
    all_y = Locations[:, :, 1].flatten()
    all_z = Locations[:, :, 2].flatten()

    max_range = max(all_x.max() - all_x.min(),
                    all_y.max() - all_y.min(),
                    all_z.max() - all_z.min()) / 2

    mid_x = (all_x.max() + all_x.min()) / 2
    mid_y = (all_y.max() + all_y.min()) / 2
    mid_z = (all_z.max() + all_z.min()) / 2

    # Set axis limits
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Gravity Simulation')

    # Create a scatter plot for each body
    colors = plt.cm.rainbow(np.linspace(0, 1, n_bodies))
    scatters = []
    trails = []

    for i in range(n_bodies):
        # Current position (will be updated)
        scatter = ax.scatter([], [], [], c=[colors[i]], s=100, marker='o')
        scatters.append(scatter)

        # Trail (path history)
        trail, = ax.plot([], [], [], c=colors[i], alpha=0.3, linewidth=1)
        trails.append(trail)

    def init():
        """Initialize animation"""
        for scatter in scatters:
            scatter._offsets3d = ([], [], [])
        for trail in trails:
            trail.set_data([], [])
            trail.set_3d_properties([])
        return scatters + trails

    def update(frame):
        """Update animation at each frame"""
        frame = frame * skip_frames  # Skip frames for speed

        if frame >= n_steps:
            frame = n_steps - 1

        # Update each body
        for i in range(n_bodies):
            # Update current position
            x, y, z = Locations[frame, i]
            scatters[i]._offsets3d = ([x], [y], [z])

            # Update trail (last 500 points or less)
            trail_start = max(0, frame - 500)
            trail_x = Locations[trail_start:frame + 1, i, 0]
            trail_y = Locations[trail_start:frame + 1, i, 1]
            trail_z = Locations[trail_start:frame + 1, i, 2]

            trails[i].set_data(trail_x, trail_y)
            trails[i].set_3d_properties(trail_z)

        # Update title with current time
        ax.set_title(f'Gravity Simulation - Time: {frame * dt:.2f}s')

        return scatters + trails

    # Create animation
    n_frames = n_steps // skip_frames
    anim = FuncAnimation(fig, update, frames=n_frames,
                         init_func=init, blit=False,
                         interval=20, repeat=True)

    plt.show()

    return anim


# Run the animation
anim = animate_simulation(Locations, dt, skip_frames=100)





