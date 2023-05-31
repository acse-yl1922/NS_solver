import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Define the file paths and names
file_prefix = 'C:/Users/yongqi/Desktop/dev/mpi-assessment/x64/Debug/out/'
file_suffix = '.dat'

p_files = ['P_{}.dat'.format(i) for i in range(101)]
u_files = ['u_{}.dat'.format(i) for i in range(101)]
v_files = ['v_{}.dat'.format(i) for i in range(101)]

# Define the animation update function
def update(frame):
    p_data = np.loadtxt(file_prefix + p_files[frame])
    u_data = np.loadtxt(file_prefix + u_files[frame])
    v_data = np.loadtxt(file_prefix + v_files[frame])
    
    plt.clf()
    plt.subplot(131)
    plt.imshow(p_data, cmap='jet')
    plt.title('P_{}.dat'.format(frame))
    plt.colorbar()

    plt.subplot(132)
    plt.imshow(u_data, cmap='jet')
    plt.title('u_{}.dat'.format(frame))
    plt.colorbar()

    plt.subplot(133)
    plt.imshow(v_data, cmap='jet')
    plt.title('v_{}.dat'.format(frame))
    plt.colorbar()

# Create the animation
fig = plt.figure(figsize=(12, 4))
ani = animation.FuncAnimation(fig, update, frames=len(p_files), interval=200)

# Save the animation as a GIF
output_path = 'C:/Users/yongqi/Desktop/dev/mpi-assessment/x64/Debug/animation.gif'
ani.save(output_path, writer='pillow', fps=5)

# Display the animation
plt.show()













