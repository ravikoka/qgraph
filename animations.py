import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation
from classes.AndersonGraph import AndersonGraph

plt.style.use('dark_background')



def animate_lattice_pdf(anderson_lattice, t_max, W):
    
    n = int(np.sqrt(anderson_lattice.num_sites))

    fig = plt.figure(figsize = (12,10))
    
    x = np.arange(0, n, 1)
    y = np.arange(0, n, 1)
    X, Y = np.meshgrid(x, y)

    ax = plt.axes(projection='3d')
    ax.set_xlabel('x', labelpad=20)
    ax.set_ylabel('y', labelpad=20)
    ax.set_zlabel('P', labelpad=20)

    def update(frame):
        ax.clear()
        psi_t = anderson_lattice.psi_at_t(frame)
        density = np.real(np.multiply(psi_t.conj(), psi_t))
        density = np.reshape(density, (n, n))
        
        ax.plot_surface(X, Y, density, cmap = plt.cm.cividis)
        ax.set_zlim(0, 1)
        
        ax.set_title(f'Wave function probability density at time {frame}')
    
    ani = FuncAnimation(fig, update, frames=range(0, t_max), interval=100)
    
    path = f'plots/lattice_animation_size_{n}_tmax_{t_max}_W_{W}.gif'
    
    #'torusanimation_n_' + str(n) + '_tmax_' + str(t_max)+ '_W_' + str(W) + '.gif'
    ani.save(path, writer='ffmpeg', fps=5)
    
    plt.show()
    