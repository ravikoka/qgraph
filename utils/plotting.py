import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

plt.style.use('dark_background')


def plot_lattice_pdf(anderson_lattice, time, fix_z_scale=True):
    '''
    Plot the PDF of the wavefunction on a 2D periodic square lattice (torus).
    
    Args
        anderson_lattice (AndersonGraph): 2D anderson lattice (ie. a torus)
        time (float): time at which the PDF is calculated. 
        fix_z_scale (bool): if True, fix the height of the z axis to 1.
    '''
    
    psi_t = anderson_lattice.psi_at_t(time)
    density = np.real(np.multiply(psi_t.conj(), psi_t))
    n = int(np.sqrt(len(density)))
    density = np.reshape(density, (n, n))

    fig = plt.figure(figsize = (12,12))
    ax = plt.axes(projection='3d')
    x = np.arange(0, n, 1)
    y = np.arange(0, n, 1)
    
    X, Y = np.meshgrid(x, y)

    ax.plot_surface(X, Y, density, cmap = plt.cm.cividis)

    if fix_z_scale:
        ax.set_zlim((0, 1))#np.max(density)))
        #ax.set_zlim((0, max(maxy, 0.2))) why 0.2???

    ax.set_title("Wave function probability density at time " + str(time))
    ax.set_xlabel('x', labelpad=20)
    ax.set_ylabel('y', labelpad=20)
    ax.set_zlabel('Probability Density', labelpad=20)

    plt.show()
