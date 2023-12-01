import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm

class AndersonLattice:
    '''
    Construct and simulate a 2D Anderson lattice with periodic boundary conditions.
    We represent all wave functions and operators in the lattice site basis. 

    Attributes
        chain (2D ndarray NxN): 
        psi0 (1D ndarray N^2): 
        eps_range (array-like):
        t_hop (float):
        store_history (boolean):
    '''
    def __init__(self, num_sites, psi_0, eps_range, t_hop, store_history=True):
        self.chain = np.zeros((num_sites, num_sites))
        self.num_sites = num_sites
        self.psi_0 = psi_0
        self.eps_range = eps_range
        self.t_hop = t_hop # hopping param
        self.binding = np.diagflat(np.random.uniform(*self.eps_range, size=(self.num_sites**2)))

        if store_history:
            self.history = []

    def _hamiltonian(self):
        '''
        Construct the hamiltonian for the Anderson tight-binding model, a matrix representation in the occupancy site basis. 
        '''
        G = nx.grid_graph(dim = (self.num_sites, self.num_sites), periodic=True)
        A = nx.to_numpy_array(G)
        hopping = -1*self.t_hop * A

        return self.binding + hopping
    
    
    def _time_evolution(self, time):
        '''
        Calculate the unitary time evolution operator for the given hamiltonian.

        Args
            time (float): time when time evolution operator is calculated

        Returns
            U(t) (ndarray of size num_sites x num_sites) 
        '''

        return expm(-1j * self._hamiltonian() * time)
    
    def solveatt(self, t):
        return self._time_evolution(t) @ self.psi_0
    
    
    def solve(self, t, nt): #t_steps):
        '''
        Calculate psi(t). 

        Args
            t (array-like): time range of form (t_initial, t_final)
            nt (int): number of time steps
        '''        

        times = np.linspace(0, t, nt)
        
        history = []

        for time in times:
            psi_t = self._time_evolution(time) @ self.psi_0
            history.append(psi_t) 
        
        return history
    
    def plotdensity(self, t, axisstabilized = False):
        '''
        Plot the probability density, aka |psi(t)|^2

        Args:
            t (float): time to plot it
        '''

        psi_t = self._time_evolution(t) @ self.psi_0
        density = np.real(np.multiply(psi_t.conj(), psi_t))
        dens2d = np.reshape(density, (self.num_sites, self.num_sites))

        fig = plt.figure(figsize = (12,12))
        plt.style.use('dark_background')
        ax = plt.axes(projection='3d')
        x = np.arange(0, self.num_sites, 1)
        y = np.arange(0, self.num_sites, 1)
        X, Y = np.meshgrid(x, y)

        surf = ax.plot_surface(X, Y, dens2d, cmap = plt.cm.cividis)
        
        # Set axes label
        ax.set_xlabel('x', labelpad=20)
        ax.set_ylabel('y', labelpad=20)
        ax.set_zlabel('P', labelpad=20)
        if axisstabilized==True:
            maxy = max(dens2d.flatten())
            ax.set_zlim((0, max(maxy, 0.2)))
        ax.set_title("Wave function probability density at time " + str(t))

        #plt.colorbar(surf)

        plt.show()
        






    