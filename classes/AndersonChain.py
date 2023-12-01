import numpy as np
from scipy.linalg import expm

class AndersonChain:
    '''
    Construct and simulate a 1D Anderson lattice with periodic boundary conditions.
    We represent all wave functions and operators in the lattice site basis. 

    Attributes
        chain (1D ndarray): 
        psi0 (1D ndarray): 
        eps_range (array-like):
        t_hop (float):
        story_history (boolean):
    '''
    def __init__(self, num_sites, psi0, eps_range, t_hop, store_history=True):
        self.chain = np.zeros((num_sites, 1))
        self.psi0 = psi0
        self.eps_range = eps_range
        self.t_hop = t_hop # hopping param

        if store_history:
            self.history = []

    def _hamiltonian(self):
        '''
        Construct the hamiltonian for the Anderson tight-binding model, a matrix representation in the occupancy site basis. 
        '''
        binding = np.diagflat(np.random.uniform(*self.eps_range, size=(len(self.chain, ))))
        
        hopping = self.t_hop * (np.roll(np.identity(len(self.chain)), 1, axis=0) + np.roll(np.identity(len(self.chain)), -1, axis=0))
        
        return binding + hopping
    
    
    def _time_evolution(self, time):
        '''
        Calculate the unitary time evolution operator for the given hamiltonian.

        Args
            time (float): time when time evolution operator is calculated

        Returns
            U(t) (ndarray of size num_sites x num_sites) 
        '''

        return expm(-1j * self._hamiltonian() * time)
    

    def solveatt(self, time):
       
       return self._time_evolution(time) @ self.psi0
    

    def solve(self, t, nt): #t_steps):
        '''
        Calculate psi(t). 

        Args
            t (array-like): time range of form (t_initial, t_final)
            nt (int): number of time steps
        '''        

        times = np.linspace(*t, nt)
        
        history = []

        for time in times:
            psi_t = self._time_evolution(time) @ self.psi0
            history.append(psi_t) 
        
        return history
    