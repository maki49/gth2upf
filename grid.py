from math import log, exp
import numpy as np

class Grid:
    def __init__(self, zmesh=1, r=None):
        self.zmesh = zmesh
        self.r = np.zeros(0)
        self.mesh = 0   # number of radial grid points
        self.dx = 0.0
        self.rab = np.zeros(0)
        self.rmax = 0.0
        self.xmin = 0.0
        if r is not None and len(r) > 1:
            self.r = r
            self.mesh = len(r)
            self.dx = log(self.r[1]/self.r[0])
        else:
            self.default_grid_cpmd2upf(self.zmesh)
        self.rab = self.r * self.dx
        self.rmax = self.r[-1]
        self.xmin = log(self.zmesh * self.r[0])
                
    # Check grid consistency
    def check_grid_consistency(self):
        for i in range(1, min(10, self.mesh)):
            expected = self.r[0] * exp(i * self.dx)
            if abs(self.r[i] - expected) > 1e-6 * expected:
                print(f"Warning: grid point {i} deviates from logarithmic")
                break
            
    # Default grid in cpmd2upf
    def default_grid_cpmd2upf(self, z:int):
        self.xmin = -7.0
        self.dx = 0.0125
        self.rmax = 100.0
        print("Using default radial grid:")
        print(f"r_i = exp(xmin+(i-1)*dx)/Z, with Z={z}, xmin={self.xmin}, dx={self.dx}, rmax={self.rmax}")
        
        self.mesh = 1 + int((log(z * self.rmax) - self.xmin) / self.dx)
        self.mesh = (self.mesh // 2) * 2 + 1  # Make odd
        self.r = np.zeros(self.mesh)
        for i in range(self.mesh):
            self.r[i] = exp(self.xmin + i * self.dx) / z        
        print(f"{self.mesh} grid points, rmax={self.r[-1]:.4f}")
    
