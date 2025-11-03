from math import log, exp
import numpy as np
from dataclasses import dataclass

@dataclass
class Grid:
    zmesh: int = 1
    r: np.ndarray = None
    mesh: int = 0   # number of radial grid points
    dx: float = 0.0
    rab: np.ndarray = None
    rmax: float = 0.0
    xmin: float = 0.0

class GridCPMD2UPF(Grid): 
    def __init__(self, zmesh=1, xmin=-7.00, rmax_init=100.0, dx=0.0125):
        super().__init__(zmesh=zmesh, rmax=rmax_init, dx=dx, xmin=xmin)
        self.rmax_init = rmax_init

        self.gen_grid(self.zmesh)

        self.rab = self.r * self.dx
        self.rmax = self.r[-1]
        self.xmin = log(self.zmesh * self.r[0])
    
    def gen_grid(self, zmesh):
        "r_i = 1/zmesh * exp(xmin + i * dx)"
        self.mesh = 1 + int((log(zmesh * self.rmax_init) - self.xmin) / self.dx)
        self.mesh = (self.mesh // 2) * 2 + 1  # Make odd
        self.r = np.zeros(self.mesh)
        for i in range(self.mesh):
            self.r[i] = exp(self.xmin + i * self.dx) / zmesh

