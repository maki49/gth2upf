from element_list import ElementList
from upf_data import UPFData, write_upf_v2, read_upf_file, rwtest
from gth_data import GTHData, read_gth_file
from grid import Grid
import numpy as np
from math import exp, log, sqrt, erf, gamma

E2 = 2.0  # Hartree to Rydberg conversion factor
def mygamma(n):
    if n<2: 
        print(f"Invalid argument {n} to mygamma = Gamma(n-1/2)")
        exit(1)
    else:
        return gamma(n-0.5)

class PSWfc:
    def __init__(self, chi=None, oc = None, lchi = None):
        self.chi = chi
        self.nwfc = 0 if chi is None else len(chi)
        self.oc = oc
        self.lchi = lchi


def gth2upf_custom(gth: GTHData, info="", 
                pswfc: PSWfc = PSWfc(),
                grid: Grid = None,
                nlcc=False, rho_atc=[])->UPFData:
    """Convert GTH data to UPF format
    input: 
        gth: GTHData object
        info: string, info section
        pswfc: PSWFC object
        grid: Grid object
        nlcc and rho_atc
    """
    if grid is None:
        grid = Grid(zmesh=gth.z)
    upf = UPFData()

    # Print info section
    upf.info = info
    
    # Get user input for parameters
    lloc = -3   #???
    rcloc = 0.0 #???
    # Set UPF header
    upf.generated = "Generated in analytical, separable form"
    upf.author = "Goedecker/Hartwigsen/Hutter/Teter"
    upf.date = "Phys.Rev.B58, 3641 (1998); B54, 1703 (1996)"
    upf.nv = "2.0.1"
    upf.comment = "Info: automatically converted from CPMD format"
    upf.psd = gth.atom_name
    upf.typ = 'NC'
    
    if gth.z > 18:
        upf.rel = 'no'
    else:
        upf.rel = 'scalar'

    upf.tvanp = False
    upf.tpawp = False
    upf.tcoulombp = False
    upf.nlcc = nlcc
    
    # Set DFT functional
    upf.dft = gth.xc

    # Set basic parameters
    upf.zp = gth.zv
    upf.etotps = 0.0
    upf.ecutrho = 0.0
    upf.ecutwfc = 0.0
    
    if gth.lmax == lloc:
        upf.lmax = gth.lmax - 1
    else:
        upf.lmax = gth.lmax
    
    upf.lloc = lloc
    upf.lmax_rho = 0
    upf.nwfc = pswfc.nwfc
    
    # Set wavefunction info
    upf.els = []
    upf.oc = []
    upf.epseu = []
    upf.lchi = []
    upf.nchi = []
    upf.rcut_chi = []
    upf.rcutus_chi = []
    
    # spdf = ['S', 'P', 'D', 'F']
    upf.oc = pswfc.oc
    upf.lchi = pswfc.lchi
    upf.rcht_chi = np.zeros(pswfc.nwfc)
    upf.rcutus_chi = np.zeros(pswfc.nwfc)
    upf.epseu = np.zeros(pswfc.nwfc)
    
    for i in range(upf.nwfc):
        upf.lchi.append(l)
        upf.rcut_chi.append(0.0)
        upf.rcutus_chi.append(0.0)
        upf.epseu.append(0.0)
    
    # Set grid
    upf.mesh = grid.mesh
    upf.dx = grid.dx
    upf.rmax = grid.rmax
    upf.xmin = grid.xmin
    upf.zmesh = gth.z
    
    upf.r = grid.r
    upf.rab = grid.rab
    
    # Core correction
    if upf.nlcc:
        upf.rho_atc = rho_atc
    else:
        upf.rho_atc = np.zeros(upf.mesh)
    
    # Local potential
    upf.rcloc = rcloc
    upf.vloc = np.zeros(upf.mesh)
    exp_loc = np.zeros(4)
    exp_loc[0:len(gth.exp_loc)] = gth.exp_loc
    for i in range(upf.mesh):
        x = upf.r[i] / gth.rc
        x2 = x**2
        upf.vloc[i] = E2 * (-upf.zp * erf(x/sqrt(2.0)) / upf.r[i] +
                            exp(-0.5 * x2) * (exp_loc[0] + x2 * (exp_loc[1] + 
                            x2 * (exp_loc[2] + x2 * exp_loc[3]))))
    
    upf.nbeta = 0
    for l in range(upf.lmax+1):
        upf.nbeta += gth.nh_nl[l]
    # Beta functions and projectors
    if upf.nbeta > 0:
        upf.els_beta = []
        upf.lll = []
        upf.kbeta = []
        upf.beta = np.zeros((upf.nbeta, upf.mesh))
        upf.dion = np.zeros((upf.nbeta, upf.nbeta))
        upf.rcut = []
        upf.rcutus = []
        
        for l in range(upf.lmax+1):
            ij = 0
            iv =0
            for i in range(gth.nh_nl[l]):
                upf.lll.append(l)
                # upf.els_beta.append(f"{i+1}{spdf[l]}")
                
                for j in range(i, gth.nh_nl[l]):
                    jv = iv + j - i
                    upf.dion[iv, jv] = gth.h[l][ij] / E2
                    ij += 1
                    if j > i:
                        upf.dion[jv, iv] = upf.dion[iv, jv]
                
                fac = sqrt(2.0 * gth.rad_nl[l]) / (gth.rad_nl[l]**(l + 2*(i+1)) * sqrt(mygamma(l + 2*(i+1))))
                
                for ir in range(upf.mesh):
                    x2 = (upf.r[ir] / gth.rad_nl[l])**2
                    upf.beta[iv][ir] = (upf.r[ir]**(l + 2*i) * 
                                        exp(-0.5 * x2) * fac * E2 * upf.r[ir])
                
                # Find kbeta
                kbeta = upf.mesh
                for ir in range(upf.mesh-1, -1, -1):
                    if abs(upf.beta[iv][ir]) > 1e-12:
                        kbeta = ir + 1
                        break
                if kbeta < 2:
                    raise ValueError(f"Zero beta function {iv}")
                elif kbeta % 2 == 0 and kbeta < upf.mesh:
                    kbeta += 1
                upf.kbeta.append(min(upf.mesh, kbeta))
                upf.rcut.append(upf.r[upf.kbeta[iv]-1])
                upf.rcutus.append(0.0)
                
                iv += 1
        
        upf.kkbeta = max(upf.kbeta) if upf.kbeta else 0
        

    # Atomic wavefunctions
    if pswfc.chi is not None:
        upf.chi = pswfc.chi
    else:
        upf.chi = np.zeros((upf.mesh, upf.nwfc))
    
    # Atomic density
    upf.rho_at = np.zeros(upf.mesh)
    for i in range(upf.nwfc):
        upf.rho_at += upf.oc[i] * upf.chi[:, i]**2
    
    print("Pseudopotential successfully converted")
    
    return upf

if __name__=="__main__":
    
    def compare_mesh(m1:Grid, m2:Grid):
        return (m1.nr == m2.nr and m1.dx==m2.dx and m1.rmax==m2.rmax and m1.xmin == m2.xmin) \
            and (np.allclose(m1.r, m2.r) and np.allclose(m1.rab, m2.rab))
        
    
    import sys
    gth_file_read = sys.argv[1]
    upf_file_ref = sys.argv[2]
    print("read file:", gth_file_read)
    gth = read_gth_file(gth_file_read)
    upf = gth2upf_custom(gth)
    
    upf_ref = read_upf_file(upf_file_ref)
    rwtest(upf, upf_ref)