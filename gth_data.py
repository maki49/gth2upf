from element_list import ElementList

# Literature: - S. Goedecker, M. Teter, and J. Hutter,
#               Phys. Rev. B 54, 1703 (1996)
#             - C. Hartwigsen, S. Goedecker, and J. Hutter,
#               Phys. Rev. B 58, 3641 (1998)
#             - M. Krack,
#               Theor. Chem. Acc. 114, 145 (2005)
#
# GTH-potential format:
#
# Element symbol  Name of the potential  Alias names
# n_elec(s)  n_elec(p)  n_elec(d)  ...
# r_loc   nexp_ppl        cexp_ppl(1) ... cexp_ppl(nexp_ppl)
# nprj
# r(1)    nprj_ppnl(1)    ((hprj_ppnl(1,i,j),j=i,nprj_ppnl(1)),i=1,nprj_ppnl(1))
# r(2)    nprj_ppnl(2)    ((hprj_ppnl(2,i,j),j=i,nprj_ppnl(2)),i=1,nprj_ppnl(2))
#  .       .               .
#  .       .               .
#  .       .               .
# r(nprj) nprj_ppnl(nprj) ((hprj_ppnl(nprj,i,j),j=i,nprj_ppnl(nprj)),
#                                               i=1,nprj_ppnl(nprj))
#
# n_elec   : Number of electrons for each angular momentum quantum number
#            (electronic configuration -> s p d ...)
# r_loc    : Radius for the local part defined by the Gaussian function
#            exponent alpha_erf
# nexp_ppl : Number of the local pseudopotential functions
# cexp_ppl : Coefficients of the local pseudopotential functions
# nprj     : Number of the non-local projectors => nprj = SIZE(nprj_ppnl(:))
# r        : Radius of the non-local part for angular momentum quantum number l
#            defined by the Gaussian function exponents alpha_prj_ppnl
# nprj_ppnl: Number of the non-local projectors for the angular momentum
#            quantum number l
# hprj_ppnl: Coefficients of the non-local projector functions
class GTHData:
    def __init__(self):
        self.atom_name = ""
        self.z = 0
        self.xc = ""
        self.zv = 0
        self.nelec_of_l = []
        self.lmax = 0
        self.rc = 0.0    # the characteristic radius parameter for the local part
        self.nexp_loc = 0 # number of exponential terms in the local potential
        self.exp_loc = [] # coefficients of the local pseudopotential functions
        self.nprj = 0    # number of angular momentum channels fo projectors
        self.rad_nl = [] # radius of the non-local part for each angular momentum channel l
        self.nh_nl = []  # number of non-local projectors for each angular momentum channel l
        self.h = []      # coefficients of the non-local projector functions for each angular momentum channel l   

def read_gth(lines) -> GTHData:
    """Read GTH pseudopotential of CP2K format"""
    gth_data=GTHData()
    i=0
    while i < len(lines):
        line = lines[i].strip()
        
        if i==0: # the first line
            ls = line.replace('-', ' ').split()
            gth_data.atom_name = ls[0]
            gth_data.z=ElementList().name2index(ls[0])
            gth_data.xc=ls[2]
            gth_data.zv=int(ls[3][1:])
        
        elif i==1:
            ls=line.split()
            gth_data.nelec_of_l=[int(s) for s in ls]
            gth_data.lmax=len(ls)-1
            
        elif i==2:
            ls=line.split()
            gth_data.rc = float(ls[0])    # the characteristic radius parameter for the local part
            gth_data.nexp_loc=int(ls[1])  # number of exponential terms in the local potential
            gth_data.exp_loc=[float(s) for s in ls[2:]]
            
        elif i==3:
            gth_data.nprj=int(line.strip()) #Number of angular momentum channels fo projectors
            lines=lines[4:]
            for l in range(gth_data.nprj):
                ls=lines[0].split()
                gth_data.rad_nl.append(float(ls[0]))
                nh = int(ls[1])
                gth_data.nh_nl.append(nh)
                gth_data.h.append([ float(s) for s in ls[2:] ])
                if nh > 1 :
                    for h in range(nh-1):
                        gth_data.h[l].extend([ float(s) for s in lines[h+1].split() ])
                assert(len(gth_data.h[l])==(nh+1)*nh/2)
                lines=lines[max(1,nh):]
                
        else: 
            print(f"Error: unexpected line at index {i}: {line}")
            continue
        i +=1
    return gth_data
        
def read_gth_file(filename: str) -> GTHData:
    """Read GTH pseudopotential of CP2K format"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    return read_gth(lines)


def write_gth_file(gth_data: GTHData, filename: str):
    """Write GTH pseudopotential of CP2K format"""
    with open(filename, 'w') as f:
        f.write(f"{gth_data.atom_name} GTH-{gth_data.xc}-q{gth_data.zv}\n")
        f.write(" ".join([str(s) for s in gth_data.nelec_of_l]) + "\n")
        f.write(f"{gth_data.rc} {gth_data.nexp_loc} " + " ".join([str(s) for s in gth_data.exp_loc]) + "\n")
        f.write(f"{gth_data.nprj}\n")
        for l in range(gth_data.nprj):
            f.write(f"{gth_data.rad_nl[l]} {gth_data.nh_nl[l]}\t")
            nh_now = gth_data.nh_nl[l]
            h_now = gth_data.h[l]
            while nh_now > 0:
                f.write("\t".join([str(s) for s in h_now[:nh_now]]) + "\n")
                h_now = h_now[nh_now:]
                nh_now-=1
                    
                    
if __name__ == "__main__":
    """For read-and-write (RW) test, usage: python gth_data.py <gth_file_in> <gth_file_out>"""
    import sys, os
    rfile = sys.argv[1]
    wfile = sys.argv[2]
    print("read file:", rfile)
    gth_ref = read_gth_file(rfile)
    write_gth_file(gth_ref, wfile)
    assert(os.path.exists(wfile))
    gth_wr = read_gth_file(wfile)
    # diff=os.popen(f'diff {rfile} {wfile}').read()
    if (gth_ref.__dict__==gth_wr.__dict__):
        print('RW test passed')
    else:
        print(f'RW test faild, ref=', gth_ref.__dict__)
        print(f'result=', gth_wr.__dict__)
    os.remove(wfile)