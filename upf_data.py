class UPFData:
    """Class to hold UPF pseudopotential data"""
    
    def __init__(self):
        # Info
        self.info=""
        # Header
        self.generated = ""
        self.author = ""
        self.date = ""
        self.nv = "2.0.1"
        self.comment = ""
        self.psd = ""
        self.rel = ""
        self.typ = ""
        self.tvanp = False
        self.tpawp = False
        self.tcoulombp = False
        self.nlcc = False
        self.dft = ""
        
        # Basic parameters
        self.zp = 0.0
        self.etotps = 0.0
        self.ecutrho = 0.0
        self.ecutwfc = 0.0
        self.lmax = 0
        self.lloc = 0
        self.lmax_rho = 0
        self.nwfc = 0
        
        # Wavefunction info
        self.els = []
        self.oc = []
        self.epseu = []
        self.lchi = []
        self.nchi = []
        self.rcut_chi = []
        self.rcutus_chi = []
        
        # Grid
        self.mesh = 0
        self.dx = 0.0
        self.rmax = 0.0
        self.xmin = 0.0
        self.zmesh = 0.0
        self.rab = None
        self.r = None
        
        # Potentials and densities
        self.rho_atc = None
        self.rcloc = 0.0
        self.vloc = None
        self.vnl = None
        
        # Beta functions
        self.nbeta = 0
        self.els_beta = []
        self.lll = []
        self.kbeta = []
        self.kkbeta = 0
        self.beta = None
        self.dion = None
        self.rcut = []
        self.rcutus = []
        
        # Atomic wavefunctions and density
        self.chi = None
        self.rho_at = None

import re
import numpy as np

def read_upf_file(filename):
    """Read UPF pseudopotential file and return UPFData object
    
    Args:
        filename (str): Path to UPF file
        
    Returns:
        UPFData: Parsed UPF data
    """
    upf_data = UPFData()
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Parse UPF version
    version_match = re.search(r'<UPF version="([^"]+)">', content)
    if version_match:
        upf_data.nv = version_match.group(1)
        
    # Parse PP_INFO section
    info_match = re.search(r'<PP_INFO[^>]*>(.*?)</PP_INFO>', content, re.DOTALL)
    if info_match:
        upf_data.info =  info_match.group(1).strip()
        
    
    # Parse PP_HEADER section
     # [^>] matches any character except >
     # re.DOTALL makes . match newlines as well
    header_match = re.search(r'<PP_HEADER([^>]*)>', content, re.DOTALL)
    if header_match:
        header_attrs = header_match.group(1)
        
        # Extract attributes using regex
        attr_patterns = {
            'generated': r'generated="([^"]*)"',
            'author': r'author="([^"]*)"',
            'date': r'date="([^"]*)"',
            'comment': r'comment="([^"]*)"',
            'element': r'element="([^"]*)"',
            'pseudo_type': r'pseudo_type="([^"]*)"',
            'relativistic': r'relativistic="([^"]*)"',
            'is_ultrasoft': r'is_ultrasoft="([^"]*)"',
            'is_paw': r'is_paw="([^"]*)"',
            'is_coulomb': r'is_coulomb="([^"]*)"',
            'has_so': r'has_so="([^"]*)"',
            'has_wfc': r'has_wfc="([^"]*)"',
            'has_gipaw': r'has_gipaw="([^"]*)"',
            'paw_as_gipaw': r'paw_as_gipaw="([^"]*)"',
            'core_correction': r'core_correction="([^"]*)"',
            'functional': r'functional="([^"]*)"',
            'z_valence': r'z_valence="([^"]*)"',
            'total_psenergy': r'total_psenergy="([^"]*)"',
            'wfc_cutoff': r'wfc_cutoff="([^"]*)"',
            'rho_cutoff': r'rho_cutoff="([^"]*)"',
            'l_max': r'l_max="([^"]*)"',
            'l_max_rho': r'l_max_rho="([^"]*)"',
            'l_local': r'l_local="([^"]*)"',
            'mesh_size': r'mesh_size="([^"]*)"',
            'number_of_wfc': r'number_of_wfc="([^"]*)"',
            'number_of_proj': r'number_of_proj="([^"]*)"'
        }
        
        for attr, pattern in attr_patterns.items():
            match = re.search(pattern, header_attrs)
            if match:
                value = match.group(1).strip()
                if attr == 'generated':
                    upf_data.generated = value
                elif attr == 'author':
                    upf_data.author = value
                elif attr == 'date':
                    upf_data.date = value
                elif attr == 'comment':
                    upf_data.comment = value
                elif attr == 'element':
                    upf_data.psd = value.strip()
                elif attr == 'pseudo_type':
                    upf_data.typ = value
                elif attr == 'relativistic':
                    upf_data.rel = value
                elif attr == 'is_ultrasoft':
                    upf_data.tvanp = value.upper() == 'T'
                elif attr == 'is_paw':
                    upf_data.tpawp = value.upper() == 'T'
                elif attr == 'is_coulomb':
                    upf_data.tcoulombp = value.upper() == 'T'
                elif attr == 'core_correction':
                    upf_data.nlcc = value.upper() == 'T'
                elif attr == 'functional':
                    upf_data.dft = value
                elif attr == 'z_valence':
                    upf_data.zp = float(value)
                elif attr == 'total_psenergy':
                    upf_data.etotps = float(value)
                elif attr == 'wfc_cutoff':
                    upf_data.ecutwfc = float(value)
                elif attr == 'rho_cutoff':
                    upf_data.ecutrho = float(value)
                elif attr == 'l_max':
                    upf_data.lmax = int(value)
                elif attr == 'l_max_rho':
                    upf_data.lmax_rho = int(value)
                elif attr == 'l_local':
                    upf_data.lloc = int(value)
                elif attr == 'mesh_size':
                    upf_data.mesh = int(value)
                elif attr == 'number_of_wfc':
                    upf_data.nwfc = int(value)
                elif attr == 'number_of_proj':
                    upf_data.nbeta = int(value)
    
    # Parse PP_MESH section
    mesh_match = re.search(r'<PP_MESH([^>]*)>', content)
    if mesh_match:
        mesh_attrs = mesh_match.group(1)
        
        dx_match = re.search(r'dx="([^"]+)"', mesh_attrs)
        if dx_match:
            upf_data.dx = float(dx_match.group(1))

        xmin_match = re.search(r'xmin="([^"]+)"', mesh_attrs)
        if xmin_match:
            upf_data.xmin = float(xmin_match.group(1))
            
        rmax_match = re.search(r'rmax="([^"]+)"', mesh_attrs)
        if rmax_match:
            upf_data.rmax = float(rmax_match.group(1))
            
        zmesh_match = re.search(r'zmesh="([^"]+)"', mesh_attrs)
        if zmesh_match:
            upf_data.zmesh = float(zmesh_match.group(1))
    
    # Parse PP_R section (radial grid)
    r_match = re.search(r'<PP_R[^>]*>([^<]+)</PP_R>', content)
    if r_match:
        r_data = r_match.group(1).strip()
        r_values = []
        for line in r_data.split('\n'):
            if line.strip():
                r_values.extend([float(x) for x in line.split()])
        upf_data.r = np.array(r_values)

    # Parse PP_RAB section (dr/di)
    rab_match = re.search(r'<PP_RAB[^>]*>([^<]+)</PP_RAB>', content)
    if rab_match:
        rab_data = rab_match.group(1).strip()
        rab_values = []
        for line in rab_data.split('\n'):
            if line.strip():
                rab_values.extend([float(x) for x in line.split()])
        upf_data.rab = np.array(rab_values)

    
    # Parse PP_LOCAL section (local potential)
    local_match = re.search(r'<PP_LOCAL[^>]*>([^<]+)</PP_LOCAL>', content)
    if local_match:
        local_data = local_match.group(1).strip()
        local_values = []
        for line in local_data.split('\n'):
            if line.strip():
                local_values.extend([float(x) for x in line.split()])
        upf_data.vloc = np.array(local_values)
    
    # Parse PP_NONLOCAL section (beta functions and D matrix)
    nonlocal_match = re.search(r'<PP_NONLOCAL>(.*?)</PP_NONLOCAL>', content, re.DOTALL)
    if nonlocal_match:
        nonlocal_content = nonlocal_match.group(1)
        # Parse angular momentums
        l_matches = re.findall(r'angular_momentum="([^"]+)"', nonlocal_content)
        if l_matches:
            upf_data.lll = np.array([int(l) for l in l_matches])
            assert(len(upf_data.lll) == upf_data.nbeta)
            assert(max(upf_data.lll) <= upf_data.lmax)
        # Parse beta functions
        beta_matches = re.findall(r'<PP_BETA[^>]*>([^<]+)</PP_BETA[^>]*>', nonlocal_content, re.DOTALL)
        if beta_matches:
            upf_data.beta = []
            for beta_data in beta_matches:
                beta_values = []
                for line in beta_data.strip().split('\n'):
                    if line.strip():
                        beta_values.extend([float(x) for x in line.split()])
                upf_data.beta.append(np.array(beta_values))
        # Parse D matrix
        dij_match = re.search(r'<PP_DIJ[^>]*>([^<]+)</PP_DIJ>', nonlocal_content, re.DOTALL)
        if dij_match:
            dij_data = dij_match.group(1).strip()
            dij_values = []
            for line in dij_data.split('\n'):
                if line.strip():
                    dij_values.extend([float(x) for x in line.split()])
            # Reshape to square matrix
            n = int(np.sqrt(len(dij_values)))
            upf_data.dion = np.array(dij_values).reshape(n, n)
    if upf_data.nbeta > 0:
        upf_data.kbeta=np.ones(upf_data.nbeta, dtype=int)*upf_data.mesh
        
    # Parse PP_PSWFC section (pseudo wavefunctions)
    pswfc_match = re.search(r'<PP_PSWFC>(.*?)</PP_PSWFC>', content, re.DOTALL)
    if pswfc_match:
        pswfc_content = pswfc_match.group(1)
        # Parse angular momentums
        l_matches = re.findall(r' l="([^"]+)"', pswfc_content)  # avoid "label="
        if l_matches:
            upf_data.lchi = np.array([int(l) for l in l_matches])
            assert(len(upf_data.lchi) == upf_data.nwfc)
            # assert(max(upf_data.lchi) <= upf_data.lmax)
        # Parse occupations
        oc_matches = re.findall(r'occupation="([^"]+)"', pswfc_content)
        if oc_matches:
            upf_data.oc = np.array([float(oc) for oc in oc_matches])
            assert(len(upf_data.oc) == upf_data.nwfc)
        # Parse pseudoenergies
        epseu_matches = re.findall(r'pseudo_energy="([^"]+)"', pswfc_content)
        if epseu_matches:
            upf_data.epseu = np.array([float(epseu) for epseu in epseu_matches])
            assert(len(upf_data.epseu) == upf_data.nwfc)
        # Parse individual wavefunctions
        wfc_matches = re.findall(r'<PP_CHI[^>]*>([^<]+)</PP_CHI[^>]*>', pswfc_content, re.DOTALL)
        if wfc_matches:
            upf_data.chi = []
            for wfc_data in wfc_matches:
                wfc_values = []
                for line in wfc_data.strip().split('\n'):
                    if line.strip():
                        wfc_values.extend([float(x) for x in line.split()])
                upf_data.chi.append(np.array(wfc_values))
    
    # Parse PP_RHOATOM section (atomic charge density)
    rhoatom_match = re.search(r'<PP_RHOATOM[^>]*>([^<]+)</PP_RHOATOM>', content, re.DOTALL)
    if rhoatom_match:
        rhoatom_data = rhoatom_match.group(1).strip()
        rhoatom_values = []
        for line in rhoatom_data.split('\n'):
            if line.strip():
                rhoatom_values.extend([float(x) for x in line.split()])
        upf_data.rho_at = np.array(rhoatom_values)
    
    # Parse PP_NLCC section (nonlinear core correction) if present
    nlcc_match = re.search(r'<PP_NLCC[^>]*>([^<]+)</PP_NLCC>', content, re.DOTALL)
    if nlcc_match:
        nlcc_data = nlcc_match.group(1).strip()
        nlcc_values = []
        for line in nlcc_data.split('\n'):
            if line.strip():
                nlcc_values.extend([float(x) for x in line.split()])
        upf_data.rho_atc = np.array(nlcc_values)
    
    return upf_data


def write_upf_v2(upf: UPFData, filename: str):
    """Write UPF v2 format file"""
    with open(filename, 'w') as f:
        # Header
        f.write('<UPF version="2.0.1">\n')
        
        # PP_INFO
        f.write('  <PP_INFO>\n')
        f.write(upf.info)
        f.write('  </PP_INFO>\n')
        
        # PP_HEADER
        f.write('  <PP_HEADER\n')
        f.write(f'    generated="{upf.generated}"\n')
        f.write(f'    author="{upf.author}"\n')
        f.write(f'    date="{upf.date}"\n')
        f.write(f'    comment="{upf.comment}"\n')
        f.write(f'    element="{upf.psd}"\n')
        f.write(f'    pseudo_type="{upf.typ}"\n')
        f.write(f'    relativistic="{upf.rel}"\n')
        f.write(f'    is_ultrasoft="{str(upf.tvanp).lower()}"\n')
        f.write(f'    is_paw="{str(upf.tpawp).lower()}"\n')
        f.write(f'    is_coulomb="{str(upf.tcoulombp).lower()}"\n')
        f.write(f'    has_so="false"\n')
        f.write(f'    has_wfc="true"\n')
        f.write(f'    has_gipaw="false"\n')
        f.write(f'    paw_as_gipaw="false"\n')
        f.write(f'    core_correction="{str(upf.nlcc).lower()}"\n')
        f.write(f'    functional="{upf.dft}"\n')
        f.write(f'    z_valence="{upf.zp:.10f}"\n')
        f.write(f'    total_psenergy="{upf.etotps:.10f}"\n')
        f.write(f'    wfc_cutoff="{upf.ecutwfc:.10f}"\n')
        f.write(f'    rho_cutoff="{upf.ecutrho:.10f}"\n')
        f.write(f'    l_max="{upf.lmax}"\n')
        f.write(f'    l_max_rho="{upf.lmax_rho}"\n')
        f.write(f'    l_local="{upf.lloc}"\n')
        f.write(f'    mesh_size="{upf.mesh}"\n')
        f.write(f'    number_of_wfc="{upf.nwfc}"\n')
        f.write(f'    number_of_proj="{upf.nbeta}"/>\n')
        
        # PP_MESH
        f.write(f'  <PP_MESH dx="{upf.dx:.16e}" mesh="{upf.mesh}" xmin="{upf.xmin:.16e}" rmax="{upf.rmax:.16e}" zmesh="{upf.zmesh:.16e}">\n')
        
        # PP_R
        f.write('    <PP_R type="real" size="{}" columns="4">\n'.format(upf.mesh))
        for i in range(0, upf.mesh, 4):
            line = '      '
            for j in range(4):
                if i+j < upf.mesh:
                    line += f'{upf.r[i+j]:18.11e} '
            f.write(line.rstrip() + '\n')
        f.write('    </PP_R>\n')
        
        # PP_RAB
        f.write('    <PP_RAB type="real" size="{}" columns="4">\n'.format(upf.mesh))
        for i in range(0, upf.mesh, 4):
            line = '      '
            for j in range(4):
                if i+j < upf.mesh:
                    line += f'{upf.rab[i+j]:18.11e} '
            f.write(line.rstrip() + '\n')
        f.write('    </PP_RAB>\n')
        
        f.write('  </PP_MESH>\n')
        
        # PP_NLCC (if present)
        if upf.nlcc:
            f.write('  <PP_NLCC type="real" size="{}" columns="4">\n'.format(upf.mesh))
            for i in range(0, upf.mesh, 4):
                line = '    '
                for j in range(4):
                    if i+j < upf.mesh:
                        line += f'{upf.rho_atc[i+j]:18.11e} '
                f.write(line.rstrip() + '\n')
            f.write('  </PP_NLCC>\n')
        
        # PP_LOCAL
        f.write('  <PP_LOCAL type="real" size="{}" columns="4">\n'.format(upf.mesh))
        for i in range(0, upf.mesh, 4):
            line = '    '
            for j in range(4):
                if i+j < upf.mesh:
                    line += f'{upf.vloc[i+j]:18.11e} '
            f.write(line.rstrip() + '\n')
        f.write('  </PP_LOCAL>\n')
        
        # PP_NONLOCAL (if present)
        if upf.nbeta > 0:
            f.write('  <PP_NONLOCAL>\n')
            
            for i in range(upf.nbeta):
                f.write(f'    <PP_BETA.{i+1} type="real" size="{upf.kbeta[i]}" columns="4" ')
                f.write(f'angular_momentum="{upf.lll[i]}" cutoff_radius_index="{upf.kbeta[i]}" ')
                # f.write(f'cutoff_radius="{upf.rcut[i]:.10f}" ultrasoft_cutoff_radius="{upf.rcutus[i]:.10f}"')
                f.write('>\n')
                
                for ir in range(0, upf.kbeta[i], 4):
                    line = '      '
                    for j in range(4):
                        if ir+j < upf.kbeta[i]:
                            line += f'{upf.beta[i][ir+j]:18.11e} '
                    f.write(line.rstrip() + '\n')
                
                f.write(f'    </PP_BETA.{i+1}>\n')
            
            f.write('    <PP_DIJ type="real" size="{}" columns="4">\n'.format(upf.nbeta * upf.nbeta))
            for i in range(upf.nbeta):
                for j in range(0, upf.nbeta, 4):
                    line = '      '
                    for k in range(4):
                        if j+k < upf.nbeta:
                            line += f'{upf.dion[i, j+k]:18.11e} '
                    f.write(line.rstrip() + '\n')
            f.write('    </PP_DIJ>\n')
            
            f.write('  </PP_NONLOCAL>\n')
        
        # PP_PSWFC
        if upf.nwfc > 0:
            f.write('  <PP_PSWFC>\n')
            
            for i in range(upf.nwfc):
                f.write(f'    <PP_CHI.{i+1} type="real" size="{upf.mesh}" columns="4" ')
                # f.write(f'label="{upf.els[i]}"')
                f.write(f'l="{upf.lchi[i]}" occupation="{upf.oc[i]:.10f}" ')
                # f.write(f'n="{upf.nchi[i]}"')
                f.write(f' pseudo_energy="{upf.epseu[i]:.10f}"')
                # f.write(f' cutoff_radius="{upf.rcut_chi[i]:.10f}" ultrasoft_cutoff_radius="{upf.rcutus_chi[i]:.10f}"')
                f.write('>\n')
                
                for ir in range(0, upf.mesh, 4):
                    line = '      '
                    for j in range(4):
                        if ir+j < upf.mesh:
                            line += f'{upf.chi[i][ir+j]:18.11e} '
                    f.write(line.rstrip() + '\n')
                
                f.write(f'    </PP_CHI.{i+1}>\n')
            
            f.write('  </PP_PSWFC>\n')
        
        # PP_RHOATOM
        f.write('  <PP_RHOATOM type="real" size="{}" columns="4">\n'.format(upf.mesh))
        for i in range(0, upf.mesh, 4):
            line = '    '
            for j in range(4):
                if i+j < upf.mesh:
                    line += f'{upf.rho_at[i+j]:18.11e} '
            f.write(line.rstrip() + '\n')
        f.write('  </PP_RHOATOM>\n')
        
        f.write('</UPF>\n')

def main():
    """Main program"""
    if len(sys.argv) != 2:
        print("Usage: python reproduce_grid_cpmd2upf.py <cpmd_file>")
        sys.exit(1)
    
    filein = sys.argv[1]
    
    if not os.path.exists(filein):
        print(f"Error: file {filein} not found")
        sys.exit(2)
    
    print(f"Reading CPMD file: {filein}")
    
    try:
        # Read CPMD file
        cpmd_data = read_cpmd(filein)
        
        # Convert to UPF format
        upf_data = convert_cpmd(cpmd_data)
        
        # Write UPF file
        fileout = filein + '.UPF'
        print(f"Output PP file in UPF format: {fileout}")
        write_upf_v2(fileout, upf_data)
        
        print("Pseudopotential successfully written")
        print("Please review the content of the PP_INFO fields")
        print("*** Please TEST BEFORE USING !!! ***")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(3)

import numpy as np


def compare_upf_objects(upf1, upf2, rtol=1e-6, atol=1e-6):
    dict1, dict2 = upf1.__dict__, upf2.__dict__
    
    if set(dict1.keys()) != set(dict2.keys()):
        return False
    
    for key in dict1.keys():
        val1, val2 = dict1[key], dict2[key]
        if isinstance(val1, np.ndarray) or isinstance(val2, np.ndarray) \
            or isinstance(val1, list) or isinstance(val2, list):
            if not np.allclose(np.asarray(val1), np.asarray(val2), rtol=rtol, atol=atol):
                return False
        elif(val1!=val2):
            return False
    return True
        
if __name__ == '__main__':
    """for read-and-write (RW) test, usage: python upf_data.py <upf_file_in> <upf_file_out>"""
    import sys, os
    rfile = sys.argv[1]
    wfile = sys.argv[2]
    print("read file:", rfile)
    upf_ref = read_upf_file(rfile)
    # print("ref=", upf_ref.__dict__)
    write_upf_v2(upf_ref, wfile)
    assert(os.path.exists(wfile))
    upf_wr = read_upf_file(wfile)

    try:
        objects_equal = (upf_ref.__dict__ == upf_wr.__dict__)
        if hasattr(objects_equal, 'all'):
            objects_equal = objects_equal.all()
    except ValueError:
        # Handle array comparison manually
        objects_equal = compare_upf_objects(upf_ref, upf_wr)
    if objects_equal:
        print('RW test passed')
    else:
        print(f'RW test faild, ref=', upf_ref.__dict__)
        print(f'result=', upf_wr.__dict__)
    os.remove(wfile)