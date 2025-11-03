# Author: Jincheng Yu <pimetamon@gmail.com>
import numpy as np
import numpy
from math import log, exp
import sys
from gth_tools import parse_gth_pp, V_loc, p_il

def eval_gth2upf(element, gthdata, zmesh=None, 
            xmin=-7.00, rmax_init=100.0, dx=0.0125, upf_file=None):
    """
    Evaluate GTH PPs on grids and generate UPF file for QE.

    Args:
    element (str): name of the element
    gthdata (str): GTH data, e.g.
        '''
            C GTH-PBE-q4 GTH-PBE
        2    2
         0.33847124    2    -8.80367398     1.33921085
        2
         0.30257575    1     9.62248665
         0.29150694    0
        '''
    zmesh (int): atomic number
    xmin (double): parameter used to generate grids, 
        ri = 1/z * exp(xmin + i * dx)
    rmax_init (double): parameter used to generate grids
    dx (double): parameter used to generate grids
    """
    gth = parse_gth_pp(gthdata)
    if upf_file == None:
        upf_file = f'%s-%s-q%s.UPF' % (element, gth.xc_label, gth.q)

    # ==> Step 1: generate grids <======================================
    if zmesh == None:
        from pyscf.data import elements
        zmesh = elements.charge(element)
    grids = GridCPMD2UPF(zmesh=zmesh, xmin=xmin, rmax_init=rmax_init, dx=dx)

    # ==> Step 2: evaluate local part of the PP <=======================
    vloc = V_loc(grids.r, Z_ion=gth.local.Z_ion, r_loc=gth.local.r_loc,
                 C1=gth.local.C[0], C2=gth.local.C[1], 
                 C3=gth.local.C[2], C4=gth.local.C[3]) * 2 # A.U. to Ry

    # ==> Step 3: evaluate projectors <=================================
    proj_list = []
    for ch in gth.channels:
        for i in range(ch.m):
            proj_list.append(
                (p_il(grids.r, ch.l, i+1, ch.r_l) * grids.r * 2, ch.l)
                )

    # ==> Step 4: write UPF file <======================================
    # ====> 4.1: print info <===========================================
    with open(upf_file, 'w+') as f:
        f.write('<UPF version="2.0.1">\n')
        f.write('<PP_INFO>\n')
        f.write('Converted from CP2K GTH format\n')
        f.write('<PP_INPUTFILE>\n')
        f.write(gthdata)
        f.write('\n')
        f.write('</PP_INPUTFILE>\n')
        f.write('</PP_INFO>\n')
        f.write('<PP_HEADER generated="Generated in analytical, separable form"\n')
        f.write('author="Goedecker/Hartwigsen/Hutter/Teter"\n')
        f.write('date="Phys.Rev.B58, 3641 (1998); B54, 1703 (1996)"\n')
        f.write('comment="PP_CHI.X and PP_RHOATOM sections deleted"\n')
        f.write('element="%s"\n' % element)
        f.write('pseudo_type="NC"\n')
        f.write('relativistic="no"\n')
        f.write('is_ultrasoft="F"\n')
        f.write('is_paw="F"\n')
        f.write('is_coulomb="F"\n')
        f.write('has_so="F"\n')
        f.write('has_wfc="F"\n')
        f.write('has_gipaw="F"\n')
        f.write('paw_as_gipaw="F"\n')
        f.write('core_correction="F"\n')
        f.write('functional="%s"\n' % gth.xc_label)
        f.write('z_valence="%f"\n' % float(gth.q))
        f.write('total_psenergy="0.000000000000000E+000"\n')
        f.write('wfc_cutoff="0.000000000000000E+000"\n')
        f.write('rho_cutoff="0.000000000000000E+000"\n')
        f.write('l_max="%d"\n' % len(gth.channels))
        f.write('l_max_rho="0"\n')
        f.write('l_local="-3"\n')
        f.write('mesh_size="%d"\n' % grids.r.shape[0])
        f.write('number_of_wfc="0"\n')
        f.write('number_of_proj="%d"/>\n' % len(proj_list))

        f.write('<PP_MESH\n') 
        f.write(f'dx="{dx:.16e}"\n') 
        f.write(f'mesh="{grids.r.shape[0]}"\n') 
        f.write(f'xmin="{xmin:.16e}"\n') 
        f.write(f'rmax="{grids.rmax:.16e}"\n')
        f.write(f'zmesh="{zmesh:.16e}">\n')
    # ====> 4.2: print grids <==========================================
        f.write(f'<PP_R type="real" size="{grids.r.shape[0]}" columns="4">\n')
        data = grids.r
        for i, val in enumerate(data, start=1):
            f.write(f"{val:20.16e} ")
            if i % 4 == 0:
                f.write("\n")
        f.write('\n</PP_R>\n')
    # ====> 4.2: print weights <========================================
        f.write(f'<PP_RAB type="real" size="{grids.r.shape[0]}" columns="4">\n')
        data = grids.rab
        for i, val in enumerate(data, start=1):
            f.write(f"{val:20.16e} ")
            if i % 4 == 0:
                f.write("\n")
        f.write('\n</PP_RAB>\n')
        f.write('</PP_MESH>\n')
    # ====> 4.2: print V_loc <==========================================
        f.write(f'<PP_LOCAL type="real" size="{grids.r.shape[0]}" columns="4">\n')
        data = vloc
        for i, val in enumerate(data, start=1):
            f.write(f"{val:20.16e} ")
            if i % 4 == 0:
                f.write("\n")
        f.write('\n</PP_LOCAL>\n')
    # ====> 4.3: print projectors <=====================================
        f.write('<PP_NONLOCAL>\n')
        for idxProj, proj in enumerate(proj_list, start=1):
            f.write(f'<PP_BETA.%d ' % idxProj)
            f.write(f'type="real" size="{grids.r.shape[0]}" columns="4" ')
            f.write(f'angular_momentum="{proj[1]}" ')
            f.write(f'cutoff_radius_index="{grids.r.shape[0]}">\n')
            data = proj[0]
            for i, val in enumerate(data, start=1):
                f.write(f"{val:20.16e} ")
                if i % 4 == 0:
                    f.write("\n")
            f.write(f'\n</PP_BETA.%d>\n' % idxProj)
    # ====> 4.4 print hij <=============================================
        f.write(f'<PP_DIJ type="real" size="%d" columns="4">\n' 
              % len(proj_list)**2)
        hij, _ = combine_channels(gth.channels)
        data = hij.ravel() * 0.5
        for i, val in enumerate(data, start=1):
            f.write(f"{val:20.16e} ")
            if i % 4 == 0:
                f.write("\n")
        if len(proj_list) % 2 == 0:
            f.write(f'</PP_DIJ>\n')
        else:
            f.write(f'\n</PP_DIJ>\n')
        f.write(f'</PP_NONLOCAL>\n')

        f.write('<PP_PSWFC>\n')
        f.write('</PP_PSWFC>\n')
        f.write('</UPF>\n')

def _packed_to_symm(hij_packed, m):
    """
    Convert a packed upper-triangular list of h_ij values into
    a full symmetric m x m matrix.

    Packing order is:
        (0,0), (0,1), ... (0,m-1),
        (1,1), (1,2), ...,
        ...
        (m-1,m-1)
    """
    need = m * (m + 1) // 2
    if len(hij_packed) != need:
        raise ValueError(f"Expected {need} coeffs for m={m}, got {len(hij_packed)}")

    M = np.zeros((m, m), dtype=float)
    k = 0
    for i in range(m):
        for j in range(i, m):
            M[i, j] = hij_packed[k]
            M[j, i] = hij_packed[k]
            k += 1
    return M

def combine_channels(channels, sort_by_l=False):
    """
    Build a single block-diagonal coupling matrix from all channels.

    Parameters
    ----------
    channels : list of Channel
        Each Channel has l, r_l, m, h_ij (packed).
    sort_by_l : bool
        If True, channels are ordered by l value (then by input order).
        If False, keep original order.

    Returns
    -------
    H_global : 2D numpy array
        Block-diagonal matrix with one block per channel.
    index_map : list of tuples
        index_map[k] = (l, ch_index, local_i)
        - l          : angular momentum for this channel
        - ch_index   : index of this channel (in sorted order)
        - local_i    : which projector inside that channel (0..m-1)
    """
    if sort_by_l:
        # stable sort by (l, original_index)
        channels = sorted(list(enumerate(channels)), key=lambda p: (p[1].l, p[0]))
    else:
        channels = list(enumerate(channels))

    blocks = []
    index_map = []
    for ch_sorted_idx, ch in channels:
        H_l = _packed_to_symm(ch.h_ij, ch.m)
        blocks.append(H_l)
        for i in range(ch.m):
            index_map.append((ch.l, ch_sorted_idx, i))

    # Assemble block-diagonal matrix
    total = sum(b.shape[0] for b in blocks)
    H_global = np.zeros((total, total), dtype=float)
    offset = 0
    for B in blocks:
        m = B.shape[0]
        H_global[offset:offset+m, offset:offset+m] = B
        offset += m

    return H_global, index_map


def combine_channels_soc(channels, soc_channels):
    """
    Build a single block-diagonal coupling matrix from all channels,
    including SOC channels.

    Parameters
    ----------
    channels : list of Channel
        Each Channel has l, r_l, m, h_ij (packed).
    soc_channels : list of Channel
        Each Channel has l, r_l, m, k_ij (packed).

    Returns
    -------
    H_global : 2D numpy array
        Block-diagonal matrix with one block per channel.
    """
    # Combine non-SOC and SOC channels
    channels = list(channels)
    soc_channels = list(soc_channels)

    n0 = len([ch for ch in channels if ch.l==0])   # #number of l=0 channels
    blocks = []
    for idx, ch in enumerate(channels):
        if ch.l == 0:
            D_l = _packed_to_symm(ch.h_ij, ch.m)
        else:
            H_l = _packed_to_symm(ch.h_ij, ch.m)
            idx_soc = idx-n0
            K_l = _packed_to_symm(soc_channels[idx_soc].h_ij, soc_channels[idx_soc].m)
            D_l = np.zeros((2*ch.m, 2*ch.m), dtype=float)
            vnlp = H_l + 0.5 * K_l * ch.l
            vnlm = H_l - 0.5 * K_l * (ch.l + 1)
            for i in range(ch.m):
                ii = 2*i
                for j in range(ch.m):
                    jj = 2*j
                    D_l[ii, jj]     = vnlm[i, j]
                    D_l[ii+1, jj+1] = vnlp[i, j]
        blocks.append(D_l)

    # Assemble block-diagonal matrix
    total = sum(b.shape[0] for b in blocks)
    H_global = np.zeros((total, total), dtype=float)
    offset = 0
    for B in blocks:
        m = B.shape[0]
        H_global[offset:offset+m, offset:offset+m] = B
        offset += m

    return H_global