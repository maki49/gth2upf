import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import re

def parse_key(content, key):
    re_dict = {
        'r': r'<PP_R[^>]*>([\s\S]*?)</PP_R>',
        'rab': r'<PP_RAB[^>]*>([\s\S]*?)</PP_RAB>',
        'local': r'<PP_LOCAL[^>]*>([\s\S]*?)</PP_LOCAL>',
        'rho': r'<PP_RHOATOM[^>]*>([\s\S]*?)</PP_RHOATOM>',
        'chi1': r'<PP_CHI\.1[^>]*>([\s\S]*?)</PP_CHI\.1>',
        'chi2': r'<PP_CHI\.2[^>]*>([\s\S]*?)</PP_CHI\.2>',
    }
    match = re.search(re_dict[key], content)
    if match:
        str = match.group(1)
        return np.array([float(x) for x in str.split()])    #string to array
    print(f'key {key} not found, content=', content)
    return None

def parse_upf_data(file_path, keys):
    """parse the part of key in upf file"""
    with open(file_path, 'r') as f:
        content = f.read()
    return {k: parse_key(content, k) for k in keys}

def cal_integral(r, rab, local):
    assert len(r) == len(rab) == len(local)
    return np.trapz(rab * local, r)

def compare_r(file_list, xmax=None, ymax=None):
    plt.figure(figsize=(10, 6))
    for f in file_list:
        r = parse_upf_data(f, ['r'])['r']
        plt.plot(r, '-*', label=f)
    plt.xlabel('mesh index', fontsize=14)
    plt.ylabel('PP_R[i] (Bohr)', fontsize=14)
    plt.title(f'R', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    if xmax is not None and ymax is not None:
        plt.xlim(0, xmax)    #a.u.
        plt.ylim(0, ymax)
    plt.tight_layout()
    plt.savefig('upf_comparison_r.png', 
                dpi=300, bbox_inches='tight')

def compare_rab(file_list, xmax=None, ymax=None):
    plt.figure(figsize=(10, 6))
    for f in file_list:
        r = parse_upf_data(f, ['r'])['r']
        rab = parse_upf_data(f, ['rab'])['rab']
        plt.plot(r, rab, '-*', label=f)
    plt.xlabel('PP_R (Bohr)', fontsize=14)
    plt.ylabel('PP_RAB (a.u.)', fontsize=14)
    plt.title(f'RAB', fontsize=16)
    if xmax is not None and ymax is not None:
        plt.xlim(0, xmax)    #a.u.
        plt.ylim(0, ymax)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'upf_comparison_rab.png', 
                dpi=300, bbox_inches='tight')
        
def compare_local(file_list):
    plt.figure(figsize=(10, 6))
    for f in  file_list:
        data = parse_upf_data(f, ['r', 'rab', 'local'])
        r, rab, local = data['r'], data['rab'], data['local']
        print(f"{f} integral:", cal_integral(r, rab, local))
        plt.plot(r, local, label=f)
    plt.xlabel('PP_R (Bohr)', fontsize=14)
    plt.ylabel('PP_LOCAL (a.u.)', fontsize=14)
    plt.title('LOCAL', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('upf_comparison_local.png', 
                dpi=300, bbox_inches='tight')
    
def compare_rho(file_list):
    plt.figure(figsize=(10, 6))
    for f in  file_list:
        data = parse_upf_data(f, ['r', 'rho'])
        r, rho = data['r'], data['rho']
        plt.plot(r, rho, label=f)
    plt.xlabel('PP_R (Bohr)', fontsize=14)
    plt.ylabel('PP_RHOATOM (a.u.)', fontsize=14)
    plt.title('RHOATOM', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('upf_comparison_rho.png', 
                dpi=300, bbox_inches='tight')
    
def compare_chi(file_list, id=1):
    plt.figure(figsize=(10, 6))
    for f in  file_list:
        data = parse_upf_data(f, ['r', f'chi{id}'])
        r, chi = data['r'], data[f'chi{id}']
        plt.plot(r, chi, label=f"{f}")
    plt.xlabel('PP_R (Bohr)', fontsize=14)
    plt.ylabel(f'PP_CHI.{id} (a.u.)', fontsize=14)
    plt.title(f'CHI.{id}', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'upf_comparison_chi{id}.png', 
                dpi=300, bbox_inches='tight')
    
    
def compare_dr(file_list, xmax = None, ymax = None):  # r[i]-r[i-1]
    plt.figure(figsize=(10, 6))
    for f in file_list:
        r = parse_upf_data(f, ['r'])['r']
        plt.plot((r[1:]-r[:-1]), label=f)
    plt.xlabel('mesh index', fontsize=14)
    plt.ylabel('PP_R[i]-PP_R[i-1] (Bohr)', fontsize=14)
    plt.title(f'dR', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    if xmax is not None and ymax is not None:
        plt.xlim(0, xmax)    #a.u.
        plt.ylim(0, ymax)
    plt.tight_layout()
    plt.savefig('upf_comparison_dr.png', 
                dpi=300, bbox_inches='tight')

def compare_logdr(file_list):  # r[i]-r[i-1]
    plt.figure(figsize=(10, 6))
    for f in file_list:
        r = parse_upf_data(f, ['r'])['r']
        plt.semilogy((r[1:]-r[:-1]), label=f)
    plt.xlabel('mesh index', fontsize=14)
    plt.ylabel('PP_R[i]-PP_R[i-1] (Bohr)', fontsize=14)
    plt.title(f'dR\n{f}', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('upf_comparison_logdr.png', 
                dpi=300, bbox_inches='tight')
    
def fit_sinh(file):
    import sinh_fitting
    r = parse_upf_data(file, ['r'])['r']
    dr=r[1:]-r[:-1]
    sinh_fitting.fit_and_plot(r[:-1], dr, file_name=file.split('.')[0])

if __name__ == '__main__':
    
    task = 3
    
    match task:
        case 0:
            # task 0: math fitting
            file_list=['C-GTH-PBE-40au-1.upf',
            'C-GTH-PBE-100au-1.upf',
            'C.pbe-hgh.UPF']
            compare_logdr(file_list)
            fit_sinh(file_list[0])
            fit_sinh(file_list[1])
        case 1:
            # task 1: change the number of grid points
            file_list=['C-GTH-PBE-40au-1.upf',
                        'C-GTH-PBE-100au-1.upf',
                        'C.pbe-hgh.UPF']
            compare_r(file_list)
            compare_rab(file_list)
            compare_local(file_list)
        case 2:
            # task 2: change the quadrature type
            file_list=['C-GTH-PBE-GC_LOG-1.upf',
                        'C-GTH-PBE-GC_SIMPLE-1.upf',
                        'C-GTH-PBE-GC_TRANSFORMED-1.upf',
                        'C.pbe-hgh.UPF']
            compare_dr(file_list, xmax=1000, ymax=5)
            compare_r(file_list, xmax=1000, ymax = 50)
            compare_rab(file_list, xmax=100, ymax=10)
        case 3:
            # task 3: reproduce C.pbe-hgh.UPF
            #    compare r and rab of the result and test/ref_upf/C.pbe-hgh.UPF 
            import sys, os
            sys.path.append('..')
            from gth2upf import main
            main('../examples/C.json')  
            res_file = 'C-GTH-PBE-1.upf'
            file_list=[res_file, '../test/ref_upf/C.pbe-hgh.UPF']
            
            compare_r(file_list)
            compare_rab(file_list)
            
            file_list= [res_file, '../test/ref_upf/C-GTH-PBE-1.upf']
            # file_list= ['../test/ref_upf/C.pbe-hgh.UPF', res_file, '../test/ref_upf/C-GTH-PBE-1.upf']
            compare_local(file_list)
            compare_rho(file_list)
            compare_chi(file_list, 1)
            compare_chi(file_list, 2)
            
            os.remove(res_file)
            os.remove('C.inp')
            os.remove('cp2k.out')