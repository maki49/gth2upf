import os
import re
import json
import sys
def grep_cp2k_basis_molopt(data):
    element = data['element']
    basis_type = data.setdefault('basis_type', 'DZVP')
    basis_file = os.path.join(data['cp2k_path'], 'data', 'BASIS_MOLOPT')
    is_sr=data.setdefault('short_range', False)
    short_range=""
    if is_sr:
        short_range="SR-"
    pattern=fr'^\s*{element}\s*{basis_type}-MOLOPT-{short_range}GTH.*'
    content=[]
    found=False
    with open(basis_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if re.match(pattern, line):
            found=True
            content.append(line)
            continue
        if found:
            if not re.match(r'^\s*\d', line):
                break
            print("line: ", line)
            content.append(line)

    writefile=f"{element}-{basis_type}-MOLOPOT-{short_range}GTH"
    with open(writefile, 'w') as f:
        f.writelines(content)
    return writefile, short_range

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python genorb.py cp2k_basis.json")
        exit(1)
    file_path = sys.argv[1]
    with open(file_path, 'r') as file:
        data = json.load(file)
    writefile, short_range = grep_cp2k_basis_molopt(data)
    orb_submodule_pth=f"{os.path.dirname(os.path.abspath(__file__))}/deps/gaussian_orbital_for_ABACUS/cp2k2abacus_gaussian.py"
    os.system(f"python {orb_submodule_pth} {writefile} {writefile} {data['basis_rcut']} MOLOPT-{short_range}GTH")