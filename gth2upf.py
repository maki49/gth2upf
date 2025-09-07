import json
import sys
import os
import re
import elec_config
def find_pseudo_content(gth_path, element, xc, valence):
    print("gth_path: ", gth_path)
    search_str = f"GTH-{xc.upper()}-q{valence}"
    pattern = fr'^{element}.*{search_str}.*'
    with open(gth_path, 'r') as f:
        lines = f.readlines()
    content = []
    found = False
    for line in lines:
        if re.match(pattern, line):
            found = True
        if found:
            if line.strip() == "#":
                break
            content.append(line)
    if not found:
        raise ValueError(f"Could not find pseudo potential for {element} with xc={xc} and valence {valence}")
    return ''.join(content[1:]) #remove the matched header line
  
def gen_cp2k_input(data):
    gth_path=os.path.join(data['cp2k_path'], 'data', 'GTH_POTENTIALS')
    element = data['element']
    xc= data.setdefault('xc', 'PBE')
    valence= data['valence']
    pseudo_content = find_pseudo_content(gth_path, element, xc, valence)
    prefix= data.setdefault('prefix', f"{element}-GTH-{xc}-q{valence}")
    core, valence = elec_config.get_core_valence(element, valence)
    
    quadrature = data.get('quadrature', 'GC_LOG')
    cp2k_quadrature = ['GC_LOG', 'GC_SIMPLE', 'GC_TRANSFORMED']
    if quadrature.upper() not in cp2k_quadrature:
        quadrature = 'GC_LOG'
        print(f"CP2K only supports QUADRATURE={cp2k_quadrature}. Using GC_LOG to generate the CP2K input file.")
    if quadrature.upper() == 'CPMD2UPF_DEFAULT':
        print("Then the grid of the CP2K-generated file will be replaced by the default grid used by cpmd2upf().")
        
    grid_points = data.setdefault('grid_points', 400)
    str_core =  f"CORE [{core[0]}]"
    if len(core) > 1:
        for item in core[1:]:
            str_core += f" {item}"
      
    str_valence ="CORE "+' '.join(valence)
    
    cp2k_inp_content=f"""
&GLOBAL
  PROJECT {element}
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT {element}
  ELECTRON_CONFIGURATION {str_valence}
  {str_core}
  &METHOD
     METHOD_TYPE  KOHN-SHAM
     &XC
       &XC_FUNCTIONAL {data['xc']}
       &END XC_FUNCTIONAL
     &END XC
  &END METHOD
  &PRINT
    &ANALYZE_BASIS
         OVERLAP_CONDITION_NUMBER T
         COMPLETENESS T
    &END ANALYZE_BASIS
    &UPF_FILE
       FILENAME ./{prefix}
    &END
  &END
  &PP_BASIS
     BASIS_TYPE GEOMETRICAL_GTO
     QUADRATURE {quadrature.upper()}
     GRID_POINTS {grid_points}
  &END PP_BASIS
  &POTENTIAL
    PSEUDO_TYPE GTH
    &GTH_POTENTIAL
    {pseudo_content}
    &END
  &END POTENTIAL
&END ATOM
"""
    with open(f"{element}.inp", 'w') as f:
        f.write(cp2k_inp_content)
        return f"{element}.inp"

def run_cp2k(cp2k_dir, inp, out="cp2k.out"):
    os.system(f"{os.path.join(cp2k_dir, 'exe', 'local', 'cp2k.popt')} -o {out} {inp}")

def postprocess(data):
    outfile=data['prefix']+"-1.upf"
    os.system(f"sed -i 's/functional=\"DFT\"/functional=\"{data['xc'].upper()}\"/' {outfile}")

def gen_single_upf(data):
    inp = gen_cp2k_input(data)
    run_cp2k(data['cp2k_path'], inp)
    postprocess(data)
    
def main(file_path):
    with open(file_path, 'r') as file:
        params = json.load(file)
    gen_single_upf(params)

    from upf_data import read_upf_file, write_upf_v2
    # postprocess
    upf_gen = read_upf_file(params['prefix']+"-1.upf")
    upf = upf_gen.correct_lmax().correct_zmesh()
    if params['quadrature'].upper() == 'CPMD2UPF_DEFAULT':
        upf = upf.replace_grid()
    write_upf_v2(upf, params['prefix']+"-1.upf")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gth2upf.py <params.json>")
        sys.exit(1)
    main(sys.argv[1])