import re
electron_configurations = {
    'H': ['1s1'],
    'He': ['1s2'],
    'Li': ['He', '2s1'],
    'Be': ['He', '2s2'],
    'B': ['He', '2s2', '2p1'],
    'C': ['He', '2s2', '2p2'],
    'N': ['He', '2s2', '2p3'],
    'O': ['He', '2s2', '2p4'],
    'F': ['He', '2s2', '2p5'],
    'Ne': ['He', '2s2', '2p6'],
    'Na': ['Ne', '3s1'],
    'Mg': ['Ne', '3s2'],
    'Al': ['Ne', '3s2', '3p1'],
    'Si': ['Ne', '3s2', '3p2'],
    'P': ['Ne', '3s2', '3p3'],
    'S': ['Ne', '3s2', '3p4'],
    'Cl': ['Ne', '3s2', '3p5'],
    'Ar': ['Ne', '3s2', '3p6'],
    'K': ['Ar', '4s1'],
    'Ca': ['Ar', '4s2'],
    'Sc': ['Ar', '4s2', '3d1'],
    'Ti': ['Ar', '4s2', '3d2'],
    'V': ['Ar', '4s2', '3d3'],
    'Cr': ['Ar', '4s1', '3d5'],
    'Mn': ['Ar', '4s2', '3d5'],
    'Fe': ['Ar', '4s2', '3d6'],
    'Co': ['Ar', '4s2', '3d7'],
    'Ni': ['Ar', '4s2', '3d8'],
    # 'Cu': ['Ar', '4s1', '3d10'],  # Free Cu atom, following the aufbau principle
    'Cu': ['Ar', '4s2', '3d9'],  # For ppotential generation, use 4s2 3d9
    'Zn': ['Ar', '4s2', '3d10'],
    'Ga': ['Ar', '4s2', '3d10', '4p1'],
    'Ge': ['Ar', '4s2', '3d10', '4p2'],
    'As': ['Ar', '4s2', '3d10', '4p3'],
    'Se': ['Ar', '4s2', '3d10', '4p4'],
    'Br': ['Ar', '4s2', '3d10', '4p5'],
    'Kr': ['Ar', '4s2', '3d10', '4p6'],
    'Rb': ['Kr', '5s1'],
    'Sr': ['Kr', '5s2'],
    'Y': ['Kr', '5s2', '4d1'],
    'Zr': ['Kr', '5s2', '4d2'],
    'Nb': ['Kr', '5s1', '4d4'],
    'Mo': ['Kr', '5s1', '4d5'],
    'Tc': ['Kr', '5s1', '4d5'],
    'Ru': ['Kr', '5s1', '4d7'],
    'Rh': ['Kr', '5s1', '4d8'],
    'Pd': ['Kr', '5s0', '4d10'],
    # 'Ag': ['Kr', '5s1', '4d10'],
    'Ag': ['Kr', '5s2', '4d9'],
    'Cd': ['Kr', '5s2', '4d10'],
    'In': ['Kr', '5s2', '4d10', '5p1'],
    'Sn': ['Kr', '5s2', '4d10', '5p2'],
    'Sb': ['Kr', '5s2', '4d10', '5p3'],
    'Te': ['Kr', '5s2', '4d10', '5p4'],
    'I': ['Kr', '5s2', '4d10', '5p5'],
    'Xe': ['Kr', '5s2', '4d10', '5p6'],
    'Cs': ['Xe', '6s1'],
    'Ba': ['Xe', '6s2'],
    'La': ['Xe', '6s2', '5d1'],
    'Ce': ['Xe', '6s2', '4f1', '5d1'],
    'Pr': ['Xe', '6s2', '4f3', '5d1'],
    'Nd': ['Xe', '6s2', '4f4', '5d1'],
    'Pm': ['Xe', '6s2', '4f5', '5d1'],
    'Sm': ['Xe', '6s2', '4f6', '5d1'],
    'Eu': ['Xe', '6s2', '4f7'],
    'Gd': ['Xe', '6s2', '4f7', '5d1'],
    'Tb': ['Xe', '6s2', '4f9', '5d1'],
    'Dy': ['Xe', '6s2', '4f10', '5d1'],
    'Ho': ['Xe', '6s2', '4f11', '5d1'],
    'Er': ['Xe', '6s2', '4f12', '5d1'],
    'Tm': ['Xe', '6s2', '4f13', '5d1'],
    'Yb': ['Xe', '6s2', '4f14'],
    'Lu': ['Xe', '6s2', '4f14', '5d1'],
    'Hf': ['Xe', '6s2', '4f14', '5d2'],
    'Ta': ['Xe', '6s2', '4f14', '5d3'],
    'W': ['Xe', '6s2', '4f14', '5d4'],
    'Re': ['Xe', '6s2', '4f14', '5d5'],
    'Os': ['Xe', '6s2', '4f14', '5d6'],
    'Ir': ['Xe', '6s2', '4f14', '5d7'],
    'Pt': ['Xe', '6s1', '4f14', '5d9'],
    # 'Au': ['Xe', '6s1', '4f14', '5d10'],
    'Au': ['Xe', '6s2', '4f14', '5d9'],
    'Hg': ['Xe', '6s2', '4f14', '5d10'],
    'Tl': ['Xe', '6s2', '4f14', '5d10', '6p1'],
    'Pb': ['Xe', '6s2', '4f14', '5d10', '6p2'],
    'Bi': ['Xe', '6s2', '4f14', '5d10', '6p3'],
    'Po': ['Xe', '6s2', '4f14', '5d10', '6p4'],
    'At': ['Xe', '6s2', '4f14', '5d10', '6p5'],
    'Rn': ['Xe', '6s2', '4f14', '5d10', '6p6'],
    'Fr': ['Rn', '7s1'],
    'Ra': ['Rn', '7s2'],
    # Add more elements as needed
}

def count_valence_num(lst):
    valence_num = 0
    for item in lst:
        if re.match(r'\d+', item):
            valence_num += int(re.search(r'\d+$', item).group())  #`\d+$` match the end number, and .group() extract the matched number
    return valence_num

# 新增：按 n 升序、壳层类型 s<p<d<f 排序
_L_ORDER = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}
_SHELL_RE = re.compile(r'^(\d+)([spdfg])(\d+)$')

def _shell_key(s: str):
    m = _SHELL_RE.match(s)
    if not m:
        return (-1, -1)  # 非壳层字符串（如惰性气体符号）放到最前面
    n = int(m.group(1))
    l = _L_ORDER[m.group(2)]
    return (n, l)

def sort_valence(shells):
    # 只对形如 "4f14" 的条目排序
    return sorted([s for s in shells], key=_shell_key)

def expand_element(element, depth=1):
    """Expand core notation to full electron configuration."""
    config = electron_configurations[element]
    if depth == 0:
        return config
    core = config[0]
    val = config[1:]
    config_core = expand_element(core, depth - 1)
    expanded = config_core + val
    return expanded

def expand_core(config):
    return expand_element(config[0], depth=0) + config[1:]

def get_core_valence(element, nval_requested):
    config = sort_valence(electron_configurations.get(element))
    core = config[0]
    val = config[1:]
    nval = count_valence_num(val)
    if element == 'H' or element == 'He':
        return ["none"], electron_configurations[element]
    elif nval == nval_requested:   #treat the core as another element to analyze
        return [core], val
    else:
        config = sort_valence(expand_core(config))  #expand the core
        val = []
        for str_try in reversed(config):
            val.insert(0, str_try)
            config.pop()
            nval_try = count_valence_num(val)
            if(nval_try == nval_requested):
                return config, val
    raise ValueError(f"Could not find core and valence for {element} with valence {nval}")


if __name__ == "__main__":
    test=2
    match test:
        case 0: # expand_element
            print(expand_element('Au', depth=1))
        case 1: # sort_valence
            print(sort_valence(expand_element('Au', depth=1)))
        case 2: # get_core_valence
            cases = [
                ('Au', 19),
                ('W', 14),
                ('W', 28)
            ]
            for element, valence in cases:
                core, valence_shells = get_core_valence(element, valence)
                print(f"(Element: {element}, Valence: {valence}) => Core: {core}, Valence Shells: {valence_shells}")