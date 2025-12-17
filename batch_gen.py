import json
import os

# Configuration
cp2k_path = "/mnt/work/hack/cp2k/exe/local/cp2k.psmp" 
gth_path = "/mnt/work/hack/cp2k/data/GTH_POTENTIALS"
xc = "PBE"
result_path = "results/"

# input:dict[element: str, z_valence: int]
input={
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "Mg": 10,
    "Al": 3,   # not having 13 
    "P": 5,
    "Ga": 3,
    "As": 5,
    "Sb": 5,   # not having 15
    "Cd": 12,   # not having 20
}

def batchgen_cp2k_default_grid(input):
    """Batch generating PP in upf format from GTH potentials with CP2K's default grid settings.
    Args:
        input (dict[element: str, z_valence: int]): Dictionary with element symbols as keys and their valence electron counts as values.
    """
    for element, z_valence in input.items():
        # 1. generate json file
        prefix=f"{element}-GTH-{xc}-cp2k-default-grid"
        dict={
            "cp2k_path": cp2k_path,
            "gth_path": gth_path,
            "element": element,
            "xc": xc,
            "valence": z_valence,
            "prefix": prefix,
        }
        json_file_name = f"{element}-default.json"
        with open(json_file_name, "w") as f:
            json.dump(dict, f, indent=4)
            
        # 2. run gth2upf
        os.system(f"python gth2upf.py {json_file_name}")
        
        # 3. copy results to result_path
        upf_raw_name = f"{prefix}-1.upf"
        upf_path = os.path.join(result_path, f"{element}-GTH-{xc}-q{z_valence}-cp2k-default-grid.upf")
        os.system(f"mv {upf_raw_name} {upf_path}")
        print(f"Generated UPF for {element}: {upf_path}")
        
        
if __name__ == "__main__":
    batchgen_cp2k_default_grid(input)