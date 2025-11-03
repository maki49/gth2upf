import sys
sys.path.append('..')
from gth2upf import gen_single_upf
from upf_tools import parser

# input_origin =    {
#     "cp2k_path": "/mnt/work/hack/cp2k/exe/local/cp2k.psmp",
#     "gth_path": "/mnt/work/hack/cp2k/data/GTH_SOC_POTENTIALS",
#     "element": "C",
#     "xc": "PBE",
#     "valence": 4,
#     "prefix": "C-GTH-SOC-PBE",
# }
# input_eval =    {
#     "cp2k_path": "/mnt/work/hack/cp2k/exe/local/cp2k.psmp",
#     "gth_path": "/mnt/work/hack/cp2k/data/GTH_SOC_POTENTIALS",
#     "element": "C",
#     "xc": "PBE",
#     "valence": 4,
#     "prefix": "C-GTH-SOC-PBE-eval",
#     "grid_points": 500,
#     "quadrature": "cpmd2upf_default",
# }

input_origin =    {
    "cp2k_path": "/mnt/work/hack/cp2k/exe/local/cp2k.psmp",
    "gth_path": "/mnt/work/hack/cp2k/data/GTH_SOC_POTENTIALS",
    "element": "W",
    "xc": "PBE",
    "valence": 14,
    "prefix": "W-GTH-SOC-PBE",
}
input_eval =    {
    "cp2k_path": "/mnt/work/hack/cp2k/exe/local/cp2k.psmp",
    "gth_path": "/mnt/work/hack/cp2k/data/GTH_SOC_POTENTIALS",
    "element": "W",
    "xc": "PBE",
    "valence": 14,
    "prefix": "W-GTH-SOC-PBE-eval",
    "grid_points": 500,
    "quadrature": "cpmd2upf_default",
}

if __name__ == "__main__":
    upf_origin = gen_single_upf(input_origin)
    upf_eval = gen_single_upf(input_eval)
    
    nodes = [
        'PP_MESH/PP_R',
        'PP_MESH/PP_RAB',
        'PP_LOCAL',
        'PP_NONLOCAL/PP_BETA.3',
        'PP_NONLOCAL/PP_BETA.9',
    ]
    data_origin = parser(upf_origin, nodes)
    data_eval = parser(upf_eval, nodes)

    # plot rab, local, beta.3 vs. r
    r = data_origin['PP_MESH/PP_R']
    r_eval = data_eval['PP_MESH/PP_R']
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 8))
    plt.subplot(4,1,1)
    plt.plot(r, label='origin')
    plt.plot(r_eval, label='eval', linestyle='--')
    plt.ylabel('R')
    plt.legend()
    plt.subplot(4,1,2)
    plt.plot(r, data_origin['PP_MESH/PP_RAB'], label='origin')
    plt.plot(r_eval, data_eval['PP_MESH/PP_RAB'], label='eval', linestyle='--')
    plt.ylabel('RAB')
    plt.legend()
    plt.subplot(4,1,3)
    plt.plot(r, data_origin['PP_LOCAL'], label='origin')
    plt.plot(r_eval, data_eval['PP_LOCAL'], label='eval', linestyle='--')
    plt.ylabel('LOCAL')
    plt.subplot(4,1,4)
    # plt.plot(r, data_origin['PP_NONLOCAL/PP_BETA.3'], label='origin')
    # plt.plot(r_eval, data_eval['PP_NONLOCAL/PP_BETA.3'], label='eval', linestyle='--')
    # plt.ylabel('BETA.3')
    plt.plot(r, data_origin['PP_NONLOCAL/PP_BETA.9'], label='origin')
    plt.plot(r_eval, data_eval['PP_NONLOCAL/PP_BETA.9'], label='eval', linestyle='--')
    plt.ylabel('BETA.9')
    plt.legend()

    plt.savefig('compare_eval_vs_origin.png', dpi=300)