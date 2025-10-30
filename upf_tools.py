# Author: Jincheng Yu <pimetamon@gmail.com>
import numpy as np

def parser(finp, node_names):
    """
    Example: 
    if __name__ == '__main__':
        upf = 'C.pbe-hgh.UPF'
        node_names = [
        'PP_MESH/PP_R',
        'PP_MESH/PP_RAB',
        'PP_LOCAL',
        'PP_NONLOCAL/PP_BETA.1',
        'PP_PSWFC/PP_CHI.1',
        'PP_PSWFC/PP_CHI.2',
        'PP_RHOATOM'
        ]

        data_all = parser(upf, node_names)
        pp_loc = data_all['PP_LOCAL']
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(finp)
    root = tree.getroot()

    def read_node(root, node_name):
        try:
            node = root.find(node_name)
        except:
            return None
        return np.asarray(list(map(float, [x for line in node.text.split('\n')
                                           for x in line.split()])))

    res = {node_name:read_node(root, node_name) for node_name in node_names}

    return res

