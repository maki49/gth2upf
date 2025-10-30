import sys
sys.path.append("../..")
import eval_gth2upf

gthdata = """
    C GTH-PBE-q4 GTH-PBE
    2    2
     0.33847124    2    -8.80367398     1.33921085
    2
     0.30257575    1     9.62248665
     0.29150694    0
    """

eval_gth2upf.eval_gth2upf('C', gthdata)

