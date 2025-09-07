# UPF RW tests
python ../upf_data.py ref_upf/C.pbe-hgh.UPF C.res.UPF
python ../upf_data.py ref_upf/C-GTH-PBE-1.upf C.res.UPF
# GTH RW tests
python ../gth_data.py ref_gth/C_gth C.res.gth
python ../gth_data.py ref_gth/Au_gth Au.res.gth

# GTH to UPF tests
python ../gth2upf_custom.py ref_gth/C_gth ref_upf/C.pbe-hgh.UPF
