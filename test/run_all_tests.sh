# UPF RW tests
python ../upf_data.py ref_upf/C.pbe-hgh.UPF C.res.UPF
python ../upf_data.py ref_upf/C-GTH-PBE-1.upf C.res.UPF
python ../upf_data.py ref_upf/Al_FR.upf Al.res.UPF  # soc


# GTH RW tests
python ../gth_data.py ref_gth/C_gth C.res.gth
python ../gth_data.py ref_gth/Au_gth Au.res.gth
python ../gth_data.py ref_gth/W-GTH-PBE-q14-SOC W.res.gth # soc

# GTH to UPF tests
python ../gth2upf_custom.py ref_gth/C_gth ref_upf/C.pbe-hgh.UPF # to be fixed, and soc cases to be added
python ../gth2upf_custom.py ref_gth/Al_gth ref_upf/Al_FR.upf # soc
python ../gth2upf_custom.py ref_gth/W-GTH-PBE-q14-SOC ref_upf/W-GTH-PBE-q14-SOC.upf # soc

#eval
python ../debug/plot_check_eval.py 
