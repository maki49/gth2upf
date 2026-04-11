from pathlib import Path

from upf_data import read_upf_file
from upf_data import write_upf_v2


ROOT = Path(__file__).resolve().parents[1]


def test_parse_spin_orbit_metadata():
    upf = read_upf_file(ROOT / "results/W-GTH-SOC-PBE-q14-cp2k-default-grid.upf")

    assert upf.has_so is True
    assert upf.jjj == [
        0.5,
        0.5,
        0.5,
        0.5,
        1.5,
        0.5,
        1.5,
        0.5,
        1.5,
        1.5,
        2.5,
        1.5,
        2.5,
    ]
    assert upf.jchi == [0.5, 0.5, 0.5, 1.5, 1.5, 2.5]


def test_normalize_soc_wavefunction_occupations():
    w_upf = read_upf_file(ROOT / "results/W-GTH-SOC-PBE-q14-cp2k-default-grid.upf")
    s_upf = read_upf_file(ROOT / "results/S-GTH-SOC-PBE-q6-cp2k-default-grid.upf")

    w_upf.canonicalize_soc_wfc_order()
    s_upf.canonicalize_soc_wfc_order()
    w_upf.normalize_soc_wfc_occupations()
    s_upf.normalize_soc_wfc_occupations()

    assert w_upf.jchi == [0.5, 0.5, 1.5, 0.5, 2.5, 1.5]
    assert s_upf.jchi == [0.5, 1.5, 0.5]

    assert w_upf.oc.tolist() == [2.0, 2.0, 4.0, 2.0, 2.4, 1.6]
    assert s_upf.oc.tolist() == [2.0, 8.0 / 3.0, 4.0 / 3.0]

    tmp_file = ROOT / "test" / "soc_w_tmp.upf"
    try:
        write_upf_v2(w_upf, tmp_file)
        roundtrip = read_upf_file(tmp_file)
        assert roundtrip.jchi == [0.5, 0.5, 1.5, 0.5, 2.5, 1.5]
        assert roundtrip.oc.tolist() == [2.0, 2.0, 4.0, 2.0, 2.4, 1.6]
    finally:
        if tmp_file.exists():
            tmp_file.unlink()


if __name__ == "__main__":
    test_parse_spin_orbit_metadata()
    test_normalize_soc_wavefunction_occupations()
    print("ok")
