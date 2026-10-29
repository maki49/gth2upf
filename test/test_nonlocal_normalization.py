"""Check GTH nonlocal operators independently of projector normalization."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from eval_gth2upf import combine_channels, combine_channels_soc
from gth_tools import parse_gth_pp, p_il
from upf_data import read_upf_file, write_upf_v2

ROOT = Path(__file__).resolve().parents[1]


def check_distributed_dion_matches_cp2k(element, valence, soc):
    stem = f"{element}-GTH{soc}-PBE-q{valence}"
    evaluated = read_upf_file(ROOT / "results" / f"{stem}.upf")
    reference, = (ROOT / "results").glob(f"{stem}-cp2k-default-grid*.upf")
    cp2k = read_upf_file(reference)
    np.testing.assert_allclose(evaluated.dion, cp2k.dion, rtol=1e-10, atol=1e-10)


def check_grid_replacement_preserves_nonlocal_operator(gth_file, upf_file, tmp_path):
    content = (ROOT / gth_file).read_text()
    gth = parse_gth_pp(content)
    upf = read_upf_file(ROOT / upf_file)
    original_dion = upf.dion.copy()
    h = (combine_channels_soc(gth.channels, gth.soc_channels) if upf.has_so
         else combine_channels(gth.channels)[0])
    upf.replace_grid_by_eval(content)
    np.testing.assert_allclose(upf.dion, original_dion, rtol=1e-10, atol=1e-10)

    # A smooth radial probe checks the actual operator: in Ry it must be
    # 2 * <u|r*p_i> h_ij <r*p_j|u>, including SOC and off-diagonal terms.
    radial = []
    for channel in gth.channels:
        for i in range(channel.m):
            p = upf.grid.r * p_il(upf.grid.r, channel.l, i + 1, channel.r_l)
            radial.append(p)
            if upf.has_so and channel.l >= 1:
                radial.append(p)
    probe = upf.grid.r * np.exp(-upf.grid.r)
    weights = probe * upf.grid.rab
    overlaps = np.asarray(radial) @ weights
    beta_overlaps = np.asarray(upf.beta) @ weights
    np.testing.assert_allclose(beta_overlaps @ upf.dion @ beta_overlaps,
                               2 * overlaps @ h @ overlaps, rtol=1e-12)

    output = tmp_path / "roundtrip.upf"
    write_upf_v2(upf, output)
    np.testing.assert_allclose(read_upf_file(output).dion, h / 2, atol=1e-10)


class TestNonlocalNormalization(unittest.TestCase):
    def test_distributed_upfs(self):
        for element, valence in [("S", 6), ("Se", 6), ("W", 14)]:
            for soc in ["", "-SOC"]:
                with self.subTest(element=element, soc=soc):
                    check_distributed_dion_matches_cp2k(element, valence, soc)

    def test_grid_replacement(self):
        for gth_file, upf_file in [
            ("test/ref_gth/C_gth", "results/C-GTH-PBE-q4-cp2k-default-grid.upf"),
            ("test/ref_gth/W-GTH-PBE-q14-SOC", "results/W-GTH-SOC-PBE-q14-cp2k-default-grid.upf"),
        ]:
            with self.subTest(gth_file=gth_file), tempfile.TemporaryDirectory() as directory:
                check_grid_replacement_preserves_nonlocal_operator(gth_file, upf_file, Path(directory))


if __name__ == "__main__":
    unittest.main()
