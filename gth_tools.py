# Author: Jincheng Yu <pimetamon@gmail.com>
import re
import sys
import numpy as np
from math import lgamma
from dataclasses import dataclass
from typing import Dict, List, Optional
import xml.etree.ElementTree as ET

def parse_gth_pp(text):
    lines = [ln.rstrip() for ln in text.strip().splitlines() if ln.strip()]
    if len(lines) < 3:
        raise ValueError("Too short for a GTH PP block.")

    # Header:  "Cu GTH-PBE-q11 GTH-PBE"
    m = re.match(r"^([A-Z][a-z]?)\s+(\S*?)-q(\d+)(?:\s+\S+)?$", lines[0])
    if not m:
        raise ValueError(f"Bad header: {lines[0]}")
    element, xc_label, q_hdr = m.group(1), m.group(2), int(m.group(3))

    # --- Local part line(s) ---
    # Expect: line 1 is metadata (varies by generator), line 2 begins local line
    # Format: r_loc, nC, then nC coefficients (which may continue on following lines)
    lp0 = lines[2].split()
    if len(lp0) < 2:
        raise ValueError(f"Local line too short: {lines[2]}")
    r_loc = float(lp0[0])
    nC = int(float(lp0[1]))

    # Collect C coefficients, possibly spanning multiple continuation lines
    C_tokens = lp0[2:]
    extra_used = 0
    while len(C_tokens) < nC:
        extra_used += 1
        if 2 + extra_used >= len(lines):
            raise ValueError("Unexpected EOF while reading local C coefficients")
        C_tokens.extend(lines[2 + extra_used].split())

    Cs = [float(x) for x in C_tokens[:nC]]
    # Pad to 4 (GTH local uses up to 4 coefficients)
    while len(Cs) < 4:
        Cs.append(0.0)

    local = LocalPart(Z_ion=float(q_hdr), r_loc=r_loc, C=Cs[:4])

    # --- Determine start index of nonlocal section ---
    idx = 3 + extra_used  # next unread line after the local block

    # --- Number of l-channels (n_l) ---
    # If the next line is a single integer token, use it. Otherwise, fall back to the
    # previous heuristic (last token of line 1), which works for many GTH files.
    def _looks_int_token(s: str) -> bool:
        try:
            v = float(s)
            return abs(v - int(v)) < 1e-12
        except Exception:
            return False

    n_l = None
    if idx < len(lines):
        tks = lines[idx].split()
        if len(tks) == 1 and _looks_int_token(tks[0]):
            n_l = int(float(tks[0]))
            idx += 1

    if n_l is None:
        # Fallback: read from the last token of line 1
        tks = lines[1].split()
        if not tks or not _looks_int_token(tks[-1]):
            raise ValueError("Cannot determine number of l-channels (n_l).")
        n_l = int(float(tks[-1]))

    # --- Nonlocal channels ---
    channels: List[Channel] = []
    for l in range(n_l):
        if idx >= len(lines):
            raise ValueError(f"Unexpected EOF before channel l={l}")

        # First line of channel l: r_l, m, then some (or none) of the h_ij
        head = lines[idx].split()
        if len(head) < 2:
            raise ValueError(f"Channel header too short at line: {lines[idx]}")

        r_l = float(head[0])
        m_proj = int(float(head[1]))
        need = m_proj * (m_proj + 1) // 2  # number of symmetric h_ij (upper triangle)

        coeff_tokens = head[2:]
        used = 0
        while len(coeff_tokens) < need:
            used += 1
            if idx + used >= len(lines):
                raise ValueError(f"Unexpected EOF while reading h_ij for l={l}")
            coeff_tokens.extend(lines[idx + used].split())

        coeffs = [float(x) for x in coeff_tokens[:need]]
        if len(coeffs) != need:
            raise ValueError(f"l={l}: expected {need} h_ij, got {len(coeffs)}")

        channels.append(Channel(l=l, r_l=r_l, m=m_proj, h_ij=coeffs))
        idx += 1 + used  # advance past this block (first line + continuation lines)

    return GTHPP(element, xc_label, q_hdr, local, channels)



def V_loc(r, Z_ion=4.0, r_loc=0.33847124, C1=-8.80367398, 
          C2=1.33921085, C3=0.0, C4=0.0):
    """
    Evaluate
        V_loc(r) = -Z_ion/r * erf(r/(sqrt(2)*r_loc))
                   + exp(-0.5*(r/r_loc)**2) * [ C1 + C2*(r/r_loc)**2
                                               + C3*(r/r_loc)**4 + C4*(r/r_loc)**6 ]
    Supports scalar or array r. Uses the finite r->0 limit at r == 0.

    Parameters
    ----------
    r : float or array_like
        Radius (same units as r_loc).
    Z_ion : float
    r_loc : float
    C1, C2, C3, C4 : floats

    Returns
    -------
    ndarray or float
        V_loc(r) with the same shape as r.
    """
    r = np.asarray(r, dtype=float)
    x = r / r_loc

    # erf: prefer numpy's ufunc if present; 
    # otherwise fall back to a vectorized math.erf
    try:
        erf = np.erf  # available in modern NumPy (np.special.erf)
    except AttributeError:
        from math import erf as _erf
        erf = np.vectorize(_erf, otypes=[float])

    term_poly = np.exp(-0.5 * x**2) * (C1 + C2 * x**2 + C3 * x**4 + C4 * x**6)

    # General expression (will be overridden at r==0 by the limit below)
    term_erf = -Z_ion / r * erf(x / np.sqrt(2))

    # Finite limit at r -> 0 for the first term: 
    #  -Z_ion * sqrt(2)/(sqrt(pi) * r_loc)
    term0_limit = -Z_ion * np.sqrt(2.0) / (np.sqrt(np.pi) * r_loc)

    out = np.where(r == 0.0, term0_limit + C1, term_erf + term_poly)
    return out

def p_il(r, l, i, r_l):
    """GTH radial projector p_i^(l)(r), normalized."""
    if l < 0 or i < 1 or r_l <= 0:
        raise ValueError("Require l>=0, i>=1, r_l>0.")
    r = np.asarray(r, dtype=float)
    power_r = l + 2*(i-1)
    alpha = l + (4*i - 1)/2.0
    log_norm = 0.5*np.log(2.0) - alpha*np.log(r_l) - 0.5*lgamma(alpha)
    prefac = np.exp(log_norm)
    return prefac * r**power_r * np.exp(-0.5*(r/r_l)**2)

def pp_beta_sum_l(gth, r_mesh: np.ndarray, l: int, i_max: int = 3):
    """
    gth: parsed GTHPP object (from the parser we wrote)
    r_mesh: 1D array of radii
    l: angular momentum channel to evaluate
    i_max: highest projector index to include (default 3)
    Returns: 1D array with PP_BETA^(l)(r) = sum_i p_i^(l)(r)
    """
    # find the channel with this l
    ch = next((c for c in gth.channels if c.l == l), None)
    if ch is None or ch.m == 0:
        return np.zeros_like(r_mesh, dtype=float)
    upto = min(i_max, ch.m)
    out = np.zeros_like(r_mesh, dtype=float)
    for i in range(1, upto+1):
        out += p_il(r_mesh, l, i, ch.r_l)
    return out

def normalize_p_il(r, p_vals, measure="r2", weight=None, return_norm=False):
    """
    Numerically normalize a radial projector on a grid.

    Parameters
    ----------
    r : (Nr,) array
        Radial grid (monotonic, same units as the projector).
    p_vals : (Nr,) array
        Values of p_i^{(l)}(r) on the grid.
    measure : {"r2", "plain"} or None
        - "r2": weight w(r) = r^2  (the 3D integral measure)  [default]
        - "plain": weight w(r) = 1 (i.e., ∫ p^2 dr = 1)
        - None: use custom 'weight' array passed below.
    weight : (Nr,) array or None
        Custom weight w(r); if given, overrides 'measure'.
    return_norm : bool
        If True, also return the computed norm before scaling.

    Returns
    -------
    p_norm : (Nr,) array
        Normalized projector on the same grid.
    (optional) norm : float
        sqrt( ∫ p^2 w(r) dr ).
    """
    r = np.asarray(r, dtype=float)
    p = np.asarray(p_vals, dtype=float)

    if weight is None:
        if measure == "r2":
            w = r**2
        elif measure == "plain":
            w = np.ones_like(r)
        else:
            raise ValueError("measure must be 'r2', 'plain', or provide weight=...")
    else:
        w = np.asarray(weight, dtype=float)

    # numerical norm
    norm_sq = 0.25 * np.sum(p**2 * w)
    #norm_sq = sum(p**2 * r**2 * w)
    if norm_sq <= 0:
        raise ValueError("Non-positive norm; check inputs.")
    norm = np.sqrt(norm_sq)
    p_out = p / norm
    return (p_out, norm) if return_norm else p_out


@dataclass
class LocalPart:
    Z_ion: float
    r_loc: float
    C: List[float]  # [C1..C4]

@dataclass
class Channel:
    l: int
    r_l: float
    m: int
    h_ij: List[float]

@dataclass
class GTHPP:
    element: str
    xc_label: str
    q: int
    local: LocalPart
    channels: List[Channel]

def _is_int_line(s: str) -> bool:
    t = s.split()
    if len(t) != 1:
        return False
    try:
        int(float(t[0])); return True
    except ValueError:
        return False

