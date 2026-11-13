#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Empirical fallback.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Audited empirical permutation fallback for a non-finite Beta approximation.

The original permutation maxima and global conditional cutoff are retained.
Floating point ties within 1e-6 in r-squared are included conservatively.
This is a fallback P value, never a fitted Beta P value.
"""

import numpy as np
from scipy.stats import t


def fallback(r2_perm, r2_nominal, dof, threshold, tolerance=1e-6):
    values = np.asarray(r2_perm, dtype=np.float64)
    observed = float(r2_nominal)
    assert values.ndim == 1 and len(values) > 0 and np.isfinite(values).all()
    assert np.isfinite(observed) and dof > 0
    assert values.min() >= -tolerance and values.max() <= 1 + tolerance
    assert -tolerance <= observed <= 1 + tolerance
    assert threshold >= 1 / (
        len(values) + 1
    ), "More permutations required: cutoff below empirical resolution"
    obs = np.clip(observed, 0, 1)
    exceed = int((np.clip(values, 0, 1) >= obs - tolerance).sum())
    p = (exceed + 1) / (len(values) + 1)
    nominal = 0.0 if obs == 1 else float(2 * t.sf(np.sqrt(dof * obs / (1 - obs)), dof))
    metadata = dict(
        method="empirical_fallback",
        nperm=len(values),
        exceed=exceed,
        tie_tolerance=tolerance,
        r2_nominal_raw=observed,
        r2_perm_min=float(values.min()),
        r2_perm_max=float(values.max()),
        p_conditional=p,
        threshold=float(threshold),
    )
    return (p, np.nan, np.nan, float(dof), nominal), metadata
