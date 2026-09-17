"""Parity thresholds vs Fortran CALMET goldens (relative RMSE).

Core fields (U, V, ZI, USTAR) are gated. Thresholds are intentionally in the
1e-3 .. ~1e-1 band: layer-1 winds and PBL for obs are tight; full 3-D winds
for obs/obs_model are looser pending a full DIAGNO/OA port. See PROGRESS.md.
"""

# relative RMSE = sqrt(mean(((pred-ref)/(|ref|+eps))^2))
THRESH = {
    "obs": {
        "U": 0.10,
        "V": 0.15,
        "U_lev1": 0.01,
        "V_lev1": 0.01,
        "ZI": 0.02,
        "USTAR": 0.01,
        "SPD": 0.05,
    },
    "noobs": {
        "U": 0.05,
        "V": 0.05,
        "U_lev1": 0.05,
        "V_lev1": 0.05,
        "ZI": 0.10,
        "USTAR": 0.05,
        "SPD": 0.05,
    },
    "obs_model": {
        "U": 0.20,
        "V": 0.10,
        "U_lev1": 0.25,  # OA blend; surface pulled toward model IGF
        "V_lev1": 0.15,
        "ZI": 0.12,
        "USTAR": 0.08,
        "SPD": 0.15,
    },
}
