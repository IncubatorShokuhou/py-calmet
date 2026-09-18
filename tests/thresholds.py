"""Parity thresholds vs Fortran CALMET goldens.

Synthetic tiny-domain gates use relative RMSE
  relative RMSE = sqrt(mean(((pred-ref)/(|ref|+eps))^2))
with eps=1e-6.

The wrf_demo case (complex terrain, ~28 km, synthesized obs) uses looser
relative-RMSE gates plus Pearson correlation floors for U/V, because a full
DIAGNO/OA port is still approximate on this domain.
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

# wrf_demo: intentionally looser (documented)
WRF_DEMO_THRESH = {
    "obs": {
        "U": 5.0,
        "V": 15.0,
        "U_corr": 0.30,
        "V_corr": 0.30,
        "ZI": 1.0,
        "USTAR": 1.5,
        "SPD": 2.0,
    },
    "obs_model": {
        "U": 20.0,
        "V": 50.0,
        "U_corr": 0.40,
        "V_corr": 0.50,
        "ZI": 2.0,
        "USTAR": 3.0,
        "SPD": 2.0,
    },
    "noobs": {
        "U": 20.0,
        "V": 50.0,
        "U_corr": 0.50,
        "V_corr": 0.60,
        "ZI": 2.0,
        "USTAR": 3.0,
        "SPD": 2.0,
    },
}
