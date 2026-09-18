"""Parity thresholds vs Fortran CALMET goldens.

Synthetic tiny-domain gates use relative RMSE
  relative RMSE = sqrt(mean(((pred-ref)/(|ref|+eps))^2))
with eps=1e-6.

The wrf_demo case (complex terrain, ~28 km) uses floored relative RMSE
(floor=0.5 m/s on U/V) plus Pearson correlation floors, because raw
relative RMSE is dominated by near-zero reference cells even when
absolute errors are ~0.5–1.3 m/s.
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

# wrf_demo: floored relative RMSE (see WRF_DEMO_UV_FLOOR) + correlation
WRF_DEMO_UV_FLOOR = 0.5  # m/s
WRF_DEMO_THRESH = {
    "obs": {
        "U": 1.2,
        "V": 2.0,
        "U_corr": 0.40,
        "V_corr": 0.55,
        "ZI": 0.50,
        "USTAR": 0.80,
        "SPD": 1.5,
    },
    "obs_model": {
        "U": 0.80,
        "V": 2.0,
        "U_corr": 0.90,
        "V_corr": 0.90,
        "ZI": 1.0,
        "USTAR": 1.5,
        "SPD": 1.0,
    },
    "noobs": {
        "U": 0.80,
        "V": 2.0,
        "U_corr": 0.90,
        "V_corr": 0.90,
        "ZI": 1.0,
        "USTAR": 1.5,
        "SPD": 1.0,
    },
}

# daytime_zi: Maul–Carson path (simplified gamma=DPTMIN; growth faster than
# Fortran sounding-based lapse). Gate on QSW, ZI correlation, convective flags.
DAYTIME_ZI_THRESH = {
    "QSW": 0.05,
    "ZI_corr_day": 0.90,
    "ZI_night": 0.15,
    "U": 0.05,
    "V": 0.05,
}
