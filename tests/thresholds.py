"""Parity thresholds vs Fortran CALMET goldens — wrf_demo (real WRF) only.

The wrf_demo case (NCAR Katrina wrfout mountain window + GEO.DAT terrain,
~28 km) uses floored relative RMSE (floor=0.5 m/s on U/V) plus Pearson
correlation floors, because raw relative RMSE is dominated by near-zero
reference cells even when absolute errors are ~0.5–1.3 m/s.

Synthetic ``small_domain`` / ``daytime_zi`` gates were removed; this file
is the sole Fortran-compare threshold table.
"""

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
