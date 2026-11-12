"""Where each input dataset lives and how to cite it.

Method and provenance: docs/METHODS.md. Copernicus Marine dataset IDs and variable names were read
from the live catalogue with ``copernicusmarine.describe`` on 2026-11-09.
"""

GREP = dict(
    title="Global Ocean Ensemble Physics Reanalysis (GREP), GLOBAL_MULTIYEAR_PHY_ENS_001_031",
    members_dataset="cmems_mod_glo_phy-all_my_0.25deg_P1D-m",
    # The three reanalyses in GREP (QUID 001-031 Issue 1.3, Table 1). FOAM was removed (same QUID, section V).
    members={"glor": "GLORYS2V4 (Mercator Ocean)", "cglo": "C-GLORS (CMCC)", "oras": "ORAS5 (ECMWF)"},
    doi="10.48670/moi-00024",
    citation="Global Ocean Ensemble Physics Reanalysis. E.U. Copernicus Marine Service Information (CMEMS). "
             "Marine Data Store (MDS). DOI: 10.48670/moi-00024",
)

ARGO = dict(
    title="Argo float data and metadata from the Global Data Assembly Centre (Argo GDAC)",
    doi="10.17882/42182",
    citation="Argo (2000). Argo float data and metadata from Global Data Assembly Centre (Argo GDAC). SEANOE. "
             "doi:10.17882/42182",
    acknowledgement="These data were collected and made freely available by the International Argo Program and "
                    "the national programs that contribute to it (https://argo.ucsd.edu, "
                    "https://www.ocean-ops.org). The Argo Program is part of the Global Ocean Observing System.",
    # Temperature accuracy after delayed-mode adjustment (Wong et al. 2020, Front. Mar. Sci. 7:700, abstract).
    temperature_accuracy_K=0.002,
)

RG = dict(
    title="Roemmich-Gilson Argo Climatology (temperature), 1 degree monthly objective analysis",
    url="https://sio-argo.ucsd.edu/RG_Climatology.html",
    citation="Roemmich, D. and J. Gilson (2009). The 2004-2008 mean and annual cycle of temperature, salinity, and "
             "steric height in the global ocean from the Argo Program. Prog. Oceanogr., 82, 81-100. "
             "doi:10.1016/j.pocean.2009.03.004",
)
