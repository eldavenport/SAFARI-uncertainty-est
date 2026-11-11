"""Where each input dataset lives and how to cite it.

Every entry is documented in docs/PRODUCT_CATALOGUE.md. Copernicus Marine dataset IDs and variable names were read
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

GLORYS12 = dict(
    title="Global Ocean Physics Reanalysis GLORYS12V1, GLOBAL_MULTIYEAR_PHY_001_030",
    dataset="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    doi="10.48670/moi-00021",
    citation="Lellouche, J.-M., et al. (2021). The Copernicus Global 1/12 degree Oceanic and Sea Ice GLORYS12 "
             "Reanalysis. Front. Earth Sci., 9, 698876. doi:10.3389/feart.2021.698876",
)

OSTIA = dict(
    title="OSTIA reprocessed sea surface temperature analysis, SST_GLO_SST_L4_REP_OBSERVATIONS_010_011",
    dataset="METOFFICE-GLO-SST-L4-REP-OBS-SST",
    variable="analysis_error",
    doi="10.48670/moi-00168",
    citation="Good, S., et al. (2020). The current configuration of the OSTIA system for operational production "
             "of foundation sea surface temperature and ice concentration analyses. Remote Sens., 12, 720. "
             "doi:10.3390/rs12040720",
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

WOA23 = dict(
    title="World Ocean Atlas 2023, temperature",
    url="https://www.ncei.noaa.gov/thredds-ocean/fileServer/woa23/DATA/temperature/netcdf/{period}/{res}/"
        "woa23_{period}_t{month:02d}_{grid}.nc",
    # File-name grid codes (WOA23 Product Documentation, Table 2): 04 = 0.25 degree, 01 = 1 degree.
    grid_codes={"0.25": "04", "1.00": "01"},
    doi="10.25923/54bh-1613",
    citation="Locarnini, R. A., et al. (2024). World Ocean Atlas 2023, Volume 1: Temperature. NOAA Atlas NESDIS 89. "
             "doi:10.25923/54bh-1613",
)
