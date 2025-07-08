
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Data used in “Ringed seal *(Pusa hispida)* haul-out behavior and emergence timing in the Bering, Chukchi, and Beaufort seas”

This repository contains:

- [:page_facing_up:
  ringed_ho_data_clean.csv](/ringed_ho_data_clean.csv): cleaned daily
  haul-out and environmental data ready to be analyzed in Hidden Markov
  Models

## Haul-out data sources

The haul-out data herein are from previous tagging projects in the
Bering, Chukchi, and Beaufort seas conducted by the National Oceanic and
Atmospheric Administration (NOAA; Kelly et al. 2010a), the Alaska
Department of Fish and Game (ADF&G; Crawford et al. 2012, Quakenbush et
al. 2019), and the North Slope Borough (NSB; Von Duyke et al. 2020,
unpubl. data). Each row represents one day of haul-out data for an
individual seal, including its daily estimated location and attached
environmental covariates. Please refer to the “data_source” column to
see which agency the haul-out data for each seal is associated with.

<div style="text-indent: -40px; padding-left: 40px;">

Crawford JA, Frost KJ, Quakenbush LT, Whiting A (2012) Different habitat
use strategies by subadult and adult ringed seals *(Phoca hispida)* in
the Bering and Chukchi seas. Polar Biol 35:241–255.

Kelly BP, Badajos OH, Kunnasranta M, Moran JR, Martinez-Bakker M,
Wartzok D, Boveng P (2010a) Seasonal home ranges and fidelity to
breeding sites among ringed seals. Polar Biol 33:1095–1109.

Quakenbush L, Bryan A, Crawford J, Olnes J (2020) Biological monitoring
of ringed seals in the Bering and Chukchi seas. Final Report to National
Oceanographic and Atmospheric Association, Award Number: NA16NMF4720079.
Fairbanks, AK. 33 pp.

Von Duyke AL, Douglas DC, Herreman JK, Crawford JA (2020) Ringed seal
*(Pusa hispida)* seasonal movements, diving, and haul-out behavior in
the Beaufort, Chukchi, and Bering Seas (2011-2017). Ecol Evol
10:5595–5616.

</div>

## Column descriptors

### Seal haul-out data

- speno_yr: unique identifier for each combination of individual seal
  and year. E.g., “RS07-12-F_2008” refers to haul-out data from seal
  RS07-12-F transmitted during February-June 2008.

- age_class: “Adult” or “Subadult”

- sex: “Female” or “Male”

- data_source: the agency that deployed the satellite-linked transmitter
  for that individual seal. “NOAA”, “ADFG”, or “NSB”.

- daily_prop_ho: daily proportion hauled out (0 to 1)

- peak_ho_hr: peak haul-out hour, in radians (-pi to pi)

- doy: day of year that haul-out data was collected (32 to 166, i.e., 1
  February to 15 June)

- year: year that haul-out data was collected (2005 to 2021)

- lat: latitude of seal’s daily estimated location

- lon: longitude of seal’s daily estimated location

### Environmental Covariates

- daylength: daylength for the seal’s latitude and doy

- iceconc: daily 25-km sea-ice concentration from National Snow & Ice
  Data Center SSM/I passive microwave data
  (<https://nsidc.org/data/nsidc-0051>; Cavalieri et al. 1996, updated
  yearly)

- air2m: daily 32-km mean air temperature (°C at 2 m above sea level)
  from the North American Regional Reanalysis (NARR) dataset
  (<https://psl.noaa.gov/data/gridded/data.narr.html>; Mesinger et
  al. 2006)

- rolling_temp: 7-day rolling average of air2m within each 32-km grid
  cell

- atdd: accumulated thawing degree days, i.e., number of days since 1
  April with a mean daily temperature of ≥0 °C

- CM_index: continuous melt index, derived from continuous melt onset
  dates from the NASA Cryosphere database
  (<https://earth.gsfc.nasa.gov/cryo/data/arctic-sea-ice-melt>; Markus
  et al. 2009)

- EM_index: early melt index, derived from early melt onset dates from
  the NASA Cryosphere database
  (<https://earth.gsfc.nasa.gov/cryo/data/arctic-sea-ice-melt>; Markus
  et al. 2009)

### Acknowledgements

We thank the Inuit communities, Brendan Kelly, and other field
researchers who contributed to ringed seal tagging in the original
studies. Tagging efforts were conducted under NMFS Research Permits Nos.
358–1787-01, 15324, and 20466 issued to ADF&G; 782-1694-00 issued to
NOAA; and Scientific License Nos. SLE04/05-328 and SLE–05/06-322 from
the Dept. of Fisheries and Oceans (DFO), Canada. Animal Care and Use
Protocols were issued by the DFO (UFWI-ACC-2004-2005-001U) and the ADF&G
Animal Care and Use Committee (Nos. 06–16, 2010-13R, 2014-03, 2016-23,
0027-2017-27, 0027-2018-29, and 0027-2019-41). Funding for tagging was
provided by the Bureau of Ocean Energy Mgmt., the Office of Naval
Research, the North Slope Borough (NSB)-Shell Baseline Studies Program,
Shell Exploration and Production Co., National Fish and Wildlife
Foundation with funding from Conoco Phillips, the Native Village of
Kotzebue, NOAA, NMFS Alaska Region, and a U.S. Dept. of the Interior,
Tribal Wildlife Grant for Federally Recognized Tribes.
