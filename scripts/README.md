
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Analysis scripts used in “Ringed seal *(Pusa hispida)* haul-out behavior and emergence timing in the Bering, Chukchi, and Beaufort seas”

## Draft Manuscript Under Active Development

Please note this manuscript associated with this repository is still in
peer review and not final. Changes to code and results are still
possible.

## Contents

This repository contains:

- [:page_facing_up:
  Demographic_comparisons.R](/scripts/Demographic_comparisons.R): code
  to compare seasonal patterns in haul-out behavior by sex & age,
  including generating manuscript figures 4 and S1.

- [:page_facing_up: HMMs_Adults.R](/scripts/HMMs_Adults.R): code to fit
  Hidden Markov Models (HMMs) to haul-out data from adult ringed seals,
  estimate emergence dates for individual seals, and generate the
  corresponding figures (figs 5,6,S2) and summary statistics reported in
  the manuscript. The model fitting code in this script is heavily based
  on the manual and associated vignettes for the momentuHMM package
  (<https://cran.r-project.org/web/packages/momentuHMM/index.html> and
  <https://github.com/bmcclintock/momentuHMM>)

- [:page_facing_up: HMMs_Subadults.R](/scripts/HMMs_Subadults.R): code
  to fit Hidden Markov Models (HMMs) to haul-out data from subadult
  ringed seals, estimate emergence dates for individual seals, and
  generate the corresponding figures (figs 7,S3) and summary statistics
  reported in the manuscript.

- [:page_facing_up:
  LOO_emergence_dates.R](/scripts/LOO_emergence_dates.R): code to run a
  leave-one-out sensitivy analysis on ringed seal emergence dates
  estimated from the best HMM for each age class. Also generates
  manuscript figures S4 & S5.
