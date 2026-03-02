# Analyzing Spatial Correlation Patterns of Water Quality Variables in Fluvial Networks
Authors: N. Pronello, S. Castiglia, V. Frontuto, N. Golini, R. Ignaccolo, L. Ippoliti

The following repository containes the R code and the data useful to analyze spatial correlation patterns of water quality variables in fluvial networks.

📁 Folder Structure
R/

This folder contains two R scripts:

1) Functions4covAnalysis.R
This file includes a set of functions useful for:

implementing cross-validation procedures,

computing Betti number 0 curves,

calculating the Kolmogorov–Smirnov (KS) distance between curves.

2) to_run.R
This script allows the user to:

estimate covariance-valued functions over a network,

compute the corresponding Betti curves,

cluster the resulting Betti curves.

data/

This folder contains supplementary files used to:

1) graphically represent the domain of interest,

2) provide covariance matrices useful for generating the synthetic data used in to_run.R.
