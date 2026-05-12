# Data and Code for "What Influences the Field Goal Attempts of Professional Players? Analysis of Basketball Shot Charts via Log Gaussian Cox Processes with Spatially Varying Coefficients"

## Overview

This archive contains the datasets and R code used in the real data analysis of basketball shot chart data for three NBA players: Stephen Curry (2014–15 season), LeBron James (2014–15 season), and Michael Jordan (2001–02 season). The analysis applies a spatially varying coefficient log-Gaussian Cox process (JSVLGCP) model to study how shot intensity patterns depend on game location (home vs. away) and opponent strength (strong vs. weak).

---

## Directory Structure

```
Supp_Data/
├── README.md
├── Code/
│   ├── basket_func_same.R           # GP eigenfunction helper functions
│   ├── SVC-LGCP.R                   # Core JSVLGCP model fitting function
│   ├── curry_homelevel.R            # MCMC model fitting script for Curry
│   ├── james_homelevel.R            # MCMC model fitting script for James
│   ├── jordan_homelevel.R           # MCMC model fitting script for Jordan
│   ├── ModelComparison.ipynb        # Model comparison via p-thinning (Section 5)
│   ├── ReplotManuscript.ipynb       # Manuscript figure generation (Section 5)
│   └── Figures/                     # Generated figures from ReplotManuscript.ipynb 
└── Data/
    ├── del_curry.txt                # Shot data for Curry 
    ├── curry_covariates.txt         # Game-level covariates for Curry
    ├── del_james.txt                # Shot data for James 
    ├── james_covariates.txt         # Game-level covariates for James
    ├── jordan_data.txt              # Shot data for Jordan 
    ├── jordan_covariates.txt        # Game-level covariates for Jordan
    ├── curry_homelevel_theta.RData  # Posterior MCMC samples for Curry
    ├── james_homelevel_theta.RData  # Posterior MCMC samples for James
    ├── jordan_homelevel_theta.RData # Posterior MCMC samples for Jordan
    ├── RR_significant_curry.RData   # Credible interval significance indicators (Curry)
    ├── RR_significant_james.RData   # Credible interval significance indicators (James)
    └── RR_significant_jordan.RData  # Credible interval significance indicators (Jordan)

```

---

## Data Description

### Shot Data

#### `del_curry.txt` / `del_james.txt`
Shot-level data for Stephen Curry and LeBron James used as direct model input.

| Variable   | Description |
|------------|-------------|
| `Game`     | Sequential game index (1 to number of games in season; links to covariates file) |
| `Missmade` | Shot outcome: 1 = made, 0 = missed |
| `Xdist`    | Horizontal coordinate on court (pixels; origin at basket center) |
| `Ydist`    | Vertical coordinate on court (pixels; origin at baseline) |
| `Dist`     | Distance from basket in feet |
| `shot`     | Shot type description (e.g., "Jump Shot") |

Before model fitting, coordinates are normalized to [−1, 1] × [−1, 1]: `Xdist` is divided by 250; `Ydist` is divided by 235 and shifted by −39/47.

#### `jordan_data.txt`
Shot-level data for Michael Jordan used as direct model input, compiled from individual game records.

| Variable   | Description |
|------------|-------------|
| `Game`     | Sequential game index |
| `Shots`    | Shot index within game |
| `Missmade` | Shot outcome: 1 = made, 0 = missed |
| `Xdist`    | Horizontal coordinate (pixels from basket center) |
| `Ydist`    | Vertical coordinate (pixels from baseline) |
| `Dist`     | Distance from basket in feet |
| `Quarter`  | Quarter of game (1–4) |

Before model fitting, shots with `Dist` > 28 or `Dist` < 1 are excluded; `Xdist` is divided by 300 and shifted by −1; `Ydist` is divided by 235 and shifted by −1, mapping coordinates to [−1, 1] × [−1, 1].

---

### Covariate Data

#### `curry_covariates.txt` / `james_covariates.txt` / `jordan_covariates.txt`
Game-level covariate files. Each row corresponds to one game.

| Variable | Description |
|----------|-------------|
| `game`   | Sequential game index (matches `Game` in shot data) |
| `home`   | Game location: 1 = home, 2 = away |
| `level`  | Opponent strength: 1 = strong opponent, 2 = weak opponent |


### Posterior Samples and Analysis Results

#### `curry_homelevel_theta.RData` / `james_homelevel_theta.RData` / `jordan_homelevel_theta.RData`
Posterior MCMC samples from the fitted JSVLGCP model for each player. Each file contains the array `Theta.result` of dimension 15,000 × L × 7, where L = 15 is the number of GP basis functions and 7 corresponds to the intercept plus six covariate-specific coefficient processes (home, away, level for missed shots; home, away, level for made shots). The last 5,000 iterations (rows 10,001–15,000) are used as posterior samples after burn-in.

#### `RR_significant_curry.RData` / `RR_significant_james.RData` / `RR_significant_jordan.RData`
Pointwise credible interval significance indicators for the spatial relative risk maps, computed from `Theta.result`.


---

## Code Description

### Helper Functions

#### `basket_func_same.R`
Core R functions implementing the Gaussian process (GP) eigenfunction approximation used in the JSVLGCP model. 

#### `SVC-LGCP.R`
Wraps the JSVLGCP MCMC sampler into the function `fit_SVCLGCP(shot_data, covar)`, used by `ModelComparison.ipynb` to fit the model on p-thinned training data across replications. Returns a list containing posterior samples (`Theta.result`) and diagnostic quantities.

---

### Model Fitting Scripts

#### `curry_homelevel.R` / `james_homelevel.R` / `jordan_homelevel.R`
Standalone MCMC scripts implementing the JSVLGCP model for each player.

---

### Analysis Notebooks

#### `ModelComparison.ipynb`
Jupyter notebook (R kernel) implementing the p-thinning model comparison reported in Section 5.


#### `ReplotManuscript.ipynb`
Jupyter notebook (R kernel) that reproduces all figures in Section 5 of the manuscript from the fitted JSVLGCP results. All required files are included in `Data/`. 

---

### Software Requirements

- **R** (≥ 4.0) with necessary packages.
- **JupyterLab** with an R kernel (via `IRkernel`) for running the `.ipynb` notebooks
