# Aquaculture with Stochastic Feeding Costs
This repository is complementing [TODO](https://www.arxiv.org/) for investigating the effect of stochastic feeding costs to optimal harvesting rules in aquaculture farms. In particular, we chose salmon farms for our investigation.

## Installation
There is no installation required, please download the code and run the main files. For Python you need all the necessary dependencies to run Tensorflow 2.x. There is also a Jupyter notebook available for testing the Python code in Colab.

[![Open In Colab TODO](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/googlecolab/colabtools/blob/main/Python/colab.ipynb)

## Code and Usage
The code is structured as follows: 

For the data preprocessing and calibration of the commodity models to the market data we used Matlab with its (Global) Optimization Toolbox. 

We provide the code for estimating the optimal harvesting times with a Least-Square Monte-Carlo (LSMC) regression in both languages [Matlab](Matlab/) and [Python](Python/). 
Everything concerning Deep Learning is only implemented in Python using Tensorflow 2.10.

### Python
Optimal Stopping:
1. Run [main](Python/main.m) with your custom model parameters to find the optimal harvesting time. Also returns a comparison of stopping rules with deterministic and stochastic feeding costs using a Deep Neural Network and a pathwise comparison.
2. Deep Optimal Stopping (Optional):
    1. Run [DeepOS_IndepFeeding](Python/DeepOS_IndepFeeding.py) to use solve the optimal stopping problem in the case of stochastic feeding costs.
    2. Run [DeepOS_ConstFeeding](Python/DeepOS_ConstFeeding.py) to use solve the optimal stopping problem in the case of constant feeding costs (different to deterministic feeding costs in this paper) to compare to 
    [Ewald et al. (2016)](https://doi.org/10.1093/ajae/aaw052).

### Matlab
Data Preprocessing:
1. Run [preprocessingRaw_salmon](Matlab/Data/preprocessingRaw_salmon.m) to preprocess salmon data.
2. Run [preprocessingRaw_soy](Matlab/Data/preprocessingRaw_soy.m) to preprocess soy data.
3. (Optinal) Run [paramCheck](Matlab/Data/paramCheck.m) to compare parameter estimates from Kalman filter and [Cortazar and Schwartz (2003)](https://doi.org/10.1016/S0140-9883(02)00096-8). For this you need to calibrate the models first and save the unobservable spots from the algorithms. Examples are provided.

Calibration:
- [Cortazar and Schwartz (2003)](https://doi.org/10.1016/S0140-9883(02)00096-8): 
    1. Run [mainSalmon](Matlab/CortazarSchwartz/mainSalmon.m) to calibrate the salmon model.
    2. Run [mainSoy](Matlab/CortazarSchwartz/mainSalmon.m) to calibrate the soy model.
- Kalman Filter:
    1. Run [mainSalmon](Matlab/Kalman/mainSalmon.m) to calibrate the salmon model.
    2. Run [mainSoy](Matlab/Kalman/mainSalmon.m) to calibrate the soy model.

Optimal Stopping:
1. Run [mainLSMC](Matlab/RealOption/mainLSMC.m) with your custom model parameters to find the optimal harvesting time. Also returns a pathwise comparison of stopping rules with deterministic and stochastic feeding costs.

## Algorithms
- Commodity Calibration
    1. Kalman filter [Schwartz (1997)](https://doi.org/10.1111/j.1540-6261.1997.tb02721.x) 
    2. Nested Least Squares [Cortazar and Schwartz (2003)](https://doi.org/10.1016/S0140-9883(02)00096-8) 
- Optimal Stopping
    1. Least Square Monte Carlo (LSMC) [Longstaff and Schwartz (2001)](https://doi.org/10.1093/rfs/14.1.113) 
    2. Deep Optimal Stopping [Becker et al. (2021)](https://doi.org/10.1017/S0956792521000073) 
- Decision Boundary
    1. Deep Classifier [DeepDecision](Python/DeepDecision.py)

# Data
The data provided in this repository may not be used for anything else, than testing the algorithms in this repository. The license for the code does not extend to the datasets.

The most recent data concerning salmon futures can be downloaded from [fishpool.eu](https://fishpool.eu/forward-price-history/) while agreeing to their terms of usage.

For historical soybean futures we only provide a preprocessed sample based on market data to test the code. Results may differ from the values presented in the paper.