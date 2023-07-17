# AquacultureStochasticFeeding
This repository is complementing [TODO](https://www.arxiv.org/) for investigating the effect of stochastic feeding costs to optimal harvesting rules in aquaculture farms.

# Installation
There is no installation required, please download the code and run the main files. For Python you need all the necessary dependencies to run Tensorflow 2.x. There is also a Jupyter notebook available for testing the Python code in Colab.

[![Open In Colab TODO](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/googlecolab/colabtools/blob/master/notebooks/colab-github-demo.ipynb)

# Code
The code is structured as follows: 

For the data preprocessing and calibration of the commodity models to the market data we used Matlab with its (Global) Optimization Toolbox. 

We provide the code for estimating the optimal harvesting times with a Least-Square Monte-Carlo (LSMC) regression in both languages [Matlab](Matlab/) and [Python](Python/main.py). 
Everything concerning Deep Learning is only implemented in Python.

## Python
Structure:

## Matlab
Structure:

## Customization

# Algorithms
- Commodity Calibration
    1. Kalman filter [Schwartz (1997)](https://doi.org/10.1111/j.1540-6261.1997.tb02721.x) 
    2. Nested Least Squares [Cortazar and Schwartz (2003)](https://doi.org/10.1016/S0140-9883(02)00096-8) 
- Optimal Stopping
    1. Least Square Monte Carlo (LSMC) [Longstaff and Schwartz (2001)](https://doi.org/10.1093/rfs/14.1.113) 
    2. Deep Optimal Stopping [Becker et al. (2021)](https://doi.org/10.1017/S0956792521000073) 
- Decision Boundary
    1. Deep Classifier [tensorflow.org](https://www.tensorflow.org/tutorials/keras/classification) 

# Data
The most recent data concerning salmon futures can be downloaded from [fishpool.eu](https://fishpool.eu/forward-price-history/). 

Historical soybean futures are not available publicly to the best of our knowledge and we only provide a preprocessed sample to test the code.  