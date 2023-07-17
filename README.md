# AquacultureStochasticFeeding
This repository is complementing [TODO](https://www.arxiv.org/) for investigating the effect of stochastic feeding costs to optimal harvesting rules in aquaculture farms.

# Installation
There is no installation required, please download the code and run the main files. There is also a Jupyter notebook available for testing the Python code in colab.

[![Open In Colab TODO](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/googlecolab/colabtools/blob/master/notebooks/colab-github-demo.ipynb)

# Code
The code is structured as follows: 

For the data preprocessing and calibration of the commodity models to the market data we used Matlab with its (Global) Optimization Toolbox. 

We provide the code for estimating the optimal harvesting times with a Least-Square Monte-Carlo (LSMC) regression in both languages [Matlab](../blob/main/Matlab/) and [Python](../blob/main/Python/). 
Everything concerning Deep Learning is only implemented in Python.

## Python
Structure:

## Matlab
Structure:

## Customization

# Algorithms
- Commodity Calibration
    1. Kalman filter [TODO](https://www.arxiv.org/) 
    2. Nested Least Squares [TODO](https://www.arxiv.org/) 
- Optimal Stopping
    1. Least Square Monte Carlo (LSMC) [TODO](https://www.arxiv.org/) 
    2. Deep Optimal Stopping [TODO](https://www.arxiv.org/) 
- Decision Boundary
    1. Deep Classifier [TODO](https://www.arxiv.org/) 

# Data
The most recent data concerning salmon futures can be downloaded from [fishpool.eu](https://fishpool.eu/forward-price-history/). 

Historical soybean futures are not available publicly to the best of our knowledge and we only provide a preprocessed sample to test the code.  