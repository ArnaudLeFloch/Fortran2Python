# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:44:18 2025

@author: lefl_ar
"""

import numpy as np
from scipy.linalg import toeplitz

def optimalfilter(p_main, p_ref, L):
    # number of samples
    Ne = len(p_main)
    
    # cross-correlation vector (g) of main and reference signal
    g_cc = np.correlate(p_main, p_ref, mode='full')
    
    # cutting off cross-correlation at L-1
    g = g_cc[Ne-1:Ne+L-1]
    
    # autocorrelation vector (r) of reference signal
    r_ac = np.correlate(p_ref, p_ref, mode='full')
    
    # cutting off autocorrelation at L-1
    r = r_ac[Ne-1:Ne+L-1]
    
    # creating Toeplitz matrix
    # RR = np.linalg.toeplitz(r)
    RR = toeplitz(r)
    
    # calculation of filter coefficients f (optimization)
    f = np.linalg.inv(RR).dot(g)
    
    # FIR filtering of the noise (filter function)
    p_noise = np.zeros(Ne)
    p_noise[:L] = 0  # first part of the signal needed for prediction
    for n in range(L+1, Ne+1):
        sum_val = 0
        for i in range(L):
            sum_val += f[i] * p_ref[n-i-1]
        p_noise[n-1] = sum_val

    # subtraction of the noise from the signal
    p_filtered = p_main - p_noise

    # cutting off the first unfiltered part (set to zero)
    p_filtered[:L] = 0

    # # adding the mean of the original signal
    # p_filtered += np.mean(p_main)
    return p_filtered, f
    # return p_filtered, p_noise



# Key Changes and Notes:
# Cross-correlation and autocorrelation:
# MATLAB uses xcorr() and toeplitz(). In Python, we use np.correlate() for cross-correlation and autocorrelation. The mode='full' in np.correlate() corresponds to MATLAB's default behavior of xcorr().
# Toeplitz Matrix:
# The toeplitz() function from numpy.linalg is used to construct the matrix in Python.
# Matrix Inversion:
# np.linalg.inv() is used for matrix inversion, similar to MATLAB's inverse operator.
# Array Indexing:
# Python arrays are 0-indexed, so adjustments are made when iterating through arrays (n-i-1 instead of n-i in the filtering loop).
# Noise and Filtering:
# The filtering loop accumulates the prediction of the noise. The p_noise array is initialized with zeros, and the filtered signal is calculated by subtracting the noise from the main signal.
# Mean Adjustment:
# At the end, the mean of the original signal is added back to the filtered signal, as in the original MATLAB code.
# This Python version should behave the same as your MATLAB code, provided the input signals p_main and p_ref are NumPy arrays. Let me know if you need further adjustments!