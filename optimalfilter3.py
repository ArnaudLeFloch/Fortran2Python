# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 17:43:34 2025

@author: lefl_ar
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 12:17:15 2025

@author: lefl_ar
"""

import numpy as np
from scipy.linalg import toeplitz
from concurrent.futures import ThreadPoolExecutor
from scipy import signal

def optimalfilter3(p_main, p_ref, L):
    # number of samples
    Ne = len(p_main)

    # cross-correlation vector (g) of main and reference signal
    g_cc = signal.fftconvolve(p_main, p_ref[::-1], mode='full')  
    # g_cc =     np.correlate(p_main, p_ref, mode='full') # longuet


    # cutting off cross-correlation at L-1
    g = g_cc[Ne-1:Ne+L-1]

    # autocorrelation vector (r) of reference signal
    r_ac = signal.fftconvolve(p_ref, p_ref[::-1], mode='full')  
    # r_ac = np.correlate(p_ref, p_ref, mode='full') # longuet aussi

    # cutting off autocorrelation at L-1
    r = r_ac[Ne-1:Ne+L-1]

    # creating Toeplitz matrix from the autocorrelation vector
    RR = toeplitz(r)

    # calculation of filter coefficients f (Wiener filter optimization)
    f = np.linalg.inv(RR).dot(g)  
    # f01 = np.linalg.inv(RR).dot(g)  
    # f02 =np.dot (np.linalg.inv(RR),g)
    
    
    # FIR filtering of the noise (filter function)
    p_noise = np.zeros(Ne)
    p_noise[:L] = 0 # first part of the signal needed for prediction

    # Function to compute p_noise[n] for each n (parallelizable part)
    def compute_p_noise(n):
        sum_val = 0
        for i in range(L):
            sum_val += f[i] * p_ref[n - i - 1]
        return n, sum_val

    # Using ThreadPoolExecutor to parallelize the computation of p_noise
    with ThreadPoolExecutor(max_workers=16) as executor:
        results = executor.map(compute_p_noise, range(L+1, Ne+1), timeout=None, chunksize=1)  # LONG

    # Assigning results to p_noise array
    for n, sum_val in results:
        p_noise[n-1] = sum_val

    # Subtraction of the noise from the signal
    p_filtered = p_main - p_noise

    # Cutting off the first unfiltered part (set to zero)
    p_filtered[:L] = 0

    return p_filtered #, p_noise
