#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  22 11:57:02 2022

@author: Aidan Pavao for Massachusetts Host-Microbiome Center
 - Logistic function and variations
 - Use python 3.8+ for best results
 - See /nmr-cdiff/venv/requirements.txt for dependencies

Copyright 2022 Massachusetts Host-Microbiome Center

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

import numpy as np

def logi(x, L, k, x0):
    """Evaluate the logistic function at point x."""
    return L*np.reciprocal(1. + np.exp(-1*k*(x - x0)))

def logi_c(x, L, k, x0, C):
    """Evaluate the logistic function at point x with vertical shift."""
    return logi(x, L, k, x0) + C

def logi_sum(x, L1, k1, x01, L2, k2, x02):
    """Evaluate the sum of two logistic functions at point x."""
    return logi(x, L1, k1, x01) + logi(x, L2, k2, x02)

def logi_flex(x, *argv):
    """Evaluate the (sum of) logistic function at point x."""
    L = argv[0]
    k = argv[1]
    x0 = argv[2]
    sol = logi(x, L, k, x0)
    if len(argv) == 3:
        return sol # typical solution
    elif len(argv) == 4:
        return sol + argv[3] # + C
    elif len(argv) == 5:
        return sol - logi_flex(x, L, *argv[3:]) # L2 = -L1
    else:
        return sol + logi_flex(x, *argv[3:])  # multiple curves

def ddx_logi(x, L, k, x0):
    """Evaluate the derivative of the logistic function at point x."""
    return k*L*np.exp(-k*(x - x0))/(1. + np.exp(-1*k*(x - x0)))**2

def ddx_logi_flex(x, *argv):
    """Evaluate the derivative of the (sum of) logistic function at point x."""
    L = argv[0]
    k = argv[1]
    x0 = argv[2]
    sol = ddx_logi(x, L, k, x0)
    if len(argv) <= 4:
        return sol  # typical solution
    elif len(argv) == 5:
        return sol - ddx_logi_flex(x, L, *argv[3:]) # L2 = -L1
    else:
        return sol + ddx_logi_flex(x, *argv[3:])  # multiple curves

def inv_logi_c(x, L, k, x0, C):
    """Evaluate the inverse of the logistic function at point x."""
    return -1/k*np.log(L/(x - C) - 1) + x0

def dydL(x, L, k, x0, C=0):
    """Logistic function partial derivative with respect to L."""
    return 1/(1 + np.exp(-k*(x - x0)))

def dydk(x, L, k, x0, C=0):
    """Logistic function partial derivative with respect to k."""
    return L*(x - x0)*np.exp(-k*(x - x0))/(1. + np.exp(-k*(x - x0)))**2

def dydx0(x, L, k, x0, C=0):
    """Logistic function partial derivative with respect to x0."""
    return -L*k*np.exp(-k*(x - x0))/(1. + np.exp(-k*(x - x0)))**2

def dydC(x, L, k, x0, C=0):
    """Logistic function partial derivative with respect to C."""
    return 1

def err_logi(x, L, k, x0, C, del_L, del_k, del_x0, del_C):
    """Evaluate standard error of logistic function at point x.

    Parameters:
    L, k, x0, C -- logistic function coefficients
    del_L, del_k, del_x0, del_C -- standard error in logistic coefficients
    """
    coeffs = [L, k, x0, C]
    L_term = (dydL(x, *coeffs)*del_L)**2
    k_term = (dydk(x, *coeffs)*del_k)**2
    x0_term = (dydx0(x, *coeffs)*del_x0)**2
    C_term = (dydC(x, *coeffs)*del_C)**2
    return np.sqrt(L_term + k_term + x0_term + C_term)

def err_logi_flex(x, *argv):
    """Evaluate standard error of logistic function at point x.
    Handles case where the coefficient C is not provided (assume 0).
    """
    if len(argv) == 6:
        return err_logi(x, *argv[:4], 0, *argv[4:], 0)
    elif len(argv) == 8:
        return err_logi(x, *argv)
    else:
        return None

def dzdL(x, L, k, x0):
    """Partial derivative of logistic functn. derivative with respect to L."""
    return ddx_logi(x, L, k, x0)/L

def dzdk(x, L, k, x0):
    """Partial derivative of logistic functn. derivative with respect to k."""
    da = 2*(x - x0)*np.exp(-k*(x - x0))/(1 + np.exp(-k*(x - x0)))
    db = x0 - x
    dc = 1/k
    return ddx_logi(x, L, k, x0)*(da + db + dc)

def dzdx0(x, L, k, x0):
    """Partial derivative of logistic functn. derivative with respect to x0."""
    da = 2*k*np.exp(-k*(x - x0))/(1 + np.exp(-k*(x - x0)))
    db = -k
    return ddx_logi(x, L, k, x0)*(da + db)

def err_ddx_logi(x, L, k, x0, C, del_L, del_k, del_x0, del_C):
    """Evaluate standard error of logistic function derivative at point x.
    Ignores coefficient C.

    Parameters:
    L, k, x0, C -- logistic function coefficients
    del_L, del_k, del_x0, del_C -- standard error in logistic coefficients
    """
    coeffs = [L, k, x0]
    L_term = (dzdL(x, *coeffs)*del_L)**2
    k_term = (dzdk(x, *coeffs)*del_k)**2
    x0_term = (dzdx0(x, *coeffs)*del_x0)**2
    return np.sqrt(L_term + k_term + x0_term)

def err_ddx_logi_flex(x, *argv):
    """Evaluate standard error of logistic function derivative at point x.
    Handles case where the coefficient C is not provided (assume 0).
    """
    if len(argv) == 6:
        return err_ddx_logi(x, *argv[:4], 0, *argv[4:], 0)
    elif len(argv) == 8:
        return err_ddx_logi(x, *argv)
    else:
        return None
