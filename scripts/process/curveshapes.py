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

class LogisticSet:
    """Class to store coefficients for multiple logistic equations.
    Equations will be evaluated independently and solutions averaged (with
    error propagation).
    """
    def __init__(self):
        """Init set -- curves and errors are lists of coefficients/errs."""
        self._curves = []
        self._errors = []
        self._scalar = 1. # scale factor for logistic
        self._time = None # time scale to evaluate
        self._val = None # evaluated function over time series
        self._err = None # evaluated errors over time series
        self._dval = None # evaluated derivative over time series
        self._derr = None # evaluated derivative error over time series

    def add_curve(self, params, errors):
        self._curves.append(params.tolist())
        self._errors.append(errors.tolist())

    def _eval(self, x, eval_fn, err_fn, avg_first=True):
        """Internal evaluate with eval_fn and error function err_fn."""
        if avg_first:
            val, err = self._eval_avg_first(x, eval_fn, err_fn)
        else:
            val, err = self._eval_avg_last(x, eval_fn, err_fn)
        return self._scalar*val, err

    def _eval_avg_last(self, x, eval_fn, err_fn):
        """Internal evaluate with eval_fn and error function err_fn.
        Evaluates each curve in self._curves independently, and averages the
        result.
        """
        eval_sol = np.array([eval_fn(x, *ps) for ps in self._curves])
        eval_err = np.array([err_fn(x, *ps, *es) for ps, es in zip(self._curves, self._errors)])
        avg_sol = np.sum(eval_sol, axis=0)/eval_sol.shape[0]
        avg_err = np.sqrt(np.sum(np.square(eval_err), axis=0))/eval_err.shape[0]
        return avg_sol, avg_err

    def _eval_avg_first(self, x, eval_fn, err_fn):
        """Internal evaluate with eval_fn and error function err_fn.
        Averages each parameter accross all curves, then evaluates using the
        average coefficient values.
        """
        avg_ps = np.average(np.array(self._curves), axis=0)
        avg_es = rss(np.array(self._errors), axis=0)
        avg_sol = eval_fn(x, *avg_ps)
        avg_err = err_fn(x, *avg_ps, *avg_es)
        return avg_sol, avg_err

    def eval(self, ts):
        """Evaluate and store function and derivative over time series ts."""
        self._time = ts
        self._val, self._err = self._eval(ts, logi_flex, err_logi_flex)
        self._dval, self._derr = self._eval(ts, ddx_logi_flex, err_ddx_logi_flex)
        return self._val, self._err, self._dval, self._derr

    def get_sol(self, t):
        """Get solution at time t if stored, else calculate it."""
        idx = np.nonzero(self._time == t)
        if idx[0].size == 0:
            val, err = self._eval(t, logi_flex, err_logi_flex)
            dval, derr = self._eval(t, ddx_logi_flex, err_ddx_logi_flex)
            return val, err, dval, derr
        else:
            i = idx[0][0]
            return self._val[i], self._err[i], self._dval[i], self._derr[i]

    def get_bounds(self, t):
        """Get 95% confint bounds at time t if stored, else calculate them."""
        val, err, dval, derr = self.get_sol(t)
        exch_l = dval - 1.96*derr
        exch_u = dval + 1.96*derr
        signal_l = val - 1.96*err
        signal_u = val + 1.96*err
        return exch_l, exch_u, signal_l, signal_u

    def halfmax(self):
        """Return timepoint of maximum flux (logistic half-max)."""
        if sum([curve[1] for curve in self._curves]) < 0:
            return self._time[np.argmin(self._dval)]
        else:
            return self._time[np.argmax(self._dval)]

    def tshift(self, t):
        """Shift all curves' half-max (x0) to t."""
        for param_set in self._curves:
            param_set[2] = t
        self.eval(self._time)

    def set_runcount(self, runcount):
        """Calculate scalar from total run count."""
        self._scalar = len(self._curves)/runcount

    def display(self):
        print(self._curves)
        print(self._errors)

    def display_avg_coeffs(self, param_names, prefix="\t"):
        """Display average coefficients p/m error. Pass in coefficient names."""
        avg_ps = np.average(np.array(self._curves), axis=0)
        avg_es = rss(np.array(self._errors), axis=0)
        for i, par in enumerate(param_names):
            print(f"{prefix}{par}: {avg_ps[i]:.3f} ± {avg_es[i]:.3f}")

def rss(args, axis=None):
    """Calculate the square root of the sum of squares."""
    return np.sqrt(np.sum(np.array(args)**2, axis=axis))

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
