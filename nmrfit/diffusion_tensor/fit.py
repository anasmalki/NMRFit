import numpy as np
import pandas as pd
from .nmr import relax, relax2
import lmfit

### Fit functions difftens
def minfunction_difftens(par, rates, err, NH_bond_vectors, fields):
        pred = np.array(list(relax(par, NH_bond_vectors, fields)))
        pred = pred[:,1]/pred[:,0]
        # df = pd.DataFrame([pred, rates, err]).dropna(axis=1)
        # pred = np.array(df.iloc[0])
        data = np.array(rates)
        data = data[:,1]/data[:,0] #R2/R1
        
        stderr = np.array(err)
        stderr = stderr[:,1]+stderr[:,0]
        
        return ((data - pred)/stderr) 



def fit_difftens(data_x, data_y, NH_bond_vectors, fields, par, method, nb_iter, func2min=minfunction_difftens): #data_x are the rates, data_y the errors
    minner= lmfit.Minimizer(func2min, par, fcn_args=(data_x,data_y, NH_bond_vectors, fields), nan_policy='propagate', scale_covar=False) 
    result = minner.minimize(method=method)
    for i in range(nb_iter):
        minner= lmfit.Minimizer(func2min, result.params, fcn_args=(data_x,data_y, NH_bond_vectors, fields), nan_policy='propagate', scale_covar=False)
        r = minner.minimize(method=method)
        if r.redchi < result.redchi:
            result = r
    return result

### Intmol = internal mobility so pseudo modelfree with anisotropic diffusion
def minfunction_intmol(par, rates, err, NH_bond_vectors, fields):
        pred = np.array(list(relax2(par, NH_bond_vectors, fields))).flatten()
        df = pd.DataFrame([pred, rates, err]).dropna(axis=1)
        pred = np.array(df.iloc[0])
        data = np.array(df.iloc[1])
        stderr = np.array(df.iloc[2])
        
        return ((data - pred)/stderr)


def fit_intmol(data_x, data_y, NH_bond_vectors, fields, par, method, nb_iter, func2min=minfunction_intmol): #data_x are the rates, data_y the errors
    minner= lmfit.Minimizer(func2min, par, fcn_args=(data_x,data_y, NH_bond_vectors, fields)) 
    result = minner.minimize(method=method)
    for i in range(nb_iter):
        minner= lmfit.Minimizer(func2min, result.params, fcn_args=(data_x,data_y, NH_bond_vectors, fields))
        r = minner.minimize(method=method)
        if r.redchi < result.redchi:
            result = r
    return result