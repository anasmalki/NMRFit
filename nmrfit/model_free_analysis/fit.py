import numpy as np
import pandas as pd
import lmfit
from .nmr import relax

### Fit functions
def minfunction_MF(par, rates, err, fields):
        pred = np.array(list(relax(par, fields))).flatten()
        df = pd.DataFrame([pred, rates, err]).dropna(axis=1)
        pred = np.array(df.iloc[0])
        data = np.array(df.iloc[1])
        stderr = np.array(df.iloc[2])

        # return ((data - pred)/stderr)**2
        return ((data - pred)/stderr)


def fit_MF(data_x, data_y, fields, par, method, nb_iter, func2min=minfunction_MF): #data_x are the rates, data_y the errors
    minner= lmfit.Minimizer(func2min, par, fcn_args=(data_x,data_y, fields)) 
    result = minner.minimize(method=method)
    for i in range(nb_iter):
        minner= lmfit.Minimizer(func2min, result.params, fcn_args=(data_x,data_y, fields))
        r = minner.minimize(method=method)
        if r.redchi < result.redchi:
            result = r
    return result

