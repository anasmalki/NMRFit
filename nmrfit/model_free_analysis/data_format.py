import numpy as np
import pandas as pd
from .nmr import relax, J
from ..common_functions import g15N, g1H, MHz2T

### Input functions
def complete_input_J(input):
    while len(input) < 3:
        input.append(0)
    return input


def prepare_input(data, fields, number):
    dd = []
    for i in range(len(data)):
        df = data[i].copy()
        df.insert(loc=0, column="field", value=fields[i])
        df = df.drop(columns=["T1", "T2", "R2/R1", "E_R2/R1"])
        dd.append(df)
    if number == 0:
        return pd.concat(dd)
    if number != 0:
        df = pd.concat(dd)
        return df[df["num"]==number]

def make_inputs(data, number):
    rr = []
    ee = []
    if number == 0:
        for i in range(len(data)):
                rr.append(data[i]["R1"])
                rr.append(data[i]["R2"])
                ee.append(data[i]["E_R1"])
                ee.append(data[i]["E_R2"])
                if "NOE" in data[i]:
                    rr.append(data[i]["NOE"])
                    ee.append(data[i]["E_NOE"])
                else:
                    rr.append([np.nan])
                    ee.append([np.nan])
                if "etaXY" in data[i]:
                    rr.append(data[i]["etaXY"])
                    ee.append(data[i]["E_etaXY"])
                else:
                    rr.append([np.nan])
                    ee.append([np.nan])
    if number != 0:
        for i in range(len(data)):
            rr.append(data[i][data[i]["num"]==number]["R1"])
            rr.append(data[i][data[i]["num"]==number]["R2"])
            ee.append(data[i][data[i]["num"]==number]["E_R1"])
            ee.append(data[i][data[i]["num"]==number]["E_R2"])
            if "NOE" in data[i]:
                rr.append(data[i][data[i]["num"]==number]["NOE"])
                ee.append(data[i][data[i]["num"]==number]["E_NOE"])
            else:
                rr.append([np.nan])
                ee.append([np.nan])
            if "etaXY" in data[i]:
                rr.append(data[i][data[i]["num"]==number]["etaXY"])
                ee.append(data[i][data[i]["num"]==number]["E_etaXY"])
            else:
                rr.append([np.nan])
                ee.append([np.nan])
    yield np.array([np.concatenate(rr).ravel(), np.concatenate(ee).ravel()])
    
def generate_MC(rates, error, nmc):
    rr = []
    for k in range(nmc):
        rr.append(np.random.normal(rates, error))
    return rr

def tidy_result_error(input_data, fit_out_error, nmodel):
    mc = []
    for k in range(nmodel):
        ff = []
        for i in range(len(input_data)):
            ff.append(pd.DataFrame(fit_out_error[i][k]).std())
        ff = pd.DataFrame(ff)
        ff.columns = ["E_tm", "E_ts", "E_tf", "E_S2s", "E_S2f", "E_S3", "E_Rex", "E_theta"]
        ff.insert(0, "num", input_data["num"].values)
        mc.append(pd.DataFrame(ff))
    return mc #shape (nmodel, nres, params)

### Function to clean fit data 
def tidy_result_stat(nmodel, input_data, fit_out_stat):
    rc = []
    for k in range(nmodel):
        ff = []
        for i in range(len(input_data["num"])):
            a = fit_out_stat[k::nmodel][i]["redchi"][0]
            b = fit_out_stat[k::nmodel][i]["chisqr"][0]
            c = fit_out_stat[k::nmodel][i]["AIC"][0]
            d = fit_out_stat[k::nmodel][i]["BIC"][0]
            # rr.append()
            ff.append([a, b, c, d])
        ff = pd.DataFrame(ff)
        ff.columns = ["redchi", "chisqr", "AIC", "BIC"]
        ff.insert(0, "num", input_data["num"].values)
        rc.append(ff)
    return rc

### Functions for rates outputs
def prep_result_rates(nmodel, input_data, fit_result, fields):
    rc = []
    for k in range(nmodel):
        rr = []
        for i in range(len(input_data["num"])):
            a = pd.DataFrame(list(relax(fit_result[k::nmodel][i], fields)))
            a.columns = ["R1", "R2", "NOE", "etaXY"]
            rr.append(a)
        rc.append(rr)
    # df = pd.DataFrame(rc).T
    # df.columns = ["model_1", "model_2", "model_3", "model_4", "model_5", "model_6", "model_7", "model_8", "model_10", "model_11", "model_12", "model_13"]
    # df.insert(0, "num", input_data["num"])
    return rc



def tidy_result_rates(nmodel, input_data, res_rate_prep, fields):
    rc = []
    for k in range(nmodel):
        ff = []
        for n in range(len(fields)):
            r1 = [res_rate_prep[k][i].iloc[n][0] for i in range(len(res_rate_prep[0]))]
            r2 = [res_rate_prep[k][i].iloc[n][1] for i in range(len(res_rate_prep[0]))]
            noe = [res_rate_prep[k][i].iloc[n][2] for i in range(len(res_rate_prep[0]))]
            etaxy = [res_rate_prep[k][i].iloc[n][3] for i in range(len(res_rate_prep[0]))]
            df=pd.DataFrame([r1, r2, noe, etaxy]).T
            f = fields[n]
            df.columns = [f"R1_pred_{f}", f"R2_pred_{f}", f"NOE_pred_{f}", f"etaXY_pred_{f}"]
            df.insert(0, "num", input_data["num"].values)
            ff.append(df)
        rc.append(ff) # (#model, field, rates)
    return rc

def tidy_result_rates_error(nmodel, fit_rates_error, num_res):
    rc = []
    for k in range(nmodel):
        ff = []
        for i in range(len(fit_rates_error)):
            e = pd.DataFrame(fit_rates_error[i][k])
            ff.append(e.std())
        f = pd.DataFrame(ff)
        f.insert(0, "num", num_res)
        rc.append(f)
    return rc

### multi_index ou non ?
# res_stat = tidy_result_redchi(nmodel, hd_600_350_293, fit_stat)
# d = {'model_1':res_stat[0], 'model_2':res_stat[1]}
#
# pd.concat(d, axis=1, keys=d.keys())["model_1"]


# res_rchi = tidy_result_redchi(nmodel, hd_600_350_293)

## Function to tidy results Params and Rex 
def tidy_result_rex(nmodel, input_data, fit_result):
    rc = []
    for k in range(nmodel):
        rr = []
        for i in range(len(input_data["num"])):
            a = fit_result[k::nmodel][i]["Rex"][0]
            rr.append(a)
        rc.append(rr)
    # df = pd.DataFrame(rc).T
    # df.columns = ["model_1", "model_2", "model_3", "model_4", "model_5", "model_6", "model_7", "model_8", "model_10", "model_11", "model_12", "model_13"]
    # df.insert(0, "num", input_data["num"])
    return rc



def tidy_result_params(nmodel, input_data, fit_result):
    rc = []
    for k in range(nmodel):
        rr = []
        for i in range(len(input_data["num"])):
            tm = fit_result[k::nmodel][i]["tm"][0]
            ts = fit_result[k::nmodel][i]["ts"][0]
            tf = fit_result[k::nmodel][i]["tf"][0]
            S2s = fit_result[k::nmodel][i]["S2s"][0]
            S2f = fit_result[k::nmodel][i]["S2f"][0]
            S3 = fit_result[k::nmodel][i]["S3"][0]
            Rex = fit_result[k::nmodel][i]["Rex"][0]
            theta = fit_result[k::nmodel][i]["theta"][0]
            model = fit_result[k::nmodel][i]["model"][0]
            rr.append([tm, ts, tf, S2s, S2f, S3, Rex, theta, model])
        df = pd.DataFrame(rr)
        df.columns = ["tm", "ts", "tf", "S2s", "S2f", "S3", "Rex", "theta", "model"]
        df.insert(loc=0, column="num", value=input_data["num"].values)
        for col in df.columns:
            if len(df[col].unique()) == 1:
                df.drop(col,inplace=True,axis=1)
        rc.append(df)
    return rc

def calculate_rate_error(input_data, rate_error, fit_rate_error, fields):
    rc = []
    for k in range(len(rate_error[0])): #model
        bb = []
        for i in range(len(rate_error)): #res           
            n = pd.DataFrame(fit_rate_error[i][k]).std()
            bb.append(n)
            ff = []
            for n in range(len(fields)):
                f = fields[n]
                ff.append((f"R1_pred_{f}", f"R2_pred_{f}", f"NOE_pred_{f}", f"etaXY_pred_{f}"))
            df = pd.DataFrame(bb)
            df.columns = [np.array(ff).ravel()]
            df.insert(0, 'num', input_data['num'])
        rc.append(df)
    return rc

def compile_results(nmodel, relax_data, fit_stat, fit_result, fields):
    res_stat = tidy_result_stat(nmodel, relax_data, fit_stat)
    res_rate_prep = prep_result_rates(nmodel, relax_data, fit_result, fields)
    res_rates = tidy_result_rates(nmodel, relax_data, res_rate_prep, fields)
    res_rex = tidy_result_rex(nmodel, relax_data, fit_result)
    res_params = tidy_result_params(nmodel, relax_data, fit_result)
    if nmodel > 1:
        dd = []
        for n in range(nmodel):
            df = res_params[n].copy()
            df["Rex"] = res_rex[n]
            for k in range(len(fields)):
                df = df.merge(res_rates[n][k], on='num', how='left')
            df = df.merge(res_stat[n], on='num', how='left')
            dd.append(df)
        return dd
    if nmodel == 1:
        df = res_params[0].copy()
        df["Rex"] = res_rex[0]
        for k in range(len(fields)):
            df = df.merge(res_rates[0][k], on='num', how='left')
        df = df.merge(res_stat[0], on='num', how='left')
        return df
    
### SDF functions


def calculate_SDF(fit_result, fields, model):
    ff = []
    if 'S3' not in fit_result.columns:
        fit_result['S3'] = 0
    if 'ts' not in fit_result.columns:
        fit_result['ts'] = 0
    for i in range(len(fit_result)):
        jj = []
        data = fit_result.iloc[i]
        taus = [data['tm'], data['tm'], data['tf']]
        Ss = [data['S2s'], data['S2f'], data['S3']]
        j0 = J(0, taus, Ss, model)
        jj.append(j0)
        for f in fields:
            wn = -g15N*(MHz2T(f))
            wh = -g1H*(MHz2T(f))
            jwn = J(wn, taus, Ss, model)
            jwh = J(wh, taus, Ss, model)
            j0870wh = J(0.870*wh, taus, Ss, model)
            j0921wh = J(0.921*wh, taus, Ss, model)
            j0955wh = J(0.955*wh, taus, Ss, model)                   
            jj = jj + [jwn, jwh, j0870wh, j0921wh, j0955wh]
        ff.append(jj)
    ff = pd.DataFrame(ff)
    col = ['J0']
    for k in range(len(fields)):
        f = fields[k]
        col = col + [f'JwN_{f}', f'JwH_{f}', f'J0870wH_{f}', f'J0921wH_{f}', f'J0955wH_{f}']
    ff.columns = col
    ff.insert(0, "num", fit_result["num"].values)
    return ff

def calculate_SDF_error(input_data, fit_resultMC_error, fields, model):
    ff = []
    for i in range(len(input_data)):
        jj = []
        for m in range(len(fit_resultMC_error[0][0])):
            bb = []
            data = fit_resultMC_error[i][0][m]
            tm, ts, tf, s2s, s2f, s3, rex, theta = data
            taus = [tm, ts, tf]
            Ss = [s2s, s2f, s3]
            j0 = J(0, taus, Ss, model)
            bb.append(j0)
            for f in fields:
                wn = -g15N*(MHz2T(f))
                wh = -g1H*(MHz2T(f))
                jwn = J(wn, taus, Ss, model)
                jwh = J(wh, taus, Ss, model)
                j0870wh = J(0.870*wh, taus, Ss, model)
                j0921wh = J(0.921*wh, taus, Ss, model)
                j0955wh = J(0.955*wh, taus, Ss, model)                      
                bb = bb + [jwn, jwh, j0870wh, j0921wh, j0955wh]
            jj.append(bb)
        ff.append(pd.DataFrame(jj).std())
    ff = pd.DataFrame(ff)
    col = ['E_J0']
    for k in range(len(fields)):
        f = fields[k]
        col = col + [f'E_JwN_{f}', f'E_JwH_{f}', f'E_J0870wH_{f}', f'E_J0921wH_{f}', f'E_J0955wH_{f}']
    ff.columns = col
    ff.insert(0, "num", input_data["num"].values)
    return ff

def add_error_to_results(fit_result, param_error, rate_error, fields):
    if 'tm' in fit_result.columns:
        idx = fit_result.columns.get_loc('tm')
        fit_result.insert(idx+1, 'E_tm' ,param_error['E_tm'])
    if 'ts' in fit_result.columns:
        idx = fit_result.columns.get_loc('ts')
        fit_result.insert(idx+1, 'E_ts' ,param_error['E_ts'])
    if 'tf' in fit_result.columns:        
        idx = fit_result.columns.get_loc('tf')
        fit_result.insert(idx+1, 'E_tf' ,param_error['E_tf'])
    if 'S2s' in fit_result.columns:
        idx = fit_result.columns.get_loc('S2s')
        fit_result.insert(idx+1, 'E_S2s' ,param_error['E_S2s'])       
    if 'S2f' in fit_result.columns:
        idx = fit_result.columns.get_loc('S2f')
        fit_result.insert(idx+1, 'E_S2f' ,param_error['E_S2f'])         
    if 'S3' in fit_result.columns:
        idx = fit_result.columns.get_loc('S3')
        fit_result.insert(idx+1, 'E_S3' ,param_error['E_S3'])           
    if 'Rex' in fit_result.columns:
        idx = fit_result.columns.get_loc('Rex')
        fit_result.insert(idx+1, 'E_Rex' ,param_error['E_Rex'])           
    if 'theta' in fit_result.columns:    
        idx = fit_result.columns.get_loc('theta')
        fit_result.insert(idx+1, 'E_theta' ,param_error['E_theta'])
        
    for n in range(len(fields)):
        f = fields[n]
        g = f'pred_{f}'
        if f'R1_{g}' in fit_result.columns:
            idx = fit_result.columns.get_loc(f'R1_{g}')
            fit_result.insert(idx+1, f'E_R1_{g}', rate_error[f'R1_{g}'])          
        if f'R2_{g}' in fit_result.columns:
            idx = fit_result.columns.get_loc(f'R2_{g}')
            fit_result.insert(idx+1, f'E_R2_{g}', rate_error[f'R2_{g}'])               
        if f'NOE_{g}' in fit_result.columns:
            idx = fit_result.columns.get_loc(f'NOE_{g}')
            fit_result.insert(idx+1, f'E_NOE_{g}', rate_error[f'NOE_{g}'])   
        if f'etaXY_{g}' in fit_result.columns:
            idx = fit_result.columns.get_loc(f'etaXY_{g}')
            fit_result.insert(idx+1, f'E_etaXY_{g}', rate_error[f'etaXY_{g}'])   
                                    
    return fit_result
