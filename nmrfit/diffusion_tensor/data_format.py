import numpy as np
import pandas as pd

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
                    rr.append(np.nan)
                    ee.append(np.nan)
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
                rr.append(np.nan)
                ee.append(np.nan)
            if "etaXY" in data[i]:
                rr.append(data[i][data[i]["num"]==number]["etaXY"])
                ee.append(data[i][data[i]["num"]==number]["E_etaXY"])
            else:
                rr.append([np.nan])
                ee.append([np.nan])
    yield np.array([np.concatenate(rr).ravel(), np.concatenate(ee).ravel()])
    
def make_inputs_difftens(data):
    rr = pd.concat(data)[["R1", "R2"]].to_numpy()
    ee = pd.concat(data)[["E_R1", "E_R2"]].to_numpy()
    return rr, ee

def make_inputs_difftens2(data):
    rr = pd.concat(data)[["R1", "R2", "NOE"]].to_numpy()
    ee = pd.concat(data)[["E_R1", "E_R2", "E_NOE"]].to_numpy()
    
    return rr, ee

# def get_data_from_multifit(fit, preinput, NH_bond_vectors, fields):
#     res_fit = []
#     if fit.params["model"] in [10,11,12]:
#      for i in range(len(preinput)):
#          num = preinput["num"].iloc[i]
#          r1r2 = np.array(list(relax2(fit.params, i, a1NH_bond_vectors, fields))).ravel()
#         #  r2r1 = r1r2[1]/r1r2[0]
#          S2 = fit.params[f"S2_{i}"].value
#          ti = fit.params[f"ti_{i}"].value
#          rex = fit.params[f"Rex_{i}"].value
#          res_fit.append((num, r1r2[0], r1r2[1], S2, ti, rex))
#     df = pd.DataFrame(res_fit, columns=["num", "R1", "R2", "S2", "ti", "Rex"])
#     return df