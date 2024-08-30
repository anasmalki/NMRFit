import numpy as np
from tqdm import tqdm
import pybroom as br
from .params import make_params, make_params_random
from .data_format import make_inputs, compile_results, extract_rate_pred, make_inputs_MC, generate_MC
from .fit import fit_MF
from .nmr import relax


def fit_1_residue(relaxation_input, fields, residue, model, tslow, tfast, theta, method, nb_iter):
    a, b = np.array(list(make_inputs(relaxation_input, residue))).reshape(len(fields),len(fields)*4)
    paramsMF = make_params(model=model, tslow=tslow, tfast=tfast, theta=theta)
    m = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
    for n in range(nb_iter):
        paramsMF = make_params_random(model=model, tslow=tslow, tfast=tfast, theta=theta)
        m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25) 
        if m2.chisqr < m.chisqr:
            m = m2              
    return m

def batch_fit(relaxation_input, fields, models, tslow, tfast, theta, method, nb_iter):
    fit_stat = []
    fit_result = []
    for k in tqdm(relaxation_input[0]["num"]): 
        a, b = np.array(list(make_inputs(relaxation_input, k))).reshape(len(fields),len(fields)*4)
        for i in models:
            paramsMF = make_params(model=i, tslow=tslow, tfast=tfast, theta=theta)
            m = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
            for n in range(nb_iter):
                paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25) 
                if m2.chisqr < m.chisqr:
                    if i in [11,12]:
                        if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                           m = m2
                        else:
                            count = 0
                            while (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) !=1 and count < 50:
                                paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                                m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
                                count +=1
                            if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                               if m2.chisqr < m.chisqr:
                                    m = m2 
                    else:
                        m = m2                                                  
                if m2.chisqr > m.chisqr:
                    if i in [11,12]:
                        if (m.params["S2s"].value + m.params["S2f"].value + m.params["S3"].value) != 1:
                            count = 0
                            while (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) !=1 and count < 50:
                                paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                                m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
                                count +=1
                            if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                                m = m2            
            fit_stat.append(br.glance(m))
            fit_result.append(br.tidy(m).set_index('name').T)
    return fit_stat, fit_result

def fit_MC(relaxation_input, fields, compiled_results, models, tslow, tfast, theta, method, nb_iter, nb_MC): ### The MC normal distribution is generated around the predicted rates and not around the experimental rates, using the experimental error.
    pred = extract_rate_pred(compiled_results, fields)
    fit_stat = []
    fit_result =[]
    for k in tqdm(relaxation_input[0]["num"]):
        a, b = np.array(list(make_inputs(relaxation_input, k))).reshape(len(fields),len(fields)*4)
        a = np.array(list(make_inputs_MC(pred, k)))[0]
        nan_positions = np.isnan(b)
        a[nan_positions] = np.nan
        mc = generate_MC(a,b,nb_MC)
        bb =[]
        cc = []
        for i in models:
            paramsMF = make_params(model=i, tslow=tslow, tfast=tfast, theta=theta)
            ff = []
            rr =[]
            for v in range(len(mc)):
                a = mc[v]
                m = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
                for n in range(nb_iter):
                    paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                    m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25) 
                    if m2.chisqr < m.chisqr:
                        if i in [11,12]:
                            if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                                m = m2
                            else:
                                count = 0
                                while (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) !=1 and count < 50:
                                    paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                                    m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
                                    count +=1
                                if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                                    if m2.chisqr < m.chisqr:
                                        m = m2 
                        else:
                            m = m2                                                  
                    if m2.chisqr > m.chisqr:
                        if i in [11,12]:
                            if (m.params["S2s"].value + m.params["S2f"].value + m.params["S3"].value) != 1:
                                count = 0
                                while (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) !=1 and count < 50:
                                    paramsMF = make_params_random(model=i, tslow=tslow, tfast=tfast, theta=theta)
                                    m2 = fit_MF(data_x=a, data_y=b, fields=fields, par=paramsMF, method=method, nb_iter=25)
                                    count +=1
                                if (m2.params["S2s"].value + m2.params["S2f"].value + m2.params["S3"].value) == 1:
                                    m = m2                   
                ar = np.array([m.params["tm"].value, m.params["ts"].value, m.params["tf"].value, m.params["S2s"].value, m.params["S2f"].value, m.params["S3"].value, m.params["Rex"].value, m.params["theta"].value])
                br = np.array(list(relax(m.params, fields))).ravel().reshape(len(fields)*4)
                ff.append(ar) # shape(mc,params)
                rr.append(br)
            bb.append(ff) # shape(nbmodel, mc, params)
            cc.append(rr)
        fit_stat.append(bb)  # shape(res, nbmodel, mc, params)
        fit_result.append(cc)
    fit_stat = np.array(fit_stat).ravel().reshape(len(relaxation_input[0]["num"]), len(models), nb_MC, 8) ## 7 is the number of params: the 3 taus, the 3 S², Rex and theta
    return fit_stat, fit_result
