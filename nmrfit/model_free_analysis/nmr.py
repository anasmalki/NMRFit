import numpy as np
from ..common_functions import g1H, g15N, dCSA, d, P2, MHz2T

def J(w, taus, Ss, model):
    tm, ts, tf = taus
    S2s, S2f, S3 = Ss
    if model==1 or model==2:
        return (S2s)*((tm/(1+(w*tm)**2)))
    if model==3 or model==4:
        return (((S2s*tm)/(1+((w*tm)**2)))+(((1-S2s)*tf)/(1+((w*tf)**2))))
    # if model==5:
    #     ts = params["ts"].value
    #     S2f = params["S2f"].value
    #     tau = ((ts*tm)/(ts+tm))
    #     return S2f*(((S2s*S2f*tm)/(1+((w*tm)**2)))+(((1-(S2s*S2f))*tau)/(1+((w*tau)**2))))
    if model==5 or model==6: ## S² = S²f*S²s
        tau_s = ((ts*tm)/(ts+tm))
        tau_f = ((tf*tm)/(tf+tm))
        return (((S2s*S2f*tm)/(1+((w*tm)**2)))+(((1-S2s)*tau_s)/(1+((w*tau_s)**2)))+(((S2f-(S2s*S2f))*tau_f)/(1+((w*tau_f)**2))))
    if model==7 or model==8:
        tau = ((ts*tf)/(ts+tf))
        return (((S2s*S2f*tm)/(1+((w*tm)**2)))+(((S2f-(S2s*S2f))*tau)/(1+((w*tau)**2))))
    if model==9 or model==10:
        return (S2f*tf/(1+(w*tf)**2))+(S2s*tm/(1+(w*tm)**2))
    if model==11 or model==12:
        return (S3*ts/(1+(w*ts)**2))+(S2f*tf/(1+(w*tf)**2))+(S2s*tm/(1+(w*tm)**2))
    
### Relax function
def relax(par, fields):
    model = int(par["model"].value)
    taus = [par["tm"].value, par["ts"].value, par["tf"].value]
    Ss = [par["S2s"].value, par["S2f"].value, par["S3"].value]
    Rex = par["Rex"].value
    if par["theta"].value is None:
        theta == 22  # noqa: F821
    else:
        theta = par["theta"].value
    for i in fields:
        fieldT = MHz2T(i)
        wh = -g1H * fieldT
        wn = -g15N * fieldT
        csa= (wn*dCSA)
        # wh = params["wh"].value
        # wn = params["wn"].value

        R1 = ((d**2)/10)*(J(wh-wn, taus, Ss, model)+(3*J(wn, taus, Ss, model))+(6*J(wh+wn, taus, Ss, model)))+((2/15)*((csa**2))*J(wn, taus, Ss, model))
        R2 = (Rex*(i**2)) + (((d**2)/20)*((4*J(0, taus, Ss, model))+J(wh-wn, taus, Ss, model)+(3*J(wn, taus, Ss, model))+(6*J(wh, taus, Ss, model))+(6*J(wh+wn, taus, Ss, model)))+(((csa**2)/45)*((4*J(0, taus, Ss, model))+(3*J(wn, taus, Ss, model)))))
        NOE = 1.0 +(d**2)/(10*R1)*(g1H/g15N)*(6*J(wh+wn, taus, Ss, model)-(J(wh-wn, taus, Ss, model)))
        # sNH = (1/10)*(d**2)*((6*J(wh+wn, taus, Ss, model))-(J(wh-wn,taus)))
        # etaZ = (1/15)*csa*d*(P2(np.cos(theta)))*(6*J(wn, taus, Ss, model))
        etaXY = (1/15)*csa*d*(P2(np.cos(theta)))*(4*(J(0, taus, Ss, model)+3*J(wn, taus, Ss, model)))

        yield R1, R2, NOE, etaXY

def relax_MC(par, fields):
    model = par[0]
    taus = par[1]
    Ss = par[2]
    Rex = par[3]
    theta = par[4]
    for i in fields:
        fieldT = MHz2T(i)
        wh = -g1H * fieldT
        wn = -g15N * fieldT
        csa= (wn*dCSA)
        # wh = params["wh"].value
        # wn = params["wn"].value

        R1 = ((d**2)/10)*(J(wh-wn, taus, Ss, model)+(3*J(wn, taus, Ss, model))+(6*J(wh+wn, taus, Ss, model)))+((2/15)*((csa**2))*J(wn, taus, Ss, model))
        R2 = (Rex*(i**2)) + (((d**2)/20)*((4*J(0, taus, Ss, model))+J(wh-wn, taus, Ss, model)+(3*J(wn, taus, Ss, model))+(6*J(wh, taus, Ss, model))+(6*J(wh+wn, taus, Ss, model)))+(((csa**2)/45)*((4*J(0, taus, Ss, model))+(3*J(wn, taus, Ss, model)))))
        NOE = 1.0 +(d**2)/(10*R1)*(g1H/g15N)*(6*J(wh+wn, taus, Ss, model)-(J(wh-wn, taus, Ss, model)))
        # sNH = (1/10)*(d**2)*((6*J(wh+wn, taus, Ss, model))-(J(wh-wn,taus)))
        # etaZ = (1/15)*csa*d*(P2(np.cos(22)))*(6*J(wn, taus, Ss, model))
        etaXY = (1/15)*csa*d*(P2(np.cos(theta)))*(4*(J(0, taus, Ss, model)+3*J(wn, taus, Ss, model)))

        yield R1, R2, NOE, etaXY