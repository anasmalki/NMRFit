import pandas as pd
import numpy as np
from ..common_functions import g1H, g15N, dCSA, d, P2, MHz2T

def J0870wh(R1, NOE):
    return (2*g15N*R1*(-1+NOE))/(g1H*(d**2))

def J_09XXwH(e, j0870wh):
    return ((0.870/e)**2)*j0870wh

def Jwn(R1, j0921wh, c):
    return (R1-((7*(d**2))/10)*j0921wh)/(3*(d**2/10)+(2*c**2)/15)

def J0(R2, jwn, j0955wh, c):
    return (R2 - ((3*(d**2)/20)+(3*(c**2)/45))*(jwn) - (13*((d**2)/20))*j0955wh)/((4*(d**2)/20)+4*((c**2)/45))

def J0_etaXY(etaXY, theta, jwn, c):
    return (etaXY - ((3*c*d*P2(np.cos(theta))*jwn*(1/15))))/(c*d*P2(np.cos(theta))*4*(1/15))


def calculate_SDM(relax_data, field, nb_MC, theta=None):
    
    wn = -g15N*MHz2T(field)
    c = dCSA*wn
    
    dd = []
    for i in relax_data['num']:
        data = relax_data[relax_data['num']==i]
        R1 = data['R1']
        NOE = data['NOE']
        R2 = data['R2']

        j0870wh = J0870wh(R1, NOE).values[0]
        j0921wh = J_09XXwH(0.921, j0870wh)
        j0955wh = J_09XXwH(0.955, j0870wh)
        jwn = Jwn(R1, j0921wh, c).values[0]
        j0 = J0(R2, jwn, j0955wh, c).values[0]
        
        R1_MC = np.random.normal(data["R1"],data["E_R1"], nb_MC)
        R2_MC = np.random.normal(data["R2"],data["E_R2"], nb_MC)
        NOE_MC = np.random.normal(data["NOE"],data["E_NOE"], nb_MC)
        
        e_j0870wh = np.std(J0870wh(R1_MC, NOE_MC), ddof=1)
        e_j0921wh = np.std(J_09XXwH(0.921, J0870wh(R1_MC, NOE_MC)), ddof=1)
        e_j0955wh = np.std(J_09XXwH(0.955, J0870wh(R1_MC, NOE_MC)), ddof=1)
        e_jwn = np.std(Jwn(R1_MC, J_09XXwH(0.921, J0870wh(R1_MC, NOE_MC)), c), ddof=1)
        e_j0 = np.std(J0(R2_MC, Jwn(R1_MC, J_09XXwH(0.921, J0870wh(R1_MC, NOE_MC)), c), J_09XXwH(0.955, J0870wh(R1_MC, NOE_MC)), c), ddof=1)
        
        if 'etaXY' in relax_data.columns:
            etaXY = data['etaXY']
            j0_etaxy = J0_etaXY(etaXY, theta, jwn, c).values[0]
            etaXY_MC = np.random.normal(data["etaXY"],data["E_etaXY"], nb_MC)
            e_j0_etaXY = np.std(J0_etaXY(etaXY_MC, theta, Jwn(R1_MC, J_09XXwH(0.921, J0870wh(R1_MC, NOE_MC)), c), c), ddof=1)
            
            dd.append([i, j0, e_j0, j0_etaxy, e_j0_etaXY, jwn, e_jwn, j0870wh, e_j0870wh, j0921wh, e_j0921wh, j0955wh, e_j0955wh])
            
        else:
            dd.append([i, j0, e_j0, jwn, e_jwn, j0870wh, e_j0870wh, j0921wh, e_j0921wh, j0955wh, e_j0955wh])
        
    df = pd.DataFrame(dd)
    if len(df.columns) == 13:
        df.columns = ['num' ,'J0', 'E_J0', 'J0_etaXY', 'E_J0_etaXY', 'JwN', 'E_JwN', 'J0870wH', 'E_J0870wH', 'J0921wH', 'E_J0921wH', 'J0955wH', 'E_J0955wH']
    if len(df.columns) == 11:
        df.columns = ['num' ,'J0', 'E_J0', 'JwN', 'E_JwN', 'J0870wH', 'E_J0870wH', 'J0921wH', 'E_J0921wH', 'J0955wH', 'E_J0955wH']
    
    return df