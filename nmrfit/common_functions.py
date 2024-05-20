import numpy as np
import pandas as pd
from scipy.constants import h, mu_0

g1H = 2.6752219 * 1e8 # rad s-1 T-1
g15N = -2.7126 * 1e7 # rad s-1 T-1
dCSA = -172.0 * 1e-6 #CSA, from Kroenke et al JACS 1999
rNH = 1.023 * 1e-10 #Average NH bond length, from Yao et al JACS 2008, in m
d = mu_0*h*g1H*g15N/(8*(np.pi**2)*(rNH**3)) #dipolar coupling

def P2(x):
    return (3*x*x - 1)/2

def MHz2T(field_MHz):
    return field_MHz*2*np.pi*1e6/g1H# MHz to Tesla

def import_data(path, prot, b0_n, conc_n, temp_n):
    df = pd.read_csv(f"{path}\input_{prot}_{b0_n}_{conc_n}_{temp_n}", header=0, delim_whitespace = True, encoding="UTF-8")
    df["T1"] = 1/df["R1"]
    df["T2"] = 1/df["R2"]
    df["R2/R1"] = df["R2"]/df["R1"]
    df["E_R2/R1"] = df["E_R2"]+df["E_R1"]
    return df

# def clean_data(datalong, datashort):
#     diff = datalong[~datalong["num"].isin(datashort["num"])].index
#     df = datalong.drop(axis=0, index=diff)
#     df = df.reset_index(drop=True)
#     return df

def clean_data(df1, df2, common_column):

    common_values = set(df1[common_column]).intersection(set(df2[common_column]))

    filtered_df1 = df1[df1[common_column].isin(common_values)]
    filtered_df2 = df2[df2[common_column].isin(common_values)]

    return filtered_df1, filtered_df2