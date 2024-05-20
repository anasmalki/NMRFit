import numpy as np
from ..common_functions import g1H, g15N, dCSA, d, P2, MHz2T
from .geometric import rotate_NH_bond_vector_individual, rotate_NH_bond_vector

### Individual anisotropic Spectral Density Functions
def J_anitop_aniDr(w, NH_bond_vectors, angles, Ds):
    #Dperp = Dxx = Dyy and Dpar = Dzz
    phi, theta, psi = angles
    x, y, z = (rotate_NH_bond_vector_individual(NH_bond_vectors, [phi, theta, psi])) # actually l, m, n

    Dxx, Dyy, Dzz = Ds
   
    Diso = (Dxx + Dzz + Dzz)/3
    L = np.sqrt(((Dxx*Dyy)+(Dxx*Dzz)+(Dyy+Dzz))/3)
    
    dx = (Dxx - Diso)/(np.sqrt((Diso**2)-(L**2)))
    dy = (Dyy - Diso)/(np.sqrt((Diso**2)-(L**2)))
    dz = (Dzz - Diso)/(np.sqrt((Diso**2)-(L**2)))
    
    t1 = 1/(4*Dxx + Dyy + Dzz)
    t2 = 1/(Dxx + 4*Dyy + Dzz)
    t3 = 1/(Dxx + Dyy + 4*Dzz)
    t4 = 1/(6*Diso + 6*(np.sqrt((Diso**2)-(L**2))))
    t5 = 1/(6*Diso - 6*(np.sqrt((Diso**2)-(L**2))))
    # for each residue i:    
    A1i = 3*(y**2)*(z**2)
    A2i = 3*(x**2)*(z**2)
    A3i = 3*(x**2)*(y**2)
    A4i = (1/4)*(3*((x**4)+(y**4)+(z**4))-1)-(1/12)*((dx*((3*(x**4))+(6*(y**2)*(z**2)))-1)+(dy*((3*(y**4))+(6*(z**2)*(x**2)))-1)+(dz*((3*(z**4))+(6*(x**2)*(y**2)))-1))
    A5i = (1/4)*(3*((x**4)+(y**4)+(z**4))-1)+(1/12)*((dx*((3*(x**4))+(6*(y**2)*(z**2)))-1)+(dy*((3*(y**4))+(6*(z**2)*(x**2)))-1)+(dz*((3*(z**4))+(6*(x**2)*(y**2)))-1))
    
    return ((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2)) + ((A4i*t4)/(1+(w*t4)**2)) + ((A5i*t5)/(1+(w*t5)**2))

def J_anitop_aniDr_oblate(w, NH_bond_vectors, angles, Ds):
    #Dperp = Dxx = Dyy and Dpar = Dzz
    phi, psi, theta = angles
    x, y, z = (rotate_NH_bond_vector_individual(NH_bond_vectors, [phi, theta, psi])) # actually l, m, n
    Dxx, Dyy, Dzz = Ds
    Dpar = Dxx
    Dperp = Dzz
       
   #  Diso = (Dpar + 2*Dperp)/3
    t1 = 1/(6*Dperp)
    t2 = 1/(5*Dperp+Dpar)
    t3 = 1/(2*Dperp + 4*Dpar)
    # for each residue i
    A1i = ((3*(np.cos(theta)**2)-1)**2)/4
    A2i = 3*(np.sin(theta)**2)*(np.cos(theta)**2)
    A3i = (3/4)*(np.sin(theta)**4)
    return ((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2))

def J_anitop_aniDr_prolate(w, NH_bond_vectors, angles, Ds):
    #Dperp = Dxx = Dyy and Dpar = Dzz
    phi, psi, theta = angles
    x, y, z = (rotate_NH_bond_vector_individual(NH_bond_vectors, [phi, theta, psi])) # actually l, m, n
    Dxx, Dyy, Dzz = Ds
    Dpar = Dzz
    Dperp = Dxx
       
   #  Diso = (Dpar + 2*Dperp)/3
    t1 = 1/(6*Dperp)
    t2 = 1/(5*Dperp+Dpar)
    t3 = 1/(2*Dperp + 4*Dpar)
    # for each residue i
    A1i = ((3*(np.cos(theta)**2)-1)**2)/4
    A2i = 3*(np.sin(theta)**2)*(np.cos(theta)**2)
    A3i = (3/4)*(np.sin(theta)**4)
    return ((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2)) 

def J(w, NH_bond_vectors, angles, Ds, S2, ti, model): 
    phi, theta, psi = angles
    Dxx, Dyy, Dzz = Ds
    S2 = S2
    ti = ti
    
    if model in [1,2,3,4,5,6]:
        x, y, z = (rotate_NH_bond_vector_individual(NH_bond_vectors, [phi, theta, psi])) # actually l, m, n
    if model in [7,8,9,10,11,12]:
        x, y, z = (rotate_NH_bond_vector(NH_bond_vectors, [phi, theta, psi])) # actually l, m, n
    
    if model in [1,2,4,5,7,8,10,11]:
        
        A1i = (1/4)*((3*(x**2)-1)**2)
        A2i = 3*(x**2)*(1-x**2)
        A3i = (3/4)*(((x**4)-1)**2)
        
        if model in [1,4,7,10]: #oblate
            Dpar = Dxx
            Dperp = Dzz
        if model in [2,5,8,11]: #prolate
            Dpar = Dzz
            Dperp = Dxx 
            
        Diso = (Dpar + 2*Dperp)/3        
        
        t1 = 1/(6*Dperp)
        t2 = 1/(5*Dperp+Dpar)
        t3 = 1/(2*Dperp + 4*Dpar)
        
        if model in [1,2,4,5]: #oblate and prolate, with or without Rex
            return ((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2))

        if model in [7,8,10,11]: #oblate and prolate, with S2 and with or without Rex
            tau = (1/(6*Diso))+ti
            return (S2*(((A1i*t1)/(1+(w*t1)**2))+((1-S2)*(tau/(1+(w+tau)**2))))) + (S2*(((A2i*t2)/(1+(w*t2)**2))+((1-S2)*(tau/(1+(w+tau)**2))))) + (S2*(((A3i*t3)/(1+(w*t3)**2))+((1-S2)*(tau/(1+(w+tau)**2)))))
    

    if model in [3,6,9,12]: ##fully anisotropic
        Diso = (Dxx + Dyy + Dzz)/3
        L = np.sqrt(((Dxx*Dyy)+(Dxx*Dzz)+(Dyy+Dzz))/3)
        
        dx = (Dxx - Diso)/(np.sqrt((Diso**2)-(L**2)))
        dy = (Dyy - Diso)/(np.sqrt((Diso**2)-(L**2)))
        dz = (Dzz - Diso)/(np.sqrt((Diso**2)-(L**2)))
        
        t1 = 1/(4*Dxx + Dyy + Dzz)
        t2 = 1/(Dxx + 4*Dyy + Dzz)
        t3 = 1/(Dxx + Dyy + 4*Dzz)
        t4 = 1/(6*Diso + 6*(np.sqrt((Diso**2)-(L**2))))
        t5 = 1/(6*Diso - 6*(np.sqrt((Diso**2)-(L**2))))
        
        # for each residue i:    
        A1i = 3*(y**2)*(z**2)
        A2i = 3*(x**2)*(z**2)
        A3i = 3*(x**2)*(y**2)
        A4i = (1/4)*(3*((x**4)+(y**4)+(z**4))-1)-(1/12)*((dx*((3*(x**4))+(6*(y**2)*(z**2)))-1)+(dy*((3*(y**4))+(6*(z**2)*(x**2)))-1)+(dz*((3*(z**4))+(6*(x**2)*(y**2)))-1))
        A5i = (1/4)*(3*((x**4)+(y**4)+(z**4))-1)+(1/12)*((dx*((3*(x**4))+(6*(y**2)*(z**2)))-1)+(dy*((3*(y**4))+(6*(z**2)*(x**2)))-1)+(dz*((3*(z**4))+(6*(x**2)*(y**2)))-1))
        
        if model in [3,6]:
            return ((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2)) + ((A4i*t4)/(1+(w*t4)**2)) + ((A5i*t5)/(1+(w*t5)**2))
        if model in [9,12]:
            tau = (1/(6*Diso))+ti
            # return S2*(((A1i*t1)/(1+(w*t1)**2)) + ((A2i*t2)/(1+(w*t2)**2)) + ((A3i*t3)/(1+(w*t3)**2)) + ((A4i*t4)/(1+(w*t4)**2)) + ((A5i*t5)/(1+(w*t5)**2))) + ((1-S2)*(tau/(1+(w+tau)**2)))
            return (S2*((A1i*t1)/(1+(w*t1)**2))+((1-S2)*(tau/(1+(w+tau)**2)))) + (S2*((A2i*t2)/(1+(w*t2)**2))+((1-S2)*(tau/(1+(w+tau)**2)))) + (S2*((A3i*t3)/(1+(w*t3)**2))+((1-S2)*(tau/(1+(w+tau)**2)))) + (S2*((A4i*t4)/(1+(w*t4)**2))+((1-S2)*(tau/(1+(w+tau)**2)))) + (S2*((A5i*t5)/(1+(w*t5)**2))+((1-S2)*(tau/(1+(w+tau)**2))))
        
        
def relax(par, NH_bond_coordinates, fields):
    model = int(par["model"].value)
    angles = [par["phi"].value, par["theta"].value, par["psi"].value]
    Ds = [par["Dxx"].value, par["Dyy"].value, par["Dzz"].value]
    Rex = par["Rex"].value
    S2 = par["S2"].value
    ti = par["ti"].value
    
    # theta_ddcsa = 22
    # else:
    #     theta_ddcsa = par["theta_ddcsa"].value
    for k in range(len(NH_bond_coordinates)):
        NH_bond_vectors = NH_bond_coordinates.iloc[k]
        for i in fields:
            fieldT = MHz2T(i)
            wh = g1H * fieldT
            wn = g15N * fieldT
            csa=(wn*dCSA)
            # wh = params["wh"].value
            # wn = params["wn"].value

            R1 = ((d**2)/10)*(J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)))+((2/15)*((csa**2))*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))
            R2 = (Rex*(i**2)) + ((d**2)/20)*((4*J(0, NH_bond_vectors, angles, Ds, S2, ti, model))+J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)))+(((csa**2)/45)*((4*J(0, NH_bond_vectors, angles, Ds, S2, ti, model))+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))))
            # NOE = 1.0 +(d**2)/(10*R1)*(g1H/g15N)*(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)-(J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)))
            # sNH = (1/10)*(d**2)*((6*J(wh+wn, angles, Ds, S2, ti, model))-(J(wh-wn,angles)))
            # etaZ = -(1/15)*csa*d*(P2(np.cos(theta_ddcsa)))*(6*J(wn, angles, Ds, S2, ti, model))
            # etaXY = -(1/15)*csa*d*(P2(np.cos(theta_ddcsa)))*(4*(J(0, NH_bond_vectors, angles, Ds, S2, ti, model)+3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model)))

            yield R1, R2
        
def relax2(par, NH_bond_coordinates, fields):
    model = int(par["model"].value)
    angles = [par["phi"].value, par["theta"].value, par["psi"].value]
    Ds = [par["Dxx"].value, par["Dyy"].value, par["Dzz"].value]
    # Rex = par[f"Rex_{res}"].value
    # S2 = par[f"S2_{res}"].value
    # ti = par[f"ti_{res}"].value
    
    Rex = par["Rex"].value
    S2 = par["S2"].value
    ti = par["ti"].value
        
    theta_ddcsa = par['theta_ddcsa'].value
    # else:
    #     theta_ddcsa = par["theta_ddcsa"].value
    # NH_bond_vectors= NH_bond_coordinates.iloc[res][["NH_x", "NH_y", "NH_z"]]
    NH_bond_vectors= NH_bond_coordinates
    for i in fields:
        fieldT = MHz2T(i)
        wh = g1H * fieldT
        wn = g15N * fieldT
        csa=(wn*dCSA)
        # wh = params["wh"].value
        # wn = params["wn"].value

        R1 = ((d**2)/10)*(J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)))+((2/15)*((csa**2))*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))
        R2 = (Rex*(i**2)) + ((d**2)/20)*((4*J(0, NH_bond_vectors, angles, Ds, S2, ti, model))+J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh, NH_bond_vectors, angles, Ds, S2, ti, model))+(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)))+(((csa**2)/45)*((4*J(0, NH_bond_vectors, angles, Ds, S2, ti, model))+(3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model))))
        NOE = 1.0 +(d**2)/(10*R1)*(g1H/g15N)*(6*J(wh+wn, NH_bond_vectors, angles, Ds, S2, ti, model)-(J(wh-wn, NH_bond_vectors, angles, Ds, S2, ti, model)))
        # sNH = (1/10)*(d**2)*((6*J(wh+wn, angles, Ds, S2, ti, model))-(J(wh-wn,angles)))
        # etaZ = -(1/15)*csa*d*(P2(np.cos(theta_ddcsa)))*(6*J(wn, angles, Ds, S2, ti, model))
        etaXY = -(1/15)*csa*d*(P2(np.cos(theta_ddcsa)))*(4*(J(0, NH_bond_vectors, angles, Ds, S2, ti, model)+3*J(wn, NH_bond_vectors, angles, Ds, S2, ti, model)))

        yield R1, R2, NOE, etaXY