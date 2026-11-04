from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import csv

GN2_R = 296.84 #J/kgK
P_0 = 3500*6894.757 #upstream pressure
T_0 = 316 #K
setP = 550*6894.757 #set pressure
S_y = 25000*6894.757
UTS = 70000*6894.757

mdot = 0.31308797 #kg/s, max mdot timestep. wont be fully accurate, need to use m dot from COPV.py after selection. COME BACK HERE
mdot1 = 0.2026508758 #kg/s, initial timestep

#upstream stag @ max mdot
T_u = 200 #K
P_u = 1108*6894.757 #Pa
#upstream stag @ first timestep
T_u1 = 316.48 #K
P_u1 = 3220*6894.757 #Pa
#downstream stag @ max mdot
T_d = 195.17 #K
P_d = setP #Pa
#downstream stag @ first timestep
T_d1 = 316 #K
P_dF = 1.1*setP #failure downstream pressure

def getfloat(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Pls try again: ")

def linesize():
    OD = getfloat("Line OD (in): ")
    thickness = getfloat("Thickness (in): ")
    ID = OD - 2*thickness
    OD *= 0.0254
    ID *= 0.0254

    area = constants.pi*((ID/2)**2)
    P_y = S_y*((OD**2)-(ID**2))/((OD**2)+(ID**2)) #yield pressure
    P_b = UTS*((OD**2)-(ID**2))/((OD**2)+(ID**2)) #burst pressure
    return area, P_y, P_b, OD, ID, thickness

def calcMach_Kappa(P0, T0, mdot, area):
    M_guess = 0.3
    rho0 = propssi("D", "P", P0, "T", T0, "nitrogen")

    k = propssi("D", "P", P0, "T", T0, "nitrogen")*(propssi("A", "P", P0, "T", T0, "nitrogen")**2)/P0
    P = P0/(1+(k-1)/2*(M_guess**2))**(k/(k-1))
    T = T0/(1+(k-1)/2*(M_guess**2))
    rho = rho0/(1+(k-1)/2*(M_guess**2))**(1/(k-1))

    velo = (mdot/rho)/area
    A = propssi("A", "T", T, "P", P, "nitrogen")
    M = velo/A
    error = abs(M_guess-M)

    for i in range(200):
        if error < 0.0001:
            break
        k = propssi("D", "P", P, "T", T, "nitrogen")*(propssi("A", "P", P, "T", T, "nitrogen")**2)/P
        P = P0/(1+(k-1)/2*(M**2))**(k/(k-1))
        T = T0/(1+(k-1)/2*(M**2))
        rho = rho0/(1+(k-1)/2*(M**2))**(1/(k-1))
    
        velo = (mdot/rho)/area
        A = propssi("A", "T", T, "P", P, "nitrogen")
        M0=M
        M = velo/A
        error = abs(M0-M)
    else:
        M = "did not converge"
    return M

while True:
    size = linesize()
    area = size[0]
    P_y = size[1]
    P_b = size[2]
    OD = size[3]
    ID = size[4]
    thickness = size[5]

    type = input("Upstream or downstream: ")
    while type.lower() != "upstream" and type.lower() != "downstream":
        type = input("Upstream or downstream: ")

    if type.lower() == "upstream":
        M = calcMach_Kappa(P_u, T_u, mdot, area)
        P_Sy = P_y/P_u
        P_Sb = P_b/P_u
        print(f"Upstream normal operation max mdot using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")

        M = calcMach_Kappa(P_u1, T_u1, mdot1, area)
        P_Sy = P_y/P_u1
        P_Sb = P_b/P_u1
        print(f"Upstream normal operation first timestep using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")
    
    elif type.lower() == "downstream":
        M = calcMach_Kappa(P_d, T_d, mdot, area)
        P_Sy = P_y/P_d
        P_Sb = P_b/P_d
        print(f"Downstream normal operation max mdot using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")

        M = calcMach_Kappa(P_d, T_d1, mdot1, area)
        P_Sy = P_y/P_d
        P_Sb = P_b/P_d
        print(f"Downstream normal operation first timestep using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")

        M = calcMach_Kappa(P_dF, T_d, mdot, area)
        P_Sy = P_y/P_dF
        P_Sb = P_b/P_dF
        print(f"Downstream failure max mdot using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")

        M = calcMach_Kappa(P_dF, T_d1, mdot1, area)
        P_Sy = P_y/P_dF
        P_Sb = P_b/P_dF
        print(f"Downstream failure first timestep using kappa: M={M:.3f}, Yield={P_Sy:.3f}, Burst={P_Sb:.3f}")
    print("----------------------------------------------------------------------------------------")



