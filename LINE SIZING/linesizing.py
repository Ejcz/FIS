from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import csv

GN2_R = 296.84 #J/kgK
P_0 = 3500*6894.757 #upstream pressure
setP = 550*6894.757 #set pressure
S_y = 25000*6894.757
UTS = 70000*6894.757

mdot_up = 0.31308797 #kg/s, max mdot timestep
T_u = 200 #K, upstream stag T
T_d = 195.17 #K
P_u = 686*6894.757 #Pa, upstream stag P
P_d = setP*6894.757 #Pa

mdot_up0 = 0.2026508758 #kg/s, initial timestep 
T_u0 = 316.48 #K
T_d0 = 306.35 #K
P_u0 = 3157*6894.757 #Pa
P_d0 = setP*6894.757
P_dF = 1.1*setP*6894.757 #failure downstream pressure

S_u = propssi("S", "T", T_u, "P", P_u, "nitrogen")
h_u = propssi("HMASS", "T", T_u, "P", P_u, "nitrogen")
S_u0 = propssi("S", "T", T_u0, "P", P_u0, "nitrogen")
h_u0 = propssi("HMASS", "T", T_u0, "P", P_u0, "nitrogen")

def size():
    OD = input("Line OD: ")
    while OD.isdigit() == False:
        OD = input("Line OD: ")
    thickness = input("Thickness: ")
    while thickness.isdigit() == False:
        thickness = input("Thickness: ")
    ID = OD - 2*thickness

    area = constants.pi*((ID/2)**2)
    P_y = S_y*(OD**2-ID**2)/(OD**2+ID**2) #yield pressure
    P_b = UTS*(OD**2-ID**2)/(OD**2+ID**2) #burst pressure
    return area, P_y, P_b

def calcMach(P, T, ):
    M_i = 0.3


    #upstream first timestep




while True:
    output = size()
    calcMach()

