from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants

#constants
CF = 2
Vdot_LOx = 0.00116 #m^3/s, max
Vdot_IPA = 0.00114 #m^3/s, max
P_t = 3548704 #Pa, postreg P
R = scipy.constants.

#initial values
V_b = 0.0131 #m^3
V_uLOx = 0.00340 #m^3
V_uIPA = 0.00307 #m^3
m_res =  #kg
m_LOx0 =  #kg
m_IPA0 =  #kg
Tb0 = 316 #K
P_postpress = 21766190.54 #Pa
mb0 = 2.8113 #kg
rho0 =  #kg/m^3
S = propssi("S", "T", Tb0, "P", P_postpress, "nitrogen") #J/kgK
time = 0 #s

#transient
def Findmused(m_used):


    return m_used


