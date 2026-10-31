<<<<<<< Updated upstream
=======
from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants


#constants
CF = 2
Vdot_LOx = 0.00115512708 #m^3/s, max
Vdot_IPA = 0.00114452214 #m^3/s, max
P_t = 3548704 #Pa, postreg P
P_tmin = 6894757.293 #Pa, min P into reg
V_b = 0.0131 #m^3
#R = scipy.constants.

T_b = 316 #K
P_b = 21766190.54 #Pa
m_b = 2.8113 #kg
m_used = 1.2 #kg
m_pres = 0.2641458904 #kg
V_uLOx = 0.00340 #m^3
V_uIPA = 0.00307 #m^3
#m_res =  #kg
#m_LOx =  #kg
#m_IPA =  #kg
rho_b = (m_b-m_used)/V_b #kg/m^3
S = propssi("S", "T", T_b, "P", P_b, "nitrogen") #J/kgK

#transient
mdots_LOx = []
mdots_IPA = []
mdots_total = []

time = 0 #s
dt = 0.1 #s
endtime = 10 #s

while (time <= endtime) & (m_b-m_used > 0) & (P_b > P_tmin):
    h_b = propssi("Hmass", "S", S, "D", rho_b, "nitrogen") #isentropic flow, specific H bottle = specific H tanks
    T_in = propssi("T", "P", P_t, "H", h_b, "nitrogen") #T going into prop tanks
    rho_in = propssi("D", "P", P_t, "T", T_in, "nitrogen") #rho GN2 going into tanks
    
    mdot_LOx = Vdot_LOx * CF * rho_in
    mdot_IPA = Vdot_IPA * CF * rho_in
    mdot_total = mdot_LOx + mdot_IPA
    mdots_LOx.append(mdot_LOx)
    mdots_IPA.append(mdot_IPA)
    mdots_total.append(mdot_total)

    m_used += mdot_total * dt
    rho_b = (m_b-m_used)/V_b #kg/m^3
    P_b = propssi("P", "S", S, "D", rho_b, "nitrogen")
    print(P_b)
    time += dt

print(f"m_used: {m_used} kg in {time}")
print(f"LOx max m dot: {max(mdots_LOx)}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt} s")
print(f"IPA max m dot: {max(mdots_IPA)}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt} s")
print(f"Total max m dot: {max(mdots_total)}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt} s")



>>>>>>> Stashed changes
