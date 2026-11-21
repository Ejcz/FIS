from CoolProp.CoolProp import PropsSI as propssi
import CoolProp.CoolProp as CP
from scipy import integrate, optimize, constants
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
import csv
import os
import io
import contextlib
from multiprocessing import Pool

#come back here: debug the end temperature or contour being weird
#figure is created inside the __main__ block at the bottom (worker processes import this file and must not make plots)
x = []
y = []
x1 = []
y1 = []
colors = []
mdots = []
mdots1 = []
m_useds = []
fill_masses = []
fill_masses1 = []
T_b1s = []
T_b1s1 = []
P_b1s = []
P_b1s1 = []
press_masses = []
press_masses1 = []
T_fs = []
T_fs1 = []
T_bfs = []
Tcolors = []
Tcolors1 = []
all_mdots = []
all_mdots1 = []
choke_mdots1 = []
choke_mdots2 = []
LOx_mdots1 = []
LOx_mdots2 = []
IPA_mdots1 = []
IPA_mdots2 = []
excesses = []
expl_times = []
part_nums = []
COPV_ODs = []
statuses = []
bweights = []
BPs = []
BTs = []
choke_times = []
times1 = []
Tins1 = []

x2 = []
y2 = []
mdots2 = []
BPs2 = []
BTs2 = []
excesses2 = []
fill_masses2 = []
press_masses2 = []
P_b1s2 = []
T_b1s2 = []
T_fs2 = []
T_bfs2 = []
Tcolors2 = []
all_mdots2 = []
expl_times2 = []
part_nums2 = []
COPV_ODs2 = []
statuses2 = []
bweights2 = []
choke_times2 = []
times2 = []
Tins2 = []

#constants/variables
Vdot_LOx = 0.001111306 #m^3/s, nominal
mdot_LOx = 1.268 #kg/s
CF = 1.5 #adjust: collapse factor
m_LOx = 10*mdot_LOx #kg, nominal
Vdot_IPA = 0.001086247 #m^3/s, nominal
mdot_IPA = 0.932 #kg/s
m_IPA = 9.82 #kg, max
P_t = 3548704 #Pa, postreg P. inconsistent? come back here
P_tmin = 6894757.293 #Pa, min P into reg
T_amb = 300 #K adjust: ambient temperature
V_uLOx = 0.00340 #m^3
V_uIPA = 0.00307 #m^3
V_utotal = 0.00699 #m^3
T_b0 = 300 #K adjust: initial bottle temperature
dt = 0.05 #s
depl_time = max(m_LOx/mdot_LOx, m_IPA/mdot_IPA) #s
P_step = 20 #psi step
V_error = 0.0001
endPressureError = 5 #psi
m_guess_prev = 2.2 #kg
V_max = 29 #L cap on volume sweep
T_rating = 215.15 #K
T_tol = 0 #K, come back here
press = "idk"
P_set = 550*6894.757 #adjust: regulator setpoint 
P_ullage_min = 3.4 * 10**6
P_ullage_max = 750*6894.757

D_orifice = 0.141 #in adjust: orifice diameter
A_orifice = ((D_orifice*0.0254/2)**2)*3.141592
C_d = 0.62 #flow coefficient, check api 520/521
A_orifice_eff = A_orifice * C_d
LID = 0.26*0.0254 #tubing ID in m
LA = (LID/2)**(2)*3.141592 #tubing area
P_loss = 0 #pressure loss across orifice

P_XTOL = 1.0 #Pa, brentq tolerance for pressure solves (default was ~1e-12 Pa, way past what matters)

#CoolProp low-level state: much faster than PropsSI (no string parsing / fluid setup per call).
#One shared state per process; read the values you need right after each update.
N2 = CP.AbstractState("HEOS", "Nitrogen")
def N2_PS(P, S): #state from pressure + entropy
    N2.update(CP.PSmass_INPUTS, P, S)
    return N2
def N2_DS(rho, S): #state from density + entropy
    N2.update(CP.DmassSmass_INPUTS, rho, S)
    return N2
def N2_HP(h, P): #state from enthalpy + pressure
    N2.update(CP.HmassP_INPUTS, h, P)
    return N2

#ullage starting state for fill_mass, constant so computed once
rho_u0 = propssi("D", "P", 101325, "T", T_amb, "nitrogen")
u_u0 = propssi("UMASS", "T", T_amb, "P", 101325, "air") #note: air and nitrogen use different CoolProp reference states, see review notes


#filename = "C:/Users/Chris/Documents/Github/FIS/COPV/GN2_COPV_Options.csv" #adjust: COPV options filename
filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "GN2_COPV_Options.csv")
def readCSV(filename, colname):
    data = []
    with open(filename, newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index(colname)
    rows = rows[1:]
    for row in rows:
        data.append(row[col])
    return data


rho_hi_fill = 2*propssi("D", "P", P_t, "T", T_amb, "nitrogen") #for fill_mass upper bound

def fill_mass(V_u, h_in):
    m_0 = V_u * rho_u0
    def f(dm):
        m_ppres = m_0 + dm
        u_f = (m_0*u_u0 + dm*h_in)/m_ppres
        N2.update(CP.DmassUmass_INPUTS, m_ppres/V_u, u_f)
        return N2.p() - P_t
    dm_hi = rho_hi_fill*V_u #just a safe upper bound, can mostly ignore
    return optimize.brentq(f, 0, dm_hi, xtol=1e-9) #kg

class ThroatCondenses(Exception):
    pass #claude help

class OrificeChokes(Exception):
    pass #claude help

def solve_sonic(S, P_b, h):
    def f(P):
        try:
            st = N2_PS(P, S) #one state update gives both A and h
            A = st.speed_sound()
            hP = st.hmass()
        except ValueError:
            return 1e6
        return max(0, 2*(h-hP))**(1/2) - A
    hi = 0.9999*P_b

    if f(hi) >= 0:
        raise ValueError("solve_sonic: flow already sonic at bottle pressure")
    for frac in (0.4, 0.2, 0.05):
        lo = frac * P_b
        if f(lo) > 0: #low guess is supersonic
            P_fill = optimize.brentq(f, lo, hi, xtol=P_XTOL)
            if abs(f(P_fill)) > 3:
                raise ThroatCondenses
            return P_fill
    raise ValueError("solve_sonic: no sonic point found")

def solve_throatP(Area, h, mdot, S, sonicP, P_b):
    def resid(P): #continuity velocity - energy velocity, one state update per guess
        st = N2_PS(P, S)
        v_cont = mdot/(st.rhomass()*Area)
        v_energy = max(0, 2*(h-st.hmass()))**(1/2)
        return v_cont - v_energy
    throat_P = optimize.brentq(resid, sonicP, 0.999999*P_b, xtol=P_XTOL)
    return throat_P, mdot/(N2_PS(throat_P, S).rhomass()*Area) #are brackets ok?

def choke(Area, S, h, P_b, sonicP = None): #pass sonicP if already solved, saves a second solve_sonic
    P = solve_sonic(S, P_b, h) if sonicP is None else sonicP
    st = N2_PS(P, S)
    return Area*st.rhomass()*st.speed_sound() #choke flow rate

def solve_recovery(LA, mdot, P_o, v_o, h_0, P_up):
    S = N2_HP(h_0, P_up).smass() #upstream conditions
    h_o = N2_PS(P_o, S).hmass() #specific enthalpy thru orifice
    v_o = (2*(h_0-h_o))**(1/2)

    rho_d = N2_HP(h_o, P_o).rhomass() #initial guess
    v_d = mdot/(rho_d*LA) #downstream velocity
    P_d0 = 0 #initial state
    error = 1e9 #placeholder
    for i in range(100): #was an unbounded while loop with 0.0001 Pa tolerance: CoolProp's own noise (~1e-3 Pa at 2500 psi) could keep it from ever stopping
        h_d = h_0 - 0.5*(v_d)**(2)
        P_d = P_o + mdot/(LA)*(v_o-v_d)
        rho_d = N2_HP(h_d, P_d).rhomass()
        v_d = mdot/(rho_d*LA)

        error = abs(P_d - P_d0)
        P_d0 = P_d
        if error < P_XTOL: #1 Pa
            return P_d
    raise ValueError(f"solve_recovery: no convergence after 100 iterations (last change {error:.3g} Pa)")


#burn start to end
def size_bottle(P_b0, m_usedguess = 2.2): #P_b0 = before prepressurization
    #constants
    m_used = 0 #kg
    m_pres = 0.2641458904 #kg, initial guess for pressurization mass
    error = 1 #placeholder

    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3
    u_b0 = N2_DS(rho_b0, S).umass() #specfic internal energy before pressurization, same every iteration

    while error > V_error:
        time = 0 #s
        m_used = 0 #kg
        mdots_LOx = [] #lists are for GN2 btw
        mdots_IPA = []
        mdots_total = []

        V_b = m_usedguess/(rho_b0-rho_end) #m^3, guess bottle volume based on m_used
        rho_b1 = rho_b0-m_pres/V_b #density after pressurization, before burn
        u_b1 = N2_DS(rho_b1, S).umass() #specfic internal energy after pressurization
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres #enthalpy of GN2 going into reg (and out of bottle) during pressurization, J/kg

        # m_presLOx = fill_mass(V_uLOx, h_in) #no line volume, need to play around with CF. kinda sketch
        # m_presIPA = fill_mass(V_uIPA, h_in)
        # m_pres = m_presLOx * 1.5 + m_presIPA
        m_pres = fill_mass(V_utotal, h_in)

        rho_b = rho_b0-(m_used+m_pres)/V_b #recalculated after properly solving m_pres
        rho_b1 = rho_b0-m_pres/V_b #this is for the start of burn btw
        P_b1 = propssi("P", "D", rho_b1, "S", S, "nitrogen")
        T_b1 = propssi("T", "P", P_b1, "D", rho_b1, "nitrogen")
        
        while (time <= depl_time):
            h_b = N2_DS(rho_b, S).hmass() #isentropic flow, specific H bottle = specific H tanks
            rho_in = N2_HP(h_b, P_t).rhomass() #rho GN2 going into tanks

            mdot_LOx = Vdot_LOx * CF * rho_in #flow of GN2 into LOx tank, kg/s does velocity decrease when paths split (ignoring friction)? come back here
            mdot_IPA = Vdot_IPA * rho_in #flow of GN2 into IPA tank, kg/s
            mdot_total = mdot_LOx + mdot_IPA

            mdots_LOx.append(mdot_LOx)
            mdots_IPA.append(mdot_IPA)
            mdots_total.append(mdot_total)

            m_used += mdot_total * dt
            rho_b = rho_b0-(m_used+m_pres)/V_b #kg/m^3

            time += dt
        error = abs(m_usedguess-m_used-m_pres)
        m_usedguess = m_used+m_pres

    P_b = propssi("P", "S", S, "D", rho_b, "nitrogen")
    T_f = propssi("T", "S", S, "D", rho_b, "nitrogen")

    if abs(P_b-P_tmin) <= endPressureError*6894.757: #5 psi
        endreached = True
    else:
        endreached = False

    # print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    # print(f"V_b: {1000*V_b:.2f}L, fill mass: {rho_b0*V_b:.3f} kg, pressurization mass: {m_pres:.3f} kg\nP_b0: {P_b0/6894.757:.0f} psi, P_b1: {P_b1/6894.757:.0f} psi, T_b1: {T_b1:.0f} K, T_f: {T_f:.0f}\nMin pressure reached: {endreached}")
    # print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    # print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    # print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    # print("---------------------------------------------------------------------------------------")
    output = {
        "V_b":V_b,
        "mdots_LOx_max":max(mdots_LOx),
        "mdots_IPA_max":max(mdots_IPA),
        "mdots_total_max":max(mdots_total),
        "endreached":endreached,
        "fill_mass":rho_b0*V_b,
        "P_b1":P_b1,
        "T_b1":T_b1,
        "m_pres":m_pres,
        "T_f":T_f,
        "mdots_total":mdots_total,
        "m_used":m_used
    }
    return output

#Function that tries bottles and appends lists for plotting
def tryP_b0(P_b0):
    global m_guess_prev
    output = size_bottle(P_b0, m_guess_prev)
    m_guess_prev = output["m_used"] + output["m_pres"]
    x.append(P_b0/6894.757293) #psi
    y.append(1000*output["V_b"]) #L
    mdots.append((output["mdots_LOx_max"], output["mdots_IPA_max"], output["mdots_total_max"]))
    if output["endreached"] == True:
        colors.append("blue")
    else:
        colors.append("red")
    fill_masses.append(output["fill_mass"])
    press_masses.append(output["m_pres"])
    P_b1s.append(output["P_b1"]/6894.757293) #psi
    T_b1s.append(output["T_b1"])
    T_fs.append(output["T_f"])
    if output["T_f"] > T_rating:
        Tcolors.append("blue")
    else:
        Tcolors.append("red")
    all_mdots.append(output["mdots_total"])
    m_useds.append(output["m_used"])

#trying and plotting COPVs
def use_bottle(P_b0, V_b): #come back here for orifice (mdot_reg < mdot_total? check reg flow curve)
    time = 0 #s
    m_used = 0 #kg
    m_pres = 0.08 #kg, initial guess
    mdots_LOx = []
    mdots_IPA = []
    mdots_total = []
    mdots_choke = []
    Ps = []
    Ts = []
    Tins = []
    times = []
    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3

    # m_presLOx = fill_mass(V_uLOx, h_in) #no line volume COME BACK HERE
    # m_presIPA = fill_mass(V_uIPA, h_in)
    # m_pres = m_presLOx * 1.5 + m_presIPA #add a check on this? come back here, kinda sketch
    m_avail = V_b*(rho_b0-rho_end)
    for i in range(50):
        m_pres = min(m_pres, m_avail*0.99999)
        rho_b1 = rho_b0-m_pres/V_b
        u_b0 = N2_DS(rho_b0, S).umass()
        u_b1 = N2_DS(rho_b1, S).umass()
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres
        m_new = fill_mass(V_utotal, h_in) #can't account for CF_pp of LOx with this
        done = abs(m_pres - m_new) < 0.00001
        m_pres = m_new
        if done:
            break

    if m_pres > m_avail: #come back here to figure out margin + m_pres
        return None 
            # "rho_margin":(m_avail-m_pres)/V_b,
            # "T_bf":propssi("T", "S", S, "D", rho_end, "nitrogen"),
            # "m_pres":m_pres
            

    m_used = 0 #kg
    rho_b = rho_b0-(m_used+m_pres)/V_b #need to recalculate these after m_pres calc
    rho_b1 = rho_b0-m_pres/V_b #these are for the start of burn btw
    P_b1 = propssi("P", "D", rho_b1, "S", S, "nitrogen")
    T_b1 = propssi("T", "P", P_b1, "D", rho_b1, "nitrogen")

    P_b = P_b1
    P_throat = 0
    t_reg = None
    t_choke = None
    min_choke_ratio = np.inf
    min_P_inreg = np.inf
    P_inreg = P_b
    P_inregs = []
    while (time < depl_time) & (P_inreg > P_tmin ): #come back here: check outputs with check valve and also line flow
        st = N2_DS(rho_b, S)
        h_b = st.hmass() #isentropic flow, specific H bottle = specific H tanks
        T_b = st.T() #temp of gn2 in bottle
        st = N2_HP(h_b, P_t)
        T_in = st.T() #T going into prop tanks. check this at the end? come back here
        rho_in = st.rhomass() #rho GN2 going into tanks

        mdot_LOx = Vdot_LOx * CF * rho_in
        mdot_IPA = Vdot_IPA * rho_in
        mdot_total = mdot_LOx + mdot_IPA

        try:
            sonicP = solve_sonic(S, P_b, h_b)
            mdot_choke = choke(A_orifice_eff, S, h_b, P_b, sonicP)
            choke_ratio = mdot_choke/mdot_total
            min_choke_ratio = min(choke_ratio, min_choke_ratio)
            if choke_ratio < 1: #orifice can't pass the demand flow: bottle fails here
                t_choke = time
                print(f"Orifice choked at {time:.3f}s and ratio {choke_ratio:.2f}")
                break
            P_throat, v_throat = solve_throatP(A_orifice_eff, h_b, mdot_total, S, sonicP, P_b)
            #deltaP = P_b - P_throat, come back here
            P_inreg = solve_recovery(LA, mdot_total, P_throat, v_throat, h_b, P_b)
        except ThroatCondenses:
            print(f"Throat condensed at {time:.3f}s and P_b {P_b/6894.757:.0f} psi")
            break
        except ValueError as e: #property/solver failure: stop this bottle instead of crashing the whole run
            print(f"Orifice solve failed at {time:.3f}s and P_b {P_b/6894.757:.0f} psi: {e}")
            break
        min_P_inreg = min(P_inreg, min_P_inreg)

        if (P_inreg < P_tmin) & (t_reg is None):
            t_reg = time

        mdots_LOx.append(mdot_LOx)
        mdots_IPA.append(mdot_IPA)
        mdots_total.append(mdot_total)
        Ps.append(P_b/6894.757)
        P_inregs.append(P_inreg/6894.757)
        Ts.append(T_b) #temp of gn2 in bottle
        Tins.append(T_in)
        times.append(time)
        mdots_choke.append(mdot_choke)

        m_used += mdot_total * dt
        rho_b = rho_b0-(m_used+m_pres)/V_b #kg/m^3
        P_b = N2_DS(rho_b, S).p()

        time += dt
    if not mdots_total:
        return None
    loss_end = (Ps[-1]-P_inregs[-1])*6894.757 #Pa, bottle P - reg inlet P at the last step (orifice loss)
    rho_end_eff = propssi("D", "S", S, "P", P_tmin + loss_end, "nitrogen") #bottle density where reg inlet hits P_tmin
    T_bf = propssi("T", "S", S, "D", max(rho_b, rho_end_eff), "nitrogen") #end of burn
    T_f = propssi("T", "S", S, "D", rho_end_eff, "nitrogen") #end of expulsion
    T_margin = T_f - T_rating
    m_excess = rho_b * V_b - rho_end_eff * V_b

    print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    print(f"V_b: {1000*V_b:.2f}L, fill mass: {rho_b0*V_b:.3f} kg, pressurization mass: {m_pres:.3f} kg\nP_b0: {P_b0/6894.757:.0f} psi, P_b1: {P_b1/6894.757:.0f} psi, T_b1: {T_b1:.0f} K, T_f: {T_f:.0f} K")
    print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    output = {
        "success":True,
        "P_margin":min_P_inreg-P_tmin,
        "rho_margin":rho_b-rho_end_eff,
        "choke_ratio":min_choke_ratio,
        "t_choke":t_choke,
        "t_reg":t_reg,
        "P_inregs":P_inregs,
        "T_margin":T_margin,
        "Ps":Ps,
        "mdots_LOx_max":max(mdots_LOx),
        "mdots_IPA_max":max(mdots_IPA),
        "mdots_total_max":max(mdots_total),
        "m_excess":m_excess,
        "fill_mass":rho_b0*V_b,
        "P_b1":P_b1,
        "T_b1":T_b1,
        "Tins":Tins,
        "m_pres":m_pres,
        "T_f":T_f,
        "T_bf":T_bf,
        "Ts":Ts,
        "mdots_choke":mdots_choke,
        "mdots_total":mdots_total,
        "mdots_LOx":mdots_LOx,
        "mdots_IPA":mdots_IPA,
        "time":time,
        "times":times,
        "t_choke":t_choke
        }
    return output

def stagnation_from_static(P, h0, mdot):
    v = 0.0
    for i in range(50):
        h = h0 - 0.5*(v**(2))
        rho = N2_HP(h, P).rhomass()
        v_new = mdot/(rho*LA)
        if abs(v - v_new) < 0.001:
            v = v_new
            break
        v = v_new
    else:
        raise ValueError("stagnation_from_static: velocity did not converge")

    h = h0 - 0.5*(v**(2))
    st = N2_HP(h, P) #static state: entropy and speed of sound
    S = st.smass()
    M = v/st.speed_sound()
    N2.update(CP.HmassSmass_INPUTS, h0, S)
    P0 = N2.p()
    return M, P0, S

def tank_state(P_guess, h_b): #loopy solve: P_tank --> rho into tank --> mdot into tank --> pressure drop and recovery --> P_tank = P_set - loss + recovery
    rho_in = N2_HP(h_b, P_guess).rhomass() #rho into tank
    mdot_LOx = Vdot_LOx * rho_in * CF
    mdot_IPA = Vdot_IPA * rho_in
    mdot = mdot_LOx + mdot_IPA
    M, P0, S = stagnation_from_static(P_set, h_b, mdot) #nonisentropic flow across reg
    sonicP = solve_sonic(S, P0, h_b)
    mdot_choke = choke(A_orifice_eff, S, h_b, P0, sonicP)
    ratio = mdot_choke/mdot
    if ratio <= 1:
        raise OrificeChokes(ratio)
    P_throat, v_throat = solve_throatP(A_orifice_eff, h_b, mdot, S, sonicP, P0)
    #deltaP = P_b - P_throat, come back here
    P_intank = solve_recovery(LA, mdot, P_throat, v_throat, h_b, P0)
    return P_intank, ratio, mdot, mdot_LOx, mdot_IPA, mdot_choke


def use_bottle1(P_b0, V_b): #come back here for orifice (mdot_reg < mdot_total? check reg flow curve)
    time = 0 #s
    m_used = 0 #kg
    m_pres = 0.08 #kg, initial guess
    mdots_LOx = []
    mdots_IPA = []
    mdots_total = []
    mdots_choke = []
    Ps = []
    Ts = []
    Tins = []
    times = []
    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3

    # m_presLOx = fill_mass(V_uLOx, h_in) #no line volume COME BACK HERE
    # m_presIPA = fill_mass(V_uIPA, h_in)
    # m_pres = m_presLOx * 1.5 + m_presIPA #add a check on this? come back here, kinda sketch
    m_avail = V_b*(rho_b0-rho_end)
    for i in range(50):
        m_pres = min(m_pres, m_avail*0.99999)
        rho_b1 = rho_b0-m_pres/V_b
        u_b0 = N2_DS(rho_b0, S).umass()
        u_b1 = N2_DS(rho_b1, S).umass()
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres
        m_new = fill_mass(V_utotal, h_in) #can't account for CF_pp of LOx with this
        done = abs(m_pres - m_new) < 0.00001
        m_pres = m_new
        if done:
            break

    if m_pres > m_avail: #come back here to figure out margin + m_pres
        return None
            # "rho_margin":(m_avail-m_pres)/V_b,
            # "T_bf":propssi("T", "S", S, "D", rho_end, "nitrogen"),
            # "m_pres":m_pres
            

    m_used = 0 #kg
    rho_b = rho_b0-(m_used+m_pres)/V_b #need to recalculate these after m_pres calc
    rho_b1 = rho_b0-m_pres/V_b #these are for the start of burn btw
    P_b1 = propssi("P", "D", rho_b1, "S", S, "nitrogen")
    T_b1 = propssi("T", "P", P_b1, "D", rho_b1, "nitrogen")

    P_b = P_b1
    t_choke = None
    min_choke_ratio = np.inf
    min_P_intank = np.inf
    P_intanks = []
    P_guess = P_set
    min_P_intank = np.inf
    max_P_intank = 0
    min_choke_ratio = np.inf
    while (time < depl_time) & (rho_b > rho_end): #come back here: check outputs with check valve and also line flow
        st = N2_DS(rho_b, S)
        h_b = st.hmass() #isentropic flow, specific H bottle = specific H tanks
        T_b = st.T() #temp of gn2 in bottle
        T_in = N2_HP(h_b, P_t).T() #T going into prop tanks. check this at the end? come back here

        try:
            for i in range(20):
                P_intank, ratio, mdot_total, mdot_LOx, mdot_IPA, mdot_choke = tank_state(P_guess, h_b) #guess P_tank of 400 psi
                if abs(P_guess - P_intank) < 10:
                    break
                P_guess = P_intank

        except OrificeChokes as e:
            t_choke = time
            min_choke_ratio = min(e.args[0], min_choke_ratio)
            print(f"Orifice choked at {time:.3f}s and ratio {e.args[0]:.2f}")
            break
        except ThroatCondenses:
            print(f"Throat condensed at {time:.3f}s and P_b {P_b/6894.757:.0f} psi")
            break
        except ValueError as e: #property/solver failure: stop this bottle instead of crashing the whole run
            print(f"Orifice solve failed at {time:.3f}s and P_b {P_b/6894.757:.0f} psi: {e}")
            break
        min_P_intank = min(P_intank, min_P_intank)

        mdots_LOx.append(mdot_LOx)
        mdots_IPA.append(mdot_IPA)
        mdots_total.append(mdot_total)
        min_choke_ratio = min(min_choke_ratio, ratio)
        min_P_intank = min(min_P_intank, P_intank)
        max_P_intank = max(max_P_intank, P_intank)
        Ps.append(P_b/6894.757)
        P_intanks.append(P_intank/6894.757)
        Ts.append(T_b) #temp of gn2 in bottle
        Tins.append(T_in)
        times.append(time)
        mdots_choke.append(mdot_choke)

        m_used += mdot_total * dt
        rho_b = rho_b0-(m_used+m_pres)/V_b #kg/m^3
        P_b = N2_DS(rho_b, S).p()
        
        time += dt
    if not mdots_total:
        return None
    T_bf = propssi("T", "S", S, "D", max(rho_b, rho_end), "nitrogen") #end of burn
    T_f = propssi("T", "S", S, "D", rho_end, "nitrogen") #end of expulsion
    T_margin = T_f - T_rating
    m_excess = rho_b * V_b - rho_end * V_b

    print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    print(f"V_b: {1000*V_b:.2f}L, fill mass: {rho_b0*V_b:.3f} kg, pressurization mass: {m_pres:.3f} kg\nP_b0: {P_b0/6894.757:.0f} psi, P_b1: {P_b1/6894.757:.0f} psi, T_b1: {T_b1:.0f} K, T_f: {T_f:.0f} K")
    print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    output = {
        "success":True,
        "P_margin":min(min_P_intank - P_ullage_min, P_ullage_max - max_P_intank), #so this becomes difference btwn P_intank and desired ullage pressure
        "rho_margin":rho_b-rho_end,
        "choke_ratio":min_choke_ratio,
        "t_choke":t_choke,
        "P_intanks":P_intanks,
        "T_margin":T_margin,
        "Ps":Ps,
        "mdots_LOx_max":max(mdots_LOx),
        "mdots_IPA_max":max(mdots_IPA),
        "mdots_total_max":max(mdots_total),
        "m_excess":m_excess,
        "fill_mass":rho_b0*V_b,
        "P_b1":P_b1,
        "T_b1":T_b1,
        "m_pres":m_pres,
        "T_f":T_f,
        "Tins":Tins,
        "Ts":Ts,
        "T_bf":T_bf,
        "mdots_total":mdots_total,
        "mdots_LOx":mdots_LOx,
        "mdots_IPA":mdots_IPA,
        "mdots_choke":mdots_choke,
        "time":time,
        "times":times,
        "t_choke":t_choke
        }
    return output


def burn_bottle(P_b0, V_b):
    time = 0 #s
    m_pres = 0.08 #kg, initial guess
    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin+P_loss, "S", S, "nitrogen") #kg/m^3
    T_end = propssi("T", "S", S, "D", rho_end, "nitrogen")
    m_avail = V_b*(rho_b0-rho_end)

    for i in range(50):
        m_pres = min(m_pres, 0.9999*m_avail)
        rho_b1 = rho_b0-m_pres/V_b
        u_b0 = N2_DS(rho_b0, S).umass()
        u_b1 = N2_DS(rho_b1, S).umass()
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres
        m_new = fill_mass(V_utotal, h_in) #can't account for CF_pp of LOx with this
        done = abs(m_pres - m_new) < 0.00001
        m_pres = m_new
        if done:
            break

    if m_pres > m_avail: #come back here to figure out margin + m_pres
        return {"rho_margin":(m_avail-m_pres)/V_b, "T_bf":propssi("T", "S", S, "D", rho_end, "nitrogen"), "m_pres":m_pres}

    m_used = 0 #kg
    rho_b = rho_b0-(m_used+m_pres)/V_b #need to recalculate these after m_pres calc
    rho_b1 = rho_b0-m_pres/V_b #these are for the start of burn btw

    while time <= depl_time: #come back here: check outputs with check valve and also line flow
        h_b = N2_DS(rho_b, S).hmass() #isentropic flow, specific H bottle = specific H tanks
        rho_in = N2_HP(h_b, P_t).rhomass() #rho GN2 going into tanks
        
        mdot_LOx = Vdot_LOx * CF * rho_in
        mdot_IPA = Vdot_IPA * rho_in
        mdot_total = mdot_LOx + mdot_IPA

        rho_b -= mdot_total*dt/V_b #kg/m^3
        if rho_b <= rho_end:
            left = depl_time - time
            return {
                "rho_margin":rho_b-rho_end-(mdot_total*left/V_b),
                "T_bf":T_end,
                "m_pres":m_pres
            }
        time += dt

    output = {
        "rho_margin":rho_b-rho_end,
        "T_bf":propssi("T", "S", S, "D", max(rho_end, rho_b), "nitrogen"),
        "m_pres":m_pres
        }
    return output


def endT(P_fill):
    S = propssi("S", "P", P_fill, "T", T_b0, "nitrogen")
    return propssi("T", "P", P_tmin+P_loss, "S", S, "nitrogen")

P_lo = P_tmin*1.1
P_hi = 5300*6894.757

#come back here, make part of loop
if (endT(P_lo) > (T_rating + T_tol)) & (endT(P_hi) < (T_rating + T_tol)): #if lower bound meets temp req and upper bound doesnt, solve for greatest upper bound
    P_fillmax = optimize.brentq(lambda P: endT(P) - (T_rating + T_tol), P_lo, P_hi, xtol = 0.1*6894.757)
else:
    P_fillmax = P_hi #if not, P_lo always < P_hi and P_lo basically always fulfills the temperature requirement so P_hi also fulfills temperature requirement within a realistic range (we're not going above 6 ksi)

def solve_fill(V_b, P_rating):
    hi = min(P_rating, P_fillmax) #something's sketchy? come back here
    margin_lo = use_bottle(P_lo, V_b)["P_margin"]
    margin_hi = use_bottle(hi, V_b)["P_margin"]
    if margin_hi < 0:
        return hi, "cooked" #kinda sketch, come back here
    elif margin_lo >= 0:
        return P_lo, "all good"
    P_fill = optimize.brentq(lambda P: use_bottle(P, V_b)["P_margin"], P_lo, hi)
    return P_fill, "ok"

def solve_fill1(V_b, P_rating):
    at_rating = burn_bottle(P_rating, V_b) #used twice below, only run it once
    if at_rating["rho_margin"] < 0:
        return P_rating, "cooked" #if bottle can't last depl time at max GN2
    elif burn_bottle(P_lo, V_b)["rho_margin"] >= 0:
        P_fill = P_lo #if bottle always makes depl time
    else: #find pressure that just makes depl time
        P_fill = optimize.brentq(lambda P: burn_bottle(P, V_b)["rho_margin"], P_lo, P_rating, xtol=P_XTOL)

    if at_rating["T_bf"] < T_rating:
        return P_rating, "cooked"
    if burn_bottle(P_fill, V_b)["T_bf"] < T_rating:
        P_fill = optimize.brentq(lambda P: burn_bottle(P, V_b)["T_bf"] - T_rating, P_fill, P_rating, xtol=P_XTOL)
    if endT(P_fill) >= T_rating:
        return P_fill, "ok"
    return P_fill, "ok at burn end"

def solve_fillloss(V_b, P_rating): #fill pressure for orifice UPSTREAM of the reg (orifice loss eats into reg inlet pressure)
    global P_loss
    P_loss = 0 #resetting P_loss for each bottle
    for i in range(6):
        P_fill, status = solve_fill1(V_b, P_rating)
        out = use_bottle(P_fill, V_b)
        if out is None or status == "cooked":
            break
        loss_new = (out["Ps"][-1] - out["P_inregs"][-1])*6894.757 #bottle pressure - P_inreg gives pressure loss across orifice at last timestep where it's worst
        if abs(loss_new - P_loss) < 6894.757: #1 psi error for loss
            break
        P_loss = loss_new
    return P_fill, status, out #out = use_bottle(P_fill) already, so it doesn't have to be rerun for plotting

def solve_fill_downstream(V_b, P_rating): #fill pressure for orifice DOWNSTREAM of the reg: nothing between bottle and reg inlet
    global P_loss
    P_loss = 0
    return solve_fill1(V_b, P_rating)

def evaluate_bottle(task): #everything for one bottle; runs in a worker process
    P_rating, V_b, partnum, OD, man, weight = task
    log = io.StringIO() #worker prints are collected and printed in order by the main process
    with contextlib.redirect_stdout(log):
        print("---------------------------------------------------------------------------------------")
        print(f"{man} {partnum}\n")
        print("Orifice upstream:")
        P_fill, status, out_up = solve_fillloss(V_b, P_rating)
        print("\nOrifice downstream:")
        P_fill_down, status_down = solve_fill_downstream(V_b, P_rating)
        out_down = use_bottle1(P_fill_down, V_b)
        print(f"\n")
    return {"V_b":V_b, "partnum":partnum, "OD":OD, "man":man, "weight":weight,
            "P_fill":P_fill, "status":status, "out_up":out_up,
            "P_fill_down":P_fill_down, "status_down":status_down, "out_down":out_down,
            "log":log.getvalue()}

def usebottle(P_fill, V_b, partnum, OD, man, status, weight, output):
    if output is None:
        print(f"{man} {partnum} skipped")
        return
    if (status != "cooked") & (output["t_choke"] is not None):
        status = f"{status}, chokes at {output["t_choke"]:.2f}s"
    x1.append(P_fill/6894.757293) #psi, come back here after solver
    y1.append(1000*V_b) #L
    mdots1.append((output["mdots_LOx_max"], output["mdots_IPA_max"], output["mdots_total_max"]))
    BPs.append(output["Ps"])
    excesses.append(output["m_excess"])
    fill_masses1.append(output["fill_mass"])
    press_masses1.append(output["m_pres"])
    P_b1s1.append(output["P_b1"]/6894.757293) #psi
    T_b1s1.append(output["T_b1"])
    T_fs1.append(output["T_f"])
    T_bfs.append(output["T_bf"])
    if (output["T_f"] >= (T_rating+T_tol)) & (output["time"] >= depl_time) & (output["choke_ratio"] > 1) & (output["P_margin"] > 0):
        Tcolors1.append("lightgreen")
    else:
        Tcolors1.append("salmon")
    all_mdots1.append(output["mdots_total"])
    expl_times.append(output["time"])
    part_nums.append(f"{man} " + partnum)
    COPV_ODs.append(OD)
    statuses.append(status)
    bweights.append(weight)
    choke_times.append(output["t_choke"])
    LOx_mdots1.append(output["mdots_LOx"])
    IPA_mdots1.append(output["mdots_IPA"])
    Tins1.append(output["Tins"])
    BTs.append(output["Ts"])
    times1.append(output["times"])
    choke_mdots1.append(output["mdots_choke"])

def usebottle1(P_fill, V_b, partnum, OD, man, status, weight, output):
    if output is None:
        print(f"{man} {partnum} skipped")
        return
    if (status != "cooked") & (output["t_choke"] is not None):
        status = f"{status}, chokes at {output["t_choke"]:.2f}s"
    x2.append(P_fill/6894.757293) #psi, come back here after solver
    y2.append(1000*V_b) #L
    mdots2.append((output["mdots_LOx_max"], output["mdots_IPA_max"], output["mdots_total_max"]))
    BPs2.append(output["Ps"])
    excesses2.append(output["m_excess"])
    fill_masses2.append(output["fill_mass"])
    press_masses2.append(output["m_pres"])
    P_b1s2.append(output["P_b1"]/6894.757293) #psi
    T_b1s2.append(output["T_b1"])
    T_fs2.append(output["T_f"])
    T_bfs2.append(output["T_bf"])
    if (output["T_f"] >= (T_rating+T_tol)) & (output["time"] >= depl_time) & (output["choke_ratio"] > 1) & (output["P_margin"] > 0):
        Tcolors2.append("lightgreen")
    else:
        Tcolors2.append("salmon")
    all_mdots2.append(output["mdots_total"])
    expl_times2.append(output["time"])
    part_nums2.append(f"{man} " + partnum)
    COPV_ODs2.append(OD)
    statuses2.append(status)
    bweights2.append(weight)
    choke_times2.append(output["t_choke"])
    LOx_mdots2.append(output["mdots_LOx"])
    IPA_mdots2.append(output["mdots_IPA"])
    Tins2.append(output["Tins"])
    BTs2.append(output["Ts"])
    times2.append(output["times"])
    choke_mdots2.append(output["mdots_choke"])

#Main run. Guarded so the worker processes (which import this file) don't rerun the sweep, plots and CSVs
if __name__ == "__main__":
    plt.close()
    fig, ax = plt.subplots(2, 2, figsize=(8,8), layout = "constrained")
    ax[0,0].set_xlabel("Fill pressure (psi)")
    ax[0,0].set_ylabel("Bottle volume (L)")
    ax[0,0].grid(True)
    ax[1,0].set_ylabel("End temperature (K)")
    ax[1,0].set_xlabel("Bottle volume (L)")
    ax[1,0].grid(True)
    ax[0,1].set_xlabel("Rated pressure (psi)")
    ax[0,1].set_ylabel("Bottle volume (L)")
    ax[0,1].grid(True)
    ax[1,1].set_xlabel("Rated pressure (psi)")
    ax[1,1].set_ylabel("Bottle volume (L)")

    V_candidates = readCSV(filename, "V (L)")
    P_candidates = readCSV(filename, "Rated pressure (psi)")
    partnums = readCSV(filename, "Part #")
    ODs = readCSV(filename, "Diameter (in)")
    mans = readCSV(filename, "Company")
    weights = readCSV(filename, "Weight (kg)")

    tasks = []
    for P_psi, V_L, partnum, OD, man, weight in zip(P_candidates, V_candidates, partnums, ODs, mans, weights, strict=True):
        if float(V_L) > V_max:
            continue
        V_b = float(V_L)/1000
        P_psi = P_psi.replace(",", "")
        P_psi = float(P_psi)*6894.757
        tasks.append((P_psi, V_b, partnum, float(OD), man, weight))

    #bottles are independent, so each one is solved in its own worker process (one per CPU core, minus one for the sweep)
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as pool:
        pending = pool.map_async(evaluate_bottle, tasks)

        #ideal-bottle sweep runs in this process while the workers do the bottles
        P_bsweep = 2000*6894.757293 #Pa, full fill pressure
        while P_bsweep <= 4700*6894.757293: #4700 psi, sweep end
            tryP_b0(P_bsweep)
            P_bsweep += P_step*6894.757293

        results = pending.get()

    for r in results: #same order as the CSV
        print(r["log"], end="")
        usebottle(r["P_fill"], r["V_b"], r["partnum"], r["OD"], r["man"], r["status"], r["weight"], r["out_up"])
        usebottle1(r["P_fill_down"], r["V_b"], r["partnum"], r["OD"], r["man"], r["status_down"], r["weight"], r["out_down"])


    points = sorted(zip(x, y, colors, mdots, fill_masses, P_b1s, T_b1s, press_masses, T_fs, Tcolors, all_mdots, m_useds), key = lambda p: p[0])
    sorted_x, sorted_y, sorted_colors, sorted_mdots, sorted_fill_masses, sorted_Pb1s, sorted_Tb1s, sorted_press_masses, sorted_Tfs, sorted_Tcolors, sorted_all_mdots, sorted_m_useds = zip(*points)
    ax[0,0].plot(sorted_x, sorted_y, color = "gray", linestyle = "-")
    ax[1,0].plot(sorted_y, sorted_Tfs, color = "gray", linestyle = "-")
    points1 = sorted(zip(x1, y1, mdots1, fill_masses1, P_b1s1, T_b1s1, press_masses1, T_fs1, Tcolors1, all_mdots1, excesses, expl_times, part_nums, COPV_ODs, statuses, bweights, BPs, T_bfs, choke_times), key = lambda p: p[12])
    sorted_x1, sorted_y1, sorted_mdots1, sorted_fill_masses1, sorted_P_b1s1, sorted_Tb1s1, sorted_press_masses1, sorted_T_fs1, sorted_Tcolors1, sorted_all_mdots1, sorted_excesses, sorted_expl_times, sorted_part_nums, sorted_COPV_ODs, sorted_statuses, sorted_weights, sorted_BPs, sorted_T_bfs, sorted_choke_times = zip(*points1, strict = True)
    if x2:
        points2 = sorted(zip(x2, y2, mdots2, fill_masses2, P_b1s2, T_b1s2, press_masses2, T_fs2, Tcolors2, all_mdots2, excesses2, expl_times2, part_nums2, COPV_ODs2, statuses2, bweights2, BPs2, T_bfs2, choke_times2), key = lambda p: p[12])
        sorted_x2, sorted_y2, sorted_mdots2, sorted_fill_masses2, sorted_P_b1s2, sorted_Tb1s2, sorted_press_masses2, sorted_T_fs2, sorted_Tcolors2, sorted_all_mdots2, sorted_excesses2, sorted_expl_times2, sorted_part_nums2, sorted_COPV_ODs2, sorted_statuses2, sorted_weights2, sorted_BPs2, sorted_T_bfs2, sorted_choke_times2 = zip(*points2, strict = True)


    csvname = f"{T_amb}K_upstream{D_orifice}in.csv"
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    os.makedirs(output_dir, exist_ok=True)
    csvname = os.path.join(output_dir, csvname)
    with open(csvname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Part_num", "OD_in", "P_b0_psi", "P_b1_psi", "V_b_L", "T_b1_K", "T_bf_K", "T_f_K", "fill_mass_kg", "press_mass_kg", "mdot_LOx_max_kg_s", "mdot_IPA_max_kg_s",
                         "mdot_total_max_kg_s", "all_mdots_kg_s", "everything_ok", "expl_time_s", "weight_kg", "Ps_psi", "choke_time_s"]) #T_sk?
        for px, vy, col, (lox, ipa, tot), fm, P1, T1, pm, T_f, amd, ext, pn, od, kg, bps, T_bf, CT in zip(sorted_x1, sorted_y1, sorted_Tcolors1, sorted_mdots1, sorted_fill_masses1, sorted_P_b1s1, sorted_Tb1s1, sorted_press_masses1, sorted_T_fs1, sorted_all_mdots1, sorted_expl_times, sorted_part_nums, sorted_COPV_ODs, sorted_weights, sorted_BPs, sorted_T_bfs, sorted_choke_times):
            writer.writerow([pn, od, px, P1, vy, T1, T_bf, T_f, fm, pm, lox, ipa, tot, amd, col == "lightgreen", ext, kg, bps, CT])

    if x2:
        csvname1 = f"{T_amb}K_downstream{D_orifice}in.csv"
        output_dir1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
        os.makedirs(output_dir1, exist_ok=True)
        csvname1 = os.path.join(output_dir1, csvname1)
        with open(csvname1, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Part_num", "OD_in", "P_b0_psi", "P_b1_psi", "V_b_L", "T_b1_K", "T_bf_K", "T_f_K", "fill_mass_kg", "press_mass_kg", "mdot_LOx_max_kg_s", "mdot_IPA_max_kg_s",
                            "mdot_total_max_kg_s", "all_mdots_kg_s", "everything_ok", "expl_time_s", "weight_kg", "Ps_psi", "choke_time_s"])
            for px, vy, col, (lox, ipa, tot), fm, P1, T1, pm, T_f, amd, ext, pn, od, kg, bps, T_bf, CT in zip(sorted_x2, sorted_y2, sorted_Tcolors2, sorted_mdots2, sorted_fill_masses2, sorted_P_b1s2, sorted_Tb1s2, sorted_press_masses2, sorted_T_fs2, sorted_all_mdots2, sorted_expl_times2, sorted_part_nums2, sorted_COPV_ODs2, sorted_weights2, sorted_BPs2, sorted_T_bfs2, sorted_choke_times2):
                writer.writerow([pn, od, px, P1, vy, T1, T_bf, T_f, fm, pm, lox, ipa, tot, amd, col == "lightgreen", ext, kg, bps, CT])

    def find_part(names, maker = "ALSAFE"): #index of the first bottle from this maker (lists are in CSV order minus skipped bottles)
        for k, name in enumerate(names):
            if name.upper().startswith(maker):
                return k
        print(f"No {maker} bottle in results, transient CSV skipped")
        return None

    k_up = find_part(part_nums)
    k_down = find_part(part_nums2) if x2 else None

    if k_up is not None:
        csvname2 = f"{T_amb}K_Alsafe_upstream{D_orifice}in_transient.csv"
        output_dir2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
        os.makedirs(output_dir2, exist_ok=True)
        csvname2 = os.path.join(output_dir2, csvname2)
        with open(csvname2, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["time_s", "P_b_psi", "T_b_K", "T_in_K", "mdot_total_kg_s", "mdot_choke_kg_s", "mdot_LOx_kg_s", "mdot_IPA_kg_s"])
            for t, Pb, Tb, Tin, tot, lox, ipa, choke_md in zip(times1[k_up], BPs[k_up], BTs[k_up], Tins1[k_up], all_mdots1[k_up], LOx_mdots1[k_up], IPA_mdots1[k_up], choke_mdots1[k_up]):
                writer.writerow([t, Pb, Tb, Tin, tot, choke_md, lox, ipa])

    if k_down is not None:
        csvname3 = f"{T_amb}K_Alsafe_downstream{D_orifice}in_transient.csv"
        output_dir3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
        os.makedirs(output_dir3, exist_ok=True)
        csvname3 = os.path.join(output_dir3, csvname3)
        with open(csvname3, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["time_s", "P_b_psi", "T_b_K", "T_in_K", "mdot_total_kg_s", "mdot_choke_kg_s", "mdot_LOx_kg_s", "mdot_IPA_kg_s"])
            for t, Pb, Tb, Tin, tot, lox, ipa, choke_md in zip(times2[k_down], BPs2[k_down], BTs2[k_down], Tins2[k_down], all_mdots2[k_down], LOx_mdots2[k_down], IPA_mdots2[k_down], choke_mdots2[k_down]):
                writer.writerow([t, Pb, Tb, Tin, tot, choke_md, lox, ipa])

    sc = ax[0,0].scatter(sorted_x, sorted_y, c = sorted_colors, zorder = 3)
    sc1 = ax[1,0].scatter(sorted_y, sorted_Tfs, c= sorted_Tcolors, zorder = 3)

    m_total = [mu + mp for mu, mp in zip(sorted_m_useds, sorted_press_masses)]
    P_grid = np.linspace(min(sorted_x)*0.95, max(sorted_x)*1.25, 80)
    V_grid = np.linspace(min(sorted_y)*0.95, V_max, 80)
    T_grid = np.full((len(V_grid), len(P_grid)), np.nan)
    for i, P in enumerate(P_grid):
        S = propssi("S", "T", T_b0, "P", 6894.757*P, "nitrogen")
        rho_0 = propssi("D", "T", T_b0, "P", 6894.757*P, "nitrogen")
        rho_min = propssi("D", "S", S, "P", P_tmin, "nitrogen")
        m_tot = np.interp(P, sorted_x, m_total)
        for j, V in enumerate(V_grid):
            rho_end = max(rho_0 - m_tot/(V/1000), rho_min)
            T_grid[j,i] = N2_DS(rho_end, S).T()

    vmin, vmax = np.nanmin(T_grid), np.nanmax(T_grid)
    levels = np.linspace(vmin, vmax, 20)

    cf = ax[0,1].contourf(P_grid, V_grid, T_grid, levels=levels, cmap="viridis", zorder=0)
    ax[0,1].contour(P_grid, V_grid, T_grid, levels=[T_rating], colors="red", linewidths=1.5, zorder=2)
    fig.colorbar(cf, ax=ax[0,1], label = "T_f (K)")
    ax[0,1].plot(sorted_x, sorted_y, color = "gray", linestyle = "-", zorder = 3)
    sc2 = ax[0,1].scatter(sorted_x1, sorted_y1, c=sorted_Tcolors1, edgecolor = "white", linewidth = 1, s = 60, zorder = 4)
    ax[0,1].set_ylim(V_grid[0], V_max)

    if x2:
        cf1 = ax[1,1].contourf(P_grid, V_grid, T_grid, levels=levels, cmap="viridis", zorder=0)
        ax[1,1].contour(P_grid, V_grid, T_grid, levels=[T_rating], colors="red", linewidths=1.5, zorder=2)
        fig.colorbar(cf, ax=ax[1,1], label = "T_f (K)")
        ax[1,1].plot(sorted_x, sorted_y, color = "gray", linestyle = "-", zorder = 3)
        sc3 = ax[1,1].scatter(sorted_x2, sorted_y2, c=sorted_Tcolors2, edgecolor = "white", linewidth = 1, s = 60, zorder = 4)
        ax[1,1].set_ylim(V_grid[0], V_max)

    fig.suptitle(f"GN2 COPV blowdown sizing {T_amb}K T_amb", fontsize = 13, fontweight = "bold")
    ax[0,0].set_title(f"Bottle volume vs fill pressure for 10 second burn\nblue = P_tmin reached, red = P_tmin not reached", fontsize = 10)
    ax[1,0].set_title(f"Burn end temperature vs bottle volume\nblue = above COPV rated temp, red = below COPV rated temp", fontsize = 10)
    ax[0,1].set_title(f"COPV options on temperature map orifice upstream\nred line = {T_rating}K, gray line = 10s burn, green = ok, salmon = not ok", fontsize = 10)
    if x2:
        ax[1,1].set_title(f"COPV options on temperature map orifice downstream")

    # params = (
    #     f"Run parameters\n"
    #     f"T_amb         {T_amb:.2f}K\n"
    #     f"T_b0          {T_b0:.2f}K\n"
    #     f"depl_time     {depl_time:.2f}s\n"
    #     f"P_t           {P_t/6894.757:.0f}psi\n"
    #     f"P_tmin        {P_tmin/6894.757:.0f}psi\n"
    #     f"LOx CF        {CF}\n" #come back here
    #     f"V_utotal      {V_utotal*1000:.2f}L\n"
    # )
    # ax[1,1].text(0.02, 0.98, params, transform=ax[1,1].transAxes, va="top", ha="left", family="monospace", fontsize=10)

    cursor1 = mplcursors.cursor(sc, hover=mplcursors.HoverMode.Transient)
    @cursor1.connect("add")
    def on_add(sel):
        LOx, IPA, total = sorted_mdots[sel.index]
        P = sorted_x[sel.index]
        P1 = sorted_Pb1s[sel.index]
        V = sorted_y[sel.index]
        T1 = sorted_Tb1s[sel.index]
        TF = sorted_Tfs[sel.index]
        FM = sorted_fill_masses[sel.index]
        PM = sorted_press_masses[sel.index]
        sel.annotation.set_text(f"{V:.2f} L tank\n{FM:.3f} kg fill mass\n{PM:.3f} kg press mass\n{P:.0f} psi full\n{P1:.0f} psi burn start\n{T1:.0f}K burn start\n{TF:.0f}K burn end\n{LOx:.3f} kg/s LOx\n{IPA:.3f} kg/s IPA\n{total:.3f} kg/s total")
        sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)

    cursor2 = mplcursors.cursor(sc1, hover=mplcursors.HoverMode.Transient)
    @cursor2.connect("add")
    def on_add(sel):
        V = sorted_y[sel.index]
        TF = sorted_Tfs[sel.index]
        sel.annotation.set_text(f"{V:.2f} L tank\n{TF:.3f}K final temp")
        sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)

    cursor3 = mplcursors.cursor(sc2, hover=mplcursors.HoverMode.Transient)
    @cursor3.connect("add")
    def on_add(sel):
        name = sorted_part_nums[sel.index]
        OD = sorted_COPV_ODs[sel.index]
        LOx, IPA, total = sorted_mdots1[sel.index]
        P = sorted_x1[sel.index]
        V = sorted_y1[sel.index]
        TF = sorted_T_fs1[sel.index]
        EX = sorted_excesses[sel.index]
        T = sorted_expl_times[sel.index]
        stat = sorted_statuses[sel.index]
        Tbf = sorted_T_bfs[sel.index]
        sel.annotation.set_text(f"{name}\n{OD:.3f} in OD\n{V:.2f} L tank\nstatus: {stat}\n{P:.0f} psi full\n{Tbf:.0f}K burn end\n{TF:.0f}K expulsion end\n{LOx:.3f} kg/s LOx\n{IPA:.3f} kg/s IPA\n{total:.3f} kg/s total\n{T:.3f}s expulsion time\n{EX:.3f}kg excess mass")
        sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)

    if x2:
        cursor4 = mplcursors.cursor(sc3, hover=mplcursors.HoverMode.Transient)
        @cursor4.connect("add")
        def on_add(sel):
            name = sorted_part_nums2[sel.index]
            OD = sorted_COPV_ODs2[sel.index]
            LOx, IPA, total = sorted_mdots2[sel.index]
            P = sorted_x2[sel.index]
            V = sorted_y2[sel.index]
            TF = sorted_T_fs2[sel.index]
            EX = sorted_excesses2[sel.index]
            T = sorted_expl_times2[sel.index]
            stat = sorted_statuses2[sel.index]
            Tbf = sorted_T_bfs2[sel.index]
            sel.annotation.set_text(f"{name}\n{OD:.3f} in OD\n{V:.2f} L tank\nstatus: {stat}\n{P:.0f} psi full\n{Tbf:.0f}K burn end\n{TF:.0f}K expulsion end\n{LOx:.3f} kg/s LOx\n{IPA:.3f} kg/s IPA\n{total:.3f} kg/s total\n{T:.3f}s expulsion time\n{EX:.3f}kg excess mass")
            sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)

    plt.show()
