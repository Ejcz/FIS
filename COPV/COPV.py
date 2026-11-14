from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
import csv
import os

#come back here: debug the end temperature or contour being weird
plt.close()
x = []
y = []
x1 = []
y1 = []
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
Tcolors = []
Tcolors1 = []
all_mdots = []
all_mdots1 = []
excesses = []
expl_times = []
part_nums = []
COPV_ODs = []
statuses = []

#constants/variables
Vdot_LOx = 0.001111306 #m^3/s, nominal
mdot_LOx = 1.268 #kg/s
CF = 2
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


def fill_mass(V_u, h_in):
    m_0 = V_u * propssi("D", "P", 101325, "T", T_amb, "nitrogen")
    u_u0 = propssi("UMASS", "T", T_amb, "P", 101325, "air")
    def f(dm):
        m_ppres = m_0 + dm
        u_f = (m_0*u_u0 + dm*h_in)/m_ppres
        return propssi("P", "D", m_ppres/V_u, "UMASS", u_f, "nitrogen") - P_t
    dm_hi = 2*propssi("D", "P", P_t, "T", T_amb, "nitrogen")*V_u #just a safe upper bound, can mostly ignore
    return optimize.brentq(f, 0, dm_hi)

#burn start to end
def size_bottle(P_b0, m_usedguess = 2.2): #P_b0 = before prepressurization
    #constants
    m_used = 0 #kg
    m_pres = 0.2641458904 #kg, initial guess for pressurization mass
    error = 1 #placeholder

    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3

    while error > V_error:
        time = 0 #s
        m_used = 0 #kg
        mdots_LOx = [] #lists are for GN2 btw
        mdots_IPA = []
        mdots_total = []
        
        V_b = m_usedguess/(rho_b0-rho_end) #m^3, guess bottle volume based on m_used
        rho_b1 = rho_b0-m_pres/V_b #density after pressurization, before burn
        u_b0 = propssi("UMASS", "S", S, "D", rho_b0, "nitrogen") #specfic internal energy before pressurization
        u_b1 = propssi("UMASS", "S", S, "D", rho_b1, "nitrogen") #specfic internal energy after pressurization
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
            h_b = propssi("HMASS", "S", S, "D", rho_b, "nitrogen") #isentropic flow, specific H bottle = specific H tanks
            rho_in = propssi("D", "P", P_t, "HMASS", h_b, "nitrogen") #rho GN2 going into tanks
            
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
def use_bottle(P_b0, V_b):
    time = 0 #s
    m_used = 0 #kg
    m_pres = 0.08 #kg, initial guess
    mdots_LOx = []
    mdots_IPA = []
    mdots_total = []
    S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
    rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
    rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3

    # m_presLOx = fill_mass(V_uLOx, h_in) #no line volume COME BACK HERE
    # m_presIPA = fill_mass(V_uIPA, h_in)
    # m_pres = m_presLOx * 1.5 + m_presIPA #add a check on this? come back here, kinda sketch
    m_avail = V_b*(rho_b0-rho_end)
    for i in range(50):
        if m_avail < m_pres:
            return {"rho_margin":(m_avail-m_pres)/V_b, "T_margin":propssi("T", "S", S, "P", P_tmin, "nitrogen")-T_rating}
        rho_b1 = rho_b0-m_pres/V_b
        u_b0 = propssi("UMASS", "S", S, "D", rho_b0, "nitrogen")
        u_b1 = propssi("UMASS", "S", S, "D", rho_b1, "nitrogen")
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres
        m_new = fill_mass(V_utotal, h_in) #can't account for CF_pp of LOx with this
        done = abs(m_pres - m_new) < 0.00001
        m_pres = m_new
        if done:
            break

    if m_pres > m_avail:
        return {"rho_margin":(m_avail-m_pres)/V_b, "T_margin":propssi("T", "S", S, "P", P_b0, "nitrogen")-T_rating}

    m_used = 0 #kg
    rho_b = rho_b0-(m_used+m_pres)/V_b #need to recalculate these after m_pres calc
    rho_b1 = rho_b0-m_pres/V_b #these are for the start of burn btw
    P_b1 = propssi("P", "D", rho_b1, "S", S, "nitrogen")
    T_b1 = propssi("T", "P", P_b1, "D", rho_b1, "nitrogen")

    P_b = P_b1
    while time <= depl_time: #come back here: check outputs with check valve and also line flow
        if P_b < P_tmin:
            break
        h_b = propssi("HMASS", "S", S, "D", rho_b, "nitrogen") #isentropic flow, specific H bottle = specific H tanks
        T_in = propssi("T", "P", P_t, "Hmass", h_b, "nitrogen") #T going into prop tanks. check this at the end? come back here
        rho_in = propssi("D", "P", P_t, "T", T_in, "nitrogen") #rho GN2 going into tanks
        
        mdot_LOx = Vdot_LOx * CF * rho_in
        mdot_IPA = Vdot_IPA * rho_in
        mdot_total = mdot_LOx + mdot_IPA
        mdots_LOx.append(mdot_LOx)
        mdots_IPA.append(mdot_IPA)
        mdots_total.append(mdot_total)

        m_used += mdot_total * dt
        rho_b = rho_b0-(m_used+m_pres)/V_b #kg/m^3
        P_b = propssi("P", "S", S, "D", rho_b, "nitrogen")
        time += dt
    T_f = propssi("T", "S", S, "D", rho_b, "nitrogen")
    T_margin = T_f - T_rating
    m_excess = rho_b * V_b - rho_end * V_b

    print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    print(f"V_b: {1000*V_b:.2f}L, fill mass: {rho_b0*V_b:.3f} kg, pressurization mass: {m_pres:.3f} kg\nP_b0: {P_b0/6894.757:.0f} psi, P_b1: {P_b1/6894.757:.0f} psi, T_b1: {T_b1:.0f} K, T_f: {T_f:.0f}")
    print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    output = {
        "rho_margin":rho_b-rho_end,
        "T_margin":T_margin,
        "mdots_LOx_max":max(mdots_LOx),
        "mdots_IPA_max":max(mdots_IPA),
        "mdots_total_max":max(mdots_total),
        "m_excess":m_excess,
        "fill_mass":rho_b0*V_b,
        "P_b1":P_b1,
        "T_b1":T_b1,
        "m_pres":m_pres,
        "T_f":T_f,
        "mdots_total":mdots_total,
        "time":time
        }
    return output

def endT(P_fill):
    S = propssi("S", "P", P_fill, "T", T_b0, "nitrogen")
    return propssi("T", "P", P_tmin, "S", S, "nitrogen")

P_lo = P_tmin*1.1
P_hi = 5300*6894.757

if (endT(P_lo) > (T_rating + T_tol)) & (endT(P_hi) < (T_rating + T_tol)): #if lower bound meets temp req and upper bound doesnt, solve for greatest upper bound
    P_fillmax = optimize.brentq(lambda P: endT(P) - (T_rating + T_tol), P_lo, P_hi, xtol = 0.1*6894.757)
else:
    P_fillmax = P_hi #if not, P_lo always < P_hi and P_lo basically always fulfills the temperature requirement so P_hi also fulfills temperature requirement within a realistic range (we're not going above 6 ksi)

def solve_fill(V_b, P_rating):
    hi = min(P_rating, P_fillmax) #something's sketchy? come back here
    margin_lo = use_bottle(P_lo, V_b)["rho_margin"]
    margin_hi = use_bottle(hi, V_b)["rho_margin"]
    if margin_hi < 0:
        return hi, "cooked" #kinda sketch, come back here
    elif margin_lo > 0:
        return P_lo, "all good"
    P_fill = optimize.brentq(lambda P: use_bottle(P, V_b)["rho_margin"], P_lo, hi)
    return P_fill, "ok"

def solve_fill1(V_b, P_rating):
    hi = P_rating #something's sketchy? come back here
    margin_lo = use_bottle(P_lo, V_b)["T_margin"]
    margin_hi = use_bottle(hi, V_b)["T_margin"]
    if margin_hi >= 0:
        return hi, "cooked" #kinda sketch, come back here
    elif margin_lo < 0:
        return P_lo, "all good"
    P_fill = optimize.brentq(lambda P: use_bottle(P, V_b)["T_margin"], P_lo, hi)
    return P_fill, "ok"

def usebottle(P_fill, V_b, partnum, OD, man):
    output = use_bottle(P_fill, V_b)
    x1.append(P_fill/6894.757293) #psi, come back here after solver
    y1.append(1000*V_b) #L
    mdots1.append((output["mdots_LOx_max"], output["mdots_IPA_max"], output["mdots_total_max"]))
    excesses.append(output["m_excess"])
    fill_masses1.append(output["fill_mass"])
    press_masses1.append(output["m_pres"])
    P_b1s1.append(output["P_b1"]/6894.757293) #psi
    T_b1s1.append(output["T_b1"])
    T_fs1.append(output["T_f"])
    if (output["T_f"] >= (T_rating+T_tol)) & (output["time"] >= depl_time):
        Tcolors1.append("lightgreen")
    else:
        Tcolors1.append("salmon")
    all_mdots1.append(output["mdots_total"])
    expl_times.append(output["time"])
    part_nums.append(f"{man} " + partnum)
    COPV_ODs.append(OD)


P_bsweep = 2000*6894.757293 #Pa, full fill pressure
while P_bsweep <= 4700*6894.757293: #4700 psi, sweep end
    tryP_b0(P_bsweep)
    P_bsweep += P_step*6894.757293

V_candidates = readCSV(filename, "V (L)")
P_candidates = readCSV(filename, "Rated pressure (psi)")
partnums = readCSV(filename, "Part #")
ODs = readCSV(filename, "Diameter (in)")
mans = readCSV(filename, "Company")

for P_psi, V_L, partnum, OD, man in zip(P_candidates, V_candidates, partnums, ODs, mans, strict=True):
    if float(V_L) > V_max:
        continue
    V_b = float(V_L)/1000
    P_psi = P_psi.replace(",", "")
    P_psi = float(P_psi)*6894.757
    print("---------------------------------------------------------------------------------------")
    print(f"{man} {partnum}\n")
    P_fill, status = solve_fill1(V_b, P_psi) #adjust: solver or solver 1
    print(f"\n")
    statuses.append(status)
    usebottle(P_fill, V_b, partnum, float(OD), man)


points = sorted(zip(x, y, colors, mdots, fill_masses, P_b1s, T_b1s, press_masses, T_fs, Tcolors, all_mdots, m_useds), key = lambda p: p[0])
sorted_x, sorted_y, sorted_colors, sorted_mdots, sorted_fill_masses, sorted_Pb1s, sorted_Tb1s, sorted_press_masses, sorted_Tfs, sorted_Tcolors, sorted_all_mdots, sorted_m_useds = zip(*points)
ax[0,0].plot(sorted_x, sorted_y, color = "gray", linestyle = "-")
ax[1,0].plot(sorted_y, sorted_Tfs, color = "gray", linestyle = "-")
points1 = sorted(zip(x1, y1, mdots1, fill_masses1, P_b1s1, T_b1s1, press_masses1, T_fs1, Tcolors1, all_mdots1, excesses, expl_times, part_nums, COPV_ODs, statuses), key = lambda p: p[12])
sorted_x1, sorted_y1, sorted_mdots1, sorted_fill_masses1, sorted_P_b1s1, sorted_Tb1s1, sorted_press_masses1, sorted_T_fs1, sorted_Tcolors1, sorted_all_mdots1, sorted_excesses, sorted_expl_times, sorted_part_nums, sorted_COPV_ODs, sorted_statuses = zip(*points1)

csvname = f"{T_amb}K_COPVs_transient.csv"
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
os.makedirs(output_dir, exist_ok=True)
csvname = os.path.join(output_dir, csvname)
with open(csvname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Part_num", "OD_in", "P_b0_psi", "P_b1_psi", "V_b_L", "T_b1_K", "T_f_K", "fill_mass_kg", "press_mass_kg", "mdot_LOx_max_kg_s", "mdot_IPA_max_kg_s",
                     "mdot_total_max_kg_s", "all_mdots_kg_s", "T_f_ok", "expl_time_s"])
    for px, vy, col, (lox, ipa, tot), fm, P1, T1, pm, T_f, amd, ext, pn, od in zip(sorted_x1, sorted_y1, sorted_Tcolors1, sorted_mdots1, sorted_fill_masses1, sorted_P_b1s1, sorted_Tb1s1, sorted_press_masses1, sorted_T_fs1, sorted_all_mdots1, sorted_expl_times, sorted_part_nums, sorted_COPV_ODs):
        writer.writerow([pn, od, px, P1, vy, T1, T_f, fm, pm, lox, ipa, tot, amd, col == "lightgreen", ext])

sc = ax[0,0].scatter(sorted_x, sorted_y, c = sorted_colors, zorder = 3)
sc1 = ax[1,0].scatter(sorted_y, sorted_Tfs, c= sorted_Tcolors, zorder = 3)

m_total = [mu + mp for mu, mp in zip(sorted_m_useds, sorted_press_masses)]
P_grid = np.linspace(min(sorted_x)*0.95, max(sorted_x)*1.05, 80)
V_grid = np.linspace(min(sorted_y)*0.95, V_max, 80)
T_grid = np.full((len(V_grid), len(P_grid)), np.nan)
for i, P in enumerate(P_grid):
    S = propssi("S", "T", T_b0, "P", 6894.757*P, "nitrogen")
    rho_0 = propssi("D", "T", T_b0, "P", 6894.757*P, "nitrogen")
    rho_min = propssi("D", "S", S, "P", P_tmin, "nitrogen")
    m_tot = np.interp(P, sorted_x, m_total)
    for j, V in enumerate(V_grid):
        rho_end = max(rho_0 - m_tot/(V/1000), rho_min)
        T_grid[j,i] = propssi("T", "D", rho_end, "S", S, "nitrogen")

vmin, vmax = np.nanmin(T_grid), np.nanmax(T_grid)
levels = np.linspace(vmin, vmax, 20)

cf = ax[0,1].contourf(P_grid, V_grid, T_grid, levels=levels, cmap="viridis", zorder=0)
ax[0,1].contour(P_grid, V_grid, T_grid, levels=[T_rating], colors="red", linewidths=1.5, zorder=2)
fig.colorbar(cf, ax=ax[0,1], label = "T_f (K)")
ax[0,1].plot(sorted_x, sorted_y, color = "gray", linestyle = "-", zorder = 3)
sc2 = ax[0,1].scatter(sorted_x1, sorted_y1, c=sorted_Tcolors1, edgecolor = "white", linewidth = 1, s = 60, zorder = 4)
ax[0,1].set_ylim(V_grid[0], V_max)

fig.suptitle(f"GN2 COPV blowdown sizing {T_amb}K T_amb", fontsize = 13, fontweight = "bold")
ax[0,0].set_title(f"Bottle volume vs fill pressure for 10 second burn\nblue = P_tmin reached, red = P_tmin not reached", fontsize = 10)
ax[1,0].set_title(f"Burn end temperature vs bottle volume\nblue = above COPV rated temp, red = below COPV rated temp", fontsize = 10)
ax[0,1].set_title(f"COPV options on temperature map\nred line = {T_rating}K, gray line = 10s burn, green = ok, salmon = not ok", fontsize = 10)
ax[1,1].axis("off")

params = (
    f"Run parameters\n"
    f"T_amb         {T_amb:.2f}K\n"
    f"T_b0          {T_b0:.2f}K\n"
    f"depl_time     {depl_time:.2f}s\n"
    f"P_t           {P_t/6894.757:.0f}psi\n"
    f"P_tmin        {P_tmin/6894.757:.0f}psi\n"
    f"LOx CF        {CF}\n" #come back here
    f"V_utotal      {V_utotal*1000:.2f}L\n"
)

cooked_bs = (
    f"cooked:\n"
)
ok_bs = (
    f"ok:\n"
)
allgood_bs = (
    f"all good:\n"
)

cookeds = []
oks = []
allgoods = []
for s in range(len(sorted_statuses)):
    if sorted_statuses[s] == "cooked":
        cookeds.append(sorted_part_nums[s])
    elif sorted_statuses[s] == "allgood":
        allgoods.append(sorted_part_nums[s])
    elif sorted_statuses[s] == "ok":
        oks.append(sorted_part_nums[s])

cookedi = 0
oki = 0
for c in cookeds:
    cooked_bs += f"{c}\n"
for o in oks:
    ok_bs += f"{o}\n"
for a in allgoods:
    allgood_bs += f"{a}\n"


ax[1,1].text(0.02, 0.98, params, transform=ax[1,1].transAxes, va="top", ha="left", family="monospace", fontsize=10)
#ax[1,1].text(0.6, 0.98, cooked_bs, transform=ax[1,1].transAxes, va="top", ha="left", family="monospace", fontsize=10)
#ax[1,1].text(0.02, 0.30, ok_bs, transform=ax[1,1].transAxes, va="top", ha="left", family="monospace", fontsize=10)
#ax[1,1].text(0.02, 0.02, allgood_bs, transform=ax[1,1].transAxes, va="top", ha="left", family="monospace", fontsize=10)

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
    sel.annotation.set_text(f"{name}\n{OD:.3f} in OD\n{V:.2f} L tank\nstatus: {stat}\n{P:.0f} psi full\n{TF:.0f}K expulsion end\n{LOx:.3f} kg/s LOx\n{IPA:.3f} kg/s IPA\n{total:.3f} kg/s total\n{T:.3f}s expulsion time\n{EX:.3f}kg excess mass")
    sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)


plt.show()
