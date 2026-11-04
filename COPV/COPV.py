from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import matplotlib.pyplot as plt
import mplcursors
import csv

plt.close()
x = []
y = []
fig, ax = plt.subplots(2, figsize=(8,8))
ax[0].set_xlabel("Fill pressure (psi)")
ax[0].set_ylabel("Bottle volume (L)")
ax[0].grid(True)
ax[1].set_ylabel("End temperature (K)")
ax[1].set_xlabel("Bottle volume (L)")
ax[1].grid(True)
colors = []
mdots = []
fill_masses = []
T_b1s = []
P_b1s = []
press_masses = []
T_fs = []
Tcolors = []
all_mdots = []

def size_bottle(P_b0): #P_b0 = before prepressurization
    #constants
    CF = 2
    Vdot_LOx = 0.00115512708 #m^3/s, max
    mdot_LOx = 1.268 #kg/s
    m_LOx = 12.677951 #kg
    Vdot_IPA = 0.00114452214 #m^3/s, max
    mdot_IPA = 0.932 #kg/s
    m_IPA = 9.3 #kg
    P_t = 3548704 #Pa, postreg P
    P_tmin = 6894757.293 #Pa, min P into reg
    T_amb = 316 #K
    m_used = 0 #kg
    #P_b0 = 21766190.54 #Pa
    m_pres = 0.2641458904 #kg, may depend weakly on V_b
    V_uLOx = 0.00340 #m^3
    V_uIPA = 0.00307 #m^3
    V_utotal = 0.00699 #m^3
    T_b0 = 316 #K
    DL_ar = 0.2924845212 #diameter/length aspect ratio, will be a slight overestimate bc not enough info to find ellipsoidal surface area
    SA_Vbig = 14.03922603 #m^2/m^3 come back here
    SA_Vsmall = 16.99784605 #m^2/m^3
    m_COPVbig = 8.24 #kg, is an average
    m_COPVsmall = 5.965 #kg, is an average

    error = 1 #placeholder
    m_usedguess = 2.2 #kg
    dt = 0.05 #s
    depl_time = max(m_LOx/mdot_LOx, m_IPA/mdot_IPA) #s


    while error > 0.00001:
        time = 0 #s
        m_used = 0 #kg
        T_w = 300 #k
        mdots_LOx = []
        mdots_IPA = []
        mdots_total = []
        S = propssi("S", "T", T_b0, "P", P_b0, "nitrogen") #J/kgK, flow assumed to be isentropic so specific S stays constant
        rho_b0 = propssi("D", "P", P_b0, "T", T_b0, "nitrogen") #kg/m^3, flow from bottle assumed to be adiabatic
        rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen") #kg/m^3
        V_b = m_usedguess/(rho_b0-rho_end) #m^3, guess bottle volume based on m_used
        if V_b < 0.013:
            SA = SA_Vsmall*V_b
            m_COPV = m_COPVsmall
        elif V_b > 0.013:
            SA = SA_Vbig*V_b
            m_COPV = m_COPVbig
        C_p = 950 #J/kgK, is a ballpark. come back here
        #C_p = (liner mass/total mass)(C_p liner) + (wrap mass/total mass)(C_p wrap)
        k_w = 5*SA/(m_COPV*C_p) #s^-1, newton's law of heating/cooling. 5 W/m^2K convective heat transfer coefficient ballpark, https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html#gsc.tab=0, come back here
        rho_b1 = rho_b0-m_pres/V_b
        u_b0 = propssi("UMASS", "S", S, "D", rho_b0, "nitrogen")
        u_b1 = propssi("UMASS", "S", S, "D", rho_b1, "nitrogen")
        h_in = V_b*(rho_b0*u_b0-rho_b1*u_b1)/m_pres

        def fill_mass(V_u, h_in):
            m_0 = V_u * propssi("D", "P", 101325, "T", T_amb, "nitrogen")
            u_u0 = propssi("UMASS", "T", T_amb, "P", 101325, "air")
            def f(dm):
                m_ppres = m_0 + dm
                u_f = (m_0*u_u0 + dm*h_in)/m_ppres
                return propssi("P", "D", m_ppres/V_u, "UMASS", u_f, "nitrogen") - P_t
            dm_hi = 2*propssi("D", "P", P_t, "T", T_amb, "nitrogen")*V_u #just a safe upper bound, can mostly ignore
            return optimize.brentq(f, 0, dm_hi)

        # m_presLOx = fill_mass(V_uLOx, h_in) #no line volume COME BACK HERE
        # m_presIPA = fill_mass(V_uIPA, h_in)
        # m_pres = m_presLOx * 1.5 + m_presIPA #add a check on this? come back here, kinda sketch
        m_pres = fill_mass(V_utotal, h_in) #can't account for CF_pp of LOx with this

        rho_b = rho_b0-(m_used+m_pres)/V_b #need to recalculate these after m_pres calc
        rho_b1 = rho_b0-m_pres/V_b #these are for the start of burn btw
        P_b1 = propssi("P", "D", rho_b1, "S", S, "nitrogen")
        T_b1 = propssi("T", "P", P_b1, "D", rho_b1, "nitrogen")
        
        while (time <= depl_time):
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
            T_b = propssi("T", "S", S, "D", rho_b, "nitrogen")
            #k_GN2 = 25*SA/(m_b*propssi("CVMASS", "T", T_b, "P", P_b, "nitrogen"))
            #T_w += (k_w*(T_amb-T_w)+k_GN2*(T_w-T_b))*dt #this assumes that GN2 temperature isnt affected by conduction from the wall which is inaccurate but conservative
            time += dt
        error = abs(m_usedguess-m_used-m_pres)
        m_usedguess = m_used+m_pres

    if abs(P_b-P_tmin) <= 5*6894.757: #5 psi
        endreached = True
    else:
        endreached = False
    T_f = propssi("T", "S", S, "D", rho_b, "nitrogen")

    print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    print(f"V_b: {1000*V_b:.2f}L, T_w: {T_w:.0f}K, fill mass: {rho_b0*V_b:.3f} kg, pressurization mass: {m_pres:.3f} kg\nP_b0: {P_b0/6894.757:.0f} psi, P_b1: {P_b1/6894.757:.0f} psi, T_b1: {T_b1:.0f} K, T_f: {T_f:.0f}\nMin pressure reached: {endreached}")
    print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    print("---------------------------------------------------------------------------------------")
    return V_b, max(mdots_LOx), max(mdots_IPA), max(mdots_total), endreached, rho_b0*V_b, P_b1, T_b1, m_pres, T_f, T_w, mdots_total


def tryP_b0(P_b0):
    output = size_bottle(P_b0)
    x.append(P_b0/6894.757293) #psi
    y.append(1000*output[0]) #L
    mdots.append((output[1], output[2], output[3]))
    if output[4] == True:
        colors.append("blue")
    else:
        colors.append("red")
    fill_masses.append(output[5])
    press_masses.append(output[8])
    P_b1s.append(output[6]/6894.757293) #psi
    T_b1s.append(output[7])
    T_fs.append(output[9])
    if output[9] > 215.15:
        Tcolors.append("blue")
    else:
        Tcolors.append("red")
    all_mdots.append(output[11])

P_bstart = 3000*6894.757293 #Pa, full fill pressure
while P_bstart <= 6500*6894.757293: #6500 psi, sweep end
    tryP_b0(P_bstart)
    P_bstart += 20*6894.757293 #20 psi

points = sorted(zip(x, y, colors, mdots, fill_masses, P_b1s, T_b1s, press_masses, T_fs, Tcolors, all_mdots), key = lambda p: p[0])
sorted_x, sorted_y, sorted_colors, sorted_mdots, sorted_fill_masses, sorted_Pb1s, sorted_Tb1s, sorted_press_masses, sorted_Tfs, sorted_Tcolors, sorted_all_mdots = zip(*points)
ax[0].plot(sorted_x, sorted_y, color = "gray", linestyle = "-")
ax[1].plot(sorted_y, sorted_Tfs, color = "gray", linestyle = "-") #is this right

csvname = input("CSV name: ")
if not csvname.endswith(".csv"):
    csvname += ".csv"
with open(csvname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["P_b0_psi", "P_b1_psi", "V_b_L", "T_b1_K", "T_f_K", "fill_mass_kg", "press_mass_kg", "mdot_LOx_max_kg_s", "mdot_IPA_max_kg_s",
                     "mdot_total_max_kg_s", "all_mdots_kg_s", "min_pressure_reached", "T_f_ok"])
    for px, vy, col, (lox, ipa, tot), fm, P1, T1, pm, T_f, Tc, amd in zip(sorted_x, sorted_y, sorted_colors, sorted_mdots, sorted_fill_masses, sorted_Pb1s, sorted_Tb1s, sorted_press_masses, sorted_Tfs, sorted_Tcolors, sorted_all_mdots):
        writer.writerow([px, P1, vy, T1, T_f, fm, pm, lox, ipa, tot, amd, col == "blue", Tc == "blue"])

sc = ax[0].scatter(sorted_x, sorted_y, c = sorted_colors, zorder = 3)
sc2 = ax[1].scatter(sorted_y, sorted_Tfs, c= sorted_Tcolors, zorder = 3)
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

cursor2 = mplcursors.cursor(sc2, hover=mplcursors.HoverMode.Transient)
@cursor2.connect("add")
def on_add(sel):
    V = sorted_y[sel.index]
    TF = sorted_Tfs[sel.index]
    sel.annotation.set_text(f"{V:.2f} L tank\n{TF:.3f}K final temp")
    sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)

plt.show()
