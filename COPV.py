from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import matplotlib.pyplot as plt
import mplcursors
import csv

plt.close()
x = []
y = []
fig, ax = plt.subplots()
ax.set_xlabel("P_b0 (psi)")
ax.set_ylabel("V_b (L)")
ax.grid(True)
colors = []
mdots = []

def size_bottle(P_b0):
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
    #R = scipy.constants.
    T_amb = 316 #K
    T_b = 316 #K
    #P_b0 = 21766190.54 #Pa
    m_b0 = 2.8113311163 #kg
    m_used = 0 #kg
    m_pres = 0.2641458904 #kg
    V_uLOx = 0.00340 #m^3
    V_uIPA = 0.00307 #m^3

    error = 1 #placeholder
    m_usedguess = 1.7 #kg
    dt = 0.05 #s
    depl_time = max(m_LOx/mdot_LOx, m_IPA/mdot_IPA) #s


    while error > 0.0001:
        time = 0 #s
        m_used = 0 #kg
        mdots_LOx = []
        mdots_IPA = []
        mdots_total = []
        S = propssi("S", "T", T_b, "P", P_b0, "nitrogen") #J/kgK
        rho_b0 = propssi("D", "P", P_b0, "T", T_b, "nitrogen") #flow from bottle assumed to be adiabatic
        rho_end = propssi("D", "P", P_tmin, "S", S, "nitrogen")
        V_b = m_usedguess/(rho_b0-rho_end) #m^3, guess bottle volume based on m_used
        m_b = V_b * rho_b0
        
        while (time <= depl_time):
            rho_b = (m_b-m_used)/V_b #kg/m^3

            h_b = propssi("HMASS", "S", S, "D", rho_b, "nitrogen") #isentropic flow, specific H bottle = specific H tanks
            T_in = propssi("T", "P", P_t, "Hmass", h_b, "nitrogen") #T going into prop tanks
            rho_in = propssi("D", "P", P_t, "T", T_in, "nitrogen") #rho GN2 going into tanks
            
            mdot_LOx = Vdot_LOx * CF * rho_in
            mdot_IPA = Vdot_IPA * rho_in
            mdot_total = mdot_LOx + mdot_IPA
            mdots_LOx.append(mdot_LOx)
            mdots_IPA.append(mdot_IPA)
            mdots_total.append(mdot_total)

            m_used += mdot_total * dt
            rho_b = (m_b-m_used)/V_b #kg/m^3
            P_b = propssi("P", "S", S, "D", rho_b, "nitrogen")
            time += dt
        error = abs(m_usedguess-m_used)
        m_usedguess = m_used

    if abs(P_b-P_tmin) <= 5*6894.757: #5 psi
        endreached = True
    else:
        endreached = False
    print(f"m_used: {m_used:.3f} kg in {time:.3f}s")
    print(f"V_b: {1000*V_b:.3f}L, m_b: {m_b:.3f}kg, P_b0: {P_b0:.3f} Pa, Min pressure reached: {endreached}")
    print(f"LOx max m dot: {max(mdots_LOx):.3f}kg/s @ {(mdots_LOx.index(max(mdots_LOx))+1)*dt:.3f}s")
    print(f"IPA max m dot: {max(mdots_IPA):.3f}kg/s @ {(mdots_IPA.index(max(mdots_IPA))+1)*dt:.3f}s")
    print(f"Total max m dot: {max(mdots_total):.3f}kg/s @ {(mdots_total.index(max(mdots_total))+1)*dt:.3f}s")
    print("---------------------------------------------------------------------------------------")
    return V_b, max(mdots_LOx), max(mdots_IPA), max(mdots_total), endreached


def tryP_b0(P_b0):
    output = size_bottle(P_b0)
    x.append(P_b0/6894.757293) #psi
    y.append(1000*output[0]) #L
    mdots.append((output[1], output[2], output[3]))
    if output[4] == True:
        colors.append("blue")
    else:
        colors.append("red")


#execute
P_bstart = 1.379e+7 #Pa
while P_bstart < 3.4474e+7: #5000 psi
    tryP_b0(P_bstart)
    P_bstart += 20*6894.757293 #20 psi

points = sorted(zip(x, y, colors, mdots))
sorted_x, sorted_y, sorted_colors, sorted_mdots = zip(*points)
ax.plot(sorted_x, sorted_y, color = "gray", linestyle = "-")

csvname = input("CSV name: ")
if not csvname.endswith(".csv"):
    csvname += ".csv"
with open(csvname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["P_b0_psi", "V_b_L", "mdot_LOx_max_kg_s", "mdot_IPA_max_kg_s",
                     "mdot_total_max_kg_s", "min_pressure_reached"])
    for px, vy, col, (lox, ipa, tot) in zip(sorted_x, sorted_y, sorted_colors, sorted_mdots):
        writer.writerow([px, vy, lox, ipa, tot, col == "blue"])

sc = ax.scatter(sorted_x, sorted_y, c = sorted_colors, zorder = 3,)
cursor = mplcursors.cursor(sc, hover=True)
@cursor.connect("add")
def on_add(sel):
    LOx, IPA, total = sorted_mdots[sel.index]
    P = sorted_x[sel.index]
    V = sorted_y[sel.index]
    sel.annotation.set_text(f"{P:.3f} psi\n{V:.3f} L\n{LOx:.3f} kg/s LOx\n{IPA:.3f} kg/s IPA\n{total:.3f} kg/s total")
    sel.annotation.get_bbox_patch().set(fc="lightgreen", alpha=1)
plt.show()
