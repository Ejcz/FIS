from CoolProp.CoolProp import PropsSI as propssi
from scipy import constants
from scipy.optimize import brentq
import csv
import ast


P_set = 550*6894.757 #reg setpoint
P_f = 1.1*P_set #overpressure, less conservative than P_set tho
dt = 0.05 #s

S_y = 25000*6894.757 #yield strength
UTS = 70000*6894.757 #UTS

#adjust: dimensions in inches
OD_up = 0.5
t_up = 0.12
ID_up = OD_up-2*t_up
OD_down = 0.5
t_down = 0.028
ID_down = OD_down-2*t_down

rho_oks = []
COPVs = []
machs = []

def readCSV(filename, colname):
    data = []
    with open(filename, newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index(colname)
    rows = rows[1:]
    for row in rows:
        data.append(row[col])
    return data

def size(OD, ID):
    OD *= 0.0254
    ID *= 0.0254
    area = constants.pi*((ID/2)**2)

    P_y = S_y*((OD**2)-(ID**2))/((OD**2)+(ID**2)) #yield pressure
    P_b = UTS*((OD**2)-(ID**2))/((OD**2)+(ID**2)) #burst pressure

    return area, P_y, P_b

def calcMach(P0, h0, S0, mdot, area):
    def residual(P):
        h = propssi("HMASS", "S", S0, "P", P, "nitrogen")
        rho = propssi("D", "S", S0, "P", P, "nitrogen")
        return rho*max(2*(h0-h), 0)**(1/2)*area-mdot
    P = brentq(residual, 0.5*P0, P0*(0.999999)) #upper bound of sweep is this to avoid negative sqrt for velocity during rounding
    h = propssi("HMASS", "S", S0, "P", P, "nitrogen")
    A = propssi("A", "P", P, "S", S0, "nitrogen")
    v = (2*(h0-h))**(1/2)
    M = v/A
    return M


filedir = "C:/Users/Chris/Documents/Github/FIS/COPV/" #adjust: CSV file directory
filename = input("Filename: ")
if not filename.endswith(".csv"):
    filename += ".csv"
filepath = filedir + filename

mdots_raw = readCSV(filepath, "all_mdots_kg_s")
T_b1s = readCSV(filepath, "T_b1_K") #come back here
P_b0s = readCSV(filepath, "P_b0_psi")
P_b1s = readCSV(filepath, "P_b1_psi")
V_bs = readCSV(filepath, "V_b_L")
m_fills = readCSV(filepath, "fill_mass_kg")
m_presses = readCSV(filepath, "press_mass_kg")
COPVs = readCSV(filepath, "Part_num")

up_area, up_P_y, up_P_b = size(OD_up, ID_up)
down_area, down_P_y, down_P_b = size(OD_down, ID_down)

for i in range(len(COPVs)):
    up_machs = []
    down_machs = []

    COPV = COPVs[i]
    mdots = ast.literal_eval(mdots_raw[i])
    T_b1 = float(T_b1s[i])
    P_b0 = float(P_b0s[i])*6894.757
    P_b1 = float(P_b1s[i])*6894.757
    V_b = float(V_bs[i])/1000
    S_b0 = propssi("S", "P", P_b1, "T", T_b1, "nitrogen") #come back here
    m_b = float(m_fills[i]) - float(m_presses[i])

    rhocheck = propssi("D", "P", P_b1, "T", T_b1, "nitrogen")
    if abs((rhocheck*V_b-m_b)/(rhocheck*V_b)) < 0.03:
        print("rho ok")
    else:
        print("rho not ok")

    for mdot in mdots:
        #upstream
        rho_b = m_b/V_b
        P_bot = propssi("P", "S", S_b0, "D", rho_b, "nitrogen")
        h_0 = propssi("HMASS", "S", S_b0, "D", rho_b, "nitrogen") #flow across reg is nonisentropic but is isenthalpic
        up_machs.append(calcMach(P_bot, h_0, S_b0, mdot, up_area))

        #downstream
        S_down = propssi("S", "HMASS", h_0, "P", P_set, "nitrogen") #no fanno flow btw. come back here
        down_machs.append(calcMach(P_set, h_0, S_down, mdot, down_area))
        #rho_d = propssi("D", "P", P_set, "T", , "nitrogen")
        
        m_b -= mdot*dt
    machs.append([max(up_machs),max(down_machs)])

if len(filename) == 21:
    csvname = f"{filename[0:3]}K_{filename[-16:-11]}s_{OD_up}upOD_{t_up}upt_{OD_down}downOD_{t_down}downt.csv"
elif len(filename) == 20:
    csvname = f"{filename[0:3]}K_{filename[-15:-11]}s_{OD_up}upOD_{t_up}upt_{OD_down}downOD_{t_down}downt.csv"
with open(csvname, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Part_num", "up_mach", "down_mach"])
    for i in range(len(COPVs)): #come back here
        writer.writerow([COPVs[i], machs[i][0], machs[i][1]])


