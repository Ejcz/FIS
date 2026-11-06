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

def getfloat(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Pls try again: ")

def readCSV(filename, row, colname):
    with open(filename, newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index(colname)
    return rows[row][col]

def size():
    OD = getfloat("Line OD (in): ")
    thickness = getfloat("Thickness (in): ")
    ID = OD - 2*thickness
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


filedir = input("File directory: ")
if not filedir.endswith("/"):
    filedir += "/"
filename = filedir + input("Filename: ")
#filename = "C:/Users/Chris/Documents/Github/FIS/COPV/" + input("Filename: ")
if not filename.endswith(".csv"):
    filename += ".csv"

while True:
    row = input("Row #: ")
    while row.isdigit() == False:
        row = input("Row #: ")
    row = int(row)-1
    mdots = ast.literal_eval(readCSV(filename, row, "all_mdots_kg_s"))
    T_b1 = float(readCSV(filename, row, "T_b1_K"))
    P_b0 = float(readCSV(filename, row, "P_b0_psi"))*6894.757
    P_b1 = float(readCSV(filename, row, "P_b1_psi"))*6894.757
    V_b = float(readCSV(filename, row, "V_b_L"))/1000
    S_b0 = propssi("S", "P", P_b1, "T", T_b1, "nitrogen")
    m_b = float(readCSV(filename, row, "fill_mass_kg")) - float(readCSV(filename, row, "press_mass_kg"))
    
    rhocheck = propssi("D", "P", P_b1, "T", T_b1, "nitrogen")
    if abs((rhocheck*V_b-m_b)/(rhocheck*V_b)) < 0.03:
        print("rho ok")
    else:
        print("rho not ok")

    area, P_y, P_b = size()
    upmachs = []
    downmachs = []
    
    for mdot in mdots:
        #upstream
        rho_b = m_b/V_b
        P_bot = propssi("P", "S", S_b0, "D", rho_b, "nitrogen")
        h_0 = propssi("HMASS", "S", S_b0, "D", rho_b, "nitrogen") #flow across reg is nonisentropic but is isenthalpic
        upmachs.append(calcMach(P_bot, h_0, S_b0, mdot, area))

        #downstream
        S_down = propssi("S", "HMASS", h_0, "P", P_set, "nitrogen") #no fanno flow btw. come back here
        downmachs.append(calcMach(P_set, h_0, S_down, mdot, area))
        #rho_d = propssi("D", "P", P_set, "T", , "nitrogen")
        
        m_b -= mdot*dt
    print(f"Upstream: {max(upmachs):.3f} mach #, {P_y/P_b0:.2f} SF_y, {P_b/P_b0:.2f} SF_b\nDownstream: {max(downmachs):.3f} mach #, {P_y/P_f:.2f} SF_y, {P_b/P_f:.2f} SF_b")

