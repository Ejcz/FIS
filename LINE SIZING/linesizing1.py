from CoolProp.CoolProp import PropsSI as propssi
from scipy import constants
from scipy.optimize import brentq
import csv
import ast

def readCSV(filename, row, colname):
    with open(filename, newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index(colname)
    return rows[row][col]

row = 28
filename = "C:/Users/Chris/Documents/GitHub/FIS/COPV/316K_ambient.csv"
mdots = ast.literal_eval(readCSV(filename, row, "all_mdots_kg_s")) #COME BACK HERE, ACTUALLY ADD TO COPV
m_b = float(readCSV(filename, row, "fill_mass_kg")) - float(readCSV(filename, row, "press_mass_kg"))
T_b0 = float(readCSV(filename, row, "T_b1_K"))
P_b0 = float(readCSV(filename, row, "P_b1_psi"))*6894.757
V_b = float(readCSV(filename, row, "V_b_L"))/1000
S_b0 = propssi("S", "P", P_b0, "T", T_b0, "nitrogen") #post prepressurization, from COPV.py
P_set = 550*6894.757 #reg setpoint
P_f = 1.1*P_set
dt = 0.05 #s

def getfloat(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Pls try again: ")

def size():
    OD = getfloat("Line OD (in): ")
    thickness = getfloat("Thickness (in): ")
    ID = OD - 2*thickness
    OD *= 0.0254
    ID *= 0.0254
    area = constants.pi*((ID/2)**2)
    return area

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


rhocheck = propssi("D", "P", P_b0, "T", T_b0, "nitrogen")
if abs((rhocheck*V_b-m_b)/(rhocheck*V_b)) < 0.03:
    print("ok")
else:
    print("bad")

area = size()
machs = []
for mdot in mdots:
    #upstream
    rho_b = m_b/V_b
    P_b = propssi("P", "S", S_b0, "D", rho_b, "nitrogen")
    h_0 = propssi("HMASS", "S", S_b0, "D", rho_b, "nitrogen") #only applies upstream
    machs.append(calcMach(P_b, h_0, S_b0, mdot, area))

    #downstream
    #rho_d = propssi("D", "P", P_set, "T", , "nitrogen")
    

    m_b -= mdot*dt
print(f"{max(machs):.3f}")

