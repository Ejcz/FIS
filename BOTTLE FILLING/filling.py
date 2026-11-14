from CoolProp.CoolProp import PropsSI as propssi
from scipy import integrate, optimize, constants
import numpy as np
import matplotlib.pyplot as plt
import mplcursors
import csv
import os

CF = 2
P_t = 3548704 #Pa, postreg P. inconsistent? come back here
P_tmin = 6894757.293 #Pa, min P into reg
T_amb = 300 #K adjust: ambient temperature
ISA_B1 = 2.17 # m^2 internal surface area of linde #1. come back here
m_b1 = 2.5985663 #kg full mass linde #1
m_B1 = 18 # kg metal mass linde #2
v_b1 = 0.017 # m^3
ISA_B2 = 1.655*((2.29/2)**2) # m^2 of linde #2. come back here
m_b2 = 15.6260454 #kg gas mass linde #2
m_B2 = 60 #kg metal mass linde #2. come back here
v_b2 = 0.043 # m^3

R_gn2 = 296.8390789 # kg/JK
C_316 = 477 # J/kgK, both lindes made of 316 SS
convect_WtoA = 5 #W/m^2K wall to ambient
convect_GtoW = 150 #W/m^2K gn2 to wall
dt = 0.05 #s


filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "300K_COPVs_transient.csv") #adjust: filename for ambient T
def readCSV(filename, colname):
    data = []
    with open(filename, newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index(colname)
    rows = rows[1:]
    for row in rows:
        data.append(row[col])
    return data

partnums = readCSV(filename, "Part_num")
v_COPVs = readCSV(filename, "V_b_L")
m_COPVs = readCSV(filename, "weight_kg")
P_fills = readCSV(filename, "P_b0_psi")
fill_masses = readCSV(filename, "fill_mass_kg")







