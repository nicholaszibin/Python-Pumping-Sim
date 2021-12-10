import numpy as np  
import matplotlib.pyplot as plt
from pumpingSim import pumpingSim
from pumpingSim import Flow_to_RPM

#Step 1: Create pump and system curves
PCurve_Flow = np.array([0, 164, 328, 492, 656, 761, 822, 952, 1039]) * 0.0864 #Convert MLD to L/s
PCurve_dkPa = np.array([16, 15.7, 15.5, 13.9, 11.8, 10.2, 9.3, 6.7, 4.4]) * 9.81 #Convert m to kPa
PCurve_kW   = np.array([49, 56, 71, 81, 89, 92, 90, 81, 67])
Sys_curve   = np.full((9), 9*9.81) #A flat system curve of 9 kPa 

#Step 2: Define a reasonable RPM operating range. Note that if there is no intersection with the VFD speed, pump curve, or
#system curve, the function below will not work.
RPM_Max = 650
RPM_Min = 550
rpm_range = np.arange(RPM_Min, RPM_Max, 1)

#Step 3: Generate Curves that relate rpm to flow and rpm to kW
Curve_rpm_to_flow, Curve_rpm_to_kW = pumpingSim(PCurve_Flow, PCurve_dkPa, PCurve_kW, Sys_curve, rpm_range)
Curve_flow_to_rpm = Flow_to_RPM(PCurve_Flow, PCurve_dkPa, Sys_curve, rpm_range)

#Sample Plot:
plt.plot(Curve_rpm_to_flow(rpm_range))
plt.grid(), plt.ylabel('Flow (L/s)'), plt.xlabel('Pump RPM');
