#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This function simulates a variable speed pumping system. First, you enter in the pump and system curve information. 
The function then calculates a curve relating rpm to flow and rpm to kW.
This is useful for estimating the flow of variable speed pumping systems
when no flow meter is available.

Example 1:

#Step 1: Create pump and system curves
PCurve_Flow = np.array([0, 164.14, 328.28, 492.61, 656.81, 761.88, 821.02, 952.51, 1039.33]) * 0.0864 #Convert MLD to L/s
PCurve_dkPa   = np.array([16.03, 15.73, 15.51, 13.87, 11.82, 10.21, 9.26, 6.70, 4.42]) * 9.81 #Convert m to kPa
PCurve_kW = np.array([48.83, 55.83, 71.34, 81.34, 88.57, 91.63, 90.22, 80.55, 66.88])
PCurve_eff = np.array([0, 0.4528, 0.6989, 0.8225, 0.8583, 0.8313, 0.8254, 0.7758, 0.6722])
Sys_curve = np.full((9), dP_meas) #A flat system of 9 kPa 

#Step 2.:Define a reasonable RPM operating range. Note that if there is no intersection with the VFD speed, pump curve, or
#system curve, the function below will not work.
RPM_Max = 670
RPM_Min = 515
rpm_range = np.arange(RPM_Min, RPM_Max, 1)

#Step 3. Generate Curves that relate rpm to flow and rpm to kW
Curve_rpm_to_flow, Curve_rpm_to_kW = pumpingSim(PCurve_Flow, PCurve_dkPa, PCurve_kW, Sys_curve, rpm_range)
Curve_flow_to_rpm = Flow_to_RPM(PCurve_MLD, PCurve_dkPa, Sys_curve, rpm_range)

@author: Nicholas Zibin
@date: Created on Fri Sep 11 15:41:46 2020
@license: MIT
"""

# Program to interpolate pump curve and system curve 
from __future__ import division 
import numpy as np

# The inperolation function is copied from:
#https://coderedirect.com/questions/501614/find-the-intersection-of-two-curves-given-by-x-y-data-with-high-precision-in
def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc

def pumpingSim(PCurve_Flow, PCurve_dkPa, PCurve_kW, Sys_curve, RPM_range):
    '''
    This function simulates a variable speed pumping system with a static lift 
    head curve. First, you enter in the pump and system curve information. 
    The function then calculates a curve relating rpm to flow and rpm to kW.
    This is useful for estimating the flow of variable speed pumping systems
    when no flow meter is available.
        
    Args:
        PCurve_Flow -> Flow data points from pump curve
        PCurve_dkPa -> Diff. Pressure data points from pump curve  
        PCurve_kW   -> Power data points from pump curve  
        Sys_curve   -> Data points from system curve  
        RPM_range   -> RPM at which you want to calculate kW
        
    Returns:
        Polynomial curve of rpm and flow data (Curve_rpm_to_flow)
        Polynomial curve of rpm and kW data (Curve_rpm_to_kW)
        
    '''
    Max_RPM = np.max(RPM_range)
    vfd_range = RPM_range/Max_RPM 
    
    #Calculate flow as a function of pump RPM
    flow_arr = np.array([])
    for vfd_spd in vfd_range:
        flow, yc = interpolated_intercept(PCurve_Flow * vfd_spd, PCurve_dkPa * vfd_spd**2, Sys_curve)
        flow_arr = np.append(flow_arr,flow)

    Curve_rpm_to_flow = np.poly1d(np.polyfit(vfd_range * Max_RPM, flow_arr, 2))

    #Calculate power as a function of pump RPM
    kW_arr  = np.array([])
    for vfd_spd in vfd_range:        
        Curve_Flow_to_kW = np.poly1d(np.polyfit(PCurve_Flow * vfd_spd, PCurve_kW * vfd_spd**3, 3)) #Establish new pump curve at the VFD speed using affinity laws
        flow = Curve_rpm_to_flow(vfd_spd * Max_RPM)
        kW = Curve_Flow_to_kW(flow)
        kW_arr  = np.append(kW_arr, kW)
        
    Curve_rpm_to_kW = np.poly1d(np.polyfit(vfd_range * Max_RPM, kW_arr, 3))
    
    return Curve_rpm_to_flow, Curve_rpm_to_kW


def Flow_to_RPM(PCurve_Flow, PCurve_dkPa, Sys_curve, RPM_range):
    '''This function works similar to pumpingSim but instead creates a function
        of flow to rpm.
    '''
    
    Max_RPM = np.max(RPM_range)
    vfd_range = RPM_range/Max_RPM 
    
    #Calculate flow as a function of pump RPM
    flow_arr = np.array([])
    for vfd_spd in vfd_range:
        flow, yc = interpolated_intercept(PCurve_Flow * vfd_spd, PCurve_dkPa * vfd_spd**2, Sys_curve)
        flow_arr = np.append(flow_arr, flow)

    Curve_flow_to_rpm = np.poly1d(np.polyfit(flow_arr, RPM_range, 2))
    
    return Curve_flow_to_rpm
