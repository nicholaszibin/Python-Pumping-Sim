# Python-Pumping-Sim

This function simulates a variable speed pumping system to calculate flow and power as a function of rpm This is useful for estimating the flow of variable speed pumping systems when no flow meter is available. This is particularly useful when the system curve is mainly static-head dominated rather than friction dominated. The function can also estimate rpm as a function of flow.

The required inputs to the function are the pump curve, system curve, and pump RPM range to be evaluated. The function outputs are polynomical curves relating:
- RPM to Flow
- RPM to Power
- Flow to RPM

See sample-pumpingSim.py for an example of using this.
