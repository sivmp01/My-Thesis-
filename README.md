Building Demand-Side Management & Thermal Storage Optimization

This repository contains MATLAB code and data developed for a master’s thesis on Demand-Side Management (DSM), Thermal Storage, and Reserve Market Participation (mFRR) using a multi-storey building model.

The project focuses on:

1. Reducing electricity cost

2. Peak shaving

3. Smart heating schedules

4. Reserve bidding (Case 3)

5. Risk-aware decisions using T-CVaR

**Repository Structure**

_Main MATLAB Codes_

1. Case_1.m — DSM without storage

2. Case_2.m — DSM + thermal storage

3. Case_3_Normal.m — DSM + storage + reserve (risk-neutral)

4. EachHour_Case3_TCVaR_Reserve.m — DSM + storage + reserve + T-CVaR

5. input_file.m - input data

6. Tin_Day.m - the thermal model code

7. compute_TMF.m — Computes thermal mass factor function

8. Building_RC_no_heat.m — RC model without heating function code 

9. Building_RC_Thermostat.m — RC model with HVAC thermostat

10. LoadWeather.m — Loads weather data

11. Price_data_filteration.m - Loads price 


_Plot Scripts_

These files contain codes that executes the graphs of the results

1. PlotCASE1.m — Plots Case 1

2. PlotCASE2.m — Plots Case 2

3. PlotCASE3TCVaR.m — Plots Case 3 T-CVaR

_Data Files (.mat)_

These are the files of my execution

1. JanFeb_Prices.mat — Electricity prices

2. weather_hourly.mat — Outdoor temperature

3. RCinput.mat — Thermal parameters

4. thermal.mat — Building thermal properties

6. Baseline_out.mat, Case1_out.mat, Case2_out.mat, Case3_TCVaR.mat — Stored results

**Code Exection Structure**
1. Execute input_file.mat
2. Execute Tin_Day.m
3. Then the Cases and the plot files 


**_If you would like to read my thesis report kindly search the title in divaportal.se - or my name in the same portal_**
