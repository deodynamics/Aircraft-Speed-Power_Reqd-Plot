# Aircraft-Speed-Power_Reqd-Plot
A short script that calculates aicraft's thrust horsepower required for a given aicraft speed.

Description:
This script is used to calculate the Thrust Horsepower required for an aircraft to produce enough thrust to overcome drag (and maintain equilibrium) at a given flight speed.
 
Sea level values are first calculated and then V vs THP pairs for different altitudes are then derived from the sea leve values, which is calculated from the density changes (through the introduction of density ratio).

CLI structure:
"python graph_speed_vs_thp_reqd.py W S C_L_max C_D_o k alt_min alt_max alt_interval save_df_as_csv save_graph_as_png"

Example:
python graph_speed_vs_thp_reqd.py 1423.15 117 1.52 0.027 0.049 0 25000 5000 True True

Variable definitions:
    W --> Gross weight of aircraft in lbs
    S --> Wing Planform area in sq.ft
    C_L_max --> Maximum Lift Coefficient of the wing
    C_D_o --> Drag coefficient at zero lift of aircraft
    k --> the factor in C_D_i = k * C_L**2. k is calculated from (pi*e*A)**-1
    alt_min --> minimum altitude in ft. 0 corresponds to Sea Level
    alt_max --> maximum altitude in ft.
    alt_interval --> interval for all altitudes of interest. 
       Note: make sure that only 7 altitude values will be formed from these due to code limitations.
    save_df_as_csv (optional) --> if you want to save the calculated table as csv
    save_graph_as_png (optional) --> if you want to save plot as png
