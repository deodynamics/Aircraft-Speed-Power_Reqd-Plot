
###############################################################################
#
# File name: graph_speed_vs_thp_reqd.py
# Version: 1.0
# Date: 02 - 11 - 2023
# Coded by: Deomary Angelo Franco
#
# Description:
# This script is used to calculate the Thrust Horsepower required for an
# aircraft to produce enough thrust to overcome drag (and maintain equilibrium)
# at a given flight speed.
# 
# Sea level values are first calculated and then V vs THP pairs for different
# altitudes are then derived from the sea leve values, which is calculated
# from the density changes (through the introduction of density ratio).
#
# CLI structure:
# "python graph_speed_vs_thp_reqd.py W S C_L_max C_D_o k alt_min alt_max
# alt_interval save_df_as_csv save_graph_as_png"
#
# Example:
# python graph_speed_vs_thp_reqd.py 1423.15 117 1.52 0.027 0.049 0 25000
# 5000 True True
#
# Variable definitions:
#
# W --> Gross weight of aircraft in lbs
# S --> Wing Planform area in sq.ft
# C_L_max --> Maximum Lift Coefficient of the wing
# C_D_o --> Drag coefficient at zero lift of aircraft
# k --> the factor in C_D_i = k * C_L**2. k is calculated from (pi*e*A)**-1
# alt_min --> minimum altitude in ft. 0 corresponds to Sea Level
# alt_max --> maximum altitude in ft.
# alt_interval --> interval for all altitudes of interest. 
#   Note: make sure that only 7 altitude values will be formed from these due
#   to code limitations.basically limited by available colors. HEHEHE!
# save_df_as_csv (optional) --> if you want to save the calculated table as csv
# save_graph_as_png (optional) --> if you want to save plot as png
#
###############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys

# numpy, pandas, and matplotlib must be installed if this script is to be used.
# sys is already a built-in and no need for pip install.

# Some constants used for calculation. These are sea level values:
rho_o = 0.002377 #density at sea level [slug/ft^3]
a = -0.003566 #lapse rate [R/ft]
T_o = 519 #temperature at sea level [Rankine]

# Some coversion values.
v_fps_to_mph = 15 / 22 # ft/s --> mph
v_mph_to_fps = 22 / 15 # mph --> ft/s
power_ftlb_to_hp = 1 / 375 # ft-lb/s --> hp


def sigma_as_f_of_h(h):
    # Calculate sigma density ratio as function of altitude h.
    return (1 + a*h/T_o)**4.26

def col_for_array(altitudes):
    # Ready the 2-level column for the power reqd table.

    # To map different altitudes to V and THP_reqd. 
    # If altitude is 0, it is mapped to sea level and other more parameters
    # since these are considered for the initial calculations
    col_dict = {}
    for h in altitudes:
        if h == 0.:
            col_dict['Sea Level'] = [
                'V_o', 'C_L', 'C_D', 'T_reqd', 'THP_reqd'
            ]
        else:
            col_dict[f'{str(int(h))} ft'] = ['V', 'THP_reqd']
    
    # A list of tuples e.g. (10_000 ft, V) to be used to cosntruct the multi-
    # level column for the dataframe.
    col_tuple = [
        (alt, param) for alt, params in col_dict.items() for param in params
    ]

    # Constructs the multi-index column.
    col = pd.MultiIndex.from_tuples(col_tuple, names=['altitude', 'param'])
    return col

def main_table(W, S, C_L_max, C_D_o, k, altitudes):
    # Creates the dataframe of data containing some aerodynamic properties.
    # But important data shows here speed corresponding to the thrust
    # horsepower reqd.

    # If there are missing values given in the commandline, it will raise a
    # a value error.
    if any(val is None for val in (W, S, C_L_max, C_D_o, k, altitudes)):
        raise ValueError
    
    else:
        # Calculate the stalling speed at sea level.
        V_s = np.sqrt((2*W)/(rho_o*S*C_L_max)) * v_fps_to_mph

        # The array containing speeds of 50 up to 250 mph. 21 is chosen as
        # the number of elements so that the interval will be 10 mph.
        # The stalling speed will be appended at the end and then the array is 
        # sorted.
        V_o__c = np.sort(
            np.append(
                np.linspace(50, 250, 21),
                V_s
            )
        )
        # Array calculating the lift coefficient from the different speeds.
        C_L__c = 2 * W / (rho_o * S * (v_mph_to_fps * V_o__c ) ** 2)
        # Drag coefficient calculated from the drag polar.
        C_D__c = C_D_o + k * C_L__c ** 2
        # Thrust required calculated from the ratio based from the fact that
        # equilibrium is assumed and maintained.
        T_reqd__c = C_D__c * W / C_L__c
        # Thrust Horsepower (or simply power) required to produce enough thrust
        # required to overcome the drag.
        THP_reqd__c = T_reqd__c * V_o__c * power_ftlb_to_hp

        # Create an empty array where speed and thp pair will be stored here.
        V_and_THP_array = np.array([])
        # Loop through the different values of altitude.
        for h in altitudes:
            # If altitude is zero, it will be renamed as sea level and reshaped
            # to become a column vector using np.newaxis.
            if h == 0:
                V_and_THP_array = THP_reqd__c[:, np.newaxis]

            # For all other altitudes, just rename them to a more readable
            # form. Solve for the corresponding V and THP based from the
            # decrease in density from the altitude increase. They will then 
            # be stacked to create a 2-dimensional ndarray and later
            # transformed to a dataframe.
            else:
                # Density ratio
                sigma = sigma_as_f_of_h(h)
                # Solve for speed corresponding to the same CL but different
                # altitude (hence the sigma or density ratio).
                V__c = V_o__c / np.sqrt(sigma)
                # Solve for the THP corresponding the the same CL but different
                # altitude (hence the sigma or density ratio).
                THP_reqd_h__c = THP_reqd__c / np.sqrt(sigma)

                V_and_THP_array = np.hstack(
                    (
                        V_and_THP_array,
                        V__c[:, np.newaxis],
                        THP_reqd_h__c[:, np.newaxis]
                    )
                )
        # Stack horizontally the V and THP pairs to the previously calculated 
        # sea level parameter columns.
        power_reqd_array = np.hstack(
            (
                V_o__c[:, np.newaxis],
                C_L__c[:, np.newaxis],
                C_D__c[:, np.newaxis],
                T_reqd__c[:, np.newaxis],
                V_and_THP_array
            )
        )

        # Transforms the array into a dataframe where the multi-level columns
        #  are created from our custom function.
        power_reqd_df = pd.DataFrame(
            power_reqd_array,
            columns=col_for_array(altitudes)
        )

        return power_reqd_df

def plot_V_vs_THP_reqd(W, S, C_L_max, C_D_o, k, altitudes, save_df_as_csv=True,
                       save_graph_as_png=True):
    # Plots the V and THP pair for all altitudes based from the created
    # dataframe.

    # Calculates the stall speed one more time.
    V_s = np.sqrt((2*W)/(rho_o*S*C_L_max)) * v_fps_to_mph

    # Calls the previous function to create the dataframe.
    power_reqd_df = main_table(W, S, C_L_max, C_D_o, k, altitudes)
    
    # A list containining the colors that will be randomly selected for the V
    # and THP pair for each altitude. But due to only having 7 available
    # colors, this limits the number of altitudes we can plot.
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    # Shuffles the colors.
    np.random.shuffle(colors)
    # The size of our figure.
    plt.figure(figsize=(12, 12))

    # Gets the maximum value of THP and rounds it to the nearest hundreds. This
    # becomes the ymax for our figure which is shown below.
    thp_reqd_max = (
        power_reqd_df
        .loc[:,power_reqd_df.columns.get_level_values(1) == 'THP_reqd']
        .max(axis=None)
    )
    ymax = int(np.ceil(thp_reqd_max/100) * 100.)

    # Gets the maximum value of speed and rounds it to the nearest hundreds.
    # This becomes the xmax for our figure which is shown below.
    v_max = (
        power_reqd_df
        .loc[:,power_reqd_df.columns.get_level_values(1) == 'V']
        .max(axis=None)
    )
    xmax = int(np.ceil(v_max/100) * 100.)

    # Explicitly call the maximum for the axis, depending on the values of our 
    # speed and THP pair.
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)

    # Since the stall speed was just added to the array and that array was
    # sorted, we need to get its index. The index is used to get the
    # corresponding THP. This will be used to plot a dot to show
    # the stalling values for V and THP.
    stall_speed_index = (
        power_reqd_df
        .loc[power_reqd_df[('Sea Level', 'V_o')] == V_s]
        .index
    )

    # To get the THP for stalling speed. The stall_speed_index was used here.
    THP_reqd_for_stall = (
        power_reqd_df
        .loc[stall_speed_index, ('Sea Level', 'THP_reqd')]
        .to_numpy()[0]
    )

    # Loop through the different altitudes then plot them.
    for h in altitudes:

        # This will be used as parameters for the plots.
        # Gets a color for the colors list.
        color = np.random.choice(colors)

        # Then once a color is selected, it will be removed from the list
        # so it will not be used anymore.
        colors.remove(color)

        # Format for the V and THP curve.
        line_format = f'{color}-'

        # Format for the point corresponding the stall point.
        point_format = f'{color}o'


        if h == 0:
            # Values used for the plot. This is defined since it is different
            # from other altitudes where their label is just the number itself
            # and uses V instead of V_o (which is used here).
            alt_label = 'Sea Level'
            x = power_reqd_df[(alt_label, 'V_o')]

            # The dot or point for the stall point.
            plt.plot(V_s, THP_reqd_for_stall, point_format)

            # Dashdotted vertical and horizontal lines connecting the stall
            # point to the x and y axis.
            plt.vlines(
                x=V_s,
                ymin=0,
                ymax=THP_reqd_for_stall,
                colors=color,
                linestyles='dashdot',
            )
            plt.hlines(
                y=THP_reqd_for_stall,
                xmin=0,
                xmax=V_s,
                colors=color,
                linestyles='dashdot',
            )

        else:
            # The labels are called out here.
            alt_label = f'{str(int(h))} ft'
            x = power_reqd_df[(alt_label, 'V')]

        # sSince both will be using  this as the y value hence coded outside
        # the if-block.
        y = power_reqd_df[(alt_label, 'THP_reqd')]
        # This code plots the curve for all values of altitude.
        plt.plot(x, y, line_format, label=alt_label)

    # Additional stuff to plot.
    plt.legend(loc=2)
    plt.xlabel('Flight Speed, V(mph)')
    plt.ylabel('Thrust Horsepower Required, $THP_{{reqd}}$ (hp)')

    # Label is broken down in this way to avoid being coded too long.
    stall_label = 'Stall point: '
    v_label = f'\n$V_s$ = {V_s:.2f} mph '
    thp_label = f'$THP_{{reqd}}$ = {THP_reqd_for_stall:.2f} hp'
    plt.text(
        x=V_s + 10, 
        y=THP_reqd_for_stall, 
        s=stall_label + v_label + thp_label,
        va='top',
        ha='left'
    )
    
    plt.title(
        'Flight Speed vs Thrust Horsepower Required at different altitudes'
    )


    # 'save_df_as_csv' and 'save_graph_as_pdf' are optional and defaulted as
    # True, hence if not explicitly added to the arguments for the command
    # line, they will be saved as csv and png.
    # If False are used for the last two arguments, none will be saved.
    # Saved to the same directory as the python file
    if save_df_as_csv:
        power_reqd_df.to_csv('Speed vs THP data.csv')

    if save_graph_as_png:
        plt.savefig('Speed vs THP image.png', format='png', dpi='figure')
    
    # To show the plot. Command line is unusable unless the plot is closed
    # manually.
    plt.show()


def main():

    # When running python graph_speed_vs_the_reqd.py <arguments>
    # 10 arguments (only 8 are required since the last 2 are optional) are 
    # needed since they are parameters for needed for the calculations.
    W = float(sys.argv[1])
    S = float(sys.argv[2])
    C_L_max = float(sys.argv[3])
    C_D_o = float(sys.argv[4])
    k = float(sys.argv[5])

    # For the altitude, make sure that only 7 altitude values are produced
    # since the usable color strings are only 7.
    alt_min = float(sys.argv[6])
    alt_max = float(sys.argv[7])
    alt_interval = float(sys.argv[8])
    save_df_as_csv = bool(sys.argv[9])
    save_graph_as_png = bool(sys.argv[10])
    
    # A ndarray based from the provided arguments.
    altitudes = np.arange(alt_min, alt_max + alt_interval, alt_interval)

    plot_V_vs_THP_reqd(W, S, C_L_max, C_D_o, k, altitudes, save_df_as_csv,
                       save_graph_as_png)

if __name__ == '__main__':
    main()
