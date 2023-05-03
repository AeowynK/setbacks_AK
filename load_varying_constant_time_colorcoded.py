'''
Script for calculating spatial distribution of thermal drawdown with multiple boreholes
'''
import sys

sys.path.insert(0, '../modules')

import matplotlib.pyplot as plt
import numpy as np
import math
import GroundLoop_superposition as groundloop

def plot_load_series(loads, delT):

    fig, ax = plt.subplots()

    ax.plot(loads, delT)
    #fig.colorbar(cmesh, ax=ax)  # , location='top')
    plt.xlabel('load, 0 - 7')
    plt.ylabel('Thermal drawdown [$^\circ C$]')
    plt.show()



def plot_loadcondmap(x, y, delT):

    # generate 2 2d grids for the load & k (thermal conductivity) bounds
    loads, ks = np.meshgrid(np.linspace(0, 10, 25), np.linspace(1.5, 5, 25))

    # Define the minimum and maximum temperature values for the plot
    temp_min = 0
    temp_max = 2

    # Define the fixed color range for the plot
    cmap = plt.cm.RdYlGn_r
    norm = plt.Normalize(temp_min, temp_max)
    cmap.set_over('lightgrey')


##    # impact thresholds
##    thresholds = [0, 2, 4]
##
##    # define the colormap object
##    cmap = plt.cm.get_cmap('BuPu', len(thresholds)-1)


    # create the pcolormesh plot
    fig, ax = plt.subplots()
    # define the temp levels 
    levels = np.linspace(1, 2, 3)
    c = ax.pcolormesh(x, y, delT, cmap=cmap, norm=norm)
    # vmin=min(thresholds), vmax=max(thresholds))
    # colors from: https://matplotlib.org/stable/tutorials/colors/colormaps.html

    # an alternative method for coloring in plot, using discrete colors for intervals
##    levels = [0, 1, 2, 3, 4]
##    origin = 'lower'
##    c = ax.contourf(x, y, delT, levels,
##                       colors=('g', 'y', 'y', 'r', 'r'),
##                       origin=origin,
##                       extend='both',
##                       alpha=0.5)

    l = ax.contour(x, y, delT, colors = 'Black', linewidths=1, levels = levels)

    fmt = {}
    strs = ['negligible impact', 'possible impact', 'request technical review']
    for label, s in zip(l.levels, strs):
        fmt[label] = s

    ax.clabel(l, l.levels, inline=True, fmt=fmt, fontsize = 10) #levels = levels)
    ax.set_title('9 boreholes: Temp. change in °C at \n12m distance over 15 years & \ninter-borehole spacing of 6m', fontsize = 17)

    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])

    #create the colorbar and label it 
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label("temperature change,°C", fontsize = 15)

    # Create axis labels
    plt.ylabel("net annual ground load, W/m", fontsize = 15)
    plt.xlabel("thermal conductivity, W/(m⋅K)", fontsize = 15)

    ax.annotate('neighboring system experiences \n'
                'excessive thermal drawdown \n'
                'sustained operation unlikely',
                xy=(2, 8), xytext=(2, 8))
    #show plot
    plt.show()



if __name__ == "__main__":
    '''
    script used for generating time series at a single location 
    '''
    loads = np.linspace(0, 10, 25)
    ks = np.linspace(1.5, 6, 25)
    times = [4.734*(10**8)]
    # 15 years in seconds 
    delT = []
    

    i = 0
    j = 0 
    for load in loads:
        j = 0 
        delT_inner = []
        
        for k in ks:
            
            params = groundloop.Data(gw=5.e-17, k=k, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)
            # If nx_obs or ny_obs are not equal to 1, then grid is constructed.  Otherwise, a single location
            config = groundloop.Configuration(nx_obs=1, ny_obs=1, B=6, nx_b=3, ny_b=3, rb=0.07, x_obs=24, y_obs=3)
            # borehole brid starts at the origin and proceeds into first quadrant
            # Call funtion so that 'result' is what is 'returned'.  In this case, two arrays, one with times the other with drawdowns
            X, Y, s = groundloop.glhe_groundwater_model(times, params, config, load)
            
            delT_inner.append(s[0][0])

            i = i + 1

        delT.append(delT_inner)
        j = j + 1
        
    #plot_load_series(loads, delT)
    plot_loadcondmap(ks, loads, delT)
