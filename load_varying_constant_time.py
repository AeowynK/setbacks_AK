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

    # fig = plt.gcf()
    # fig.set_size_inches(9, 5.5)
    # imagefile = '../../ModelOutput.png'
    # plt.savefig(imagefile, dpi=300)


def plot_loadcondmap(X, Y, s)

    # generate 2 2d grids for the load & k (thermal conductivity) bounds
    loads, ks = np.meshgrid(np.linspace(0, 10, 25), np.linspace(1.5, 5, 25))

    z = (1 - x / 2. + x ** 5 + y ** 3) * np.exp(-x ** 2 - y ** 2)
    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    z = z[:-1, :-1]
    z_min, z_max = -np.abs(z).max(), np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
    ax.set_title('load and thermal conductivity')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    plt.xlabel("load factor, unitless")
    plt.ylabel("thermal conductivity, W/m")

    plt.show()



if __name__ == "__main__":
    '''
    script used for generating time series at a single location 
    '''
    times = [10**9] 
    delT = []
    # loads = np.linspace(0, 7, num=8)
    loads, k = np.meshgrid(np.linspace(0, 10, 25), np.linspace(1.5, 5, 25))

    i = 0
    j = 0 
    for load in loads:
        j = 0 
        i = i + 1 
        for k in ks:
            j = j + 1
            params = groundloop.Data(gw=5.e-17, k=k, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)
        # If nx_obs or ny_obs are not equal to 1, then grid is constructed.  Otherwise, a single location
            config = groundloop.Configuration(nx_obs=1, ny_obs=1, B=3, nx_b=1, ny_b=2, rb=0.07, x_obs=-3, y_obs=0)
        # Call funtion so that 'result' is what is 'returned'.  In this case, two arrays, one with times the other with drawdowns
            X, Y, s = groundloop.glhe_groundwater_model(times, params, config, load)

            delT[i][j].append(s[0][0])

    
    plot_load_series(loads, delT)
