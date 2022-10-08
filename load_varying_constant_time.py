'''
Script for calculating spatial distribution of thermal drawdown with multiple boreholes
'''
import sys

sys.path.insert(0, '../modules')

import matplotlib.pyplot as plt
import numpy as np
import math
import GroundLoop_superposition as groundloop

def plot_load_series(times, X, Y, s):

    x = X[0]
    y = Y[0]

    delT = []
    for i in range(len(s)):
        delT.append(s[i][0])
    z = np.asarray(delT)
    print(x, y, z)


    fig, ax = plt.subplots()

    ax.plot(loads, z)
    #fig.colorbar(cmesh, ax=ax)  # , location='top')
    plt.xlabel('load, 0 - 7')
    plt.ylabel('Thermal drawdown [$^\circ C$]')
    plt.show()

    # fig = plt.gcf()
    # fig.set_size_inches(9, 5.5)
    # imagefile = '../../ModelOutput.png'
    # plt.savefig(imagefile, dpi=300)


if __name__ == "__main__":
    '''
    script used for generating time series at a single location 
    '''
    times = [10**9] 

    loads = np.linspace(0, 7, num=8) 
    
    params = groundloop.Data(gw=5.e-17, k=1.5, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)
    config = groundloop.Configuration(nx_obs=1, ny_obs=1, B=3, nx_b=1, ny_b=2, rb=0.07, x_obs=-3, y_obs=0)
    # Call funtion so that 'result' is what is 'returned'.  In this case, two arrays, one with times the other with drawdowns
    X, Y, s = groundloop.glhe_groundwater_model(times, params, config, loads)


    plot_load_series(loads, X, Y, s)
