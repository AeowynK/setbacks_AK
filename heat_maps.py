'''
Script for calculating spatial distribution of thermal drawdown with multiple boreholes
'''
import sys

sys.path.insert(0, '../modules')

import matplotlib.pyplot as plt
import numpy as np
import math
import GroundLoop_superposition as groundloop

def plot_heatmap(X, Y, s):
    x = np.asarray(X)
    y = np.asarray(Y)
    z = np.asarray(s[-1])

    nx = int(math.sqrt(len(X)))
    ny = int(math.sqrt(len(Y)))

    d = []
    e = []
    f = []
    count = 0
    for i in range(nx):
        a = []
        b = []
        c = []
        for j in range(ny):
            a.append(x[count])
            b.append(y[count])
            c.append(z[count])
            count = count + 1
        d.append(a)
        e.append(b)
        f.append(c)
    zmin, zmax = z.min(), z.max()
    levels = np.linspace(zmin - 0.1, zmax + 0.1, num=25)
    fig, ax = plt.subplots()

    contours = ax.contour(d, e, f, levels=np.asarray([1.0]))
    cmesh = ax.pcolormesh(d, e, f, cmap='jet', shading='auto')
    cbar = fig.colorbar(cmesh, ax=ax)  # , location='top')
    cbar.set_label("temperature change,°C", fontsize = 15)
    plt.title('Heat map of subsurface temperature change [°C] over \n50 years in a 3 x 3 borehole field with groundwater flow', fontsize = 15)
    plt.xlabel(r'$x [m]$')
    plt.ylabel(r'$y [m]$')
    plt.show()

    # fig = plt.gcf()
    # fig.set_size_inches(9, 5.5)
    # imagefile = '../../ModelOutput.png'
    # plt.savefig(imagefile, dpi=300)


if __name__ == "__main__":
    '''
    script used for generating color contour plots 
    '''
    #times = np.linspace(1, 10**8, num=5)
    times = [1.578*(10**9)]
    # 50 years in seconds 
    load = 7.0
    params = groundloop.Data(gw=5.e-7, k=1.5, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)
    #if nx_obs or ny_obs are greater than 1, then x_obs and y_obs are not used but need to be set.
    config = groundloop.Configuration(nx_obs=24, ny_obs=24, B=3, nx_b=3, ny_b=3, rb=0.07, x_obs=0, y_obs=0)
    # Call funtion so that 'result' is what is 'returned'.  In this case, two arrays, one with times the other with drawdowns
    X, Y, s = groundloop.glhe_groundwater_model(times, params, config, load)

    plot_heatmap(X, Y, s)
