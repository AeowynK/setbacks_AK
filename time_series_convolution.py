import datetime
import math

import matplotlib.pyplot as plt
import numpy as np

from GroundLoop_superposition import gwModels, Data


def plot_time_series(times, X, Y, s):

    x = X[0]
    y = Y[0]

    delT = []
    for i in range(len(s)):
        delT.append(s[i][0])
    z = np.asarray(delT)
    print(x, y, z)

    times_days = times/(86400*365)

    fig, ax = plt.subplots()

    ax.plot(times_days, z)
    #fig.colorbar(cmesh, ax=ax)  # , location='top')
    plt.xlabel('time [years]')
    plt.ylabel('Thermal drawdown [$^\circ C$]')
    plt.show()

    # fig = plt.gcf()
    # fig.set_size_inches(9, 5.5)
    # imagefile = '../../ModelOutput.png'
    # plt.savefig(imagefile, dpi=300)

def get_convolution(forcing, response, duration):
    W = np.asarray(response)
    # Block Response is first difference of well function
    Br = np.diff(W, n=1)
    # Insert first element of W into Block Response
    Br = np.insert(Br, 0, W[0])
    # Perform convolution
    j = np.convolve(forcing, Br)
    return j[0:duration:1]

def calculate_time_series(times, q, glhe_gw):

    '''
    Compute simulated values of ground loop temperature
    '''
    theta = []
    # quad_error=[]
    # integral=[]

    for t in times:
        delT = glhe_gw.Tmfls(glhe_gw.H/2, t)
        # integral.append(delT[0])
        # quad_error.append(delT[1])
        theta.append(delT)
    W = np.asarray(theta)
    Br = np.diff(W, n=1)  # Block Response is first difference of well function
    Br = np.insert(Br, 0, W[0])  # insert first element of W into Block Response
    j = np.convolve(q, Br)  # Perform convolution
    result = j[0:len(q):1]
    return result

if __name__ == "__main__":
    '''
    Establish aquifer parameters, dimensionless parameters, and initialize models
    '''
    times = np.arange(1, 15*365*86400, 30*86400)
    # input years desired here as multiplier
    q = 5 + 10*np.sin(np.pi*times/(180*86400))
    # specifies unbalanced base load and range

    params = Data(gw=5.e-17, k=2.5, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)

    for r in [1.5, 6, 9]:

        glhe_gw = gwModels(x=r, y=0, h=params.H, vt=params.vt, a=params.a, k=params.k)
        s = calculate_time_series(times,q, glhe_gw)
        years = times/(365*86400)
        plt.plot(years, s, label= 'r='+str(r)+'m')

    plt.legend()
    plt.title('Thermal drawdown at different distances for \n5 W/m unbalanced load with +/- 10 W/m operational range', fontsize = 14)
    plt.ylabel("temperature change,°C", fontsize = 14)
    plt.xlabel("time [years]", fontsize = 14)
    plt.show()
