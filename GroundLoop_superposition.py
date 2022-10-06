'''
Script for calculating superposition of thermal responses for the moving finite line source
model of Molina Giraldo in GLHE_groundwater.py

Many values are hard coded for this initial script.

'''
import sys

sys.path.insert(0, '../modules')

import datetime
import pandas
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import special
import math
from scipy import integrate


class Data:

    def __init__(self, gw, k, ps, cs, pw, cw, n, to, H):
        """Aquifer attributes are assigned or calculated from the input

        In addition to the attributes that are passed directly as parameters
        the following parameters are passed as additional arguments.

        These parameters are not bound to the class as they may vary.

        Args:
            gw: groundwater flow rate (specific discharge)
            k:  thermal conductivity
            ps: density of solids
            cs: specific heat capacity of solids
            pw: density of water
            cw: specific heat capacity of water
            n:  porosity
            H:  borehole depth
        """

        self.k = k
        self.gw = gw
        self.ps = ps
        self.cs = cs
        self.pw = pw
        self.cw = cw
        self.n = n
        self.pc = n * pw * cw + (1 - n) * ps * cs
        self.a = k / self.pc
        self.vt = (gw * pw * cw) / self.pc
        self.To = to
        self.H = H

class Configuration:
    """configuraiton of borefield and observations
    Attributes:
        nx_obs = number of observation location in x direction
        ny_obs = number of observations in y-direction
        B = borehole spacing [m]
        nx_b = number of boreholes in x-direction
        ny_b = number of boreholes in y-direction
        rb = radius of borehole [m]

    """
    def __init__(self, nx_obs, ny_obs, B, nx_b, ny_b, rb, x_obs, y_obs):
        self.nx_obs = nx_obs
        self.ny_obs = ny_obs
        self.nx_b = nx_b
        self.ny_b = ny_b
        self.B = B
        self.rb = rb
        self.x_obs = x_obs
        self.y_obs = y_obs


class Aquifer:
    """Parameters assigned as static characteristics of the aquifer and GLHE

    Attributes:
        k  =  thermal conductivity
        gw  =  groundwater flow rate (specific discharge)
        pc  =  effective heat capacity of medium
        a  =  thermal diffusivity
        vt  =  dimensionless velocity
        To  =  undisturbed temperature of aquifer

    """

    def __init__(self, gw, k, ps, cs, pw, cw, n, to):
        """Aquifer attributes are assigned or calculated from the input

        In addition to the attributes that are passed directly as parameters
        the following parameters are passed as additional arguments.

        These parameters are not bound to the class as they may vary.

        Args:
            ps: density of solids
            cs: specific heat capacity of solids
            pw: density of water
            cw: specific heat capacity of water
            n:  porosity
            H:  borehole depth
        """

        self.k = k
        self.gw = gw
        self.pc = n * pw * cw + (1 - n) * ps * cs
        self.a = k / self.pc
        self.vt = (gw * pw * cw) / self.pc
        self.To = to


class gwModels:
    """A collection of analytical and seminanalytical solutions

    All models use dimensionless parameters as defined in Molina-Giraldo (2011)

    Attributes:
        H  =  borehole depth
        r  =  radial distance away from borehole
        Rprime  =  dimensionless distance
        a  =  thermal diffusivity
        vt  =  dimensionless velocity
        pe  =  Peclet Number
        phi  =  azimuthal direction of location relative to groundwater flow

    Dimensionless parameters that include vertical location (z) and/or time (t)
    are passed to models.
    """

    def __init__(self, x, y, h, vt, a, k):
        """ Assignment and calculation of attributes based on iput

        NOTE that DB has been added as an attribute so it can be determined
        if integration on theta is necessary.  Result will be integrated on
        theta if abs(r - DB/2)< r

        Also NOTE that DB is input in units of inches
        """

        self.H = h
        self.r = math.sqrt(x * x + y * y)
        self.Rprime = math.sqrt(x * x + y * y) / h
        self.a = a
        self.k = k
        self.vt = vt
        self.pe = (vt * h) / a
        if (x >= 0) and (y > 0):
            self.phi = math.atan(x/y)
        elif (x >= 0) and (y < 0):
            self.phi = 1/2*math.pi - math.atan(y/x)
        elif x <= 0 and y < 0:
            self.phi = math.pi - math.atan(x/y)
        elif (x <= 0) and (y > 0):
            self.phi = 3/2*math.pi - math.atan(y/x)
        elif (x >= 0) and (y == 0):
            self.phi = 0
        elif (x <= 0) and (y == 0):
            self.phi = math.pi

        self.A = math.exp((self.pe / 2.0) * self.Rprime * math.cos(self.phi))

    '''
    Integrands
    '''

    def f_mfls(self, z, t):
        """"
        Integrand for moving finite line source model of Molina-Giraldo (2011)
        equation 17.
        """
        fo = self.a * t / (self.H * self.H)
        ulim = 1
        llim = 0
        integrand = []
        xvalues = np.linspace(llim, ulim, num=2000)
        for Zprime in xvalues:
            r = math.sqrt(self.Rprime ** 2 + ((z / self.H) - Zprime) ** 2)
            arg1 = (self.pe * r) / 2
            arg2 = (r - (self.pe * fo)) / (2 * math.sqrt(fo))
            arg3 = (r + (self.pe * fo)) / (2 * math.sqrt(fo))
            y = (1 / (4 * r)) * (math.exp(-arg1) * math.erfc(arg2) +
                                 math.exp(arg1) * math.erfc(arg3))
            integrand.append(y)
        return xvalues, integrand

    def f_mils(self, t):
        """
        Integrand for moving finite line source model of Molina-Giraldo (2011)
        equation 11 and Sutton et al (2003) equation 11.
        """
        llim = (self.r ** 2) / (4.0 * (self.a * t))
        ulim = 5  # can't us np.inf as upper limit
        integrand = []
        xvalues = np.logspace(llim, ulim, num=80)
        for eta in xvalues:
            f = (1 / eta) * math.exp(-eta - (self.pe ** 2 * self.Rprime ** 2) / (16 * eta))
            integrand.append(f)
        return xvalues, integrand

    '''
    Dimensionless temperature calculations
    '''

    def Tmfls(self, z, t):
        """
        Moving finite line source model transient solution from Molina-
        Giraldo (2011), eqn 16

        The horizontal location of interest (R') is set as class attribute,
        the vertical location (z) is provided as an argument.

        5.22.17: Problems arise with scipy.integrate.quad when the integrand
        is very narrow and very steep.  This occurs for small values of time
        and distance around the value z/H. To rectify this, an option is used
        to inform the algorithm to focus on this location that is idenfied with
        variables peak1 (range [0,1]) and peak2 (range [-1,0].
        """

        if self.pe < 600.:
            fo = self.a * t / (self.H * self.H)

            def y(zprime):
                r = np.sqrt(self.Rprime ** 2 + ((z / self.H) - zprime) ** 2)
                arg1 = (self.pe * r) / 2
                arg2 = (r - (self.pe * fo)) / (2 * math.sqrt(fo))
                arg3 = (r + (self.pe * fo)) / (2 * math.sqrt(fo))
                y = (1 / (4 * r)) * (math.exp(-arg1) * math.erfc(arg2) +
                                     math.exp(arg1) * math.erfc(arg3))
                return y

            if abs(z / self.H) < 1:
                peak1 = [z / self.H]
                peak2 = [-z / self.H]
                first = integrate.quad(lambda Zprime: y(Zprime), 0.0, 1.0, points=peak1)
                second = integrate.quad(lambda Zprime: y(Zprime), -1.0, 0.0, points=peak2)
            else:
                first = integrate.quad(lambda Zprime: y(Zprime), 0.0, 1.0)
                second = integrate.quad(lambda Zprime: y(Zprime), -1.0, 0.0)
            result = first[0] - second[0]

            return 1.0/(4.0*math.pi*self.k)*2.0 * self.A * result
        else:
            print ('warning: large Peclet number, not computing Tmfls ', self.pe)
            return

    def Tmflss(self, z):  # noqa
        """
        Moving finite line source model steady state solution from Molina-
        Giraldo (2011), eqn 18

        The horizontal location of interest (R') is set as class attribute,
        the vertical location (z) is provided as an argument.
        """
        if self.pe < 600.:
            # Fo = self.a*t/(self.H*self.H)
            ulim = 1
            llim = 1E-60
            ulim2 = 1E-60
            llim2 = -1

        def y(zprime):
            r = np.sqrt(self.Rprime ** 2 + ((z / self.H) - zprime) ** 2)
            arg = (-self.pe / 2) * r
            y = (1 / r) * math.exp(arg)
            return y

        first = integrate.quad(lambda Zprime: y(Zprime), llim, ulim)
        second = integrate.quad(lambda Zprime: y(Zprime), llim2, ulim2)
        result = first[0] - second[0]
        return 1.0/(4.0*math.pi*self.k)*self.A * result

    def Tmils(self, t):  # noqa
        """
        Moving infinite line source model transient solution from Molina-
        Giraldo eqn 11 and Sutton et al eqn 11.

        The limits of integration are from Sutton et al (2003) and differ from
        those given by Molina-Giraldo (2011).  It has been found that the
        two solutions are equal for most cases and when the Sutton equation
        is more computationally stable for very large values of time.
        """
        fo = (self.a * t) / (self.H * self.H)
        llim = (self.Rprime ** 2) / (4.0 * fo)
        ulim = np.inf

        def y(eta):
            yv = (1 / eta) * math.exp(-eta - (self.pe ** 2 * self.Rprime ** 2) / (16 * eta))
            return yv

        result = integrate.quad(lambda eta: y(eta), llim, ulim)
        return 1.0/(4.0*math.pi*self.k)*self.A * result[0]

    def Tmilss(self):  # noqa
        """
        Steady state moving infinite line source model from Molina-Giraldo
        eqn 12 and Sutton et al (2003) eqns 16 and 19
        """
        arg = (self.pe / 2) * self.Rprime
        besselk = special.kn(0, arg)
        return 1.0/(4.0*math.pi*self.k)*2.0 * self.A * besselk

    def Theis(self, t):  # noqa
        u = (self.r * self.r) / (4.0 * self.a * t)
        w = scipy.special.expn(1, u)
        print(u,w)
        return 1.0/(4.0*math.pi*self.k)*w



def log_interp(zz, xx, yy):
    logz = np.log10(zz)
    logx = np.log10(xx)
    logy = np.log10(yy)
    return np.power(10.0, np.interp(logz, logx, logy))


def glhe_groundwater_model(times, params, config, load):

    """
    Calculates thermal response at some specified location (usually borehole
    wall) to heatflow in to and out of ground loop. Temporal data can be at
    any resolution as long as units are correct (seconds, W/m, and Celcius).
    """

    x_values = []
    y_values = []
    Tmg_values = []
    num_g_calcs = 1

    # Establish aquifer parameters, dimensionless parameters,
    # and initialize models

    aquifer = Aquifer(params.gw, params.k, params.ps, params.cs,
                      params.pw, params.cw, params.n, params.To)

    x_bore = np.arange(0, config.nx_b*config.B, config.B).tolist()
    y_bore = np.arange(0, config.ny_b*config.B, config.B).tolist()

    bore_grid = [[row, col] for row in x_bore for col in y_bore]

    if (config.nx_obs > 1) or (config.ny_obs > 1):
        x_obs = np.arange(-config.B, config.nx_b*config.B, config.nx_b*config.B/config.nx_obs).tolist()
        y_obs = np.arange(-config.B, config.ny_b*config.B, config.ny_b*config.B/config.ny_obs).tolist()
        obs_grid = [[row, col] for row in x_obs for col in y_obs]
    else:
        obs_grid = [[config.x_obs, config.y_obs]]

    s = []

    for t in times:
        g = []
        X = []
        Y = []
        for x1, y1 in obs_grid:
            theta_loc = []
            x_row = []
            y_col = []
            for x2, y2 in bore_grid:
                x = x1 - x2
                y = y1 - y2
                if abs(x) < config.rb:
                    x = config.rb
                if abs(y) < config.rb:
                    y = config.rb
                #Initialize models class with new locaitons and set value of z to the midpoint
                glhe_gw = gwModels(x, y, params.H, aquifer.vt, aquifer.a, aquifer.k)
                z = glhe_gw.H/2

                # Compute simulated values of ground loop temperature for each location
                delT = glhe_gw.Tmfls(z, t)  # calls the desired function
                # Append delT values to theta array for each of the borehole locations at time t
                theta_loc.append(delT)

            X.append(x1)
            Y.append(y1)
            g.append(np.asarray(theta_loc).sum())   # sum the delT values over boreholes locations for time t and store in g

        #before moving the next time, store time series for each observation s
        s.append(g)
        

    return X, Y, s


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
    fig.colorbar(cmesh, ax=ax) #, location='top')
    plt.xlabel(r'$x [m]$')
    plt.ylabel(r'$y [m]$')
    plt.show()

    # fig = plt.gcf()
    # fig.set_size_inches(9, 5.5)
        # imagefile = '../../ModelOutput.png'
        # plt.savefig(imagefile, dpi=300)

if __name__ == "__main__":
    '''
    script used for testing and debugging of classes and functions
    '''
    # params = Data(gw=5.e-7, k=1.5, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100)
    #config = Configuration(nx_obs=24, ny_obs=24, B=3, nx_b=3, ny_b=3, rb=0.07)

    # Call funtion so that 'result' is what is 'returned'.  In this case, two arrays, one with times the other with drawdowns
    # X, Y, s = glhe_groundwater_model(times, params, config, load)

    # plot_heatmap(X, Y, s)

   



