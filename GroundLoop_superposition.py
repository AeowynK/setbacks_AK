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

    def __init__(self, gw, k, ps, cs, pw, cw, n, to, H, phi):
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
        self.phi = phi
        self.H = H

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

    def __init__(self, x, y, h, vt, a, k, phi):
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
        self.phi = phi
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


def glhe_groundwater_model(params, x_locs, y_locs, times):

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

    x_row = []
    y_row = []
    delT_loc = []
    times = range[0, 1.578*(10**9)] #50 years span in seconds 

    #Initialize models class and set value of z to the midpoint
    
    glhe_gw = gwModels(x, y, params.H, aquifer.vt, aquifer.a, aquifer.k, params.phi)
    z = glhe_gw.H/2
            
    for t in times:

        delT = []  # create an empty array to store the delT values 
        theta = []
        Tmg_row = []
        
    #g = np.asarray(theta)

        for x, y in zip(x_locs, y_locs):


            # Compute simulated values of ground loop temperature
            

            #tsparse = np.logspace(math.log10(times[0]), math.log10(times[-1]), num=num_g_calcs)

            #delT = glhe_gw.Tmilss(z)

            delT_loc.append(theta) 
            x_row.append(x)
            y_row.append(y)
            #g.append(detT)            

            load = 7.0 # degree to which load is unbalanced
    
            B = 6.0   # in meters; CSA standard is 3m to property lines
    
            x_locs = [-3*B, -2*B, -1*B, -3*B, -2*B, -1*B, -3*B, -2*B, -1*B]
            y_locs = [-B, -B, -B, 0, 0, 0, B, B, B]
    
            s = np.asarray(result[2]).sum()    # sum the delT values
            temp = s*load    # multiply the sum by the load to get the change in temp
            

        delT = glhe_gw.Tmfls(z, t)  # calls the desired function 
        theta.append(delT)    # appends delT values to theta array   
        delT.append(delT)   # create array of summed delT values from each time
        
    return (x_row, y_row, delT_loc)


if __name__ == "__main__":
    
        params = Data(gw=5.e-16, k=1.5, ps=2650, cs=880, pw=1016, cw=3850, n=.1, to=0, H=100, phi=0.)
        
        # Call funtion so that 'result' is what is 'returned'
        result = glhe_groundwater_model(params, x_locs, y_locs, delT)
        
        #print(result)
    
   

##  want to get Q/k


