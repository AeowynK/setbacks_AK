# MODFLOW 2005:
# use flopy to simulate thermal drawdown in geothermal borehole heat exchanger(s)


import numpy as np
import flopy
import flopy.utils as fpu
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import pandas as pd

# assigning model info: [metric]

Lx = 50.
Ly = 50.
# size of model, in meters
ztop = 0.0
zbot = -1.0
# thickness, b, meters
nlay = 1
nrow = 200
ncol = 200
# number of rows and columns
delr = Lx / ncol
delc = Ly / nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
hk = 2.5
# thermal conductivity, units of W/mK
# sy = 0.1 # specific yield, don't need bc confined
ss = 2.0e6
# specific storage in J/m^3*K
laytyp = 0
# zero for confined

# creating the BAS package variables
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
# starting head set at zero. head represents thermal reservoir in this model
strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)

'''
get data
+'''

data = pd.read_csv('./data/varying_load_periodic.csv')
Q = np.asarray(data['load [W/m]'].fillna(float(0)))

'''
set up stress period times
'''
# define stress periods
nper = len(data)
nstp = np.ones(nper, dtype=int)*5
perlen = np.ones(nper)*3600

# multiplier on ts
tsmult = 1.2

# steady state condition - False for transient
steady = np.zeros(nper, dtype=bool)


# creating the static flopy objects (not time-dependent):
modelname = "./superposition/varying_load_test"
mf = flopy.modflow.Modflow(modelname, exe_name="./mflow/mf2005")
# executable called
dis = flopy.modflow.ModflowDis(
    mf,
    nlay,
    nrow,
    ncol,
    delr=delr,
    delc=delc,
    tsmult=1.2, # time step multiplier 
    top=ztop,
    botm=botm[1:],
    nper=nper,
    perlen=perlen,
    nstp=nstp,
    steady=steady,
)
# discretization of parameters, tsmult = time step multiplier
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
lpf = flopy.modflow.ModflowLpf(
    mf, hk=hk, ss=ss, laytyp=laytyp, ipakcb=53)

ghb_t = 0.0
ghb_boundary = []
for il in range(nlay):
    ghb_c = hk * (ghb_t - zbot) * delc
    for ir in range(nrow):
        #ghb_boundary.append([il, ir, 0, ghb_t, ghb_c])
        ghb_boundary.append([il, ir, ncol - 1, ghb_t, ghb_c])
    for ic in range(ncol):
        #ghb_boundary.append([il, 0, ic, ghb_t, ghb_c])
        ghb_boundary.append([il, nrow-1, ic, ghb_t, ghb_c])

print("Adding ", len(ghb_boundary), "GHBs for stress period 1.")
stress_period_data = {0: ghb_boundary, 1: ghb_boundary}

ghb = flopy.modflow.ModflowGhb(mf, stress_period_data=stress_period_data)

pcg = flopy.modflow.ModflowPcg(mf)

# create the well package object, of type: flopy.modflow.ModflowWel()
# Remember to use zero-based layer, row, column indices!

#well info
x_wellcells = [50]
y_wellcells = [50]
well_layers = [0]

#well package for w/m energy input/extraction

def BH_StressPeriod_Dict(Well_layers, Well_rows, Well_columns, Q):
    sp_list = []
    for count, Q_val in enumerate(Q, 0):
        sp_list.append(count)
        sp_temp = []
        for layer in Well_layers:
            for row in Well_rows:
                for column in Well_columns:
                    sp_temp.append((layer, row, column, Q_val))
        sp_list.append(sp_temp)

    it = iter(sp_list)
    sp_dict = dict(zip(it, it))
    return sp_dict

sp_dict = BH_StressPeriod_Dict(well_layers, x_wellcells, y_wellcells, Q)


wel = flopy.modflow.ModflowWel(mf, stress_period_data=sp_dict)

# create the output control package, of type: flopy.modflow.ModflowOc()
stress_period_data = {}
for kper in range(nper):
    for kstp in range(nstp[kper]):
        stress_period_data[(kper, kstp)] = [
            "save head",
            "save drawdown",
            "save budget",
            "print budget",
        ]
oc = flopy.modflow.ModflowOc(
    mf, stress_period_data=stress_period_data, compact=True
)

# run the model using the run_model method

# Write the model input files
mf.write_input()

# Run the model
success, mfoutput = mf.run_model(silent=True, pause=False)
if not success:
    raise Exception("MODFLOW did not terminate normally.")

# post-processing results
# read temps from the MODFLOW binary output file, using flopy.utils.binaryfile()

# create temp plot from model outputs

# Create the headfile and budget file objects
headobj = bf.HeadFile(modelname + ".hds")
times = headobj.get_times()
cbb = bf.CellBudgetFile(modelname + ".cbc")

# Plot the temperature change versus time
idx = (0, int(nrow / 4) - 1, int(ncol / 4) - 1)
ts = headobj.get_ts(idx)
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1)
ttl = "Temperature change in °C at ({0},{1},{2}) for a thermal conductitivity of \n2.5 [W/m*K] and varying load [W/m]".format(idx[0] + 1, idx[1] + 1, idx[2] + 1)
ax.set_title(ttl)
ax.set_xlabel("time [s]")
ax.set_ylabel("Temperature Change[°C]")
ax.plot(ts[:, 0], ts[:, 1], "bo-")
plt.show()

# post-processing results

# extracting the heads
hds = bf.HeadFile(f"{modelname}.hds")
times = hds.get_times()
head = hds.get_data(totim=times[-1])

# contouring the heads
extent = (delr / 2.0, Lx - delr / 2.0, Ly - delc / 2.0, delc / 2.0)
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
CS = ax.contour(head[0, :, :], levels=np.arange(0, 1, .01), extent=extent)
ax.clabel(CS, CS.levels, inline=True, fontsize=8)
# extracting the cell-by-cell flows
#cbb = bf.CellBudgetFile(f"{modelname}.cbc")
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[159])[0]
fff = cbb.get_data(text="FLOW FRONT FACE", totim=times[159])[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
    (frf, fff, None), mf, head)
#plt.plot(qx[0][0])
#plt.ylim(-1, 1)
#plt.show()

# creating the plot/figure
# fig = plt.figure(figsize=(6, 6))
# ax = fig.add_subplot(1, 1, 1, aspect="equal")
# modelmap = flopy.plot.PlotMapView(model=mf, layer=0, ax=ax)
# qm = modelmap.plot_ibound()
# lc = modelmap.plot_grid()
# cs = modelmap.contour_array(head, levels=np.linspace(0, 10, 11))
# quiver = modelmap.plot_vector(qx, qy)

# remember to add flopy contour map of temperature (color contour or line contour)
plt.show()


