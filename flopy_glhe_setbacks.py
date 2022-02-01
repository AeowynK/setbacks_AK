# MODFLOW 2005:
# use flopy to simulate thermal drawdown in geothermal borehole heat exchanger(s)


import numpy as np
import flopy
import flopy.utils as fpu
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf


# assigning model info:

Lx = 20.0
Ly = 20.0
# size of model, in meters
ztop = 0.0
zbot = -50.0
# thickness, b
nlay = 1
nrow = 100
ncol = 100
# number of rows and columns
delr = Lx / ncol
delc = Ly / nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
hk = 2.0
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

# define stress periods
nper = 1
# just one stress period

perlen = [182 * 86400]
# period length, 182 days in seconds
nstp = [50]
# number of steps
steady = [False]
# steady state condition

# creating the static flopy objects (not time-dependent):
modelname = "heat_flow"
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

pcg = flopy.modflow.ModflowPcg(mf)

# create the well package object, of type: flopy.modflow.ModflowWel()
# Remember to use zero-based layer, row, column indices!

pumping_rate = -500.0
# thermal 'discharge' of 10 W/m multiplied by thickness, b
wel_sp1 = [[0, nrow / 2 - 1, ncol / 2 - 1, pumping_rate], [0, nrow / 4, ncol / 4, pumping_rate]]
stress_period_data = {0: wel_sp1}
wel = flopy.modflow.ModflowWel(mf, stress_period_data=stress_period_data)

# create the output control package, of type: flopy.modflow.ModflowOc()
stress_period_data = {}
for kper in range(nper):
    for kstp in range(nstp[kper]):
        stress_period_data[(kper, kstp)] = [
            "save head",
            "save drawdown",
            "save budget",
            "print head",
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
idx = (0, int(nrow / 2) - 1, int(ncol / 2) + 0)
ts = headobj.get_ts(idx)
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
ttl = "Temperature ({0},{1},{2})".format(idx[0] + 1, idx[1] + 1, idx[2] + 1)
ax.set_title(ttl)
ax.set_xlabel("time [s]")
ax.set_ylabel("Temp Change[C]")
ax.plot(ts[:, 0], ts[:, 1], "bo-")
plt.show()

# post-processing results

# extracting the heads
hds = bf.HeadFile(f"{modelname}.hds")
times = hds.get_times()
head = hds.get_data(totim=times[-1])

# contouring the heads
extent = (delr / 2.0, Lx - delr / 2.0, Ly - delc / 2.0, delc / 2.0)
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
ax.contour(head[0, :, :], levels=np.arange(-1, 0, .1), extent=extent)

# extracting the cell-by-cell flows
#cbb = bf.CellBudgetFile(f"{modelname}.cbc")
kstpkper_list = cbb.get_kstpkper()
#frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[-1])[0]
#fff = cbb.get_data(text="FLOW FRONT FACE", totim=times[-1])[0]
#qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
#    (frf, fff, None), mf, head
#)

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


