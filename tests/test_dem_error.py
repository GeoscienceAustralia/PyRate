import numpy as np
from scipy.interpolate import griddata
import math

# read base.par list
ifg_dates = []
baseline_list = '/g/data/dg9/INSAR_ANALYSIS/MEXICO/S1/UNIT-TEST-DATA/CROP-A/lists-and-config-files/from-geotiff/' \
                'pyrate_baseline_files.list'

with open(baseline_list, 'r') as f:
    for line in f.readlines():
        basepar = line.split('/')[8]
        ifg_dates.append(basepar[:17])

print(ifg_dates)

# interpolation location (azimuth and range coordinate of crop A pixel (50,29)
az0 = 2636.9912
rg0 = 103.4759

# round azimuth and range coordinates to closest step (500 for az, 200 for rg)
azstep = 500
rgstep = 200
az1 = azstep * math.floor(az0/azstep)
rg1 = rgstep * math.floor(rg0/rgstep)
az2 = azstep * math.ceil(az0/azstep)
rg2 = rgstep * math.ceil(rg0/rgstep)

# four coordinates for bi-linear interpolation
teststr1 = str(az1).rjust(6, ' ') + str(rg1).rjust(7, ' ')
teststr2 = str(az1).rjust(6, ' ') + str(rg2).rjust(7, ' ')
teststr3 = str(az2).rjust(6, ' ') + str(rg1).rjust(7, ' ')
teststr4 = str(az2).rjust(6, ' ') + str(rg2).rjust(7, ' ')

# loop through all corresponding bperp.par files in base.par list
int_path = '/g/data/dg9/INSAR_ANALYSIS/MEXICO/S1/GAMMA/T005A/INT/'
bperp_int = np.empty(shape=(len(ifg_dates)))
i = 0
for ifg in ifg_dates:
    # read Mexico city bperp file
    bperpfile = int_path + ifg + '/' + ifg + '_VV_8rlks_bperp.par'
    with open(bperpfile, 'r') as f:
        for line in f.readlines():
            if teststr1 in line:
                bperp1 = line.split()[7]
            if teststr2 in line:
                bperp2 = line.split()[7]
            if teststr3 in line:
                bperp3 = line.split()[7]
            if teststr4 in line:
                bperp4 = line.split()[7]

    # setup numpy array for bi-linear interpolation
    n = np.array([(az1, rg1, bperp1),
                  (az1, rg2, bperp2),
                  (az2, rg1, bperp3),
                  (az2, rg2, bperp4)])
    # interpolate using scipy function "griddata"
    bperp_int[i] = griddata(n[:,0:2], n[:,2], [(az0, rg0)], method='linear')
    i += 1
print(bperp_int)


# the calculation of Bperp values in PyRate is done on-the-fly in module dem_error.py
# extract az, rg and corresponding bperp from PyRate for pixel (50,29), e.g. using the following lines of code:
 #  print(az_parts[50, 29])
 #  print(rg_parts[50, 29])
 #  print(bperp[:, 50, 29])
 #-> resulting in the following output:
#2636.9912
#103.47585
#[  33.48592183    3.44669685  -75.37369399  -26.88597679  -33.25298942
# -108.84360354    3.74075472   -3.19700977  -14.58390611   10.06920291
#  -51.12649599   -5.74544068  -17.36872483  -30.4772929     7.12691256
#  -37.68943916  -73.14248882  -11.45674522  -24.64851804   12.69928323
#  -32.16248418  -20.86746046   61.514626     48.30928659  -13.17640207
#   24.28126177  -36.84111057  -20.5870326    77.8291117    -8.66115426]
pyrate_bperp = np.array([33.48592183, 3.44669685, -75.37369399, -26.88597679, -33.25298942, -108.84360354, 3.74075472, \
                         -3.19700977, -14.58390611, 10.06920291, -51.12649599, -5.74544068, -17.36872483, -30.4772929, \
                         7.12691256, -37.68943916, -73.14248882, -11.45674522, -24.64851804, 12.69928323, -32.16248418,\
                         -20.86746046, 61.514626, 48.30928659, -13.17640207, 24.28126177, -36.84111057, -20.5870326, \
                         77.8291117, -8.66115426])
print(pyrate_bperp)

# compare the GAMMA and PyRate Bperp estimates
diff = bperp_int - pyrate_bperp
print(1000*np.mean(diff)) # mean difference in mm
print(1000*np.max(diff)) # max difference in mm