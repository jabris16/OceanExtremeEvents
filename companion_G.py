'''

Companion file.

Set parameters, load eddyTracking and heatwave data.

Jamie Atkins

'''

#****************************************************************

# LIBRARIES

import numpy as np
from netCDF4 import Dataset
import heatwave_functions
import eddySupp_functions

#****************************************************************

# PARAMETERS (adjust as necessary)

pctile = 90 # percentile threshold of choice (90% in line with Hobday et al., 2016)
num_yearsCLIM = 5 # time length [years] of the climatological mean/threshold period
resolution = 0.25 # model spatial resolution
num_years = 23 # number of years of analysis
days_in_year = 365 # number of days in a year

run = 'run' # name of run, should match that from eddyTracking
nc_filename = 'dataset.nc' # filename of temperature NetCDF file

kelvin_convert = 273.15 # Kelvin to Celcius conversion factor; switch off if temperature dataset is in degrees C already

#****************************************************************

# CONSTANTS

# eddyTracking output (from Eric Oliver code)
path = './'
data = np.load(path + 'eddy_track_' + run + '.npz') # 'encoding' may be required if switching between Python 2. and Python 3.
eddies_tracked = data['eddies']

# temperature data
z = path + nc_filename
fileobj = Dataset(z)
temperature = fileobj.variables['temperature'][0:num_years * days_in_year,0,:,:]  - kelvin_convert
lon = fileobj.variables['longitude'][:]
lat = fileobj.variables['latitude'][:]
t = (fileobj.variables['time'][0:num_years * days_in_year] - 376956) / 24 # [days]

#****************************************************************
# execute functions (use/remove as necessary)

# general heatwaves
clim_thresh, clim_mean = heatwave_functions.clim_calcs(temperature, num_yearsCLIM, pctile)
heatwaves = heatwave_functions.mhw_metrics(temperature, clim_thresh, clim_mean, lat, lon, t)

# eddyHeatwaves

eddies, eddies_a, eddies_c = eddySupp_functions.eddy_census_calc(eddies_tracked, temperature, lon, lat, clim_mean)
sst_absolute_a, sst_anom_a, amp_a, scale_a, rot_velocity_a = eddySupp_functions.eddy_plotready(eddies_a)
