'''

Companion file.

Set parameters, load eddyTracking and heatwave data.

Jamie Atkins

'''

#****************************************************************

# LIBRARIES

import numpy as np
from netCDF4 import Dataset
import heatwave_functions_G
import eddySuppli_functions_G

#****************************************************************

# PARAMETERS (adjust as necessary)

pctile = 90 # percentile threshold of choice (90% in line with Hobday et al., 2016)
num_yearsCLIM = 5 # time length [years] of the climatological mean/threshold period
resolution = 0.25 # model spatial resolution
num_years = 23 # number of years of analysis
days_in_year = 365 # number of days in a year

path = './' # working directory route 
run = 'MetO' # name of run, should match that from eddyTracking
nc_filename = 'dataset.nc' # filename of temperature NetCDF file
nc_filename_gridarea = 'gridarea.nc' # filename of grid cell area NetCDF file

kelvin_convert = 273.15 # Kelvin to Celcius conversion factor; switch off if temperature dataset is in degrees C already

#****************************************************************

# DATA LOAD

# eddyTracking output (from Eric Oliver code)
data = np.load(path + 'eddy_track_' + run + '.npz') # 'encoding' may be required if switching between Python 2. and Python 3.
eddies_tracked = data['eddies']

# temperature data
z = path + nc_filename
fileobj = Dataset(z)
temperature = fileobj.variables['temperature'][0:num_years * days_in_year,0,:,:]  - kelvin_convert
lon = fileobj.variables['longitude'][:]
lat = fileobj.variables['latitude'][:]
t = (fileobj.variables['time'][0:num_years * days_in_year] - 376956) / 24 # [days]

# grid cell area
### N.B. CDO 'gridarea' function is useful for creating a NetCDF file of grid cell areas
fileobj = Dataset(nc_filename_gridarea)
grid_area = fileobj.variables['cell_area'][:]

#****************************************************************
# execute functions (use/remove as necessary)

# general heatwaves
clim_thresh, clim_mean = heatwave_functions_G.clim_calcs(temperature, num_yearsCLIM, pctile)
heatwaves = heatwave_functions_G.mhw_metrics(temperature, clim_thresh, clim_mean, lat, lon, t)
area, area_yearly = heatwave_functions_G.mhw_area(heatwaves, grid_area, num_years, t)

# eddyHeatwaves

eddies, eddies_a, eddies_c = eddySuppli_functions_G.eddy_census_calc(eddies_tracked, temperature, lon, lat, clim_mean, resolution)
sst_absolute_a, sst_anom_a, amp_a, scale_a, rot_velocity_a = eddySuppli_functions_G.eddy_plotready(eddies_a)
