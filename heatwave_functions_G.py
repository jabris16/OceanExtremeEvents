'''

Functions for general heatwave analysis.

Jamie Atkins

'''

#****************************************************************

import numpy as np
import scipy.ndimage as ndimage

#****************************************************************

###############################################

### MEAN/THRESHOLD CLIMATOLOGY CALCULATIONS ###

###############################################

def clim_calcs(field, num_yearsCLIM, pctile):
       
    ''' 

    Heatwave analysis.
    
    function to calculate the mean and specified threshold value for each grid cell across each day of a year over a given climatological period, where: 
       
    INPUT:
    field = daily SST dataset in array format of shape e.g. temperature[time,lat,lon] where time is a multiple of 365 (i.e. complete years)
    num_yearsCLIM = time length [years] of the climatological mean/threshold period, first x years of the 'field' dataset for example
    pctile = specified percentile to compute
    
    OUTPUT:
    clim_mean = 3D array (of shape [365, len(lat), len(lon)] of mean SST for each grid cell averaged for each day of the year across the climatological period
    clim_thresh = 3D array (of shape [365, len(lat), len(lon)] of threshold percentile SST for each grid cell calulated for each day of the year across the climatological period
    
    *** TO FOLLOW: UPDATE TO PROCESS LEAP YEARS, I.E. SO 'FIELD' DOES NOT NEED TO BE A MULTIPLE OF 365
        UPDATE TO INCLUDE OPTION USE WEEKLY, MONTHLY DATA AS WELL AS DAILY ***
        
        
    '''
    
    print('Starting clim_calcs():')
    
    days_in_year = 365
    array_shape = np.shape(field[0:365,:,:])
    t = num_yearsCLIM * days_in_year 
    clim_thresh = np.zeros(array_shape)
    clim_mean = np.zeros(array_shape)
    t_len = np.array(range(t))
    t_store = np.zeros((days_in_year, int(t / days_in_year)), dtype = object)

    # subset timesteps for upcoming calculations
    for i in range(days_in_year):
        for j in range(num_yearsCLIM):
            if i == 0:
                t_store[0][j] = t_len[j * days_in_year]
            else:
                t_store[i][:] = t_store[0] + int(i)

    # calculate percentile and mean values
    for x in range(len(field[0,:,0])):
        print(str(x) + ' of ' + str(len(field[0,:,0])-1))
        for y in range(len(field[0,0,:])):
            for i in range(days_in_year):
                sst_int = []
                for j in range(num_yearsCLIM):
                    sst_int.append(field[t_store[i][j],x,y])
                sst_array = np.array(sst_int)
                sst_thresh = np.nanpercentile(sst_array, pctile)
                sst_mean = np.nanmean(sst_array)
                clim_thresh[i,x,y] = sst_thresh
                clim_mean[i,x,y] = sst_mean
    
    return clim_thresh, clim_mean

###############################

### MARINE HEATWAVE METRICS ###

###############################

def mhw_metrics(field, clim_thresh, clim_mean, lat, lon, t):

    '''

    Calculating Marine Heatwave Metrics.

    Using some functions that are based on those by Eric Oliver in marineHeatwaves (https://github.com/ecjoliver/marineHeatWaves/blob/master/marineHeatWaves.py)
    Expands on Eric Oliver code by applying the methods to all grid cells within a given climatological region.

    INPUT:
    field = SST dataset in array format of shape e.g. temperature[time,lat,lon] where time is a multiple of 365 (i.e. complete years)
    baseline_thresh = SST array of threshold values at a certain percentile for each time step and grid cell , calculated in clim_calc function
    baseline_mean = as above but mean values
    lat = array of latitude values corresponding to 'field'
    lon = array of longitude values corresponding to 'field'
    t = array of time values corresponding to 'field' 
    
    OUTPUT:
    heatwaves = dataset of marine heatwave metrics
    *** For more detailed outputs information see Eric Oliver marineHeatwaves 'outputs' (https://github.com/ecjoliver/marineHeatWaves/blob/master/marineHeatWaves.py) 
    
    '''
    
    print('Starting mhw_metrics():')

    
    heatwaves = [] # create empty heatwaves list
    days_in_year = 365 # number of days in a year
    num_years = int(len(field) / days_in_year) # number of years in analysis period
    # repeat clim_mean for number of years in climatological period so compatible with 'field'
    baseline_mean_int = [clim_mean,] * num_years
    baseline_mean = np.reshape(baseline_mean_int, (num_years * days_in_year, len(lat), len(lon))) 
    # repeat clim_thresh for number of years in climatological period so compatible with 'field'
    baseline_thresh_int = [clim_thresh,] * num_years
    baseline_thresh = np.reshape(baseline_thresh_int, (num_years * days_in_year, len(lat), len(lon)))

    for x in range(len(lat)):
        print(str(x) + ' of ' + str(len(lat)-1))
        for y in range(len(lon)):
            mhw = {}
            mhw['lat'] = []
            mhw['lon'] = []
            mhw['lat_index'] = []
            mhw['lon_index'] = []
            mhw['time_start'] = [] # [index]
            mhw['time_end'] = [] # [index]
            mhw['n_events'] = []
            mhw['time_peak'] = [] # [index]
            mhw['duration'] = [] # [days]
            mhw['duration_moderate'] = [] # [days]
            mhw['duration_strong'] = [] # [days]
            mhw['duration_severe'] = [] # [days]
            mhw['duration_extreme'] = [] # [days]
            mhw['intensity_max'] = [] # [deg C]
            mhw['intensity_mean'] = [] # [deg C]
            mhw['intensity_var'] = [] # [deg C]
            mhw['intensity_cumulative'] = [] # [deg C]
            mhw['intensity_max_relThresh'] = [] # [deg C]
            mhw['intensity_mean_relThresh'] = [] # [deg C]
            mhw['intensity_var_relThresh'] = [] # [deg C]
            mhw['intensity_cumulative_relThresh'] = [] # [deg C]
            mhw['intensity_max_abs'] = [] # [deg C]
            mhw['intensity_mean_abs'] = [] # [deg C]
            mhw['intensity_var_abs'] = [] # [deg C]
            mhw['intensity_cumulative_abs'] = [] # [deg C]
            mhw['category'] = []
            mhw['rate_onset'] = [] # [deg C / day]
            mhw['rate_decline'] = [] # [deg C / day]

            # find where temp exceeds threshold
            exceed = field[:,x,y] - baseline_thresh[:,x,y]
            exceed[exceed>0] = True
            exceed[exceed<=0] = False
            # label the events that exceed threshold
            events, n_events = ndimage.label(exceed)

            # find all heatwaves where duration exceeds the min_duration
            min_duration = 5 # in line with Hobday et al. (2016)
            for ev in range(1,n_events+1):
                event_duration = (events == ev).sum()
                if event_duration < min_duration:
                    continue
                mhw['time_start'].append(t[np.where(events == ev)[0][0]])
                mhw['time_end'].append(t[np.where(events == ev)[0][-1]])

            # calculate heatwave metrics
            mhw['n_events'] = len(mhw['time_start'])
            categories = np.array(['Moderate', 'Strong', 'Severe', 'Extreme'])
            for ev in range(mhw['n_events']):
                tt_start = np.where(t==mhw['time_start'][ev])[0][0]
                tt_end = np.where(t==mhw['time_end'][ev])[0][0]
                temp_mhw = field[tt_start:tt_end+1,x,y]
                thresh_mhw = baseline_thresh[tt_start:tt_end+1,x,y]
                seas_mhw = baseline_mean[tt_start:tt_end+1,x,y]
                mhw_relSeas = temp_mhw - seas_mhw
                mhw_relThresh = temp_mhw - thresh_mhw
                mhw_relThreshNorm = (temp_mhw - thresh_mhw) / (thresh_mhw - seas_mhw)
                mhw_abs = temp_mhw
                # lat and lon
                mhw['lat_index'] = x
                mhw['lon_index'] = y
                mhw['lat'] = lat[x]
                mhw['lon'] = lon[y]
                # Find peak
                tt_peak = np.argmax(mhw_relSeas)
                mhw['time_peak'].append(tt_start + tt_peak)
                # MHW Duration
                mhw['duration'].append(len(mhw_relSeas))
                # MHW Intensity metrics
                mhw['intensity_max'].append(mhw_relSeas[tt_peak])
                mhw['intensity_mean'].append(mhw_relSeas.mean())
                mhw['intensity_var'].append(np.sqrt(mhw_relSeas.var()))
                mhw['intensity_cumulative'].append(mhw_relSeas.sum())
                mhw['intensity_max_relThresh'].append(mhw_relThresh[tt_peak])
                mhw['intensity_mean_relThresh'].append(mhw_relThresh.mean())
                mhw['intensity_var_relThresh'].append(np.sqrt(mhw_relThresh.var()))
                mhw['intensity_cumulative_relThresh'].append(mhw_relThresh.sum())
                mhw['intensity_max_abs'].append(mhw_abs[tt_peak])
                mhw['intensity_mean_abs'].append(mhw_abs.mean())
                mhw['intensity_var_abs'].append(np.sqrt(mhw_abs.var()))
                mhw['intensity_cumulative_abs'].append(mhw_abs.sum())
                # Fix categories
                tt_peakCat = np.argmax(mhw_relThreshNorm)
                cats = np.floor(1. + mhw_relThreshNorm)
                mhw['category'].append(categories[np.min([cats[tt_peakCat], 4]).astype(int) - 1])
                mhw['duration_moderate'].append(np.sum(cats == 1.))
                mhw['duration_strong'].append(np.sum(cats == 2.))
                mhw['duration_severe'].append(np.sum(cats == 3.))
                mhw['duration_extreme'].append(np.sum(cats >= 4.))

            # ajoin to overall large dataset
            heatwaves.append(mhw)

    return heatwaves

###############################

 ### CUMULATIVE MHW AREA ###

###############################
 
def mhw_area(heatwaves, grid_area, years, t):

    '''
    
    Function for calculating the total cumulative area (km2) experiencing MHW conditions per year (MHW threshold as defined by Hobday et al. 2016).

    INPUTS:
    heatwaves = dataset of MHW metrics, from 'heatwave_functions_G.py'
    grid_area = NetCDF file of the area (km2) of each grid cell within the analysis domain the format (lat, lon)
        N.B. this can be prduced with CDO 'gridarea' function and an existing NetCDF file containing data across the analysis domain
    years = time length of analysis period [years]
    t = array of timesteps in analysis period
    
    OUPUTS:
    area = array of total area (km2) covered by MHW conditions at each timestep across the analysis period
    area_yearly = array of cumulative area covered by MHW conditions in yearly bins across the analysis period

    '''
    
    # CONSTANTS and LISTS
    m2_to_km2_convert = 1000000 # conversion factor between m2 and km2
    days_in_year = 365 # number of days in a year
    t_len = years * 365
    grid_timestore = []
    area_store = []
    time_area = []
    area_sum = []
    area_yearly = []
    
    # CALCULATE MHW AREA AT EACH TIME STEP
    
    # assign the timestep of each heatwave incidence to a list
    for grid in range(len(heatwaves)):
        time_store = []
        for ev in range(len(heatwaves[grid]['duration'])):
            time_store.append(range(int(heatwaves[grid]['time_start'][ev]),int(heatwaves[grid]['time_end'][ev])))
        grid_timestore.append(time_store)
    
    # assign the area of each grid cell to a list in format to match len(heatwave) which is total number of cells in analysis region
    for x in range(len(grid_area[:,0])):
        for y in range(len(grid_area[0,:])):
            area_store.append(grid_area[x,y])
    
    # assign the area of each grid cell which is experiencing heatwave conditions at a given timestep to a list of len(total number of timesteps) to the index of the timestep in which the conditions are occuring
    # N.B. this runs quite slow...
    for i in range(t_len):
        print(str(i) + ' of ' + str(t_len - 1))
        areas = []
        for grid in range(len(grid_timestore)):
            for ev in range(len(grid_timestore[grid])):
                for x in range(len(grid_timestore[grid][ev])):
                    if grid_timestore[grid][ev][x] == t[i]:
                        areas.append(area_store[grid])
                    else:
                        pass
        time_area.append(areas)
    
    # sum the areas of each grid cell experiencing heatwave conditions at a given timestep
    # N.B. each index of area_sum corresponds to a time step, i.e area_sum[0] = 1st timestep
    for i in range(t_len):
        area_add = np.sum(time_area[i])
        area_sum.append(area_add)
    
    area_km = area_sum / m2_to_km2_convert # convert to km2

    # calculate yearly cumulative area (i.e. cumulative area in yearly bins)
    for i in range(years):
        if i == 0:
            area_yearly.append(np.sum(area_km[i * days_in_year:days_in_year - 1]))
        else:
            area_yearly.append(np.sum(area_km[i * days_in_year:i * (days_in_year -1) + days_in_year]))
    area_yearly_np = np.array(area_yearly)
    area_yearly = area_yearly_np

    # SAVE OUTPUT
    # each timestep
    area = np.array(area_km)
    np.savez('area_timestep', area = area)  
    # yearly bins
    np.savez('area_yearly', area_yearly = area_yearly)
    
    return area, area_yearly
