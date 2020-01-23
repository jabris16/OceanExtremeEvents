'''

Functions for eddy-specific heatwave analysis.

Jamie Atkins

'''

#****************************************************************

import numpy as np
import math

#****************************************************************

###################

### EDDY CENSUS ###

###################

def eddy_census_calc(eddies_tracked, field, lon, lat, clim_mean, spat_res):
    
    ''' 
    *** companion function to Eric Oliver eddyTracking; requires output from eddyTracking ***
       
    Calculate and store  additional eddy metrics to supplement those from eddyTracking.
 
    INPUT:
    eddies_tracked = array of tracked eddies from eddyTracking (N.B. tracked eddies must come from a temporal period that does not exceed the temporal period of 'field', otherwise must subset 'eddies_tracked to match 'field')
    field = SST dataset in array format of shape e.g. temperature[time,lat,lon] where time is a multiple of 365 (i.e. complete years)
    lat = array of latitude values corresponding to 'field'
    lon = array of longitude values corresponding to 'field'
    clim_mean = mean SST value for each grid cell across each day of the year
    spat_res = model spatial resolution (adjust as necessary) [degrees]; ideally 0.25 (eddy-permitting) or higher
    
    OUTPUT:
    eddies =  all metrics
    eddies_a = subset of 'eddies' with only anticyclonic eddies
    eddies_c = subset of 'eddies' with only cyclonic eddies
    
        METRICS:
        'time' = time index of eddy feature
        'type' = anticyclonic or cyclonic eddy
        'age' = total age of eddy [days]
        'sst_absolute' = absolute eddy grid cell SST value [deg C]
        'sst_anom' = eddy grid cell SST anomaly (absolute SST - mean for grid cell for each day of the year) [deg C]
        'amp' = eddy amplitude [cm]
        'scale' = eddy radius [km]
        'rot_velocity' = rotational velocity of eddy feature [cm/s]
        'translation_speed' = translation speed of eddy feature [cm/s] (age must exceed 1 day; no value for first index as translation speed (c) must be calculated across subsequent eddy features/time steps)
        'nonlin' = nonlinerity value of eddy feature (age must exceed 1 day, no value for first index as c must be calculated across subsequent eddy features/time steps)
        'lat' = latitude of eddy feature [deg]
        'lon' = longitude of eddy feature [deg]
        'lat_index' = array index of latitude
        'lon_index' = array index of longitude
            
    '''

    print('Starting eddy_census_calc()...')

    res_factor = int(1 / spat_res)
    grav = 981 # gravitational constant [cm/s]
    earth_radius = 637813700 # [cm]
    earth_rot = 0.000072921 # Earth rotation rate [rad/s]
    time_24hrs = 86400 # number of seconds in a day [seconds]
    days_in_year = 365 # number of days in a year [days]
    num_years = int(len(field) / days_in_year) # number of years in analysis period [years]
    radians_convert = np.pi/180 # radians to degrees conversion factor
    
    def gridcell_round(x):
    
        ''' function to round/assign an eddy feature to a [spat_res] degree grid cell '''
        
        return math.floor(x * res_factor) / res_factor
    
    eddies = []

    for ed in range(len(eddies_tracked)):
        eddies_grid = {}
        eddies_grid['lon_index'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lat_index'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lon'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['lat'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['time'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['sst'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['age'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['type'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['amp'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['scale'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['rot_velocity'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['translation_speed'] = np.zeros(len(eddies_tracked[ed]['lon']))
        eddies_grid['nonlin'] = np.zeros(len(eddies_tracked[ed]['lon']))
        
        for x in range(len(eddies_tracked[ed]['lon'])):
            # index values of lat/lon [for purposes of eddy-SST assigning]
            eddies_grid['lon_index'][x] = list(lon).index(gridcell_round(eddies_tracked[ed]['lon'][x]))
            eddies_grid['lat_index'][x] = list(lat).index(gridcell_round(eddies_tracked[ed]['lat'][x]))
            eddies_grid['lat'][x] = eddies_tracked[ed]['lat'][x]
            eddies_grid['lon'][x] = eddies_tracked[ed]['lon'][x]

            # standard metrics from eddyTracking 
            eddies_grid['time'][x] = eddies_tracked[ed]['time'][x]
            eddies_grid['age'][x] = eddies_tracked[ed]['age']
            eddies_grid['amp'][x] = eddies_tracked[ed]['amp'][x]
            eddies_grid['scale'][x] = eddies_tracked[ed]['scale'][x]

            # SST values
            eddies_grid['sst'][x] = field[int(eddies_grid['time'][x] - 1),int(eddies_grid['lat_index'][x]),int(eddies_grid['lon_index'][x])] # time -1 to correspond to timestep in 'field'
            
            # assign eddy type; 0 == cyclonic, 1 == anticyclonic
            if eddies_tracked[ed]['type'] == 'anticyclonic':
                eddies_grid['type'][x] = 1
            elif eddies_tracked[ed]['type'] == 'cyclonic':
                eddies_grid['type'][x] = 0
        
            # rotational speed (U) calculation [useful metric + component of nonlinearity calculation]
            coriolis = 2 * earth_rot * np.sin(eddies_tracked[ed]['lat'][x] * radians_convert)
            scalar = grav / coriolis
            amp_scale = eddies_tracked[ed]['amp'][x] / (eddies_tracked[ed]['scale'][x] * 100000)
            rot_velocity = scalar * amp_scale * 100
            eddies_grid['rot_velocity'][x] = rot_velocity
            
        for x in range(1,len(eddies_tracked[ed]['lon'])):
            # translation speed (c) calculation  [useful metric + component of nonlinearity calculation]
            if eddies_tracked[ed]['age'] >= 2: # the eddy must have at least 2 track points to calc c
                angles = np.sin(eddies_tracked[ed]['lat'][x]) * np.sin(eddies_tracked[ed]['lat'][x-1]) + \
                    np.cos(eddies_tracked[ed]['lon'][x] - eddies_tracked[ed]['lon'][x-1]) * \
                    np.cos(eddies_tracked[ed]['lat'][x]) * np.cos(eddies_tracked[ed]['lat'][x-1])
                distance = earth_radius * np.arccos(angles) * radians_convert
                speed = distance / time_24hrs
                eddies_grid['translation_speed'][x] = speed
                
            # calculate nonlinearity (U/c)
                nonlin = rot_velocity / speed
                eddies_grid['nonlin'][x] = nonlin
                    
            else:
                pass

        eddies.append(eddies_grid)

    eddies_a = []
    eddies_c = []
    
    baseline_mean_int = [clim_mean,] * num_years
    baseline_mean = np.reshape(baseline_mean_int, (num_years * days_in_year, len(lat), len(lon)))
    
    # further split into anticylonic and cyclonic datasets
    for ed in range(len(eddies)):
        eddy_int = {}
        eddy_int['time'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['sst_absolute'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['sst_anomaly'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['amp'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['scale'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['rot_velocity'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['lat'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['lon'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['rot_velocity'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['translation_speed'] = np.zeros(len(eddies[ed]['time']))
        eddy_int['nonlin'] = np.zeros(len(eddies[ed]['time']))
        if eddies[ed]['type'][0] == 1: # anticyclonic
            for x in range(len(eddies[ed]['time'])):
                eddy_int['time'][x] = eddies[ed]['time'][x]
                eddy_int['sst_absolute'][x] = eddies[ed]['sst'][x]
                eddy_int['sst_anomaly'][x] = eddies[ed]['sst'][x] - baseline_mean[int(eddies[ed]['time'][x]),int(eddies[ed]['lat_index'][x]),int(eddies[ed]['lon_index'][x])]
                eddy_int['amp'][x] = eddies[ed]['amp'][x] * 100 # convert to cm
                eddy_int['scale'][x] = eddies[ed]['scale'][x]
                eddy_int['rot_velocity'][x] = eddies[ed]['rot_velocity'][x]
                eddy_int['lat'][x] = eddies[ed]['lat'][x]
                eddy_int['lon'][x] = eddies[ed]['lon'][x]
                eddy_int['rot_velocity'][x] = eddies[ed]['rot_velocity'][x]
                eddy_int['translation_speed'][x] = eddies[ed]['translation_speed'][x]
                eddy_int['nonlin'][x] =  eddies[ed]['nonlin'][x]
            eddies_a.append(eddy_int)
        elif eddies[ed]['type'][0] == 0: # cyclonic
            for x in range(len(eddies[ed]['time'])):
                eddy_int['time'][x] = eddies[ed]['time'][x]
                eddy_int['sst_absolute'][x] = eddies[ed]['sst'][x]
                eddy_int['sst_anomaly'][x] = eddies[ed]['sst'][x] - baseline_mean[int(eddies[ed]['time'][x]),int(eddies[ed]['lat_index'][x]),int(eddies[ed]['lon_index'][x])]
                eddy_int['amp'][x] = eddies[ed]['amp'][x] * 100 # to convert to cm
                eddy_int['scale'][x] = eddies[ed]['scale'][x]
                eddy_int['rot_velocity'][x] = eddies[ed]['rot_velocity'][x]
                eddy_int['lat'][x] = eddies[ed]['lat'][x]
                eddy_int['lon'][x] = eddies[ed]['lon'][x]
                eddy_int['rot_velocity'][x] = eddies[ed]['rot_velocity'][x]
                eddy_int['translation_speed'][x] = eddies[ed]['translation_speed'][x]
                eddy_int['nonlin'][x] =  eddies[ed]['nonlin'][x]
            eddies_c.append(eddy_int)
    
    return eddies, eddies_a, eddies_c

###################

### EDDY PLOTS ###

###################

def eddy_plotready(eddyCHOICE):
    
    ''' 
    
    Create lists ready to be easily plotted using matplotlib.
        
    INPUT:
    eddyCHOICE = eddies_a or eddies_c from eddy_census_calc
    
    OUTPUT:
    lists of absolute SST, SST anomaly, Eddy amplitude, Eddy scale/radius and Rotational velocity
    
    '''
    
    print('Starting eddy_plotready():')

    sst_absolute = []
    sst_anom = []
    amp = []
    scale = []
    rot_velocity = []
    
    for ed in range(len(eddyCHOICE)):
        for x in range(len(eddyCHOICE[ed]['time'])):
            sst_absolute.append(float(eddyCHOICE[ed]['sst_absolute'][x]))
            sst_anom.append(float(eddyCHOICE[ed]['sst_anomaly'][x]))
            amp.append(float(eddyCHOICE[ed]['amp'][x]))
            scale.append(float(eddyCHOICE[ed]['scale'][x]))
            rot_velocity.append(float(eddyCHOICE[ed]['rot_velocity'][x]))

    return sst_absolute, sst_anom, amp, scale, rot_velocity
