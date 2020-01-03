'''

IN PROGRESS.

General heatwave plots.

Jamie Atkins

'''
#****************************************************************

import numpy as np
from matplotlib import pyplot as plt

#****************************************************************

# FUNCTIONS

def hw_histograms(heatwaves, years):
    
    '''
    MHW metric plots. Using output from 'heatwave_functions_G.py'.
    
    INPUT:
    heatwaves = dataset of MHW metrics, from 'heatwave_functions_G.py'
    years = number of years in the analysis period
    
    OUTPUTS:
    Plots of MHW metrics
    'Duration' = histogram of MHW duration distribution across the analysis period
    'Max Intensity' = histogram of maximum SST anomaly (relative to MHW threshold) during MHW events distribution across the analysis period
    'Category' = histogram of MHW category distribution across the analysis period (categories based on Hobday et al., 2016, Oceanography)
    'Incidence' = MHW count in each year of the analysis period (number of events counted in yearly bins)
    
    '''
    
    duration = []
    max_intensity = []
    time = []
    category = []
        
    # create list version of different metrics
    for grid in range(len(heatwaves)):
        for ev in range(len(heatwaves[grid]['duration'])):
            duration.append(heatwaves[grid]['duration'][ev])
            max_intensity.append(heatwaves[grid]['intensity_max_relThresh'][ev])
            time.append(heatwaves[grid]['time_start'][ev] / 365)
            category.append(heatwaves[grid]['category'][ev])
    
    # set up figure
    fig, axs = plt.subplots(1, 4, sharey = False, tight_layout = True)

    # a) duration
    
    axs[0].hist(duration, weights = np.ones(len(duration)) / len(duration) * 100, histtype = 'stepfilled', color = 'r')
    axs[0].set_xlabel('Duration (days)')
    axs[0].set_ylabel('% MHW events')
    axs[0].set_title('A) Duration Distribution', loc = 'left') 
    
    # b) max intensity
    
    axs[1].hist(max_intensity, weights = np.ones(len(max_intensity)) / len(max_intensity) * 100, histtype = 'stepfilled', color = 'r') 
    axs[1].set_xlabel('Maximum Intensity ($^\circ$C)')
    axs[1].set_ylabel('% MHW events')
    axs[1].set_title('B) Maximum Intensity Distribution', loc = 'left')
    
   
    # c) category distribution
    
    axs[2].hist(category, weights = np.ones(len(category)) / len(category) * 100, histtype = 'stepfilled', color = 'r')
    axs[2].set_xlabel('MHW Category')
    axs[2].set_ylabel('% MHW events')
    axs[2].set_title('C) Category Distribution', loc = 'left')
    
    # d) incidence (yearly bins)
    
    n, x, _ = plt.hist(time, bins = np.linspace(0, years, years), histtype = 'stepfilled', color = 'w')
    bin_centers = 0.5 * (x[1:] + x[:-1])
    axs[3].plot(bin_centers, n, color = 'r')
    axs[3].set_xlabel('Time since 1 Jan 1993 (years)')
    axs[3].set_ylabel('Number of MHW events')
    axs[3].set_title('D) Incidence over time', loc = 'left')
