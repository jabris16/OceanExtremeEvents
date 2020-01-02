'''

STILL IN PROGRESS.

General heatwave plots.

Jamie Atkins

'''
#****************************************************************

import numpy as np
from matplotlib import pyplot as plt

#****************************************************************

# FUNCTIONS

def hw_histograms(heatwaves):
    
    duration = []
    max_intensity = []
    time = []
    category = []
        
    # create list version of different metrics
    for grid in range(len(heatwaves)):
        for ev in range(len(heatwaves[grid]['duration'])):
            duration.append(heatwaves[grid]['duration'][ev])
            max_intensity.append(heatwaves[grid]['intensity_max_relThresh'][ev])
            time.append(heatwaves[grid]['time_start'][ev])
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
    axs[2].set_ylabel('Number of MHW events')
    axs[2].set_title('C) Category Distribution', loc = 'left')
     
    # d) incidence
    
    axs[3].hist(time, histtype = 'step', color = 'r')
    axs[3].set_xlabel('Time (days)')
    axs[3].set_ylabel('Number of MHW events')
    axs[3].set_title('D) Incidence over time', loc = 'left')

