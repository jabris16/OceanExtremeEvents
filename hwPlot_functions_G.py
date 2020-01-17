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
    
def area_plot(area_yearly):
    
    '''
    
    Plot cumulative area (km2) experiencing MHW conditions in each yearly bin. 
    
    INPUTS:
    area_yearly = array of cumulative area covered by MHW conditions in yearly bins across the analysis period, calculated in 'heatwave_functions_G.py'

    OUTPUT:
    Plot of cumulative area vs. time [years]
    
    '''

    area_yearly = area_yearly * 1000 # multiply by 1000 to convert to units of (1000s of km2) to make the plot clearer
    years = len(area_yearly) # number of years in analysis period
    
    plt.plot(range(years),area_yearly, 'r', linewidth = 2)
    plt.xlabel('Time since 1 Jan 1993 (years)')
    plt.ylabel('Cumulative area (thousand $km^2$)')
    plt.ticklabel_format(axis = 'y', style = 'sci')
    plt.title('Cumulative area per year experiencing MHW conditions')
    plt.show()

def hw_histogramsComparison(heatwaves, cut_off, num_years):
    
    '''
    
    MHW metric plots for the beginning of the analysis period vs. end [ e.g. first 5 years vs. last 5 years ]. 

    INPUTS:
    heatwaves = dataset of MHW metrics, from 'heatwave_functions_G.py'
    cut_off = number referring to the time length of the subsets for comparison at the start and end of the analysis period [ e.g. cut_off = 5 then first 5 years and last 5 years ] [years] 
    num_years = total number of years in analysis period, must be integer [years]

    OUTPUTS:
    Plots of MHW metrics
    'Duration' = histograms of MHW duration distribution for start and end of analysis period
    'Max Intensity' = histograms of maximum SST anomaly (relative to MHW threshold) during MHW events distribution for start and end of analysis period
    'Category' = histograms of MHW category (categories based on Hobday et al., 2016, Oceanography) distribution for the start and end of analysis period 
    
    '''

    days_in_year = 365
    duration_start = []
    duration_end = []
    max_intensity_start = []
    max_intensity_end = []
    category_start = []
    category_end = []

    for grid in range(len(heatwaves)):
        for ev in range(len(heatwaves[grid]['duration'])):
            if heatwaves[grid]['time_start'][ev] < cut_off * days_in_year:
                duration_start.append(heatwaves[grid]['duration'][ev])
                max_intensity_start.append(heatwaves[grid]['intensity_max_relThresh'][ev])
                category_start.append(heatwaves[grid]['category'][ev])
            elif heatwaves[grid]['time_start'][ev] > (num_years * days_in_year) - (cut_off * days_in_year):
                duration_end.append(heatwaves[grid]['duration'][ev])
                max_intensity_end.append(heatwaves[grid]['intensity_max_relThresh'][ev])
                category_end.append(heatwaves[grid]['category'][ev])
            else:
                pass

    # set up axes
    
    fig, axs = plt.subplots(2, 3, sharey = False, tight_layout = True)

    # a) duration
    
    axs[0].hist(duration_start, weights = np.ones(len(duration_start)) / len(duration_start) * 100, histtype = 'stepfilled', color = 'r')
    axs[0].set_xlabel('Duration (days)')
    axs[0].set_ylabel('% MHW events')
    axs[0].set_title('A) Duration Distribution', loc = 'left') 
    
    # b) max intensity
    
    axs[1].hist(max_intensity_start, weights = np.ones(len(max_intensity_start)) / len(max_intensity_start) * 100, histtype = 'stepfilled', color = 'r') 
    axs[1].set_xlabel('Maximum Intensity ($^\circ$C)')
    axs[1].set_ylabel('% MHW events')
    axs[1].set_title('B) Maximum Intensity Distribution', loc = 'left')
    
    # c) category distribution
    
    axs[2].hist(category_start, weights = np.ones(len(category_start)) / len(category_start) * 100, histtype = 'stepfilled', color = 'r')
    axs[2].set_xlabel('MHW Category')
    axs[2].set_ylabel('% MHW events')
    axs[2].set_title('C) Category Distribution', loc = 'left')

