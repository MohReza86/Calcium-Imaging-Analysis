''' @uthor: Mohammadreza Baghery '''

# Importing required libraries
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import statistics as stat
from scipy.signal import find_peaks, peak_widths, peak_prominences
import csv

path = ''

# Setting constants
std_above_mean = 1          # standard deviations above mean as threshold to detect peaks 
peak_dist = 10000           # Required minimal horizontal distance in samples between neighbouring peaks
avg_window = 5000           # size of moving average box for smoothing
local_min_window = 80000    # size of moving local minima window for detrending
local_avg_window = 50000    # size of moving local average window for baseline detection

# reading the csv file from a specified path
df = pd.read_csv(path, header=1)

# removing rows contaning Not a Number (NaN) from desired columns
df.dropna(subset=['Time(s)', 'AIn-1 - Dem (AOut-1)', 'AIn-2 - Dem (AOut-2)'], inplace=True)

# the time column, getting every other data point
time = df['Time(s)']
time = np.array(time)[::2]

# length of each time stamps in seconds 
point_dur  = time[1] - time[0]

# the GCaMP column; getting every other data point
GCaMP = df['AIn-1 - Dem (AOut-1)']
GCaMP = np.array(GCaMP)[::2]

# the mCherry column; getting every other data point
mCherry = df['AIn-2 - Dem (AOut-2)']
mCherry = np.array(mCherry)[::2]


# subtracting a local minimum value for each time point, 
# defined as the minimum value within a radius of given data points of the smoothed trace
def local_min_win(signal, window_size):

    '''
    signal: array of signal
    window_size: moving local minima window size

    Returns the de-trended trace 
    '''

    local_min = []
    new_signal = []
    i = 0
    while i < len(signal) - window_size + 1:
        this_window = signal[i : i + window_size]
        window_min = min(this_window) 
        local_min.append(window_min)
        for v in this_window:
            new = v - window_min
            new_signal.append(new)
        i += window_size

    return new_signal

detrended_GCaMP = local_min_win(GCaMP, local_min_window)
detrended_mCherry = local_min_win(mCherry, local_min_window)


# revising the time list to match the sigbal list
time_rev = time[0:len(detrended_GCaMP)]


def local_avg_win(signal, window_size):

    '''
    signal: array of signal
    window_size: moving local minima window size

    Returns baseline as the mean of an event-free portion of trace and standard deviation of 
    that event-free portion
    '''

    min_wins = []
    i = 0
    while i < len(signal) - window_size + 1:
        this_window = list(signal[i : i + window_size])
        min_wins.append(this_window)
        i += window_size

    min_win = min(min_wins)
    baseline = stat.mean(min_win)
    #std = stat.stdev(min_win)

    return(baseline)


# baseline as the mean of the event-free portion of trace 
baseline = local_avg_win(detrended_GCaMP, local_min_window)


# dF/F = (F(t) - F0)/F0 
# F0 denotes the minimal averaged fluorescence within 8s windows throughout the measurement
GCaMP_delta = [(i-baseline)/baseline for i in detrended_GCaMP]
mCherry_delta = [(i-baseline)/baseline for i in detrended_mCherry]


# smoothing using a moving average of given points, scaled by a given value
def smooth(y, box_pts):

    '''
    y: array of signal
    box_pts: size of the moving window in data points

    Returns smoothed trace
    '''

    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    scaling_value = 1.2    # scaling after smoothing to match original data
    return (y_smooth * scaling_value) 


smoot_GCaMP = smooth(GCaMP_delta, avg_window)
smoot_mCherry = smooth(mCherry_delta, avg_window)


def standard_deviation(signal, window_size):

    min_wins = []
    i = 0
    while i < len(signal) - window_size + 1:
        this_window = list(signal[i : i + window_size])
        min_wins.append(this_window)
        i += window_size

    min_win = min(min_wins)
    std = stat.stdev(min_win)
    return(std)

# standard deviation of the event-free portion of trace 

std = standard_deviation(smoot_GCaMP, local_min_window)


# removing the first 10 seconds
time_rev = time_rev[60000:]
smoot_GCaMP = smoot_GCaMP[60000:]
smoot_mCherry = smoot_mCherry[60000:]


# converting to numpy array 
GCaMP_np = np.array(smoot_GCaMP)
mCherry_np = np.array(smoot_mCherry)


# Finding peaks of certain height and distance
peak_idx, properties = find_peaks(smoot_GCaMP, 
height=(baseline + (std_above_mean * std)), distance = peak_dist)

peak_idx = peak_idx.tolist()


# coordinates of selected peaks 
peak_x = time_rev[peak_idx]
peak_y = GCaMP_np[peak_idx]


prominences = peak_prominences(smoot_GCaMP, peak_idx)[0]

contour_heights = peak_y - prominences

peak_height = peak_y - contour_heights

peak_height_idx = []
for index, value in enumerate(peak_height):
    if value >= (baseline + (std_above_mean * std)):
        peak_height_idx.append(index)


sel_peak_x = peak_x[peak_height_idx].tolist()
sel_peak_y = peak_y[peak_height_idx].tolist()


sel_peak_idx = np.searchsorted(time_rev, sel_peak_x)

# determining the duration of each peak at its half prominence 
width = peak_widths(smoot_GCaMP, sel_peak_idx, rel_height=0.5)

peak_width = width[0]


# plotting the detrended signal 
plt.plot(time_rev, GCaMP_np)
#plt.plot(time_rev, mCherry_np)

# marking the peaks with x
plt.plot(sel_peak_x, sel_peak_y, "x")

#plt.vlines(x=peak_x, ymin=contour_heights, ymax=peak_y)

# drawing a line marking standard deviations above mean as threshold to detect peaks 
plt.axhline(y= (baseline + (std_above_mean * std)), ls='--', c='g', lw=2)


# limiting the y axis 
plt.ylim(0, 10)

plt.show()


# collecting the duration of the peaks in seconds 
peak_width_sec = []
for i in peak_width:
    dur_sec = i * point_dur
    peak_width_sec.append(round(dur_sec, 2))

# rounding times of peaks
sel_peak_x_round = [i for i in sel_peak_x sel_peak_x_round.append(round(i, 2))]

# rounding amplitude of peaks
sel_peak_y_round = [i for i in sel_peak_y sel_peak_y_round.append(round(i, 2))]


index = range(1, (len(peak_width_sec)+1))


# writing the peak proerties in a csv file in the same folder as the current python script 
header = ['Event', 'Time', 'Amplitude', 'Duration (sec)'] 
filename = "GCaMP events.csv"

with open(filename, 'w') as f:
        wr = csv.writer(f)
        wr.writerow(header)

        for w in range(len(sel_peak_x)):
            wr.writerow([index[w], sel_peak_x_round[w], sel_peak_y_round[w], peak_width_sec[w]])







