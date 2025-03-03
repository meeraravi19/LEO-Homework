import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from leo_hw1 import calc_orbital_elems

# n = 2460409.06736111 - 2451545.0
# # print(n)

# L = (280.460 + (0.9857464 * n) )% 360
# #print(L)

# g = (357.526 + (0.9856003 * n)) % 360
# #print(g)

# lamb = L + (1.915 * math.degrees(math.sin(math.radians(g)))) + (0.020* math.degrees(math.sin(math.radians(2 * g))))
# #print("lambda: ", lamb)

# eps = 23.439 - (0.0000004 * n)
# #print("epsilon: ", eps)

# obiq = 23.439 - (0.0000004 * n)
# f = 180 / math.pi # to convert to degrees
# t = math.tan(eps / 2) ** 2

# alpha = lamb - (f * t * math.sin(math.radians(2 * lamb))) + ((f / 2) * (t ** 2) * math.sin(math.radians(4 * lamb)))
# # calculates in radians
# alpha_2 = math.atan(math.cos(math.radians(eps))* math.tan(math.radians(lamb)))

# delt = math.degrees(math.asin(math.sin(math.radians(eps)) * math.sin(math.radians(lamb))))

# Homework 3, Problem 2
file_path = "C:/Users/meera/Documents/Desktop/UT Junior Year (24-25)/Spring Semester/LEO/1yr.180s.ECI.txt"
data = pd.read_csv(file_path, delim_whitespace=True, header=None)

time = data[0].values  # Time in seconds or Julian Date
r = data.iloc[:, 2:5].values  # Position vectors
v = data.iloc[:, 5:8].values

def detect_revolutions(time, r):
    z = r[:, 2]  # Z-component (height above equator)
    crossings = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]  

    rev_times = time[crossings]
    return crossings, rev_times

rev_indices, rev_times = detect_revolutions(time, r)
#print(rev_indices)

mu = 398600.44 * 1000

osculating_elements = np.array([calc_orbital_elems(r[i] / 1000, v[i] / 1000, 398600.44) for i in range(len(time))])
#print(osculating_elements.size)

mean_elements = []
mean_times = []

for i in range(len(rev_indices) - 1):
    start, end = rev_indices[i], rev_indices[i + 1]
    rev_oe = osculating_elements[start:end]

    mean_oe = np.mean(rev_oe, axis=0)
    mean_elements.append(mean_oe)
    mean_times.append(np.mean(time[start:end]))

mean_elements = np.array(mean_elements)
mean_times = np.array(mean_times)

start_index = 0
end_index = np.searchsorted(time, time[0] + (86400))
time_day = time[start_index:end_index]
osculating_day = osculating_elements[start_index:end_index]
mean_mask = mean_times < time[0] + (86400)
mean_times_day = mean_times[mean_mask] / 60  # Convert to minutes
mean_elements_day = mean_elements[mean_mask]

# first partial
first_partials = []
first_partial_time = (180 * rev_indices[0]) / 2
for i in range(osculating_elements.shape[1]):
    first_partial = osculating_elements[0:rev_indices[0], i]
    first_partial_mean = np.mean(first_partial)
    print(first_partial_mean)
    first_partials.append(first_partial_mean)
#first_partial = osculating_elements[0:rev_indices[0]]
#first_partial_mean = np.mean(first_partial)


# last partial
this = 0
for i, val in enumerate(rev_indices):
    if val < end_index:
        this = i

last_partial_time = np.array(np.mean(time[rev_indices[this]:end_index]))
last_partials = []
for i in range(osculating_elements.shape[1]):
    last_partial = osculating_elements[rev_indices[this]:end_index, i]
    last_partial_mean = np.array(np.mean(last_partial))
    last_partials.append(last_partial_mean)

mean_times_day = np.insert(mean_times_day, 0, first_partial_time / 60)
mean_times_day = np.insert(mean_times_day, len(mean_times_day), 85500 / 60)
mean_elements_day = np.insert(mean_elements_day, 0, first_partials, axis=0)
mean_elements_day = np.insert(mean_elements_day, len(mean_elements_day), last_partials, axis=0)

def plot_elements(time, oe, mean_times, mean_oe, title, ylab, day):
    plt.figure(figsize=(10, 5))
    if day == True:
        plt.plot(time / 60, oe, label="Osculating", linewidth = 0.5, linestyle="-", marker="o", markersize = 3)
        plt.plot(mean_times, mean_oe, 'o-', label="Mean", color='red')
        plt.xlabel("Time (min)")
    else:
        plt.plot(mean_times / 86400, mean_oe, 'o-', label="Mean", color='red', markersize = 3.5)
        plt.xlabel("Time (days since July 31, 2023)")

    plt.ylabel(ylab)
    plt.title(title)
    plt.legend(loc = "upper right")
    plt.show()

ylabels = ["a (km)", "e (-)", "i (deg)", "RAAN (deg)", "w (deg)", "M (deg)"]

for i, title in enumerate(["Semi-Major Axis", "Eccentricity", "Inclination", "RAAN", "Arg of Perigee", "Mean Anomaly"]):
    plot_elements(time, osculating_elements[:, i], mean_times, mean_elements[:, i], title, ylabels[i], False)

for i, title in enumerate(["Semi-Major Axis", "Eccentricity", "Inclination", "RAAN", "Arg of Perigee", "Mean Anomaly"]):
    plot_elements(time_day, osculating_day[:, i], mean_times_day, mean_elements_day[:, i], title, ylabels[i], True)