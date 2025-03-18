import math
import numpy as np
from leo_hw3 import find_sun_coor
import matplotlib.pyplot as plt
from leo_hw2 import calc_gmst
from datetime import datetime, timedelta

# constants
mu_earth = 3.986004415e14
Re = 6378136.3
earth_J2 = 0.001082
omega_E = 7.292115e-5 
omega_s = 0.9856  # deg/day
solar_day_E = 86400

# Mars
mu_mars = 4.2828e13
Rm = 3396.2e3
mars_J2 = 1960.45e-6
solar_day_m = 88775.244

# Problem 1 ------------------------------------------------------------------
def find_zero_crossings(time, beta):
    sign_changes = np.where(np.diff(np.sign(beta)))[0]  # Indices where sign changes
    zero_crossing_times = time[sign_changes]  # Corresponding times
    return zero_crossing_times

def raan_dot(alt, e, i, radius, j2, mu, solar_day):
        
        a = radius + alt
        n = np.sqrt(mu / (a ** 3))

        val = (-3/2) * j2 * ((radius / a) ** 2) * n * np.cos(np.radians(i)) * (180 / np.pi) / (1 - (e ** 2))**2
        val *= solar_day
        T_kep = (2 * math.pi) / n

        #m_dot = n * ((1 - (3/4)*((Re/a)**2)) * j2 * (1 - 3*(math.cos(math.radians(i))**2)) / ((1 - (e ** 2))**(3/2)))
        m_dot = n * (1 - (j2 * (1- 3 * (math.cos(math.radians(i))**2)) * (radius / a)**2 * (3/4) / (1 - (e **2))**(3/2)))
        T_anom = (2 * math.pi) / m_dot

        omega_dot = (-3/4) * n * ((radius / a) ** 2) * j2 * (1 - (5 * math.cos(math.radians(i))**2)) / (1 - (e **2) ** (3/2))
        arg_lat_dot = m_dot + omega_dot
        T_drac = (2 * math.pi) / arg_lat_dot

        gmst_dot = (1 / (86164 / (2 * math.pi))) * (180 / math.pi)
        nodal_day = (360) / abs((val / solar_day) - gmst_dot)

        sun_cyc = ((360) / abs(val - omega_s))

        return val, T_kep, T_anom, T_drac, nodal_day, sun_cyc, math.degrees(omega_dot)

a_vals = np.array([717, 1336, 490, 782, 475, 812]) * 1000
e_vals = np.array([0.0003665, 0.0007665, 0.00005855, 0.0031394, 0.0004512, 0.0205681])
i_vals = np.array([92.0194, 66.0425, 88.9691, 98.7932, 97.4358, 49.8242])
sats = np.array(["Cryosat-2", "TOPEX/Poseidon", "GRACE-FO-1", "ERS-1", "Pelican-2", "Starlette"])

raan_ders = []
sun_cycs = []

for i, val in enumerate(a_vals):
    raan_der, T_kep, T_anom, T_drac, nodal_day, sun_cyc, o_dot = raan_dot(a_vals[i], e_vals[i], i_vals[i], Re, earth_J2, mu_earth, solar_day_E)
    raan_ders.append(raan_der)
    sun_cycs.append(sun_cyc)
    #days = np.linspace(0, sun_cyc * 4)
    days = np.arange(sun_cyc * 4)
    beta_primes = []
    for j, day in enumerate(days):
        lamb, eps, g, capL = find_sun_coor(day)
        #r_sun = 1.00014 - 0.01671*math.cos(math.radians(g)) - 0.00014*math.cos(math.radians(2 * g))
        x = math.cos(math.radians(capL))
        y = math.cos(math.radians(eps)) * math.sin(math.radians(capL))
        z = math.sin(math.radians(eps)) * math.sin(math.radians(capL))
        e_sun = np.array([x, y, z])
        #e_sun = e_sun / np.linalg.norm(e_sun)

        raan_new = (raan_der * day) % 360
        e_h_x = math.sin(math.radians(i_vals[i])) * math.sin(math.radians(raan_new))
        #e_h_x = math.degrees(e_h_x)
        e_h_y = math.sin(math.radians(i_vals[i])) * -1 * math.cos(math.radians(raan_new))
        #e_h_y = math.degrees(e_h_y)
        e_h_z = math.cos(math.radians(i_vals[i]))
        #e_h_z = math.degrees(e_h_z)
        e_h = np.array([e_h_x, e_h_y, e_h_z])

        beta_prime = 90 - math.degrees(math.acos(np.dot(e_h, e_sun)))
        beta_primes.append(beta_prime)
    
    plt.plot(days, beta_primes)
    plt.title(sats[i])
    plt.xlabel("Time (days)")
    plt.ylabel("Beta Prime (degrees)")
    plt.show()

    if sats[i] != "ERS-1" and sats[i] != "Pelican-1":
        zero_crossings = find_zero_crossings(days, beta_primes)
        #zero_crossing_intervals = np.diff(zero_crossings)
        #print(zero_crossings)

        
    #sun_cyc = ((360) / abs(raan_der - omega_s))
    # print("new sat")
    # print("T kep: ", T_kep)
    # print("T anom", T_anom)
    # print("T drac", T_drac)
    # print("Nodal Day", nodal_day)
    # print("Sun cyc", sun_cyc)

# P1d
#print(raan_ders)
# End Problem 1 ------------------------------------------------------------

# Problem 2 ------------------------------------------------------
incs = np.arange(180)
eccens_raan = np.array([0, 0.128, 0.065, 0.037, 0.023, 0.008, 0])
alts_raan = np.array([[2000, 2000], [100, 2000], [100, 1000], [100, 600], [100, 400], [100, 200], [100, 100]]) * 1000

eccens_om = np.array([0, 0.030, 0.065, 0.128, 0.183, 0.231, 0])
alts_om = np.array([[100, 100], [100, 500], [100, 1000], [100, 2000], [100, 3000], [100, 4000], [4000, 4000]]) * 1000

for i, val in enumerate(eccens_raan):
     raan_der_mars = []
     alt_p = alts_raan[i][0] + Rm
     alt_a = alts_raan[i][1] + Rm
     a = (alt_a + alt_p) / 2
     a = a - Rm
     for j, inc in enumerate(incs):
        raan_der_m, T_kep_m, T_anom_m, T_drac_m, nodal_day_m, sun_cyc_mars, o_dot_m = raan_dot(a, val, inc, Rm, mars_J2, mu_mars, solar_day_m)
        raan_der_mars.append(raan_der_m)
     
     label_text = f"e={val}, {int(alts_raan[i][0]/1000)}x{int(alts_raan[i][1]/1000)} km"
     plt.plot(incs, raan_der_mars, label = label_text)

plt.xlabel("Inclination (degrees)")
plt.ylabel("Rate of Change of RAAN (deg/day)")
plt.legend(loc = "upper left", prop={'size': 8})
plt.show()

for i, val in enumerate(eccens_om):
     om_der_mars = []
     alt_p = alts_om[i][0] + Rm
     alt_a = alts_om[i][1] + Rm
     a = (alt_a + alt_p) / 2
     a = a - Rm
     for j, inc in enumerate(incs):
        raan_der_m, T_kep_m, T_anom_m, T_drac_m, nodal_day_m, sun_cyc_mars, o_dot_m = raan_dot(a, val, inc, Rm, mars_J2, mu_mars, solar_day_m)
        om_der_mars.append(o_dot_m)
     
     label_text = f"e={val}, {int(alts_om[i][0]/1000)}x{int(alts_om[i][1]/1000)} km"
     plt.plot(incs, om_der_mars, label = label_text)

plt.xlabel("Inclination (degrees)")
plt.ylabel("Rate of Change of Argument of Perigee (deg/day)")
plt.legend(loc = "upper center", prop={'size': 12})
plt.show()
# End Problem 2 --------------------------------------------------------

# Problem 3 -------------------------------------------------------

i_trmm = np.radians(35)
a_trmm = 350 * 1000
e_trmm = 0

raan_drift_trmm, kep_tr, anom_tr, drac_tr, nodal_tr, suncyc_tr, o_dot_tr = raan_dot(a_trmm, e_trmm, e_trmm, Re, earth_J2, mu_earth, solar_day_E)
#print(raan_drift_trmm)

# need to find local hour angle at t0
jd0 = 2451200.363738
gmst0 = calc_gmst(jd0)
lamb0, eps, g, L = find_sun_coor(jd0 - 2451545.0) # degrees

long0 = -5.157 # degrees

raan0 = gmst0 + long0
#print(raan0)

ra_equin0 = 280.46 + (360/365.2425) * (-245)
#print(ra_equin0)

lha0 = raan0 - ra_equin0
#print(lha0)

year_days = np.arange(0, 365, 1)
lha_vals = []

print(suncyc_tr)

for ind, day in enumerate(year_days):
    days_past = ind - 20
    lha = lha0 + ((360 / suncyc_tr)*days_past)
    lha_vals.append(lha)

lha_vals = np.array(lha_vals)
lmt_vals = [(lha / 15 + 12) % 24 for lha in lha_vals]
# print(lmt_vals)

target_lmt = 22 + (20 / 60)  # 22.3333 hours
tolerance = 0

matching_days = [] 
crossing_days = []  # Store days where LMT crosses 22:20

for i in range(1, len(lmt_vals)): 
    # Direct match within tolerance
    if abs(lmt_vals[i] - target_lmt) <= tolerance:
        matching_days.append(year_days[i])

    # Check if LMT crosses 22:20 (sign change)
    if (lmt_vals[i] - target_lmt) * (lmt_vals[i - 1] - target_lmt) < 0:
        crossing_days.append(year_days[i])

# Convert days to actual dates
matching_dates = [datetime(1999, 1, 1) + timedelta(days=int(d) - 1) for d in matching_days]
crossing_dates = [datetime(1999, 1, 1) + timedelta(days=int(d) - 1) for d in crossing_days]

# Print the results
print("\nDays where LMT = 22:20 ±15 min:")
for date in matching_dates:
    print(date.strftime("%Y-%m-%d"))

print("\nDays where LMT crosses 22:20:")
for date in crossing_dates:
    print(date.strftime("%Y-%m-%d"))

# Plot LMT over time with markers
# plt.figure(figsize=(8, 4))
# plt.plot(year_days, lmt_vals, label="LMT at Ascending Node", linestyle="-")
# #plt.scatter(matching_days, [22.3333] * len(matching_days), color='red', label="LMT = 22:20 Matches", zorder=3)
# plt.scatter(crossing_days, [22.3333] * len(crossing_days), color='blue', label="LMT Crossings", zorder=3)
# plt.axhline(22.3333, color="gray", linestyle="dashed", label="LMT = 22:20 Reference")

# # Labels
# plt.xlabel("Day of Year (1999)")
# plt.ylabel("Local Mean Time (hours)")
# plt.title("LMT at Ascending Node Over 1999")
# plt.legend()
# plt.grid()
# plt.show()
# End Problem 3 ---------------------------------------------------------

# Problem 4 ---------------------------------------------------------------
a_grace = 500000  # 500 km altitude
e_grace = 0  # Circular orbit
i_grace = 89  # Degrees
raan0 = 0  # Assume arbitrary starting RAAN (it will change over time)

mission_years = 12
mission_days = mission_years * 365

raan_drift_gr, kep_gr, anom_gr, drac_gr, nodal_gr, suncyc_gr, o_dot_gr = raan_dot(a_grace, e_grace, i_grace, Re, earth_J2, mu_earth, solar_day_E)

days_gr = np.arange(0, mission_days, 1)
jd0 = 2.4583e6
n = jd0 - 2451545.0
beta_primes_gr = []

for i, day in enumerate(days_gr):
    lamb, eps, g, capL = find_sun_coor(n + day)
    raan_new = (raan0 + (raan_drift_gr * day)) % 360
    #r_sun = 1.00014 - 0.01671*math.cos(math.radians(g)) - 0.00014*math.cos(math.radians(2 * g))
    x = math.cos(math.radians(capL))
    y = math.cos(math.radians(eps)) * math.sin(math.radians(capL))
    z = math.sin(math.radians(eps)) * math.sin(math.radians(capL))
    e_sun = np.array([x, y, z])
        #e_sun = e_sun / np.linalg.norm(e_sun)

    # raan_new = raan0 + raan_der * day
    e_h_x = math.sin(math.radians(i_grace)) * math.sin(math.radians(raan_new))
        #e_h_x = math.degrees(e_h_x)
    e_h_y = math.sin(math.radians(i_grace)) * -1 * math.cos(math.radians(raan_new))
        #e_h_y = math.degrees(e_h_y)
    e_h_z = math.cos(math.radians(i_grace))
        #e_h_z = math.degrees(e_h_z)
    e_h = np.array([e_h_x, e_h_y, e_h_z])

    beta_prime = 90 - math.degrees(math.acos(np.dot(e_h, e_sun)))
    beta_primes_gr.append(beta_prime)

beta_primes_gr = np.array(beta_primes_gr)
plt.plot(days_gr, beta_primes_gr)
plt.xlabel("Days Since May 22, 2018")
plt.ylabel("Beta Prime (degrees)")
plt.title("Beta Prime for GRACE-FO-1")
plt.show()
# End Problem 4 -------------------------------------------------------------------------