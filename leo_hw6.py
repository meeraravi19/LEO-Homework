import numpy as np
# from leo_hw4 import raan_dot
import matplotlib.pyplot as plt

# constants
mu_earth = 3.986004415e14
Re = 6378136.3
earth_J2 = 0.001082
omega_E = 7.292115e-5 
omega_s = 0.9856  # deg/day
solar_day_E = 86400

#val, T_kep, T_anom, T_drac, nodal_day, sun_cyc, omega_dot = raan_dot(300, )

uc = np.linspace(0, 360)
lhas = np.zeros_like(uc)
td = 5672.53062
nd = 86131.0075

for i, val in enumerate(uc):
    lha = (td / nd) * val
    lhas[i] = lha

plt.plot(uc, lhas)
plt.xlabel("Mean Argument of Latitude (Uc) (degrees)")
plt.ylabel("Local Solar Hour Angle of Satellite (degrees)")
plt.title("Change in Local Solar Hour Angle vs Mean Argument of Latitude")
plt.show()
