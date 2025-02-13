import math
import numpy as np
from leo_hw1 import calc_orbital_elems, calc_pos_vel

h = 23 + (59 / 60) + (42/3600)
# print(h)

# jd0 = int(B / 4) + int(C / 12) - int(3 * E / 4) + 31 - 32075
# jd0 = (367 * 2023) - ((7 * (2023+(16/12))) / 4) + (275 * 7 / 9) + 31 + 1721013.5
jd0 = (367 * 2023 
          - math.floor((7 * (2023 + math.floor((7 + 9) / 12))) / 4) 
          + math.floor((275 * 7) / 9) 
          + 31 + 1721013.5 + (h / 24)
          - 0.5 * math.copysign(1, (100 * 2023 + 7 - 190002.5)) - 0.5)

#print(jd0)

def calc_gmst(jd0):


    jdut = jd0 + (23.995000000000000 / 24)
    # print(jdut)

    dtt = jdut - 2451545
    # print(dtt)

    t = dtt / 36525
    # print(t)

    dut = jd0 - 2451545
    # print(dut)

    gmst = (6.697375 + 0.065707485828 * dut + 1.0027379 * h + 0.0854103 * t + 0.0000258 * (t * t)) % 24
    #print(gmst)

    gmst = gmst * (360 / 24)
    #print(gmst)
    
    return gmst

# Problem 1, Part c ---------------------------------
# hnow = 16 + (27 / 60) + (38/3600)

# jd0_now = (367 * 2025 
          # - math.floor((7 * (2025 + math.floor((2 + 9) / 12))) / 4) 
          # + math.floor((275 * 2) / 9) 
          # + 31 + 1721013.5 + (hnow / 24)
          # - 0.5 * math.copysign(1, (100 * 2025 + 2 - 190002.5)) - 0.5)

# print(jd0_now)
# gmst_now = calc_gmst(jd0_now)
# print(gmst_now)
# ---------------------------------------------------

r_eci = np.array([6093861.926471250131726, -1840393.570987288141623, 2590090.032104292418808])
v_eci = np.array([-2709.870119294875167, 975.605814756081486, 7047.604436212530345])

gmst = calc_gmst(jd0)
gmst = math.radians(gmst)

r3_gmst = [[math.cos(gmst), math.sin(gmst), 0],
           [-math.sin(gmst), math.cos(gmst), 0],
           [0, 0, 1]]

r_ecef_est = np.dot(r3_gmst, r_eci)
long_est = math.atan2(r_ecef_est[1], r_ecef_est[0])
# print(math.degrees(long_est))
lat_est = math.atan2(r_ecef_est[2], math.sqrt((r_ecef_est[0] * r_ecef_est[0]) + (r_ecef_est[1] * r_ecef_est[1])))
# print(math.degrees(lat_est))

r_ecef = np.array([5259791.807165756821632, 3575614.384020817000419, 2603908.595081204082817])
long_real = math.atan2(r_ecef[1], r_ecef[0])
#print(long_real)
lat_real = math.atan2(r_ecef[2], math.sqrt(r_ecef[0] * r_ecef[0] + (r_ecef[1] * r_ecef[1])))

diff_long = math.degrees(abs(long_est - long_real))
diff_lat = math.degrees(abs(lat_est - lat_real))
# print(diff_long)
# print(diff_lat)

a, e, i, raan, s_omega, m_anom = calc_orbital_elems(r_eci / 1000, v_eci / 1000)
r, v, theta = calc_pos_vel(a, e, i, raan, s_omega, m_anom)

arg_lat = math.radians(s_omega) + math.radians(theta)
sin_phi = math.sin(math.radians(i)) * math.sin(arg_lat)
latitude = math.asin(sin_phi)

cos_lambda = math.cos(arg_lat) / math.cos(latitude)
sin_lambda = (math.sin(arg_lat) / math.cos(latitude)) * math.cos(math.radians(i))
lambda_val = math.atan2(sin_lambda, cos_lambda)

longitude = lambda_val + math.radians(raan) - gmst
#print(math.degrees(longitude))
#print(math.degrees(latitude))
# longitude = (longitude + 180) % 360 - 180

diff_longitude = math.degrees(abs(longitude - long_real))
diff_latitude = math.degrees(abs(latitude - lat_real))
# print(diff_longitude)
#print(diff_latitude)

w = [0, 0, 7.29*(10**-5)]
v_ecef_est = v_eci - np.cross(w, r_ecef_est)
# print(v_ecef_est)