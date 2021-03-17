import datetime
import time
import sgp4
from get_tle import *
from skyfield.api import load, Topos, EarthSatellite
import math
import numpy as np
import mathematics as m

TLE_FILE = "https://celestrak.com/NORAD/elements/active.txt"  # DB file to download
SAT_NAME = "NOAA 19                 "
# LK coords
LATITUDE = 55.93013
LONGITUDE = 37.51832

kep = from_strings(TLE_FILE, SAT_NAME)
print(kep)
pi = math.atan(1) * 4
mu = 2.9755363405824e15

a = (mu / ((2 * pi * kep['freq']) ** 2)) ** (1 / 3)  # большая полуось
P = a * (1 - kep['e'])  # перицентр

dt = 6.944444e-4
ndt = int(dt * kep['aver anomaly'])
big_brackets = int((kep['aver anomaly'] + 360 * (dt * kep['aver anomaly'] - ndt)) / 360)

M = kep['aver anomaly'] + big_brackets
M = m.to_rad(M)  # перевод в радианы
E = M
nu = 2 * math.atan(((1 + kep['e']) / (1 - kep['e'])) ** 0.5 * math.tan(E / 2))
r = P * (1 + kep['e']) / (1 + kep['e'] * math.cos(nu))

J2 = 1.0826267e-3
Re = 6378.137
a1 = a / Re
kep['inclination'] = m.to_rad(kep['inclination'])
d1 = 3 * J2 * (Re ** 2) * (3 * math.cos(kep['inclination']) ** 2 - 1) / \
     (4 * a1 ** 2 * (1 - kep['e'] ** 2) ** (2 / 3))
a0 = -a1 * (((((134 * d1) ** 3) / 81)) + (d1 ** 2) + (d1 * 3) - 1)
p0 = a0 * (1 - (kep['e']) ** 2)
omega = kep['periapsis'] + 360 * (
        ((3 * J2 * (Re) ** 2 * kep['freq'] * (5 * math.cos(kep['inclination']) ** 2) - 1)) * dt /
        (4 * p0 ** 2))

omega = m.to_rad(omega)
u = omega + nu  # широта для даты расчета

d_alpha = math.acos(
    (math.cos(u)) / (1 - math.sin(kep['inclination']) ** 2 * (math.sin(u) ** 2)) ** (1 / 2))  # разница долгот ...

delta_g = np.sign(math.sin(u)) * math.acos(math.cos(u) / math.cos(d_alpha))  # склонение

alpha_OMEGA = kep['node long'] + 360 * (-3) * J2 * Re ** 2 * kep['freq'] * kep['inclination'] * dt \
              / (2 * p0 ** 2)  # значение долготы восходящего узла
alpha_OMEGA = m.to_rad(alpha_OMEGA)

alpha_g = d_alpha + alpha_OMEGA

coords_in_abs = [r * math.cos(alpha_g) * math.cos(delta_g), r * math.sin(alpha_g) * math.cos(delta_g),
                 r * math.sin(delta_g)]  # координаты в АГЭСК
