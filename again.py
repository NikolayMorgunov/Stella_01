import sgp4
from get_tle import *
import math
import numpy as np
from mathematics import *
import datetime as dt
import matrix_multiplication
import new_time
import matplotlib.pyplot as plt

TLE_FILE = "https://celestrak.com/NORAD/elements/active.txt"  # DB file to download
SAT_NAME = "NOAA 19                 "
# LK coords
LATITUDE = 55.93013
LONGITUDE = 37.51832
pi = math.atan(1) * 4
mu = 2.9755363405824e15
J2 = 1.0826267e-3
Re = 6378.137
Rp = 6356.7523  # полярный радиус Земли
h = 0.2  # высота НКУ над уровнем моря

# получение кеплеровских величин
kep = from_strings(TLE_FILE, SAT_NAME)
# большая полуось
a = (mu / ((2 * pi * kep['freq']) ** 2)) ** (1 / 3)
# перицентр
P = a * (1 - kep['e'])
a1 = a / Re

kep['inclination'] = to_rad(kep['inclination'])
# какие-то параметры
d1 = 3 * J2 * (Re ** 2) * (3 * math.cos(kep['inclination']) ** 2 - 1) / \
     (4 * a1 ** 2 * (1 - kep['e'] ** 2) ** (2 / 3))
a0 = -a1 * (((((134 * d1) ** 3) / 81)) + (d1 ** 2) + (d1 * 3) - 1)
p0 = a0 * (1 - (kep['e']) ** 2)

# ввод интервала
print('Введите начало периода в формате ГГГГ ММ ДД ЧЧ ММ')
inp = input().split()
beg_dt = dt.datetime(int(inp[0]), int(inp[1]), int(inp[2]), int(inp[3]), int(inp[4]))
beg = [int(i) for i in inp]
beg.append(0)
print('Введите конец периода в формате ГГГГ ММ ДД ЧЧ ММ')
inp = input().split()
end_dt = dt.datetime(int(inp[0]), int(inp[1]), int(inp[2]), int(inp[3]), int(inp[4]))
end = [int(i) for i in inp]
end.append(0)

# Расчет элевации и азимута для каждого промежуточного времени
t = beg_dt.timestamp()
cur_time = beg.copy()
all_elev = []
all_azim = []

# настоящее время
now = dt.datetime(kep['year'], 1, 1).timestamp()
now += kep['day'] * 86400

while t <= end_dt.timestamp():
    dt = t - now

    ndt = kep['freq'] * dt
    M = kep['aver anomaly'] + 360 * (ndt - int(ndt) - int((kep['aver anomaly'] + 360 * (ndt - int(ndt)))
                                                          / 360))
    M = borders2pi(to_rad(M))
    E = M / (1 - kep['e']) / 2
    nu = 2 * math.atan2(((1 + kep['e']) / (1 - kep['e'])) ** 0.5 * math.sin(E), math.cos(E))  # в радианах
    nu = borders2pi(nu)

    # геоцентрическое расстояние
    r = P * (1 + kep['e']) / (1 + kep['e'] * math.cos(nu))

    # аргумент перицентра
    omega = kep['periapsis'] + 360 * (
            (3 * J2 * (Re) ** 2 * kep['freq'] * (5 * math.cos(kep['inclination']) ** 2 - 1)) * dt /
            (4 * p0 ** 2))  # в градусах
    omega = borders2pi(to_rad(omega))

    u = omega + nu  # широта для даты расчета
    u = borders2pi(u)
    # разница долгот ...
    d_alpha = math.acos(
        (math.cos(u)) / (1 - math.sin(kep['inclination']) ** 2 * (math.sin(u) ** 2)) ** (1 / 2))  # в радианах

    delta_g = np.sign(math.sin(u)) * math.acos(math.cos(u) / math.cos(d_alpha))  # склонение в радианах

    alpha_OMEGA = kep['node long'] + 360 * (-3) * J2 * Re ** 2 * kep['freq'] * kep['inclination'] * dt \
                  / (2 * p0 ** 2)  # значение долготы восходящего узла в градусах
    alpha_OMEGA = borders2pi(to_rad(alpha_OMEGA))

    alpha_g = d_alpha + alpha_OMEGA

    coords_in_abs = [r * math.cos(alpha_g) * math.cos(delta_g), r * math.sin(alpha_g) * math.cos(delta_g),
                     r * math.sin(delta_g)]  # координаты в АГЭСК

    year = cur_time[0]
    month = cur_time[1]
    day = cur_time[2]
    hour = cur_time[3]
    minutes = cur_time[4]
    seconds = cur_time[5]

    # Расчет местного звездного времени
    a_tilda = int((14 - month) / 12)
    m = month + 12 * a_tilda - 3
    y = year + 4800 - a_tilda
    JDN = day + int((153 * m + 2) / 5) + 365 * y + int(y / 4) - int(y / 100) + int(
        y / 400) - 32045  # номер юлианского дня
    JD = JDN + (hour - 12) / 24 + (minutes / 1440) + (seconds / 86400)  # юлианская дата
    T = (JD - 2451545) / 36525
    T0 = 6.697374558 + 2400.051336 * T + 0.000025862 * (T ** 2)  # число юлианских столетий
    GST = T0 + (hour + (minutes / 60) + (seconds / 3600)) * 1.002737909  # зведное время на меридиане Гринвича

    delta_lambda = LONGITUDE  # разность долгото НКУ и Гринвича
    LST = GST + np.sign(delta_lambda) * 24 * delta_lambda / 360  # местное звездное время
    alpha_i = 2 * pi * LST / 24

    delta_i = LATITUDE / (1 + 2 * pi / 180)  # склонение НКУ
    delta_i = to_rad(delta_i)
    r_g = (((math.cos(delta_i) / Re) ** 2) + (math.cos(delta_i) / Rp) ** 2) ** (-0.5) + h

    coords_in_NKY = [r_g * math.cos(alpha_i) * math.cos(delta_i),
                     r_g * math.sin(alpha_i) * math.cos(delta_i),
                     r_g * math.sin(delta_i)]  # координаты НКУ

    # ФИНАЛЬНЫЙ ПУНКТ

    # вектор дальности КА по отношению к НКУ
    delta_r = [(coords_in_abs[i] - coords_in_NKY[i]) for i in range(3)]

    x_s = delta_r[0]
    y_s = delta_r[1]
    z_s = delta_r[2]

    delta_r_pro = (x_s ** 2 + y_s ** 2 + z_s ** 2) ** 0.5  # проекци вектора delta_r
    cosinusy = [[-math.cos(alpha_i) * math.sin(delta_i), -math.sin(delta_i) * math.sin(alpha_i), math.cos(delta_i)],
                [math.cos(delta_i) * math.cos(alpha_i), math.cos(delta_i) * math.sin(alpha_i), math.sin(delta_i)],
                [-math.sin(alpha_i), math.cos(alpha_i), 0]]
    coords_div_dr = [x_s / delta_r_pro, y_s / delta_r_pro, z_s / delta_r_pro]

    coords_st = matrix_multiplication.multiplic(cosinusy, coords_div_dr)

    x_st = coords_st[0]
    y_st = coords_st[1]
    z_st = coords_st[2]

    elevation = math.acos(((x_st ** 2) + (z_st ** 2)) ** (1 / 2) / delta_r_pro)  # элевация
    elevation = 180 * elevation / pi

    A = math.atan(z_st / x_st)  # азимут
    A = 180 * A / pi
    if A < 0:
        A += 360
    print('theta =', elevation)
    print('phi =', A)
    print()
    t += 60
    cur_time = new_time.new_time(cur_time)

    # Реализация полярного графика
    all_elev.append(elevation)
    all_azim.append(A)
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim(bottom=90, top=0)
ax.plot(all_azim, all_elev)
fig.set_size_inches(7, 7)
plt.show()
