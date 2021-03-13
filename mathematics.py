import math


def to_rad(alpha):
    pi = math.atan(1) * 4
    alpha /= 180
    alpha *= pi
    return alpha


def borders2pi(a):
    pi = math.atan(1) * 4
    while not (0 <= a < (2 * pi)):
        if a >= (2 * pi):
            a -= 2 * pi
        elif a < 0:
            a += 2 * pi
    return a
