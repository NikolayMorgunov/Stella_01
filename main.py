import datetime
import time
import sgp4
from skyfield.api import load, Topos, EarthSatellite

TLE_FILE = "https://celestrak.com/NORAD/elements/active.txt"  # DB file to download
SAT_NAME = "NOAA 19"
# LK coords
LATITUDE = 55.93013
LONGITUDE = 37.51832
