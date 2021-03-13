import requests
import pars_tle

def from_strings(tle_file, sat_name):
    r = requests.get(tle_file, allow_redirects=True)
    open('tle.txt', 'wb').write(r.content)
    txt = open('tle.txt', 'r').read().split('\n')
    for i in range(len(txt)):
        if (txt[i] == sat_name):
            print(sat_name, txt[i + 1], txt[i + 2], sep='\n')
            tle = pars_tle.pars([sat_name, txt[i + 1], txt[i + 2]])
            return tle


