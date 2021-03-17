def pars(lst):
    ans = {}
    ans['year'] = 2000 + int(lst[1][18:20]) #текущий год
    ans['day'] = float(lst[1][20:32]) #день в году
    ans['inclination'] = float(lst[2][8:16]) #наклонение в градусах
    ans['node long'] = float(lst[2][17:25]) #долгота восходящего узла
    ans['e'] = float('.' + lst[2][26:33]) #эксцентриситет
    ans['aver anomaly'] = float(lst[2][43:51]) #средняя аномалия
    ans['periapsis'] = float(lst[2][43:51]) #аргумент перицентра
    ans['freq'] = float(lst[2][52:63]) #частота обращения витков/день
    return ans
