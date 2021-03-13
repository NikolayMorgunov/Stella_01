def new_time(time):
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    time[4] += 1
    if time[4] == 60:
        time[4] = 0
        time[3] += 1
    if time[3] == 24:
        time[3] = 0
        time[2] += 1
    if time[2] > days_in_months[time[1] - 1]:
        time[2] = 1
        time[1] += 1
    if time[1] == 13:
        time[1] = 1
        time[0] += 1
    return time
