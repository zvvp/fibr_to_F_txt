import os
import numpy as np
from numba import njit


def get_periods(lines):
    period_2 = int(lines[0].split(';')[1])
    period_1 = int(lines[1].split(';')[1])
    period = int(lines[2].split(';')[1])
    period1 = int(lines[3].split(';')[1])
    period2 = int(lines[4].split(';')[1])
    return period_2, period_1, period, period1, period2

@njit
def get_coef_cor(x: np.ndarray, y: np.ndarray) -> float:
    mean_x: float = np.mean(x)
    mean_y: float = np.mean(y)
    mean_xy: float = np.mean(x * y)
    std_x: float = np.std(x)
    std_y: float = np.std(y)
    if std_x * std_y:
        return (mean_xy - mean_x * mean_y) / (std_x * std_y)
    else:
        return 0
            
def get_S():
    try:
        os.remove("C:/EcgVar/B1.txt")
    except FileNotFoundError:
        pass
    try:
        os.remove("C:/EcgVar/F.txt")
    except FileNotFoundError:
        pass
    s = 0
    with open("C:/EcgVar/B.txt", "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if (i >= 14) and (i < len(lines) - 2):  # i > 13   i < len(lines) - 1
            if (';N' in lines[i]) and (not ';V' in lines[i - 1]) and (not ';S' in lines[i - 1]):
                periods = get_periods(lines[i - 2:i + 3])
                tf = np.array(periods)
                ref_t = np.array([tf[1], tf[1] * 0.8, tf[1] * 1.2, tf[1]])
                ref_t0 = np.array([tf[1], tf[1] * 0.65, tf[1] * 1.15, tf[1]])
                ref_t1 = np.array([tf[1], tf[1] * 0.65, tf[1] * 1.1, tf[1] * 0.65])
                ref_t3 = np.array([tf[1], tf[1] * 1.6, tf[1], tf[1] * 1.6])
                ref_tA = np.array([tf[1], tf[1] * 0.5, tf[1] * 0.5, tf[1]])
                ref_tA1 = np.array([tf[1], tf[1] * 0.65, tf[1] * 0.35, tf[1]])
                ref_tA2 = np.array([tf[1], tf[1] * 0.35, tf[1] * 0.65, tf[1]])
                max_t = np.max(tf)
                min_t = np.min(tf)
                mean_t = (np.sum(tf) - np.max(tf) - np.min(tf)) / 3
                if (tf[2] > 50) and (tf[3] > 50) and (tf[2] < mean_t * 2) and (tf[3] < mean_t * 2):
                    coef_cor = get_coef_cor(ref_t, tf[1:])
                    coef_cor0 = get_coef_cor(ref_t0, tf[1:])
                    coef_cor1 = get_coef_cor(ref_t1, tf[1:])
                    coef_cor3 = get_coef_cor(ref_t3, tf[1:])
                    coef_corA = get_coef_cor(ref_tA, tf[1:])
                    coef_corA1 = get_coef_cor(ref_tA1, tf[1:])
                    coef_corA2 = get_coef_cor(ref_tA2, tf[1:])
                    trs = 0.975  # trs = 0.985
                    if (coef_corA > 0.95) or (coef_corA1 > 0.95) or (coef_corA2 > 0.95):
                        if (tf[3] + tf[2]) < tf[1] * 1.1:
                            lines[i] = lines[i].replace(';N', ';A')
                            lines[i + 1] = lines[i + 1].replace(';N', ';A')
                    elif (coef_cor > trs) or (coef_cor0 > trs) or (coef_cor1 > trs):
                        if (tf[3] - tf[2]) > tf[1] * 0.06:
                            lines[i] = lines[i].replace(';N', ';S')
                            s += 1
                    elif coef_cor3 > trs:
                        if (tf[3] - tf[2]) > tf[1] * 0.06:
                            lines[i + 1] = lines[i + 1].replace(';N', ';S')
                            s += 1
                else:
                    lines[i] = lines[i].replace(';N', ';A')
                    lines[i + 1] = lines[i + 1].replace(';N', ';A')
    lines[6] = lines[6] + f"НЖ: {s}"
    with open("C:/EcgVar/B1.txt", "w") as f:
        for i, line in enumerate(lines):
            f.write(line)
            
def get_fname():
    dir = os.getcwd()
    for file in os.listdir(dir):
        if file.endswith(".ecg"):
            fname = os.path.join(dir, file)
            return fname
    return "Unknown"

def get_start_time(fname):
    with open(fname, "rb") as f:
        f.seek(151)
        dlmt = f.read(1)
        if dlmt == b":":
            f.seek(150)
            start_h = int(f.read(1))
            f.seek(152)
            start_m = int(f.read(2))
            f.seek(155)
            start_s = int(f.read(2))
        else:
            f.seek(150)
            start_h = int(f.read(2))
            f.seek(153)
            start_m = int(f.read(2))
            f.seek(156)
            start_s = int(f.read(2))
    return start_h, start_m, start_s

def parse_B_txt():
    r_pos = []
    intervals = []
    chars = []
    forms = []
    with open("C:/EcgVar/B1.txt", "r") as f:
        for line in f:
            if ';' in line:
                temp = int(line.split(';')[0])
                r_pos.append(temp)
                temp = int(line.split(';')[1])
                intervals.append(temp)
                temp = line.split(';')[2][0]
                chars.append(temp)
                temp = int(line.split(':')[1])
                forms.append(temp)
        r_pos = np.array(r_pos)
        intervals = np.array(intervals)
        chars = np.array(chars)
        forms = np.array(forms)
        # intervals[intervals < 50] = 50
        # intervals[intervals > 500] = 500

    return r_pos, intervals, chars, forms

def del_V_S(intervals, chars):
    len_in = len(intervals)
    out = intervals.copy()
    for i in np.arange(3, len_in - 3):
        if ('V' in chars[i]) or ('S' in chars[i]) or ('A' in chars[i]):
            # out[i:i+2] = (np.sum(intervals[i-2:i+3]) - np.sum(intervals[i:i+2])) / 3
            out[i:i + 2] = np.mean(intervals[i-3:i+4])
    return out

def get_coef_fibr(intervals, chars):
    len_in = len(intervals)
    out = np.zeros(len_in)
    for i in np.arange(2, len_in - 3):
        win_t = intervals[i - 2:i + 3].copy()
        # win_t[win_t > np.mean(win_t) * 2] = np.mean(win_t)
        diff_t = np.abs(win_t - np.roll(win_t, 1))
        diff_t = np.sort(diff_t)
        sum_diff_tf = np.sum(diff_t[:3])
        win_t = np.sort(win_t)
        mean_win_t = np.mean(win_t[1:-1])
        out[i] = (sum_diff_tf * (1 + 100000 / mean_win_t ** 2)) ** 2 * 0.0019 * 5
    return out

def get_ranges_fibr(fintervals, fcoef_fibr, r_pos):
    start = []
    stop = []
    temp = 0
    w = 1500   # 7500
    for i in range(1, len(fintervals)):
        if (fcoef_fibr[i-1] <= fintervals[i-1]) and (fcoef_fibr[i] > fintervals[i]):
            if len(start) == 0:
                start.append(r_pos[i])
                temp = r_pos[i]
                continue
            if (r_pos[i] - temp) >= w:
                start.append(r_pos[i])
                temp = r_pos[i]
            elif len(stop) > 0:
                stop.pop(-1)
                if len(stop) > 0:
                    temp = stop[-1]
        elif (fcoef_fibr[i-1] >= fintervals[i-1]) and (fcoef_fibr[i] < fintervals[i]):
            if len(stop) == 0:
                stop.append(r_pos[i])
                temp = r_pos[i]
                continue
            if (r_pos[i] - temp) >= w:
                stop.append(r_pos[i])
                temp = r_pos[i]
            elif len(start) > 0:
                start.pop(-1)
                if len(start) > 0:
                    temp = start[-1]

    return start, stop

def get_time_qrs(addr, start_time):
    s = addr * 4 // 1000
    m = s // 60
    s = s % 60
    s = s + start_time[2]
    if s >= 60:
        s = s - 60
        m = m + 1
    h = m // 60
    m = m % 60
    m = m + start_time[1]
    if m >= 60:
        m = m - 60
        h = h + 1
    d = h // 24
    h = h % 24
    h = h + start_time[0]
    if h >= 24:
        h = h - 24
        d = d + 1
    # d = d + 1
    h = 24 * d + h
    # return f"{d:02d} день {h:02d}:{m:02d}:{s:02d}"
    return h, m, s

def get_diff_time(start, stop):
    h1, m1, s1 = start
    h2, m2, s2 = stop
    if s1 <= s2:
        diff_s = s2 - s1
    else:
        diff_s = s2 + 60 - s1
        m2 -= 1
    if m1 <= m2:
        diff_m = m2 - m1
    else:
        diff_m = m2 + 60 - m1
        h2 -= 1
    if h1 <= h2:
        diff_h = h2 - h1
    else:
        diff_h = h2 + 24 - h1

    return diff_h, diff_m, diff_s