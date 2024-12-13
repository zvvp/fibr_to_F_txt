import numpy as np
from scipy.signal import medfilt, filtfilt
from functions import get_S, get_fname, get_start_time, parse_B_txt, get_coef_fibr, get_ranges_fibr, get_time_qrs, \
    get_diff_time, del_V_S, moving_average, get_number_of_peaks


def get_P(lead1, lead2, lead3, intervals, r_pos, chars):
    bl = [0.01591456, 0.03182911, 0.01591456]
    al = [1.0, -1.61277905,  0.67643727]
    bh = [0.9989957, - 0.9989957]
    ah = [1.0, - 0.9979914]
    mean_interval = np.mean(intervals) * 2.5
    out = np.zeros(len(r_pos))
    for i in range(1, len(r_pos) - 1):
        if chars[i] == 'N':
            len_pr = int(round(intervals[i] * 0.36))
            start = r_pos[i] - len_pr
            stop = r_pos[i] - 11
            fragment_l1 = lead1[start:stop]
            fragment_l1 = filtfilt(bl, al, fragment_l1)
            fragment_l1 = filtfilt(bh, ah, fragment_l1)
            fragment_l2 = lead2[start:stop]
            fragment_l2 = filtfilt(bl, al, fragment_l2)
            fragment_l2 = filtfilt(bh, ah, fragment_l2)
            fragment_l3 = lead3[start:stop]
            fragment_l3 = filtfilt(bl, al, fragment_l3)
            fragment_l3 = filtfilt(bh, ah, fragment_l3)
            n1 = get_number_of_peaks(fragment_l1)
            n2 = get_number_of_peaks(fragment_l2)
            n3 = get_number_of_peaks(fragment_l3)
            sum_n = n1 + n2 + n3

            if sum_n == 3:
                out[i] = 0.0
            elif sum_n == 2:
                out[i] = 50.0
            elif sum_n == 1:
                out[i] = mean_interval
            elif sum_n == 0:
                out[i] = mean_interval
        else:
            out[i] = 50.0
    out = moving_average(out, 111)
    return out
# def get_P(lead1, lead2, lead3, intervals, r_pos, chars):
#     bl = [0.09635722, 0.19271443, 0.09635722]
#     al = [1.0, -0.95071164, 0.33614051]
#     # bl = [0.01591456, 0.03182911, 0.01591456]
#     # al = [1.0, -1.61277905,  0.67643727]
#     bh = [0.9989957, - 0.9989957]
#     ah = [1.0, - 0.9979914]
#     out = np.ones(len(r_pos))
#     for i in range(1, len(r_pos) - 1):
#         if chars[i] == 'N':
#             len_pr = int(round(intervals[i] * 0.36))
#             start = r_pos[i] - len_pr
#             stop = r_pos[i] - 11
#             fragment_l1 = lead1[start:stop]
#             fragment_l1 = filtfilt(bl, al, fragment_l1)
#             fragment_l1 = filtfilt(bh, ah, fragment_l1)
#             fragment_l2 = lead2[start:stop]
#             fragment_l2 = filtfilt(bl, al, fragment_l2)
#             fragment_l2 = filtfilt(bh, ah, fragment_l2)
#             fragment_l3 = lead3[start:stop]
#             fragment_l3 = filtfilt(bl, al, fragment_l3)
#             fragment_l3 = filtfilt(bh, ah, fragment_l3)
#             n1 = get_number_of_peaks(fragment_l1)
#             n2 = get_number_of_peaks(fragment_l2)
#             n3 = get_number_of_peaks(fragment_l3)
#             sum_n = n1 + n2 + n3
#             if sum_n == 3:
#                 out[i] = 0.0  #  out[i] = 0.3
#             elif sum_n == 2:
#                 out[i] = 0.0  # 0.7
#             elif sum_n == 1:
#                 out[i] = 1000.0  # 1.7
#             elif sum_n == 0:
#                 out[i - 1:i + 2] = 1000.0  # 2.0
#         else:
#             out[i] = 0.0   # 0.3
#     win = int(np.mean(intervals) * 2)
#     out = moving_average(out, win)
#     return out

# get_S()
fname = get_fname()
if fname == 'Unknown':
    text = "No file *.ecg.\n"
else:
    start_time = get_start_time(fname)
    text = "\n"
    text += f"{fname}\n\n"
    text += f"start_time {start_time[0]}:{start_time[1]}:{start_time[2]}\n\n"

    lead1 = np.load("d:/Kp_01/clean_lead1.npy")
    lead2 = np.load("d:/Kp_01/clean_lead2.npy")
    lead3 = np.load("d:/Kp_01/clean_lead3.npy")
    # lead1 = np.load("clean_lead1.npy")
    # lead2 = np.load("clean_lead2.npy")
    # lead3 = np.load("clean_lead3.npy")
    get_S()
    r_pos, intervals, chars, forms = parse_B_txt()
    pzub = get_P(lead1, lead2, lead3, intervals, r_pos, chars)
    mean_pzub = np.mean(pzub)
    pzub = (pzub - mean_pzub) * 0.8 + mean_pzub
    fintervals = del_V_S(intervals, chars)
    coef_fibr = get_coef_fibr(fintervals)
    n = 51
    fintervals = medfilt(fintervals, n)
    fintervals = moving_average(fintervals, n)

    coef_fibr = medfilt(coef_fibr, n)
    coef_fibr = moving_average(coef_fibr, n)
    mean_fibr = np.mean(coef_fibr)
    coef_fibr = (coef_fibr - mean_fibr) * 0.8 + mean_fibr
    p_coef_fibr = coef_fibr * pzub / 200
    mean_p_fibr = np.mean(p_coef_fibr)
    p_coef_fibr = (p_coef_fibr - mean_p_fibr) * 0.9 + mean_p_fibr


    ranges_fibr = get_ranges_fibr(fintervals, p_coef_fibr, r_pos)

    len_f = len(ranges_fibr[0])
    sum_time = np.array([0, 0, 0])
    for i in range(len_f):
        start = get_time_qrs(ranges_fibr[0][i], start_time)
        stop = get_time_qrs(ranges_fibr[1][i], start_time)
        diff = np.array(get_diff_time(start, stop))
        sum_time += diff
        time_start = np.array(start)
        time_stop = np.array(stop)
        if time_start[0] >= 24:
            time_start[0] = time_start[0] - 24
        if time_stop[0] >= 24:
            time_stop[0] = time_stop[0] - 24
        text += f"{time_start[0]}:{time_start[1]}:{time_start[2]}  {time_stop[0]}:{time_stop[1]}:{time_stop[2]}\n"
    if sum_time[2] > 59:
        overflow_s = sum_time[2] // 60
        sum_time[1] += overflow_s
        sum_time[2] %= 60
    if sum_time[1] > 59:
        overflow_m = sum_time[1] // 60
        sum_time[0] += overflow_m
        sum_time[1] %= 60
    sum_fibr_time = f"{sum_time[0]}:{sum_time[1]}:{sum_time[2]}\n\n"
    text += f"\nВсего эпизодов {len_f}, суммарное время фибрилляции {sum_fibr_time}"
with open("C:/EcgVar/F.txt", "w") as f:
    for i, line in enumerate(text):
        f.write(line)