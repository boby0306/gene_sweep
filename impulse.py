import numpy as np
from math import log10

def calc_RT(data, decaydB = 30, start_dB = 5.0):
    times = 60 / decaydB
    # decay curve
    ir_bpf_square = data ** 2.0
    ir_bpf_square_sum = np.sum(ir_bpf_square)
    temp = ir_bpf_square_sum
    curve=[]
    for i in range(len(data)):
        temp = temp - ir_bpf_square[i]
        curve.append(temp)
    curve_dB = 10.0 * np.log10(curve)
    curve_offset = max(curve_dB)
    decay_curve = curve_dB - curve_offset

    # find regression target
    i = 0
    while decay_curve[i] > - start_dB:
        i += 1
    start_sample = i
    print(f"start:{i}")
    while decay_curve[i] > -1.0*decaydB - start_dB:
        i += 1
    end_sample = i
    print(f"end:{i}")
    regression_target = decay_curve[start_sample:end_sample]

    # linear regression for T
    x = np.linspace(start_sample, end_sample, end_sample-start_sample)
    a, b = np.polyfit(x, regression_target, 1)
    rt_sec = (-1.0*decaydB/a)*times/48000
    return rt_sec, decay_curve