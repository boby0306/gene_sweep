
import numpy as np
from scipy.signal import butter, sosfiltfilt, sosfreqz, filtfilt, tf2zpk
import matplotlib.pyplot as plt
from weighting import A_weighting

#########周波数の計算######################
by3OctFcList = (20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000)
octFcList = (31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000)
def calc_fc_true(fc):
    index = by3OctFcList.index(fc)
    fc_ture = 1000*(((10**(3/10))**(1/3))**(index-17))
    return fc_ture

def calc_f_cutoff(fc, bandwidth=1/3):
    fc = calc_fc_true(fc)
    f_lowcut = fc/((10**(3/10))**(bandwidth/2))
    f_highcut = fc*((10**(3/10))**(bandwidth/2))
    return f_lowcut, f_highcut

#########フィルタの作成######################
def mk_passFilt(lowcut, highcut, fs=48000.0, order=8):
    fnyq = 0.5 * fs
    low = lowcut / fnyq
    high = highcut / fnyq
    sos = butter(order, [low, high], btype="bandpass", output="sos")  # 発散するので2次セクション化
    z, p, k = butter(order, [low, high], btype="bandpass", output="zpk")
    IRLen = getIRLen(p)
    return sos, IRLen

def getIRLen(p, eps=1e-9):  # eps:０近似ライン（計算機イプシロン？？,1e-9or7)
    r = np.max(np.abs(p))
    approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
    return approx_impulse_len

#########フィルタ周波数特性の表示######################
def mk_graph_bandpass(sos, fs=48000.0, xlim = None, ylim = [-75, 5]):
    plt.figure()
    freq, h = sosfreqz(sos, worN=24000)
    h[0] = h[1]
    plt.subplot(2, 1, 1)
    db = 20 * np.log10(np.abs(h))
    plt.xscale("log")
    plt.plot(24000 * freq / np.pi, db)
    plt.ylim(ylim[0], ylim[1])
    plt.grid(True)
    plt.ylabel('Gain [dB]')
    plt.title('Frequency Response')
    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])
    plt.xlabel('frequency')
    plt.show()

#app分析用
def get_AEne(data, fs = 48000):
    b, a = A_weighting(fs)
    z, p, k = tf2zpk(b, a)
    IRLen = getIRLen(p)
    filtData = filtfilt(b, a, data, padlen=IRLen)
    Aene = (filtData*filtData + 1e-15).sum()
    return Aene

# low, high = calc_f_cutoff(1000, bandwidth = 1)
# sos, _ = mk_passFilt(low, high)
# low_fig, high_fig = calc_f_cutoff(1000, bandwidth = 1)
# mk_graph_bandpass(sos, xlim=[low_fig, high_fig])

# print(sos)