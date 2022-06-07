#!/usr/bin/env python3

import matplotlib.pyplot as plt

data = """
perf_1_1_1.txt_v3_0:0.07user 0.00system 1:27.73elapsed 0%CPU (0avgtext+0avgdata 11812maxresident)k
perf_1_1_1.txt_v3_1:0.07user 0.00system 1:27.71elapsed 0%CPU (0avgtext+0avgdata 12048maxresident)k
perf_1_1_1.txt_v3_2:0.07user 0.00system 1:26.97elapsed 0%CPU (0avgtext+0avgdata 12076maxresident)k
perf_12_1_12.txt_v3_0:0.06user 0.01system 0:15.91elapsed 0%CPU (0avgtext+0avgdata 12084maxresident)k
perf_12_1_12.txt_v3_1:0.06user 0.01system 0:16.02elapsed 0%CPU (0avgtext+0avgdata 11908maxresident)k
perf_12_1_12.txt_v3_2:0.07user 0.00system 0:15.67elapsed 0%CPU (0avgtext+0avgdata 12036maxresident)k
perf_16_1_16.txt_v3_0:0.07user 0.00system 0:14.82elapsed 0%CPU (0avgtext+0avgdata 12060maxresident)k
perf_16_1_16.txt_v3_1:0.07user 0.00system 0:13.64elapsed 0%CPU (0avgtext+0avgdata 11900maxresident)k
perf_20_1_20.txt_v3_0:0.07user 0.00system 0:14.00elapsed 0%CPU (0avgtext+0avgdata 11760maxresident)k
perf_20_1_20.txt_v3_1:0.07user 0.00system 0:12.71elapsed 0%CPU (0avgtext+0avgdata 12028maxresident)k
perf_20_1_20.txt_v3_3:0.07user 0.01system 0:13.63elapsed 0%CPU (0avgtext+0avgdata 12044maxresident)k
perf_2_1_2.txt_v3_0:0.07user 0.00system 1:19.82elapsed 0%CPU (0avgtext+0avgdata 12072maxresident)k
perf_2_1_2.txt_v3_1:0.06user 0.01system 1:20.03elapsed 0%CPU (0avgtext+0avgdata 12052maxresident)k
perf_2_1_2.txt_v3_2:0.07user 0.00system 1:20.21elapsed 0%CPU (0avgtext+0avgdata 12072maxresident)k
perf_24_1_24.txt_v3_0:0.07user 0.01system 0:12.36elapsed 0%CPU (0avgtext+0avgdata 12052maxresident)k
perf_24_1_24.txt_v3_1:0.07user 0.01system 0:12.56elapsed 0%CPU (0avgtext+0avgdata 12092maxresident)k
perf_28_1_28.txt_v3_0:0.07user 0.00system 0:11.56elapsed 0%CPU (0avgtext+0avgdata 12048maxresident)k
perf_28_1_28.txt_v3_1:0.07user 0.00system 0:12.21elapsed 0%CPU (0avgtext+0avgdata 12052maxresident)k
perf_32_1_32.txt_v3_0:0.07user 0.00system 0:11.41elapsed 0%CPU (0avgtext+0avgdata 12028maxresident)k
perf_32_1_32.txt_v3_1:0.07user 0.00system 0:10.76elapsed 0%CPU (0avgtext+0avgdata 12056maxresident)k
perf_36_1_36.txt_v3_0:0.07user 0.00system 0:09.66elapsed 0%CPU (0avgtext+0avgdata 12052maxresident)k
perf_36_1_36.txt_v3_1:0.06user 0.01system 0:10.80elapsed 0%CPU (0avgtext+0avgdata 11916maxresident)k
perf_40_1_40.txt_v3_0:0.07user 0.00system 0:09.81elapsed 0%CPU (0avgtext+0avgdata 11880maxresident)k
perf_4_1_4.txt_v3_0:0.07user 0.00system 0:41.24elapsed 0%CPU (0avgtext+0avgdata 11772maxresident)k
perf_4_1_4.txt_v3_1:0.07user 0.03system 0:41.26elapsed 0%CPU (0avgtext+0avgdata 12052maxresident)k
perf_4_1_4.txt_v3_2:0.08user 0.00system 0:41.17elapsed 0%CPU (0avgtext+0avgdata 12060maxresident)k
perf_44_1_44.txt_v3_0:0.06user 0.01system 0:09.30elapsed 0%CPU (0avgtext+0avgdata 11836maxresident)k
perf_48_1_48.txt_v3_0:0.07user 0.00system 0:08.88elapsed 0%CPU (0avgtext+0avgdata 11764maxresident)k
perf_48_1_48.txt_v3_3:0.07user 0.00system 0:08.80elapsed 0%CPU (0avgtext+0avgdata 11828maxresident)k
perf_8_1_8.txt_v3_0:0.07user 0.00system 0:20.81elapsed 0%CPU (0avgtext+0avgdata 12064maxresident)k
perf_8_1_8.txt_v3_1:0.07user 0.00system 0:20.07elapsed 0%CPU (0avgtext+0avgdata 12040maxresident)k
"""


def get_data(l):
    n = int(l.split("_")[1])
    t_str = l.split()[2].split("e")[0]
    t = 100 * (60 * int(t_str[0]) + int(t_str[2:4])) + int(t_str[5:7])
    return n, t


res = {}
for n, t in map(get_data, data.splitlines()[1:]):
    if n not in res:
        res[n] = []
    res[n].append(t)

vals = {n: sum(ts) / len(ts) for n, ts in res.items()}

X, Y = zip(*sorted(vals.items()))

plt.figure(figsize=(8, 6), dpi=80)
plt.xlabel("Number of processes")
plt.ylabel("Speedup")
plt.title("Strong scaling")
plt.plot(X, [Y[0] / y for y in Y])
plt.show()
