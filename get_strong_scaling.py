#!/usr/bin/env python3

import matplotlib.pyplot as plt

data = """
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
