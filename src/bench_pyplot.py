#!/usr/bin/env python3
from matplotlib import pyplot as plt

def read_file(File):
    with open(File, 'r') as f:
        data = f.read().split('[')
        intervals = list(map(int, data[1].replace("\n", "").replace("]", "").split(",")))
        values = list(map(float, data[2].replace("\n", "").replace("]", "").split(",")))
        return intervals, values

plt.figure(0)
f_names = ["mult_erl", "mult_dnif", "mult_cnif"]
for f in f_names:
    i,v = read_file("../benchRes/"+ f + ".txt")
    plt.plot(i,v, label=f)
    
plt.title("Running time.")
plt.xlabel("Vector size")
plt.ylabel("Running time [ms]")
plt.legend()
plt.savefig("../benchRes/test_run.png")