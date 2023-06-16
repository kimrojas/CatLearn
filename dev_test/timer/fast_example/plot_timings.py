import numpy as np
import matplotlib.pyplot as plt
import sys

# data = np.genfromtxt("info.txt", delimiter=",")

datadict = {
    "Train": {"mean": [], "std": []},
    "ML NEB Opt": {"mean": [], "std": []},
    "Eval": {"mean": [], "std": []},
}

threads = list(range(1, 9))
for thread in threads:
    file = f"info-{thread}.txt"
    data = np.genfromtxt(file, delimiter=",")

    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    datadict["Train"]["mean"].append(mean[0])
    datadict["Train"]["std"].append(std[0])
    datadict["ML NEB Opt"]["mean"].append(mean[1])
    datadict["ML NEB Opt"]["std"].append(std[1])
    datadict["Eval"]["mean"].append(mean[2])
    datadict["Eval"]["std"].append(std[2])

fig, ax = plt.subplots(1, 3, figsize=(10, 5))
ax[0].errorbar(
    threads,
    datadict["Train"]["mean"],
    yerr=datadict["Train"]["std"],
    fmt="-o",
    capsize=5,
    label="Train",
)
ax[1].errorbar(
    threads,
    datadict["ML NEB Opt"]["mean"],
    yerr=datadict["ML NEB Opt"]["std"],
    fmt="-o",
    capsize=5,
    label="ML NEB Opt",
)
ax[2].errorbar(
    threads,
    datadict["Eval"]["mean"],
    yerr=datadict["Eval"]["std"],
    fmt="-o",
    capsize=5,
    label="Eval",
)

for i, k in enumerate(datadict.keys()):
    ax[i].set_xlabel("Threads")
    ax[i].set_ylabel("Time (s)")
    ax[i].set_title(k)

fig.tight_layout()
plt.show()
fig.savefig("timings-threads.png")


sys.exit()
