#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt

STATE_LABELS = [1, 2, 3, 15, 20, 30]

p = argparse.ArgumentParser()
p.add_argument("--n", type=int, required=True)
p.add_argument("--spectrum", default="../data/spectrum.csv")
p.add_argument("--vectors", default="../data/eigenvectors.csv")
args = p.parse_args()


spec = np.loadtxt(args.spectrum)
k, e_num, e_exact = spec[:, 0], spec[:, 1], spec[:, 2]

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(k, e_num, "k.", markersize=4, label="numeric")
ax.plot(k, e_exact, "r.", markersize=4, label="exact")
ax.set_xlabel("k")
ax.set_ylabel(r"$E_k$")
ax.set_title(f"number of nodes: n = {args.n}")
ax.legend()
fig.tight_layout()
fig.savefig(f"spectrum_n{args.n}.png", dpi=150)
plt.close(fig)

data = np.loadtxt(args.vectors)
x = data[:, 0]

fig, axes = plt.subplots(3, 2, figsize=(10, 9))
for s, (label, ax) in enumerate(zip(STATE_LABELS, axes.flat)):
    psi_num = data[:, 1 + 2*s]
    psi_exact = data[:, 2 + 2*s]
    ax.plot(x, psi_exact, "k-", label=rf"$\psi_{{{label},exact}}$")
    ax.plot(x, psi_num, "r.", markersize=4, label=rf"$\psi_{{{label},num}}$")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\psi$")
    ax.legend(loc="upper right")
fig.suptitle(f"n = {args.n}")
fig.tight_layout()
fig.savefig(f"states_n{args.n}.png", dpi=150)
plt.close(fig)