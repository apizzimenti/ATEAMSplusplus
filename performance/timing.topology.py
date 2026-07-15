
import numpy as np
import matplotlib.pyplot as plt

EXECS = ["persistence"]
SCALES = [2, 4, 5, 8, 11, 16]
FIELDS = [2, 3]
DIMENSIONS = [2, 4]
ATTEMPTS=100

UNIT = 1

labels = {
	1e0: "1µs",
	1e1: "10µs",
	1e2: "100µs",
	1e3: "1ms",
	1e6: "1s",
	3e7: "30s",
	6e7: "1m",
	9e8: "30m",
	3.6e9: "1h",
	4.32e10: "12h",
	8.64e10: "1d",
	6.048e11: "1w"
}

lw = 2/3

flierprops = dict(
	marker="x",
	markersize=4
)

medianprops = dict(
	linewidth=lw
)

boxprops = dict(
	linewidth=lw
)

whiskerprops = dict(
	linewidth=lw
)

capprops = dict(
	linewidth=lw
)

plotprops = dict(
	widths=0.35,
	whis=(1,99),
	showfliers=False,
	capprops=capprops,
	flierprops=flierprops,
	medianprops=medianprops,
	boxprops=boxprops,
	whiskerprops=whiskerprops
)

cellcounts = {
	2: {
		2: 8,
		4: 32,
		5: 50,
		8: 128,
		11: 242,
		16: 512,
		22: 968,
		32: 2048
	},
	4: {
		2: 96,
		4: 1536,
		5: 3750,
		8: 24576,
		11: 87846,
		16: 393216,
		22: 1405536,
		32: 6291456
	}
}

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})

for EXEC in EXECS:
	for FIELD in FIELDS:
		for DIMENSION in DIMENSIONS:
			boxes = [
				np.loadtxt(f"timing/{EXEC}.{str(SCALE).zfill(2)}.{DIMENSION}.{FIELD}.{ATTEMPTS}.data", dtype=float)/UNIT
				for SCALE in SCALES
			]

			lo, hi = boxes[0].min(), boxes[-1].max()

			fig, ax = plt.subplots(figsize=(4,4))
			boxplot = ax.boxplot(boxes, **plotprops)

			ax.set_yscale("log")
			ax.set_ylabel("$\log($microseconds$)$", fontsize=8)
			ax.set_ylim(lo, hi)

			ax.set_xlabel("cells\n(linear scale)", fontsize=8)
			# ax.set_xticks([cellcounts[DIMENSION][SCALE] for SCALE in SCALES])
			ax.set_xticklabels([f"{cellcounts[DIMENSION][SCALE]}\n({SCALE})" for SCALE in SCALES], fontsize=6)
			ax.set_xlim(0,len(boxes)+1)

			for tick, label in labels.items():
				if lo <= tick <= hi:
					ax.text(0.05, tick, label, fontsize=6, ha="left", alpha=1/2, va="bottom")
					ax.hlines(tick, xmin=0, xmax=len(boxes)+1, color="k", alpha=1/4, zorder=-1000)
			plt.savefig(f"./timing/{EXEC}.{DIMENSION}.{FIELD}.{ATTEMPTS}.jpeg", dpi=300, bbox_inches="tight")

			plt.clf()
			plt.close()
