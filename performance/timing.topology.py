
import numpy as np
import matplotlib.pyplot as plt

HOST = "pangolin"
EXECS = ["persistence"]
FIELDS = [2, 3]
DIMENSIONS = [2]

MODES = ["SERIAL", "PARALLEL"]
COLORS = ["#DDAA33", "#BB5566"]

OFFSET = 0.2
POSITIONS = [-OFFSET, OFFSET]

SCALES = {
	2: [2, 4, 5, 8, 11, 16, 22, 32, 45, 64, 90],
	4: [2, 4, 5, 8]
}

UNIT = 1

WIDTHS = [512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]

cellcounts = {
	2: { 2: 16, 4: 64, 5: 100, 8: 256, 11: 484, 16: 1024, 22: 1936, 32: 4096, 45: 8100, 64: 16384, 90: 32400, 128: 65536, 181: 131044, 256: 262144, 362: 524176, 512: 1048576,  },
	4: { 2: 256, 4: 4096, 5: 10000, 8: 65536, 11: 234256, 16: 1048576, 22: 3748096, 32: 16777216,  }
}

density = {
	2: { 2: 0.222222, 4: 0.0555556, 5: 0.0355556, 8: 0.0138889, 11: 0.00734619, 16: 0.00347222, 22: 0.00183655, 32: 0.000868056, 45: 0.000438957, 64: 0.000217014, 90: 0.000109739, 128: 5.42535e-05, 181: 2.71325e-05, 256: 1.35634e-05, 362: 6.78313e-06, 512: 3.39084e-06,  },
	4: { 2: 0.0177778, 4: 0.00111111, 5: 0.000455111, 8: 6.94444e-05, 11: 1.94279e-05, 16: 4.34028e-06, 22: 1.21425e-06, 32: 2.71267e-07,  }
}

labels = {
	1e0: "1$\mu$s",
	1e1: "10$\mu$s",
	1e2: "100$\mu$s",
	1e3: "1ms",
	# 0.5e6: r"$\nicefrac 12$s",
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
	linewidth=0,
	color="none"
)

boxprops = dict(
	linewidth=0
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
	whiskerprops=whiskerprops
)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman",
	"text.latex.preamble": r"\usepackage{nicefrac}"
})

for EXEC in EXECS:
	for FIELD in FIELDS:
		for DIMENSION in DIMENSIONS:
			fig, ax = plt.subplots(figsize=(4,4))

			PREFIX = f"{HOST}.{EXEC}.{FIELD}.{DIMENSION}"

			lo = np.inf
			hi = 0
			
			for POSITION, COLOR, MODE in zip(POSITIONS, COLORS, MODES):
				boxes = [
					np.loadtxt(f"timing/{PREFIX}.{str(SCALE).zfill(3)}.{MODE}.data")
					for SCALE in SCALES[DIMENSION]
				]

				if boxes[0].min() < lo: lo = boxes[0].min()
				if boxes[-1].max() > hi: hi = boxes[-1].max()

				boxplot = ax.boxplot(
					boxes,
					positions=[k+POSITION for k in range(len(boxes))],
					patch_artist=True,
					boxprops={"facecolor": COLOR, **boxprops},
					**plotprops
				)

				for j, box in enumerate(boxes):
					ax.scatter(j+POSITION, np.median(box), marker="D", s=2, color="k" , zorder=1000)

			ax.set_yscale("log")
			ax.set_ylabel(r"$\log_{10}(\mu\mathrm{s})$", fontsize=8)
			ax.set_ylim(lo, hi)

			ax.set_xlabel("cells (linear scale)", fontsize=8)
			ax.set_xticks(range(len(boxes)))
			ax.set_xticklabels(
				[
					f"{cellcounts[DIMENSION][SCALE]:,} ({SCALE})"
					for SCALE in SCALES[DIMENSION]
				],
				fontsize=6,
				rotation=315,
				ha="left"
			)

			ax.set_xlim(-1,len(boxes)+1)

			for tick, label in labels.items():
				if lo <= tick <= hi:
					ax.text(0.05, tick, label, fontsize=6, ha="left", alpha=1/2, va="bottom")
					ax.hlines(tick, xmin=-1, xmax=len(boxes)+1, color="k", alpha=1/4, zorder=-1000)
			plt.savefig(f"./timing/{PREFIX}.jpeg", dpi=300, bbox_inches="tight")

			plt.clf()
			plt.close()


## A plot of matrix density vs time?
