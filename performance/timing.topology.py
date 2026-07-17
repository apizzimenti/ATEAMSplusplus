
import numpy as np
import matplotlib.pyplot as plt

HOST = "pangolin"
EXECS = ["persistence"]
FIELDS = [2, 3]
DIMENSIONS = [2, 4]
ATTEMPTS=51

SCALES = {
	2: [2, 4, 5, 8, 11, 16, 22, 32, 45, 64, 90, 128, 181, 256, 362, 512],
	4: [2, 4, 5, 8]
}

UNIT = 1

cellcounts = {
	2: { 2: 16, 4: 64, 5: 100, 8: 256, 11: 484, 16: 1024, 22: 1936, 32: 4096, 45: 8100, 64: 16384, 90: 32400, 128: 65536,  },
	4: { 2: 256, 4: 4096, 5: 10000, 8: 65536, 11: 234256, 16: 1048576, 22: 3748096, 32: 16777216,  }
}

density = {
	2: { 2: 0.222222, 4: 0.0555556, 5: 0.0355556, 8: 0.0138889, 11: 0.00734619, 16: 0.00347222, 22: 0.00183655, 32: 0.000868056, 45: 0.000438957, 64: 0.000217014, 90: 0.000109739, 128: 5.42535e-05,  },
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

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman",
	"text.latex.preamble": r"\usepackage{nicefrac}"
})

for EXEC in EXECS:
	for FIELD in FIELDS:
		for DIMENSION in DIMENSIONS:
			boxes = [
				np.loadtxt(f"timing/{HOST}.{EXEC}.{str(SCALE).zfill(3)}.{DIMENSION}.{FIELD}.{ATTEMPTS}.data", dtype=float)[:-1]/UNIT
				for SCALE in SCALES[DIMENSION]
			]

			lo, hi = boxes[0].min(), boxes[-1].max()

			fig, ax = plt.subplots(figsize=(4,4))
			boxplot = ax.boxplot(
				boxes,
				# positions=[cellcounts[DIMENSION][SCALE] for SCALE in SCALES[DIMENSION]],
				**plotprops
			)

			ax.set_yscale("log")
			ax.set_ylabel(r"$\log_{10}(\mu\mathrm{s})$", fontsize=8)
			ax.set_ylim(lo, hi)

			# ax.set_xlabel("linear scale\n(cells)\n(density)", fontsize=8)
			# ax.set_xticks([cellcounts[DIMENSION][SCALE] for SCALE in SCALES])
			# ax.set_xticklabels(
			# 	[
			# 		f"{SCALE}\n({cellcounts[DIMENSION][SCALE]:,})\n({density[DIMENSION][SCALE]*100:.3f}\%)"
			# 		for SCALE in SCALES
			# 	],
			# 	fontsize=6
			# )

			ax.set_xlabel("cells (linear scale)", fontsize=8)

			ax.set_xticklabels(
				[
					f"{cellcounts[DIMENSION][SCALE]:,} ({SCALE})"
					for SCALE in SCALES[DIMENSION]
				],
				fontsize=6,
				rotation=315,
				ha="left"
			)

			ax.set_xlim(0,len(boxes)+1)
			# ax.set_xscale("log", base=np.sqrt(2))

			for tick, label in labels.items():
				if lo <= tick <= hi:
					ax.text(0.05, tick, label, fontsize=6, ha="left", alpha=1/2, va="bottom")
					ax.hlines(tick, xmin=0, xmax=len(boxes)+1, color="k", alpha=1/4, zorder=-1000)
			plt.savefig(f"./timing/{EXEC}.{DIMENSION}.{FIELD}.{ATTEMPTS}.jpeg", dpi=300, bbox_inches="tight")

			plt.clf()
			plt.close()


## A plot of matrix density vs time?
