
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

from config import CONFIG

metadata = CONFIG.metadata
persistence = CONFIG.topics.persistence
colors = CONFIG.colors

# Specify available data.
persistence.trials = 64

defaults = persistence.plots.defaults

# Initialize the plot.
metadata.initialize(plt)

for host in persistence.hosts:
	for field in persistence.fields:
		for dimension in persistence.dimensions:

			fig, ax = plt.subplots(**persistence.plots.defaults.subplots)
			# positions = [-0.3, -0.15, 0, 0.15, 0.3]
			positions = [-0.3, -0.1, 0.1, 0.3]

			xticks = np.array([])
			# labeled = [False, False, False, False, False]
			labeled = [False]*len(positions)

			for i, (color, pos, name) in enumerate(zip(persistence.colors, positions, persistence.computing)):
				try: data = persistence.data(host, name)
				except: break

				subset = data[(data.DIMENSION == dimension) & (data.FIELD == field)]

				boxes = [
					subset[subset.SCALE == scale].TTC for scale in subset.SCALE.unique()
				]

				# Plot boxes.
				if not len(boxes): continue
				if len(boxes) >= len(xticks): xticks = subset.SCALE.unique()

				posns = np.arange(1, len(boxes)+1) + pos
				ax.boxplot(boxes, positions=posns, **persistence.plots.defaults.boxplot.props(color))

				if not all(labeled):
					ax.text(
						posns[0], np.percentile(boxes[0], 98, axis=None),
						rf".~~{name}",
						fontsize=7,
						rotation=90,
						ha="center",
						va="bottom"
					)
					labeled[i] = True

				# Plot medians.
				medians = np.array([
					np.median(box) for box in boxes
				])

				ax.scatter(posns, medians, **persistence.plots.defaults.boxplot.medians(color))

			ticks = range(1,len(xticks)+1)
			ax.set_xticks(ticks)
			ax.set_xticklabels([f"{t}\n({metadata.cells[dimension][t]:,})" for t in xticks], fontsize=7)
			ax.set_xlim(1/2, len(xticks)+1/2)

			ax.set_title(rf"$\mathbb T^{{{dimension}}}_k$, $\mathbb Z/{{{field}}}\mathbb Z$", fontsize=8)
			# ax.set_title(rf"{host}, $\mathbb T^{{{dimension}}}_k$, $\mathbb Z/{{{field}}}\mathbb Z$", fontsize=8)

			CONFIG._defaults.yaxis.logTime(ax)

			plt.savefig(f"./timing/{host}.{dimension}.{field}.jpeg", **CONFIG._defaults.savefig)
			plt.close()
			plt.clf()

