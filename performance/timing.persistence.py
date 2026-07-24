
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
persistence.trials = 100

defaults = persistence.plots.defaults
# timeByOverlap = persistence.plots.timeByOverlap

# Initialize the plot.
metadata.initialize(plt)

for host in persistence.hosts:
	standard = persistence.data(host, "standard")
	twist = persistence.data(host, "twist")

	scales = standard.SCALE.unique()
	dimensions = standard.DIMENSION.unique()
	fields = standard.FIELD.unique()

	# sgrouped = standard.groupby(["SCALE", "DIMENSION"])
	# tgrouped = twist.groupby(["SCALE", "DIMENSION"])

	for field in fields:
		for dimension in dimensions:
			sgroup = standard[(standard.DIMENSION == dimension) & (standard.FIELD == field)]
			tgroup = twist[(twist.DIMENSION == dimension) & (twist.FIELD == field)]

			sboxes = [
				sgroup[sgroup.SCALE == scale].TTC for scale in sgroup.SCALE.unique()
			]

			tboxes = [
				tgroup[tgroup.SCALE == scale].TTC for scale in tgroup.SCALE.unique()
			]

			if len(tboxes) != len(sboxes): continue
				
			fig, ax = plt.subplots(figsize=(2.5,4))

			OFFSET = 0.15
			pos = np.arange(1, len(sboxes)+1)
			L = pos-OFFSET
			R = pos+OFFSET

			ax.boxplot(sboxes, positions=L, **persistence.plots.defaults.boxplot.props(colors.tol.highcontrast.red))
			ax.boxplot(tboxes, positions=R, **persistence.plots.defaults.boxplot.props(colors.tol.highcontrast.yellow))

			ax.set_xticks(range(1,len(sboxes)+1))
			ax.set_xticklabels(sgroup.SCALE.unique())
			ax.set_title(rf"{host}, $\mathbb Z^{{{dimension}}}$, $\mathbb Z/{{{field}}}\mathbb Z$", fontsize=8)

			CONFIG._defaults.yaxis.logTime(ax)

			plt.savefig(f"./timing/{host}.{dimension}.{field}.jpeg", **CONFIG._defaults.savefig)
			plt.close()
			plt.clf()

