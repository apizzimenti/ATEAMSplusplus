
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

from config import CONFIG

metadata = CONFIG.metadata
persistence = CONFIG.topics.persistence

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

	# sgrouped = standard.groupby(["SCALE", "DIMENSION"])
	# tgrouped = twist.groupby(["SCALE", "DIMENSION"])

	for dimension in dimensions:
		sgroup = standard[standard.DIMENSION == dimension]
		tgroup = twist[twist.DIMENSION == dimension]

		sboxes = [
			sgroup[sgroup.SCALE == scale].TTC for scale in sgroup.SCALE.unique()
		]

		tboxes = [
			tgroup[tgroup.SCALE == scale].TTC for scale in tgroup.SCALE.unique()
		]
			
		fig, ax = plt.subplots()

		ax.boxplot(sboxes, widths=0.2, positions=[k-0.25 for k in range(1, len(sboxes)+1)], patch_artist=True, boxprops=dict(facecolor="k"))
		ax.boxplot(tboxes, widths=0.2, positions=[k+0.25 for k in range(1, len(tboxes)+1)])

		ax.set_xticks(range(1,len(sboxes)+1))
		ax.set_xticklabels(sgroup.SCALE.unique())
		ax.set_title(f"dim {dimension}")

		ax.set_yscale("log")

		plt.savefig(f"./timing/{host}.{dimension}.jpeg")
		# sys.exit(1)

