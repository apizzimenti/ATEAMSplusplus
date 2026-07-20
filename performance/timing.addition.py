
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

from config import CONFIG

metadata = CONFIG.metadata

boxplots = CONFIG.plots.boxplots
addition = boxplots.addition

metadata.initialize(plt)

# First, within each host and length, plot density against time.
densityByTime = addition.subplots.densityByTime

for host in metadata.hosts:
	data = boxplots.addition.data(host)
	base = data.BASE[0]

	lengths = data.LENGTH.unique()

	for length in lengths:
		subset = data[data.LENGTH == length]
		grouped = subset.groupby("POWER")

		entries = np.array([int(base**power*length) for power, _ in grouped])
		boxes = [group.TTC for _, group in grouped]

		# Viewing only the last 12-16 boxes.
		window = 16
		rightward = 1
		ll = len(boxes)-window-rightward
		ul = len(boxes)-rightward

		viewboxes = boxes[ll:ul]
		viewentries = entries[ll:ul]

		# Boxplot; plot medians.
		fig, ax = plt.subplots(**boxplots.defaults.subplots)
		ax.boxplot(
			viewboxes,
			**boxplots.defaults.plotprops
		)

		boxplots.defaults.medians(ax, viewboxes, ax.get_xticks())

		# Configure horizontal axis.
		boxplots.defaults.xaxis.boxedLabeled(
			ax,
			[f"{e:,}" for e in viewentries]
		)

		# Configure vertical axis.
		ax.set_yscale("log")
		boxplots.defaults.yaxis.boxedYlim(ax, viewboxes)
		boxplots.defaults.yaxis.logTime(ax, viewboxes)

		# Plot some guide ranges.
		X = np.arange(0, window+2)
		rescale = 1/100

		Y = entries[ll-1:ul+1]*rescale
		wiggle = np.sqrt(10)

		# ax.plot(X, Y, lw=1/2, color="k", zorder=-10000)

		ax.fill_between(
			X, Y*(1/wiggle), Y*wiggle,
			alpha=1/10,
			color="k",
			zorder=-10000,
			edgecolor="none"
		)


		plt.savefig(
			densityByTime.out(host, length),
			**densityByTime.savefig
		)
		plt.close()
		plt.clf()

