
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

from config import CONFIG

metadata = CONFIG.metadata
addition = CONFIG.topics.addition

defaults = addition.plots.defaults
timeByOverlap = addition.plots.timeByOverlap

# Initialize the plot.
metadata.initialize(plt)

for host in ["pangolin"]:
	for computing in metadata.computing:
		data = addition.data(host, computing)

		lengths = data.N.unique()

		for length in lengths:
			subset = data[data.N == length]

			# time-by-density scatterplot.
			fig, ax = plt.subplots(**defaults.subplots)
			ax.scatter(subset.OVERLAP, subset.TTC, **timeByOverlap.scatter)

			# Prepare horizontal and vertical axes.
			timeByOverlap.xaxis(ax)
			timeByOverlap.yaxis(ax)

			# Perform a linear least squares fit.
			popt, pcov = curve_fit(timeByOverlap.lsq.f, subset.OVERLAP, subset.TTC)
			m, b = popt

			X = np.linspace(-1, 2, 1000)
			Y = m*X + b
			ax.plot(X, Y, **timeByOverlap.lsq.plot)

			xlo, xhi = ax.get_xlim()
			ax.text(xhi, m*xhi+b, rf"$\mathrm{{O}}(d)$", **timeByOverlap.lsq.text)

			plt.savefig(
				timeByOverlap.out(host, length, computing),
				**timeByOverlap.savefig
			)
			plt.close()
			plt.clf()
