
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

from config import CONFIG

metadata = CONFIG.metadata
addition = CONFIG.topics.addition

defaults = addition.plots.defaults
timeByDensity = addition.plots.timeByDensity

# Initialize the plot.
metadata.initialize(plt)

for host in metadata.hosts:
	data = addition.data(host)

	lengths = data.N.unique()

	for length in lengths:
		subset = data[data.N == length]
		print(host, length, len(subset))

		# time-by-density scatterplot.
		fig, ax = plt.subplots(**defaults.subplots)
		ax.scatter(subset.OVERLAP, subset.TTC, **timeByDensity.scatter)

		# Prepare horizontal and vertical axes.
		timeByDensity.xaxis(ax)
		timeByDensity.yaxis(ax)

		# Perform a linear least squares fit.
		popt, pcov = curve_fit(timeByDensity.lsq.f, subset.OVERLAP, subset.TTC)
		m, b = popt

		X = np.linspace(-1, 2, 1000)
		Y = m*X + b
		ax.plot(X, Y, **timeByDensity.lsq.plot)

		xlo, xhi = ax.get_xlim()
		ax.text(xhi, m*xhi+b, rf"$\mathrm{{O}}(d)$", **timeByDensity.lsq.text)

		plt.savefig(
			timeByDensity.out(host, length),
			**timeByDensity.savefig
		)
		plt.close()
		plt.clf()

