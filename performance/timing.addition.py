
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

from config import CONFIG

metadata = CONFIG.metadata
addition = CONFIG.topics.addition

# Specify available data.
addition.trials = 100

defaults = addition.plots.defaults
timeByOverlap = addition.plots.timeByOverlap

# Initialize the plot.
metadata.initialize(plt)

for host in addition.hosts:
	for computing in addition.computing:
		data = addition.data(host, computing)

		lengths = data.N.unique()
		corecounts = data.CORES.unique()

		for cores in corecounts:
			for length in lengths:
				toosmall = False

				subset = data[(data.N == length) & ((data.CORES == cores) ^ (data.CORES == 1))]
				if len(subset) < 1: toosmall = True; break

				# time-by-density scatterplot.
				fig, ax = plt.subplots(**defaults.subplots)
				ax.scatter(subset.OVERLAP*subset.L, subset.TTC, **timeByOverlap.scatter)

				# Perform a linear least squares fit.
				popt, pcov = curve_fit(timeByOverlap.lsq.f, subset.OVERLAP*subset.L, subset.TTC)
				z, k = popt

				X = np.linspace(1, subset.L.max(), 1000)
				Y = timeByOverlap.lsq.f(X, z, k)

				ax.plot(X, Y, **timeByOverlap.lsq.plot)

				ax.annotate(
					rf"${z:.2f}x^{{{k:.2f}}}$",
					xy=(subset.L.max(), timeByOverlap.lsq.f(subset.L.max()*1.05, z, k)),
					xycoords="data",
					**timeByOverlap.lsq.text
				)

				ax.annotate(
					f"$N=$ {int(addition.trials):,}\n$L = $ {int(length):,}",
					xy=(0.995, 0.01),
					xycoords="axes fraction",
					ha="right",
					va="bottom",
					fontsize=6,
					alpha=1/2,
					bbox=dict(fc="w",lw=0,boxstyle="square,pad=0.15")
				)

				# Prepare horizontal and vertical axes.
				timeByOverlap.xaxis(ax, 0, length)
				timeByOverlap.yaxis(ax)

				if not toosmall:
					plt.savefig(
						timeByOverlap.out(host, length, computing, cores),
						**timeByOverlap.savefig
					)
				
				plt.close()
				plt.clf()


# Now, plot the datasets against one another.
overlaid = addition.plots.timeByOverlapOverlaid


for host in addition.hosts:
	parallel = addition.data(host, "parallel")
	serial = addition.data(host, "serial")

	# These should be the same for both datasets.
	lengths = parallel.N.unique()
	corecounts = parallel.CORES.unique()

	for length in lengths:
		fig, ax = plt.subplots(**defaults.subplots)

		# Plot the serial data once.
		subsetserial = serial[(serial.N == length)]

		scatterargs = overlaid.scatter
		scatterargs.update(dict(color=CONFIG.colors.tol.highcontrast.yellow))
		ax.scatter(subsetserial.OVERLAP*length, subsetserial.TTC, **scatterargs)

		popt, pcov = curve_fit(overlaid.lsq.f, subsetserial.OVERLAP*length, subsetserial.TTC)
		z, k = popt

		X = np.linspace(0, length, 1000)
		Y = overlaid.lsq.f(X, z, k)

		ax.plot(X, Y, **overlaid.lsq.plot(CONFIG.colors.tol.highcontrast.yellow))

		ax.annotate(
			f"1 core",
			xy=(length, overlaid.lsq.f(length*0.85, z, k)),
			xycoords="data",
			color=CONFIG.colors.tol.highcontrast.yellow,
			ha="right",
			va="top",
			fontsize=6,
		)

		# Plot the data for all the core values we have.
		subsetparallel = parallel[(parallel.N == length)]
		corecounts = subsetparallel.CORES.unique()

		for cores in corecounts:
			subsubset = subsetparallel[subsetparallel.CORES == cores]
			scatterargs = overlaid.scatter
			scatterargs.update(dict(color=CONFIG.colors.tol.highcontrast.blue))
			
			ax.scatter(subsubset.OVERLAP*length, subsubset.TTC, **scatterargs)

			popt, pcov = curve_fit(overlaid.lsq.f, subsubset.OVERLAP*length, subsubset.TTC)
			z, k = popt

			X = np.linspace(0, length, 1000)
			Y = overlaid.lsq.f(X, z, k)

			ax.plot(X, Y, **overlaid.lsq.plot(CONFIG.colors.tol.highcontrast.blue))

			ax.annotate(
				f"{cores} core{'s' if cores > 1 else ''}",
				xy=(length, overlaid.lsq.f(length*1.05, z, k)),
				xycoords="data",
				color=CONFIG.colors.tol.highcontrast.blue,
				ha="right",
				va="top",
				fontsize=4,
				# bbox=dict(fc="w",lw=0,boxstyle="square,pad=0")
			)

		# Prepare axes.
		overlaid.xaxis(ax, 0, length)
		overlaid.yaxis(ax)

		ax.annotate(
			f"$N = $ {int(length):,}",
			xy=(0.995, 0.02),
			xycoords="axes fraction",
			ha="right",
			va="bottom",
			fontsize=6,
			alpha=1/2,
			bbox=dict(fc="w",lw=0,boxstyle="square,pad=0.15")
		)

		plt.savefig(
			overlaid.out(host, length),
			**overlaid.savefig
		)
		
		plt.close()
		plt.clf()
