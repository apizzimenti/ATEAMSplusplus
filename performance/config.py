
from sklearn.utils import Bunch
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np

CONFIG = Bunch()
CONFIG.topics = Bunch()

################################################################################
## METADATA ####################################################################
################################################################################
metadata = Bunch()

metadata.hosts = ["pangolin", "meglTower"]
metadata.computing = ["serial", "parallel"]

metadata.initialize = lambda p: p.rcParams.update({
	"text.usetex": True,
	"font.family": "Times New Roman",
	"text.latex.preamble": r"\usepackage{nicefrac}"
})

# Cell counts by dimension and scale.
metadata.cells = {
	2: {
		2: 16,
		4: 64,
		5: 100,
		8: 256,
		11: 484,
		16: 1024,
		22: 1936,
		32: 4096,
		45: 8100,
		64: 16384,
		90: 32400,
		128: 65536,
		181: 131044,
		256: 262144,
		362: 524176,
		512: 1048576
	},
	4: {
		2: 256,
		4: 4096,
		5: 10000,
		8: 65536,
		11: 234256,
		16: 1048576,
		22: 3748096,
		32: 16777216
	}
}

# Boundary matrix density by dimension and scale.
metadata.density = {
	2: {
		2: 0.222222,
		4: 0.0555556,
		5: 0.0355556,
		8: 0.0138889,
		11: 0.00734619,
		16: 0.00347222,
		22: 0.00183655,
		32: 0.000868056,
		45: 0.000438957,
		64: 0.000217014,
		90: 0.000109739,
		128: 5.42535e-05,
		181: 2.71325e-05,
		256: 1.35634e-05,
		362: 6.78313e-06,
		512: 3.39084e-06
	},
	4: {
		2: 0.0177778,
		4: 0.00111111,
		5: 0.000455111,
		8: 6.94444e-05,
		11: 1.94279e-05,
		16: 4.34028e-06,
		22: 1.21425e-06,
		32: 2.71267e-07
	}
}

metadata.logTimeLabels = {
	1e0: "1$\mu$s",
	1e1: "10$\mu$s",
	1e2: "100$\mu$s",
	1e3: "1ms",
	# 0.5e6: r"$\nicefrac 12$s",
	1e4: "10ms",
	1e5: "100ms",
	1e6: "1s",
	3e7: "30s",
	6e7: "1m",
	9e8: "30m",
	3.6e9: "1h",
	4.32e10: "12h",
	8.64e10: "1d",
	6.048e11: "1w"
}

CONFIG.metadata = metadata

################################################################################
## COLORS ######################################################################
################################################################################

colors = Bunch()

tol = Bunch()

tol.bright = Bunch(
	blue="#4477AA",
	cyan="#66CCEE",
	green="#228833",
	yellow="#CCBB44",
	red="#EE6677",
	purple="#AA3377",
	grey="#BBBBBB"
)

tol.highcontrast = Bunch(
	yellow="#DDAA33",
	red="#BB5566",
	blue="#004488"
)

tol.vibrant = Bunch(
	blue="#0077BB",
	cyan="#33BBEE",
	teal="#009988", 
	orange="#EE7733",
	red="#CC3311",
	magenta="#EE3377",
	grey="#BBBBBB"
)

tol.muted = Bunch(
	indigo="#332288",
	cyan="#88CC33",
	teal="#44AA99",
	green="#117733",
	olive="#999933",
	sand="#DDCC77",
	rose="#CC6677",
	wine="#882255",
	purple="#AA4499"
)

colors.tol = tol
CONFIG.colors = colors


################################################################################
## PLOT DEFAULTS ###############################################################
################################################################################

CONFIG._defaults = Bunch()
CONFIG._defaults.yaxis = Bunch()
CONFIG._defaults.xaxis = Bunch()

# Saved figure properties.
CONFIG._defaults.savefig = dict(
	dpi=600,
	bbox_inches="tight"
)

# Fallback scatter plots.
CONFIG._defaults.scatter = dict(
	marker="o",
	color="k",
	s=2,
	linewidths=0
)

# Log-scaled time defaults.
CONFIG._defaults.yaxis.logTime = Bunch()

CONFIG._defaults.yaxis.logTime.textprops = dict(
	fontsize=6,
	ha="left",
	alpha=1/2,
	va="bottom"
)

CONFIG._defaults.yaxis.logTime.lineprops = dict(
	color="k",
	alpha=1/4,
	zorder=-1000
)

CONFIG._defaults.yaxis.logTime = \
	lambda \
		ax, \
		textprops=CONFIG._defaults.yaxis.logTime.textprops, \
		lineprops=CONFIG._defaults.yaxis.logTime.lineprops: \
		_defaultLogTimeVerticalAxis(ax, textprops, lineprops)



################################################################################
## BOXPLOTS ####################################################################
################################################################################

# boxplots = Bunch()
# defaults = Bunch()

# boxplots.defaults = defaults
# CONFIG.plots.boxplots = boxplots

# defaults.offsets = Bunch(
# 	two=lambda w: [-w, w],
# 	three=lambda w: [-2/3*w, 0, 2/3*w],
# 	four=lambda w: [-w, -w/2, w/2, w]
# )

# defaults.subplots = dict(
# 	figsize=(4,4)
# )

# defaults.lw = 2/3

# defaults.flierprops = dict(
# 	marker="x",
# 	markersize=4
# )

# defaults.medianprops = dict(
# 	linewidth=0,
# 	color="none"
# )

# defaults.boxprops = dict(
# 	linewidth=0,
# 	facecolor=f"#{tol.highcontrast.yellow}"
# )

# defaults.whiskerprops = dict(
# 	linewidth=defaults.lw
# )

# defaults.capprops = dict(
# 	linewidth=defaults.lw
# )

# defaults.plotprops = dict(
# 	widths=0.35,
# 	whis=(1,99),
# 	showfliers=False,
# 	patch_artist=True,
# 	boxprops=defaults.boxprops,
# 	capprops=defaults.capprops,
# 	flierprops=defaults.flierprops,
# 	medianprops=defaults.medianprops,
# 	whiskerprops=defaults.whiskerprops
# )

# defaults.savefig = dict(
# 	dpi=600,
# 	bbox_inches="tight"
# )

# defaults.medianscatter = dict(
# 	marker="D",
# 	s=2,
# 	color="k"
# )

# defaults.xaxis = Bunch()

# defaults.xaxis.fontsize = 6
# defaults.xaxis.rotation = -35
# defaults.xaxis.labelpad = -0.5
# defaults.xaxis.ha = "left"

# defaults.xaxis.labelprops = dict(
# 	fontsize=defaults.xaxis.fontsize,
# 	rotation=defaults.xaxis.rotation,
# 	ha=defaults.xaxis.ha
# )

# defaults.xaxis.tickprops = dict(
# 	pad=defaults.xaxis.labelpad
# )


# def _defaultLabeledXaxis(ax, labels, labelprops, tickprops):
# 	ax.set_xlim(1/2, len(labels)+1/2)
# 	ax.set_xticklabels(labels, **labelprops)
# 	ax.tick_params(**tickprops)

# defaults.xaxis.boxedLabeled = lambda ax, labels, labelprops=defaults.xaxis.labelprops, tickprops=defaults.xaxis.tickprops: \
# 	_defaultLabeledXaxis(ax, labels, labelprops, tickprops)


# def _defaultLabeledLogXaxis(ax, base, ticks, labels, labelprops, tickprops):
# 	ax.set_logscale("log", base=base)
# 	ax.set_xlim(1/2, labels[-1]+1/2)
# 	ax.set_xticklabels(labels, **labelprops)
# 	ax.tick_params(**tickprops)

# defaults.yaxis = Bunch()


# # Adjusts the vertical axis limits relative to a set of boxes.
# def _defaultYlim(ax, boxset):
# 	if not isinstance(boxset[0], list): boxset = [boxset]

# 	lo = min(ax.get_ylim()[0], min(boxes[0].min() for boxes in boxset))
# 	hi = min(ax.get_ylim()[1], max(boxes[-1].max() for boxes in boxset))

# 	ax.set_ylim(lo, hi)

# defaults.yaxis.boxedYlim = lambda ax, boxset: _defaultYlim(ax, boxset)


# def _defaultBoxplotMedians(ax, boxes, locs, props):
# 	for loc, box in zip(locs, boxes):
# 		ax.scatter(loc, np.median(box), **props)

# defaults.medians = lambda ax, boxes, locs, props=defaults.medianscatter: _defaultBoxplotMedians(ax, boxes, locs, props)


# def _defaultlogTimeLabels(ax, boxes, textprops, lineprops):
# 	lo, hi = ax.get_ylim()
# 	xlo, _ = ax.get_xlim()

# 	for tick, label in metadata.logTimeLabels.items():
# 		if lo <= tick <= hi:
# 			ax.text(xlo+0.075, tick, label, fontsize=6, ha="left", alpha=1/2, va="bottom")
# 			ax.hlines(tick, xmin=-1, xmax=len(boxes)+1, color="k", alpha=1/4, zorder=-1000)

# defaults._logTime = Bunch()
# defaults._logTime.textprops = dict(
# 	fontsize=6,
# 	ha="left",
# 	alpha=1/2,
# 	va="bottom"
# )

# defaults._logTime.lineprops = dict(
# 	color="k",
# 	alpha=1/4,
# 	zorder=-1000
# )

# defaults.yaxis.logTime = \
# 	lambda ax, boxes, textprops=defaults._logTime.textprops, lineprops=defaults._logTime.lineprops: \
# 		_defaultlogTimeLabels(ax, boxes, textprops, lineprops)

# ################################################################################
# ## ADDITION (BOXPLOT) ##########################################################
# ################################################################################

# addition = Bunch()
# addition.subplots = Bunch()

# boxplots.addition = addition

# addition.exec = "addition"
# addition.trials = 1000

# addition.sep = ","
# addition.prefix = lambda host: f"{host}.{addition.exec}.{addition.trials}"

# addition.columns = [
# 	"N", "lN", "rN", "OVERLAP", "TTC"
# ]

# addition.dtypes = {
# 	"N": int,
# 	"lN": int,
# 	"rN": int,
# 	"OVERLAP": float,
# 	"TTC": int
# }

# def _additionReadFile(host):
# 	pref = addition.prefix(host)
# 	data = pd.read_csv(f"./timing/{pref}.csv", sep=addition.sep, names=addition.columns)

# 	for column in addition.columns: data[column] = data[column].astype(addition.dtypes[column])

# 	return data

# addition.data = lambda host: _additionReadFile(host)

# ###############################
# ## ADDITION, DENSITY BY TIME ##
# ###############################

# densityByTime = Bunch()
# addition.subplots.densityByTime = densityByTime

# densityByTime.ext = "jpeg"
# densityByTime.out = lambda host, length: f"./timing/{addition.prefix(host)}.{length}.densityByTime.{densityByTime.ext}"
# densityByTime.savefig = defaults.savefig



################################################################################
## TOPIC: ADDITION #############################################################
################################################################################

CONFIG.topics.addition = Bunch()

CONFIG.topics.addition.exec = "addition"
CONFIG.topics.addition.computing = CONFIG.metadata.computing
CONFIG.topics.addition.trials = 1000
CONFIG.topics.addition.sep = ","
CONFIG.topics.addition.prefix = lambda host, computing: f"{host}.{CONFIG.topics.addition.exec}.{computing}.{CONFIG.topics.addition.trials}"

CONFIG.topics.addition.columns = [
	"N",
	"lN",
	"rN",
	"OVERLAP",
	"TTC"
]

CONFIG.topics.addition.dtypes = {
	"N": int,
	"lN": int,
	"rN": int,
	"OVERLAP": float,
	"TTC": int
}


CONFIG.topics.addition.data = lambda host, computing: _additionReadFile(host, computing, CONFIG.topics.addition)

################################################################################
## PLOTS #######################################################################
################################################################################

CONFIG.topics.addition.plots = Bunch()
CONFIG.topics.addition.plots.defaults = Bunch()

CONFIG.topics.addition.plots.defaults.subplots = dict(
	figsize=(5,2.5)
)

CONFIG.topics.addition.plots.defaults.scatter = CONFIG._defaults.scatter

################################################################
## ADDITION, DENSITY BY TIME ###################################
################################################################

CONFIG.topics.addition.plots.timeByOverlap = Bunch()

CONFIG.topics.addition.plots.timeByOverlap.scatter = CONFIG.topics.addition.plots.defaults.scatter


## AXES ############
def _timeByOverlapVaxis(ax):
	CONFIG._defaults.yaxis.logTime(ax)

def _timeByOverlapHaxis(ax):
	ax.set_xlim(-0.05,1.05)
	ax.set_xticks([0, 1/4, 1/2, 3/4, 1])

	ax.set_xticklabels([
		f"${t}$" for t in [r"0", r"\nicefrac 14", r"\nicefrac 12", r"\nicefrac 34", r"1"]
	])

CONFIG.topics.addition.plots.timeByOverlap.yaxis = _timeByOverlapVaxis
CONFIG.topics.addition.plots.timeByOverlap.xaxis = _timeByOverlapHaxis

## LSQ #############

CONFIG.topics.addition.plots.timeByOverlap.lsq = Bunch()
CONFIG.topics.addition.plots.timeByOverlap.lsq.f = lambda x, m, b: m*x + b
CONFIG.topics.addition.plots.timeByOverlap.lsq.plot = dict(
	c=CONFIG.colors.tol.vibrant.red,
	zorder=-10000
)

CONFIG.topics.addition.plots.timeByOverlap.lsq.text = dict(
	ha="right",
	va="bottom",
	fontsize=4,
	color=CONFIG.colors.tol.vibrant.red
)

## OUTPUT ##########
CONFIG.topics.addition.plots.timeByOverlap.ext = "jpeg"
CONFIG.topics.addition.plots.timeByOverlap.out = lambda host, length, computing: f"./timing/{CONFIG.topics.addition.prefix(host, computing)}.{length}.timeByOverlap.{CONFIG.topics.addition.plots.timeByOverlap.ext}"
CONFIG.topics.addition.plots.timeByOverlap.savefig = CONFIG._defaults.savefig


################################################################################
## MISCELLANEOUS FUNCTIONS #####################################################
################################################################################
def _additionReadFile(host, computing, topic):
	pref = topic.prefix(host, computing)
	data = pd.read_csv(f"./timing/{pref}.csv", sep=topic.sep, names=topic.columns)

	for column in topic.columns: data[column] = data[column].astype(topic.dtypes[column])

	return data


def _defaultLogTimeVerticalAxis(ax, textprops, lineprops):
	# Fiddle with axis limits; make it look a bit cleaner, perhaps?
	ax.set_yscale("log")
	ylo, yhi = ax.get_ylim()

	loglo = int(np.log10(ylo))
	loghi = int(np.log10(yhi))

	ax.set_yticks([10**k for k in range(loglo, loghi+2)])

	ylo = 10**loglo
	yhi = 10**(loghi+1)
	ax.set_ylim(ylo, yhi+10)

	xlo, xhi = ax.get_xlim()

	for tick, label in CONFIG.metadata.logTimeLabels.items():
		if ylo <= tick <= yhi:
			ax.text(xlo+0.005, tick, label, **textprops)
			ax.plot([xlo, xhi], [tick, tick], **lineprops)


