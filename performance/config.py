
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
## TOPIC: ADDITION #############################################################
################################################################################

CONFIG.topics.addition = Bunch()

CONFIG.topics.addition.exec = "addition"
CONFIG.topics.addition.computing = CONFIG.metadata.computing
CONFIG.topics.addition.hosts = CONFIG.metadata.hosts
CONFIG.topics.addition.trials = 1000
CONFIG.topics.addition.sep = ","
CONFIG.topics.addition.prefix = lambda host, computing: f"{host}.{CONFIG.topics.addition.exec}.{computing}.{CONFIG.topics.addition.trials}"

CONFIG.topics.addition.columns = [
	"N",
	"L",
	"OVERLAP",
	"TTC",
	"CORES"
]

CONFIG.topics.addition.dtypes = {
	"N": int,
	"L": int,
	"OVERLAP": float,
	"TTC": int,
	"CORES": int
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
CONFIG.topics.addition.plots.timeByOverlap.yaxis = lambda ax: _timeByOverlapVaxis(ax)
CONFIG.topics.addition.plots.timeByOverlap.xaxis = lambda ax, lo, hi: _timeByOverlapHaxis(ax, lo, hi)

## LSQ #############
CONFIG.topics.addition.plots.timeByOverlap.lsq = Bunch()
CONFIG.topics.addition.plots.timeByOverlap.lsq.f = lambda x, a, k: a*x**k
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
CONFIG.topics.addition.plots.timeByOverlap.out = lambda host, length, computing, cores: f"./timing/{CONFIG.topics.addition.prefix(host, computing)}.{length}.{cores}.timeByOverlap.{CONFIG.topics.addition.plots.timeByOverlap.ext}"
CONFIG.topics.addition.plots.timeByOverlap.savefig = CONFIG._defaults.savefig


################################################################
## ADDITION, DENSITY BY TIME, OVERLAID #########################
################################################################

CONFIG.topics.addition.plots.timeByOverlapOverlaid = Bunch()

CONFIG.topics.addition.plots.timeByOverlapOverlaid.scatter = CONFIG.topics.addition.plots.defaults.scatter

## AXES ############
CONFIG.topics.addition.plots.timeByOverlapOverlaid.yaxis = CONFIG.topics.addition.plots.timeByOverlap.yaxis
CONFIG.topics.addition.plots.timeByOverlapOverlaid.xaxis = CONFIG.topics.addition.plots.timeByOverlap.xaxis

## LSQ #############
CONFIG.topics.addition.plots.timeByOverlapOverlaid.lsq = Bunch()
CONFIG.topics.addition.plots.timeByOverlapOverlaid.lsq.f = CONFIG.topics.addition.plots.timeByOverlap.lsq.f
CONFIG.topics.addition.plots.timeByOverlapOverlaid.lsq.plot = lambda color: dict(
	c=color,
	zorder=-10000
)

CONFIG.topics.addition.plots.timeByOverlapOverlaid.lsq.text = lambda color: dict(
	ha="right",
	va="bottom",
	fontsize=4,
	color=color
)

## OUTPUT ##########
CONFIG.topics.addition.plots.timeByOverlapOverlaid.ext = "jpeg"

CONFIG.topics.addition.plots.timeByOverlapOverlaid.out = lambda host, length, computing="overlaid": \
	f"./timing/{CONFIG.topics.addition.prefix(host, computing)}.{length}.timeByOverlap.{CONFIG.topics.addition.plots.timeByOverlapOverlaid.ext}"

CONFIG.topics.addition.plots.timeByOverlapOverlaid.savefig = CONFIG._defaults.savefig



























################################################################################
## TOPIC: PERSISTENCE ##########################################################
################################################################################

CONFIG.topics.persistence = Bunch()

CONFIG.topics.persistence.exec = "persistence"
CONFIG.topics.persistence.computing = CONFIG.metadata.computing
CONFIG.topics.persistence.hosts = CONFIG.metadata.hosts
CONFIG.topics.persistence.trials = 100
CONFIG.topics.persistence.sep = ","
CONFIG.topics.persistence.prefix = lambda host, computing: f"{host}.{CONFIG.topics.persistence.exec}.{computing}.{CONFIG.topics.persistence.trials}"

CONFIG.topics.persistence.columns = [
	"SCALE",
	"DIMENSION",
	"TTC",
	"CORES"
]

CONFIG.topics.persistence.dtypes = {
	"SCALE": int,
	"DIMENSION": int,
	"TTC": int,
	"CORES": int
}


CONFIG.topics.persistence.data = lambda host, computing: _additionReadFile(host, computing, CONFIG.topics.persistence)
























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


def _timeByOverlapVaxis(ax):
	CONFIG._defaults.yaxis.logTime(ax)

def _timeByOverlapHaxis(ax, lo, hi):
	ax.set_xlim(lo, hi)
	ax.set_xticks([hi*k for k in [0, 1/4, 1/2, 3/4, 1]])

	ax.set_xticklabels([
		f"${t}$" for t in [r"0", r"\nicefrac 14", r"\nicefrac 12", r"\nicefrac 34", r"1"]
	])
