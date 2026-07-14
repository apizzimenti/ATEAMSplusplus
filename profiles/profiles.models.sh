#!/bin/zsh

MODELS=("Bernoulli" "Glauber" "InvadedCluster" "Invasion" "SwendsenWang")
SCALES=(3 4)
FIELDS=(2 3)
TOPDIMENSIONS=(2 4)
ITERATIONS=1000

EXECPREFIX="../build"

root="ATEAMS::models::Bernoulli::sample"

for MODEL in "${MODELS[@]}"; do
	for FIELD in "${FIELDS[@]}"; do
		for SCALE in "${SCALES[@]}"; do
			for TOPDIMENSION in "${TOPDIMENSIONS[@]}"; do
				# Delete existing ones first.
				PADDED=${(l(2)(0))SCALE}
				prefix="$MODEL.$FIELD.$PADDED.$TOPDIMENSION"
				rm -f profiles/reports/$prefix.*

				echo "########################################"
				echo "## $prefix"
				echo "########################################"

				# Record basic statistics.
				perf stat -o profiles/reports/$prefix.stat ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!

				# Find bottlenecks.
				perf record -g -o profiles/reports/$prefix.record ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!

				# Make a plot?
				perf script --input=profiles/reports/$prefix.record \
					| c++filt \
					| gprof2dot --format=perf --depth=3 --strip --node-thres=5 \
					| dot -Tpng -o profiles/reports/$prefix.png

				echo "\n\n"
				# perf report --input=profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.record --output=profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.record.csv &
				# wait $!
			done
		done
	done
done
