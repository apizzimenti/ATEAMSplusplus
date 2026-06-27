#!/bin/zsh

MODELS=("Bernoulli" "Glauber" "InvadedCluster" "InvasionPercolation" "SwendsenWang")
SCALES=(3 4 6 8 11 16 23)
FIELDS=(2 3 5 7)
TOPDIMENSIONS=(2 4)
ITERATIONS=100

EXECPREFIX="../build"

for MODEL in "${MODELS[@]}"; do
	for FIELD in "${FIELDS[@]}"; do
		for SCALE in "${SCALES[@]}"; do
			for TOPDIMENSION in "${TOPDIMENSIONS[@]}"; do
				# Delete existing ones first.
				PADDED=${(l(2)(0))SCALE}
				prefix="$MODEL.$FIELD.$PADDED.$TOPDIMENSION"
				rm profiles/reports/$prefix.*

				echo $prefix

				# Record basic statistics.
				perf stat -o profiles/reports/$prefix.stat ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!

				# Find bottlenecks.
				perf record -o profiles/reports/$prefix.record ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!

				echo \n\n

				# perf report --input=profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.record --output=profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.record.csv &
				# wait $!
			done
		done
	done
done
