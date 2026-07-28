#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16 23 32)
DIMENSIONS=(2 4 6)
FIELDS=(2 3 5)
TRIALS=${1:-100}
STRATEGIES=("twist" "stagger" "JIT" "standardParallel" "standard")

HOST=$(hostname -f)

# Test persistence.
for EXEC in "${EXECS[@]}"; do
	for SCALE in "${SCALES[@]}"; do
		for DIMENSION in "${DIMENSIONS[@]}"; do
			for FIELD in "${FIELDS[@]}"; do
				for STRATEGY in "${STRATEGIES[@]}"; do
					./build/timing.$EXEC $HOST $SCALE $DIMENSION $FIELD $TRIALS $STRATEGY

					PADDEDSCALE=${(l(2)( ))SCALE}
					echo "completed $HOST $PADDEDSCALE $DIMENSION $FIELD $TRIALS $STRATEGY"
					# echo $?
				done
			done
		done
	done
done