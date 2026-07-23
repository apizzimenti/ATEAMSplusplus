#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16)
DIMENSIONS=(4)
FIELDS=(3)
TRIALS=${1:-100}
STRATEGIES=("standard" "twist")
PARALLELS=(1)

HOST=$(hostname -f)

# Test persistence.
for EXEC in "${EXECS[@]}"; do
	for SCALE in "${SCALES[@]}"; do
		for DIMENSION in "${DIMENSIONS[@]}"; do
			for FIELD in "${FIELDS[@]}"; do
				for STRATEGY in "${STRATEGIES[@]}"; do
					for PARALLEL in "${PARALLELS[@]}"; do
						./build/timing.$EXEC $HOST $SCALE $DIMENSION $FIELD $TRIALS $STRATEGY $PARALLEL
						echo "completed $HOST $SCALE $DIMENSION $FIELD $TRIALS $STRATEGY $PARALLEL"
						# echo $?
					done
				done
			done
		done
	done
done