#!/bin/zsh

SCALES=(2 4 5 8 11)
DIMENSIONS=(6)
TRIALS=${1:-100}

for SCALE in "${SCALES[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			./build/performance.samples $SCALE $DIMENSION $TRIALS
			echo "completed $SCALE $DIMENSION $TRIALS"
		done
	done
done