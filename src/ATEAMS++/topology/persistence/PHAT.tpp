
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PHAT_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_PHAT_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PHAT_H
#error __FILE__ should only be included from topology/persistence/PHAT.h.
#endif

#include <phat/compute_persistence_pairs.h>
#include <phat/helpers/misc.h>

using namespace std;

namespace ATEAMS::topology::persistence {

	template <typename RingLike>
	inline vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration
	) {
		// The filtration specifies the order in which we add the cells of all
		// dimensions. Create a map that specifies to which position each cell was
		// moved. For example, if the filtration has
		//
		//				[... 12, 9, 10, 19, ...]
		//
		// that was originally
		//				[... 9, 10, 11, 12, ...]
		//
		// then the mapping should have entries
		//
		//				{... 12: 9, 9: 10, 10: 11, 19: 12, ...}

		// TODO optimize, this takes way too long
		vector<int> remapping(filtration.size(), 0);
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		PHATBoundaryMatrix boundary;
		PHATColumn column;

		vector<int> ground;
		int dim;

		boundary.set_num_cols(complex->Boundary.Flat.size());

		// this also takes an insane amount of time; no need to go all the way
		// up.
		for (int t=0; t < filtration.size(); t++) {
			ground = complex->Boundary.Flat[filtration[t]];

			for (int j=0; j < ground.size(); j++) {
				column.push_back(remapping[ground[j]]);
			}

			std::sort(column.begin(), column.end());
		
			// Determine the dimension of the cell.
			for (int d=0; d < complex->Breaks.size(); d++) {
				int lo = complex->Breaks[d][0];
				int hi = complex->Breaks[d][1];

				if (lo <= filtration[t] && filtration[t] < hi) {
					dim = d;
					break;
				}
			}

			boundary.set_dim(t, dim);
			boundary.set_col(t, column);
			column.clear();
		}

		// Compute the persistence pairs and populate a Vector to return to the
		// user. We know when the giant(/essential) cycles are born because they
		// do not appear in the list of births or deaths.
		PHATPairs pairs;
		set<int> all;
		phat::compute_persistence_pairs<PHATTwist>(pairs, boundary);
		pairs.sort();

		for (phat::index i=0; i < pairs.get_num_pairs(); i++) {
			all.insert(pairs.get_pair(i).first);
			all.insert(pairs.get_pair(i).second);
		}

		vector<int> essential;

		for (int t=0; t < filtration.size(); t++) {
			if (!all.contains(filtration[t])) essential.push_back(filtration[t]);
		}

		return essential;
	}



	template <typename RingLike>
	inline vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension
	) {
		// No sense in doing a ton of policy work for this one.
		vector<int> essential = PHAT<RingLike>(complex, filtration);
		int low = complex->Breaks[dimension][0], high = complex->Breaks[dimension][1];
		std::erase_if(essential, [filtration,low,high](int t) { return !(low <= filtration[t] && filtration[t] < high); });

		return essential;
	}
}

#endif
