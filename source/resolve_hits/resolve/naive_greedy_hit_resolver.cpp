/// \file
/// \brief The naive_greedy_hit_resolver class definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "naive_greedy_hit_resolver.hpp"

#include "common/algorithm/sort_uniq_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "resolve_hits/algo/scored_arch_proxy_fn.hpp"
#include "resolve_hits/scored_hit_arch.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

/// \brief Naively, greedily resolve hits
scored_hit_arch cath::rslv::naive_greedy_resolve_hits(const calc_hit_list &prm_hits ///< The hits to resolve
                                                      ) {
	scored_arch_proxy results;

	// Get a vector of the indices of the hits in descending order of score
	const auto hit_indices_best_to_worst = sort_build<hitidx_vec>(
		indices( prm_hits.size() ),
		[&] (const hitidx_t &x, const hitidx_t &y) {
			return ( prm_hits[ y ].get_score() < prm_hits[ x ].get_score() );
		}
	);

	// Work through the hit indices and add any that don't clash with those already added
	for (const hitidx_t &hit_index : hit_indices_best_to_worst) {
		const calc_hit &the_hit = prm_hits[ hit_index ];
		const resscr_t &score   = the_hit.get_score();
		add_hit_if_does_not_overlap( results, score, hit_index, prm_hits );
	}

	// Build a scored_hit_arch of the results
	return make_scored_hit_arch( results, prm_hits );
}
