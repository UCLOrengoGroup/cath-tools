/// \file
/// \brief The read_and_process_mgr class definitions

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

#include "read_and_process_mgr.hpp"

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/min_proj_element.hpp"
#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"
#include "cath/resolve_hits/options/spec/crh_input_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_spec.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.hpp"
#include "cath/resolve_hits/resolve/hit_resolver.hpp"
#include "cath/resolve_hits/scored_hit_arch.hpp"

#include <iostream>
#include <string>
#include <thread>
#include <utility>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::rslv::detail;

using ::std::ostream;
using ::std::string;

// /// \brief
// template <typename Rng, typename Comp, typename Proj>
// size_vec get_ranks(Rng  &&prm_range, ///<
//                    Comp &&prm_comp,  ///<
//                    Proj &&prm_proj   ///<
//                    ) {
// 	const auto comp_fn = [&] (const size_t &x, const size_t &y) {
// 		return forward<Comp>( prm_comp )(
// 			prm_proj( prm_range[ x ] ),
// 			prm_proj( prm_range[ y ] )
// 		);
// 	};

// 	const auto sorted_indices = sort_build<size_vec>(
// 		indices( size( prm_range ) ),
// 		comp_fn
// 	);

// 	size_vec ranks( size( prm_range ), 0_z );

// 	size_t ctr = 0;
// 	for (const size_t &index : indices( size( prm_range ) ) ) {
// 		const size_t &prev_sorted_index  = sorted_indices[ index - 1 ];
// 		const size_t &this_sorted_index  = sorted_indices[ index     ];
// 		if ( index == 0 || comp_fn( prev_sorted_index, this_sorted_index ) ) {
// 			ctr = index + 1;
// 		}
// 		ranks[ this_sorted_index ] = ctr;
// 	}

// 	return ranks;
// }

// /// \todo Make this smarter re joint-positions
// size_t lowest_rank_of_soln_hits(const calc_hit_list &prm_hits, ///<
//                                 const hit_arch &prm_soln  ///<
//                                 ) {
// 	const size_vec hit_ranks = get_ranks( prm_hits, std::greater{}, [] (const calc_hit &x) { return x.get_score(); } );

// 	return min_proj(
// 		prm_soln,
// 		std::less<>{},
// 		[&] (const calc_hit &x) {
// 			const size_t hit_idx = *boost::find_if(
// 				indices( prm_hits.size() ),
// 				[&] (const size_t &y) {
// 					return ( prm_hits[ y ] == x );
// 				}
// 			);
// 			return hit_ranks[ hit_idx ];
// 		}
// 	);
// }

/// \brief Make a read_and_process_mgr from the specified hits_processor, crh_filter_spec and crh_input_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor   &prm_hits_processor,   ///< The hits_processor that should process the hits
                                                           const crh_filter_spec  &prm_filter_spec,      ///< The crh_filter_spec to specify how hits should be filtered
                                                           const crh_score_spec   &prm_crh_score_spec,   ///< The score spec to apply to incoming hits
                                                           const crh_segment_spec &prm_crh_segment_spec, ///< The segment spec to apply to incoming hits
                                                           const crh_input_spec   &prm_input_spec        ///< The crh_input_spec to specify how hits should be read in
                                                           ) {
	return read_and_process_mgr{
		hits_processor_list{ prm_crh_score_spec, prm_crh_segment_spec, { prm_hits_processor.clone() } },
		prm_filter_spec,
		prm_input_spec.get_input_hits_are_grouped()
	};
}

/// \brief Make a read_and_process_mgr from the specified hits_processor and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor &prm_hits_processor, ///< The hits_processor that should process the hits
                                                           const crh_spec       &prm_spec            ///< The crh_spec to specify behaviour
                                                           ) {
	return make_read_and_process_mgr(
		prm_hits_processor,
		prm_spec.get_filter_spec(),
		prm_spec.get_score_spec(),
		prm_spec.get_segment_spec(),
		prm_spec.get_input_spec()
	);
}

/// \brief Make a read_and_process_mgr from the specified hits_processor_list and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor_list &prm_hits_processors, ///< The hits_processor that should process the hits
                                                           const crh_spec            &prm_spec             ///< The crh_spec to specify behaviour
                                                           ) {
	return read_and_process_mgr{
		prm_hits_processors,
		prm_spec.get_filter_spec(),
		prm_spec.get_input_spec().get_input_hits_are_grouped()
	};
}

/// \brief Make a read_and_process_mgr from the specified ostream and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(ofstream_list  &prm_ofstreams, ///< The ofstream_list to which the read_and_process_mgr's hits_processors should write
                                                           const crh_spec &prm_spec       ///< The crh_spec to specify what to do
                                                           ) {
	return make_read_and_process_mgr(
		make_hits_processors(
			prm_ofstreams,
			prm_spec.get_single_output_spec(),
			prm_spec.get_output_spec(),
			prm_spec.get_score_spec(),
			prm_spec.get_segment_spec(),
			prm_spec.get_html_spec()
		),
		prm_spec
	);
}
