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

#include "common/boost_addenda/range/min_proj_element.hpp"
#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/hit_resolver.hpp"
#include "resolve_hits/options/spec/crh_input_spec.hpp"
#include "resolve_hits/options/spec/crh_spec.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.hpp"
#include "resolve_hits/scored_hit_arch.hpp"

#include <iostream>
#include <string>
#include <thread>
#include <utility>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::irange;
using boost::size;
using std::forward;
using std::future;
using std::launch;
using std::move;
using std::ostream;
using std::ref;
using std::string;

constexpr bool read_and_process_mgr::DEFAULT_INPUT_HITS_ARE_GROUPED;

// /// \brief
// template <typename Rng, typename Comp, typename Proj>
// size_vec get_ranks(Rng  &&arg_range, ///<
//                    Comp &&arg_comp,  ///<
//                    Proj &&arg_proj   ///<
//                    ) {
// 	const auto comp_fn = [&] (const size_t &x, const size_t &y) {
// 		return forward<Comp>( arg_comp )(
// 			arg_proj( arg_range[ x ] ),
// 			arg_proj( arg_range[ y ] )
// 		);
// 	};

// 	const auto sorted_indices = sort_build<size_vec>(
// 		irange( 0_z, size( arg_range ) ),
// 		comp_fn
// 	);

// 	size_vec ranks( size( arg_range ), 0_z );

// 	size_t ctr = 0;
// 	for (const size_t &index : irange( 0_z, size( arg_range ) ) ) {
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
// size_t lowest_rank_of_soln_hits(const calc_hit_list &arg_hits, ///<
//                                 const hit_arch &arg_soln  ///<
//                                 ) {
// 	const size_vec hit_ranks = get_ranks( arg_hits, std::greater{}, [] (const calc_hit &x) { return x.get_score(); } );

// 	return min_proj(
// 		arg_soln,
// 		std::less<>{},
// 		[&] (const calc_hit &x) {
// 			const size_t hit_idx = *boost::find_if(
// 				irange( 0_z, arg_hits.size() ),
// 				[&] (const size_t &y) {
// 					return ( arg_hits[ y ] == x );
// 				}
// 			);
// 			return hit_ranks[ hit_idx ];
// 		}
// 	);
// }

/// \brief Make a read_and_process_mgr from the specified hits_processor, crh_filter_spec and crh_input_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor   &arg_hits_processor,   ///< The hits_processor that should process the hits
                                                           const crh_filter_spec  &arg_filter_spec,      ///< The crh_filter_spec to specify how hits should be filtered
                                                           const crh_score_spec   &arg_crh_score_spec,   ///< The score spec to apply to incoming hits
                                                           const crh_segment_spec &arg_crh_segment_spec, ///< The segment spec to apply to incoming hits
                                                           const crh_input_spec   &arg_input_spec        ///< The crh_input_spec to specify how hits should be read in
                                                           ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return read_and_process_mgr{
		hits_processor_list{ arg_crh_score_spec, arg_crh_segment_spec, { arg_hits_processor.clone() } },
		arg_filter_spec,
		arg_input_spec.get_input_hits_are_grouped()
	};
}

/// \brief Make a read_and_process_mgr from the specified hits_processor and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor &arg_hits_processor, ///< The hits_processor that should process the hits
                                                           const crh_spec       &arg_spec            ///< The crh_spec to specify behaviour
                                                           ) {
	return make_read_and_process_mgr(
		arg_hits_processor,
		arg_spec.get_filter_spec(),
		arg_spec.get_score_spec(),
		arg_spec.get_segment_spec(),
		arg_spec.get_input_spec()
	);
}

/// \brief Make a read_and_process_mgr from the specified hits_processor_list and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(const hits_processor_list &arg_hits_processors, ///< The hits_processor that should process the hits
                                                           const crh_spec            &arg_spec             ///< The crh_spec to specify behaviour
                                                           ) {
	return read_and_process_mgr{
		arg_hits_processors,
		arg_spec.get_filter_spec(),
		arg_spec.get_input_spec().get_input_hits_are_grouped()
	};
}

/// \brief Make a read_and_process_mgr from the specified ostream and crh_spec
///
/// \relates read_and_process_mgr
read_and_process_mgr cath::rslv::make_read_and_process_mgr(ostream        &arg_ostream, ///< The ostream to which the results should be written
                                                           const crh_spec &arg_spec     ///< The crh_spec to specify what to do
                                                           ) {
	return make_read_and_process_mgr(
		make_hits_processors(
			arg_ostream,
			arg_spec.get_single_output_spec(),
			arg_spec.get_score_spec(),
			arg_spec.get_segment_spec(),
			arg_spec.get_html_spec()
		),
		arg_spec
	);
}
