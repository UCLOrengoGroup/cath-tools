/// \file
/// \brief The hit_seg class definitions

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

#include "hit_seg.h"

#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/range/adaptor/adjacented.h"
#include "exception/invalid_argument_exception.h"

#include <iostream>

using namespace cath::common;
using namespace cath::rslv;

using boost::range::adjacent_find;
using boost::range::sort;
using std::ostream;
using std::pair;
using std::string;

/// \brief TODOCUMENT
void cath::rslv::start_sort_hit_segs(hit_seg_vec &arg_hit_segs ///< TODOCUMENT
                                     ) {
	sort(
		arg_hit_segs,
		hit_seg::get_hit_seg_start_less()
	);
}

/// \brief TODOCUMENT
hit_seg_vec cath::rslv::start_sort_hit_segs_copy(hit_seg_vec arg_hit_segs ///< TODOCUMENT
                                                 ) {
	start_sort_hit_segs( arg_hit_segs );
	return arg_hit_segs;
}

/// \brief TODOCUMENT
hit_seg_vec cath::rslv::make_fragments_of_segments(hit_seg_vec arg_segments ///< TODOCUMENT
                                                   ) {
	start_sort_hit_segs( arg_segments );
	return make_fragments_of_start_sorted_segments( arg_segments );
}

/// \brief TODOCUMENT
bool cath::rslv::segments_are_start_sorted_and_non_overlapping(const hit_seg_vec &arg_segments ///< TODOCUMENT
                                                               ) {
	const auto seg_pair_are_invalid = [] (const hit_seg &x, const hit_seg &y) {
		return ( y.get_start_arrow() < x.get_stop_arrow() );
	};
	return ( adjacent_find( arg_segments, seg_pair_are_invalid ) == common::cend( arg_segments ) );
}

/// \brief TODOCUMENT
hit_seg_vec cath::rslv::make_fragments_of_start_sorted_segments(const hit_seg_vec &arg_segments ///< TODOCUMENT
                                                                ) {
	if ( arg_segments.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot make fragments from empty vector of segments"));
	}

#ifndef NDEBUG
	if ( ! segments_are_start_sorted_and_non_overlapping( arg_segments ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot make fragments from segments with misordered/overlapping invalid segments (start > stop)"));
	}
#endif

	return transform_build<hit_seg_vec>(
		arg_segments | adjacented,
		[] (const pair<hit_seg, hit_seg> &x) {
			return hit_seg{
				x.first .get_stop_arrow(),
				x.second.get_start_arrow()
			};
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates hit_seg
string cath::rslv::to_string(const hit_seg &arg_hit_seg ///< TODOCUMENT
                             ) {
	return "hit_seg["
		+ ::std::to_string( get_start_res_index( arg_hit_seg ) )
		+ "-"
		+ ::std::to_string( get_stop_res_index ( arg_hit_seg ) )
		+ "]";
}

/// \brief TODOCUMENT
///
/// \relates hit_seg
ostream & cath::rslv::operator<<(ostream       &arg_ostream, ///< TODOCUMENT
                                 const hit_seg &arg_hit_seg  ///< TODOCUMENT
                                 ) {
	arg_ostream << to_string( arg_hit_seg );
	return arg_ostream;
}


// /// \brief TODOCUMENT
// res_arrow_res_arrow_pair_vec get_hit_segs(const hit_vec &arg_hits ///< TODOCUMENT
//                                           ) {
// 	res_arrow_res_arrow_pair_vec segments;
// 	for (const hit &the_hit : arg_hits) {
// 		for (const auto &x : irange( 0_z, the_hit.get_num_segments() ) ) {
// 			hit_segs.emplace_back(
// 				the_hit.get_start_arrow_of_segment( seg_ctr ),
// 				the_hit.get_stop_arrow_of_segment ( seg_ctr )
// 			)
// 		}
// 	}
// 	return segments;
// }

// /// \brief TODOCUMENT
// res_arrow_res_arrow_pair_vec get_start_sorted_hit_segs(const hit_vec &arg_hits ///< TODOCUMENT
//                                                        ) {
// 	auto segments = get_hit_segs( arg_hits );
// 	sort( segments, [] (const res_arrow_res_arrow_pair_vec &x, const res_arrow_res_arrow_pair_vec &y) { return ( x.first < y.first); } );
// 	return segments;
// }

// /// \brief TODOCUMENT
// residx_residx_pair_vec get_start_sorted_segments_of_hits(const hit_vec &arg_hits ///< TODOCUMENT
//                                                          ) {
// 	residx_residx_pair_vec segments;
// 	for (const hit &the_hit : arg_hits) {
// 		for (const auto &seg_ctr : irange( 0_z, the_hit.get_num_segments() ) ) {
// 			segments.emplace_back(
// 				the_hit.get_start_of_segment( seg_ctr ),
// 				the_hit.get_stop_of_segment( seg_ctr )
// 			);
// 		}
// 	}
// 	sort( segments, [] (const residx_residx_pair &x, const residx_residx_pair &y) { return ( x.first < y.first); } );
// 	return segments;
// }


// {  }
// res_before( )

// { { 30, 35 } }
// res_before( )

// { { 30, 35 }, { 40, 45 } }
// res_before( )


// /// \brief TODOCUMENT
// res_arrow_res_arrow_pair_vec get_free_space_of_fitting_hits_in_region(const hit_vec   &arg_hits,         ///< TODOCUMENT
//                                                                       const res_arrow &arg_region_start, ///< TODOCUMENT
//                                                                       const res_arrow &arg_region_stop   ///< TODOCUMENT
//                                                                       ) {
// 	const auto arg_region_start = 
// 	if ( ! arg_hits.empty() ) {
// 		BOOST_THROW_EXCEPTION(invalid_argument_exception("Arghh"));
// 	}
// 	return { { arg_region_start, arg_region_stop } };
// 	// }
// 	// el
// 	// const auto start_sorted_hit_segs = get_start_sorted_hit_segs( arg_hits );
// }
