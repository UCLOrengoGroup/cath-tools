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

/// \brief Make the (sorted) fragments between the specified segments
///
/// \relates hit_seg
hit_seg_vec cath::rslv::make_fragments_of_segments(hit_seg_vec arg_segments ///< The segments from which the fragments should be made
                                                   ) {
	start_sort_hit_segs( arg_segments );
	return make_fragments_of_start_sorted_segments( arg_segments );
}

/// \brief Return whether the specified segments are sorted by their starts and non-overlapping
///
/// \relates hit_seg
bool cath::rslv::segments_are_start_sorted_and_non_overlapping(const hit_seg_vec &arg_segments ///< The segments to query
                                                               ) {
	const auto seg_pair_are_invalid = [] (const hit_seg &x, const hit_seg &y) {
		return ( y.get_start_arrow() < x.get_stop_arrow() );
	};
	return ( adjacent_find( arg_segments, seg_pair_are_invalid ) == common::cend( arg_segments ) );
}

/// \brief Make the fragments between the specified segments
///
/// \pre The segments must be pre-sorted by their starts
///
/// \relates hit_seg
hit_seg_vec cath::rslv::make_fragments_of_start_sorted_segments(const hit_seg_vec &arg_segments ///< The segments from which the fragments must be made, pre-sorted by their starts
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

/// \brief Return whether the midpoint of the first specified hit_seg is less than that of the second
///
/// \relates hit_seg
bool cath::rslv::midpoint_less(const hit_seg &arg_hit_seg_lhs, ///< The first  hit_seg to compare
                               const hit_seg &arg_hit_seg_rhs  ///< The second hit_seg to compare
                               ) {
	const auto midpoint_x2_fn = [] (const hit_seg &x) {
		return x.get_start_arrow().get_index() + x.get_stop_arrow().get_index();
	};
	return midpoint_x2_fn( arg_hit_seg_lhs ) < midpoint_x2_fn( arg_hit_seg_rhs );
}

/// \brief Generate a string describing the specified hit_seg
///
/// \relates hit_seg
string cath::rslv::to_simple_string(const hit_seg &arg_hit_seg ///< The hit_seg to describe
                                    ) {
	return
		  ::std::to_string( get_start_res_index( arg_hit_seg ) )
		+ "-"
		+ ::std::to_string( get_stop_res_index ( arg_hit_seg ) );
}

/// \brief Generate a string describing the specified hit_seg
///
/// \relates hit_seg
string cath::rslv::to_string(const hit_seg &arg_hit_seg ///< The hit_seg to describe
                             ) {
	return "hit_seg[" + to_simple_string( arg_hit_seg ) + "]";
}

/// \brief Insert a description of the specified hit_seg into the specified ostream
///
/// \relates hit_seg
ostream & cath::rslv::operator<<(ostream       &arg_ostream, ///< The ostream into which the description of the hit_seg should be inserted
                                 const hit_seg &arg_hit_seg  ///< The hit_seg to describe
                                 ) {
	arg_ostream << to_string( arg_hit_seg );
	return arg_ostream;
}

