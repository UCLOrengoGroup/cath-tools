/// \file
/// \brief The seq_seg class definitions

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

#include "seq_seg.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/exception/invalid_argument_exception.hpp"

#include <iostream>

using namespace cath::common;
using namespace cath::seq;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::range::adjacent_find;
using boost::range::sort;
using std::ostream;
using std::pair;
using std::string;

/// \brief Get a vector of the entries that are present in the specified vector of optional seq_segs
///
/// \relates seq_seg
seq_seg_vec cath::seq::get_present_segments(const seq_seg_opt_vec &arg_segments ///< The vector of optional seq_segs to query
                                            ) {
	return transform_build<seq_seg_vec>(
		arg_segments
			| filtered( [&] (const seq_seg_opt &x) {
				return static_cast<bool>( x );
			} ),
		[&] (const seq_seg_opt &x) { return *x; }
	);
}

/// \brief Get a simple string describing the specified segments (eg "1-3,5-9")
///
/// \relates seq_seg
string cath::seq::get_segments_string(const seq_seg_vec &arg_segments ///< The segments to describe
                                      ) {
	return join(
		arg_segments
			| transformed( [&] (const seq_seg &x) {
				return to_simple_string( x );
			} ),
		","
	);
}

/// \brief Make the (sorted) fragments between the specified segments
///
/// \relates seq_seg
seq_seg_vec cath::seq::make_fragments_of_segments(seq_seg_vec arg_segments ///< The segments from which the fragments should be made
                                                  ) {
	start_sort_seq_segs( arg_segments );
	return make_fragments_of_start_sorted_segments( arg_segments );
}

/// \brief Return whether the specified segments are sorted by their starts and non-overlapping
///
/// \relates seq_seg
bool cath::seq::segments_are_start_sorted_and_non_overlapping(const seq_seg_vec &arg_segments ///< The segments to query
                                                              ) {
	const auto seg_pair_are_invalid = [] (const seq_seg &x, const seq_seg &y) {
		return ( y.get_start_arrow() < x.get_stop_arrow() );
	};
	return ( adjacent_find( arg_segments, seg_pair_are_invalid ) == common::cend( arg_segments ) );
}

/// \brief Make the fragments between the specified segments
///
/// \pre The segments must be pre-sorted by their starts
///
/// \relates seq_seg
seq_seg_vec cath::seq::make_fragments_of_start_sorted_segments(const seq_seg_vec &arg_segments ///< The segments from which the fragments must be made, pre-sorted by their starts
                                                               ) {
	if ( arg_segments.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot make fragments from empty vector of segments"));
	}

#ifndef NDEBUG
	if ( ! segments_are_start_sorted_and_non_overlapping( arg_segments ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot make fragments from segments with misordered/overlapping invalid segments (start > stop)"));
	}
#endif

	return transform_build<seq_seg_vec>(
		arg_segments | adjacented,
		[] (const pair<seq_seg, seq_seg> &x) {
			return seq_seg{
				x.first .get_stop_arrow(),
				x.second.get_start_arrow()
			};
		}
	);
}

/// \brief Return whether the midpoint of the first specified seq_seg is less than that of the second
///
/// This cuts out the halving operation, which isn't actually needed here
///
/// \relates seq_seg
bool cath::seq::midpoint_less(const seq_seg &arg_seq_seg_lhs, ///< The first  seq_seg to compare
                              const seq_seg &arg_seq_seg_rhs  ///< The second seq_seg to compare
                              ) {
	const auto midpoint_x2_fn = [] (const seq_seg &x) {
		return x.get_start_arrow().get_index() + x.get_stop_arrow().get_index();
	};
	return midpoint_x2_fn( arg_seq_seg_lhs ) < midpoint_x2_fn( arg_seq_seg_rhs );
}

/// \brief Generate a string describing the specified seq_seg
///
/// \relates seq_seg
string cath::seq::to_simple_string(const seq_seg &arg_seq_seg ///< The seq_seg to describe
                                   ) {
	return
		  ::std::to_string( get_start_res_index( arg_seq_seg ) )
		+ "-"
		+ ::std::to_string( get_stop_res_index ( arg_seq_seg ) );
}

/// \brief Generate a string describing the specified seq_seg_opt
///
/// \relates seq_seg_opt
string cath::seq::to_simple_string(const seq_seg_opt &arg_seq_seg_opt ///< The seq_seg_opt to describe
                                   ) {
	return arg_seq_seg_opt
		? to_simple_string( *arg_seq_seg_opt )
		: "absent";
}

/// \brief Generate a string describing a segment between the two specified positions
///
/// \relates seq_seg
string cath::seq::to_simple_seg_string(const seq_arrow &arg_start_arrow, ///< The start position
                                       const seq_arrow &arg_stop_arrow   ///< The stop  position
                                       ) {
	return to_simple_string( seq_seg{ arg_start_arrow, arg_stop_arrow } );
}

/// \brief Generate a string describing the specified seq_seg
///
/// \relates seq_seg
string cath::seq::to_string(const seq_seg &arg_seq_seg ///< The seq_seg to describe
                            ) {
	return "seq_seg[" + to_simple_string( arg_seq_seg ) + "]";
}

/// \brief Insert a description of the specified seq_seg into the specified ostream
///
/// \relates seq_seg
ostream & cath::seq::operator<<(ostream       &arg_ostream, ///< The ostream into which the description of the seq_seg should be inserted
                                const seq_seg &arg_seq_seg  ///< The seq_seg to describe
                                ) {
	arg_ostream << to_string( arg_seq_seg );
	return arg_ostream;
}

