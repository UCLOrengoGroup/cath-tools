/// \file
/// \brief The calc_hit class definitions

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

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "resolve_hits/calc_hit.hpp"

#include <string>

using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::algorithm::join;
using boost::irange;
using std::ostream;
using std::string;

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates calc_hit
string cath::rslv::get_segments_string(const calc_hit &arg_hit ///< The calc_hit whose segments should be described
                                       ) {
	return join(
		boost::irange( 0_z, get_num_segments( arg_hit ) )
			| boost::adaptors::transformed( [&] (const size_t &x) {
				return ::std::to_string( get_start_res_index_of_segment( arg_hit, x ) )
				     + "-"
				     + ::std::to_string( get_stop_res_index_of_segment ( arg_hit, x ) );
			} ),
		","
	);
}

/// \brief Generate a string describing the specified calc_hit
///
/// \relates calc_hit
string cath::rslv::to_string(const calc_hit &arg_hit ///< The calc_hit to describe
                             ) {
	return "calc_hit[" + get_segments_string( arg_hit ) + "]";
}

/// \brief Insert a description of the specified calc_hit into the specified ostream
///
/// \relates calc_hit
ostream & cath::rslv::operator<<(ostream        &arg_os,      ///< The ostream into which the description should be inserted
                                 const calc_hit &arg_calc_hit ///< The calc_hit to describe
                                 ) {
	arg_os << to_string( arg_calc_hit );
	return arg_os;
}

/// \brief Return whether the two specified hits are identical
///
/// \relates calc_hit
bool cath::rslv::operator==(const calc_hit &arg_lhs, ///< The first  calc_hit to compare
                            const calc_hit &arg_rhs  ///< The second calc_hit to compare
                            ) {
	return (
		( get_num_segments( arg_lhs ) == get_num_segments( arg_rhs ) )
		&&
		( arg_lhs.get_score()         == arg_rhs.get_score()         )
		&&
		( arg_lhs.get_label_idx()     == arg_rhs.get_label_idx()     )
		&&
		all_of(
			irange( 0_z, get_num_segments( arg_lhs ) ),
			[&] (const size_t &seg_ctr) {
				return (
					get_start_arrow_of_segment( arg_lhs, seg_ctr ) == get_start_arrow_of_segment( arg_rhs, seg_ctr )
					&&
					get_stop_arrow_of_segment ( arg_lhs, seg_ctr ) == get_stop_arrow_of_segment ( arg_rhs, seg_ctr )
				);
			}
		)
	);
}
