/// \file
/// \brief The hit class definitions

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

#include <boost/algorithm/string/join.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "resolve_hits/hit.h"

#include <string>

using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::irange;
using std::ostream;
using std::string;

// template <size_t T> class SD;
// SD< sizeof( hit         ) > sizeof_hit;
// SD< sizeof( res_arrow   ) > sizeof_res_arrow;
// SD< sizeof( resscr_t    ) > sizeof_resscr_t;
// SD< sizeof( hitidx_t    ) > sizeof_hitidx_t;
// SD< sizeof( hit_seg_vec ) > sizeof_hit_seg_vec;

/// \brief Generate a string describing the segements of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates hit
string cath::rslv::get_segments_string(const hit &arg_hit ///< The hit whose segments should be described
                                       ) {
	return join(
		boost::irange( 0_z, arg_hit.get_num_segments() )
			| boost::adaptors::transformed( [&] (const size_t &x) {
				return ::std::to_string( get_start_res_index_of_segment( arg_hit, x ) )
				     + "-"
				     + ::std::to_string( get_stop_res_index_of_segment ( arg_hit, x ) );
			} ),
		","
	);
}

/// \brief Generate a string describing the specified hit
///
/// \relates hit
string cath::rslv::to_string(const hit               &arg_hit,        ///< The hit to describe
                             const str_vec           &arg_hit_labels, ///< The list of labels correspding to the hit
                             const hit_output_format &arg_format,     ///< The format in which to generate the output
                             const string            &arg_prefix      ///< Any prefix that should come before the hit in hit_output_format::JON
                             ) {
	if ( arg_format != hit_output_format::JON && ! arg_prefix.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot specify prefix for any hit format other than JON"));
	}

	switch ( arg_format ) {
		case( hit_output_format::JON ) : {
			return arg_prefix
				+ ( arg_prefix.empty() ? ""s : " "s )
				+ arg_hit.get_label( arg_hit_labels )
				+ " "
				+ ::std::to_string( arg_hit.get_score() )
				+ " "
				+ get_segments_string( arg_hit );
		}
		case( hit_output_format::CLASS ) : {
			return "hit["
				+ get_segments_string( arg_hit )
				+ "; score: "
				+ ::std::to_string( arg_hit.get_score() )
				+ "; label: \""
				+ arg_hit.get_label( arg_hit_labels )
				+ "\"]";
		}
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of hit_output_format not recognised in to_string() for hit"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Return whether the two specified hits are identical
///
/// \relates hit
bool cath::rslv::operator==(const hit &arg_lhs, ///< The first  hit to compare
                            const hit &arg_rhs  ///< The second hit to compare
                            ) {
	const bool similar = (
		arg_lhs.get_num_segments() == arg_rhs.get_num_segments()
		&&
		arg_lhs.get_score()        == arg_rhs.get_score()
		&&
		arg_lhs.get_label_idx()    == arg_rhs.get_label_idx()
	);
	if ( ! similar ) {
		return false;
	}
	for (const size_t seg_ctr : irange( 0_z, arg_lhs.get_num_segments() ) ) {
		const bool segs_differ = (
			arg_lhs.get_start_arrow_of_segment( seg_ctr ) != arg_rhs.get_start_arrow_of_segment( seg_ctr )
			||
			arg_lhs.get_stop_arrow_of_segment ( seg_ctr ) != arg_rhs.get_stop_arrow_of_segment ( seg_ctr )
		);
		if ( segs_differ ) {
			return false;
		}
	}
	return true;
}
