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

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::irange;
using std::ostream;
using std::string;

/// \brief TODOCUMENT
///
/// \relates hit
string cath::rslv::get_segments_string(const hit &arg_hit ///< TODOCUMENT
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

/// \brief TODOCUMENT
///
/// \relates hit
string cath::rslv::to_string(const hit               &arg_hit,   ///< TODOCUMENT
                             const hit_output_format &arg_format ///< TODOCUMENT
                             ) {
	if ( arg_format == hit_output_format::JON ) {
		return arg_hit.get_label()
			+ " "
			+ ::std::to_string( arg_hit.get_score() )
			+ " "
			+ get_segments_string( arg_hit );
	}
	else {
		return "hit["
			+ get_segments_string( arg_hit )
			+ "; score: "
			+ ::std::to_string( arg_hit.get_score() )
			+ "; label: \""
			+ arg_hit.get_label()
			+ "\"]";
	}
}

/// \brief TODOCUMENT
///
/// \relates hit
ostream & cath::rslv::operator<<(ostream   &arg_ostream, ///< TODOCUMENT
                                 const hit &arg_hit      ///< TODOCUMENT
                                 ) {
	arg_ostream << to_string( arg_hit );
	return arg_ostream;
}


/// \brief TODOCUMENT
///
/// \relates hit
bool cath::rslv::operator==(const hit &arg_lhs, ///< TODOCUMENT
                            const hit &arg_rhs  ///< TODOCUMENT
                            ) {
	const bool similar = (
		arg_lhs.get_num_segments() == arg_rhs.get_num_segments()
		&&
		arg_lhs.get_score()        == arg_rhs.get_score()
		&&
		arg_lhs.get_label()        == arg_rhs.get_label()
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
