/// \file
/// \brief The full_hit class definitions

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

#include "full_hit.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>

using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::algorithm::join;
using boost::format;
using boost::optional;
using std::ostream;
using std::string;

/// \brief Generate a formatted string for the specified score of the specified type with the specified number of significant figures (roughly)
std::string cath::rslv::get_score_string(const double         &arg_score,      ///< The score to represent in a string
                                         const hit_score_type &arg_score_type, ///< The type of score to represent
                                         const size_t         &arg_num_figures ///< The number of significant figures (roughly)
                                         ) {
	switch ( arg_score_type ) {
		case ( hit_score_type::FULL_EVALUE ) : { return ( format( "%." + ::std::to_string( arg_num_figures ) + "e" ) % arg_score ).str(); }
		case ( hit_score_type::BITSCORE    ) : { return ( format( "%." + ::std::to_string( arg_num_figures ) + "g" ) % arg_score ).str(); }
		case ( hit_score_type::CRH_SCORE   ) : { return ( format( "%." + ::std::to_string( arg_num_figures ) + "g" ) % arg_score ).str(); }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of hit_score_type not recognised whilst getting score string"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Generate a formatted string for the specified's full_hit's native score with the specified number of significant figures (roughly)
///
/// \relates full_hit
string cath::rslv::get_score_string(const full_hit &arg_full_hit,   ///< The full_hit containing the score to represent in a string
                                    const size_t   &arg_num_figures ///< The number of significant figures (roughly)
                                    ) {
	return get_score_string(
		arg_full_hit.get_score(),
		arg_full_hit.get_score_type(),
		arg_num_figures
	);
}

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
string cath::rslv::get_segments_string(const full_hit            &arg_full_hit,     ///< The full_hit whose segments should be described
                                       const optional<trim_spec> &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                       ) {
	return get_segments_string( arg_full_hit.get_segments(), arg_trim_spec_opt );
}

/// \brief Generate a string describing the specified full_hit
///
/// \relates full_hit
string cath::rslv::to_string(const full_hit            &arg_full_hit,     ///< The full_hit to describe
                             const hit_output_format   &arg_format,       ///< The format in which to generate the output
                             const string              &arg_prefix,       ///< A prefix string, typically used to put the query_id at the front. (Any non-empty string will have a space appended.)
                             const optional<trim_spec> &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                             ) {
	if ( arg_format != hit_output_format::JON && ! arg_prefix.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot specify prefix for any full_hit format other than JON"));
	}

	switch ( arg_format ) {
		case( hit_output_format::JON ) : {
			return arg_prefix
				+ ( arg_prefix.empty() ? ""s : " "s )
				+ arg_full_hit.get_label()
				+ " "
				+ get_score_string( arg_full_hit, 6 )
				+ " "
				+ get_segments_string( arg_full_hit, arg_trim_spec_opt );
		}
		case( hit_output_format::CLASS ) : {
			return "full_hit["
				+ get_segments_string( arg_full_hit, arg_trim_spec_opt )
				+ "; score: "
				+ get_score_string( arg_full_hit, 6 )
				+ "; label: \""
				+ arg_full_hit.get_label()
				+ "\"]";
		}
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of full_hit_output_format not recognised in to_string() for full_hit"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}


