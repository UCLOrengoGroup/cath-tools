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

#include "full_hit.hpp"

#include <boost/format.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::common;
using namespace ::cath::rslv;

using ::boost::format;
using ::std::string;

/// \brief Generate a formatted string for the specified score of the specified type with the specified number of significant figures (roughly)
std::string cath::rslv::get_score_string(const double         &prm_score,      ///< The score to represent in a string
                                         const hit_score_type &prm_score_type, ///< The type of score to represent
                                         const size_t         &prm_num_figures ///< The number of significant figures (roughly)
                                         ) {
	switch ( prm_score_type ) {
		case ( hit_score_type::FULL_EVALUE ) : { return ( format( "%." + ::std::to_string( prm_num_figures ) + "e" ) % prm_score ).str(); }
		case ( hit_score_type::BITSCORE    ) : { return ( format( "%." + ::std::to_string( prm_num_figures ) + "g" ) % prm_score ).str(); }
		case ( hit_score_type::CRH_SCORE   ) : { return ( format( "%." + ::std::to_string( prm_num_figures ) + "g" ) % prm_score ).str(); }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of hit_score_type not recognised whilst getting score string"));
}

/// \brief Generate a formatted string for the specified full_hit's native score with the specified number of significant figures (roughly)
///
/// \relates full_hit
string cath::rslv::get_score_string(const full_hit &prm_full_hit,   ///< The full_hit containing the score to represent in a string
                                    const size_t   &prm_num_figures ///< The number of significant figures (roughly)
                                    ) {
	return get_score_string(
		prm_full_hit.get_score(),
		prm_full_hit.get_score_type(),
		prm_num_figures
	);
}

/// \brief Getter for the name to use when outputting the prefix (ie query-id)
string full_hit::get_prefix_name() {
	return "query-id";
}

/// \brief Getter for the name to use when outputting the label (ie match-id)
string full_hit::get_label_name() {
	return "match-id";
}

/// \brief Getter for the name to use when outputting the resolved boundaries
string full_hit::get_resolved_name() {
	return "resolved";
}

/// \brief Getter for the name to use when outputting the score
string full_hit::get_score_name() {
	return "score";
}

/// \brief Getter for the name to use when outputting the score type
string full_hit::get_score_type_name() {
	return "score-type";
}

/// \brief Getter for the name to use when outputting the boundaries
string full_hit::get_segments_name() {
	return "boundaries";
}

/// \brief Getter for the name to use when outputting the trimmed boundaries
string full_hit::get_trimmed_name() {
	return "trimmed";
}
