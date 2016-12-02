/// \file
/// \brief The cath_id_score_category definitions

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

#include "hmmer_hmmsearch_domtblout.hpp"

#include <boost/utility/string_ref.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "exception/out_of_range_exception.hpp"

#include <regex>

using namespace cath::rslv;

using boost::string_ref;
using std::ostream;
using std::regex;
using std::string;

/// \brief Get the cath_id_score_category for the specified ID
///        (or return cath_id_score_category::NORMAL if `!arg_apply_cath_policies` )
///
/// \relates cath_id_score_category
cath_id_score_category cath::rslv::cath_score_category_of_id(const string_ref &arg_id,                 ///< The ID to examine
                                                             const bool       &arg_apply_cath_policies ///< Whether to actually examine arg_id, rather than just returning cath_id_score_category::NORMAL
                                                             ) {
	if ( ! arg_apply_cath_policies ) {
		return cath_id_score_category::NORMAL;
	}

	static const regex  dc_regex        { R"(^dc_\w{32}$)" };
	static const regex  round_regex     { R"(_round_\d+$)" };
	static const string round_one_suffix{ "_round_1" };
	static const string dc_prefix_suffix{ "dc_" };

	if ( arg_id.ends_with( round_one_suffix ) ) {
		return cath_id_score_category::NORMAL;
	}
	if ( arg_id.length() == 35 && arg_id.starts_with( dc_prefix_suffix ) ) {
		if (   regex_search( common::cbegin( arg_id ), common::cend( arg_id ), dc_regex        ) ) {
			return cath_id_score_category::DC_TYPE;
		}
	}
	if ( regex_search( common::cbegin( arg_id ), common::cend( arg_id ), round_regex     ) ) {
		return cath_id_score_category::LATER_ROUND;
	}
	return cath_id_score_category::NORMAL;
}

/// \brief Generate a string describing the specified cath_id_score_category
///
/// \relates cath_id_score_category
string cath::rslv::to_string(const cath_id_score_category &arg_cath_id_score_category ///< The cath_id_score_category to describe
                             ) {
	switch ( arg_cath_id_score_category ) {
		case ( cath_id_score_category::NORMAL      ) : { return "cath_id_score_category::NORMAL"      ; }
		case ( cath_id_score_category::DC_TYPE     ) : { return "cath_id_score_category::DC_TYPE"     ; }
		case ( cath_id_score_category::LATER_ROUND ) : { return "cath_id_score_category::LATER_ROUND" ; }
		default : {
			BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("Value of cath_id_score_category not recognised whilst converting to_string()"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Insert a description of the specified cath_id_score_category into the specified ostream
///
/// \relates cath_id_score_category
ostream & cath::rslv::operator<<(ostream                      &arg_os,                    ///< The ostream into which the description should be inserted
                                 const cath_id_score_category &arg_cath_id_score_category ///< The cath_id_score_category to describe
                                 ) {
	arg_os << to_string( arg_cath_id_score_category );
	return arg_os;
}