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

#include "cath_id_score_category.hpp"

#include "cath/common/cpp20/starts_with_ends_with_contains.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

#include <regex>

using namespace ::cath::common;
using namespace ::cath::rslv;

using ::std::ostream;
using ::std::regex;
using ::std::string;
using ::std::string_view;

/// \brief Get the cath_id_score_category for the specified ID
///        (or return cath_id_score_category::NORMAL if `!prm_apply_cath_policies` )
///
/// \relates cath_id_score_category
cath_id_score_category cath::rslv::cath_score_category_of_id(const string_view &prm_id,                 ///< The ID to examine
                                                             const bool        &prm_apply_cath_policies ///< Whether to actually examine prm_id, rather than just returning cath_id_score_category::NORMAL
                                                             ) {
	if ( prm_apply_cath_policies ) {
		static const regex  dc_regex        { R"(^dc_\w{32}$)" };
		constexpr string_view dc_prefix_suffix = "dc_";


		if ( prm_id.length() == 35 && starts_with( prm_id, dc_prefix_suffix ) ) {
			if ( regex_search( cbegin( prm_id ), cend( prm_id ), dc_regex ) ) {
				return cath_id_score_category::DC_TYPE;
			}
		}
	}
	return cath_id_score_category::NORMAL;
}

/// \brief Generate a string describing the specified cath_id_score_category
///
/// \relates cath_id_score_category
string cath::rslv::to_string(const cath_id_score_category &prm_cath_id_score_category ///< The cath_id_score_category to describe
                             ) {
	switch ( prm_cath_id_score_category ) {
		case ( cath_id_score_category::NORMAL      ) : { return "cath_id_score_category::NORMAL"      ; }
		case ( cath_id_score_category::DC_TYPE     ) : { return "cath_id_score_category::DC_TYPE"     ; }
	}
	BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("Value of cath_id_score_category not recognised whilst converting to_string()"));
}

/// \brief Insert a description of the specified cath_id_score_category into the specified ostream
///
/// \relates cath_id_score_category
ostream & cath::rslv::operator<<(ostream                      &prm_os,                    ///< The ostream into which the description should be inserted
                                 const cath_id_score_category &prm_cath_id_score_category ///< The cath_id_score_category to describe
                                 ) {
	prm_os << to_string( prm_cath_id_score_category );
	return prm_os;
}