/// \file
/// \brief The link_dirn class definitions

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

#include "link_dirn.hpp"

#include <boost/algorithm/string/case_conv.hpp>

#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/program_options/validator.hpp"

#include <map>
#include <string>

using namespace ::cath::clust;
using namespace ::cath::clust::detail;
using namespace ::cath::common;

using ::boost::algorithm::to_upper;
using ::boost::any;
using ::std::istream;
using ::std::map;
using ::std::string;

/// \brief Getter for a map from name to link_dirn
///
/// \relates link_dirn
map<string, link_dirn> link_dirn_by_name::get() {
	using ::std::to_string;

	map<string, link_dirn> result;
	for (const auto &x : all_link_dirns) {
		result.emplace( to_string( x ), x );
	}
	return result;
}

/// \brief Generate a string describing the specified link_dirn (ie its name)
///
/// \relates link_dirn
string cath::clust::to_string(const link_dirn &prm_format_tag ///< The link_dirn to describe in a string
                              ) {
	switch ( prm_format_tag ) {
		case ( link_dirn::DISSIMILARITY ) : { return "DISTANCE" ; }
		case ( link_dirn::STRENGTH      ) : { return "STRENGTH" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("link_dirn value not recognised"));
}

/// \brief Extract into the specified link_dirn from the specified stream
///
/// \relates link_dirn
istream & cath::clust::operator>>(istream   &prm_is,       ///< The stream from which the link_dirn should be extracted
                                  link_dirn &prm_link_dirn ///< The link_dirn to populate from the specified stream
                                  ) {
	string input_string;
	prm_is >> input_string;
	to_upper( input_string );

	const auto all_link_dirns_by_name = link_dirn_by_name::get();
	prm_link_dirn = all_link_dirns_by_name.at( input_string );
	return prm_is;
}

/// \brief Generate a string containing a description of the specified link_dirn
///
/// \relates link_dirn
string cath::clust::description_of_link_dirn(const link_dirn &prm_format_tag ///< The link_dirn to describe
                                             ) {
	switch ( prm_format_tag ) {
		case ( link_dirn::DISSIMILARITY ) : { return "A higher value means the corresponding two entries are more distant"   ; }
		case ( link_dirn::STRENGTH      ) : { return "A higher value means the corresponding tow entries are more connected" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("link_dirn value not recognised"));
}

/// \brief Provide Boost program_options validation for link_dirn
///
/// \relates link_dirn
void cath::clust::validate(any           &prm_value,         ///< The value to populate
                           const str_vec &prm_value_strings, ///< The string values to validate
                           link_dirn *, int) {
	prm_value = lex_castable_validator<link_dirn>::perform_validate( prm_value, prm_value_strings );
}
