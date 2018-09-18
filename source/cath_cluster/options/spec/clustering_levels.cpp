/// \file
/// \brief The clustering_levels class definitions

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

#include "clustering_levels.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/program_options/validator.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::algorithm::is_any_of;
using ::boost::any;
using ::boost::lexical_cast;
using ::std::istream;
using ::std::string;

/// \brief Extract into the specified clustering_levels from the specified stream
///
/// \relates clustering_levels
std::istream & cath::clust::operator>>(istream           &prm_is,               ///< The stream from which the clustering_levels should be extracted
                                       clustering_levels &prm_clustering_levels ///< The clustering_levels to populate from the specified stream
                                       ) {
	string input_string;
	prm_is >> input_string;

	prm_clustering_levels.levels = transform_build<strength_vec>(
		split_build<str_vec>( input_string, is_any_of( "," ) ),
		[] (const string &x) { return lexical_cast<strength>( x ); }
	);

	return prm_is;
}

/// \brief Provide Boost program_options validation for clustering_levels
///
/// \relates clustering_levels
void cath::clust::validate(any           &prm_value,         ///< The value to populate
                           const str_vec &prm_value_strings, ///< The string values to validate
                           clustering_levels *, int) {
	prm_value = lex_castable_validator<clustering_levels>::perform_validate( prm_value, prm_value_strings );
}
