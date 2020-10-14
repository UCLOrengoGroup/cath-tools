/// \file
/// \brief The value_list_scaling class definitions

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

#include "value_list_scaling.hpp"

#include <iostream>
#include <limits>
#include <string>

using namespace ::cath::score;
using namespace ::std;

constexpr double value_list_scaling::BAD_SCALED_VALUE;

/// \brief Simple to_string() overload for value_list_scaling
///
/// \relates value_list_scaling
string cath::score::to_string(const value_list_scaling &prm_scaling ///< The value_list_scaling to be output as a string
                              ) {
	return "value_list_scaling[ ( "
		+ ::std::to_string( prm_scaling.get_multiplier() )
		+ " * x ) + "
		+ ::std::to_string( prm_scaling.get_constant() )
		+ " ]";
}

/// \brief Simple insertion operator for value_list_scaling
///
/// \relates value_list_scaling
ostream & cath::score::operator<<(ostream                  &prm_os,     ///< The ostream to which the value_list_scaling should be output
                                  const value_list_scaling &prm_scaling ///< The value_list_scaling to output
                                  ) {
	prm_os << to_string( prm_scaling );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates value_list_scaling
void cath::score::scale_value(const value_list_scaling &prm_scaling, ///< TODOCUMENT
                              double                   &prm_value    ///< TODOCUMENT
                              ) {
	const auto &multiplier = prm_scaling.get_multiplier();

	// If the value is the worst_possible_value, then set it to BAD_SCALED_VALUE
	if ( ( multiplier < 0.0 && prm_value == numeric_limits<double>::max()    )
	     ||
	     ( multiplier > 0.0 && prm_value == numeric_limits<double>::lowest() ) ) {
		prm_value = value_list_scaling::BAD_SCALED_VALUE;
	}
	// Otherwise, perform the normal scaling
	else {
		prm_value *= multiplier;
		prm_value += prm_scaling.get_constant();
	}
}

/// \brief TODOCUMENT
///
/// \relates value_list_scaling
double cath::score::scale_value_copy(const value_list_scaling &prm_scaling, ///< TODOCUMENT
                                     double                    prm_value    ///< TODOCUMENT
                                     ) {
	scale_value( prm_scaling, prm_value );
	return prm_value;
}
