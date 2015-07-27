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

#include "value_list_scaling.h"

using namespace cath::score;
// using namespace std;

/// \brief TODOCUMENT
value_list_scaling::value_list_scaling(const double &arg_multiplier, ///< TODOCUMENT
                                       const double &arg_constant    ///< TODOCUMENT
                                       ) : multiplier ( arg_multiplier ),
                                           constant   ( arg_constant   ) {
}

/// \brief TODOCUMENT
const double & value_list_scaling::get_multiplier() const {
	return multiplier;
}

/// \brief TODOCUMENT
const double & value_list_scaling::get_constant() const {
	return constant;
}

/// \brief TODOCUMENT
void cath::score::scale_value(const value_list_scaling &arg_scaling, ///< TODOCUMENT
                              double                   &arg_value    ///< TODOCUMENT
                              ) {
	arg_value *= arg_scaling.get_multiplier();
	arg_value += arg_scaling.get_constant();
}

/// \brief TODOCUMENT
double cath::score::scale_value_copy(const value_list_scaling &arg_scaling, ///< TODOCUMENT
                                     double                    arg_value    ///< TODOCUMENT
                                     ) {
	scale_value( arg_scaling, arg_value );
	return arg_value;
}
