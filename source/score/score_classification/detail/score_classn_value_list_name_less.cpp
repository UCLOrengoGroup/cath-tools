/// \file
/// \brief The score_classn_value_list_name_less class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "score_classn_value_list_name_less.h"

#include "score/score_classification/score_classn_value_list.h"

using namespace cath::score::detail;
using namespace std;

/// \brief TODOCUMENT
bool score_classn_value_list_name_less::operator()(const score_classn_value_list &arg_score_classn_value_list_a, ///< The first  score_classn_value_list to compare
                                                   const score_classn_value_list &arg_score_classn_value_list_b  ///< The second score_classn_value_list to compare
                                                   ) const {
	return arg_score_classn_value_list_a.get_name() < arg_score_classn_value_list_b.get_name();
}

/// \brief TODOCUMENT
bool score_classn_value_list_name_less::operator()(const score_classn_value_list &arg_score_classn_value_list_a, ///< The first  score_classn_value_list to compare
                                                   const string                  &arg_name_b                     ///< The second score_classn_value_list to compare
                                                   ) const {
	return arg_score_classn_value_list_a.get_name() < arg_name_b;
}
