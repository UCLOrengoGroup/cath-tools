/// \file
/// \brief The make_clone definitions

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include "length_getter_make_clone.hpp"

#include "cath/score/length_getter/protein_only_length_getter.hpp"
#include "cath/score/length_getter/sym_protein_only_length_getter.hpp"

using namespace cath::score;
using namespace std;

/// \brief Specialisation for getting a protein_only_length_getter clone from a protein_only_length_getter
///        object.
///
/// clone() is already used by length_getter::clone()
template <>
unique_ptr<protein_only_length_getter> cath::common::detail::make_clone(const protein_only_length_getter &prm_value ///< The value to be cloned
                                                                        ) {
	return prm_value.protein_only_clone();
}

/// \brief Specialisation for getting a sym_protein_only_length_getter clone from a sym_protein_only_length_getter
///        object.
///
/// clone() is already used by length_getter::clone()
template <>
unique_ptr<sym_protein_only_length_getter> cath::common::detail::make_clone(const sym_protein_only_length_getter &prm_value ///< The value to be cloned
                                                                            ) {
	return prm_value.sym_protein_only_clone();
}
