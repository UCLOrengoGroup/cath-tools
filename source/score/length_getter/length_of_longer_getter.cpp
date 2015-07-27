/// \file
/// \brief The length_of_longer_getter class definitions

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

#include "length_of_longer_getter.h"

#include "common/clone/make_uptr_clone.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"

using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<sym_protein_only_length_getter> length_of_longer_getter::do_sym_protein_only_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
size_t length_of_longer_getter::do_get_length(const protein   &arg_protein_a,     ///< TODOCUMENT
                                              const protein   &arg_protein_b      ///< TODOCUMENT
                                              ) const {
	const size_t length_a = arg_protein_a.get_length();
	const size_t length_b = arg_protein_b.get_length();
	return length_a >= length_b ? length_a : length_b;
}

/// \brief TODOCUMENT
length_getter_category length_of_longer_getter::do_get_length_getter_category() const {
	return length_getter_category::LONGER;
}

/// \brief TODOCUMENT
string length_of_longer_getter::do_get_choice_adjective() const {
	return "longer";
}

/// \brief TODOCUMENT
bool length_of_longer_getter::do_less_than_with_same_dynamic_type(const length_getter &/*arg_length_getter*/ ///< TODOCUMENT
                                                                   ) const {
	return false;
}
