/// \file
/// \brief The geometric_mean_length_getter class definitions

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

#include "geometric_mean_length_getter.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"

using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::numeric_cast;

/// \brief A standard do_clone method.
unique_ptr<sym_protein_only_length_getter> geometric_mean_length_getter::do_sym_protein_only_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
size_t geometric_mean_length_getter::do_get_length(const protein   &prm_protein_a,     ///< TODOCUMENT
                                                   const protein   &prm_protein_b      ///< TODOCUMENT
                                                   ) const {
	const size_t length_a       = prm_protein_a.get_length();
	const size_t length_b       = prm_protein_b.get_length();
	const double geom_mean_doub = sqrt( numeric_cast<double>( length_a * length_b ) );
	return numeric_cast<size_t>( round( geom_mean_doub ) );
}

/// \brief TODOCUMENT
length_getter_category geometric_mean_length_getter::do_get_length_getter_category() const {
	return length_getter_category::OTHER;
}

/// \brief TODOCUMENT
string geometric_mean_length_getter::do_get_choice_adjective() const {
	return "geometric_mean";
}

/// \brief TODOCUMENT
bool geometric_mean_length_getter::do_less_than_with_same_dynamic_type(const length_getter &/*prm_length_getter*/ ///< TODOCUMENT
                                                                       ) const {
	return false;
}
