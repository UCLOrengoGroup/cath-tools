/// \file
/// \brief The num_aligned_length_getter class definitions

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

#include "num_aligned_length_getter.h"

#include <boost/logic/tribool.hpp>

#include "alignment/common_atom_selection_policy/common_atom_selection_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"
#include "common/clone/make_uptr_clone.h"
#include "structure/geometry/coord.h"

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<length_getter> num_aligned_length_getter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
tribool num_aligned_length_getter::do_higher_is_better() const {
	return true;
}

/// \brief TODOCUMENT
size_t num_aligned_length_getter::do_get_length(const alignment &arg_alignment, ///< TODOCUMENT
                                                const protein   &arg_protein_a, ///< TODOCUMENT
                                                const protein   &arg_protein_b  ///< TODOCUMENT
                                                ) const {
	const pair<coord_list, coord_list> common_coords = common_coord_handler.get_common_coords(
		arg_alignment,
		arg_protein_a,
		arg_protein_b
	);
	return common_coords.first.size();
}

/// \brief TODOCUMENT
length_getter_category num_aligned_length_getter::do_get_length_getter_category() const {
	return length_getter_category::OTHER;
}

/// \brief TODOCUMENT
string num_aligned_length_getter::do_id_name() const {
	return "num_aligned_residues";
}

/// \brief TODOCUMENT
str_bool_pair_vec num_aligned_length_getter::do_short_name_suffixes() const {
	return common_coord_handler.short_name_suffixes();
}

/// \brief TODOCUMENT
string num_aligned_length_getter::do_long_name() const {
	return "Number of aligned residues" + common_coord_handler.long_suffix_string();
}

/// \brief TODOCUMENT
string num_aligned_length_getter::do_description() const {
	return { "The number of pairs of residues that have been aligned together"
	            + common_coord_handler.description_brackets_string()
	            + ". Note that the fact that a pair of residues has been aligned"
	            + " doesn't necessarily mean that they have been deemed highly similar/equivalent." };
}

///// \brief TODOCUMENT
//string num_aligned_length_getter::do_short_suffix_string() const {
//	return common_coord_handler.short_suffix_string();
//}

///// \brief TODOCUMENT
//string num_aligned_length_getter::do_long_suffix_string() const {
//	return common_coord_handler.long_suffix_string();
//}

/// \brief TODOCUMENT
const string num_aligned_length_getter::do_description_brackets_string() const {
	return common_coord_handler.long_suffix_string();
}

/// \brief TODOCUMENT
bool num_aligned_length_getter::do_less_than_with_same_dynamic_type(const length_getter &arg_length_getter ///< TODOCUMENT
                                                                    ) const {
	const num_aligned_length_getter &casted_length_getter = dynamic_cast<const num_aligned_length_getter &>( arg_length_getter );
	return ( *this < casted_length_getter );
}

/// \brief TODOCUMENT
const score_common_coord_handler & num_aligned_length_getter::get_common_coord_handler() const {
	return common_coord_handler;
}

/// \brief Ctor for num_aligned_length_getter
num_aligned_length_getter::num_aligned_length_getter(const common_residue_selection_policy &arg_comm_res_seln_pol ///< TODOCUMENT
                                                     ) : common_coord_handler(
                                                         	arg_comm_res_seln_pol,
                                                         	*make_default_common_atom_selection_policy()
                                                         ) {
}

/// \brief TODOCUMENT
///
/// \relates num_aligned_length_getter
bool cath::score::operator<(const num_aligned_length_getter &arg_num_aligned_length_getter_a, ///< TODOCUMENT
                            const num_aligned_length_getter &arg_num_aligned_length_getter_b  ///< TODOCUMENT
                            ) {
	return (
		arg_num_aligned_length_getter_a.get_common_coord_handler()
		<
		arg_num_aligned_length_getter_b.get_common_coord_handler()
	);

}
