/// \file
/// \brief The protein_only_length_getter class definitions

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

#include "protein_only_length_getter.hpp"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/score/length_getter/length_of_longer_getter.hpp"
#include "cath/score/length_getter/length_of_shorter_getter.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::indeterminate;
using boost::numeric_cast;
using boost::ptr_vector;
using boost::tribool;

/// \brief TODOCUMENT
unique_ptr<length_getter> protein_only_length_getter::do_clone() const {
	return { protein_only_clone() };
}

/// \brief TODOCUMENT
tribool protein_only_length_getter::do_higher_is_better() const {
	return indeterminate;
}

/// \brief TODOCUMENT
size_t protein_only_length_getter::do_get_length(const alignment &/*prm_alignment*/, ///< TODOCUMENT
			                                     const protein   &prm_protein_a,     ///< TODOCUMENT
			                                     const protein   &prm_protein_b      ///< TODOCUMENT
			                                     ) const {
	return do_get_length(
		prm_protein_a,
		prm_protein_b
	);
}

/// \brief TODOCUMENT
string protein_only_length_getter::do_id_name() const {
	return do_get_choice_adjective() + "_protein_length";
}

/// \brief TODOCUMENT
str_bool_pair_vec protein_only_length_getter::do_short_name_suffixes() const {
	return {};
}

/// \brief TODOCUMENT
string protein_only_length_getter::do_long_name() const {
	return do_get_choice_adjective() + " length";
}

/// \brief TODOCUMENT
string protein_only_length_getter::do_description() const {
	return "The " + do_get_choice_adjective() + " protein length";
}

/// \brief TODOCUMENT
const string protein_only_length_getter::do_description_brackets_string() const {
	return " (over the " + do_get_choice_adjective() + " protein length)";
}

/// \brief TODOCUMENT
unique_ptr<protein_only_length_getter> protein_only_length_getter::protein_only_clone() const {
	return do_protein_only_clone();
}

/// \brief TODOCUMENT
size_t protein_only_length_getter::get_prot_only_length(const protein &prm_protein_a, ///< TODOCUMENT
                                                        const protein &prm_protein_b  ///< TODOCUMENT
                                                        ) const {
	return do_get_length( prm_protein_a, prm_protein_b );
}

/// \brief TODOCUMENT
string protein_only_length_getter::get_choice_adjective() const {
	return do_get_choice_adjective();
}

/// \brief TODOCUMENT
ptr_vector<protein_only_length_getter> cath::score::get_all_protein_only_length_getters() {
	ptr_vector<protein_only_length_getter> all_protein_only_length_getters;

	const ptr_vector<sym_protein_only_length_getter> all_sym_protein_only_length_getters = get_all_sym_protein_only_length_getters();
	for (const sym_protein_only_length_getter &getter : all_sym_protein_only_length_getters) {
		push_back( all_protein_only_length_getters, getter.protein_only_clone() );
	}
	return all_protein_only_length_getters;
}

/// \brief TODOCUMENT
///
/// \relates protein_only_length_getter
score_value cath::score::get_length_score(const protein_only_length_getter &prm_length_getter, ///< TODOCUMENT
                                          const protein                    &prm_protein_a,     ///< TODOCUMENT
                                          const protein                    &prm_protein_b      ///< TODOCUMENT
                                          ) {
	return numeric_cast<score_value>(
		prm_length_getter.get_prot_only_length( prm_protein_a, prm_protein_b )
	);
}
