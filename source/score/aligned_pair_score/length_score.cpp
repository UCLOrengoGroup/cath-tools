/// \file
/// \brief The length_score class definitions

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

#include "length_score.hpp"

#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/serialization/export.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

BOOST_CLASS_EXPORT(length_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> length_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that more aligned residues is generally better
tribool length_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the number of common residues defined by this alignment
///
/// This uses the score_common_coord_handler (and hence its policies)
score_value length_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
                                       const protein   &arg_protein_a, ///< The protein associated with the first  half of the alignment
                                       const protein   &arg_protein_b  ///< The protein associated with the second half of the alignment
                                       ) const {
	return get_length_score( *length_getter_ptr, arg_alignment, arg_protein_a, arg_protein_b );
}

/// \brief Concrete implementation that describes what this score means
string length_score::do_description() const {
	return length_getter_ptr->description();
}

/// \brief TODOCUMENT
string length_score::do_id_name() const {
	return length_getter_ptr->id_name();
}

/// \brief TODOCUMENT
str_bool_pair_vec length_score::do_short_name_suffixes() const {
	return length_getter_ptr->short_name_suffixes();
}


/// \brief Concrete implementation providing long name
string length_score::do_long_name() const {
	return length_getter_ptr->long_name();
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> length_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                           ) const {
//	cerr << "Should build a length_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool length_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                       ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief TODOCUMENT
length_score::length_score(const length_getter &arg_length_getter ///< TODOCUMENT
                           ) : length_getter_ptr( arg_length_getter.clone() ) {
}

/// \brief TODOCUMENT
const length_getter & length_score::get_length_getter() const {
	return *length_getter_ptr;
}

/// \brief Pass-through method to provide public access to the score_common_coord_handler's description_brackets_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string length_score::description_brackets_string() const {
	return length_getter_ptr->description_brackets_string();
}

/// \brief TODOCUMENT
///
/// \relates length_score
bool cath::score::operator<(const length_score &arg_length_score_a, ///< TODOCUMENT
                            const length_score &arg_length_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_length_score_a, arg_length_score_b );
	the_helper.register_comparison_field( &length_score::get_length_getter );
	return final_less_than_result( the_helper );
}

