/// \file
/// \brief The length_getter class definitions

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

#include "length_getter.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>

#include "cath/alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/score/detail/score_name_helper.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"
#include "cath/score/length_getter/protein_only_length_getter.hpp"

#include <cassert>
#include <iostream> // ***** TEMPORARY *****

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::cath::score::detail;
using namespace ::std;

using ::boost::algorithm::any_of;
using ::boost::assign::ptr_push_back;
using ::boost::numeric_cast;
using ::boost::ptr_vector;
using ::boost::range::join;
using ::boost::tribool;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<length_getter> length_getter::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
tribool length_getter::higher_is_better() const {
	return do_higher_is_better();
}

/// \brief TODOCUMENT
size_t length_getter::get_length(const alignment &prm_alignment, ///< TODOCUMENT
                                 const protein   &prm_protein_a, ///< TODOCUMENT
                                 const protein   &prm_protein_b  ///< TODOCUMENT
                                 ) const {
	return do_get_length( prm_alignment, prm_protein_a, prm_protein_b );
}

/// \brief TODOCUMENT
length_getter_category length_getter::get_length_getter_category() const {
	return do_get_length_getter_category();
}

/// \brief TODOCUMENT
string length_getter::human_friendly_short_name() const {
	return score_name_helper::human_friendly_short_name( id_name(), short_name_suffixes() );
}

/// \brief TODOCUMENT
string length_getter::full_short_name() const {
	return score_name_helper::full_short_name( id_name(), short_name_suffixes() );
}

/// \brief TODOCUMENT
string length_getter::id_name() const {
	return do_id_name();
}

/// \brief TODOCUMENT
str_bool_pair_vec length_getter::short_name_suffixes() const {
	return do_short_name_suffixes();
}

/// \brief TODOCUMENT
string length_getter::long_name() const {
	return do_long_name();
}

/// \brief TODOCUMENT
string length_getter::description() const {
	return do_description();
}

/// \brief TODOCUMENT
string length_getter::description_brackets_string() const {
	return do_description_brackets_string();
}

/// \brief An NVI pass-through to the concrete class's do_less_than_with_same_dynamic_type(),
///        which defines the less-than operator when the argument's known to have the same dynamic type
bool length_getter::less_than_with_same_dynamic_type(const length_getter &prm_length_getter ///< TODOCUMENT
                                                     ) const {
	assert( typeid( *this ) == typeid( prm_length_getter ) );
	return do_less_than_with_same_dynamic_type( prm_length_getter );
}

/// \brief TODOCUMENT
str_bool_pair_vec cath::score::length_getter_as_short_name_suffixes(const length_getter &prm_length_getter,                    ///< TODOCUMENT
                                                                    const bool          &prm_include_id_name_in_human_friendly ///< TODOCUMENT
                                                                    ) {
	const auto short_name_suffixes            = prm_length_getter.short_name_suffixes();
	const bool has_suffixes_in_human_friendly = any_of( short_name_suffixes, [] (const str_bool_pair &x) { return x.second; } );
	const bool include_id_name                = has_suffixes_in_human_friendly || prm_include_id_name_in_human_friendly;
//	const str_bool_pair_vec id_name_suffixes               = { make_pair( prm_length_getter.id_name(), include_id_name ) };
	const auto id_name_suffixes               = { make_pair( prm_length_getter.id_name(), include_id_name ) };

	return copy_build<str_bool_pair_vec>( join(
		id_name_suffixes,
		short_name_suffixes
	) );
}

/// \brief TODOCUMENT
str_bool_pair_vec cath::score::length_getter_as_short_name_suffixes(const length_getter              &prm_length_getter,                          ///< TODOCUMENT
                                                                    const length_getter_category_vec &prm_category_to_exclude_from_human_friendly ///< TODOCUMENT
																	) {
	return length_getter_as_short_name_suffixes(
		prm_length_getter,
		! contains( prm_category_to_exclude_from_human_friendly, prm_length_getter.get_length_getter_category() )
	);
}

/// \brief TODOCUMENT
///
/// \relates length_getter
ptr_vector<length_getter> cath::score::get_all_length_getters() {
	ptr_vector<length_getter> all_length_getters;

	const ptr_vector<protein_only_length_getter> all_protein_only_length_getters = get_all_protein_only_length_getters();
	for (const protein_only_length_getter &getter : all_protein_only_length_getters) {
		push_back( all_length_getters, getter.clone() );
	}
	ptr_push_back< num_aligned_length_getter >( all_length_getters )( );

	return all_length_getters;
}

/// \brief TODOCUMENT
///
/// \relates length_getter
score_value cath::score::get_length_score(const length_getter &prm_length_getter, ///< TODOCUMENT
                                          const alignment     &prm_alignment,     ///< TODOCUMENT
                                          const protein       &prm_protein_a,     ///< TODOCUMENT
                                          const protein       &prm_protein_b      ///< TODOCUMENT
                                          ) {
	return numeric_cast<score_value>(
		prm_length_getter.get_length(
			prm_alignment,
			prm_protein_a,
			prm_protein_b
		)
	);
}
