/// \file
/// \brief The alignment_split_list class definitions

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

#include "alignment_split_list.hpp"

#include "alignment/alignment.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"

#include <iostream>

using namespace cath;
using namespace cath::align::detail;
using namespace cath::common;

/// \brief TODOCUMENT
void alignment_split_list::insert(const alignment_split &arg_alignment_split ///< TODOCUMENT
                                  ) {
	splits.insert( arg_alignment_split );
}

/// \brief TODOCUMENT
alignment_split_list::iterator alignment_split_list::begin() {
	return std::begin( splits );
}

/// \brief TODOCUMENT
alignment_split_list::iterator alignment_split_list::end() {
	return std::end( splits );
}

/// \brief TODOCUMENT
alignment_split_list::const_iterator alignment_split_list::begin() const {
	return common::cbegin( splits );
}

/// \brief TODOCUMENT
alignment_split_list::const_iterator alignment_split_list::end() const {
	return common::cend( splits );
}

/// \brief TODOCUMENT
alignment_split_list cath::align::detail::make_list_of_alignment_split(const alignment &arg_alignment, ///< TODOCUMENT
                                                                       const size_vec  &arg_indices    ///< TODOCUMENT
                                                                       ) {
	alignment_split_list new_alignment_splits;
	new_alignment_splits.insert( make_alignment_split(
		arg_indices,
		arg_alignment.num_entries()
	) );
	return new_alignment_splits;
}

/// \brief TODOCUMENT
alignment_split_list cath::align::detail::get_all_single_alignment_splits(const alignment &arg_alignment ///< TODOCUMENT
                                                                          ) {
	const size_t num_entries = arg_alignment.num_entries();
	alignment_split_list new_alignment_splits;
	for (const size_t &aln_entry : indices( num_entries ) ) {
		const alignment_split single_split = make_single_alignment_split( aln_entry, num_entries );
		if ( is_valid_split( single_split) ) {
			new_alignment_splits.insert( get_least_version( single_split ) );
		}
	}
	return new_alignment_splits;
}

/// \brief TODOCUMENT
alignment_split_list cath::align::detail::get_preexisting_alignment_splits(const alignment &arg_alignment ///< TODOCUMENT
                                                                           ) {
	const size_t num_entries = arg_alignment.num_entries();
	const size_t aln_length  = arg_alignment.length();
	alignment_split_list new_alignment_splits;
	for (const size_t &aln_index : indices( aln_length ) ) {
		const size_vec present_positions = entries_present_at_index( arg_alignment, aln_index );
		const alignment_split multi_split = make_alignment_split( present_positions, num_entries );
		if ( is_valid_split( multi_split) ) {
			new_alignment_splits.insert( get_least_version( multi_split ) );
		}
	}
	return new_alignment_splits;
}

/// \brief TODOCUMENT
alignment_split_list cath::align::detail::get_standard_alignment_splits(const alignment &arg_alignment ///< TODOCUMENT
                                                                        ) {
	return add_alignment_splits_copy(
		get_all_single_alignment_splits ( arg_alignment ),
		get_preexisting_alignment_splits( arg_alignment )
	);
}

/// \brief TODOCUMENT
void cath::align::detail::add_alignment_splits(alignment_split_list       &arg_alignment_split,      ///< TODOCUMENT
                                               const alignment_split_list &arg_other_alignment_split ///< TODOCUMENT
                                               ) {
	for (const alignment_split &the_split : arg_other_alignment_split) {
		arg_alignment_split.insert( the_split );
	}
}

/// \brief TODOCUMENT
alignment_split_list cath::align::detail::add_alignment_splits_copy(alignment_split_list        arg_alignment_split,      ///< TODOCUMENT
                                                                    const alignment_split_list &arg_other_alignment_split ///< TODOCUMENT
                                                                    ) {
	add_alignment_splits( arg_alignment_split, arg_other_alignment_split );
	return arg_alignment_split;
}
