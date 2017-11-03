/// \file
/// \brief The alignment_split class definitions

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

#include "alignment_split.hpp"

#include "common/algorithm/contains.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/size_t_literal.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"

using namespace cath;
using namespace cath::align::detail;
using namespace cath::common;

/// \brief Ctor for alignment_split
alignment_split::alignment_split(const size_t &arg_num_entries ///< TODOCUMENT
                                 ) : num_entries( arg_num_entries ) {
}

/// \brief TODOCUMENT
size_t alignment_split::get_num_entries() const {
	return num_entries;
}

/// \brief TODOCUMENT
void alignment_split::add_first_half_entry(const size_t &arg_first_half_entry ///< TODOCUMENT
                                           ) {
	if ( arg_first_half_entry >= num_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add first_half_entry that is out of range given number of entries"));
	}
	first_half_entries.insert( arg_first_half_entry );
}

/// \brief TODOCUMENT
size_t alignment_split::get_num_first_half_entries() const {
	return first_half_entries.size();
}

/// \brief TODOCUMENT
const size_set & alignment_split::get_first_half_entries() const {
	return first_half_entries;
}

/// \brief TODOCUMENT
alignment_split::const_iterator alignment_split::begin() const {
	return common::cbegin( first_half_entries );
}

/// \brief TODOCUMENT
alignment_split::const_iterator alignment_split::end() const {
	return common::cend( first_half_entries );
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
bool cath::align::detail::operator<(const alignment_split &arg_split_a, ///< TODOCUMENT
                                    const alignment_split &arg_split_b  ///< TODOCUMENT
                                    ) {
	if ( arg_split_a.get_num_entries() != arg_split_b.get_num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot compare alignment_splits for different numbers of entries"));
	}
	const size_t num_first_half_entries_a = arg_split_a.get_num_first_half_entries();
	const size_t num_first_half_entries_b = arg_split_b.get_num_first_half_entries();
	if ( num_first_half_entries_a != num_first_half_entries_b ) {
		return ( num_first_half_entries_a < num_first_half_entries_b );
	}
	return ( arg_split_a.get_first_half_entries() < arg_split_b.get_first_half_entries() );
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
size_set cath::align::detail::entries_of_alignment_split_half(const alignment_split      &arg_split,     ///< TODOCUMENT
                                                              const alignment_split_half &arg_split_half ///< TODOCUMENT
                                                              ) {
	switch ( arg_split_half ) {
		case ( alignment_split_half::FIRST  ) : {
			return arg_split.get_first_half_entries();
		}
		case ( alignment_split_half::SECOND ) : {
			return make_opposite_version( arg_split ).get_first_half_entries();
		}
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("Split half is neither first half nor second half"));
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
alignment_split cath::align::detail::make_opposite_version(const alignment_split &arg_alignment_split ///< TODOCUMENT
                                                           ) {
	const size_t    num_entries        = arg_alignment_split.get_num_entries();
	const size_set &first_half_entries = arg_alignment_split.get_first_half_entries();
	alignment_split new_split( num_entries );
	for (const size_t &new_first_half_ctr : indices( num_entries ) ) {
		if ( ! contains( first_half_entries, new_first_half_ctr ) ) {
			new_split.add_first_half_entry( new_first_half_ctr );
		}
	}
	return new_split;
}


/// \brief TODOCUMENT
///
/// \relates alignment_split
alignment_split cath::align::detail::get_least_version(const alignment_split &arg_alignment_split ///< TODOCUMENT
                                                       ) {
	if ( ! is_valid_split(arg_alignment_split) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get least version for invalid arg_alignment_split"));
	}

	const size_t num_entries        = arg_alignment_split.get_num_entries();
	const size_t num_in_first_half  = arg_alignment_split.get_num_first_half_entries();
	const size_t num_in_second_half = num_entries - num_in_first_half;

	// If the second half is smaller then return that
	if ( num_in_second_half < num_in_first_half ) {
		return make_opposite_version( arg_alignment_split );
	}

	// If the two halves are of equal size and the 0 is in the second half, not the first, then return the second half
	if ( num_in_second_half == num_in_first_half ) {
		const auto &first_half_entries = arg_alignment_split.get_first_half_entries();
		if ( ! contains( first_half_entries, 0_z ) ) {
			return make_opposite_version( arg_alignment_split );
		}
	}

	// Otherwise return this first half
	return arg_alignment_split;
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
bool cath::align::detail::is_valid_split(const alignment_split &arg_alignment_split ///< TODOCUMENT
                                         ) {
	const size_t num_entries        = arg_alignment_split.get_num_entries();
	const size_t num_in_first_half  = arg_alignment_split.get_num_first_half_entries();
	const size_t num_in_second_half = num_entries - num_in_first_half;
	return ( num_in_first_half > 0 && num_in_second_half > 0 );
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
alignment_split cath::align::detail::make_single_alignment_split(const size_t &arg_first_half_entry_index, ///< TODOCUMENT
                                                                 const size_t &arg_num_entries             ///< TODOCUMENT
                                                                 ) {
	alignment_split new_split( arg_num_entries );
	new_split.add_first_half_entry( arg_first_half_entry_index );
	return new_split;
}

/// \brief TODOCUMENT
///
/// \relates alignment_split
alignment_split cath::align::detail::make_alignment_split(const size_vec &arg_first_half_entry_indices, ///< TODOCUMENT
                                                          const size_t   &arg_num_entries               ///< TODOCUMENT
                                                          ) {
	alignment_split new_split( arg_num_entries );
	for (const size_t &first_half_entry_index : arg_first_half_entry_indices) {
		new_split.add_first_half_entry( first_half_entry_index );
	}
	return new_split;
}

///// \brief TODOCUMENT
/////
///// \relates alignment_split
//size_t cath::align::detail::get_split_length(const alignment            &arg_alignment,       ///< TODOCUMENT
//                                             const protein_list         &arg_proteins,        ///< TODOCUMENT
//                                             const alignment_split      &arg_alignment_split, ///< TODOCUMENT
//                                             const alignment_split_half &arg_split_half       ///< TODOCUMENT
//                                             ) {
//	// * IF THIS HALF ONLY HAS ONE ENTRY, THEN BE SURE TO TAKE LENGTH FROM PROTEIN *
//}
