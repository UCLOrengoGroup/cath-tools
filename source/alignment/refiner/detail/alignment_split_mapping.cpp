/// \file
/// \brief The alignment_split_mapping class definitions

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

#include "alignment_split_mapping.h"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/lower_bound.hpp>

#include "alignment/alignment.h"
#include "alignment/alignment_row.h"
#include "alignment/pair_alignment.h"
#include "alignment/refiner/detail/alignment_split.h"
#include "common/c++14/cbegin_cend.h"
#include "exception/invalid_argument_exception.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"

#include <iostream> // ***** TEMPORARY *****
#include <iterator>

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;
using boost::range::lower_bound;

/// \brief Ctor for alignment_split_mapping
///
/// \todo Consider splitting part of this out into a function that creates a copy of the alignment with any missing residues reinserted
alignment_split_mapping::alignment_split_mapping(const alignment &arg_alignment,      ///< TODOCUMENT
                                                 const size_set  &arg_entries,        ///< TODOCUMENT
                                                 const size_vec  &arg_correct_lengths ///< TODOCUMENT
                                                 ) : orig_length           ( arg_alignment.length()                                     ),
                                                     orig_num_entries      ( arg_alignment.num_entries()                                ),
                                                     did_insert_entries    ( false                                                      ),
                                                     local_aln             ( arg_entries.size()                                         ),
                                                     idx_of_orig_aln_idx   ( orig_length                                                ),
                                                     orig_aln_entries      ( common::cbegin( arg_entries ), common::cend( arg_entries ) ),
                                                     index_of_pdb_res_index( arg_alignment.num_entries()                                ) {
	// Sanity check the input list of entries
	if ( orig_aln_entries.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("List of entries is empty"));
	}
	if ( orig_aln_entries.back() >= orig_num_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Entries out of range of the entries in the alignment"));
	}

	// Loop along the length of the original alignment
	for (size_t orig_aln_index = 0; orig_aln_index < orig_length; ++orig_aln_index) {

		const alignment_row local_row = get_row_of_entries_of_alignment( arg_alignment, orig_aln_entries, orig_aln_index );

		// Only do anything if at least one of the specified entries is present at this index
		if ( any_entries_present( local_row ) ) {

			// For each entry that's present at this index,
			// insert any missing residues into the local_aln, recording where each has gone
			for (size_t entry = 0; entry < orig_aln_entries.size(); ++entry) {
				// Grab some details for this entry:
				//  * a reference to the relevant entry of index_of_pdb_res_index
				//  * the equivalent entry in the original alignment
				//  * the position for this index/entry in the original alignment
				size_vec          &pdb_res_indices = index_of_pdb_res_index[ entry ];
				const size_t      &orig_aln_entry  = orig_aln_entries      [ entry ];
				const opt_aln_posn position        = arg_alignment.position_of_entry_of_index( orig_aln_entry, orig_aln_index );

				// If there is something present here in the original alignment then add any missing residues
				if ( position ) {
					// Add any missing residues...
					while ( pdb_res_indices.size() < *position ) {
						did_insert_entries = true;
						pdb_res_indices.push_back( local_aln.length() );
						append_row_with_single_value( local_aln, entry, pdb_res_indices.size() - 1 );
					}
				}
			}

			// For each entry that's present at this index,
			// record where the residue for this index's row will go
			for (size_t entry = 0; entry < orig_aln_entries.size(); ++entry) {
				// Grab some details for this entry:
				//  * a reference to the relevant entry of index_of_pdb_res_index
				//  * the equivalent entry in the original alignment
				//  * the position for this index/entry in the original alignment
				size_vec          &pdb_res_indices = index_of_pdb_res_index[ entry ];
				const size_t      &orig_aln_entry  = orig_aln_entries      [ entry ];
				const opt_aln_posn position        = arg_alignment.position_of_entry_of_index( orig_aln_entry, orig_aln_index );

				// If there is something present here in the original alignment then store the residue for this index will go
				if ( position ) {
					// ...and store the residue for this index will go
					pdb_res_indices.push_back( local_aln.length() );
				}
			}

			// Store where this alignment row is placed in the local_aln and then put it there
			idx_of_orig_aln_idx[ orig_aln_index ] = local_aln.length();
			append_row( local_aln, local_row );
		}
	}

	// Append any extra residues that are missing at the end of the original alignment
	for (size_t entry = 0; entry < orig_aln_entries.size(); ++entry) {
		// Grab some details for this entry:
		//  * a reference to the relevant entry of index_of_pdb_res_index
		//  * the equivalent entry in the original alignment
		//  * the correct length for this entry
		size_vec          &pdb_res_indices = index_of_pdb_res_index[ entry ];
		const size_t      &orig_aln_entry  = orig_aln_entries      [ entry ];
		const size_t      &correct_length  = arg_correct_lengths   [ orig_aln_entry ];

		// Add any missing residues...
		while ( pdb_res_indices.size() < correct_length ) {
			did_insert_entries = true;
			pdb_res_indices.push_back( local_aln.length() );
			append_row_with_single_value( local_aln, entry, pdb_res_indices.size() - 1 );
		}
	}


	// Very basic sanity checks
	assert( orig_aln_entries.size()    == local_aln.num_entries() );
	assert( idx_of_orig_aln_idx.size() == orig_length             );
}

/// \brief TODOCUMENT
size_t alignment_split_mapping::orig_aln_length() const {
	return orig_length;
}

/// \brief TODOCUMENT
size_t alignment_split_mapping::orig_aln_num_entries() const {
	return orig_num_entries;
}

/// \brief TODOCUMENT
bool alignment_split_mapping::inserted_entries() const {
	return did_insert_entries;
}

/// \brief TODOCUMENT
size_t alignment_split_mapping::length() const {
	return local_aln.length();
}

/// \brief TODOCUMENT
size_t alignment_split_mapping::num_entries() const {
	return local_aln.num_entries();
}

/// \brief TODOCUMENT
opt_size alignment_split_mapping::index_of_orig_aln_index(const size_t &arg_orig_aln_index ///< TODOCUMENT
                                                          ) const {
	if ( arg_orig_aln_index >= orig_aln_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Index is too large for the original alignment"));
	}
	return idx_of_orig_aln_idx[ arg_orig_aln_index ];
//	if ( ! index ) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception(
//			"Index "
//			+ lexical_cast<string>( arg_orig_aln_index )
//			+ " in original alignment doesn't appear in this alignment_split_mapping"
//		));
//	}
//	return *index;
}

/// \brief TODOCUMENT
opt_size alignment_split_mapping::entry_of_orig_aln_entry(const size_t &arg_orig_aln_entry ///< TODOCUMENT
                                                          ) const {

	if ( arg_orig_aln_entry >= orig_aln_num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Entry is too large for the original alignment"));
	}
	const auto find_itr = lower_bound(
		orig_aln_entries,
		arg_orig_aln_entry
	);
	if ( find_itr == common::cend( orig_aln_entries ) || *find_itr != arg_orig_aln_entry ) {
		return none;
//		BOOST_THROW_EXCEPTION(invalid_argument_exception("Entry in original alignment doesn't appear in this alignment_split_mapping"));
	}
	return numeric_cast<size_t>( distance( common::cbegin( orig_aln_entries ), find_itr ) );
}

/// \brief TODOCUMENT
size_t alignment_split_mapping::orig_aln_entry_of_entry(const size_t &arg_orig_aln_entry ///< TODOCUMENT
                                                        ) const {
	if ( arg_orig_aln_entry > orig_aln_entries.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Entry is out of range"));
	}
	return orig_aln_entries[ arg_orig_aln_entry ];
}

/// \brief TODOCUMENT
opt_aln_posn alignment_split_mapping::position_of_entry_of_index(const size_t &arg_entry, ///< TODOCUMENT
                                                                 const size_t &arg_index  ///< TODOCUMENT
                                                                 ) const {
	return local_aln.position_of_entry_of_index( arg_entry, arg_index );
}

///// \brief TODOCUMENT
/////
///// \relates alignment_split_mapping
//ostream & cath::align::detail::operator<<(ostream                       &arg_os,     ///< TODOCUMENT
//                                          const alignment_split_mapping &arg_mapping ///< TODOCUMENT
//                                          ) {
//	arg_os << "alignment_split_mapping[";
//	arg_os << orig_aln_length() const;
//	size_t orig_aln_num_entries() const;
//	bool inserted_entries() const;
//
//	size_t length() const;
//	size_t num_entries() const;
//
//	size_t index_of_orig_aln_index(const size_t &) const;
//	size_t entry_of_orig_aln_entry(const size_t &) const;
//	size_t orig_aln_entry_of_entry(const size_t &) const;
//
//	opt_aln_posn position_of_entry_of_index(const size_t &,
//	                                        const size_t &) const;
//
//	size_t index_of_protein_index(const size_t &,
//	                              const size_t &) const;
//	arg_os << "]";
//	return arg_os;
//}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
size_vec cath::align::detail::get_orig_entries(const alignment_split_mapping &arg_aln_mapping ///< TODOCUMENT
                                               ) {
	const size_t num_entries = arg_aln_mapping.num_entries();
	size_vec orig_entries;
	orig_entries.reserve( num_entries );
	for (size_t entry = 0; entry < num_entries; ++entry) {
		orig_entries.push_back( arg_aln_mapping.orig_aln_entry_of_entry( entry ) );
	}
	return orig_entries;
}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
alignment_row cath::align::detail::get_row_of_alignment(const alignment_split_mapping &arg_aln_mapping, ///< TODOCUMENT
                                                        const size_t                  &arg_index        ///< TODOCUMENT
                                                        ) {
	const size_t num_entries = arg_aln_mapping.num_entries();
	opt_aln_posn_vec positions;
	positions.reserve( num_entries );
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const opt_aln_posn position = arg_aln_mapping.position_of_entry_of_index( entry_ctr, arg_index );
		positions.push_back( position );
	}
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
bool cath::align::detail::has_position_of_entry_of_index(const alignment_split_mapping &arg_aln_mapping, ///< TODOCUMENT
                                                         const size_t                  &arg_entry,       ///< TODOCUMENT
                                                         const size_t                  &arg_index        ///< TODOCUMENT
                                                         ) {
	return static_cast<bool>( arg_aln_mapping.position_of_entry_of_index( arg_entry, arg_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
size_t cath::align::detail::get_position_of_entry_of_index(const alignment_split_mapping &arg_aln_mapping, ///< TODOCUMENT
                                                           const size_t                  &arg_entry,       ///< TODOCUMENT
                                                           const size_t                  &arg_index        ///< TODOCUMENT
                                                           ) {
	return *arg_aln_mapping.position_of_entry_of_index( arg_entry, arg_index );
}

///// \brief TODOCUMENT
/////
///// \relates alignment_split_mapping
//
//size_vec cath::align::detail::present_orig_aln_entries_of_orig_aln_index(const alignment_split_mapping &arg_aln_mapping,   ///< TODOCUMENT
//                                                                         const size_t                  &arg_orig_aln_index ///< TODOCUMENT
//                                                                         ) {
//	const opt_size index = arg_aln_mapping.index_of_orig_aln_index( arg_orig_aln_index );
//	return index ? present_orig_aln_entries_of_index( arg_aln_mapping, *index )
//	             : size_vec();
//}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
size_vec cath::align::detail::present_orig_aln_entries_of_index(const alignment_split_mapping &arg_aln_mapping, ///< TODOCUMENT
                                                                const size_t                  &arg_index        ///< TODOCUMENT
                                                                ) {
	const size_t orig_aln_num_entries = arg_aln_mapping.orig_aln_num_entries();
	size_vec present_positions;
	for (size_t orig_aln_entry_ctr = 0; orig_aln_entry_ctr < orig_aln_num_entries; ++orig_aln_entry_ctr) {
		const opt_size entry = arg_aln_mapping.entry_of_orig_aln_entry( orig_aln_entry_ctr );
		if ( entry && has_position_of_entry_of_index( arg_aln_mapping, *entry, arg_index ) ) {
			present_positions.push_back( orig_aln_entry_ctr );
		}
	}
	return present_positions;
}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
alignment_split_mapping cath::align::detail::make_alignment_split_mapping(const alignment            &arg_alignment,            ///< TODOCUMENT
                                                                          const alignment_split      &arg_alignment_split,      ///< TODOCUMENT
                                                                          const alignment_split_half &arg_alignment_split_half, ///< TODOCUMENT
                                                                          const size_vec             &arg_correct_lengths       ///< TODOCUMENT
                                                                          ) {
	return alignment_split_mapping(
		arg_alignment,
		entries_of_alignment_split_half( arg_alignment_split, arg_alignment_split_half ),
		arg_correct_lengths
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_split_mapping
alignment cath::align::detail::build_alignment(const alignment               &arg_inter_mapping_alignment, ///< TODOCUMENT
                                               const alignment_split_mapping &arg_mapping_a,               ///< TODOCUMENT
                                               const alignment_split_mapping &arg_mapping_b                ///< TODOCUMENT
                                               ) {
	check_alignment_is_a_pair( arg_inter_mapping_alignment );

	const size_t inter_mapping_aln_length = arg_inter_mapping_alignment.length();

	const size_vec entries_a     = get_orig_entries( arg_mapping_a );
	const size_vec entries_b     = get_orig_entries( arg_mapping_b );
	const size_t   num_entries_a = arg_mapping_a.num_entries();
	const size_t   num_entries_b = arg_mapping_b.num_entries();

	alignment new_alignment( num_entries_a + num_entries_b );
	new_alignment.reserve( inter_mapping_aln_length );

	for (size_t inter_mapping_index = 0; inter_mapping_index < inter_mapping_aln_length; ++inter_mapping_index) {
		const opt_aln_posn a_position = a_position_of_index( arg_inter_mapping_alignment, inter_mapping_index );
		const opt_aln_posn b_position = b_position_of_index( arg_inter_mapping_alignment, inter_mapping_index );
		if ( ! a_position && ! b_position) {
			BOOST_THROW_EXCEPTION(out_of_range("Inter-mapping alignment has an empty row"));
		}
		const alignment_row aln_row_a = a_position ? get_row_of_alignment( arg_mapping_a, *a_position )
		                                           : make_empty_aln_row( num_entries_a );
		const alignment_row aln_row_b = b_position ? get_row_of_alignment( arg_mapping_b, *b_position )
		                                           : make_empty_aln_row( num_entries_b );

		append_row( new_alignment, weave( aln_row_a, entries_a, aln_row_b, entries_b ) );
	}
	return new_alignment;
}
