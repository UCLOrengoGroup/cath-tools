/// \file
/// \brief The multi_align_builder class definitions

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

#include "multi_align_builder.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>

#include "alignment/detail/multi_align_group.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "alignment/io/outputter/horiz_align_outputter.hpp" // ***** TEMPORARY *****
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/invert_permutation.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::algorithm::any_of;
using boost::lexical_cast;

/// \brief Return the index of the group containing the specified entry
size_t multi_align_builder::find_group_of_entry(const size_t &arg_index ///< The index of the entry to be located
                                                ) const {
	return group_index_of_entry[ arg_index ];
}

/// \brief Update the group_index_of_entry to reflect the members of the group with the specified group
void multi_align_builder::update_group_index_of_entry(const size_t &arg_index_of_enlarged_group ///< The index of the group for which group_index_of_entry should be updated
                                                      ) {
	const multi_align_group &enlarged_group = groups[arg_index_of_enlarged_group];
	const size_vec          &groups_entries = enlarged_group.get_entries();
	for (const size_t &groups_entry : groups_entries) {
		group_index_of_entry[groups_entry] = arg_index_of_enlarged_group;
	}
}

/// \brief Ctor for multi_align_builder
multi_align_builder::multi_align_builder(const size_t &arg_num_entries ///< The number of entries to be joined together by alignments
                                         ) {
	groups.reserve              ( arg_num_entries );
	group_index_of_entry.reserve( arg_num_entries );
	// Populate the entries of group_index_of_index with their own indices
	for (const size_t &entry_ctr : indices( arg_num_entries ) ) {
		groups.push_back              ( multi_align_group( entry_ctr) );
		group_index_of_entry.push_back(                    entry_ctr  );
	}
}

/// \brief Get a list of indices of groups that are currently active
size_set multi_align_builder::get_active_groups() const {
	return size_set(
		common::cbegin( group_index_of_entry ),
		common::cend  ( group_index_of_entry )
	);
}

/// \brief Return a const reference to the group with the specified index
const multi_align_group & multi_align_builder::get_group_of_index(const size_t &arg_index ///< The index of the group to be returned
                                                                  ) const {
	return groups[arg_index];
}

/// \brief Add an alignment branch to join two entries that currently exist in separate groups
void multi_align_builder::add_alignment_branch(const size_t       &arg_entry_a,     ///< The index of the first  entry to be joined
                                               const size_t       &arg_entry_b,     ///< The index of the second entry to be joined
                                               const alignment    &arg_alignment,   ///< The alignment between the two entries (with entries in the same order)
                                               const protein_list &/*arg_proteins*/     ///< The PDBs associated with the entries whose alignment is being built
                                               ) {
	// Find the groups currently containing the two specified entries
	const size_t group_index_a = find_group_of_entry( arg_entry_a );
	const size_t group_index_b = find_group_of_entry( arg_entry_b );

	// If the two entries are already in the same group then throw an exception
	if ( group_index_a == group_index_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot add alignment between entry "
			+ lexical_cast<string>( arg_entry_a )
			+ " and "
			+ lexical_cast<string>( arg_entry_b )
			+ " because they are already aligned together"
		));
	}

	// Get (non-const, const) references to the two groups
	multi_align_group       &group_a = groups[ group_index_a ];
	const multi_align_group &group_b = groups[ group_index_b ];

	// Add entry_b to group_a:
	//  (ie glue the joining alignment into group_a (by identifying the left half with entry_a) )
	glue_in_alignment( group_a, arg_alignment, arg_entry_a, arg_entry_b );

	// Add the rest of group_b to group_a:
	//  (ie glue group_b into group_a by identifying entry_b from the two groups)
	group_a.glue_in_copy_of_group( group_b, arg_entry_b );

	// // Refine the join between the two groups
	// group_a.refine_join(
	// 	the_refiner,
	// 	make_subset_protein_list( arg_proteins, group_a.get_entries() ),
	// 	gap_penalty( 50, 0 )
	// );

	// // Refine the alignment
	// group_a.refine_alignment(
	// 	the_refiner,
	// 	make_subset_protein_list( arg_proteins, group_a.get_entries() ),
	// 	gap_penalty( 50, 0 )
	// );

	// Update all the group indices of the entries in group_a
	update_group_index_of_entry( group_index_a );

//	// ***** TEMPORARY *****
//	cerr << "DEBUG: Added alignment branch between " << arg_entry_a << " and " << arg_entry_b;
//	cerr << endl;
//	cerr << "...and the new status is..." << endl;
//	cerr << *this;
//	cerr << endl;
}

/// \brief Get the final alignment from the multi_align_builder
///
/// \pre The multi_align_builder must not be empty
///      else an invalid_argument_exception will be thrown
///
/// \pre All entries must have been built into one active group (ie `get_active_groups().size() == 1` )
///      else an invalid_argument_exception will be thrown
alignment multi_align_builder::get_alignment() const {
	// Sanity check that this isn't empty
	if ( group_index_of_entry.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get alignment of multi_align_builder with no entries"));
	}

	// Grab the group of the first entry
	const size_t group_index = group_index_of_entry.front();

	// Sanity check all elements of group_index_of_entry are group_index
	if ( any_of( group_index_of_entry, [&] (const size_t &x) { return ( x != group_index ); } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get alignment of multi_align_builder that has not been joined into one group"));
	}

	// Get the entries and the alignment from the single active group,
	// then use entries to permute the alignment back to sequential ordering
	const size_vec  &the_temp_entries   = groups[ group_index_of_entry[ group_index ] ].get_entries();
	const alignment &the_temp_alignment = groups[ group_index_of_entry[ group_index ] ].get_alignment();
	const alignment  permuted_alignment = make_permuted_alignment( the_temp_alignment, the_temp_entries );

//	cerr << "Alignment before permutation is : " << horiz_align_outputter( the_temp_alignment ) << endl;
//	cerr << "Alignment after  permutation is : " << horiz_align_outputter( permuted_alignment ) << endl;
//	cerr << "Permutation is                  :";
//	for (const size_t &permutation_value : the_temp_entries) {
//		cerr << " " << permutation_value;
//	}
//	cerr << endl;

//	cerr << "Final alignment is : " << horiz_align_outputter( permuted_alignment ) << endl;

	// Return the result
	return permuted_alignment;
}

/// \brief Simple insertion operator for multi_align_builder
///
/// \relates multi_align_builder
ostream & cath::align::detail::operator<<(ostream                   &arg_os,                 ///< The ostream to which the multi_align_builder should be output
                                          const multi_align_builder &arg_multi_align_builder ///< The multi_align_builder to output
                                          ) {
	const size_set active_groups = arg_multi_align_builder.get_active_groups();
	arg_os << "multi_align_builder[";
	arg_os << active_groups.size();
	arg_os << " active groups : ";
	for ( const size_t &group_index : active_groups ) {
		arg_os << "\n\t";
		arg_os << arg_multi_align_builder.get_group_of_index( group_index );
	}
	arg_os << ( ( ! active_groups.empty() ) ? "\n" : "" );
	arg_os << "]";
	return arg_os;
}

/// \brief Convenience function for adding an alignment branch to a multi_align_builder with a size_size_alignment_tuple
///
/// \relates multi_align_builder
void cath::align::detail::add_alignment_branch(multi_align_builder             &arg_multi_align_builder, ///< The multi_align_builder to which the branch should be added
                                               const size_size_alignment_tuple &arg_branch_alignment,    ///< The details of the branch to be added
                                               const protein_list              &arg_proteins_list        ///< The PDBs associated
                                               ) {
	arg_multi_align_builder.add_alignment_branch(
		get<0>( arg_branch_alignment ),
		get<1>( arg_branch_alignment ),
		get<2>( arg_branch_alignment ),
		arg_proteins_list
	);
}
