/// \file
/// \brief The multi_align_group class definitions

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

#include "multi_align_group.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/remove_copy.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

#include "cath/alignment/alignment_action.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/alignment/io/outputter/horiz_align_outputter.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/common/boost_addenda/range/adaptor/lexical_casted.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::align::gap;
using namespace ::cath::common;

using ::boost::lexical_cast;
using ::boost::numeric_cast;
using ::boost::range::find;
using ::boost::range::max_element;
using ::boost::range::push_back;
using ::boost::range::remove_copy;
using ::std::ostream;
using ::std::string;

/// \brief TODOCUMENT
void multi_align_group::refine_join(alignment_refiner  &prm_alignment_refiner, ///< TODOCUMENT
                                    const protein_list &prm_proteins,          ///< TODOCUMENT
                                    const gap_penalty  &prm_gap_penalty,       ///< TODOCUMENT
                                    const size_vec     &prm_join_point         ///< TODOCUMENT
                                    ) {
	if ( entries.empty() || prm_join_point.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot refine join with empty entries or empty join point"));
	}
	// std::cerr << "Alignment has        : " << the_alignment.num_entries()                             << " entries\n" << std::flush;
	// const size_vec &entries_cr = entries;
	// std::cerr << "Entries are          : " << join( entries_cr        | lexical_casted<string>(), ", " ) <<         "\n" << std::flush;
	// std::cerr << "Join point is        : " << join( prm_join_point    | lexical_casted<string>(), ", " ) <<         "\n" << std::flush;
	size_opt_vec index_of_entry( *max_element( entries ) + 1 );
	for (const size_t &entry_index : indices( entries.size() ) ) {
		index_of_entry[ entries[ entry_index ] ] = entry_index;
	}
	const auto mapped_join_point = transform_build<size_vec>(
		prm_join_point,
		[&] (const size_t &x) { return *index_of_entry[ x ]; }
	);

	// std::cerr << "Mapped join point is : " << join( mapped_join_point | lexical_casted<string>(), ", " ) <<         "\n" << std::flush;

	// cerr << "Alignment to be refined (at join point between "
		// << join( mapped_join_point | lexical_casted<string>(), ", " )
		// << " and the rest) : \n";
	// write_alignment_as_fasta_alignment( cerr, the_alignment, prm_proteins );
	// cerr << "\n";

	the_alignment = prm_alignment_refiner.iterate_join( the_alignment, prm_proteins, prm_gap_penalty, mapped_join_point );

	// cerr << "Refined alignment (at join point between "
		// << join( mapped_join_point | lexical_casted<string>(), ", " )
		// << " and the rest) : \n";
	// write_alignment_as_fasta_alignment( cerr, the_alignment, prm_proteins );
	// cerr << "\n";
}

/// \brief Ctor for multi_align_group that creates a single entry of the specified index
multi_align_group::multi_align_group(const size_t &prm_index ///< The index of the single entry in this multi_align_group
                                     ) : the_alignment( 1            ),
                                         entries      ( 1, prm_index ) {
}

/// \brief Ctor for multi_align_group that creates a pair of entries with the specified indices, joined by the specified pair-alignment
multi_align_group::multi_align_group(alignment     prm_pair_alignment, ///< The pair-alignment joining the first and second entry (the indices of which are given below)
                                     const size_t &prm_entry_a,        ///< The index of the first  entry in this multi_align_group
                                     const size_t &prm_entry_b         ///< The index of the second entry in this multi_align_group
                                     ) : the_alignment{ std::move( prm_pair_alignment ) },
                                         entries      { prm_entry_a, prm_entry_b        } {
	check_alignment_is_a_pair( the_alignment );
}

/// \brief Getter for the current alignment of the members of this group
const alignment & multi_align_group::get_alignment() const {
	return the_alignment;
}

/// \brief Getter for the entries that this group currently contains
const size_vec & multi_align_group::get_entries() const {
	return entries;
}

/// \brief Glue (a copy of) another multi_align_group into this one by identifying the single entry that is shared between the two groups
void multi_align_group::glue_in_copy_of_group(const multi_align_group &prm_other_group,       ///< The other multi_align_group to glue into this one
                                              const size_t            &prm_entry_to_identify ///< The single entry that's common to the two multi_align_groups
                                              ) {
	// Grab some details of the two entries and alignments
	const size_t     index_of_entry_in_this = get_index_of_entry( *this,           prm_entry_to_identify );
	const size_t     index_of_entry_in_that = get_index_of_entry( prm_other_group, prm_entry_to_identify );
	const alignment &this_alignment         =                 get_alignment();
	const alignment &that_alignment         = prm_other_group.get_alignment();
	const size_vec  &that_entries           = prm_other_group.get_entries();

//	cerr << "Attempting to glue groups by identifying :\n";
//
//	cerr << "Element " << prm_entry_to_identify;
//	cerr << " index_of_entry_in_this : " << index_of_entry_in_this << ", [";
//	for (const size_t &entry : entries) {
//		cerr << " " << entry;
//	}
//	cerr << " ] this alignment :\n";
//	cerr << horiz_align_outputter( this_alignment );
//	cerr << "\n";
//
//	cerr << "Element " << prm_entry_to_identify;
//	cerr << " index_of_entry_in_that : " << index_of_entry_in_this << ", [";
//	for (const size_t &entry : that_entries) {
//		cerr << " " << entry;
//	}
//	cerr << " ] that alignment :\n";
//	cerr << horiz_align_outputter( that_alignment );
//	cerr << "\n";
//
//	cerr << flush;

	alignment updated_alignment = glue_two_alignments(
		this_alignment,
		index_of_entry_in_this,
		that_alignment,
		index_of_entry_in_that
	);

	// Create a vector of the new entries but not the prm_entry_to_identify
	//
	/// \todo Write a remove_copy_build() wrapper for Boost Range's remove_copy() and change this to:
	///       const auto new_entries = remove_copy_build<size_vec>( entry_from_that, prm_entry_to_identify );
	size_vec new_entries;
	new_entries.reserve( that_entries.size() - 1 );
	remove_copy(
		that_entries,
		back_inserter( new_entries ),
		prm_entry_to_identify
	);

	// Update the alignment and then append the new entries to the back of entries
	the_alignment = updated_alignment;
	push_back( entries, new_entries );
}

/// \brief TODOCUMENT
void multi_align_group::refine_alignment(alignment_refiner  &prm_alignment_refiner, ///< TODOCUMENT
                                         const protein_list &prm_proteins,          ///< TODOCUMENT
                                         const gap_penalty  &prm_gap_penalty        ///< TODOCUMENT
                                         ) {
	// cerr << "Alignment to be refined: \n";
	// write_alignment_as_fasta_alignment( cerr, the_alignment, prm_proteins );
	// cerr << "\n" << flush;

	the_alignment = prm_alignment_refiner.iterate( the_alignment, prm_proteins, prm_gap_penalty );

	// cerr << "Refined alignment: \n";
	// write_alignment_as_fasta_alignment( cerr, the_alignment, prm_proteins );
	// cerr << "\n" << flush;
}

/// \brief Simple insertion operator for multi_align_group
///
/// \relates multi_align_group
ostream & cath::align::detail::operator<<(ostream                 &prm_os,               ///< The ostream to which the multi_align_group should be output
                                          const multi_align_group &prm_multi_align_group ///< The multi_align_group to output
                                          ) {
	const alignment &prm_alignment = prm_multi_align_group.get_alignment();
	const size_vec  &prm_entries   = prm_multi_align_group.get_entries();

	prm_os << "multi_align_group[entries:";
	for (const size_t &entry_num : prm_entries) {
		prm_os << " " << entry_num;
	}
	prm_os << " : ";
	prm_os << horiz_align_outputter( prm_alignment );
	prm_os << "\n]";
	return prm_os;
}

/// \brief Get the index of a particular entry within a multi_align_group
size_t cath::align::detail::get_index_of_entry(const multi_align_group &prm_align_group, ///< The multi_align_group to query
                                               const size_t            &prm_entry        ///< The entry for which the index is required
                                               ) {
	// Grab a cref to the list of entries in the group
	const size_vec &indices = prm_align_group.get_entries();

	// Find prm_entry within the indices
	const auto find_itr = find( indices, prm_entry );

	// If prm_entry wasn't found then throw an exception
	if ( find_itr == cend( indices ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to find entry "
			+ lexical_cast<string>( prm_entry )
			+ " in multi_align_group"
		));
	}

	// Else return the index of the found entry
	return numeric_cast<size_t>( distance( cbegin( indices ), find_itr ) );
}

/// \brief Join another entry to the specified multi_align_group using a pair-alignment to one of the group's existing entries
void cath::align::detail::glue_in_alignment(multi_align_group &prm_align_group,       ///< The group to which an entry should be added
                                            const alignment   &prm_alignment,         ///< The pair-alignment between an existing entry and a new entry
                                            const size_t      &prm_entry_to_identify, ///< The existing entry to identify with the first entry of the pair alignment
                                            const size_t      &prm_other_entry        ///< The entry to add
                                            ) {
	prm_align_group.glue_in_copy_of_group(
		multi_align_group(
			prm_alignment,
			prm_entry_to_identify,
			prm_other_entry
		),
		prm_entry_to_identify
	);
}

