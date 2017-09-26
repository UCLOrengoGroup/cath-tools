/// \file
/// \brief The alignment_acquirer class definitions

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

#include "alignment_acquirer.hpp"

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/residue_name_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "common/clone/check_uptr_clone_against_this.hpp"
#include "common/logger.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "exception/runtime_error_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "options/options_block/alignment_input_options_block.hpp"
#include "options/options_block/alignment_input_spec.hpp"
#include "superposition/superpose_orderer.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

// using boost::lexical_cast;

constexpr size_t alignment_acquirer::MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<alignment_acquirer> alignment_acquirer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> alignment_acquirer::get_alignment_and_spanning_tree(const pdb_list &arg_pdbs
                                                                                        ) const {
	// Call the concrete class's implementation of do_get_alignment_and_orderer() and grab the resulting alignment and superpose_orderer
	const pair<alignment, superpose_orderer> alignment_and_orderer = do_get_alignment_and_orderer(arg_pdbs);
	const alignment         &new_alignment = alignment_and_orderer.first;
	const superpose_orderer &my_orderer    = alignment_and_orderer.second;

	// Check that both are of the correct size
	const size_t num_pdbs = arg_pdbs.size();
	if ( new_alignment.num_entries() != num_pdbs ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Number of entries in alignment ("
			+ ::std::to_string( new_alignment.num_entries() )
			+ ") does not match expected number ("
			+ ::std::to_string( num_pdbs                    )
			+ ")"
		));
	}
	if ( my_orderer.get_num_items()  != num_pdbs  ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Number of entries in superpose_orderer ("
			+ ::std::to_string( my_orderer.get_num_items() )
			+ ") does not match expected number ("
			+ ::std::to_string( num_pdbs                    )
			+ ")"
		));
	}

	// Try to construct a spanning tree
	// Catch any failures and report them sensibly
	size_size_pair_vec spanning_tree;
	try {
//		spanning_tree = get_spanning_tree( my_orderer );
		spanning_tree = get_spanning_tree_ordered_by_desc_score( my_orderer );
	}
	/// \todo Make the condition of this error message come from the method by which the orderer was populated
	catch (const invalid_argument_exception &err) {
		logger::log_and_exit(
			logger::return_code::INSUFFICIENT_RESIDUE_NAME_OVERLAPS,
			"Cannot construct a tree connecting all PDBs with pairs overlapping by at least " + ::std::to_string( alignment_acquirer::MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR ) + " residues"
		);
	}

	// If the size of the spanning tree isn't one less than the number of PDBs then throw a wobbly
	//
	/// \todo Also check the size of the alignment matches?
	if ( spanning_tree.size() + 1 != num_pdbs ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The spanning tree does not correctly cover the number of PDBs"));
	}

	// Return the resulting alignment and spanning tree
	return { new_alignment, spanning_tree };
}

/// \brief Construct suitable alignment_acquirer objects implied by the specified alignment_input_spec
///
/// NOTE: Keep this code in sync with get_num_acquirers()
///
/// \relates alignment_input_spec
uptr_vec<alignment_acquirer> cath::align::get_alignment_acquirers(const alignment_input_spec &arg_alignment_input_spec ///< The alignment_input_spec to query
                                                                  ) {
	uptr_vec<alignment_acquirer> alignment_acquirers;

	// If the alignment is to be created by reading a legacy CORA alignment file, do that
	// If the alignment is to be created by aligning using residue names,         do that
	// If the alignment is to be created by reading a FASTA alignment file,       do that
	// If the alignment is to be created by reading a legacy SSAP alignment file, do that
	// If the alignment is to be created by reading a list of SSAP scores,        do that
	if ( ! arg_alignment_input_spec.get_cora_alignment_file().empty()  ) {
		alignment_acquirers.push_back( make_unique< cora_aln_file_alignment_acquirer    >( arg_alignment_input_spec.get_cora_alignment_file()  ) );
	}
	if (   arg_alignment_input_spec.get_residue_name_align()           ) {
		alignment_acquirers.push_back( make_unique< residue_name_alignment_acquirer     >(                                                     ) );
	}
	if ( ! arg_alignment_input_spec.get_fasta_alignment_file().empty() ) {
		alignment_acquirers.push_back( make_unique< fasta_aln_file_alignment_acquirer   >( arg_alignment_input_spec.get_fasta_alignment_file() ) );
	}
	if ( ! arg_alignment_input_spec.get_ssap_alignment_file().empty()  ) {
		alignment_acquirers.push_back( make_unique< ssap_aln_file_alignment_acquirer    >( arg_alignment_input_spec.get_ssap_alignment_file()  ) );
	}
	if ( ! arg_alignment_input_spec.get_ssap_scores_file().empty()     ) {
		alignment_acquirers.push_back( make_unique< ssap_scores_file_alignment_acquirer >( arg_alignment_input_spec.get_ssap_scores_file()     ) );
	}

	if ( alignment_acquirers.size() != get_num_acquirers( arg_alignment_input_spec ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"The number of alignment acquirers "     + ::std::to_string( alignment_acquirers.size() )
			+ " doesn't match the expected number (" + ::std::to_string( get_num_acquirers( arg_alignment_input_spec ) ) + ")"
		));
	}

	return alignment_acquirers;
}

/// \brief Construct suitable alignment_acquirer objects implied by the specified alignment_input_options_block
///
/// \relates alignment_input_options_block
uptr_vec<alignment_acquirer> cath::align::get_alignment_acquirers(const alignment_input_options_block &arg_alignment_input_options_block ///< The alignment_input_options_block to query
                                                                  ) {
	return get_alignment_acquirers( arg_alignment_input_options_block.get_alignment_input_spec() );
}

/// \brief Construct the single alignment_acquirer implied by the specified alignment_input_spec
///        (or throw an invalid_argument_exception if fewer/more are implied)
///
/// \relates alignment_input_spec
unique_ptr<alignment_acquirer> cath::align::get_alignment_acquirer(const alignment_input_spec &arg_alignment_input_spec ///< The alignment_input_spec to query
                                                                   ) {
	const auto alignment_acquirers = get_alignment_acquirers( arg_alignment_input_spec );
	if ( alignment_acquirers.size() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to get alignment_acquirer failed because the number of alignment_acquirers isn't one"));
	}
	return alignment_acquirers.front()->clone();
}
