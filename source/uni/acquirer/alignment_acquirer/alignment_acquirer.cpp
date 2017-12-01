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
#include "acquirer/alignment_acquirer/do_the_ssaps_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/residue_name_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.hpp"
#include "acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "alignment/options_block/alignment_input_options_block.hpp"
#include "alignment/options_block/alignment_input_spec.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/clone/check_uptr_clone_against_this.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/logger.hpp"
#include "file/strucs_context.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;

using std::make_unique;
using std::pair;
using std::unique_ptr;

constexpr size_t alignment_acquirer::MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<alignment_acquirer> alignment_acquirer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
///
/// Acquire an alignment for the specified PDBs and return it along with a spanning tree that can be
/// used for superposing the corresponding structures
///
/// The order of the edges in the spanning tree doesn't matter and the scores are discarded
/// because none of that matters for superposing
pair<alignment, size_size_pair_vec> alignment_acquirer::get_alignment_and_spanning_tree(const strucs_context &arg_strucs_context, ///< The structures to which the alignment should correspond
                                                                                        const align_refining &arg_align_refining  ///< How much refining should be done to the alignment
                                                                                        ) const {
	using std::to_string;

	// Call the concrete class's implementation of do_get_alignment_and_orderer() and grab the resulting alignment and spanning_tree
	const pair<alignment, size_size_pair_vec> alignment_and_orderer = do_get_alignment_and_spanning_tree( arg_strucs_context, arg_align_refining );

	const size_t num_alignment_entries = alignment_and_orderer.first.num_entries();
	const size_t num_span_tree_entries = alignment_and_orderer.second.size();

	// Check that both are of the correct size
	const size_t num_pdbs = size( arg_strucs_context );
	if ( num_alignment_entries != num_pdbs ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Number of entries in alignment ("
			+ to_string( num_alignment_entries )
			+ ") does not match expected number ("
			+ to_string( num_pdbs              )
			+ ")"
		));
	}
	if ( ( num_span_tree_entries > 0 || num_pdbs > 0 ) && num_span_tree_entries + 1 != num_pdbs ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Number of entries in spanning tree between aligned entries ("
			+ to_string( num_span_tree_entries                                         )
			+ ") does not match expected number ("
			+ to_string( num_pdbs                                                      )
			+ ") - perhaps there aren't enough PDBs with pairs overlapping by at least "
			+ to_string( alignment_acquirer::MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR )
			+ " residues"
		));
	}

	// Return the resulting alignment and spanning tree
	return alignment_and_orderer;
}

/// \brief Construct suitable alignment_acquirer objects implied by the specified alignment_input_spec
///
/// NOTE: Keep this code in sync with get_num_acquirers()
///
/// \relates alignment_input_spec
uptr_vec<alignment_acquirer> cath::align::get_alignment_acquirers(const alignment_input_spec &arg_alignment_input_spec ///< The alignment_input_spec to query
                                                                  ) {
	using std::to_string;
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
	if ( arg_alignment_input_spec.get_do_the_ssaps_dir() ) {
		alignment_acquirers.push_back( make_unique< do_the_ssaps_alignment_acquirer     >( *arg_alignment_input_spec.get_do_the_ssaps_dir()    ) );
	}

	if ( alignment_acquirers.size() != get_num_acquirers( arg_alignment_input_spec ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"The number of alignment acquirers "     + to_string( alignment_acquirers.size() )
			+ " doesn't match the expected number (" + to_string( get_num_acquirers( arg_alignment_input_spec ) ) + ")"
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

	// If no alignment_acquirer has been specified then use a do_the_ssaps_alignment_acquirer
	if ( alignment_acquirers.empty() ) {
		return make_unique< do_the_ssaps_alignment_acquirer >();
	}

	if ( alignment_acquirers.size() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to get alignment_acquirer failed because the number of alignment_acquirers isn't one"));
	}
	return alignment_acquirers.front()->clone();
}
