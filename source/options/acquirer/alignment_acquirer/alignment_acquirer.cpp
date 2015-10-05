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

#include "alignment_acquirer.h"

#include <boost/lexical_cast.hpp>

#include "alignment/alignment.h"
#include "common/clone/check_uptr_clone_against_this.h"
#include "common/logger.h"
#include "exception/invalid_argument_exception.h"
#include "exception/runtime_error_exception.h"
#include "superposition/superpose_orderer.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/pdbs_acquirer/pdbs_acquirer.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::lexical_cast;

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
			+ lexical_cast<string>( new_alignment.num_entries() )
			+ ") does not match expected number ("
			+ lexical_cast<string>( num_pdbs                    )
			+ ")"
		));
	}
	if ( my_orderer.get_num_items()  != num_pdbs  ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Number of entries in superpose_orderer ("
			+ lexical_cast<string>( my_orderer.get_num_items() )
			+ ") does not match expected number ("
			+ lexical_cast<string>( num_pdbs                    )
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
			"Cannot construct a tree connecting all PDBs with pairs overlapping by at least " + lexical_cast<string>(alignment_acquirer::MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR) + " residues"
		);
	}

	// If the size of the spanning tree isn't one less than the number of PDBs then throw a wobbly
	//
	/// \todo Also check the size of the alignment matches?
	if ( spanning_tree.size() + 1 != num_pdbs ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The spanning tree does not correctly cover the number of PDBs"));
	}

	// Return the resulting alignment and spanning tree
	return make_pair(new_alignment, spanning_tree);
}
