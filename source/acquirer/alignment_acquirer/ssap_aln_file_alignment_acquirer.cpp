/// \file
/// \brief The ssap_aln_file_alignment_acquirer class definitions

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

#include "ssap_aln_file_alignment_acquirer.hpp"

#include "alignment/alignment.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/residue_score/residue_scorer.hpp" // ***** TEMPORARY *****
#include "common/clone/make_uptr_clone.hpp"
#include "common/file/open_fstream.hpp"
#include "common/logger.hpp"
#include "exception/runtime_error_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/protein/protein.hpp" // ***** TEMPORARY *****
#include "structure/protein/protein_list.hpp" // ***** TEMPORARY *****
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp" // ***** TEMPORARY *****
#include "structure/protein/sec_struc_planar_angles.hpp" // ***** TEMPORARY *****
#include "superposition/superpose_orderer.hpp"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;

using std::cerr;
using std::ifstream;
using std::make_pair;
using std::pair;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> ssap_aln_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pair<alignment, superpose_orderer> ssap_aln_file_alignment_acquirer::do_get_alignment_and_orderer(const pdb_list &arg_pdbs ///< TODOCUMENT
                                                                                                  ) const {
	const size_t num_pdbs = arg_pdbs.size();
	if ( num_pdbs != 2 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Superposing with a SSAP alignment requires exactly two PDBs"));
	}
	ifstream my_aln_stream;
	open_ifstream(my_aln_stream, get_ssap_alignment_file());
	alignment new_alignment( num_pdbs );
	superpose_orderer my_orderer( num_pdbs );
	try {
		new_alignment = read_alignment_from_cath_ssap_legacy_format(
			my_aln_stream,
			arg_pdbs[alignment::PAIR_A_IDX],
			arg_pdbs[alignment::PAIR_B_IDX],
			cerr
		);
		my_aln_stream.close();
		// Set an arbitrary score to ensure that the spanning tree will connect the two entries
		my_orderer.set_score(alignment::PAIR_B_IDX, alignment::PAIR_A_IDX, 0.0);
	}
	catch (const runtime_error_exception &er) {
		logger::log_and_exit(
			logger::return_code::UNABLE_TO_LOAD_SSAP_LEGACY_ALIGNMENT,
			string("Unable to load SSAP legacy alignment file (a frequent cause of this is that the PDB files don't correspond to the alignment or are in the wrong order). ") +
			"Details: \"" + er.what() + "\""
		);
	}

	 const protein_list proteins_of_pdbs = build_protein_list_of_pdb_list( arg_pdbs );
	 score_alignment( residue_scorer(), new_alignment, proteins_of_pdbs );

	// Return the results
	return make_pair( new_alignment, my_orderer );
}

/// \brief Ctor for ssap_aln_file_alignment_acquirer
ssap_aln_file_alignment_acquirer::ssap_aln_file_alignment_acquirer(const path &arg_ssap_alignment_file ///< TODOCUMENT
                                                                   ) : ssap_alignment_file( arg_ssap_alignment_file ) {
}

/// \brief TODOCUMENT
path ssap_aln_file_alignment_acquirer::get_ssap_alignment_file() const {
	return ssap_alignment_file;
}
