/// \file
/// \brief The fasta_aln_file_alignment_acquirer class definitions

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

#include "fasta_aln_file_alignment_acquirer.hpp"

#include "alignment/alignment.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/file/open_fstream.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "superposition/superpose_orderer.hpp"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> fasta_aln_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pair<alignment, superpose_orderer> fasta_aln_file_alignment_acquirer::do_get_alignment_and_orderer(const pdb_list &arg_pdbs ///< TODOCUMENT
                                                                                                   ) const {
	// Construct an alignment from the FASTA alignment file
	const protein_list proteins_of_pdbs = build_protein_list_of_pdb_list( arg_pdbs );
	const alignment new_alignment = read_alignment_from_fasta_file( get_fasta_alignment_file(), proteins_of_pdbs, cerr );

	// Construct a superpose_orderer and set arbitrary scores to ensure that the spanning tree will connect the entries
	const size_t num_pdbs      = arg_pdbs.size();
	superpose_orderer my_orderer( num_pdbs );
	for (size_t link_ctr = 1; link_ctr < num_pdbs; ++link_ctr) {
		my_orderer.set_score(link_ctr, 0, 0.0);
	}

	const alignment scored_new_alignment = score_alignment_copy( residue_scorer(), new_alignment, proteins_of_pdbs );

	// Return the results
	return make_pair( scored_new_alignment, my_orderer );
}

/// \brief Ctor for fasta_aln_file_alignment_acquirer
fasta_aln_file_alignment_acquirer::fasta_aln_file_alignment_acquirer(const path &arg_fasta_alignment_file ///< TODOCUMENT
                                                                     ) : fasta_alignment_file(arg_fasta_alignment_file) {
}

/// \brief TODOCUMENT
path fasta_aln_file_alignment_acquirer::get_fasta_alignment_file() const {
	return fasta_alignment_file;
}
