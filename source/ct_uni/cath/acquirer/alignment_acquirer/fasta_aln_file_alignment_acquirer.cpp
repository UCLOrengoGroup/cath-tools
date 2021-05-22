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

#include <filesystem>
#include <fstream>
#include <utility>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/alignment/residue_score/residue_scorer.hpp"
#include "cath/common/boost_addenda/graph/spanning_tree.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::std;

using ::std::filesystem::path;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> fasta_aln_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Specify that this does not require backbone-complete input
bool fasta_aln_file_alignment_acquirer::do_requires_backbone_complete_input() const {
	return false;
}

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> fasta_aln_file_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context  &prm_strucs_context, ///< TODOCUMENT
                                                                                                          const ostream_ref_opt &prm_ostream         ///< An (optional reference_wrapper of an) ostream to which warnings/errors should be written
                                                                                                          ) const {
	// Construct an alignment from the FASTA alignment file
	const alignment new_alignment = read_alignment_from_fasta_file(
		get_fasta_alignment_file(),
		prm_strucs_context.get_pdbs(),
		get_domain_or_specified_or_from_acq_names( prm_strucs_context.get_name_sets() ),
		prm_ostream ? prm_ostream->get() : cerr
	);

	const protein_list proteins_of_pdbs     = build_protein_list( prm_strucs_context );
	const alignment    scored_new_alignment = score_alignment_copy( residue_scorer(), new_alignment, proteins_of_pdbs );

	// Return the results
	return make_pair(
		scored_new_alignment,
		make_simple_unweighted_spanning_tree( proteins_of_pdbs.size() )
	);
}

/// \brief Ctor for fasta_aln_file_alignment_acquirer
fasta_aln_file_alignment_acquirer::fasta_aln_file_alignment_acquirer(path prm_fasta_alignment_file ///< TODOCUMENT
                                                                     ) : fasta_alignment_file( std::move( prm_fasta_alignment_file ) ) {
}

/// \brief TODOCUMENT
path fasta_aln_file_alignment_acquirer::get_fasta_alignment_file() const {
	return fasta_alignment_file;
}
