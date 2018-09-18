/// \file
/// \brief The cora_aln_file_alignment_acquirer class definitions

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

#include "cora_aln_file_alignment_acquirer.hpp"

#include "alignment/alignment.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "common/boost_addenda/graph/spanning_tree.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/file/open_fstream.hpp"
#include "file/strucs_context.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <fstream>

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;

using boost::filesystem::path;
using std::cerr;
using std::ifstream;
using std::make_pair;
using std::pair;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<alignment_acquirer> cora_aln_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Specify that this does require backbone-complete input
bool cora_aln_file_alignment_acquirer::do_requires_backbone_complete_input() const {
	return true;
}

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> cora_aln_file_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context &prm_strucs_context ///< TODOCUMENT
                                                                                                         ) const {
	// Construct an alignment from the CORA alignment file
	const auto   &the_pdbs = prm_strucs_context.get_pdbs();
	const size_t  num_pdbs = the_pdbs.size();
	ifstream my_aln_stream;
	open_ifstream(my_aln_stream, get_cora_alignment_file());
	const alignment     new_alignment = read_alignment_from_cath_cora_legacy_format( my_aln_stream, the_pdbs, ostream_ref{ cerr } );
	my_aln_stream.close();

	const protein_list proteins_of_pdbs = build_protein_list( prm_strucs_context );
	const alignment scored_new_alignment = score_alignment_copy( residue_scorer(), new_alignment, proteins_of_pdbs );

	// Return the results with a simple spanning tree that will connect the entries
	return make_pair(
		scored_new_alignment,
		make_simple_unweighted_spanning_tree( num_pdbs )
	);
}

/// \brief Ctor for cora_aln_file_alignment_acquirer
cora_aln_file_alignment_acquirer::cora_aln_file_alignment_acquirer(const path     &prm_cora_alignment_file ///< TODOCUMENT
                                                                   ) : cora_alignment_file(prm_cora_alignment_file) {
}

/// \brief TODOCUMENT
path cora_aln_file_alignment_acquirer::get_cora_alignment_file() const {
	return cora_alignment_file;
}
