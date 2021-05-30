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

#include <filesystem>
#include <fstream>
#include <utility>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/alignment/residue_score/residue_scorer.hpp" // ***** TEMPORARY *****
#include "cath/common/boost_addenda/graph/spanning_tree.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/logger.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/structure/protein/protein.hpp" // ***** TEMPORARY *****
#include "cath/structure/protein/protein_list.hpp" // ***** TEMPORARY *****
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp" // ***** TEMPORARY *****
#include "cath/structure/protein/sec_struc_planar_angles.hpp" // ***** TEMPORARY *****

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::std::filesystem::path;
using ::std::ifstream;
using ::std::make_pair;
using ::std::pair;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> ssap_aln_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Specify that this does require backbone-complete input
bool ssap_aln_file_alignment_acquirer::do_requires_backbone_complete_input() const {
	return true;
}

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> ssap_aln_file_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context  &prm_strucs_context, ///< TODOCUMENT
                                                                                                         const ostream_ref_opt &prm_ostream         ///< An (optional reference_wrapper of an) ostream to which warnings/errors should be written
                                                                                                         ) const {
	const auto   &the_pdbs = prm_strucs_context.get_pdbs();
	const size_t  num_pdbs = the_pdbs.size();
	if ( num_pdbs != 2 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Superposing with a SSAP alignment requires exactly two PDBs"));
	}
	ifstream my_aln_stream = open_ifstream( get_ssap_alignment_file() );
	alignment new_alignment( num_pdbs );
	try {
		new_alignment = read_alignment_from_cath_ssap_legacy_format(
			my_aln_stream,
			the_pdbs[alignment::PAIR_A_IDX],
			the_pdbs[alignment::PAIR_B_IDX],
			prm_ostream
		);
		my_aln_stream.close();
	}
	catch (const runtime_error_exception &er) {
		logger::log_and_exit(
			logger::return_code::UNABLE_TO_LOAD_SSAP_LEGACY_ALIGNMENT,
			string("Unable to load SSAP legacy alignment file (a frequent cause of this is that the PDB files don't correspond to the alignment or are in the wrong order). ") +
			"Details: \"" + er.what() + "\""
		);
	}

	const protein_list proteins_of_pdbs = build_protein_list( prm_strucs_context );
	score_alignment( residue_scorer(), new_alignment, proteins_of_pdbs );

	// Return the results
	return make_pair(
		new_alignment,
		make_simple_unweighted_spanning_tree( num_pdbs )
	);
}

/// \brief Ctor for ssap_aln_file_alignment_acquirer
ssap_aln_file_alignment_acquirer::ssap_aln_file_alignment_acquirer(path prm_ssap_alignment_file ///< TODOCUMENT
                                                                   ) : ssap_alignment_file( std::move( prm_ssap_alignment_file ) ) {
}

/// \brief TODOCUMENT
path ssap_aln_file_alignment_acquirer::get_ssap_alignment_file() const {
	return ssap_alignment_file;
}
