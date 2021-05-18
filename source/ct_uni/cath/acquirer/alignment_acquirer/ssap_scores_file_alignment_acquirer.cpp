/// \file
/// \brief The ssap_scores_file_alignment_acquirer class definitions

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

#include "ssap_scores_file_alignment_acquirer.hpp"

#include <filesystem>
#include <fstream>
#include <utility>

#include <boost/range/adaptor/transformed.hpp>

#include <spdlog/spdlog.h>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/alignment_action.hpp"
#include "cath/alignment/alignment_action.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/alignment/io/outputter/horiz_align_outputter.hpp" /// *** TEMPORARY? ***
#include "cath/alignment/residue_score/residue_scorer.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/graph/spanning_tree.hpp"
#include "cath/common/boost_addenda/graph/spanning_tree.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/file/name_set/name_set_list.hpp"
#include "cath/file/pdb/protein_info.hpp"
#include "cath/file/ssap_scores_file/ssap_scores_file.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"               /// *** TEMPORARY? ***
#include "cath/structure/protein/sec_struc_planar_angles.hpp" /// *** TEMPORARY? ***

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::std::cerr;
using ::std::filesystem::path;
using ::std::make_pair;
using ::std::pair;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> ssap_scores_file_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Specify that this does require backbone-complete input
bool ssap_scores_file_alignment_acquirer::do_requires_backbone_complete_input() const {
	return true;
}

/// \brief Get the aln_glue_style that should be used to implement the specified align_refining
static inline constexpr aln_glue_style aln_glue_style_of_align_refining(const align_refining &prm_align_refining ///< How much refining should be done to the alignment
                                                                        ) {
	return
		( prm_align_refining == align_refining::NO    ) ? aln_glue_style::SIMPLY                           :
		( prm_align_refining == align_refining::LIGHT ) ? aln_glue_style::INCREMENTALLY_WITH_PAIR_REFINING :
		                                                  aln_glue_style::WITH_HEAVY_REFINING;
}

static_assert( aln_glue_style_of_align_refining( align_refining::NO    ) == aln_glue_style::SIMPLY,                           "" );
static_assert( aln_glue_style_of_align_refining( align_refining::LIGHT ) == aln_glue_style::INCREMENTALLY_WITH_PAIR_REFINING, "" );
static_assert( aln_glue_style_of_align_refining( align_refining::HEAVY ) == aln_glue_style::WITH_HEAVY_REFINING,              "" );

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> ssap_scores_file_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context &prm_strucs_context, ///< TODOCUMENT
                                                                                                            const align_refining &prm_align_refining  ///< How much refining should be done to the alignment
                                                                                                            ) const {
	// Parse the SSAP scores file
	const pdb_list &the_pdbs         = prm_strucs_context.get_pdbs();
	const size_t    num_pdbs         = the_pdbs.size();
	const path      ssaps_filename   = get_ssap_scores_file();
	const auto      ssap_scores_data = ssap_scores_file::parse_ssap_scores_file( ssaps_filename );
	const str_vec  &names            = ssap_scores_data.first;
	const auto     &scores           = ssap_scores_data.second;

	if ( names.size() != num_pdbs ) {
		if ( !names.empty() && num_pdbs != 1 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"The number of PDBs is "
				+ ::std::to_string( num_pdbs         )
				+ ", which doesn't match the "
				+ ::std::to_string( names.size() )
				+ " structures required for combining with the SSAP scores file \""
				+ ssaps_filename.string()
				+ "\""
			));
		}
	}

	// Construct the new alignment
	const auto aln_and_spantree = build_multi_alignment(
		the_pdbs,
		names,
		scores,
		ssaps_filename.parent_path(),
		aln_glue_style_of_align_refining( prm_align_refining ),
		ostream_ref{ cerr }
	);
	const alignment          &new_alignment = aln_and_spantree.first;
	const size_size_pair_vec &spanning_tree = aln_and_spantree.second;

	// TODOCUMENT
	if ( names.empty() ) {
		// Return the results
		return make_pair( new_alignment, spanning_tree );
	}

	// ::spdlog::warn( "About to attempt to build protein list using data that's been read from ssaps_filename (with {} pdbs and {} names)",
	//                 num_pdbs,
	//                 names.size() );

	const protein_list proteins_of_pdbs     = build_protein_list_of_pdb_list_and_names(
		prm_strucs_context.get_pdbs(),
		build_name_set_list( names )
	);
	const alignment    scored_new_alignment = score_alignment_copy( residue_scorer(), new_alignment, proteins_of_pdbs );

//	cerr << "Did generate alignment : \n";
//	cerr << horiz_align_outputter( scored_new_alignment ) << endl;
//	write_alignment_as_fasta_alignment( cerr, scored_new_alignment, build_protein_list_of_pdb_list( prm_strucs_context ) );
//	cerr << endl;

	// Return the results
	return make_pair( scored_new_alignment, spanning_tree );
}

/// \brief Ctor for ssap_scores_file_alignment_acquirer
ssap_scores_file_alignment_acquirer::ssap_scores_file_alignment_acquirer(path prm_ssap_scores_filename ///< TODOCUMENT
                                                                         ) : ssap_scores_filename { std::move( prm_ssap_scores_filename ) } {
}

/// \brief TODOCUMENT
path ssap_scores_file_alignment_acquirer::get_ssap_scores_file() const {
	return ssap_scores_filename;
}

/// \brief Build an alignment between the specified PDBs & names using the specified scores and directory of SSAP alignments
pair<alignment, size_size_pair_vec> cath::align::build_multi_alignment(const pdb_list               &prm_pdbs,           ///< The PDBs to be aligned
                                                                       const str_vec                &prm_names,          ///< The names of the structures to be aligned
                                                                       const size_size_doub_tpl_vec &prm_scores,         ///< The SSAP scores between the structures
                                                                       const path                   &prm_alignments_dir, ///< The directory containing alignments for the structures
                                                                       const aln_glue_style         &prm_aln_glue_style, ///< The approach that should be used for glueing alignments together
                                                                       const ostream_ref_opt        &prm_ostream         ///< An (optional reference_wrapper of an) ostream to which warnings/errors should be written
                                                                       ) {
	const protein_list prots = build_protein_list_of_pdb_list( prm_pdbs );
	auto aln_and_spantree = build_alignment(
		prots,
		prm_scores,
		prm_aln_glue_style,

		// Lambda to define how to get the alignment corresponding to the pair of proteins
		// corresponding to the specified indices
		[&] (const size_t  &prm_index_a, //< The index of the first  protein for which the alignment is required
		     const size_t  &prm_index_b  //< The index of the second protein for which the alignment is required
		     ) {
			return read_alignment_from_cath_ssap_legacy_format(
				prm_alignments_dir / ( prm_names[ prm_index_a ] + prm_names[ prm_index_b ] + ".list" ),
				prots[ prm_index_a ],
				prots[ prm_index_b ],
				prm_ostream
			);
		}
	);

	return {
		std::move( aln_and_spantree.first ),
		get_edges_of_spanning_tree( aln_and_spantree.second )
	};
}
