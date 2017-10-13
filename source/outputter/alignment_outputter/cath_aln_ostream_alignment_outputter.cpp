/// \file
/// \brief The cath_aln_ostream_alignment_outputter class definitions

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

#include "cath_aln_ostream_alignment_outputter.hpp"

#include "alignment/alignment_context.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/pair_alignment.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "chopping/region/region.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"
#include "score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::score;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_outputter> cath_aln_ostream_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void cath_aln_ostream_alignment_outputter::do_output_alignment(const alignment_context &arg_alignment_context, ///< TODOCUMENT
                                                               ostream                 &arg_ostream            ///< TODOCUMENT
                                                               ) const {
	check_alignment_is_a_pair( arg_alignment_context.get_alignment() );

	const alignment     &the_alignment = arg_alignment_context.get_alignment();
	const name_set_list &name_sets     = get_name_sets( arg_alignment_context );
	const pdb_list      &pdbs          = get_pdbs     ( arg_alignment_context );
	const protein_list   proteins      = build_protein_list_of_pdb_list_and_names( pdbs, name_sets );
	const alignment      scored_aln    = score_alignment_copy( residue_scorer(), the_alignment, proteins );

//	const aligned_pair_score_list the_score_list = make_full_aligned_pair_score_list();
//	cerr << "make_full_aligned_pair_score_list() returns : " << endl;
//	for (const aligned_pair_score &the_score : the_score_list) {
//		cerr << "\t" << the_score.human_friendly_short_name() << endl;
//	}
//	cerr << endl;
//	cerr << endl;

	const protein &protein_a = proteins[0];
	const protein &protein_b = proteins[1];
	const aligned_pair_score_value_list the_scores_and_values = make_aligned_pair_score_value_list(
		make_full_aligned_pair_score_list(),
		scored_aln,
		protein_a,
		protein_b
	);
	arg_ostream << score_value_list_json_outputter( the_scores_and_values );
}

/// \brief TODOCUMENT
bool cath_aln_ostream_alignment_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Get a name for this alignment_outputter
string cath_aln_ostream_alignment_outputter::do_get_name() const {
	return "CATH";
}
