/// \file
/// \brief The residue_scorer class definitions

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

#include "residue_scorer.h"

#include <boost/filesystem/path.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>

#include "alignment/alignment.h"
#include "alignment/io/alignment_io.h"
#include "alignment/residue_score/alignment_residue_scores.h"
#include "ssap/context_res.h"
#include "ssap/ssap.h"
#include "structure/entry_querier/residue_querier.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"

#include <algorithm>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::range::set_intersection;

/// \brief TODOCUMENT
alignment_residue_scores residue_scorer::get_alignment_residue_scores(const alignment    &arg_alignment, ///< TODOCUMENT
                                                                      const protein_list &arg_proteins   ///< TODOCUMENT
                                                                      ) const {
	// Grab the num_entries and length of the alignment and sanity check that there are at least two entries
	const size_t num_entries = arg_alignment.num_entries();
	const size_t length      = arg_alignment.length();
	if ( num_entries <= 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot score alignment with fewer than two entries"));
	}

	float_score_vec_vec numerators  ( num_entries, float_score_vec( length, 0.0 ) );
	float_score_vec_vec denominators( num_entries, float_score_vec( length, 0.0 ) );

//	cerr << "Attempting to residue score : " << endl;
//	cerr << horiz_align_outputter( arg_alignment ) << endl;

	/// Loop down the length of the alignment
	for (size_t from_index = 0; from_index < length; ++from_index) {
		const size_vec from_entries = present_positions_of_index( arg_alignment, from_index );
		if ( from_entries.size() <= 1 ) {
			continue;
		}

		for (size_t to_index = 0; to_index < length; ++to_index) {
			const size_vec to_entries = present_positions_of_index( arg_alignment, to_index );
			if ( to_entries.size() <= 1 ) {
				continue;
			}

			size_vec common_entries;
			common_entries.reserve( min( from_entries.size(), to_entries.size() ) );
			set_intersection(
				from_entries,
				to_entries,
				back_inserter( common_entries )
			);

//			fprintf( stderr, "after common_entries\n" );
//			fprintf( stderr, "from_entries has %ld entries\n", from_entries.size() );
//			fprintf( stderr, "to_entries has %ld entries\n", to_entries.size() );
//			fprintf( stderr, "common_entries has %ld entries\n", common_entries.size() );

//			cerr << "From index\t" << from_index << " to index\t" << to_index << endl;
//			cerr << "\t\tCommon entries are :";
//			for (const size_t &common_entry : common_entries) {
//				cerr << "\t" << common_entry;
//			}
//			cerr << endl;

			// For some reason, using the following commented code lead to a crash on the release build on orengobuild64.
			// The crash didn't occur on debug/relwithdebinfo builds on ob64 or any build on bsmlx62 (Ubuntu).
			// I couldn't see any other problems when using the code without BOOST_FOREACH() (eg valgrind ran clean).
			// When using the BOOST_FOREACH() code, the crash disappeared on adding a bunch of of debug statements.
			//
			//
			//			for (const size_t &common_entry_a : common_entries) {
			////				fprintf(stderr, "common_entry_a   : %ld\n", common_entry_a);
			//				for (const size_t &common_entry_b : common_entries) {
			//					fprintf(stderr, "common_entry_b : %ld\n", common_entry_b);
			//					if ( common_entry_a < common_entry_b ) {
			////						fprintf(stderr, "yes\n");

			const size_t common_entries_size = common_entries.size();
			for (size_t comm_ent_ctr_a = 0; comm_ent_ctr_a < common_entries_size; ++comm_ent_ctr_a) {
//				fprintf(stderr, "comm_ent_ctr_a : %ld\n", comm_ent_ctr_a);
				for (size_t comm_ent_ctr_b = 0; comm_ent_ctr_b < common_entries_size; ++comm_ent_ctr_b) {
//					fprintf(stderr, "comm_ent_ctr_b : %ld\n", comm_ent_ctr_b);
					if ( comm_ent_ctr_a < comm_ent_ctr_b ) {
//						fprintf(stderr, "comm_ent_ctr_a <  comm_ent_ctr_b\n" );

						const size_t &common_entry_a = common_entries[comm_ent_ctr_a];
						const size_t &common_entry_b = common_entries[comm_ent_ctr_b];

//						cerr << "\t\tEntry\t" << common_entry_a << " and\t" << common_entry_b << endl;

						const protein       &protein_a   = arg_proteins[ common_entry_a ];
						const protein       &protein_b   = arg_proteins[ common_entry_b ];

						const aln_posn_type  posn_a_from = get_position_of_entry_of_index( arg_alignment, common_entry_a, from_index );
						const aln_posn_type  posn_b_from = get_position_of_entry_of_index( arg_alignment, common_entry_b, from_index );

						const aln_posn_type  posn_a_to   = get_position_of_entry_of_index( arg_alignment, common_entry_a, to_index   );
						const aln_posn_type  posn_b_to   = get_position_of_entry_of_index( arg_alignment, common_entry_b, to_index   );

						const residue       &res_a_from  = protein_a.get_residue_ref_of_index( posn_a_from );
						const residue       &res_b_from  = protein_b.get_residue_ref_of_index( posn_b_from );
						const residue       &res_a_to    = protein_a.get_residue_ref_of_index( posn_a_to   );
						const residue       &res_b_to    = protein_b.get_residue_ref_of_index( posn_b_to   );

						const float_score_type max_score = residue_querier::RESIDUE_A_VALUE / residue_querier::RESIDUE_B_VALUE;
						const float_score_type score     = context_res(
							res_a_from,
							res_b_from,
							res_a_to,
							res_b_to
						);
						numerators  [ common_entry_a ][ from_index ] += score;
						numerators  [ common_entry_b ][ from_index ] += score;
						denominators[ common_entry_a ][ from_index ] += max_score;
						denominators[ common_entry_b ][ from_index ] += max_score;

						numerators  [ common_entry_a ][ to_index   ] += score;
						numerators  [ common_entry_b ][ to_index   ] += score;
						denominators[ common_entry_a ][ to_index   ] += max_score;
						denominators[ common_entry_b ][ to_index   ] += max_score;
					}
				}
			}
		}
	}

	opt_score_vec_vec scores( num_entries, opt_score_vec( length ) );
	for (size_t index_ctr = 0; index_ctr < length; ++index_ctr) {
		for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
			if ( has_position_of_entry_of_index( arg_alignment, entry_ctr, index_ctr ) ) {
				const float_score_type &numerator   = numerators  [ entry_ctr ][ index_ctr ];
				const float_score_type &denominator = denominators[ entry_ctr ][ index_ctr ];
				scores[ entry_ctr ][ index_ctr ] = ( denominator != 0.0 ) ? ( numerator / denominator ) : 0.0;
			}
		}
	}

	return make_alignment_residue_scores( arg_alignment, scores );
}


/// \brief TODOCUMENT
///
/// \relates residue_scorer
void cath::align::score_alignment(const residue_scorer &arg_residue_scorer, ///< TODOCUMENT
                                  alignment            &arg_alignment,      ///< TODOCUMENT
                                  const protein_list   &arg_proteins        ///< TODOCUMENT
                                  ) {
	const alignment_residue_scores the_scores = arg_residue_scorer.get_alignment_residue_scores(
		arg_alignment,
		arg_proteins
	);
	arg_alignment.set_scores( the_scores );
}

/// \brief TODOCUMENT
///
/// \relates residue_scorer
alignment cath::align::score_alignment_copy(const residue_scorer &arg_residue_scorer, ///< TODOCUMENT
                                            alignment             arg_alignment,      ///< TODOCUMENT
                                            const protein_list   &arg_proteins        ///< TODOCUMENT
                                            ) {
	score_alignment( arg_residue_scorer, arg_alignment, arg_proteins );
	return arg_alignment;
}

/// \brief TODOCUMENT
///
/// \relates residue_scorer
alignment cath::align::read_and_rescore_fasta_alignment(const path           &arg_fasta_aln_file, ///< TODOCUMENT
                                                        const protein_list   &arg_proteins,       ///< TODOCUMENT
                                                        const residue_scorer &arg_residue_scorer, ///< TODOCUMENT
                                                        ostream              &arg_stderr          ///< TODOCUMENT
                                                        ) {
	return score_alignment_copy(
		arg_residue_scorer,
		read_alignment_from_fasta_file(
			arg_fasta_aln_file,
			arg_proteins,
			arg_stderr
		),
		arg_proteins
	);
}
