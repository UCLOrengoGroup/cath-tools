/// \file
/// \brief The view_cache_index class definitions

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

#include "view_cache_index.h"

#include <boost/filesystem/path.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.h"
#include "alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.h"
#include "alignment/dyn_prog_align/detail/matrix_plotter/matrix_plot.h"
#include "alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.h"
#include "alignment/dyn_prog_align/std_dyn_prog_aligner.h"
#include "alignment/gap/gap_penalty.h"
#include "alignment/io/alignment_io.h"
#include "alignment/pair_alignment.h"
#include "alignment/refiner/alignment_refiner.h"
#include "common/difference.h"
#include "common/file/simple_file_read_write.h"
#include "ssap/context_res.h"
#include "structure/geometry/coord.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "structure/view_cache/index/detail/vcie_match_criteria.h"
#include "structure/view_cache/index/quad_find_action.h"
#include "structure/view_cache/index/quad_find_action_check.h"

#include <chrono>
#include <iostream>

using namespace boost::filesystem;
using namespace boost::math::constants;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::index;
using namespace cath::index::detail;
using namespace std;
using namespace std::chrono;

using boost::lexical_cast;
using boost::numeric_cast;

/// \brief Ctor for view_cache_index
view_cache_index::view_cache_index(const dim_tuple &arg_dim_defaults ///< TODOCUMENT
                                   ) : dim_defaults( arg_dim_defaults        ),
                                       the_index   ( dim_defaults.get_head() ) {
}

/// \brief TODOCUMENT
void view_cache_index::store(const view_cache_index_entry &arg_entry ///< TODOCUMENT
                             ) {
	the_index.store(
		arg_entry,
		dim_defaults.get_tail()
	);
}

///// \brief TODOCUMENT
/////
///// \relates view_cache_index
//bool cath::index_contains_value(const view_cache_index &arg_view_cache_index, ///< TODOCUMENT
//                                const coord            &arg_coord             ///< TODOCUMENT
//                                ) {
//	const opt_vcie_vec_cref entries = arg_view_cache_index.get_entries( arg_coord );
//	/// \todo Make this neater with ranges/lambdas/algorithms etc
//	if ( entries ) {
//		for (const view_cache_index_entry &the_entry : *entries) {
//			if ( the_entry.get_view() == arg_coord ) {
//				return true;
//			}
//		}
//	}
//	return false;
//}

/// \brief TODOCUMENT
///	
/// \relates view_cache_index
view_cache_index cath::index::build_view_cache_index(const double              &arg_xyz_cell_width,       ///< TODOCUMENT
                                                     const angle_type          &arg_phi_angle_cell_width, ///< TODOCUMENT
                                                     const angle_type          &/*arg_psi_angle_cell_width*/, ///< TODOCUMENT
                                                     const protein             &arg_protein,              ///< TODOCUMENT
                                                     const vcie_match_criteria &arg_criteria              ///< TODOCUMENT
                                                     ) {
	view_cache_index new_view_cache_index( boost::make_tuple(
		view_cache_index_dim_dirn           (                          ),
		view_cache_index_dim_linear_from_phi( arg_phi_angle_cell_width ),
		// view_cache_index_dim_linear_from_psi( arg_psi_angle_cell_width ),
		// view_cache_index_dim_linear_to_phi  ( arg_phi_angle_cell_width ),
		// view_cache_index_dim_linear_to_psi  ( arg_psi_angle_cell_width ),
		view_cache_index_dim_linear_x       ( numeric_cast<view_base_type>( arg_xyz_cell_width ) ),
		view_cache_index_dim_linear_y       ( numeric_cast<view_base_type>( arg_xyz_cell_width ) ),
		view_cache_index_dim_linear_z       ( numeric_cast<view_base_type>( arg_xyz_cell_width ) ),
		view_cache_index_tail               (                          )
	) );
	const size_t num_residues = arg_protein.get_length();
	for (size_t from_ctr = 0; from_ctr < num_residues; ++from_ctr) {
		for (size_t to_ctr = 0; to_ctr < num_residues; ++to_ctr) {
			const view_cache_index_entry the_entry = make_view_cache_index_entry( arg_protein, from_ctr, to_ctr );
			if ( arg_criteria( the_entry ) ) {
				new_view_cache_index.store( the_entry );
			}
		}
	}
//	cerr << "At finish, contains_it? is : " << boolalpha << index_contains_value( new_view_cache_index, coord_of_interest );
	return new_view_cache_index;
}

/// \brief TODOCUMENT
///
/// \relates view_cache_index
high_resolution_clock::duration cath::index::process_quads_indexed(const protein             &arg_protein_a,           ///< TODOCUMENT
                                                                   const protein             &arg_protein_b,           ///< TODOCUMENT
                                                                   const double              &arg_xyz_cell_size,       ///< TODOCUMENT
                                                                   const angle_type          &arg_phi_angle_cell_size, ///< TODOCUMENT
                                                                   const angle_type          &arg_psi_angle_cell_size, ///< TODOCUMENT
                                                                   const vcie_match_criteria &arg_criteria,            ///< TODOCUMENT
                                                                   quad_find_action_check    &arg_action               ///< TODOCUMENT
                                                                   ) {
	const protein_list proteins = make_protein_list( { arg_protein_a, arg_protein_b } );

//	const auto build_start_time   = high_resolution_clock::now();

	const view_cache_index view_cache_index_a = build_view_cache_index(
		arg_xyz_cell_size,
		arg_phi_angle_cell_size,
		arg_psi_angle_cell_size,
		arg_protein_a,
		arg_criteria
	);
	const view_cache_index view_cache_index_b = build_view_cache_index(
		arg_xyz_cell_size,
		arg_phi_angle_cell_size,
		arg_psi_angle_cell_size,
		arg_protein_b,
		arg_criteria
	);

	const auto scan_start_time   = high_resolution_clock::now();

	// cerr << "Before leaves" << endl;
	// view_cache_index_a.perform_action_on_all_match_at_leaves( view_cache_index_b, arg_criteria, arg_action );
	// cerr << "After leaves" << endl;

	cerr << "Before nodes" << endl;
	view_cache_index_a.perform_action_on_all_match_at_nodes ( view_cache_index_b, arg_criteria, arg_action );
	cerr << "After nodes" << endl;
	
	// cerr << "Error: currently doing both of perform_action_on_all_match_at_leaves and perform_action_on_all_match_at_nodes" << endl;

	const auto scan_complete_time = high_resolution_clock::now();

	return scan_complete_time - scan_start_time;
}



//										const bool in_alignment = (
//											alignment_contains_pair( arg_alignment, entry.get_from_index(), match.get_from_index() )
//											&&
//											alignment_contains_pair( arg_alignment, entry.get_to_index(),   match.get_to_index()   )
//										);
//
//										total_scores[ vcie_a.get_from_index() ][ vcie_b.get_from_index() ] += score;
//										total_scores[ vcie_a.get_to_index()   ][ vcie_b.get_to_index()   ] += score;

//	const score_alignment_pair score_and_alignment = std_dyn_prog_aligner().align(
//		new_matrix_dyn_prog_score_source( total_scores, num_entries_a, num_entries_b ),
//		the_gap_penalty,
//		get_window_width_for_full_matrix( num_entries_a, num_entries_b )
//	);
//
//	const auto                align_complete_time  = high_resolution_clock::now();
//
//	const alignment &scan_alignment    = score_and_alignment.second;
//	const alignment  refined_alignment = alignment_refiner().iterate(
//		scan_alignment,
//		make_protein_list( { arg_protein_a, arg_protein_b } ),
//		the_gap_penalty
//	);
//
//	const auto                refine_complete_time  = high_resolution_clock::now();

//	const duration build_durn  = scan_start_time      - build_start_time;
//	const duration scan_durn   = scan_complete_time   - scan_start_time;
//	const duration align_durn  = align_complete_time  - scan_complete_time;
//	const duration refine_durn = refine_complete_time - align_complete_time;

//	const double build_durn_in_seconds  = numeric_cast<double>( build_durn.total_nanoseconds()  ) / 1000000000.0;
//	const double scan_durn_in_seconds   = numeric_cast<double>( scan_durn.total_nanoseconds()   ) / 1000000000.0;
//	const double align_durn_in_seconds  = numeric_cast<double>( align_durn.total_nanoseconds()  ) / 1000000000.0;
//	const double refine_durn_in_seconds = numeric_cast<double>( refine_durn.total_nanoseconds() ) / 1000000000.0;
//	cerr << "build_durn : "    << right << setw( 9 ) <<       build_durn_in_seconds  << "s";
//	cerr << " ("               << right << setw( 9 ) << 1.0 / build_durn_in_seconds  << "/s)";
//	cerr << "scan_durn : "     << right << setw( 9 ) <<       scan_durn_in_seconds   << "s";
//	cerr << " ("               << right << setw( 9 ) << 1.0 / scan_durn_in_seconds   << "/s)";
//	cerr << ", align_durn : "  << right << setw( 9 ) <<       align_durn_in_seconds  << "s";
//	cerr << " ("               << right << setw( 9 ) << 1.0 / align_durn_in_seconds  << "/s)";
//	cerr << ", refine_durn : " << right << setw( 9 ) <<       refine_durn_in_seconds << "s";
//	cerr << " ("               << right << setw( 9 ) << 1.0 / refine_durn_in_seconds << "/s)";
//	cerr << "\n";
//
//	using doub_doub_doub_pair_pair = pair<double, doub_doub_pair>;
//	for (const doub_doub_doub_pair_pair &histogram_entry : histogram_data) {
//		cerr        << histogram_entry.first;
//		cerr << " " << histogram_entry.second.first;
//		cerr << " " << histogram_entry.second.second;
//		cerr << " " << numeric_cast<double>( histogram_entry.second.first  ) / numeric_cast<double>( num_in_alignment );
//		cerr << " " << numeric_cast<double>( histogram_entry.second.second ) / numeric_cast<double>( num_not_in_alignment );
//
//		cerr << "\n";
//	}
//
//	write_alignment_as_fasta_alignment( "refined_alignment.fa", refined_alignment, proteins );
//
////	cerr << scan_complete_time
////	cerr << total_score << endl;
////	cerr      << "Found " << right << setw( 9) << num_matches << ", which took " << search_durn.total_microseconds() << " microseconds" << endl;
////	ostringstream dummy_oss;
////	dummy_oss << total_score << "\n";
////	dummy_oss << "Found " << right << setw( 9) << num_matches << ", which took " << search_durn.total_microseconds() << " microseconds\n";

//	const new_matrix_dyn_prog_score_source scorer( total_scores, num_entries_a, num_entries_b );
//	matrix_plot<gnuplot_matrix_plotter>( "view_cache_index_matrix", scorer );



/// \brief TODOCUMENT
///
/// \relates view_cache_index
high_resolution_clock::duration cath::index::process_quads_complete(const protein             &arg_protein_a, ///< TODOCUMENT
                                                                    const protein             &arg_protein_b, ///< TODOCUMENT
                                                                    const double              &/*arg_cell_size*/, ///< TODOCUMENT
                                                                    const vcie_match_criteria &arg_criteria,  ///< TODOCUMENT
                                                                    quad_find_action          &arg_action     ///< TODOCUMENT
                                                                    ) {
	const auto num_entries_a   = arg_protein_a.get_length();
	const auto num_entries_b   = arg_protein_b.get_length();
	float_score_vec_vec full_scores_all( num_entries_a, float_score_vec( num_entries_b, 0 ) );

	const auto scan_start_time = high_resolution_clock::now();

	for (size_t from_a = 0; from_a < num_entries_a; ++from_a) {
		for (size_t to_a = 0; to_a < num_entries_a; ++to_a) {
			if ( from_a != to_a ) {
				const size_size_pair indices_a( from_a, to_a );
				for (size_t from_b = 0; from_b < num_entries_b; ++from_b) {
					for (size_t to_b = 0; to_b < num_entries_b; ++to_b) {
						if ( from_b != to_b ) {
							const size_size_pair indices_b( from_b, to_b );
							if ( arg_criteria( indices_a, indices_b, arg_protein_a, arg_protein_b ) ) {

//								const double sq_dist = squared_distance( indices_a, indices_b, arg_protein_a, arg_protein_b );
//								const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
//
//								std::cerr << "compl_match ";
//								std::cerr << std::right << std::setw( 4 ) << indices_a.first;
//								std::cerr << "";
//								std::cerr << std::right << std::setw( 4 ) << indices_a.second;
//								std::cerr << " ";
//								std::cerr << std::right << std::setw( 4 ) << indices_b.first;
//								std::cerr << " ";
//								std::cerr << std::right << std::setw( 4 ) << indices_b.second;
//								std::cerr << " " << score;
//								std::cerr << std::endl;


								arg_action( indices_a, indices_b );
							}
						}
					}
				}
			}
		}
	}
	const auto scan_complete_time = high_resolution_clock::now();

//	const new_matrix_dyn_prog_score_source scorer_all( full_scores_all, num_entries_a, num_entries_b );
//	matrix_plot<gnuplot_matrix_plotter>( "full_scores_all", scorer_all );

	return scan_complete_time - scan_start_time;
}



//	using view_comp_set = tuple<double, double, double, double, double, double>;
//	using bool_bool_pair = pair<bool, bool>;
//	using view_comp_record = pair<view_comp_set, bool_bool_pair>;
//	using view_comp_record_vec = vector<view_comp_record>;
//
//  view_comp_record_vec records;
//
//	// In the loop...
// const bool             index_dirn_diff  = ! same_direction          ( indices_a, indices_b );
// const size_t           min_index_diff   = min_index_difference      ( indices_a, indices_b );
// const double           sq_dist          = squared_distance          ( indices_a, indices_b, arg_protein_a, arg_protein_b );
// const angle            frame_diff_angle = angle_between_frame_angles( indices_a, indices_b, arg_protein_a, arg_protein_b );
// const angle            max_phi_diff     = max_phi_angle_difference  ( indices_a, indices_b, arg_protein_a, arg_protein_b );
// const angle            max_psi_diff     = max_psi_angle_difference  ( indices_a, indices_b, arg_protein_a, arg_protein_b );
//	records.push_back(
//		make_pair(
//			make_tuple(
//				sq_dist,
//				( index_dirn_diff ? 1.0 : 0.0 ),
//				0.0 - numeric_cast<double>( min_index_diff ),
//				angle_in_degrees( frame_diff_angle ),
//				angle_in_degrees( max_phi_diff ),
//				angle_in_degrees( max_psi_diff )
//			),
//			make_pair(
//				alignment_contains_pair( arg_alignment, from_a, from_b ),
//				alignment_contains_pair( arg_alignment, to_a,   to_b   )
//			)
//		)
//	);
////	sq_dist, index_dirn_diff, min_index_distance, frame_diff_angle, min_phi_diff, min_psi_diff
//
//	str_vec num_output_strings;
//	str_vec val_output_strings;
////	for (double dist_co = 2.0; dist_co <= 6.5; dist_co += 2.0) {
//	for (double dist_co = 4.0; dist_co <= 40.0; dist_co += 4.0) {
////		for (double dirn_co = 1.0; dirn_co <= 1.0; dirn_co += 1.0) {
//		for (double dirn_co = 0.0; dirn_co < 1.5; dirn_co += 1.0) {
////			for (double index_dist_co = -1.0; index_dist_co <= -1.0; index_dist_co += 1.0) {
//			for (double index_dist_co = -31.0; index_dist_co <= -1.0; index_dist_co += 5.0) {
////				for (double frame_ang_co = 180.0; frame_ang_co < 185.0; frame_ang_co += 22.5) {
//				for (double frame_ang_co = 22.5; frame_ang_co < 185.0; frame_ang_co += 22.5) {
////					for (double phi_ang_co = 180.0; phi_ang_co < 185.0; phi_ang_co += 22.5) {
//					for (double phi_ang_co = 22.5; phi_ang_co < 185.0; phi_ang_co += 22.5) {
//						cerr << "dist_co : "         << dist_co;
//						cerr << ", dirn_co : "       << dirn_co;
//						cerr << ", index_dist_co : " << index_dist_co;
//						cerr << ", frame_ang_co : "  << frame_ang_co;
//						cerr << ", phi_ang_co : "    << phi_ang_co;
//						cerr << endl;
////						for (double psi_ang_co = 180.0; psi_ang_co < 185.0; psi_ang_co += 22.5) {
//						for (double psi_ang_co = 22.5; psi_ang_co < 185.0; psi_ang_co += 22.5) {
//							size_t true_pos_num = 0;
//							size_t fals_pos_num = 0;
//							size_t fals_neg_num = 0;
////							size_t true_neg_num = 0;
//
//							double true_pos_val = 0;
//							double fals_pos_val = 0;
//							double fals_neg_val = 0;
////							double true_neg_val = 0;
//							for (const view_comp_record &record : records) {
//								const double &dist         = record.first.get<0>();
//								const double &dirn         = record.first.get<1>();
//								const double &index_dist   = record.first.get<2>();
//								const double &frame_ang    = record.first.get<3>();
//								const double &phi_ang      = record.first.get<4>();
//								const double &psi_ang      = record.first.get<5>();
//								const bool   &from_correct = record.second.first;
//								const bool   &to_correct   = record.second.second;
//
//								const double score = 10.0 / (dist + 10.0);
//
//								const bool    positive = ( dist <= dist_co && dirn <= dirn_co && index_dist <= index_dist_co && frame_ang <= frame_ang_co && phi_ang <= phi_ang_co && psi_ang <= psi_ang_co );
//								const bool    answer   = ( from_correct && to_correct );
//
////								if ( ! positive ) {
////									cerr << "Oooooooo!!!" << endl;
////									cerr << "dist_co : "         << dist_co;
////									cerr << ", dirn_co : "       << dirn_co;
////									cerr << ", index_dist_co : " << index_dist_co;
////									cerr << ", frame_ang_co : "  << frame_ang_co;
////									cerr << ", phi_ang_co : "    << phi_ang_co;
////									cerr << ", psi_ang_co : "    << phi_ang_co;
////									cerr << endl;
////
////									cerr << "dist : "         << dist;
////									cerr << ", dirn : "       << dirn;
////									cerr << ", index_dist : " << index_dist;
////									cerr << ", frame_ang : "  << frame_ang;
////									cerr << ", phi_ang : "    << phi_ang;
////									cerr << ", psi_ang : "    << phi_ang;
////									cerr << endl;
////								}
//
//								if ( positive ) {
//									if ( answer ) {
//										++true_pos_num;
//										true_pos_val += score;
//									}
//									else {
//										++fals_pos_num;
//										fals_pos_val += score;
//									}
//								}
//								else {
//									if ( answer ) {
//										++fals_neg_num;
//										fals_neg_val += score;
//									}
//									else {
////										++true_neg_num;
////										true_neg_val += score;
//									}
//								}
//							}
//
//							const double num_sensitivity = 100.0 * numeric_cast<double>( true_pos_num ) / numeric_cast<double>( true_pos_num + fals_neg_num );
//							const double num_precision   = 100.0 * numeric_cast<double>( true_pos_num ) / numeric_cast<double>( true_pos_num + fals_pos_num );
//							const double val_sensitivity = 100.0 *                       true_pos_val   /                     ( true_pos_val + fals_neg_val );
//							const double val_precision   = 100.0 *                       true_pos_val   /                     ( true_pos_val + fals_pos_val );
////							cerr << "dist_co : "         << dist_co;
////							cerr << ", dirn_co : "       << dirn_co;
////							cerr << ", index_dist_co : " << index_dist_co;
////							cerr << ", frame_ang_co : "  << frame_ang_co;
////							cerr << ", phi_ang_co : "    << phi_ang_co;
//							num_output_strings.push_back(
//								        to_string( num_sensitivity )
//								+ " " + to_string( num_precision   )
//								+ " # dist_co: "     + to_string( dist_co       )
//								+ " dirn_co: "       + to_string( dirn_co       )
//								+ " index_dist_co: " + to_string( index_dist_co )
//								+ " frame_ang_co: "  + to_string( frame_ang_co  )
//								+ " phi_ang_co: "    + to_string( phi_ang_co    )
//								+ " psi_ang_co: "    + to_string( psi_ang_co    )
//							);
//
//							val_output_strings.push_back(
//								        to_string( val_sensitivity )
//								+ " " + to_string( val_precision   )
//								+ " # dist_co: "     + to_string( dist_co       )
//								+ " dirn_co: "       + to_string( dirn_co       )
//								+ " index_dist_co: " + to_string( index_dist_co )
//								+ " frame_ang_co: "  + to_string( frame_ang_co  )
//								+ " phi_ang_co: "    + to_string( phi_ang_co    )
//								+ " psi_ang_co: "    + to_string( psi_ang_co    )
//							);
//						}
//					}
//				}
//			}
//		}
//	}
//
//	write_file( path("residue_sensitivity_precision.num.data"), num_output_strings );
//	write_file( path("residue_sensitivity_precision.val.data"), val_output_strings );


