/// \file
/// \brief The indexed_refiner class definitions

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

#include "indexed_refiner.hpp"

#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/replace.hpp>
// #include <boost/container/small_vector.hpp> // ***** small_vector was only added in Boost 1.58.0 *****
#include <boost/core/demangle.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <gnuplot-iostream.h>

#include "common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/chrono/duration_to_seconds_string.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/metaprogramming/combine_params_lists_with_template_list.hpp"
#include "common/metaprogramming/template_list.hpp"
#include "common/type_to_string.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::scan;

#include <tuple>
#include <type_traits>
#include <vector>

using boost::adaptors::transformed;
using boost::algorithm::erase_all;
using boost::algorithm::replace_all;
using boost::algorithm::replace_all_copy;
// using boost::container::small_vector;
using boost::filesystem::path;
using boost::irange;
using std::chrono::high_resolution_clock;
using std::pair;
using std::string;
using std::tuple;
using std::vector;

using cell_type_list = tuple<
                              tuple< vector      <simple_locn_index>     >
                              // tuple< small_vector<simple_locn_index,  1> >,
                              // tuple< small_vector<simple_locn_index,  2> >,
                              // tuple< small_vector<simple_locn_index,  4> >,
                              // tuple< small_vector<simple_locn_index, 10> >,
                              // tuple< small_vector<simple_locn_index, 20> >,
                              // tuple< small_vector<simple_locn_index, 40> >
                              >;

using index_template_list = template_list< 
                                           // vector_refine_index,
                                           vec_of_vectors_refine_index,
                                           // hash_refine_index,
                                           vec_of_hashes_refine_index,
                                           // lattice_refine_index,
                                           vec_of_lattices_refine_index
                                           >;

using combined_list = combine_params_lists_with_template_list_t< cell_type_list, index_template_list >;

constexpr float indexed_refiner_constants::MAX_DIST;
constexpr float indexed_refiner_constants::MAX_SQ_DIST;

namespace {

	inline void simplify_name(std::string &arg) {
		erase_all  ( arg, "cath::align::"            );
		erase_all  ( arg, "cath::scan::"             );
		erase_all  ( arg, "boost::container::"       );
		erase_all  ( arg, "std::"                    );
		erase_all  ( arg, "_refine_index"            );
		replace_all( arg, "simple_locn_index", "SLI" );
		erase_all  ( arg, ", new_allocator<SLI> "    );
		erase_all  ( arg, ", allocator<SLI>"         );
		replace_all( arg, "ul>>", ">>"               );
		erase_all  ( arg, " "                        );

		// erase_all  ( arg, ", allocator<simple_locn_index>" );
	}

	inline std::string simplify_name_copy(std::string arg
	                                      ) {
		simplify_name( arg );
		return arg;
	}

	template <typename Fn, typename... Ts>
	hrc_duration time_execution(Fn &&prm_fn,
	                            Ts &&...prm_vars
	                            ) {
		const auto start_time = high_resolution_clock::now();
		invoke( prm_fn, std::forward< Ts >( prm_vars )... );
		return high_resolution_clock::now() - start_time;
	}

	template <typename T>
	class index_type_pair_processor final {};

	template <typename... Ts>
	class index_type_pair_processor< tuple< Ts... > > final {
	private:
		template <sod Sod,
		          typename T,
		          typename U>
		void process_index_type_pair_sod(vector<pair<string, doub_doub_pair_vec>> &prm_data,
		                                 const protein                            &prm_protein_a,    ///< TODOCUMENT
		                                 const protein                            &prm_protein_b,    ///< TODOCUMENT
		                                 const alignment                          &prm_alignment ///< TODOCUMENT
		                                 ) {
			// std::cerr << "Emplacing back " << demangle( typeid( T ).name() ) << " and " << demangle( typeid( U ).name() ) << "\n";

			// const auto keyer = make_res_pair_keyer(
			// 	simple_locn_x_keyer_part{ prm_cell_size },
			// 	simple_locn_y_keyer_part{ prm_cell_size },
			// 	simple_locn_z_keyer_part{ prm_cell_size }
			// );

			// return make_pair(
			// 	store_maker<sod::DENSE, T>{}(
			// 	prm_rng,
			// 	keyer,
			// 	simple_locn_crit{ prm_max_dist * prm_max_dist }
			// );

			// store_maker<sod::DENSE, T>{}(

			// );

			doub_doub_pair_vec series_data;

			std::cerr << ( to_string( Sod )
					+ ":"
					+ simplify_name_copy( common::type_to_string<T>() ) ) << std::endl;

			// ** \todo ENSURE THE READING VARIES WITH SOD **

			constexpr size_t start_offset = ( Sod == sod::DENSE ? 3_z : 3_z );

			for (const size_t &cell_size_index : irange( start_offset, 12_z ) ) {
			// for (const size_t &cell_size_index : irange( 5_z, 6_z ) ) {

				const float cell_size = static_cast<float>( cell_size_index ) * 4.0f;
				std::cerr << "\t" << cell_size << std::endl;
				info_quantity memory_usage;
				const auto build_durn = time_execution(
					[&] {
						const T store_wrapper_a{
							std::integral_constant< sod, Sod >{},
							prm_protein_a,
							fot::FROM,
							cell_size
						};
						const T store_wrapper_b{
							std::integral_constant< sod, Sod >{},
							prm_protein_b,
							fot::FROM,
							cell_size
						};
						memory_usage = store_wrapper_a.get_info_size() + store_wrapper_b.get_info_size();



						for (const auto &aln_index : indices( prm_alignment.length() ) ) {
							if ( has_both_positions_of_index( prm_alignment, aln_index ) ) {
								const auto &a_posn = get_a_position_of_index( prm_alignment, aln_index );
								const auto &b_posn = get_b_position_of_index( prm_alignment, aln_index );
								const auto &cells_a = store_wrapper_a.the_store[ a_posn ];
								const auto &cells_b = store_wrapper_b.the_store[ b_posn ];

								for (const auto &the_cell_a : cells_a ) {
									for (const simple_locn_index &the_entry_a : the_cell_a.second ) {
										for (const simple_locn_index &candidate_b : cells_b.find_matches( store_wrapper_a.the_keyer.make_key( the_entry_a ) ) ) {
											if ( get_squared_distance( candidate_b, the_entry_a ) < indexed_refiner_constants::MAX_SQ_DIST ) {
												std::cerr << to_string( the_entry_a ) << ", " << to_string( candidate_b ) << "\n";
											}
										}
									}
								}
							}
						}
					}
				);

				// std::cerr << series_name << "\t" << cell_size << "\t" << durn_to_seconds_string( build_durn ) << "\n";

				series_data.emplace_back(
					cell_size,
					durn_to_seconds_double( build_durn )
				);
				// "size "s + std::to_string( memory_usage.value() ) + "b"
				// "build durn : " + durn_to_seconds_string( build_durn ) + " size "s + boost::lexical_cast<string>( memory_usage.value() ) + "b"
			}

			prm_data.emplace_back(
				to_string( Sod )
					+ ":"
					+ simplify_name_copy( common::type_to_string<T>() ),
				series_data
			);
		}

		template <typename T,
		          typename U>
		int process_index_type_pair(vector<pair<string, doub_doub_pair_vec>> &prm_data,      ///< TODOCUMENT
		                            const protein                            &prm_protein_a, ///< TODOCUMENT
		                            const protein                            &prm_protein_b, ///< TODOCUMENT
		                            const alignment                          &prm_alignment  ///< TODOCUMENT
		                            ) {
			process_index_type_pair_sod< sod::SPARSE, T, U >( prm_data, prm_protein_a, prm_protein_b, prm_alignment );
			process_index_type_pair_sod< sod::DENSE,  T, U >( prm_data, prm_protein_a, prm_protein_b, prm_alignment );
			return 0;
		}

	public:
		/// \brief TODOCUMENT
		vector<pair<string, doub_doub_pair_vec>> operator()(const protein   &prm_protein_a, ///< TODOCUMENT
		                                                    const protein   &prm_protein_b, ///< TODOCUMENT
		                                                    const alignment &prm_alignment  ///< TODOCUMENT
		                                                    ) {
			vector<pair<string, doub_doub_pair_vec>> data;
			const std::initializer_list<int> dummy_list = { process_index_type_pair< Ts, Ts >( data, prm_protein_a, prm_protein_b, prm_alignment )... };
			boost::ignore_unused( dummy_list );
			return data;
		}
	};
} // namespace

void cath::align::do_some_gubbins(const protein   &prm_protein_a, ///< TODOCUMENT
                                  const protein   &prm_protein_b, ///< TODOCUMENT
                                  const alignment &prm_alignment  ///< TODOCUMENT
                                  ) {
	// std::cerr << "Doing some gubbins" << std::endl;
	const auto results = index_type_pair_processor< combined_list >{}( prm_protein_a, prm_protein_b, prm_alignment );

	const path base_filename    = "/tmp/indexed_refiner_graph.eps";

	const path gnuplot_file     = replace_extension_copy( base_filename, ".gnuplot"         );
	const path eps_file         = replace_extension_copy( base_filename, ".eps"             );
	const path the_data_file    = replace_extension_copy( base_filename, ".data.txt"        );
	const path filter_data_file = replace_extension_copy( base_filename, ".filter_data.txt" );
	Gnuplot gp("tee " + gnuplot_file.string() + " | gnuplot"); // Write to an intermediate gnuplot file

	gp << "set   terminal postscript color\n";
	gp << "set   output " << eps_file << "\n";
//	gp << "set   size square\n";

//	gp << "set   xtics 0,10\n";
//	gp << "set   ytics 0,10\n";
	gp << "set   xtics font \"Helvetica,10\"\n";
	gp << "set   ytics font \"Helvetica,10\"\n";
	gp << "set   key font \"Helvetica,7\"\n";
	gp << "set   style line 11 lc rgb '#808080' lt 1\n";
	gp << "set   border 3 back ls 11\n";
	gp << "set   tics nomirror\n";
	gp << "set   style line 12 lc rgb '#808080' lt 0 lw 1\n";
	gp << "set   style line 1 lt rgb \"#000000\"\n";
	gp << "set   grid back ls 12\n";

	gp << "set   title \"Time to build a full index\"\n";
	gp << "set   xlabel \"Cell size (angstroms)\"\n";
	gp << "set   ylabel \"Time (seconds)\"\n";

	gp
		<< "plot "
		<< boost::algorithm::join(
			results
				| transformed( [&] (const pair<string, doub_doub_pair_vec> &x) {
					return "'-' with linespoints title '" + replace_all_copy( x.first, "_", "\\_" ) + "' ";
				} ),
			", "
		) << "\n";

	for (const auto &x : results) {
		gp.send1d( x.second );
	}
}


// /// \brief TODOCUMENT
// bool_aln_pair indexed_refiner::iterate_step(const alignment       &prm_alignment,       ///< TODOCUMENT
//                                             const protein_list    &prm_proteins,        ///< TODOCUMENT
//                                             // const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
//                                             const gap_penalty     &prm_gap_penalty      ///< TODOCUMENT
//                                             ) {
// 	cerr << "Will try to find ways to split alignment with " << prm_alignment.num_entries() << " entries" << endl;

// 	const size_t num_entries = prm_alignment.num_entries();
// 	if ( prm_proteins.size() != num_entries ) {
// 		BOOST_THROW_EXCEPTION(not_implemented_exception("Mismatch between number of entries in alignment and in protein list"));
// 	}

// 	return iterate_step_for_alignment_split_list(
// 		prm_alignment,
// 		prm_proteins,
// 		// prm_view_cache_list,
// 		prm_gap_penalty,
// 		get_standard_alignment_splits( prm_alignment )
// 	);
// }

// /// \brief TODOCUMENT
// bool_aln_pair indexed_refiner::iterate_step_for_alignment_split_list(const alignment            &prm_alignment,           ///< TODOCUMENT
//                                                                      const protein_list         &prm_proteins,            ///< TODOCUMENT
//                                                                      // const view_cache_list      &prm_view_cache_list,     ///< TODOCUMENT
//                                                                      const gap_penalty          &prm_gap_penalty,         ///< TODOCUMENT
//                                                                      const alignment_split_list &prm_alignment_split_list ///< TODOCUMENT
//                                                                      ) {
// 	alignment iter_aln( prm_alignment );
// 	bool inserted_residues = false;
// 	for (const alignment_split &the_split : prm_alignment_split_list) {
// 		const bool_aln_pair inserted_res_and_aln = iterate_step_for_alignment_split(
// 			iter_aln,
// 			prm_proteins,
// 			// prm_view_cache_list,
// 			prm_gap_penalty,
// 			the_split
// 		);
// 		inserted_residues = ( inserted_res_and_aln.first || inserted_residues );
// 		iter_aln          =   inserted_res_and_aln.second;
// 	}
// 	return make_pair( inserted_residues, iter_aln );
// }

// /// \brief TODOCUMENT
// bool_aln_pair indexed_refiner::iterate_step_for_alignment_split(const alignment       &prm_alignment,       ///< TODOCUMENT
//                                                                 const protein_list    &prm_proteins,        ///< TODOCUMENT
//                                                                 // const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
//                                                                 const gap_penalty     &prm_gap_penalty,     ///< TODOCUMENT
//                                                                 const alignment_split &prm_alignment_split  ///< TODOCUMENT
//                                                                 ) {
// 	const size_vec correct_lengths = get_protein_lengths( prm_proteins );

// 	cerr << "Iterating alignment with " << prm_alignment.num_entries() << " entries using split :";
// 	for (const size_t &split_member : prm_alignment_split) {
// 		cerr << " " << split_member;
// 	}
// 	cerr << endl;

// 	if ( prm_alignment_split.get_num_entries() != prm_alignment.num_entries() ) {
// 		BOOST_THROW_EXCEPTION(not_implemented_exception("Number of entries in alignment split doesn't match number in alignment"));
// 	}

// 	// bob
// 	const alignment_split_mapping mapping_a = make_alignment_split_mapping( prm_alignment, prm_alignment_split, alignment_split_half::FIRST,  correct_lengths );
// 	const alignment_split_mapping mapping_b = make_alignment_split_mapping( prm_alignment, prm_alignment_split, alignment_split_half::SECOND, correct_lengths );
// 	const bool inserted_residues = ( mapping_a.inserted_entries() || mapping_b.inserted_entries() );

// 	// Grab the lengths
// 	const size_t full_length_a     = mapping_a.length();
// 	const size_t full_length_b     = mapping_b.length();
// 	const size_t full_window_width = get_window_width_for_full_matrix( full_length_a, full_length_b );

// 	from_and_to_alignment_scores.assign( full_length_a, float_score_vec( full_length_b, 0.0 ) );
// 	// to_alignment_scores.assign  ( full_length_a, float_score_vec( full_length_b, 0.0 ) );

// //	cerr << "number of entries in alignment is " << prm_alignment.num_entries() << endl;
// //	cerr << "number of entries in half a of split is " << mapping_a.num_entries() << endl;
// //	cerr << "number of entries in half b of split is " << mapping_b.num_entries() << endl;

// 	const alignment::size_type alignment_length = prm_alignment.length();
// 	for (const size_t &aln_ctr : indices( alignment_length ) ) {
// 		const size_opt mapping_index_a = mapping_a.index_of_orig_aln_index( aln_ctr );
// 		const size_opt mapping_index_b = mapping_b.index_of_orig_aln_index( aln_ctr );
// 		if ( ! mapping_index_a || ! mapping_index_b ) {
// 			continue;
// 		}
// //		cerr << "At alignment counter " << aln_ctr << endl;

// 		const size_vec present_orig_aln_entries_a = present_orig_aln_entries_of_index( mapping_a, *mapping_index_a );
// 		const size_vec present_orig_aln_entries_b = present_orig_aln_entries_of_index( mapping_b, *mapping_index_b );

// 		for (const size_t &present_orig_aln_entry_a : present_orig_aln_entries_a) {
// 			for (const size_t &present_orig_aln_entry_b : present_orig_aln_entries_b) {
// 				const size_t        present_entry_a = * mapping_a.entry_of_orig_aln_entry( present_orig_aln_entry_a );
// 				const size_t        present_entry_b = * mapping_b.entry_of_orig_aln_entry( present_orig_aln_entry_b );
// 				const aln_posn_type a_position      = get_position_of_entry_of_index( mapping_a, present_entry_a, *mapping_index_a );
// 				const aln_posn_type b_position      = get_position_of_entry_of_index( mapping_b, present_entry_b, *mapping_index_b );

// 	//			cerr << "Rescoring based on pair " << a_position << ", " << b_position << endl;
// 				const protein &protein_a = prm_proteins[ present_orig_aln_entry_a ];
// 				const protein &protein_b = prm_proteins[ present_orig_aln_entry_b ];
// 				const size_t   length_a  = protein_a.get_length();
// 				const size_t   length_b  = protein_b.get_length();
// //				const residue &residue_a = protein_a.get_residue_ref_of_index( a_position );
// //				const residue &residue_b = protein_b.get_residue_ref_of_index( b_position );

// 				for (const size_t &res_ctr_a : indices( length_a ) ) {
// 					for (const size_t &res_ctr_b : indices( length_b ) ) {
// 						if ( res_ctr_a != a_position && res_ctr_b != b_position ) {
// 							const size_t   other_mapping_index_a = mapping_a.index_of_protein_index( present_entry_a, res_ctr_a );
// 							const size_t   other_mapping_index_b = mapping_b.index_of_protein_index( present_entry_b, res_ctr_b );
// //							const residue &other_residue_a       = protein_a.get_residue_ref_of_index( res_ctr_a );
// //							const residue &other_residue_b       = protein_b.get_residue_ref_of_index( res_ctr_b );

// 							from_and_to_alignment_scores[ other_mapping_index_a ][ other_mapping_index_b ] += get_residue_context(
// 								present_orig_aln_entry_a,
// 								present_orig_aln_entry_b,
// 								a_position,
// 								b_position,
// 								res_ctr_a,
// 								res_ctr_b
// 							);
// 							from_and_to_alignment_scores[ other_mapping_index_a ][ other_mapping_index_b ] += get_residue_context(
// 								present_orig_aln_entry_a,
// 								present_orig_aln_entry_b,
// 								res_ctr_a,
// 								res_ctr_b,
// 								a_position,
// 								b_position
// 							);
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}

// 	float_score_vec_vec avg_scores( full_length_a, float_score_vec( full_length_b, 0 ) );
// 	for (const size_t &ctr_a : indices( full_length_a ) ) {
// 		for (const size_t &ctr_b : indices( full_length_b ) ) {
// 			avg_scores[ctr_a][ctr_b] = (from_alignment_scores[ctr_a][ctr_b] + to_alignment_scores[ctr_a][ctr_b]) / 2.0;
// 		}
// 	}

// //	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot_from"), new_matrix_dyn_prog_score_source( from_alignment_scores, length_a, length_b) );
// //	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot___to"), new_matrix_dyn_prog_score_source( to_alignment_scores,   length_a, length_b) );
// //	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot__avg"), new_matrix_dyn_prog_score_source( avg_scores,            length_a, length_b) );

// 	const new_matrix_dyn_prog_score_source scorer( avg_scores, full_length_a, full_length_b );

// 	const score_alignment_pair score_and_alignment = std_dyn_prog_aligner().align( scorer, prm_gap_penalty, full_window_width );

// 	alignment new_alignment = set_empty_scores_copy(
// 		build_alignment(
// 			score_and_alignment.second,
// 			mapping_a,
// 			mapping_b
// 		)
// 	);
// 	return make_pair( inserted_residues, new_alignment );
// }

// // /// \brief TODOCUMENT
// // alignment indexed_refiner::iterate(const alignment    &prm_alignment,  ///< TODOCUMENT
// //                                    const protein_list &prm_proteins,   ///< TODOCUMENT
// //                                    const gap_penalty  &prm_gap_penalty ///< TODOCUMENT
// //                                    ) {
// // 	return iterate(
// // 		prm_alignment,
// // 		prm_proteins,
// // 		view_cache_list( prm_proteins ),
// // 		prm_gap_penalty
// // 	);
// // }

// /// \brief TODOCUMENT
// alignment indexed_refiner::iterate(const alignment       &prm_alignment,       ///< TODOCUMENT
//                                    const protein_list    &prm_proteins,        ///< TODOCUMENT
//                                    // const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
//                                    const gap_penalty     &prm_gap_penalty      ///< TODOCUMENT
//                                    ) {
// //	if (prm_proteins.size() != 2 || prm_alignment.num_entries() != 2) {
// //		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to iterate alignments of more than two structures"));
// //	}

// 	/// \todo Move to using scores to prevent loops

// 	/// \todo Ensure that if using loops, a step that fills in alignment holes is always accepted

// 	size_t iter_ctr = 0;
// 	alignment prev_alignment( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT );
// 	alignment curr_alignment( prm_alignment );
// 	bool inserted_residues = true;
// 	while ( inserted_residues || curr_alignment != prev_alignment ) {
// 		const bool_aln_pair ins_res_and_next_aln = iterate_step( curr_alignment, prm_proteins, prm_gap_penalty );
// 		inserted_residues        = ins_res_and_next_aln.first;
// 		alignment next_alignment = ins_res_and_next_aln.second;
// 		const bool next_matches_prev= ( next_alignment == prev_alignment );
// 		swap( prev_alignment, curr_alignment );
// 		swap( curr_alignment, next_alignment );

// 		if (iter_ctr > 0) {
// 			cerr << "Refining alignment, step : " << iter_ctr << endl;

// //			// For debugging why 1fyvA00 vs 2rirA01 never stops
// //			const protein &protein_a         = prm_proteins[0];
// //			const protein &protein_b         = prm_proteins[1];
// //			const path temp_align_out_file( "temp_align." + lexical_cast<string>(iter_ctr) + ".txt" );
// //			ofstream temp_align_ofstream;
// //			open_ofstream( temp_align_ofstream, temp_align_out_file );
// //			output_alignment_to_cath_ssap_legacy_format( temp_align_ofstream, next_alignment, protein_a, protein_b );
// //			temp_align_ofstream.close();
// 		}

// 		if ( next_matches_prev ) {
// 			break;
// 		}

// 		++iter_ctr;
// 	}
// //	const protein &protein_a         = prm_proteins[0];
// //	const protein &protein_b         = prm_proteins[1];
// //	score_alignment( residue_scorer(), curr_alignment, prm_proteins );
// //	output_alignment_to_cath_ssap_legacy_format( cout, curr_alignment, protein_a, protein_b );

// 	return curr_alignment;
// }

