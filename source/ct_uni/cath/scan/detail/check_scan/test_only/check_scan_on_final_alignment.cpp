/// \file
/// \brief The check_scan_on_final_alignment class definitions

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

#include "check_scan_on_final_alignment.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/operators.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/difference.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/scan/detail/check_scan/test_only/alignment_scan_comparison.hpp"
#include "cath/scan/detail/check_scan/test_only/quad_and_rep_criteria_result.hpp"
#include "cath/scan/detail/check_scan/test_only/quad_criteria_result.hpp"
#include "cath/scan/detail/quad_criteria_are_met_by.hpp"
#include "cath/scan/detail/res_pair/functions/res_pair_core_functions.hpp"
#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/scan/detail/stride/rep_strider.hpp"
#include "cath/scan/detail/stride/roled_scan_stride.hpp"
#include "cath/scan/quad_criteria.hpp"
#include "cath/ssap/context_res.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"

//#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::scan;
using namespace ::cath::scan::detail;
using namespace ::std;

using ::boost::accumulate;
using ::boost::algorithm::join;
using ::boost::irange;
using ::boost::numeric_cast;
using ::boost::range::set_intersection;
using ::boost::range::set_symmetric_difference;

/// \brief TODOCUMENT
quad_criteria_result check_scan_on_final_alignment::criteria_result_of(const quad_criteria            &prm_criteria,   ///< TODOCUMENT
                                                                       const multi_struc_res_rep_pair &prm_query_pair, ///< TODOCUMENT
                                                                       const multi_struc_res_rep_pair &prm_index_pair  ///< TODOCUMENT
                                                                       ) {
	if ( ! are_not_violated_by( prm_criteria, prm_query_pair ) ) {
		return quad_criteria_result::QUERY_FAILS_SINGLE_CHECKS;
	}
	if ( ! are_not_violated_by( prm_criteria, prm_index_pair ) ) {
		return quad_criteria_result::INDEX_FAILS_SINGLE_CHECKS;
	}

	const auto &query_core = prm_query_pair.get_res_pair_core();
	const auto &index_core = prm_index_pair.get_res_pair_core();

	if ( squared_distance         ( query_core, index_core ) > prm_criteria.get_maximum_squared_distance()       ) {
		return quad_criteria_result::FAILS_VIEW_CHECK;
	}
	if ( max_phi_angle_difference ( query_core, index_core ) > prm_criteria.get_maximum_phi_angle_difference()   ) {
		return quad_criteria_result::FAILS_PHI_CHECK;
	}
	if ( max_psi_angle_difference ( query_core, index_core ) > prm_criteria.get_maximum_psi_angle_difference()   ) {
		return quad_criteria_result::FAILS_PSI_CHECK;
	}
	if ( distance_1_between_frames( query_core, index_core ) > prm_criteria.get_maximum_frame_angle_distance_1() ) {
		return quad_criteria_result::FAILS_FRAME_CHECK;
	}
	if ( ! are_met_by         ( prm_criteria, prm_query_pair, prm_index_pair ) ) {
		return quad_criteria_result::FAILS_QUAD_CHECKS;
	}
	return quad_criteria_result::PASS;
}

/// \brief TODOCUMENT
quad_criteria_result check_scan_on_final_alignment::criteria_result_of(const quad_criteria         &prm_criteria,   ///< TODOCUMENT
                                                                       const single_struc_res_pair &prm_query_pair, ///< TODOCUMENT
                                                                       const single_struc_res_pair &prm_index_pair  ///< TODOCUMENT
                                                                       ) {
	if ( ! are_not_violated_by( prm_criteria, prm_query_pair ) ) {
		return quad_criteria_result::QUERY_FAILS_SINGLE_CHECKS;
	}
	if ( ! are_not_violated_by( prm_criteria, prm_index_pair ) ) {
		return quad_criteria_result::INDEX_FAILS_SINGLE_CHECKS;
	}
	if ( ! are_met_by         ( prm_criteria, prm_query_pair, prm_index_pair ) ) {
		return quad_criteria_result::FAILS_QUAD_CHECKS;
	}
	return quad_criteria_result::PASS;
}

/// \brief TODOCUMENT
///
/// Provide functions that allow this function to be expressed slightly more tightly
quad_criteria_result check_scan_on_final_alignment::rep_quad_criteria_result_of(const protein            &prm_query_protein,    ///< TODOCUMENT
                                                                                const protein            &prm_index_protein,    ///< TODOCUMENT
                                                                                const quad_criteria      &prm_criteria,         ///< TODOCUMENT
                                                                                const scan_stride        &prm_scan_stride,      ///< TODOCUMENT
                                                                                const index_type         &prm_query_from_index, ///< TODOCUMENT
                                                                                const index_type         &prm_query_to_index,   ///< TODOCUMENT
                                                                                const index_type         &prm_index_from_index, ///< TODOCUMENT
                                                                                const index_type         &prm_index_to_index    ///< TODOCUMENT
                                                                                ) {
	const auto from_rep_indices = get_from_rep_of_indices( prm_scan_stride, prm_query_from_index, prm_index_from_index );
	const auto to_rep_indices   = get_to_rep_of_indices  ( prm_scan_stride, prm_query_to_index,   prm_index_to_index   );

	if ( ! from_rep_indices || ! to_rep_indices ) {
		return quad_criteria_result::HAS_NO_REP;
	}

	const auto query_from_rep_index = from_rep_indices->first;
	const auto query_to_rep_index   = to_rep_indices->first;
	const auto index_from_rep_index = from_rep_indices->second;
	const auto index_to_rep_index   = to_rep_indices->second;
	const auto query_res_pair       = make_multi_struc_res_rep_pair(
		prm_query_protein.get_residue_ref_of_index( get_index_of_rep_index( prm_scan_stride.get_query_from_strider(), query_from_rep_index ) ),
		prm_query_protein.get_residue_ref_of_index( get_index_of_rep_index( prm_scan_stride.get_query_to_strider  (), query_to_rep_index   ) ),
		0,
		query_from_rep_index,
		query_to_rep_index
	);
	const auto index_res_pair       = make_multi_struc_res_rep_pair(
		prm_index_protein.get_residue_ref_of_index( get_index_of_rep_index( prm_scan_stride.get_index_from_strider(), index_from_rep_index ) ),
		prm_index_protein.get_residue_ref_of_index( get_index_of_rep_index( prm_scan_stride.get_index_to_strider  (), index_to_rep_index   ) ),
		0,
		index_from_rep_index,
		index_to_rep_index
	);
	const auto result               = criteria_result_of( prm_criteria, query_res_pair, index_res_pair );
//	if ( result != quad_criteria_result::PASS && result != quad_criteria_result::QUERY_FAILS_SINGLE_CHECKS && result != quad_criteria_result::INDEX_FAILS_SINGLE_CHECKS ) {
//		quad_criteria_result::PASS;
		cerr << "Quad: [("  << setw(  2 ) << right << prm_query_from_index
		     << ","         << setw(  2 ) << right << prm_index_from_index
		     << ")->("      << setw(  2 ) << right << prm_query_to_index
		     << ","         << setw(  2 ) << right << prm_index_to_index
		     << ")], Rep[(" << setw(  2 ) << right << query_from_rep_index
		     << ","         << setw(  2 ) << right << index_from_rep_index
		     << ")->("      << setw(  2 ) << right << query_to_rep_index
		     << ","         << setw(  2 ) << right << index_to_rep_index
		     << ")]="       << setw( 25 ) << right << result
//		     << ",A/"                              << query_res_pair
//		     << ",B/"                              << index_res_pair
		     << "\n";
//	}

	return result;
}

/// \brief TODOCUMENT
quad_criteria_result check_scan_on_final_alignment::quad_criteria_result_of(const protein       &prm_query_protein,    ///< TODOCUMENT
                                                                            const protein       &prm_index_protein,    ///< TODOCUMENT
                                                                            const quad_criteria &prm_criteria,         ///< TODOCUMENT
                                                                            const index_type    &prm_query_from_index, ///< TODOCUMENT
                                                                            const index_type    &prm_query_to_index,   ///< TODOCUMENT
                                                                            const index_type    &prm_index_from_index, ///< TODOCUMENT
                                                                            const index_type    &prm_index_to_index    ///< TODOCUMENT
                                                                            ) {
	return criteria_result_of(
		prm_criteria,
		make_single_res_pair( prm_query_protein, prm_query_from_index, prm_query_to_index ),
		make_single_res_pair( prm_index_protein, prm_index_from_index, prm_index_to_index )
	);
}

/// \brief TODOCUMENT
quad_and_rep_criteria_result check_scan_on_final_alignment::quad_and_rep_criteria_result_of(const protein       &prm_query_protein,    ///< TODOCUMENT
                                                                                            const protein       &prm_index_protein,    ///< TODOCUMENT
                                                                                            const quad_criteria &prm_criteria,         ///< TODOCUMENT
                                                                                            const scan_stride   &prm_scan_stride,      ///< TODOCUMENT
                                                                                            const index_type    &prm_query_from_index, ///< TODOCUMENT
                                                                                            const index_type    &prm_query_to_index,   ///< TODOCUMENT
                                                                                            const index_type    &prm_index_from_index, ///< TODOCUMENT
                                                                                            const index_type    &prm_index_to_index    ///< TODOCUMENT
                                                                                            ) {

	return {
		rep_quad_criteria_result_of(
			prm_query_protein,
			prm_index_protein,
			prm_criteria,
			prm_scan_stride,
//			from_stride, from_co_stride,
//			entry_index_of_stride_rep() get_query_from_rep_of_index( prm_query_from_index, prm_scan_stride ),
			prm_query_from_index,
			prm_query_to_index,
			prm_index_from_index,
			prm_index_to_index
		),
		quad_criteria_result_of(
			prm_query_protein,
			prm_index_protein,
			prm_criteria,
			prm_query_from_index,
			prm_query_to_index,
			prm_index_from_index,
			prm_index_to_index
		)
	};
}

/// \brief TODOCUMENT
///
/// Things for investigation:
///  * raw/final score from alignment
///  * expected raw/final score
///  * top categorised reasons by raw score
///  * repeat for varying scan_stride (can just be called by nmnf fn)
alignment_scan_comparison check_scan_on_final_alignment::do_check(const alignment     &prm_alignment,  ///< TODOCUMENT
                                                                  const protein       &prm_protein_a,  ///< TODOCUMENT
                                                                  const protein       &prm_protein_b,  ///< TODOCUMENT
                                                                  const quad_criteria &prm_criteria,   ///< TODOCUMENT
                                                                  const scan_stride   &prm_scan_stride ///< TODOCUMENT
                                                                  ) const {
	const auto aln_range = indices( prm_alignment.length() );
	cerr << "SHOULD THE RANGE BE 7.0 RATHER THAN SQRT(40.0)????\n";
	return accumulate(
		cross( aln_range, aln_range ),
		alignment_scan_comparison{},
		[&] (alignment_scan_comparison x, const size_size_tpl &y) {
			const size_t aln_from_ctr   = get<0>( y );
			const size_t aln_to_ctr     = get<1>( y );
			const bool   from_alns_both = has_both_positions_of_index( prm_alignment, aln_from_ctr );
			const bool   to_alns_both   = has_both_positions_of_index( prm_alignment, aln_to_ctr   );

			if ( from_alns_both && to_alns_both ) {

				const auto a_from     = get_a_position_of_index( prm_alignment, aln_from_ctr );
				const auto b_from     = get_b_position_of_index( prm_alignment, aln_from_ctr );
				const auto a_to       = get_a_position_of_index( prm_alignment, aln_to_ctr   );
				const auto b_to       = get_b_position_of_index( prm_alignment, aln_to_ctr   );
				const bool a_included = difference( a_from, a_to ) > NUM_EXCLUDED_ON_SIDES;
				const bool b_included = difference( b_from, b_to ) > NUM_EXCLUDED_ON_SIDES;

				if ( a_included && b_included ) {
					const auto the_distance = distance_between_points(
						view_vector_of_residue_pair(
							prm_protein_a.get_residue_ref_of_index( a_from ),
							prm_protein_a.get_residue_ref_of_index( a_to   )
						),
						view_vector_of_residue_pair(
							prm_protein_b.get_residue_ref_of_index( b_from ),
							prm_protein_b.get_residue_ref_of_index( b_to   )
						)
					);
//					const auto score       = ( the_distance >= 7.0 ) ? 0.0 : ( 1.0 - the_distance / 7.0 );
					const auto score       = ( the_distance >= sqrt( 40.0 ) ) ? 0.0 : ( 1.0 - the_distance / 7.0 );
					if ( score > 0.0 ) {
						const auto scan_result = quad_and_rep_criteria_result_of(
							prm_protein_a,
							prm_protein_b,
							prm_criteria,
							prm_scan_stride,
							numeric_cast<index_type>( a_from ),
							numeric_cast<index_type>( a_to   ),
							numeric_cast<index_type>( b_from ),
							numeric_cast<index_type>( b_to   )
						);
						x += make_pair( scan_result, score );
					}
				}
			}

			return x;
		}
	);
}

/// \brief TODOCUMENT
pair<str_vec, str_vec> check_scan_on_final_alignment::get_rep_name_lists(const protein           &prm_protein,          ///< TODOCUMENT
                                                                         const roled_scan_stride &prm_roled_scan_stride ///< TODOCUMENT
                                                                         ) const {
//	const auto &num_residues        = prm_protein.get_length();
	const auto  rep_list_indices    = get_rep_index_lists( prm_roled_scan_stride, numeric_cast<index_type>( prm_protein.get_length() ) );
	const auto  residue_id_of_index = [&] (const index_type &x) { return get_pdb_residue_id_string( prm_protein.get_residue_ref_of_index( x ) ); };
	return make_pair(
		transform_build<str_vec>( rep_list_indices.first,  residue_id_of_index ),
		transform_build<str_vec>( rep_list_indices.second, residue_id_of_index )
	);
}

/// \brief TODOCUMENT
pair<index_vec, index_vec> cath::scan::detail::get_rep_index_lists(const roled_scan_stride &prm_roled_scan_stride, ///< TODOCUMENT
                                                                   const index_type        &prm_num_residues       ///< TODOCUMENT
                                                                   ) {
	const rep_strider &from_strider  = get_this_from_strider( prm_roled_scan_stride );
	const rep_strider &to_strider    = get_this_to_strider  ( prm_roled_scan_stride );
	const auto         num_from_reps = get_num_reps_of_num_residues( from_strider, prm_num_residues );
	const auto         num_to_reps   = get_num_reps_of_num_residues( to_strider,   prm_num_residues );
	return make_pair(
		transform_build<index_vec>(
			irange<res_rep_index_type>( 0, num_from_reps ),
			[&] (const res_rep_index_type &x) { return get_index_of_rep_index( from_strider, x ); }
		),
		transform_build<index_vec>(
			irange<res_rep_index_type>( 0, num_to_reps ),
			[&] (const res_rep_index_type &x) { return get_index_of_rep_index( from_strider, x ); }
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates
void cath::scan::detail::print_highlight_rep_pymol_commands(ostream                      &prm_os,      ///< TODOCUMENT
                                                            const string                 &prm_title,   ///< TODOCUMENT
                                                            const pair<str_vec, str_vec> &prm_rep_list ///< TODOCUMENT
			                                                ) {
	const str_vec sorted_from_reps = sort_copy( prm_rep_list.first  );
	const str_vec sorted_to_reps   = sort_copy( prm_rep_list.second );

	/// \todo Write a set_intersection_build and a set_symmetric_difference_build and then deploy them here
	str_vec from_and_to_reps;
	set_intersection        ( sorted_from_reps, sorted_to_reps, back_inserter( from_and_to_reps ) );

	str_vec from_xor_to_reps;
	set_symmetric_difference( sorted_from_reps, sorted_to_reps, back_inserter( from_xor_to_reps ) );

	const auto and_colour = string( "black" );
	const auto xor_colour = string( "grey"  );
	const auto and_suffix = string( "_from_and_to_reps" );
	const auto xor_suffix = string( "_from_xor_to_reps" );

	for (const auto &list_colour_and_suffix : { tie( from_and_to_reps, and_colour, and_suffix ),
	                                            tie( from_xor_to_reps, xor_colour, xor_suffix ) } ) {
		const auto &the_list      = get<0>( list_colour_and_suffix );
		const auto &colour        = get<1>( list_colour_and_suffix );
		const auto &suffix        = get<2>( list_colour_and_suffix );
		const auto selection_name = prm_title + suffix;
		const auto join_string    = "/ OR " + prm_title + "///";
		if ( ! the_list.empty() ) {

			prm_os << "select " << selection_name
			       << ", "      << prm_title
			       << "///"     << join( the_list, join_string )
			       << "/\n";
			prm_os << "colour " << colour
			       << ", "      << selection_name
			       << "\n";
			prm_os << "deselect\n";
		}
	}
	prm_os << "show_as lines, (name n,ca,c)\n";
}

/// \brief TODOCUMENT
///
/// \relates
void cath::scan::detail::print_highlight_rep_pymol_commands(ostream                             &prm_os,               ///< TODOCUMENT
                                                            const check_scan_on_final_alignment &prm_csofa,            ///< TODOCUMENT
                                                            const protein                       &prm_protein,          ///< TODOCUMENT
                                                            const roled_scan_stride             &prm_roled_scan_stride ///< TODOCUMENT
                                                            ) {
	print_highlight_rep_pymol_commands(
		prm_os,
		get_domain_or_specified_or_name_from_acq( prm_protein ),
		prm_csofa.get_rep_name_lists(
			prm_protein,
			prm_roled_scan_stride
		)
	);
}

