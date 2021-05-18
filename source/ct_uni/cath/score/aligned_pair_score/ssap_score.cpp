/// \file
/// \brief The ssap_score class definitions

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

#include "ssap_score.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>
#include <boost/serialization/export.hpp>

#include "cath/alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "cath/alignment/gap/alignment_gap.hpp"
#include "cath/alignment/gap/gap_penalty.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/less_than_helper.hpp"
#include "cath/score/length_getter/length_of_longer_getter.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"
#include "cath/ssap/context_res.hpp"

#include <iostream> // ***** TEMPORARY *****

#include "cath/ssap/distance_score_formula.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::cath::score::detail;
using namespace ::std;

using ::boost::lexical_cast;
using ::boost::numeric_cast;
using ::boost::range::join;
using ::boost::tribool;

BOOST_CLASS_EXPORT(ssap_score)

/// \brief TODOCUMENT
///
/// \todo Parameterise on atoms (CA/CB) and whether to use all residues in alignment (common_residue_selection_policy)
score_size_pair ssap_score::calculate_total_score_and_num_quads(const alignment              &prm_alignment,             ///< The pair alignment to be scored
                                                                const protein                &prm_protein_a,             ///< The protein associated with the first  half of the alignment
                                                                const protein                &prm_protein_b,             ///< The protein associated with the first  half of the alignment
                                                                const ssap_score_accuracy    &prm_accuracy,              ///< TODOCUMENT
                                                                const size_t                 &prm_num_excluded_on_sides, ///< TODOCUMENT
                                                                const distance_score_formula &prm_distance_formula       ///< TODOCUMENT
                                                                ) {
	const size_t length       = prm_alignment.length();
	const bool   use_rounding = ( prm_accuracy == ssap_score_accuracy::LOW );

	score_value the_score = 0.0;
	size_t      num_quads = 0;
	for (const size_t &from_ctr : indices( length ) ) {
		if ( has_both_positions_of_index( prm_alignment, from_ctr ) ) {
			const aln_posn_type a_from = get_a_position_of_index( prm_alignment, from_ctr );
			const aln_posn_type b_from = get_b_position_of_index( prm_alignment, from_ctr );

			for (const size_t &to_ctr : indices( length ) ) {
				if ( has_both_positions_of_index( prm_alignment, to_ctr ) ) {
					const aln_posn_type a_to = get_a_position_of_index( prm_alignment, to_ctr );
					const aln_posn_type b_to = get_b_position_of_index( prm_alignment, to_ctr );

					const bool a_included = pair_is_not_excluded( prm_num_excluded_on_sides, a_from, a_to );
					const bool b_included = pair_is_not_excluded( prm_num_excluded_on_sides, b_from, b_to );

					if ( a_included && b_included ) {
						const score_value raw_quad_score = context_res(
							prm_protein_a, prm_protein_b,
							a_from,        b_from,
							a_to,          b_to,
							use_rounding,
							prm_distance_formula
						);
						const score_value quad_score = use_rounding ? floor( raw_quad_score )
						                                            :        raw_quad_score;

//						cerr << a_from << "\t" << a_to << "\t" << b_from << "\t" << b_to << "\t" << quad_score << endl;
//						the_score += quad_score;
						the_score += ( quad_score * residue_querier::RESIDUE_B_VALUE
						                          / residue_querier::RESIDUE_A_VALUE );
						++num_quads;
					}
				}
			}
		}
	}

//	const float_score_type gap_penalty_value = gap_penalty_value_of_alignment( prm_alignment, gap_penalty( 1.0, 0.0 ) );
//	const float_score_float_score_pair open_and_extend_counts = gap_open_and_extend_counts_of_alignment( prm_alignment );
//	cerr << endl;
//	cerr << "open count         : " << open_and_extend_counts.first             << endl;
//	cerr << "extend count       : " << open_and_extend_counts.second            << endl;
//	cerr << "open_gap_penalty   : " << the_gap_penalty.get_open_gap_penalty()   << endl;
//	cerr << "extend_gap_penalty : " << the_gap_penalty.get_extend_gap_penalty() << endl;
//	cerr << "gap_penalty_value  : " << gap_penalty_value                        << endl;
//	cerr << "the_score          : " << the_score                                << endl;
//	cerr << endl;

//	cerr << "the_score is                           " << the_score << endl;
//	cerr << "get_naive_num_gaps( prm_alignment ) is " << get_naive_num_gaps( prm_alignment ) << endl;

	the_score -= numeric_cast<score_value>( get_naive_num_gaps( prm_alignment ) );
//	the_score -= gap_penalty_value_of_alignment( prm_alignment, gap_penalty( 1, 0 ) );
	return make_pair(
		max( min_total_score, the_score),
		num_quads
	);
}

/// \brief TODOCUMENT
score_value ssap_score::log_score_copy(const score_value &prm_score ///< TODOCUMENT
                                       ) {
	const double      final_score_scaling  = 1000.0;
	const score_value k                    = final_score_scaling * residue_querier::RESIDUE_A_VALUE
	                                                             / residue_querier::RESIDUE_B_VALUE;
	const score_value log_k                = log( k );
	return 1.0 + log( prm_score ) / log_k;
}

/// \brief TODOCUMENT
score_value ssap_score::simple_normalise(const score_value &prm_score,           ///< TODOCUMENT
                                         const size_t      &num_aligned_length,  ///< TODOCUMENT
                                         const size_t      &normalisation_length ///< TODOCUMENT
                                         ) {
	return prm_score * numeric_cast<score_value>( num_aligned_length   )
	                 / numeric_cast<score_value>( normalisation_length );
}

/// \brief TODOCUMENT
score_value ssap_score::complex_normalise(const score_value &prm_score,            ///< TODOCUMENT
                                          const size_t      &num_quads,            ///< TODOCUMENT
                                          const size_t      &normalisation_length, ///< TODOCUMENT
                                          const size_t      &num_excluded_on_sides ///< TODOCUMENT
                                          ) {
	const size_t num_comparable = _num_comparable_impl( num_excluded_on_sides, normalisation_length );
	return prm_score * numeric_cast<score_value>( num_quads      )
	                 / numeric_cast<score_value>( num_comparable );
}

/// \brief A standard do_clone method
unique_ptr<aligned_pair_score> ssap_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that says a higher SSAP score's (generally) better
tribool ssap_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the SSAP score of an alignment
score_value ssap_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                     const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                     const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                     ) const {

	const score_size_pair total_score_and_num_quads = calculate_total_score_and_num_quads(
		prm_alignment,
		prm_protein_a,
		prm_protein_b,
		accuracy,
		num_excluded_on_sides,
		distance_formula
	);


	const score_value &total_score = total_score_and_num_quads.first;
	const size_t      &num_quads   = total_score_and_num_quads.second;
	const score_value  mean_score  = total_score / numeric_cast<double>( num_quads );
//	const score_value  mean_score  = total_score * residue_querier::RESIDUE_B_VALUE / ( numeric_cast<double>( num_quads ) * residue_querier::RESIDUE_A_VALUE );

	const bool         pre_log     = has_pre_log ( post_processing );
	const score_value  pre_score   = pre_log ? log_score_copy( mean_score )
	                                         : mean_score;

	const bool         simple_normls        = normalisation_is_simple( post_processing );
	const size_t       num_aligned_length   = num_aligned_length_getter().get_length( prm_alignment, prm_protein_a, prm_protein_b );
	const size_t       normalisation_length = length_getter_ptr->get_length         ( prm_alignment, prm_protein_a, prm_protein_b );
	const score_value  post_score           = simple_normls ? simple_normalise ( pre_score, num_aligned_length, normalisation_length )
	                                                        : complex_normalise( pre_score, num_quads, normalisation_length, num_excluded_on_sides );

	const bool         post_log    = has_post_log( post_processing );
	const score_value  final_score = post_log ? log_score_copy( post_score )
	                                          : post_score;

//	cerr << "total_score is : " << total_score << endl;
//	cerr << "num_quads   is : " << num_quads   << endl;
//	cerr << "mean_score  is : " << mean_score  << endl;
//	cerr << "pre_score   is : " << pre_score   << endl;
//	cerr << "post_score  is : " << post_score  << endl;
//	cerr << "final_score is : " << final_score << endl;
//	cerr << "returned    is : " << 100.0 * final_score << endl;

	return 100.0 * final_score;
}

/// \brief Concrete implementation that describes what this score means
string ssap_score::do_description() const {
	return "The SSAP score"; // TODOCUMENT
//	return "SSAP score over " + length_getter_ptr->description() + the_score_handler.long_suffix_string();
}

/// \brief TODOCUMENT
string ssap_score::do_id_name() const {
	return "ssap";
}

/// \brief TODOCUMENT
str_bool_pair_vec ssap_score::do_short_name_suffixes() const {
	const auto length_getter_suffixes = length_getter_as_short_name_suffixes(
		*length_getter_ptr,
		{ length_getter_category::LONGER }
	);
	const auto coord_handler_suffixes = the_coord_handler.short_name_suffixes();

//	constexpr ssap_score_post_processing ssap_score::default_post_processing;
//	constexpr ssap_score_accuracy        ssap_score::default_accuracy;
//	constexpr size_t                     ssap_score::default_num_excluded_on_sides;
//	constexpr distance_score_formula     ssap_score::default_distance_formula;

	const str_bool_pair_vec other_suffixes = {
		{ lexical_cast<string>( post_processing ),                                  post_processing       != default_post_processing       },
		{ lexical_cast<string>( accuracy        ),                                  accuracy              != default_accuracy              },
		{ "num_excluded_on_sides:" + lexical_cast<string>( num_excluded_on_sides ), num_excluded_on_sides != default_num_excluded_on_sides },
		{ lexical_cast<string>( distance_formula ),                                 distance_formula      != default_distance_formula      },
	};

	return copy_build<str_bool_pair_vec>(
		join(
			join(
				length_getter_suffixes,
				coord_handler_suffixes
			),
			other_suffixes
		)
	);
}

/// \brief Concrete implementation providing long name
string ssap_score::do_long_name() const {
	return "SSAP score over " + length_getter_ptr->long_name() + the_coord_handler.long_suffix_string();
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string ssap_score::do_reference() const {
	return "Taylor W R, Orengo C A (1989) \"Protein Structure Alignment\", Journal of Molecular Biology 208, 1-22 PMID: 2769748";
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> ssap_score::do_build_from_short_name_spec(const string &prm_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                         ) const {
//	cerr << "Should build a ssap_score from string \"" << prm_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool ssap_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                     ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Default ctor for ssap_score (because all three arguments have default values)
ssap_score::ssap_score(const ssap_score_post_processing &prm_post_processing,      ///< (optional) How to post process the basic score once it has been calculated
                       const ssap_score_accuracy        &prm_accuracy,             ///< (optional) Whether to use high accuracy (floating point numbers; no explicit rounding) or low accuracy (ints)
                       const size_t                     &prm_num_excluded_on_sides ///< (optional) The number of residues on each side that are excluded from view calculations
                       ) : post_processing      ( prm_post_processing         ),
                           accuracy             ( prm_accuracy                ),
                           num_excluded_on_sides( prm_num_excluded_on_sides   ) {
}

/// \brief Ctor for ssap_score that allows the caller to specify the protein_only_length_getter
ssap_score::ssap_score(const length_getter              &prm_length_getter,        ///< The method for choosing which of the two proteins should be used in the normalisation
                       const ssap_score_post_processing &prm_post_processing,      ///< (optional) How to post process the basic score once it has been calculated
                       const ssap_score_accuracy        &prm_accuracy,             ///< (optional) Whether to use high accuracy (floating point numbers; no explicit rounding) or low accuracy (ints)
                       const size_t                     &prm_num_excluded_on_sides ///< (optional) The number of residues on each side that are excluded from view calculations
                       ) : length_getter_ptr    ( prm_length_getter.clone() ),
                           post_processing      ( prm_post_processing       ),
                           accuracy             ( prm_accuracy              ),
                           num_excluded_on_sides( prm_num_excluded_on_sides ) {
}

/// \brief Ctor for ssap_score that allows the caller to specify the protein_only_length_getter, common_residue_selection_policy and common_atom_selection_policy
ssap_score::ssap_score(const length_getter                   &prm_length_getter,         ///< The method for choosing which of the two proteins should be used in the normalisation
                       const common_residue_selection_policy &prm_comm_res_seln_pol,     ///< The policy to use for selecting common residues
                       const common_atom_selection_policy    &prm_comm_atom_seln_pol,    ///< The policy to use for selecting common atoms
                       const ssap_score_post_processing      &prm_post_processing,       ///< (optional) How to post process the basic score once it has been calculated
                       const ssap_score_accuracy             &prm_accuracy,              ///< (optional) Whether to use high accuracy (floating point numbers; no explicit rounding) or low accuracy (ints)
                       const size_t                          &prm_num_excluded_on_sides, ///< (optional) The number of residues on each side that are excluded from view calculations
                       const distance_score_formula          &prm_distance_formula       ///< (optional) The formula to use to score the distance between distance within views
                       ) : the_coord_handler     ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol ),
                           length_getter_ptr     ( prm_length_getter.clone() ),
                           post_processing       ( prm_post_processing       ),
                           accuracy              ( prm_accuracy              ),
                           num_excluded_on_sides ( prm_num_excluded_on_sides ),
                           distance_formula      ( prm_distance_formula      ) {
}

/// \brief TODOCUMENT
const length_getter & ssap_score::get_length_getter() const {
	return *length_getter_ptr;
}

/// \brief TODOCUMENT
const ssap_score_post_processing & ssap_score::get_post_processing() const {
	return post_processing;
}

/// \brief TODOCUMENT
const ssap_score_accuracy & ssap_score::get_accuracy() const {
	return accuracy;
}

/// \brief TODOCUMENT
const size_t & ssap_score::get_num_excluded_on_sides() const {
	return num_excluded_on_sides;
}

/// \brief TODOCUMENT
const score_common_coord_handler & ssap_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates ssap_score
bool cath::score::operator<(const ssap_score &prm_ssap_score_a, ///< TODOCUMENT
                            const ssap_score &prm_ssap_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_ssap_score_a, prm_ssap_score_b );
	the_helper.register_comparison_field( &ssap_score::get_length_getter         );
	the_helper.register_comparison_field( &ssap_score::get_post_processing       );
	the_helper.register_comparison_field( &ssap_score::get_accuracy              );
	the_helper.register_comparison_field( &ssap_score::get_num_excluded_on_sides );
	return final_less_than_result( the_helper );
}
