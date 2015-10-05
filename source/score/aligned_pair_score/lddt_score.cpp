/// \file
/// \brief The lddt_score class definitions

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

#include "lddt_score.h"

#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/range/numeric.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"
#include "common/algorithm/copy_build.h"
#include "common/clone/make_uptr_clone.h"
#include "common/difference.h"
#include "common/less_than_helper.h"
#include "common/size_t_literal.h"
#include "exception/not_implemented_exception.h"
#include "exception/out_of_range_exception.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/coord_list.h"

#include <iostream> // ***** TEMPORARY *****
#include <numeric>

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace std;

using boost::accumulate;
using boost::algorithm::none_of;
using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;
using boost::range::join;

BOOST_CLASS_EXPORT(lddt_score)

constexpr double lddt_score::THRESHOLD_FOR_ALL;

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> lddt_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool lddt_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value lddt_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
                                     const protein   &arg_protein_a, ///< The protein associated with the first  half of the alignment
                                     const protein   &arg_protein_b  ///< The protein associated with the second half of the alignment
                                     ) const {
	// As described in the paper, use a R_0 value of 15.0 angstroms
	const score_value R_0 = 15.0;

	// A threshold value that will ensure everything gets considered
	

	// Extract the common coordinates to be chosen
	const pair<coord_list_vec, coord_list_vec> common_coords_by_residue = the_coord_handler.get_common_coords_by_residue(
		arg_alignment,
		arg_protein_a,
		arg_protein_b
	);

	// Check that there are some coords
	const size_t num_common_residues = common_coords_by_residue.first.size();
	if ( num_common_residues == 0 ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"Cannot calculate "
			+ do_id_name()
			+ " for alignment from which no common residues are extracted"
		));
	}
	// If none of the residue common coord entries has size() > 0, then throw an exception
	if ( none_of( common_coords_by_residue.first, [] (const coord_list &x) { return x.size() > 0; } ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"Cannot calculate " + do_id_name() + " for alignment from which no common coords are extracted"
		));
	}

	// Construct a map from threshold values to the counts of residues pairs that satisfy each
	using doub_size_map = std::map< double, size_t >;
	using doub_size_map_val = doub_size_map::value_type;
	doub_size_map count_within_threshold = { { THRESHOLD_FOR_ALL, 0 } };
	for (const score_value &threshold_value : threshold_values) {
		count_within_threshold[ threshold_value ] = 0;
	}

	// Loop over the pairs of residues
	for (size_t res_ctr_1 = 0; res_ctr_1 < num_common_residues; ++res_ctr_1) {
		for (size_t res_ctr_2 = res_ctr_1 + 1; res_ctr_2 < num_common_residues; ++res_ctr_2) {

			const size_t &num_atoms_1 = common_coords_by_residue.first[res_ctr_1].size();
			const size_t &num_atoms_2 = common_coords_by_residue.first[res_ctr_2].size();

			for (size_t atom_ctr_1 = 0; atom_ctr_1 < num_atoms_1; ++atom_ctr_1) {
				for (size_t atom_ctr_2 = 0; atom_ctr_2 < num_atoms_2; ++atom_ctr_2) {
					const double distance_a  = distance_between_points(
						common_coords_by_residue.first[  res_ctr_1 ][ atom_ctr_1 ],
						common_coords_by_residue.first[  res_ctr_2 ][ atom_ctr_2 ]
					);
					const double distance_b  = distance_between_points(
						common_coords_by_residue.second[  res_ctr_1 ][ atom_ctr_1 ],
						common_coords_by_residue.second[  res_ctr_2 ][ atom_ctr_2 ]
					);
					const size_t num_suitable_pairs =   ( distance_a < R_0 ? 1_z : 0_z )
					                                  + ( distance_b < R_0 ? 1_z : 0_z );

					const double distance_difference = difference( distance_a, distance_b );

					for (doub_size_map_val &threshold_count_pair : count_within_threshold) {
						if (distance_difference < threshold_count_pair.first) {
							threshold_count_pair.second += num_suitable_pairs;
						}
					}
				}
			}
		}
	}

	const size_t all_count = count_within_threshold.at( THRESHOLD_FOR_ALL );
	vector<score_value> fractions;
	fractions.reserve( count_within_threshold.size() - 1 );
	for (const doub_size_map_val &threshold_count : count_within_threshold) {
		const double &threshold = threshold_count.first;
		const size_t &count     = threshold_count.second;
		if (threshold != THRESHOLD_FOR_ALL) {
			fractions.push_back(
				numeric_cast<score_value>( count ) / numeric_cast<score_value>( all_count  )
			);
		}
	}
	const score_value fraction_total = accumulate( fractions, 0.0 );
	const score_value fraction_avg = fraction_total / numeric_cast<score_value>( fractions.size() );
	return fraction_avg;
}

/// \brief Concrete implementation that describes what this score means
string lddt_score::do_description() const {
	return "The root of the mean squared deviation between equivalent atoms"
	       + the_coord_handler.description_brackets_string()
	       + ". This is measured in angstroms.";
}

/// \brief TODOCUMENT
string lddt_score::do_id_name() const {
	return "lDDT";
}

/// \brief TODOCUMENT
str_bool_pair_vec lddt_score::do_short_name_suffixes() const {
	const auto threshold_suffixes = { make_pair(
		string( "threshold[" ) + get_thresholds_summary_string( *this ) + "]",
		true
	) };
	return copy_build<str_bool_pair_vec>( join (
		threshold_suffixes,
		the_coord_handler.short_name_suffixes()
	) );
}

/// \brief Concrete implementation providing long name
string lddt_score::do_long_name() const {
	return " Local Distance Difference Test"
	       + the_coord_handler.long_suffix_string()
	       + ". This is measured in angstroms.";
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string lddt_score::do_reference() const {
	return { "Mariani V, Biasini M, Barbato A, Schwede T. "
	         "lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. "
	         "Bioinformatics. 2013 Nov 1;29(21):2722-8. doi: 10.1093/bioinformatics/btt473. "
	         "Epub 2013 Aug 27. PubMed PMID: 23986568; PubMed Central PMCID: PMC3799472." };
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> lddt_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                         ) const {
//	cerr << "Should build a lddt_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool lddt_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                     ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Initialise the actual numeric thresholds over which the results should be averages
///        based on the specified policy
void lddt_score::init_distance_threshold_values(const lddt_distance_threshold &arg_distance_threshold ///< TODOCUMENT
                                                ) {
	switch (arg_distance_threshold) {
		case ( lddt_distance_threshold::DEFAULT_AVG ) : {
			threshold_values = { 0.5, 1.0, 2.0, 4.0 };
			return;
		}
		case ( lddt_distance_threshold::HALF_A      ) : {
			threshold_values = { 0.5 };
			return;
		}
		case ( lddt_distance_threshold::ONE_A       ) : {
			threshold_values = { 1.0 };
			return;
		}
		case ( lddt_distance_threshold::TWO_A       ) : {
			threshold_values = { 2.0 };
			return;
		}
		case ( lddt_distance_threshold::FOUR_A      ) : {
			threshold_values = { 4.0 };
			return;
		}
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("Failed to handle value of lddt_distance_threshold enum"));
//	return;
}

/// \brief Default ctor for lddt_score that uses the defaults for score_common_coord_handler (selecting CA atoms for all aligned residues)
lddt_score::lddt_score() {
	init_distance_threshold_values( lddt_distance_threshold::DEFAULT_AVG );
}

/// \brief Ctor for lddt_score that allows the caller to specify the distance threshold
///        but just use sensible defaults for the common_residue_selection_policy and common_atom_selection_policy
lddt_score::lddt_score(const lddt_distance_threshold &arg_distance_threshold ///< The standard lDDT distance threshold that should be used
                       ) {
	init_distance_threshold_values( arg_distance_threshold );
}

/// \brief Ctor for lddt_score that allows the caller to specify the common_residue_selection_policy and common_atom_selection_policy
lddt_score::lddt_score(const lddt_distance_threshold         &arg_distance_threshold, ///< The standard lDDT distance threshold that should be used
                       const common_residue_selection_policy &arg_comm_res_seln_pol,  ///< The policy to use for selecting common residues
                       const common_atom_selection_policy    &arg_comm_atom_seln_pol  ///< The policy to use for selecting common atoms
                       ) : the_coord_handler( arg_comm_res_seln_pol, arg_comm_atom_seln_pol ) {
	init_distance_threshold_values( arg_distance_threshold );
}

/// \brief TODOCUMENT
const doub_vec & lddt_score::get_threshold_values() const {
	return threshold_values;
}

/// \brief TODOCUMENT
const score_common_coord_handler & lddt_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates lddt_score
bool cath::score::operator<(const lddt_score &arg_lddt_score_a, ///< TODOCUMENT
                            const lddt_score &arg_lddt_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_lddt_score_a, arg_lddt_score_b );
	the_helper.register_comparison_field( &lddt_score::get_threshold_values           );
	the_helper.register_comparison_field( &lddt_score::get_score_common_coord_handler );
	return final_less_than_result( the_helper );
}

/// \brief Get a string that summarises the lDDT threshold policy being used by the specified lddt_score
///
/// \relates lddt_score
string cath::score::get_thresholds_summary_string(const lddt_score &arg_lddt_score ///< TODOCUMENT
                                                  ) {
	const doub_vec the_thresholds = arg_lddt_score.get_threshold_values();
	if ( the_thresholds.size() == 1 ) {
		return lexical_cast<string>( the_thresholds.front() );
	}
	const doub_vec avg_vector = { 0.5, 1.0, 2.0, 4.0 };
	if ( the_thresholds == avg_vector ) {
		return "STD_MEAN";
	}

	BOOST_THROW_EXCEPTION(not_implemented_exception("lddt_score cannot handle an unrecognised thresholds list"));
	return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}

/// \brief Return all possible lddt_distance_threshold values
///
/// This can be used in the code that aligned_pair_score_list factory code for generating
/// all possible lDDT scores
vector<lddt_distance_threshold> cath::score::get_all_lddt_distance_thresholds() {
	return { lddt_distance_threshold::DEFAULT_AVG,
	         lddt_distance_threshold::HALF_A,
	         lddt_distance_threshold::ONE_A,
	         lddt_distance_threshold::TWO_A,
	         lddt_distance_threshold::FOUR_A };
}
