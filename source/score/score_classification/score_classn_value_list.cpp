/// \file
/// \brief The score_classn_value_list class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "score_classn_value_list.h"

#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/transform_build.h"
#include "common/algorithm/sort_uniq_copy.h"
#include "common/boost_addenda/range/adaptor/equal_grouped.h"
#include "common/boost_addenda/range/back.h"
#include "common/boost_addenda/range/front.h"
#include "common/boost_addenda/sorted_insert.h"
#include "common/c++14/cbegin_cend.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "score/score_classification/value_list_scaling.h"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list.h"
#include "score/true_pos_false_neg/true_false_pos_neg.h"
#include "score/true_pos_false_neg/true_false_pos_neg_list.h"

#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::adaptors::reversed;
using boost::adaptors::transformed;
using boost::const_rend;
using boost::filesystem::path;
using boost::irange;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::count_if;
using boost::range::equal;
using boost::range::find_if;
using boost::range::sort;

/// \brief Private function to sort all the values in descending order of goodness
void score_classn_value_list::sort_values() {
	sort( score_classn_values, better_than );
}

/// \brief Ctor for score_classn_value_list
///
/// For now, this is kept private to ensure clients use the factory function make_score_classn_value_list()
score_classn_value_list::score_classn_value_list(const score_classn_value_vec &arg_score_classn_values, ///< TODOCUMENT
                                                 const bool                   &arg_higher_is_better,    ///< TODOCUMENT
                                                 const string                 &arg_name                 ///< TODOCUMENT
                                                 ) : score_classn_values( arg_score_classn_values ),
                                                     better_than        ( arg_higher_is_better    ),
                                                     name               ( arg_name                ) {
	sort_values();
}

/// \brief TODOCUMENT
bool score_classn_value_list::empty() const {
	return score_classn_values.empty();
}

/// \brief TODOCUMENT
size_t score_classn_value_list::size() const {
	return score_classn_values.size();
}

/// \brief TODOCUMENT
const score_classn_value & score_classn_value_list::operator[](const size_t &arg_index ///< TODOCUMENT
                                                               ) const {
	return score_classn_values[ arg_index ];
}

/// \brief TODOCUMENT
score_classn_value_list::const_iterator score_classn_value_list::begin() const {
	return common::cbegin( score_classn_values );
}

/// \brief TODOCUMENT
score_classn_value_list::const_iterator score_classn_value_list::end() const {
	return common::cend( score_classn_values );
}

/// \brief TODOCUMENT
const string score_classn_value_list::get_name() const {
	return name;
}

/// \brief TODOCUMENT
const score_classn_value_better_value & score_classn_value_list::get_better_than() const {
	return better_than;
}

/// \brief TODOCUMENT
void score_classn_value_list::add_score_classn_value(const score_classn_value &arg_score_classn_value ///< TODOCUMENT
                                                     ) {
	sorted_insert( score_classn_values, arg_score_classn_value, better_than );
}

/// \brief Get the best score for a score_classn_value_list
///
/// \relates score_classn_value_list
double cath::score::best_score(const score_classn_value_list &arg_score_classn_value_list ///< The score_classn_value_list to be queried
                               ) {
	if ( arg_score_classn_value_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to retrieve best score from empty score_classn_value_list"));
	}
	return front( arg_score_classn_value_list ).get_score_value();
}

/// \brief Get the worst score for a score_classn_value_list
///
/// \relates score_classn_value_list
double cath::score::worst_score(const score_classn_value_list &arg_score_classn_value_list ///< The score_classn_value_list to be queried
                                ) {
	if ( arg_score_classn_value_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to retrieve worst score from empty score_classn_value_list"));
	}
	return back( arg_score_classn_value_list ).get_score_value();
}

/// \brief Get scaling to normalise the score_classn_value_list
///
/// For now, this just linearly scales to make the worst value 0 and the best value 1.
/// If there is no.
/// Of course, this could be extended to allow other approaches if there's need for it.
///
/// \relates score_classn_value_list
value_list_scaling cath::score::get_scaling(const score_classn_value_list &arg_score_classn_value_list ///< The score_classn_value_list which the scaling should be designed to normalised
                                            ) {
	if ( arg_score_classn_value_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to determine scaling for empty score_classn_value_list"));
	}
	const double the_best_score  = best_score ( arg_score_classn_value_list );
	const double the_worst_score = worst_score( arg_score_classn_value_list );
	const double worst_to_best   = the_best_score - the_worst_score;
	const double multiplier      = ( worst_to_best != 0.0 ) ? ( 1.0 / worst_to_best ) : 1.0;
	const double constant        = - ( multiplier * the_worst_score );

//	value_list_scaling test_scaling( multiplier, constant );
//	cerr << "best is "           << the_best_score;
//	cerr << ", worst is "        << the_worst_score;
//	cerr << ", multiplier is "   << multiplier;
//	cerr << ", constant is "     << constant;
//	cerr << ", scaled_best is "  << scale_value_copy( test_scaling, the_best_score  );
//	cerr << ", scaled_worst is " << scale_value_copy( test_scaling, the_worst_score );
//	cerr << endl;

	return { multiplier, constant };
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
const score_classn_value & cath::score::best_scoring_actual_positive(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                     ) {
	const auto find_itr = find_if( arg_score_classn_value_list, [] (const score_classn_value &x) { return x.get_instance_is_positive(); } );
	if ( find_itr == cath::common::cend( arg_score_classn_value_list ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return *find_itr;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
const score_classn_value & cath::score::best_scoring_actual_negative(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                     ) {
	const auto find_itr = find_if( arg_score_classn_value_list, [] (const score_classn_value &x) { return ! x.get_instance_is_positive(); } );
	if ( find_itr == cath::common::cend( arg_score_classn_value_list ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return *find_itr;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
const score_classn_value & cath::score::worst_scoring_actual_positive(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                      ) {
	const auto find_itr = find_if( arg_score_classn_value_list | reversed, [] (const score_classn_value &x) { return x.get_instance_is_positive(); } );
	if ( find_itr == const_rend( arg_score_classn_value_list ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return *find_itr;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
const score_classn_value & cath::score::worst_scoring_actual_negative(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                      ) {
	const auto find_itr = find_if( arg_score_classn_value_list | reversed, [] (const score_classn_value &x) { return ! x.get_instance_is_positive(); } );
	if ( find_itr == const_rend( arg_score_classn_value_list ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return *find_itr;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
ostream & cath::score::summarise_score_classn_value_list(ostream                       &arg_ostream,                ///< TODOCUMENT
                                                         const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
														 ) {
	arg_ostream << "###############################################" << endl;
	arg_ostream << "# Summary of score_classn_value_list"            << endl;
	arg_ostream << "###############################################" << endl;
	arg_ostream << "# Name                          : " << arg_score_classn_value_list.get_name()              << endl;
	arg_ostream << "# Higher is better              : " << boolalpha << get_higher_is_better( arg_score_classn_value_list ) << endl;
	arg_ostream << "# Size                          : " << arg_score_classn_value_list.size()                  << endl;

	if ( ! arg_score_classn_value_list.empty() ) {
		arg_ostream << "###############################################" << endl;
		arg_ostream << "# Area under ROC curve          : " << area_under_roc_curve( arg_score_classn_value_list ) << endl;
		arg_ostream << "# Best score                    : " << front( arg_score_classn_value_list ).get_score_value() << endl;
		arg_ostream << "# Worst score                   : " << back ( arg_score_classn_value_list ).get_score_value() << endl;
		arg_ostream << "# Best scoring actual positive  : " << best_scoring_actual_positive ( arg_score_classn_value_list ) << endl;
		arg_ostream << "# Best scoring actual negative  : " << best_scoring_actual_negative ( arg_score_classn_value_list ) << endl;
		arg_ostream << "# Worst scoring actual positive : " << worst_scoring_actual_positive( arg_score_classn_value_list ) << endl;
		arg_ostream << "# Worst scoring actual negative : " << worst_scoring_actual_negative( arg_score_classn_value_list ) << endl;
	}
	arg_ostream << "###############################################" << endl;
	return arg_ostream;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
str_set cath::score::get_sorted_instance_labels(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                ) {
	const auto instance_labels = get_instance_labels( arg_score_classn_value_list );
	return {
		cath::common::cbegin( instance_labels ),
		cath::common::cend  ( instance_labels )
	};
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
str_vec cath::score::get_instance_labels(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                         ) {
	return copy_build<str_vec>(
		arg_score_classn_value_list
			| transformed( [] (const score_classn_value &x) { return x.get_instance_label(); } )
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
str_size_pair_vec cath::score::get_sorted_instance_labels_and_indices(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                      ) {
	return sort_copy(
		transform_build<str_size_pair_vec>(
			irange( 0_z, arg_score_classn_value_list.size() ),
			[&] (const size_t &x) { return make_pair( arg_score_classn_value_list[ x ].get_instance_label(), x ); }
		),
		[] (const str_size_pair &x, const str_size_pair &y) { return x.first < y.first; }
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
score_classn_value_vec cath::score::get_score_classn_values_of_instance_labels(const score_classn_value_list &arg_score_classn_value_list, ///< TODOCUMENT
                                                                               const str_vec                 &arg_instance_labels          ///< TODOCUMENT,
                                                                               ) {
	const auto instance_labels_and_indices = get_sorted_instance_labels_and_indices( arg_score_classn_value_list );
	return transform_build<score_classn_value_vec>(
		arg_instance_labels,
		[&] (const string &x) {
			const auto find_itr = find_if(
				instance_labels_and_indices,
				[&] (const str_size_pair &y) { return x == y.first; }
			);
			if ( find_itr == cath::common::cend( instance_labels_and_indices ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find instance label in score_classn_value_list"));
			}
			return arg_score_classn_value_list[ find_itr->second ];
		}
	);

}

/// \brief TODOCUMENT
///
/// This isn't necessarily a massively efficient approach and it could probably be substantially improved
/// if profiling flags it as a bottleneck
///
/// \relates score_classn_value_list
bool cath::score::instance_labels_match(const score_classn_value_list &arg_score_classn_value_list_a, ///< TODOCUMENT
                                        const score_classn_value_list &arg_score_classn_value_list_b  ///< TODOCUMENT
                                        ) {
	const auto instance_labels_a = get_sorted_instance_labels( arg_score_classn_value_list_a );
	const auto instance_labels_b = get_sorted_instance_labels( arg_score_classn_value_list_b );

//	cerr << "names_a has " << instance_labels_a.size() << " entries" << endl;
//	cerr << "names_b has " << instance_labels_b.size() << " entries" << endl;

	return ( instance_labels_a.size() == arg_score_classn_value_list_a.size() )
	       &&
	       ( instance_labels_b.size() == arg_score_classn_value_list_b.size() )
	       &&
		   equal( instance_labels_a, instance_labels_b );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
const bool & cath::score::get_higher_is_better(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                               ) {
	return arg_score_classn_value_list.get_better_than().get_higher_is_better();
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
double cath::score::worst_possible_score(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                         ) {
	return cath::score::get_worst_possible_score( arg_score_classn_value_list.get_better_than() );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
score_classn_value_list cath::score::make_score_classn_value_list(const score_classn_value_vec &arg_score_classn_values, ///< TODOCUMENT
                                                                  const bool                   &arg_higher_is_better,    ///< TODOCUMENT
                                                                  const string                 &arg_name                 ///< TODOCUMENT
                                                                  ) {
	return score_classn_value_list( arg_score_classn_values, arg_higher_is_better, arg_name );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
doub_doub_pair_vec cath::score::correlated_data(const score_classn_value_list &arg_score_classn_values_a, ///< TODOCUMENT
                                                const score_classn_value_list &arg_score_classn_values_b  ///< TODOCUMENT
                                                ) {
	if ( ! instance_labels_match( arg_score_classn_values_a, arg_score_classn_values_b ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get correlated data for score_classn_value_lists with mismatching instance labels"));
	}
	const auto index_b_of_label = transform_build<str_size_map>(
		irange( 0_z, arg_score_classn_values_b.size() ),
		[&] (const size_t &x) {
			return make_pair( arg_score_classn_values_b[ x ].get_instance_label(), x );
		}
	);

	return transform_build<doub_doub_pair_vec>(
		irange( 0_z, arg_score_classn_values_a.size() ),
		[&] (const size_t &x) {
			const auto &score_classn_value_a = arg_score_classn_values_a[ x ];
			const auto &instance_label       = score_classn_value_a.get_instance_label();
			const auto &score_classn_value_b = arg_score_classn_values_b[ index_b_of_label.at( instance_label ) ];
			return make_pair(
				score_classn_value_a.get_score_value(),
				score_classn_value_b.get_score_value()
			);
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
score_classn_value_list cath::score::read_svmlight_predictions_file(const path   &arg_path, ///< TODOCUMENT
                                                                    const string &arg_name  ///< TODOCUMENT
                                                                    ) {
	constexpr bool higher_is_better = true;
	return read_score_classn_value_list(
		arg_path,
		higher_is_better,
		arg_name,
		[] (const str_vec &x) { return ( stod( x.at( 2 ) ) >= 0.0 ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
score_classn_value_list_vec cath::score::read_svmlight_predictions_files(const vector<pair<path, string>> &arg_paths_and_names ///< TODOCUMENT
                                                                         ) {
	return transform_build<score_classn_value_list_vec>(
		arg_paths_and_names,
		[] (const pair<path, string> &x) {
			return read_svmlight_predictions_file( x.first, x.second );
		}
	);
}

/// \brief TODOCUMENT
///
/// \todo Create an adaptor called something like equal_grouped that returns a range whose
///       iterators dereference to sub_ranges for groups of equivalent entries according to some
///       less_than predicate. Then, simplify this code accordingly.
///
/// \relates score_classn_value_list
named_true_false_pos_neg_list cath::score::make_named_true_false_pos_neg_list(const score_classn_value_list &arg_score_classn_value_list ///< TODOCUMENT
                                                                              ) {
	// Grab functors for determining:
	//  * whether a score_classn_value represents a positive instance
	//  * whether one score_classn_value is better_than another
	const auto reps_positive_instance = [] (const score_classn_value &x) { return x.get_instance_is_positive(); };
	const auto better_than            = arg_score_classn_value_list.get_better_than();

	// Count the number of positive examples and subtract it from the total to get the total number of negatives
	const size_t total_num           = arg_score_classn_value_list.size();
	const size_t total_num_positives = numeric_cast<size_t>( count_if( arg_score_classn_value_list, reps_positive_instance ) );
	const size_t total_num_negatives = total_num - total_num_positives;

	// Initialise running_tfpn_counts to have no positives and the correct number of true/false negatives
	auto running_tfpn_counts = true_false_pos_neg{ 0, total_num_negatives, 0, total_num_positives };
	true_false_pos_neg_vec tfpns = { running_tfpn_counts };

	// Loop over the groups of equivalent score_classn_values
	// (the unequal predicate only need test for better_than because the values are already sorted by better_than)
	for (const auto &equivalent_score_classn_values : arg_score_classn_value_list | equal_grouped( better_than ) ) {

		// Within the group of equivalent score_classn_values, count the number of entries, positives and negatives
		const size_t num_in_group  = equivalent_score_classn_values.size();
		const size_t num_positives = numeric_cast<size_t>( count_if( equivalent_score_classn_values, reps_positive_instance ) );
		const size_t num_negatives = num_in_group - num_positives;

		// Update the true_false_pos_neg counts with the number of newly predicted positives and negatives
		// and then add that new group of counts to the back of tfpns
		update_with_predicted_positives( running_tfpn_counts, num_positives, num_negatives);
		tfpns.push_back( running_tfpn_counts );
	}

	// Return a true_false_pos_neg_list constructed from the true_false_pos_neg_vec
	return { tfpns, arg_score_classn_value_list.get_name() };
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
double cath::score::area_under_curve(const score_classn_value_list &arg_tfpns,         ///< TODOCUMENT
                                     const classn_stat             &arg_classn_stat_x, ///< TODOCUMENT
                                     const classn_stat             &arg_classn_stat_y  ///< TODOCUMENT
                                     ) {
	return area_under_curve(
		make_named_true_false_pos_neg_list( arg_tfpns ),
		arg_classn_stat_x,
		arg_classn_stat_y
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_list
double cath::score::area_under_roc_curve(const score_classn_value_list &arg_tfpns ///< TODOCUMENT
                                         ) {
	return area_under_roc_curve( make_named_true_false_pos_neg_list( arg_tfpns ) );
}


//const double & get_score_value() const;
//const bool & get_instance_is_positive() const;
