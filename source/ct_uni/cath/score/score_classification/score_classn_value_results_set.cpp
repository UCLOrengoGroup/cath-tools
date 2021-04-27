/// \file
/// \brief The score_classn_value_results_set class definitions

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

#include "score_classn_value_results_set.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm_ext/is_sorted.hpp>
#include <boost/range/irange.hpp>

#include <spdlog/spdlog.h>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/random_split.hpp"
#include "cath/common/algorithm/set_difference_build.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "cath/common/boost_addenda/range/adaptor/limited.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/tribool/tribool.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"
#include "cath/score/score_classification/value_list_scaling.hpp"
#include "cath/score/true_pos_false_neg/classn_rate_stat.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_pair_series_list.hpp"
#include "cath/score/true_pos_false_neg/named_true_false_pos_neg_list.hpp"
#include "cath/score/true_pos_false_neg/named_true_false_pos_neg_list_list.hpp"
#include "cath/score/true_pos_false_neg/true_false_pos_neg.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::cath::score::detail;
using namespace ::std;

using ::boost::algorithm::join;
//using ::boost::algorithm::join; // ***** TEMPORARY *****
using ::boost::filesystem::path;
using ::boost::irange;
using ::boost::lexical_cast;
using ::boost::range::adjacent_find;
using ::boost::range::equal;
using ::boost::range::find_if;
using ::boost::range::includes;
using ::boost::range::is_sorted;
using ::boost::range::lower_bound;
using ::boost::range::sort;

/// \brief Determine whether the score_classn_value_lists is correctly sorted+uniqued on the names
bool score_classn_value_results_set::is_sorted_uniqued() const {
	// Try to find a pair with decreasing or equal names
	const auto decreasing_or_equal_pair_itr = adjacent_find(
		score_classn_value_lists,
		[] (const score_classn_value_list &x, const score_classn_value_list &y) { return ! score_classn_value_list_name_less()( y, x ); }
	);
	// Return whether none was found
	return ( decreasing_or_equal_pair_itr == common::cend( score_classn_value_lists ) );
}

/// \brief Check that score_classn_value_lists is correctly sorted+uniqued on the names and throw if not
void score_classn_value_results_set::check_is_sorted_uniqued() const {
	if ( ! is_sorted_uniqued() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("score_classn_value_results_set has violated its invariant of keeping score_classn_value_lists sorted on name"));
	}
}

/// \brief Private member function to sort the score_classn_value_lists on the names
void score_classn_value_results_set::sort_score_classn_value_lists() {
	sort( score_classn_value_lists, score_classn_value_list_name_less() );
}

/// \brief TODOCUMENT
score_classn_value_list & score_classn_value_results_set::get_score_classn_value_list_of_name(const string &prm_name ///< TODOCUMENT
                                                                                              ) {
	return get_score_classn_value_list_of_name_impl( *this, prm_name );
}

/// \brief TODOCUMENT
const score_classn_value_list & score_classn_value_results_set::get_score_classn_value_list_of_name(const string &prm_name ///< TODOCUMENT
                                                                                                    ) const {
	return get_score_classn_value_list_of_name_impl( *this, prm_name );
}

/// \brief TODOCUMENT
bool score_classn_value_results_set::empty() const {
	return score_classn_value_lists.empty();
}

/// \brief TODOCUMENT
size_t score_classn_value_results_set::size() const {
	return score_classn_value_lists.size();
}

/// \brief TODOCUMENT
const score_classn_value_list & score_classn_value_results_set::operator[](const size_t &prm_index ///< TODOCUMENT
                                                                           ) const {
	return score_classn_value_lists[ prm_index ];
}

/// \brief TODOCUMENT
score_classn_value_results_set::const_iterator score_classn_value_results_set::begin() const {
	return cath::common::cbegin( score_classn_value_lists );
}

/// \brief TODOCUMENT
score_classn_value_results_set::const_iterator score_classn_value_results_set::end() const {
	return cath::common::cend( score_classn_value_lists );
}

/// \brief TODOCUMENT
void score_classn_value_results_set::add_score_classn_value_list(const score_classn_value_list &prm_score_classn_value_list ///< TODOCUMENT
                                                                 ) {


	// If there are already entries, confirm that the new entry's list of instances matches the
	// list of instances in one of the existing entries
	if ( ! score_classn_value_lists.empty() ) {
		if ( ! instance_labels_match( prm_score_classn_value_list, score_classn_value_lists.front() ) ) {
			if ( score_classn_value_lists.front().empty() || prm_score_classn_value_list.empty() ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add score_classn_value_list with instance labels that don't match those already present in this score_classn_value_results_set (note: one has no labels at all)"));
			}
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Cannot add score_classn_value_list with instance labels that don't match those already present in this score_classn_value_results_set (example label present: \""
				+ front( score_classn_value_lists.front() ).get_instance_label()
				+ "\"; example label from score_classn_value_list to add: \""
				+ front( prm_score_classn_value_list ).get_instance_label()
				+ "\")"
			));
		}
	}

	// Find the correct place to insert the new list
	const auto insert_location_itr = lower_bound(
		score_classn_value_lists,
		prm_score_classn_value_list,
		score_classn_value_list_name_less()
	);

	// If the list wouldn't be inserted at the end, check that the existing entry in the location
	// doesn't have the same name and throw if it does
	if ( insert_location_itr != common::cend( score_classn_value_lists ) ) {
		if ( insert_location_itr->get_name() == prm_score_classn_value_list.get_name() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
		}
	}

	// Insert the new entry at the previously identified location
	score_classn_value_lists.insert( insert_location_itr, prm_score_classn_value_list );
}

/// \brief TODOCUMENT
///
/// \todo Could this be generalised (perhaps with a non-member function
///       to make it convenient for the particulars of aligned_pair_score_value_list) ?
void score_classn_value_results_set::add_aligned_pair_score_value_list(const aligned_pair_score_value_list &prm_aligned_pair_score_value_list, ///< TODOCUMENT
                                                                       const bool                          &prm_condition_is_positive,         ///< TODOCUMENT
                                                                       const string                        &prm_instance_label                 ///< TODOCUMENT
                                                                       ) {
	const size_t num_score_values = prm_aligned_pair_score_value_list.size();
	if ( empty() ) {
		for (const size_t &score_value_ctr : indices( num_score_values ) ) {
			const auto   &score            = prm_aligned_pair_score_value_list.get_aligned_pair_score_of_index( score_value_ctr );
			const string &name             = score.human_friendly_short_name();
			const bool   &higher_is_better = is_true( score.higher_is_better() );
			score_classn_value_lists.push_back( make_score_classn_value_list( {}, higher_is_better, name ) );
		}
		sort_score_classn_value_lists();
	}
	else {
		const auto    new_names      = get_all_names( prm_aligned_pair_score_value_list );
		const str_set new_names_set( common::cbegin( new_names ), common::cend( new_names ) );
		const size_t  num_series     = size();
		const auto    existing_names = transform_build<str_vec>(
			indices( num_series ),
			[&] (const size_t &x) { return get_name_of_index( *this, x ); }
		);

		if ( ! equal( existing_names, new_names_set ) ) {
//			cerr << "existing_names : " << join( existing_names, ", " ) << endl;
//			cerr << "new_names_set  : " << join( new_names_set,  ", " ) << endl;
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add aligned_pair_score_value_list because the series names don't match the existing ones"));
		}
	}

	for (const size_t &score_value_ctr : indices( num_score_values ) ) {
		const auto   &value            = prm_aligned_pair_score_value_list.get_value_of_index             ( score_value_ctr );
		const auto   &score            = prm_aligned_pair_score_value_list.get_aligned_pair_score_of_index( score_value_ctr );
		const string &name             = score.human_friendly_short_name();
		const bool   &higher_is_better = is_true( score.higher_is_better() );

		score_classn_value_list &the_list = get_score_classn_value_list_of_name( name );
		if ( get_higher_is_better( the_list ) != higher_is_better ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception( "Cannot add aligned_pair_score_value_list to score_classn_value_results_set because they conflict regarding the higher_is_better value for score " + name ));
		}

		the_list.add_score_classn_value(
			score_classn_value(
				value,
				prm_condition_is_positive,
				prm_instance_label
			)
		);
	}
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
score_classn_value_results_set cath::score::make_score_classn_value_results_set(const score_classn_value_list_vec &prm_score_classn_value_lists ///< TODOCUMENT
                                                                                ) {
	score_classn_value_results_set new_set;
	for (const auto &scv_list : prm_score_classn_value_lists) {
		add_score_classn_value_list_and_add_missing( new_set, scv_list, worst_possible_score( scv_list ) );
	}
	return new_set;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
score_classn_value_list_vec cath::score::make_score_classn_value_list_vec(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                                                                          ) {
	return copy_build<score_classn_value_list_vec>( prm_results_set );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
value_list_scaling_vec cath::score::get_value_list_scalings(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                                                            ) {
	return transform_build<value_list_scaling_vec>(
		prm_results_set,
		[] (const score_classn_value_list &x) {
			return get_scaling( x );
		}
	);
}

/// \brief TODOCUMENT
doub_doub_pair_vec cath::score::get_correlated_data(const score_classn_value_results_set &prm_results,        ///< TODOCUMENT
                                                    const string                         &prm_name_of_first, ///< TODOCUMENT
                                                    const string                         &prm_name_of_second ///< TODOCUMENT
                                                    ) {
	const auto &first_data   = prm_results.get_score_classn_value_list_of_name( prm_name_of_first  );
	const auto &second_data  = prm_results.get_score_classn_value_list_of_name( prm_name_of_second );
	return correlated_data( first_data, second_data );
}

/// \brief TODOCUMENT
void cath::score::detail::write_to_svm_light_data_files_impl(const score_classn_value_vec_vec &prm_results,             ///< TODOCUMENT
                                                             const value_list_scaling_vec     &prm_value_list_scalings, ///< TODOCUMENT
                                                             const path                       &prm_output_file,         ///< TODOCUMENT
                                                             const size_vec                   &prm_indices              ///< TODOCUMENT
                                                             ) {
	// Sanity check the inputs
	if ( ! is_sorted( prm_indices ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write SVM data for indices that aren't sorted"));
	}
	if ( contains_adjacent_match( prm_results, [] (const score_classn_value_vec &x, const score_classn_value_vec &y) { return x.size() != y.size(); } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Not all results to be written to an SVM data file are of equal size"));
	}
	if ( contains_if( prm_results, [&] (const score_classn_value_vec &x) { return x.size() <= prm_indices.back(); } ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The results to be written to an SVM data file are not big enough for the specified indices"));
	}
	if ( prm_value_list_scalings.size() != prm_results.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The number of scalings doesn't match the number of score_classn_value_vecs when attempting to write SVM data "));
	}
//	if ( ! results.empty() ) {
//		const auto first_results_size = results.front().size();
//		if ( ! all_of( results, [] (const score_classn_value_vec &x) { return x.size() == first_results_size; } ) ) {
//			BOOST_THROW_EXCEPTION(invalid_argument_exception("Not all results to be written to an SVM data file are of equal size"));
//		}
//	}
//	if ( ! prm_indices.empty() && ! results.empty() ) {
//		if ( prm_indices.back() >= prm_results.front().size() ) {
//			BOOST_THROW_EXCEPTION(invalid_argument_exception("The results to be written to an SVM data file are not big enough for the specified indices"));
//		}
//	}

	// Open an ostream for the file
	ofstream out_stream;
	open_ofstream( out_stream, prm_output_file );

	// Loop over the instances of the requested indices
	for (const size_t &index : prm_indices) {

		// Check for any mismatching entries wrt their instance label/is_positive value
		const bool mismatching_instances = contains_adjacent_match(
			prm_results,
			[&] (const score_classn_value_vec &x, const score_classn_value_vec &y) {
				return (    ( x[ index ].get_instance_is_positive() != y[ index ].get_instance_is_positive() )
						 || ( x[ index ].get_instance_label()       != y[ index ].get_instance_label()       ) );
			}
		);
		if ( mismatching_instances ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("When attempting to write SVM data, detected mismatching entries wrt their instance label/is_positive value"));
		}

		// If there are prm_results, then output a line of data
		if ( ! prm_results.empty() ) {
			const auto &first_value = prm_results.front()[ index ];

			// First, output '+1 ' if this instance is positive and '-1 ' otherwise
			out_stream << ( first_value.get_instance_is_positive() ? string( "+1 " ) : string( "-1 " ) );

			// Next output the data in the format : '1:value_1 2:value_2 [...] n:value_n'
			out_stream << join(
				transform_build<str_vec>(
					indices( prm_results.size() ),
					[&] (const size_t &x) {
						const auto &scaling      = prm_value_list_scalings[ x ];
						const auto  scaled_score = scale_value_copy( scaling, prm_results[ x ][ index ].get_score_value() );
						return lexical_cast<string>( x + 1 ) + ":" + lexical_cast<string>( scaled_score );
					}
				),
				" "
			);

			// Finally, append a comment containing the label of the instance
			out_stream << " # " << first_value.get_instance_label() << "\n";
		}
	}

	// Flush and close the output stream
	out_stream << flush;
	out_stream.close();
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
void cath::score::write_to_svm_light_data_files(const score_classn_value_results_set &prm_results_set,      ///< TODOCUMENT
                                                const path                           &prm_output_file_stem, ///< TODOCUMENT
                                                const size_t                         &prm_num_repeats,      ///< TODOCUMENT
                                                mt19937                              &prm_rng,              ///< TODOCUMENT
                                                const double                         &prm_fraction_train    ///< TODOCUMENT
                                                ) {
	const auto num_instances   = get_num_instances      ( prm_results_set );
	const auto names           = get_names              ( prm_results_set );
	const auto instance_labels = get_instance_labels    ( prm_results_set );
	const auto scalings        = get_value_list_scalings( prm_results_set );

	// Get a data structure in which all lists are sorted in the same way
	const auto results = transform_build<score_classn_value_vec_vec>(
		names,
		[&] (const string &x) {
			const auto &results_list = prm_results_set.get_score_classn_value_list_of_name( x );
			return get_score_classn_values_of_instance_labels( results_list, instance_labels );
		}
	);

	for (const size_t &repeat_ctr : irange( 1_z, prm_num_repeats + 1 ) ) {
		const auto repeat_ctr_str = lexical_cast<string>( repeat_ctr );

		const path train_file  = replace_extension_copy( prm_output_file_stem, "." + repeat_ctr_str + ".train" );
		const path test_file   = replace_extension_copy( prm_output_file_stem, "." + repeat_ctr_str + ".test"  );

		const size_vec_size_vec_pair split_indices = random_split( prm_rng, num_instances, prm_fraction_train );

		write_to_svm_light_data_files_impl( results, scalings, train_file, split_indices.first  );
		write_to_svm_light_data_files_impl( results, scalings, test_file,  split_indices.second );
	}
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
void cath::score::add_score_classn_value_list_and_add_missing(score_classn_value_results_set &prm_results_set,             ///< TODOCUMENT
                                                              score_classn_value_list         prm_score_classn_value_list, ///< TODOCUMENT
                                                              const double                   &prm_score_for_missing,       ///< TODOCUMENT
                                                              const bool                     &prm_warn_on_missing          ///< TODOCUMENT
                                                              ) {
	// If the score_classn_value_results_set is empty, just do a normal add and then return
	if ( prm_results_set.empty() ) {
		prm_results_set.add_score_classn_value_list( prm_score_classn_value_list );
		return;
	}

	// Grab sorted lists of the instance labels in the results_set and value_list
	const auto results_set_instance_labels = get_sorted_instance_labels( prm_results_set             );
	const auto value_list_labels           = get_sorted_instance_labels( prm_score_classn_value_list );

	// If value_list's instance labels aren't a subset of the results_set's instance labels, throw
	if ( ! includes( results_set_instance_labels, value_list_labels ) ) {
		const auto spurious_instance_labels = set_difference_build<str_set>( value_list_labels, results_set_instance_labels );
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot add score_classn_value_list to score_classn_value_results_set (even by adding missing instance labels) because it contains extra instance labels, eg: "
			+ join( spurious_instance_labels | limited( 6_z ), ", " )
		));
	}

	// Create a list of any result_set instance_labels that the value_list is missing
	const auto missing_instance_labels = set_difference_build<str_set>(
		results_set_instance_labels,
		value_list_labels
	);

	// If there are missing labels...
	if ( ! missing_instance_labels.empty() ) {

		// Emit a warning if warnings were requested
		if ( prm_warn_on_missing ) {
			::spdlog::warn( "Having to add {} missing instances to score_classn_value_list before adding to "
			                "score_classn_value_results_set; examples include: {}",
			                missing_instance_labels.size(),
			                join( missing_instance_labels | limited( 6_z ), ", " ) );
		}

		// For each entry in the results_set that's in the missing list, add an equivalent to the value_list
		for (const score_classn_value &the_value : front( prm_results_set ) ) {
			const auto &instance_label = the_value.get_instance_label();
			if ( contains( missing_instance_labels, instance_label ) ) {
				prm_score_classn_value_list.add_score_classn_value( score_classn_value(
					prm_score_for_missing,
					the_value.get_instance_is_positive(),
					instance_label
				) );
			}
		}
	}

	// Use the standard results_set method to add the now complete value_list
	prm_results_set.add_score_classn_value_list( prm_score_classn_value_list );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
str_set cath::score::get_sorted_instance_labels(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                                                ) {
	if ( prm_results_set.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get sorted instance labels of an empty score_classn_value_results_set"));
	}
	return get_sorted_instance_labels( front( prm_results_set ) );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
str_vec cath::score::get_instance_labels(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                                         ) {
	if ( prm_results_set.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return get_instance_labels( front( prm_results_set ) );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
size_t cath::score::get_num_instances(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                                      ) {
	if ( prm_results_set.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return ( front( prm_results_set ) ).size();
}


/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
string cath::score::get_name_of_index(const score_classn_value_results_set &prm_results_set, ///< TODOCUMENT
                                      const size_t                         &prm_index        ///< TODOCUMENT
                                      ) {
	return prm_results_set[ prm_index ].get_name();
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
str_vec cath::score::get_names(const score_classn_value_results_set &prm_results_set ///< TODOCUMENT
                               ) {
	return transform_build<str_vec>(
		indices( prm_results_set.size() ),
		[&] (const size_t &x) {
			return get_name_of_index( prm_results_set, x );
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
const bool & cath::score::get_higher_is_better_of_index(const score_classn_value_results_set &prm_results_set, ///< TODOCUMENT
                                                        const size_t                         &prm_index        ///< TODOCUMENT
                                                        ) {
	return get_higher_is_better( prm_results_set[ prm_index ] );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
const score_classn_value_list & cath::score::find_score_classn_value_list_of_name(const score_classn_value_results_set &prm_results_set, ///< TODOCUMENT
                                                                                  const string                         &prm_name         ///< TODOCUMENT
                                                                                  ) {
	const auto find_itr = find_if(
		prm_results_set,
		[&] (const score_classn_value_list &x) { return ( x.get_name() == prm_name ); }
	);
	if ( find_itr == cath::common::cend( prm_results_set ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception( "No score_classn_value_list found with specified name " + prm_name ));
	}
	return *find_itr;
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
named_true_false_pos_neg_list_list cath::score::make_named_true_false_pos_neg_list_list(const score_classn_value_results_set &prm_score_classn_value_results_set ///< TODOCUMENT
                                                                                        ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return named_true_false_pos_neg_list_list{
		transform_build<named_true_false_pos_neg_list_vec>(
			indices( prm_score_classn_value_results_set.size() ),
			[&] (const size_t &x) {
				const auto &name   = get_name_of_index( prm_score_classn_value_results_set, x );
				const auto &scores = prm_score_classn_value_results_set.get_score_classn_value_list_of_name( name );

				// \todo Put this check in get_score_classn_value_list_of_name()
				if ( name != scores.get_name() ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("Mismatching names in score_classn_value_results_set"));
				}
				return make_named_true_false_pos_neg_list( scores );
			}
		)
	};
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
named_true_false_pos_neg_list_list cath::score::make_named_true_false_pos_neg_list_list(const score_classn_value_list_vec &prm_score_classn_value_lists ///< TODOCUMENT
                                                                                        ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return named_true_false_pos_neg_list_list{
		transform_build<named_true_false_pos_neg_list_vec>(
			prm_score_classn_value_lists,
			[] (const score_classn_value_list &x) {
				return make_named_true_false_pos_neg_list( x );
			}
		)
	};
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
classn_stat_pair_series_list cath::score::make_classn_stat_pair_series_list(const score_classn_value_results_set &prm_score_classn_value_results_set, ///< TODOCUMENT
                                                                            const classn_stat                    &prm_classn_stat_a,                  ///< TODOCUMENT
                                                                            const classn_stat                    &prm_classn_stat_b                   ///< TODOCUMENT
                                                                            ) {
	return make_classn_stat_pair_series_list(
		make_named_true_false_pos_neg_list_list(
			prm_score_classn_value_results_set
		),
		prm_classn_stat_a,
		prm_classn_stat_b
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
classn_stat_pair_series_list cath::score::make_classn_stat_pair_series_list(const score_classn_value_list_vec &prm_score_classn_value_lists, ///< TODOCUMENT
                                                                            const classn_stat                 &prm_classn_stat_a,            ///< TODOCUMENT
                                                                            const classn_stat                 &prm_classn_stat_b             ///< TODOCUMENT
                                                                            ) {
	return make_classn_stat_pair_series_list(
		make_named_true_false_pos_neg_list_list(
			prm_score_classn_value_lists
		),
		prm_classn_stat_a,
		prm_classn_stat_b
	);
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
classn_stat_pair_series_list cath::score::make_roc_series_list(const score_classn_value_results_set &prm_score_classn_value_results_set ///< TODOCUMENT
                                                               ) {
	return make_standard_classn_stat_pair_series_list<roc_rates>( prm_score_classn_value_results_set );
}

/// \brief TODOCUMENT
///
/// \relates score_classn_value_results_set
classn_stat_pair_series_list cath::score::make_precision_recall_series_list(const score_classn_value_results_set &prm_score_classn_value_results_set ///< TODOCUMENT
                                                                            ) {
	return make_standard_classn_stat_pair_series_list<precision_recall_rates>( prm_score_classn_value_results_set );
}
