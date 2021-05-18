/// \file
/// \brief The filter_vs_full_score_list class definitions

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

#include "filter_vs_full_score_list.hpp"

#include <algorithm>
#include <filesystem>
#include <numeric>

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/sub_range.hpp>

#include <gnuplot-iostream.h>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "cath/common/boost_addenda/sorted_insert.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/score/true_pos_false_neg/classn_rate_stat.hpp"
#include "cath/score/true_pos_false_neg/true_false_pos_neg.hpp"
#include "cath/structure/view_cache/filter/detail/filter_vs_full_score_less.hpp"
#include "cath/structure/view_cache/filter/filter_vs_full_score.hpp"

/// \todo Write some decent unit tests for sorted_insert and then remove this code
#ifndef NDEBUG
#include <boost/algorithm/cxx11/is_sorted.hpp>
using ::boost::algorithm::is_sorted;
#endif

using namespace ::cath::common;
using namespace ::cath::index::filter;
using namespace ::cath::index::filter::detail;
using namespace ::cath::score;

using ::boost::accumulate;
using ::boost::distance;
using ::boost::numeric_cast;
using ::boost::range::lower_bound;
using ::boost::range::nth_element;
using ::boost::range::sort;
using ::boost::rational_cast;
using ::boost::sub_range;
using ::std::cend;
using ::std::cerr;
using ::std::end;
using ::std::endl;
using ::std::filesystem::path;
using ::std::make_pair;
using ::std::to_string;

/// \brief Sort the filter_vs_full_score objects by full_score (ascending) to uphold that class invariant
void filter_vs_full_score_list::sort_filter_vs_full_scores() {
	sort( filter_vs_full_scores, full_score_less() );
}

/// \brief Ctor from a vector of filter_vs_full_score objects
///
/// The vector specified in the argument does not need to be pre-sorted;
/// once this ctor has made the local copy, it will sort that copy
filter_vs_full_score_list::filter_vs_full_score_list(filter_vs_full_score_vec prm_filter_vs_full_scores ///< The vector of filter_vs_full_score objects from which to construct this filter_vs_full_score_list
                                                     ) : filter_vs_full_scores( std::move( prm_filter_vs_full_scores ) ) {
	sort_filter_vs_full_scores();
}

/// \brief Add a filter_vs_full_score to the list
///
/// This will resort the filter_vs_full_score objects once the new one has been added
void filter_vs_full_score_list::add_filter_vs_full_score(const filter_vs_full_score &prm_filter_vs_full_score ///< The filter_vs_full_score to be added to the list
                                                         ) {
	// Insert into the correct place in the sorted vector
	sorted_insert( filter_vs_full_scores, prm_filter_vs_full_score, full_score_less() );

	/// \todo Write some decent unit tests for sorted_insert and then remove this check
#ifndef NDEBUG
	if ( ! is_sorted( filter_vs_full_scores, full_score_less() ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Attempt to use sorted_insert failed"));
	}
#endif
}

/// \brief Standard size method()
size_t filter_vs_full_score_list::size() const {
	return filter_vs_full_scores.size();
}

/// \brief Standard subscript operator
const filter_vs_full_score & filter_vs_full_score_list::operator[](const size_t &prm_index ///< The index of the filter_vs_full_score to access
                                                                   ) const {
	return filter_vs_full_scores[ prm_index ];
}

/// \brief Standard const begin() method as part of making filter_vs_full_score_list into a (const) range
filter_vs_full_score_list::const_iterator filter_vs_full_score_list::begin() const {
	return cbegin( filter_vs_full_scores );
}

/// \brief Standard const end() method as part of making filter_vs_full_score_list into a (const) range
filter_vs_full_score_list::const_iterator filter_vs_full_score_list::end() const {
	return cend( filter_vs_full_scores );
}

/// \brief Find the minimum filter score that achieves the specified sensitivity re finding entries with the specified full score
///
/// \relates filter_vs_full_score_list
///
/// \todo Test that the resulting score does indeed achieve the requested sensitivity
double cath::index::filter::filter_score_full_score_with_sensitivity(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                     const double                    &prm_full_score,                ///< The score that defines which filter_vs_full_scores are "wanted" (those which have a full_score >= this score)
                                                                     const double                    &prm_sensitivity_fraction       ///< The required sensitivity that the filter score must achieve
                                                                     ) {
	// Find an iterator to the first entry with a full_score >= the specified full score in the filter_vs_full_score_list
	// (which is already sorted by full_score)
	const filter_vs_full_score_list_citr begin_of_wanted = lower_bound(
		prm_filter_vs_full_score_list,
		prm_full_score,
		full_score_less()
	);

	// Construct a sub_range for the values after that point (ie all those with a full_score >= the specified full score),
	// grab the number of such "wanted" elements and throw if it's zero
	const sub_range<filter_vs_full_score_list> wanted_range( begin_of_wanted, cend( prm_filter_vs_full_score_list ) );
	const size_t num_wanted = numeric_cast<size_t>( distance( wanted_range ) );
	if ( num_wanted == 0 ) {
		cerr << "Last entry is : " << *prev( end(prm_filter_vs_full_score_list ) ) << endl;
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot find filter_score_full_score_with_sensitivity because no filter_vs_full_scores match specified full_score "
			+ to_string( prm_full_score )
		));
	}

	// Grab the filter scores of all wanted values
	auto filter_scores_of_wanted = transform_build<doub_vec>(
		wanted_range,
		[] (const filter_vs_full_score &x) { return x.get_filter_score(); }
	);

	// Aim is to find the filter value that at least prm_sensitivity_fraction of the "wanted" filter values
	// are greater than equal to
	//
	// Start by calculating the maximum number of "wanted" elements that can be less than the final value
	const double miss_rate       = 1.0 - prm_sensitivity_fraction;
	const size_t fraction_offset = numeric_cast<size_t>( floor( miss_rate * numeric_cast<double>( num_wanted ) ) );

	// Use nth_element to find the filter score associated with the
	// entry in the fraction_offset position.
	//
	// (nth_element is a std algorithm designed specifically for precisely this sort of percentile-type
	//  calculation and it achieves it in (average) linear time (whereas sort requires linearithmic (ie O(n.log(n))) time)
	nth_element(
		filter_scores_of_wanted,
		std::begin( filter_scores_of_wanted ) + numeric_cast<ptrdiff_t>( fraction_offset )
	);
	return *( cbegin( filter_scores_of_wanted ) + numeric_cast<ptrdiff_t>( fraction_offset ) );
}

/// \brief Find the minimum filter score that achieves the specified sensitivity re finding entries with the specified full score
///
/// \returns A filter_vs_full_score containing the filter score and the original full score
///
/// \relates filter_vs_full_score_list
///
/// \todo Test that the resulting score does indeed achieve the requested sensitivity
filter_vs_full_score cath::index::filter::filter_attempt_full_score_with_sensitivity(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                                     const double                    &prm_full_score,                ///< The score that defines which filter_vs_full_scores are "wanted" (those which have a full_score >= this score)
                                                                                     const double                    &prm_sensitivity_fraction       ///< The required sensitivity that the filter score must achieve
                                                                                     ) {
	return filter_vs_full_score(
		filter_score_full_score_with_sensitivity(
			prm_filter_vs_full_score_list,
			prm_full_score,
			prm_sensitivity_fraction
		),
		prm_full_score
	);
}

/// \brief Find the minimum filter score that achieves the specified sensitivity re finding entries with the specified full score
///        and then assess the result of applying that filter on the data
///
/// \returns A true_false_pos_neg describing results of applying the calculated filter score to find entries with the specified full
///          score
///
/// \relates filter_vs_full_score_list
///
/// \todo Test that the resulting score does indeed achieve the requested sensitivity
true_false_pos_neg cath::index::filter::filter_result_full_score_with_sensitivity(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                                  const double                    &prm_full_score,                ///< The score that defines which filter_vs_full_scores are "wanted" (those which have a full_score >= this score)
                                                                                  const double                    &prm_sensitivity_fraction       ///< The required sensitivity that the filter score must achieve
                                                                                  ) {
	const true_false_pos_neg results = assess_results_on_filter_attempt(
		prm_filter_vs_full_score_list,
		filter_attempt_full_score_with_sensitivity(
			prm_filter_vs_full_score_list,
			prm_full_score,
			prm_sensitivity_fraction
		)
	);
	const double the_sensitivity = rational_cast<double>( sensitivity().calculate( results ) );
	if ( the_sensitivity < prm_sensitivity_fraction ) {

		cerr <<
			"Attempt to calculate filter for sensitivity of "
			<< to_string( prm_sensitivity_fraction )
			<< " resulted in an inadequate sensitivity of "
			<< to_string( the_sensitivity ) << endl;
//		));
//		BOOST_THROW_EXCEPTION(out_of_range_exception(
//			"Attempt to calculate filter for sensitivity of "
//			+ to_string( prm_sensitivity_fraction )
//			+ " resulted in an inadequate sensitivity of "
//			+ to_string( the_sensitivity )
//		));
	}
	return results;
}

/// \brief TODOCUMENT
///
/// \relates filter_vs_full_score_list
filter_vs_full_score_vec cath::index::filter::filter_attempts_with_sensitivity(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                               const double                    &prm_sensitivity_fraction       ///< The score that defines which filter_vs_full_scores are "wanted" (those which have a full_score >= this score)
                                                                               ) {
	filter_vs_full_score_vec results;
	for (const filter_vs_full_score &real_score : prm_filter_vs_full_score_list) {
		const double &full_score = real_score.get_full_score();
		if ( results.empty() || results.back().get_full_score() != full_score ) {
			// Calculate the result for full_score
			const filter_vs_full_score result = filter_attempt_full_score_with_sensitivity(
				prm_filter_vs_full_score_list,
				full_score,
				prm_sensitivity_fraction
			);

			// Add the result to results
			results.push_back( result );
		}
	}
	return results;
}

/// \brief Calculate the true_false_pos_neg results associated with filters calculated to achieve the specified sensitivity
///        for each of the full_scores in the filter_vs_full_score_list
///
/// \relates filter_vs_full_score_list
doub_true_false_pos_neg_pair_vec cath::index::filter::filter_results_with_sensitivity(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                                      const double                    &prm_sensitivity_fraction       ///< The score that defines which filter_vs_full_scores are "wanted" (those which have a full_score >= this score)
                                                                                      ) {
	const filter_vs_full_score_vec filter_attempts = filter_attempts_with_sensitivity(
		prm_filter_vs_full_score_list,
		prm_sensitivity_fraction
	);

	doub_true_false_pos_neg_pair_vec results;
	results.reserve( prm_filter_vs_full_score_list.size() );
	for (const filter_vs_full_score &filter_attempt : filter_attempts) {
		const double             &full_score = filter_attempt.get_full_score();
		const true_false_pos_neg  result     = assess_results_on_filter_attempt(
			prm_filter_vs_full_score_list,
			filter_attempt
		);

		// Add the result to results
		results.push_back( make_pair( full_score, result ) );
	}

	return results;
}

/// \brief Sum up the true_false_pos_neg results of assessing all the specified filter_vs_full_score_list results
///        against the specified filter attempt (which tries to identify all full scores >= some value by
///        selecting all filter scores >= some value)
///
/// \relates filter_vs_full_score_list
true_false_pos_neg cath::index::filter::assess_results_on_filter_attempt(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                                                         const filter_vs_full_score      &prm_filter_attempt             ///< The filter attempt (interpretation: this attempts to identify all entries with full score >= filter_vs_full_score's full score by selecting all entries with filter score >= filter_vs_full_score's filter score)
                                                                         ) {
	/// \todo Come C++14 (GCC >= 4.9), change the lambda to use [] (const auto &x) { ... }
	return accumulate(
		prm_filter_vs_full_score_list,
		true_false_pos_neg(),
		[&] (const true_false_pos_neg &x, const filter_vs_full_score &y) {
			return x + assess_real_scores_on_filter_attempt( y, prm_filter_attempt );
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates filter_vs_full_score_list
void cath::index::filter::gnuplot_data(const filter_vs_full_score_list &prm_filter_vs_full_score_list, ///< The real scores to be assessed
                                       const path                      &prm_output_stem,               ///< The stem of the file to generate (ie without the ".gnuplot" or ".eps" suffix)
                                       const filter_vs_full_score_list &prm_filter_attempts            ///< The filter policy to be plotted
                                       ) {
	const path gnuplot_file     = replace_extension_copy( prm_output_stem, ".gnuplot"         );
	const path eps_file         = replace_extension_copy( prm_output_stem, ".eps"             );
	const path the_data_file    = replace_extension_copy( prm_output_stem, ".data.txt"        );
	const path filter_data_file = replace_extension_copy( prm_output_stem, ".filter_data.txt" );
	Gnuplot gp("tee " + gnuplot_file.string() + " | gnuplot"); // Write to an intermediate gnuplot file

	gp << "set   terminal postscript color\n";
	gp << "set   output " << eps_file << "\n";
//	gp << "set   size square\n";

//	gp << "set   xtics 0,10\n";
//	gp << "set   ytics 0,10\n";
	gp << "set   xtics font \"Helvetica,10\"\n";
	gp << "set   ytics font \"Helvetica,10\"\n";
	gp << "set   style line 11 lc rgb '#808080' lt 1\n";
	gp << "set   border 3 back ls 11\n";
	gp << "set   tics nomirror\n";
	gp << "set   style line 12 lc rgb '#808080' lt 0 lw 1\n";
	gp << "set   style line 1 lt rgb \"#000000\"\n";
	gp << "set   grid back ls 12\n";

	gp << "set   title \"Filter scores versus real scores (SSAP scores)\"\n";
	gp << "set   xlabel \"Real score (SSAP score)\"\n";
	gp << "set   ylabel \"Filter score\"\n";

	const auto values = transform_build<doub_doub_pair_vec>(
		prm_filter_vs_full_score_list,
		[] (const filter_vs_full_score &x) {
			return make_pair( x.get_full_score(), x.get_filter_score() );
		}
	);
	const auto filters = transform_build<doub_doub_pair_vec>(
		prm_filter_attempts,
		[] (const filter_vs_full_score &x) {
			return make_pair( x.get_full_score(), x.get_filter_score() );
		}
	);

	gp << "plot "
	   << gp.file1d( values,  the_data_file.string()    ) << " with points notitle linetype rgb \"#555555\" pointtype 7 pointsize 0.1, "
//	   << gp.file1d( filters, filter_data_file.string() ) << " with points notitle linetype rgb \"#000000\" pointtype 1 pointsize 0.5"
	   << gp.file1d( filters, filter_data_file.string() ) << " with lines  notitle linestyle 1\n";
	// To check: has removing the final endl from here caused any problems?
}

/// \brief TODOCUMENT
///
/// \relates filter_vs_full_score_list
void cath::index::filter::gnuplot_classsn_stat_for_recall(const filter_vs_full_score_list &prm_data,        ///< TODOCUMENT
                                                          const path                      &prm_output_stem, ///< TODOCUMENT
                                                          const classn_stat               &prm_classn_stat, ///< TODOCUMENT
                                                          const double                    &prm_recall       ///< TODOCUMENT
                                                          ) {
	const doub_true_false_pos_neg_pair_vec results = filter_results_with_sensitivity( prm_data, prm_recall);
	gnuplot_classsn_stat_for_recall( results, prm_output_stem, prm_classn_stat );
}

/// \brief TODOCUMENT
///
/// \relates filter_vs_full_score_list
void cath::index::filter::gnuplot_classsn_stat_for_recall(const doub_true_false_pos_neg_pair_vec &prm_data,        ///< TODOCUMENT
                                                          const path                             &prm_output_stem, ///< TODOCUMENT
                                                          const classn_stat                      &prm_classn_stat  ///< TODOCUMENT
                                                          ) {
	const path gnuplot_file  = replace_extension_copy( prm_output_stem, ".gnuplot"  );
	const path eps_file      = replace_extension_copy( prm_output_stem, ".eps"      );
	const path the_data_file = replace_extension_copy( prm_output_stem, ".data.txt" );
	Gnuplot gp("tee " + gnuplot_file.string() + " | gnuplot"); // Write to an intermediate gnuplot file

	gp << "set   terminal postscript color\n";
	gp << "set   output " << eps_file << "\n";
//	gp << "set   size square\n";

//	gp << "set   xtics 0,10\n";
//	gp << "set   ytics 0,10\n";
	gp << "set   xtics font \"Helvetica,10\"\n";
	gp << "set   ytics font \"Helvetica,10\"\n";
	gp << "set   style line 11 lc rgb '#808080' lt 1\n";
	gp << "set   border 3 back ls 11\n";
	gp << "set   tics nomirror\n";
	gp << "set   style line 12 lc rgb '#808080' lt 0 lw 1\n";
	gp << "set   grid back ls 12\n";
	gp << "set   style line 1 lt rgb \"#A00000\" lw 2 pt 1\n";

	gp << "set   title \"Filter statistic versus real scores (SSAP scores)\"\n";
	gp << "set   xlabel \"Real score (SSAP score)\"\n";
	gp << "set   ylabel \"Statistic\"\n";

	const auto values = transform_build<doub_doub_pair_vec>(
		prm_data,
		[&] (const doub_true_false_pos_neg_pair &x) {
			return make_pair(
				x.first,
				rational_cast<double>( prm_classn_stat.calculate( x.second ) )
			);
		}
	);

	gp << "plot " << gp.file1d( values, the_data_file.string() ) << " with lines notitle linestyle 1\n";

	// To check: has removing the final endl from here caused any problems?
}

