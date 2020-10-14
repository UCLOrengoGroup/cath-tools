/// \map
/// \brief The overlap_frac_distn class definitions

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

#include "overlap_frac_distn.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/max_proj_element.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;

#include <numeric>
#include <tuple>

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::format;
using ::std::accumulate;
using ::std::get;
using ::std::max;
using ::std::min;
using ::std::string;

constexpr size_t overlap_frac_distn::num_dec_places;
constexpr size_t overlap_frac_distn::num_gaps;
constexpr size_t overlap_frac_distn::num_posns;

/// \brief Find the index of the bin containing the i-th overlap fraction for specified i (where the overlap fractions are in ascending order)
size_t overlap_frac_distn::find_index_of_nth(const size_t &prm_n ///< The index of the overlap fraction required
                                             ) const {
	size_t sum = 0;

	// \TODO Come C++17 and structured bindings, use here
	for (const boost::tuple<const size_t &, size_t> &value_and_index : combine( *this, indices( num_posns ) ) ) {
		const size_t &value = value_and_index.get<0>();
		const size_t &index = value_and_index.get<1>();

		if ( value > 0 ) {
			sum += value;
			if ( sum > prm_n ) {
				return index;
			}
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find n values in overlap_frac_distn"));
}

/// \brief Get the number of overlap fractions that are in the specified [begin, end) half-open range
size_t overlap_frac_distn::get_num_in_range(const double &prm_lower_fraction, ///< The lower fraction, which is inclusive
                                            const double &prm_upper_fraction  ///< The upper fraction, which is exclusive and can go over 1
                                            ) const {
	for (const double &fraction : { prm_lower_fraction, prm_upper_fraction  } ) {
		if ( ! boost::math::isfinite( fraction ) || fraction < 0.0 ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get num from overlap_frac_distn in range with invalid fraction"));
		}
	}
	const auto lower_index = static_cast<size_t>( round( get_index_of_fraction( prm_lower_fraction ) ) );
	const auto upper_index = static_cast<size_t>( round( get_index_of_fraction( prm_upper_fraction ) ) );

	return accumulate(
		next( common::cbegin( *this ), static_cast<ptrdiff_t>( min( num_posns, lower_index ) ) ),
		next( common::cbegin( *this ), static_cast<ptrdiff_t>( min( num_posns, upper_index ) ) ),
		0_z
	);
}

/// \brief Get the number at the specified fraction
size_t overlap_frac_distn::get_num_at_fraction(const double &prm_fraction ///< The fraction to query
                                               ) const {
	if ( ! boost::math::isfinite( prm_fraction ) || prm_fraction < 0.0 || prm_fraction > 1.0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get_num_at_fraction() from overlap_frac_distn with invalid fraction"));
	}
	return *next(
		common::cbegin( *this ),
		static_cast<ptrdiff_t>( round( get_index_of_fraction( prm_fraction ) ) )
	);
}

/// \brief Get the overlap fraction to be found at the specified percentile
///
/// This uses the nearest-rank percentile (that does no interpolation and just takes the value at the i-th place
/// where i = ceil( P / 100 * N )
double overlap_frac_distn::get_frac_at_percentile(const double        &prm_percentile,   ///< The percentile to look for
                                                  const zeroes_policy &prm_zeroes_policy ///< Whether to exclude any zero overlap fractions from calculations
                                                  ) const {
	const size_t zeroes_offset = ( prm_zeroes_policy == zeroes_policy::EXCLUDE )
	                             ? fraction_counts.front()
	                             : 0_z;
	const size_t num_fractions = size() - zeroes_offset;
	if ( num_fractions == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get_frac_at_percentile on an unpopulated overlap_frac_distn (after zeroes removed if requested)"));
	}

	const size_t nth = max( 1_z, static_cast<size_t>( ceil( prm_percentile / 100.0 * static_cast<double>( num_fractions ) ) ) ) - 1_z;

	return get_fraction_of_index( find_index_of_nth( nth + zeroes_offset ) );
}

/// \brief Get the number in the specified half open-closed range
///
/// \relates overlap_frac_distn
size_t cath::clust::get_num_in_open_closed_range(const overlap_frac_distn &prm_overlap_frac_distn, ///< The overlap_frac_distn to query
                                                 const double             &prm_lower_fraction,     ///< The lower fraction, which is exclusive
                                                 const double             &prm_upper_fraction      ///< The upper fraction, which is inclusive
                                                 ) {
	return prm_overlap_frac_distn.get_num_in_range( prm_lower_fraction, prm_upper_fraction )
		- prm_overlap_frac_distn.get_num_at_fraction( prm_lower_fraction )
		+ prm_overlap_frac_distn.get_num_at_fraction( prm_upper_fraction );
}

/// \brief Build a overlap_frac_distn from the specified list of overlap fractions
///
/// \relates overlap_frac_distn
overlap_frac_distn cath::clust::build_overlap_frac_distn_from_overlap_fractions(const doub_vec &prm_overlap_fractions ///< The list of overlap fractions from which to populate the new overlap_frac_distn
                                                                                ) {
	overlap_frac_distn result;
	for (const double &overlap_fraction : prm_overlap_fractions) {
		result.add_overlap_fraction( overlap_fraction );
	}
	return result;
}

/// \brief Generate percentile data for the specified overlap_frac_distn at the specified percentile points
///
/// \relates overlap_frac_distn
doub_doub_pair_vec cath::clust::percentile_data(const overlap_frac_distn                &prm_overlap_frac_distn, ///< The overlap_frac_distn to query
                                                const doub_vec                          &prm_percentiles,        ///< The percentile points to calculate values for (eg median is 50.0)
                                                const overlap_frac_distn::zeroes_policy &prm_zeroes_policy       ///< Whether to exclude any zero overlap fractions from calculations
                                                ) {
	return transform_build<doub_doub_pair_vec>(
		prm_percentiles,
		[&] (const double &percentile) -> doub_doub_pair {
			return { percentile, 100.0 * prm_overlap_frac_distn.get_frac_at_percentile( percentile, prm_zeroes_policy ) };
		}
	);
}

/// \brief Generate a Markdown table of percentile data for the specified overlap_frac_distn at the specified percentile points
///
/// \relates overlap_frac_distn
string cath::clust::percentile_markdown_table(const overlap_frac_distn                &prm_overlap_frac_distn, ///< The overlap_frac_distn to query
                                              const doub_vec                          &prm_percentiles,        ///< The percentile points to calculate values for (eg median is 50.0)
                                              const string                            &prm_percentile_title,   ///< The title for the percentiles column
                                              const string                            &prm_value_title,        ///< The title for the values column
                                              const overlap_frac_distn::zeroes_policy &prm_zeroes_policy       ///< Whether to exclude any zero overlap fractions from calculations
                                              ) {
	using ::std::to_string;

	const doub_doub_pair_vec data = percentile_data( prm_overlap_frac_distn, prm_percentiles, prm_zeroes_policy );
	return "| " + prm_percentile_title + " | " + prm_value_title + " |\n"
		+ "|" + string( 2 + prm_percentile_title.length(), '-' )
		// + "--------------------------------------------------------"
		+ "|" + string( 2 + prm_value_title.length(),      '-' )
		+ "|"
		+ join(
			data
				| transformed( [&] (const doub_doub_pair &x) {
					return
						  "\n| "
						+ ( format( R"(%)" + to_string( prm_percentile_title.length() ) + "d" ) % x.first  ).str()
						+ " |"
						+ ( format( R"(%)" + to_string( prm_value_title.length()      ) + "d" ) % x.second ).str()
						+ R"(% |)";
				} ),
			""
		);
}

/// \brief Get histogram data for the specified overlap_frac_distn given the specified
///        number of domains that mapped at 0% with nothing on the parent sequence
///
/// \relates overlap_frac_distn
str_size_doub_tpl_vec cath::clust::histogram_data(const overlap_frac_distn &prm_overlap_frac_distn,        ///< The overlap_frac_distn to query
                                                  const size_t             &prm_num_with_nothing_on_parent ///< The number of domains that mapped at 0% with nothing on the parent sequence
                                                  ) {
	using ::std::to_string;

	const size_t total       = prm_overlap_frac_distn.size();
	const size_t num_at_zero = prm_overlap_frac_distn.get_num_at_fraction( 0.0 );

	if ( prm_num_with_nothing_on_parent > num_at_zero || num_at_zero > total ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The number at zero is incompatible with the prm_num_with_nothing_on_parent or total"));
	}

	str_size_pair_vec results = { {
		str_size_pair{ R"(All)",                                        total                                        },
		             { R"(0% (with no other match on the sequence))",   prm_num_with_nothing_on_parent               },
		             { R"(0% (with some other match on the sequence))", num_at_zero - prm_num_with_nothing_on_parent },
	} };

	for (const size_t &frac_tenth : indices( 10_z ) ) {
		const size_t range_pc_begin   = 10 *       frac_tenth;
		const size_t range_pc_end     = 10 * ( 1 + frac_tenth );
		const double range_frac_begin = static_cast<double>( range_pc_begin ) / 100.0;
		const double range_frac_end   = static_cast<double>( range_pc_end   ) / 100.0;
		results.emplace_back(
			to_string( range_pc_begin ) + R"(% < x <= )" + to_string( range_pc_end ) + R"(%)",
			get_num_in_open_closed_range( prm_overlap_frac_distn, range_frac_begin, range_frac_end )
		);
	}

	return transform_build<str_size_doub_tpl_vec>(
		results,
		[&] (const str_size_pair &x) {
			return make_tuple( x.first, x.second, 100.0 * static_cast<double>( x.second ) / static_cast<double>( total ) );
		}
	);
}

/// \brief Generate a Markdown histogram string for the specified overlap_frac_distn
///
/// \relates overlap_frac_distn
string cath::clust::histogram_markdown_table(const overlap_frac_distn &prm_overlap_frac_distn,        ///< The overlap_frac_distn to query
                                             const string             &prm_range_title,               ///< The title for the range column
                                             const string             &prm_value_title,               ///< The title for the values column
                                             const string             &prm_percent_title,             ///< The title for the percent column
                                             const size_t             &prm_num_with_nothing_on_parent ///< The number of domains that mapped at 0% with nothing on the parent sequence
                                             ) {
	using ::std::to_string;

	const auto data = histogram_data( prm_overlap_frac_distn, prm_num_with_nothing_on_parent );

	const size_t longest_first_col_length = data.empty()
		? 0_z
		: max_proj(
			data,
			std::less<>{},
			[] (const str_size_doub_tpl &x) {
				return get<0>( x ).length();
			}
		);

	return
		  "| "  + ( format( R"(%)" + to_string( longest_first_col_length ) + "s"   ) % prm_range_title ).str()
		+ " | " + prm_value_title
		+ " | " + prm_percent_title
		+ " |\n"
		+ "|-"  + string( longest_first_col_length,   '-' )
		+ "-|-" + string( prm_value_title.length(),   '-' )
		+ "-|-" + string( prm_percent_title.length(), '-' )
		+ "-|"
		+ join(
			data
				| transformed( [&] (const str_size_doub_tpl &x) {
					return
						  "\n| "
						+ ( format( R"(%)" + to_string( longest_first_col_length   ) + "s"   ) % get<0>( x ) ).str()
						+ " | "
						+ ( format( R"(%)" + to_string( prm_value_title.length()   ) + "d"   ) % get<1>( x ) ).str()
						+ " |"
						+ ( format( R"(%)" + to_string( prm_percent_title.length() ) + ".3g" ) % get<2>( x ) ).str()
						+ R"(% |)";
				} ),
			""
		);
}
