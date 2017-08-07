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
#include <boost/range/irange.hpp>

#include "common/algorithm/transform_build.hpp"

using namespace cath;
using namespace cath::clust;
using namespace cath::common;

#include <numeric>

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::format;
using boost::irange;
using std::accumulate;
using std::make_pair;
using std::max;
using std::min;
using std::string;

constexpr size_t overlap_frac_distn::num_dec_places;
constexpr size_t overlap_frac_distn::num_gaps;
constexpr size_t overlap_frac_distn::num_posns;

/// \brief Find the index of the bin containing the i-th overlap fraction for specified i (where the overlap fractions are in ascending order)
size_t overlap_frac_distn::find_index_of_nth(const size_t &arg_n ///< The index of the overlap fraction required
                                             ) const {
	size_t sum = 0;

	// \TODO Come C++17 and structured bindings, use here
	for (const boost::tuple<const size_t &, size_t> &value_and_index : combine( *this, irange( 0_z, num_posns ) ) ) {
		const size_t &value = value_and_index.get<0>();
		const size_t &index = value_and_index.get<1>();

		if ( value > 0 ) {
			sum += value;
			if ( sum > arg_n ) {
				return index;
			}
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find n values in overlap_frac_distn"));
}

/// \brief Get the number of overlap fractions that are in the specified [begin, end) half-open range
size_t overlap_frac_distn::get_num_in_range(const double &arg_lower_fraction, ///< The lower fraction, which is inclusive
                                            const double &arg_upper_fraction  ///< The upper fraction, which is exclusive and can go over 1
                                            ) const {
	for (const double &fraction : { arg_lower_fraction, arg_upper_fraction  } ) {
		if ( ! boost::math::isfinite( fraction ) || fraction < 0.0 ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get num from overlap_frac_distn in range with invalid fraction"));
		}
	}
	const auto lower_index = static_cast<size_t>( round( get_index_of_fraction( arg_lower_fraction ) ) );
	const auto upper_index = static_cast<size_t>( round( get_index_of_fraction( arg_upper_fraction ) ) );

	return accumulate(
		next( common::cbegin( *this ), static_cast<ptrdiff_t>( min( num_posns, lower_index ) ) ),
		next( common::cbegin( *this ), static_cast<ptrdiff_t>( min( num_posns, upper_index ) ) ),
		0_z
	);
}

/// \brief Get the overlap fraction to be found at the specified percentile
///
/// This uses the nearest-rank percentile (that does no interpolation and just takes the value at the i-th place
/// where i = ceil( P / 100 * N )
double overlap_frac_distn::get_frac_at_percentile(const double        &arg_percentile,   ///< The percentile to look for
                                                  const zeroes_policy &arg_zeroes_policy ///< Whether to exclude any zero overlap fractions from calculations
                                                  ) const {
	const size_t zeroes_offset = ( arg_zeroes_policy == zeroes_policy::EXCLUDE )
	                             ? fraction_counts.front()
	                             : 0_z;
	const size_t num_fractions = size() - zeroes_offset;
	if ( num_fractions == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get_frac_at_percentile on an unpopulated overlap_frac_distn (after zeroes removed if requested)"));
	}

	const size_t nth = max( 1_z, static_cast<size_t>( ceil( arg_percentile / 100.0 * static_cast<double>( num_fractions ) ) ) ) - 1_z;

	return get_fraction_of_index( find_index_of_nth( nth + zeroes_offset ) );
}

/// \brief Build a overlap_frac_distn from the specified list of overlap fractions
///
/// \relates overlap_frac_distn
overlap_frac_distn cath::clust::build_overlap_frac_distn_from_overlap_fractions(const doub_vec &arg_overlap_fractions ///< The list of overlap fractions from which to populate the new overlap_frac_distn
                                                                                ) {
	overlap_frac_distn result;
	for (const double &overlap_fraction : arg_overlap_fractions) {
		result.add_overlap_fraction( overlap_fraction );
	}
	return result;
}

/// \brief Generate percentile data for the specified overlap_frac_distn at the specified percentile points
///
/// \relates overlap_frac_distn
doub_doub_pair_vec cath::clust::percentile_data(const overlap_frac_distn                &arg_overlap_frac_distn, ///< The overlap_frac_distn to query
                                                const doub_vec                          &arg_percentiles,        ///< The percentile points to calculate values for (eg median is 50.0)
                                                const overlap_frac_distn::zeroes_policy &arg_zeroes_policy       ///< Whether to exclude any zero overlap fractions from calculations
                                                ) {
	return transform_build<doub_doub_pair_vec>(
		arg_percentiles,
		[&] (const double &percentile) -> doub_doub_pair {
			return { percentile, 100.0 * arg_overlap_frac_distn.get_frac_at_percentile( percentile, arg_zeroes_policy ) };
		}
	);
}

/// \brief Generate a Markdown table of percentile data for the specified overlap_frac_distn at the specified percentile points
///
/// \relates overlap_frac_distn
string cath::clust::percentile_markdown_table(const overlap_frac_distn                &arg_overlap_frac_distn, ///< The overlap_frac_distn to query
                                              const doub_vec                          &arg_percentiles,        ///< The percentile points to calculate values for (eg median is 50.0)
                                              const string                            &arg_percentile_title,   ///< The title for the percentiles column
                                              const string                            &arg_value_title,        ///< The title for the values column
                                              const overlap_frac_distn::zeroes_policy &arg_zeroes_policy       ///< Whether to exclude any zero overlap fractions from calculations
                                              ) {
	using std::to_string;

	const doub_doub_pair_vec data = percentile_data( arg_overlap_frac_distn, arg_percentiles, arg_zeroes_policy );
	return "| " + arg_percentile_title + " | " + arg_value_title + " |\n"
		+ "|" + string( 2 + arg_percentile_title.length(), '-' )
		// + "--------------------------------------------------------"
		+ "|" + string( 2 + arg_value_title.length(),      '-' )
		+ "|"
		+ join(
			data
				| transformed( [&] (const doub_doub_pair &x) {
					return
						  "\n| "
						+ ( format( R"(%)" + to_string( arg_percentile_title.length() ) + "d" ) % x.first  ).str()
						+ " |"
						+ ( format( R"(%)" + to_string( arg_value_title.length()      ) + "d" ) % x.second ).str()
						+ R"(% |)";
				} ),
			""
		);
}
