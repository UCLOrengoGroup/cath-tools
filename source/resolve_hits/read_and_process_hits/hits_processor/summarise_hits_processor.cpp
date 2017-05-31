/// \file
/// \brief The summarise_hits_processor class definitions

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

#include "summarise_hits_processor.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include "common/boost_addenda/range/front.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "resolve_hits/full_hit_fns.hpp"
#include "resolve_hits/full_hit_list.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::rslv::detail;
using namespace std::literals::string_literals;

using boost::format;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::max_element;
using boost::range::min_element;
using std::move;
using std::ostream;
using std::right;
using std::setw;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<hits_processor> summarise_hits_processor::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void summarise_hits_processor::do_process_hits_for_query(const string           &arg_query_id,         ///< The query_protein_id string
                                                         const crh_filter_spec  &/*arg_filter_spec*/,  ///< The filter_spec to apply to the hits
                                                         const crh_score_spec   &/*arg_score_spec*/,   ///< The score spec to apply to the hits
                                                         const crh_segment_spec &/*arg_segment_spec*/, ///< The segment spec to apply to the hits
                                                         const calc_hit_list    &arg_calc_hits         ///< The hits to process
                                                         ) {
	const full_hit_list &full_hits = arg_calc_hits.get_full_hits();
	if ( ! example_query_id_and_hit && ! full_hits.empty() ) {
		example_query_id_and_hit = make_pair(
			arg_query_id,
			front( full_hits )
		);
	}

	const auto max_stop_opt = get_max_stop( full_hits );
	if ( max_stop_opt ) {
		max_stops.push_back( *max_stop_opt );
	}

	num_hits += full_hits.size();
}

/// \brief Calculate the median of an unsorted bunch of size_t values
///
/// \TODO Move this into common/
double median(size_vec &args ///< The bunch of unsorted size_t values for which to calculate the median
              ) {
	const size_t size      = args.size();
	const size_t half_size = size / 2;
	std::nth_element(
		args.begin(),
		args.begin() + numeric_cast<ptrdiff_t>( half_size ),
		args.end()
	);
	const auto halfway_arg = numeric_cast<double>( args[ half_size ] );
	if ( size % 2 == 1) {
		return halfway_arg;
	}

	std::nth_element(
		args.begin(),
		args.begin() + numeric_cast<ptrdiff_t>( half_size ) - 1,
		args.end()
	);
	return 0.5 * ( halfway_arg + numeric_cast<double>( args[ half_size - 1 ] ) );
}

/// \brief Write the HTML suffix to finish the work (if it has been started)
void summarise_hits_processor::do_finish_work() {
	if ( num_hits > 0 ) {
		for (const ostream_ref &ostream_ref : get_ostreams() ) {
			ostream_ref.get()
				<< "Summary of input data\n"
				<< "---------------------\n"
				<< " * Number of queries : " << right << setw( 6 ) << max_stops.size() << " (excludes any queries with no hits)\n"
				<< " * Number of hits    : " << right << setw( 6 ) << num_hits
				<< " (ie an average of " << ( numeric_cast<double>( num_hits ) / numeric_cast<double>( max_stops.size() ) ) << " per query)\n"
				<< " * Minimum max-stop  : " << right << setw( 6 ) << ( max_stops.empty() ? "<N/A>" : std::to_string( *min_element( max_stops ) ) ) << "\n"
				<< " * Median  max-stop  : "
				<< (
					max_stops.empty()
						? " <N/A>"
						: ( format( "%8.1f" ) % median( max_stops ) ).str()
				)
				<< "\n"
				<< " * Maximum max-stop  : " << right << setw( 6 ) << ( max_stops.empty() ? "<N/A>" : std::to_string( *max_element( max_stops ) ) ) << "\n"
				<< " * Example hit       :\n"
				<< (
					example_query_id_and_hit
					?
						  "    * Query ID : " + example_query_id_and_hit->first                         + "\n"
						+ "    * Match ID : " + example_query_id_and_hit->second.get_label()            + "\n"
						+ "    * Score    : " + get_score_string   ( example_query_id_and_hit->second ) + "\n"
						+ "    * Segments : " + get_segments_string( example_query_id_and_hit->second ) + "\n"
					:
						""
				);
		}
	}
}

/// \brief Return true: read_and_resolve_mgr should still parse hits that fail the score filter and pass them to this processor
bool summarise_hits_processor::do_wants_hits_that_fail_score_filter() const {
	return true;
}

/// \brief Ctor for the summarise_hits_processor
summarise_hits_processor::summarise_hits_processor(ref_vec<ostream> arg_ostreams ///< The ostream to which the results should be written
                                                   ) noexcept : super{ move( arg_ostreams ) } {
}
