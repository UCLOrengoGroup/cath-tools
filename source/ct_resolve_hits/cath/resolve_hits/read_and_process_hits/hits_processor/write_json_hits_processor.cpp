/// \file
/// \brief The write_json_hits_processor class definitions

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

#include "write_json_hits_processor.hpp"

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"
#include "cath/resolve_hits/resolve/hit_resolver.hpp"
#include "cath/resolve_hits/scored_hit_arch.hpp"

using namespace cath::common;
using namespace cath::rslv::detail;

using std::move;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<hits_processor> write_json_hits_processor::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void write_json_hits_processor::do_process_hits_for_query(const string           &prm_query_id,        ///< The query_protein_id string
                                                          const crh_filter_spec  &/*prm_filter_spec*/, ///< The filter_spec to apply to the hits
                                                          const crh_score_spec   &prm_score_spec,      ///< The score spec to apply to the hits
                                                          const crh_segment_spec &prm_segment_spec,    ///< The segment spec to apply to the hits
                                                          const calc_hit_list    &prm_calc_hits        ///< The hits to process
                                                          ) {
	if ( ! has_started ) {
		json_writers.start_object();
		has_started = true;
	}

	// Resolve the hits
	const auto result_hit_arch  = resolve_hits( prm_calc_hits, prm_score_spec.get_naive_greedy() );
	const auto result_full_hits = get_full_hits_of_hit_arch(
		result_hit_arch,
		prm_calc_hits.get_full_hits()
	);

	// Output the results to prm_ostream
	json_writers.write_key( prm_query_id );
	json_writers.write_raw_string( to_json_string_with_compact_fullhits( result_full_hits, prm_segment_spec, 1 ) );
}

/// \brief Do nothing to finish the batch of work
void write_json_hits_processor::do_finish_work() {
	if ( has_started && ! json_writers.is_complete() ) {
		json_writers.end_object();
	}
}

/// \brief Return false: read_and_resolve_mgr needn't parse hits that fail the score filter or pass them to this processor
bool write_json_hits_processor::do_wants_hits_that_fail_score_filter() const {
	return false;
}

/// \brief Return false: read_and_resolve_mgr may strip out strictly worse hits from the data; they aren't required
bool write_json_hits_processor::do_requires_strictly_worse_hits() const {
	return false;
}

/// \brief Ctor for write_json_hits_processor
write_json_hits_processor::write_json_hits_processor(ref_vec<ostream> prm_ostreams ///< The ostream to which the results should be written
                                                     ) noexcept : super { move( prm_ostreams ) } {
	for (const ostream_ref &ostream_ref : get_ostreams() ) {
		json_writers.emplace_back( ostream_ref.get() );
	}
}


/// \brief Copy ctor for write_json_hits_processor
write_json_hits_processor::write_json_hits_processor(const write_json_hits_processor &prm_rhs ///< The other write_json_hits_processor from which to copy construct
                                                     ) : super       { prm_rhs             },
                                                         json_writers{ get_ostreams()      },
                                                         has_started { prm_rhs.has_started } {
	if ( has_started && ! json_writers.is_complete() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Unable to copy construct from write_json_hits_processor that's in-process of writing"));
	}
}

/// \brief Move ctor for write_json_hits_processor
write_json_hits_processor::write_json_hits_processor(write_json_hits_processor &&prm_rhs ///< The other write_json_hits_processor from which to move construct
                                                     ) : super       { move( prm_rhs )     },
                                                         json_writers{ get_ostreams()      },
                                                         has_started { prm_rhs.has_started } {
	if ( has_started && ! json_writers.is_complete() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Unable to copy construct from write_json_hits_processor that's in-process of writing"));
	}
}

