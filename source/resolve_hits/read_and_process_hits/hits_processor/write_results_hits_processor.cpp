/// \file
/// \brief The write_results_hits_processor class definitions

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

#include "write_results_hits_processor.hpp"

#include "common/clone/make_uptr_clone.hpp"
#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/hit_resolver.hpp"
#include "resolve_hits/scored_hit_arch.hpp"

using namespace cath::common;
using namespace cath::rslv::detail;

using std::move;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<hits_processor> write_results_hits_processor::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void write_results_hits_processor::do_process_hits_for_query(const string          &arg_query_id,    ///< The query_protein_id string
                                                             const crh_filter_spec &arg_filter_spec, ///< The filter_spec to apply to the hits
                                                             full_hit_list         &arg_full_hits    ///< The hits to process
                                                             ) {
	// Build a calc_hit_list of the hits and labels
	const calc_hit_list the_calc_hit_list{ move( arg_full_hits ), get_score_spec(), get_segment_spec(), arg_filter_spec };
	arg_full_hits = full_hit_list{};

	// Resolve the hits
	const auto result_hit_arch  = resolve_hits( the_calc_hit_list );
	const auto result_full_hits = get_full_hits_of_hit_arch(
		result_hit_arch,
		the_calc_hit_list.get_full_hits()
	);

	// Output the results to arg_ostream
	get_ostream() << to_output_string(
		result_full_hits,
		get_segment_spec(),
		hit_output_format::JON,
		arg_query_id,
		boundary_output
	);
}

/// \brief Do nothing to finish the batch of work
void write_results_hits_processor::do_finish_work() {
}

/// \brief Return true: read_and_resolve_mgr needn't parse hits that fail the score filter or pass them to this processor
bool write_results_hits_processor::do_parse_hits_that_fail_score_filter() const {
	return false;
}

/// \brief Ctor for write_results_hits_processor
write_results_hits_processor::write_results_hits_processor(ostream                   &arg_ostream,        ///< The ostream to which the results should be written
                                                           const crh_score_spec      &arg_score_spec,     ///< The score_spec to apply to hits
                                                           const crh_segment_spec    &arg_segment_spec,   ///< The segment_spec to apply to hits
                                                           const hit_boundary_output &arg_boundary_output ///< Whether to trim the boundaries before outputting them
                                                           ) noexcept : super{ arg_ostream, arg_score_spec, arg_segment_spec },
                                                                        boundary_output{ arg_boundary_output } {
}

