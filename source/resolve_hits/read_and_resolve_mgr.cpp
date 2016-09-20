/// \file
/// \brief The read_and_resolve_mgr class definitions

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

#include "read_and_resolve_mgr.h"

#include "exception/invalid_argument_exception.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/hit_resolver.h"
#include "resolve_hits/scored_hit_arch.h"

#include <iostream>
#include <string>
#include <thread>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;

using std::cerr;
using std::cout;
using std::future;
using std::launch;
using std::ostream;
using std::string;

constexpr bool read_and_resolve_mgr::DEFAULT_INPUT_HITS_ARE_PRESORTED;

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void read_and_resolve_mgr::process_query_id(const string &arg_query_id, ///< The query_protein_id string
                                            hit_vec      &arg_hits,     ///< The hits to process
                                            str_vec      &arg_labels,   ///< The labels corresponding to the hits to be processed
                                            ostream      &arg_ostream   ///< The ostream to which the results should be written
                                            ) {
	// Build a hit_list of the hits and labels
	const hit_list the_hit_list{ move( arg_hits ), move( arg_labels ) };
	arg_hits   = hit_vec{};
	arg_labels = str_vec{};

	// Resolve the hits
	const auto best_result = resolve_hits( the_hit_list );

	// Output the results to arg_ostream
	arg_ostream << to_string(
		best_result.get_arch(),
		the_hit_list.get_labels(),
		hit_output_format::JON,
		arg_query_id
	);
}

/// \brief Trigger asynchronous processing of the data corresponding to the specified protein_query_id
void read_and_resolve_mgr::trigger_async_process_query_id(const string &arg_query_id ///< THe protein_query_id
                                                          ) {
	// Wait until any previously triggered asynchronous processing is complete
	wait_for_any_active_work();

	// Mark to_be_erased_query_id with this query ID so it's possible to detect if there's any
	// attempt to modify the data for this query ID while it's being processed
	to_be_erased_query_id = arg_query_id;

	// Start asynchronous processing of the work to process these hits and output the results
	resolve_future = async(
		launch::async,
		process_query_id,
		arg_query_id,
		ref( hit_list_by_query_id[ arg_query_id ].first  ),
		ref( hit_list_by_query_id[ arg_query_id ].second ),
		ref( output_stream.get()                         )
	);
}
