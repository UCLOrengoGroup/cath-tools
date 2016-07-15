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

/// \brief Sort & resolve the hits data and output the results to the specified ostream
///
/// This is the work that's performed asynchronously in a separate thread once the data has been completely loaded
void process_hit_list(const string   &arg_query_id, ///< The query_id of the data
                      const hit_vec &&arg_hits,     ///< The hits
                      const str_vec &&arg_labels,   ///< The corresponding labels
                      ostream        &arg_ostream   ///< The ostream to which the results should be output
                      ) {
	// Build a hit_list of the hits and labels
	const hit_list the_hit_list{ move( arg_hits ), move( arg_labels ) };

	// Resolve the hits
	const auto best_result = resolve_hits( the_hit_list );

	// Output the results to arg_ostream
	arg_ostream << to_string( best_result.get_arch(), the_hit_list.get_labels(), hit_output_format::JON, arg_query_id );
}

/// \brief Record that the data that's been loaded under the current query_id is now complete and ready for processing
///
/// This: shifts the data into the computation_ data members
void read_and_resolve_mgr::complete() {
	// Check that any previously triggered computation is complete
	if ( resolve_future.valid() ) {
		resolve_future.wait();
	}

	// Shift the data into the computation_... data members so that
	// the computations won't conflict with the loading of new data
	//
	// (actually, this is probably overkill because std::async() copies values
	//  but not 100% familiar and better to be safe and this can all be done with fast moves here)
	computation_query_id = query_id;
	computation_hits     = move( the_hits );
	computation_labels   = move( labels   );

	// Start asynchronous processing of the work to process these hits and output the results
	resolve_future = async(
		launch::async,
		process_hit_list,
		computation_query_id,
		move( computation_hits    ),
		move( computation_labels  ),
		ref ( output_stream.get() )
	);

	// Re initialise the data_members to prepare for more data to be loaded in
	query_id.clear();
	the_hits = hit_vec{};
	labels   = str_vec{};
	active   = false;
}

/// \brief Check that any previously triggered computation is complete
void read_and_resolve_mgr::final_wait() const {
	if ( resolve_future.valid() ) {
		resolve_future.wait();
	}
}
