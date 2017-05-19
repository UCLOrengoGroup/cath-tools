/// \file
/// \brief The write_html_hits_processor class definitions

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

#include "write_html_hits_processor.hpp"

#include "common/clone/make_uptr_clone.hpp"
#include "resolve_hits/full_hit.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/html_output/resolve_hits_html_outputter.hpp"

#include <memory>

using namespace cath::common;
using namespace cath::rslv::detail;

using std::move;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<hits_processor> write_html_hits_processor::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void write_html_hits_processor::do_process_hits_for_query(const string           &arg_query_id,     ///< The query_protein_id string
                                                          const crh_filter_spec  &arg_filter_spec,  ///< The filter_spec to apply to the hits
                                                          const crh_score_spec   &arg_score_spec,   ///< The score spec to apply to the hits
                                                          const crh_segment_spec &arg_segment_spec, ///< The segment spec to apply to the hits
                                                          const calc_hit_list    &arg_calc_hits     ///< The hits to process
                                                          ) {
	// If the prefix hasn't already been printed, then do so and record
	if ( ! printed_prefix ) {
		if ( ! html_spec.get_restrict_html_within_body() ) {
			get_ostream() << resolve_hits_html_outputter::html_prefix();
		}
		printed_prefix = true;
	}

	// Output the HTML for this query and its hits
	get_ostream() << resolve_hits_html_outputter::output_html(
		arg_query_id,
		arg_calc_hits,
		arg_score_spec,
		arg_segment_spec,
		html_spec,
		false,
		arg_filter_spec,
		batch_counter
	);
	++batch_counter;
}

/// \brief Write the HTML suffix to finish the work (if it has been started)
void write_html_hits_processor::do_finish_work() {
	if ( printed_prefix ) {
		get_ostream() << resolve_hits_html_outputter::html_key();
		if ( ! html_spec.get_restrict_html_within_body() ) {
			get_ostream() << resolve_hits_html_outputter::html_suffix();
		}
	}
}

/// \brief Return true: read_and_resolve_mgr should still parse hits that fail the score filter and pass them to this processor
bool write_html_hits_processor::do_wants_hits_that_fail_score_filter() const {
	return true;
}

/// \brief Ctor for the write_html_hits_processor
write_html_hits_processor::write_html_hits_processor(ostream       &arg_ostream,  ///< The ostream to which the results should be written
                                                     crh_html_spec  arg_html_spec ///< The specification for how to render the HTML
                                                     ) noexcept : super    { arg_ostream,               },
                                                                  html_spec{ std::move( arg_html_spec ) } {
}
