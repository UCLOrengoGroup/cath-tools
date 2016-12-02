/// \file
/// \brief The hits_processor class definitions

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

#include "hits_processor.hpp"

#include "resolve_hits/options/spec/crh_output_spec.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/write_html_hits_processor.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/write_results_hits_processor.hpp"

using namespace cath::rslv::detail;

using std::make_unique;
using std::ostream;
using std::unique_ptr;

/// \brief Make the hits_processor implied by the specified spec object and the specified ostream
unique_ptr<hits_processor> cath::rslv::detail::make_hits_processor(ostream                &arg_ostream,     ///< The ostream to which the new hits_processor should write
                                                                   const crh_output_spec  &arg_output_spec, ///< The crh_output_spec defining the type of hits_processor to make
                                                                   const crh_score_spec   &arg_score_spec,  ///< The crh_score_spec defining the type of hits_processor to make
                                                                   const crh_segment_spec &arg_segment_spec ///< The crh_segment_spec defining the type of hits_processor to make
                                                                   ) {
	if ( arg_output_spec.get_generate_html_output() || arg_output_spec.get_restrict_html_within_body() ) {
		return make_unique<write_html_hits_processor>(
			arg_ostream,
			arg_score_spec,
			arg_segment_spec,
			arg_output_spec.get_restrict_html_within_body()
		);
	}
	else {
		return make_unique<write_results_hits_processor>(
			arg_ostream,
			arg_score_spec,
			arg_segment_spec,
			arg_output_spec.get_boundary_output()
		);
	}
}