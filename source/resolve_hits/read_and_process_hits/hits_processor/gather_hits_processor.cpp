/// \file
/// \brief The gather_hits_processor class definitions

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

#include "gather_hits_processor.hpp"

#include "common/clone/make_uptr_clone.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<hits_processor> gather_hits_processor::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Process the specified data
///
/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
void gather_hits_processor::do_process_hits_for_query(const string           &prm_query_id,         ///< The query_protein_id string
                                                      const crh_filter_spec  &/*prm_filter_spec*/,  ///< The filter_spec to apply to the hits
                                                      const crh_score_spec   &/*prm_score_spec*/,   ///< The score spec to apply to the hits
                                                      const crh_segment_spec &/*prm_segment_spec*/, ///< The segment spec to apply to the hits
                                                      const calc_hit_list    &prm_calc_hits         ///< The hits to process
                                                      ) {
	hit_lists.get().emplace_back(
		prm_query_id,
		prm_calc_hits
	);
}

/// \brief Write the HTML suffix to finish the work (if it has been started)
void gather_hits_processor::do_finish_work() {
}

/// \brief Return true: should store all hits
bool gather_hits_processor::do_wants_hits_that_fail_score_filter() const {
	return true;
}

/// \brief Return true: should store all hits
bool gather_hits_processor::do_requires_strictly_worse_hits() const {
	return true;
}

/// \brief Ctor from the data structure into which the data should be placed
gather_hits_processor::gather_hits_processor(str_calc_hit_list_pair_vec &prm_hit_lists ///< The data structure into which the data should be placed
                                             ) noexcept : hit_lists { prm_hit_lists } {
}
