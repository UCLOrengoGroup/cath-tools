/// \file
/// \brief The summarise_hits_processor class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_SUMMARISE_HITS_PROCESSOR_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_SUMMARISE_HITS_PROCESSOR_H

#include "resolve_hits/full_hit.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/hits_processor.hpp"

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief A hits_processor to summarise the input data
			class summarise_hits_processor final : public hits_processor {
			private:
				/// \brief Convenience type alias for the parent class
				using super = hits_processor;

				/// \brief For each query, record the maximum stop of the hits
				size_vec max_stops;

				/// \brief Record the number of hits
				size_t num_hits = 0;

				/// \brief Record an example query_id/full_hit pair
				str_full_hit_pair_opt example_query_id_and_hit;

				std::unique_ptr<hits_processor> do_clone() const final;

				void do_process_hits_for_query(const std::string &,
				                               const crh_filter_spec &,
				                               full_hit_list &) final;

				void do_finish_work() final;

				bool do_parse_hits_that_fail_score_filter() const final;

			public:
				explicit summarise_hits_processor(std::ostream &,
				                                  const crh_score_spec &,
				                                  const crh_segment_spec &) noexcept;
			};


		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
