/// \file
/// \brief The write_results_hits_processor class header

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

#ifndef WRITE_RESULTS_HITS_PROCESSOR_H_INCLUDED
#define WRITE_RESULTS_HITS_PROCESSOR_H_INCLUDED

#include "resolve_hits/options/spec/hit_boundary_output.h"
#include "resolve_hits/read_and_process_hits/hits_processor/hits_processor.h"

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Hits processor that writes results output to the hits_processor's ostream
			class write_results_hits_processor final : public hits_processor {
			private:
				/// \brief Convenience type alias for the parent class
				using super = hits_processor;

				/// \brief Whether to trim the boundaries before outputting them
				hit_boundary_output boundary_output;

				virtual std::unique_ptr<hits_processor> do_clone() const override final;

				virtual void do_process_hits_for_query(const std::string &,
				                                       const crh_filter_spec &,
				                                       full_hit_list &) override final;

				virtual void do_finish_work() override final;

				virtual bool do_parse_hits_that_fail_score_filter() const override final;

			public:
				explicit write_results_hits_processor(std::ostream &,
				                                      const crh_score_spec &,
				                                      const crh_segment_spec &,
				                                      const hit_boundary_output &) noexcept;
				virtual ~write_results_hits_processor() noexcept = default;
			};

		}
	}
}

#endif
