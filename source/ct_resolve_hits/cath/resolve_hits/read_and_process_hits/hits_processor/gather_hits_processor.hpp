/// \file
/// \brief The gather_hits_processor class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_GATHER_HITS_PROCESSOR_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_GATHER_HITS_PROCESSOR_HPP

#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief A hits_processor to gather the input data
			class gather_hits_processor final : public hits_processor {
			private:
				/// \brief Convenience type alias for the parent class
				using super = hits_processor;

				/// \brief A reference to the data structure into which the query IDs and
				///        associated calc_hit_lists should be placed
				std::reference_wrapper<str_calc_hit_list_pair_vec> hit_lists;

				[[nodiscard]] std::unique_ptr<hits_processor> do_clone() const final;

				void do_process_hits_for_query(const std::string &,
				                               const crh_filter_spec &,
				                               const crh_score_spec &,
				                               const crh_segment_spec &,
				                               const calc_hit_list &) final;

				void do_finish_work() final;

				[[nodiscard]] bool do_wants_hits_that_fail_score_filter() const final;

				[[nodiscard]] bool do_requires_strictly_worse_hits() const final;

			  public:
				explicit gather_hits_processor(str_calc_hit_list_pair_vec &) noexcept;
			};


		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_GATHER_HITS_PROCESSOR_HPP
