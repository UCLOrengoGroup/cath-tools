/// \file
/// \brief The write_json_hits_processor class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_WRITE_JSON_HITS_PROCESSOR_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_WRITE_JSON_HITS_PROCESSOR_HPP

#include <rapidjson/ostreamwrapper.h>

#include "cath/common/rapidjson_addenda/rapidjson_writer_list.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor.hpp"

namespace cath::rslv::detail {

	/// \brief Hits processor that writes results output to the hits_processor's ostream
	class write_json_hits_processor final : public hits_processor {
	private:
		/// \brief Convenience type alias for the parent class
		using super = hits_processor;

		/// \brief The JSON writer, which writes to an OStreamWrapper of the hits_processor's ostream
		common::rapidjson_writer_list<common::json_style::PRETTY, rapidjson::OStreamWrapper> json_writers;

		/// \brief Whether anything has been written to this yet
		bool has_started = false;

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
		explicit write_json_hits_processor(ref_vec<std::ostream>) noexcept;

		write_json_hits_processor(const write_json_hits_processor &);
		write_json_hits_processor(write_json_hits_processor &&);
		write_json_hits_processor & operator=(const write_json_hits_processor &) = delete;
		write_json_hits_processor & operator=(write_json_hits_processor &&) = delete;
	};

} // namespace cath::rslv::detail

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_WRITE_JSON_HITS_PROCESSOR_HPP
