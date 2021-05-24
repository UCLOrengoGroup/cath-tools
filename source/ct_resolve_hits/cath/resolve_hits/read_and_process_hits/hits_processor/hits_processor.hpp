/// \file
/// \brief The hits_processor class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_HPP

#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"
#include "cath/resolve_hits/options/spec/crh_score_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"

#include <functional>
#include <iosfwd>
#include <string>

// clang-format off
namespace cath::rslv { class calc_hit_list; }
namespace cath::rslv { class crh_filter_spec; }
namespace cath::rslv { class crh_html_spec; }
namespace cath::rslv { class crh_single_output_spec; }
namespace cath::rslv { class full_hit_list; }
namespace cath::rslv::detail { class hits_processor_list; }
// clang-format on

namespace cath::rslv::detail {

	/// \brief Provide and ABC interface for classes that will process hits in read_and_process_mgr
	class hits_processor {
	private:
		/// \brief (A reference_wrapper to) the output stream to which the results should be written
		ref_vec<std::ostream> output_streams;

		/// \brief Pure virtual method with which each concrete hits_processor must define how to create a clone of itself
		[[nodiscard]] virtual std::unique_ptr<hits_processor> do_clone() const = 0;

		/// \brief Pure virtual method with which each concrete hits_processor must define how it processes a new hit for a query
		virtual void do_process_hits_for_query(const std::string &,
		                                       const crh_filter_spec &,
		                                       const crh_score_spec &,
		                                       const crh_segment_spec &,
		                                       const calc_hit_list &) = 0;

		/// \brief Pure virtual method with which each concrete hits_processor must define how it finishes work
		virtual void do_finish_work() = 0;

		/// \brief Pure virtual method with which each concrete hits_processor must define whether it expects
		///        read_and_resolve_mgr to still parse hits that fail the score filter and pass them to this processor
		[[nodiscard]] virtual bool do_wants_hits_that_fail_score_filter() const = 0;

		/// \brief Pure virtual method with which each concrete hits_processor must define whether it expects
		///        to see results even if they're strictly worse than other hits in the results
		[[nodiscard]] virtual bool do_requires_strictly_worse_hits() const = 0;

	  protected:
		const ref_vec<std::ostream> & get_ostreams();

	public:
		hits_processor() = default;
		explicit hits_processor(std::ostream &) noexcept;
		explicit hits_processor(ref_vec<std::ostream>) noexcept;
		virtual ~hits_processor() noexcept = default;

		hits_processor(const hits_processor &) = default;
		hits_processor(hits_processor &&) noexcept = default;
		hits_processor & operator=(const hits_processor &) = default;
		hits_processor & operator=(hits_processor &&) noexcept = default;

		[[nodiscard]] std::unique_ptr<hits_processor> clone() const;

		void process_hits_for_query(const std::string &,
		                            const crh_filter_spec &,
		                            const crh_score_spec &,
		                            const crh_segment_spec &,
		                            full_hit_list &);
		void process_hits_for_query(const std::string &,
		                            const crh_filter_spec &,
		                            const crh_score_spec &,
		                            const crh_segment_spec &,
		                            const calc_hit_list &);
		void finish_work();
		[[nodiscard]] bool wants_hits_that_fail_score_filter() const;
		[[nodiscard]] bool requires_strictly_worse_hits() const;
	};

	/// \brief Getter for the ostreams to which results should be written
	inline const ref_vec<std::ostream> & hits_processor::get_ostreams() {
		return output_streams;
	}

	/// \brief Ctor
	inline hits_processor::hits_processor(std::ostream &prm_output_stream ///< The ostream to which results should be written
	                                      ) noexcept : output_streams { { prm_output_stream } } {
	}

	/// \brief Ctor
	inline hits_processor::hits_processor(ref_vec<std::ostream> prm_stream_refs ///< A vector of reference_wrappers to the ostreams to which results should be written
	                                      ) noexcept : output_streams{ std::move( prm_stream_refs ) } {

	}

	/// \brief Standard approach to achieving a virtual copy-ctor
	inline std::unique_ptr<hits_processor> hits_processor::clone() const {
		return common::check_uptr_clone_against_this( do_clone(), *this );
	}

	/// \brief NVI pass-through to the virtual do_process_hits_for_query() method
	inline void hits_processor::process_hits_for_query(const std::string      &prm_query_id,         ///< The query_protein_id string
	                                                   const crh_filter_spec  &prm_filter_spec,      ///< The filter spec to apply to hits
	                                                   const crh_score_spec   &prm_crh_score_spec,   ///< The score spec to apply to incoming hits
	                                                   const crh_segment_spec &prm_crh_segment_spec, ///< The segment spec to apply to incoming hits
	                                                   full_hit_list          &prm_full_hits         ///< The full hits to be processed
	                                                   ) {
		const calc_hit_list the_calc_hit_list{
			std::move( prm_full_hits ),
			prm_crh_score_spec,
			prm_crh_segment_spec,
			prm_filter_spec
		};
		prm_full_hits = full_hit_list{};
		return do_process_hits_for_query(
			prm_query_id,
			prm_filter_spec,
			prm_crh_score_spec,
			prm_crh_segment_spec,
			the_calc_hit_list
		);
	}

	/// \brief NVI pass-through to the virtual do_process_hits_for_query() method
	inline void hits_processor::process_hits_for_query(const std::string      &prm_query_id,         ///< The query_protein_id string
	                                                   const crh_filter_spec  &prm_filter_spec,      ///< The filter spec to apply to hits
	                                                   const crh_score_spec   &prm_crh_score_spec,   ///< The score spec to apply to incoming hits
	                                                   const crh_segment_spec &prm_crh_segment_spec, ///< The segment spec to apply to incoming hits
	                                                   const calc_hit_list    &prm_calc_hits         ///< The calc hits to be processed
	                                                   ) {
		return do_process_hits_for_query(
			prm_query_id,
			prm_filter_spec,
			prm_crh_score_spec,
			prm_crh_segment_spec,
			prm_calc_hits
		);
	}

	/// \brief NVI pass-through to the virtual do_finish_work() method
	inline void hits_processor::finish_work() {
		do_finish_work();
	}

	/// \brief NVI pass-through to the virtual do_wants_hits_that_fail_score_filter() method
	inline bool hits_processor::wants_hits_that_fail_score_filter() const {
		return do_wants_hits_that_fail_score_filter();
	}

	/// \brief NVI pass-through to the virtual do_requires_strictly_worse_hits() method
	inline bool hits_processor::requires_strictly_worse_hits() const {
		return do_requires_strictly_worse_hits();
	}

} // namespace cath::rslv::detail

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_HPP
