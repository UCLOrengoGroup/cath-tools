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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_H

#include "common/clone/check_uptr_clone_against_this.hpp"
#include "resolve_hits/options/spec/crh_score_spec.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"

#include <functional>
#include <iosfwd>
#include <string>

namespace cath { namespace rslv { class crh_filter_spec; } }
namespace cath { namespace rslv { class crh_output_spec; } }
namespace cath { namespace rslv { class full_hit_list; } }


namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Provide and ABC interface for classes that will process hits in read_and_process_mgr
			class hits_processor {
			private:
				/// \brief (A reference_wrapper to) the output stream to which the results should be written
				std::reference_wrapper<std::ostream> output_stream;

				/// \brief The score spec to apply to incoming hits
				crh_score_spec the_score_spec;

				/// \brief The segment spec to apply to incoming hits
				crh_segment_spec the_segment_spec;

				/// \brief Pure virtual method with which each concrete hits_processor must define how to create a clone of itself
				virtual std::unique_ptr<hits_processor> do_clone() const = 0;

				/// \brief Pure virtual method with which each concrete hits_processor must define how it processes a new hit for a query
				virtual void do_process_hits_for_query(const std::string &,
				                                       const crh_filter_spec &,
				                                       full_hit_list &) = 0;

				/// \brief Pure virtual method with which each concrete hits_processor must define how it finishes work
				virtual void do_finish_work() = 0;

				/// \brief Pure virtual method with which each concrete hits_processor must define whether it expects
				///        read_and_resolve_mgr to still parse hits that fail the score filter and pass them to this processor
				virtual bool do_parse_hits_that_fail_score_filter() const = 0;

			protected:
				std::ostream & get_ostream();
				const crh_score_spec & get_score_spec() const;
				const crh_segment_spec & get_segment_spec() const;

			public:
				explicit hits_processor(std::ostream &,
				                        const crh_score_spec &,
				                        const crh_segment_spec &) noexcept;
				virtual ~hits_processor() noexcept = default;

				hits_processor(const hits_processor &) = default;
				hits_processor(hits_processor &&) noexcept = default;
				hits_processor & operator=(const hits_processor &) = default;
				hits_processor & operator=(hits_processor &&) noexcept = default;

				std::unique_ptr<hits_processor> clone() const;

				void process_hits_for_query(const std::string &,
				                            const crh_filter_spec &,
				                            full_hit_list &);
				void finish_work();
				bool parse_hits_that_fail_score_filter() const;

			};

			/// \brief Getter for the ostream to which results should be written
			inline std::ostream & hits_processor::get_ostream() {
				return output_stream;
			}

			/// \brief Getter for the score spec to apply to incoming hits
			inline const crh_score_spec & hits_processor::get_score_spec() const {
				return the_score_spec;
			}

			/// \brief Getter for the segment spec to apply to incoming hits
			inline const crh_segment_spec & hits_processor::get_segment_spec() const {
				return the_segment_spec;
			}


			/// \brief Ctor
			inline hits_processor::hits_processor(std::ostream           &arg_output_stream,   ///< The ostream to which results should be written
			                                      const crh_score_spec   &arg_crh_score_spec,  ///< The score spec to apply to incoming hits
			                                      const crh_segment_spec &arg_crh_segment_spec ///< The segment spec to apply to incoming hits
			                                      ) noexcept : output_stream    { arg_output_stream    },
			                                                   the_score_spec   { arg_crh_score_spec   },
			                                                   the_segment_spec { arg_crh_segment_spec } {
			}

			/// \brief Standard approach to achieving a virtual copy-ctor
			inline std::unique_ptr<hits_processor> hits_processor::clone() const {
				return common::check_uptr_clone_against_this( do_clone(), *this );
			}

			/// \brief NVI pass-through to the virtual do_process_hits_for_query() method
			inline void hits_processor::process_hits_for_query(const std::string     &arg_query_id,    ///< The query_protein_id string
			                                                   const crh_filter_spec &arg_filter_spec, ///< The filter spec to apply to hits
			                                                   full_hit_list         &arg_full_hits    ///< The full hits to be processed
			                                                   ) {
				return do_process_hits_for_query(
					arg_query_id,
					arg_filter_spec,
					arg_full_hits
				);
			}

			/// \brief NVI pass-through to the virtual do_finish_work() method
			inline void hits_processor::finish_work() {
				do_finish_work();
			}

			/// \brief NVI pass-through to the virtual do_parse_hits_that_fail_score_filter() method
			inline bool hits_processor::parse_hits_that_fail_score_filter() const {
				return do_parse_hits_that_fail_score_filter();
			}

			std::unique_ptr<hits_processor> make_hits_processor(std::ostream &,
			                                                    const crh_output_spec &,
			                                                    const crh_score_spec &,
			                                                    const crh_segment_spec &);

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
