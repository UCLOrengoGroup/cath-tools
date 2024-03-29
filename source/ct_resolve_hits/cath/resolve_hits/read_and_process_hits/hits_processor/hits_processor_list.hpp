/// \file
/// \brief The hits_processor_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_LIST_HPP

#include <boost/range/adaptor/indirected.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/clone/clone_ptr.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor.hpp"

#include <initializer_list>
#include <utility>

// clang-format off
namespace cath::common { class ofstream_list; }
namespace cath::rslv { class crh_output_spec; }
// clang-format on

namespace cath::rslv::detail {

	/// \brief A list of hits_processors to process hits
	///
	/// Be careful: this has an initializer_list ctor
	class hits_processor_list final {
	private:
		/// \brief A vector of (clone_ptrs to) hits_processors
		hits_processor_clptr_vec processors;

		/// \brief The score spec to apply to incoming hits
		crh_score_spec the_score_spec;

		/// \brief The segment spec to apply to incoming hits
		crh_segment_spec the_segment_spec;

		[[nodiscard]] const crh_score_spec &  get_score_spec() const;
		[[nodiscard]] const crh_segment_spec &get_segment_spec() const;

	  public:
		/// \brief A const_iterator type alias as part of making this a range
		///
		/// Note that this pipes through boost::indirected_range so the range is
		/// over elements of type `const hits_processor &` rather than `clone_ptr<const hits_processor>`
		using const_iterator = common::range_const_iterator_t< boost::indirected_range<const hits_processor_clptr_vec> >;

		hits_processor_list() = default;

		explicit hits_processor_list(const crh_score_spec &,
		                             const crh_segment_spec & = crh_segment_spec{},
		                             std::initializer_list<hits_processor_clptr> = {});

		hits_processor_list(const hits_processor_list &) = default;
		hits_processor_list(hits_processor_list &&) = default;
		hits_processor_list & operator=(const hits_processor_list &) = default;
		hits_processor_list & operator=(hits_processor_list &&) = default;

		hits_processor_list & add_processor(const hits_processor &);
		hits_processor_list & add_processor(hits_processor_uptr);
		hits_processor_list & add_processor(hits_processor_clptr);

		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t size() const;

		const hits_processor & operator[](const size_t &) const;

		[[nodiscard]] bool wants_hits_that_fail_score_filter() const;

		[[nodiscard]] bool requires_strictly_worse_hits() const;

		void process_hits_for_query(const std::string &,
		                            const crh_filter_spec &,
		                            full_hit_list);
		void finish_work();

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	hits_processor_list make_hits_processors(common::ofstream_list &,
	                                         const crh_single_output_spec &,
	                                         const crh_output_spec &,
	                                         const crh_score_spec &,
	                                         const crh_segment_spec &,
	                                         const crh_html_spec &);


	/// \brief Process the specified full_hit_list for the specified query using the specified crh_filter_spec
	///
	/// This builds a calc_hit_list from the specified full_hit_list once and then passes it to each of the hits_processors
	inline void hits_processor_list::process_hits_for_query(const std::string     &prm_query_id,    ///< The query_protein_id string
	                                                        const crh_filter_spec &prm_filter_spec, ///< The filter spec to apply to hits
	                                                        full_hit_list          prm_full_hits    ///< The full hits to be processed
	                                                        ) {
		const calc_hit_list the_calc_hit_list{
			std::move( prm_full_hits ),
			get_score_spec(),
			get_segment_spec(),
			prm_filter_spec,
			(
				// If strictly worse hits were required in the full_hits, then it's worth pruning them out now,
				// otherwise don't (because that'd just be wasted effort)
				requires_strictly_worse_hits()
					? seg_dupl_hit_policy::PRUNE
					: seg_dupl_hit_policy::PRESERVE
			)
		};
		prm_full_hits = full_hit_list{};
		boost::for_each(
			processors,
			[&] (common::clone_ptr<hits_processor> &x) {
				x->process_hits_for_query(
					prm_query_id,
					prm_filter_spec,
					get_score_spec(),
					get_segment_spec(),
					the_calc_hit_list
				);
			}
		);
	}

	/// \brief Get each of the hits_processors in the list to finish any work they've started
	inline void hits_processor_list::finish_work() {
		boost::for_each(
			processors,
			[] (common::clone_ptr<hits_processor> &x) {
				x->finish_work();
			}
		);
	}

} // namespace cath::rslv::detail

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR_HITS_PROCESSOR_LIST_HPP
