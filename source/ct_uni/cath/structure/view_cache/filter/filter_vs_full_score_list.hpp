/// \file
/// \brief The filter_vs_full_score_list class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_LIST_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_LIST_HPP

#include <filesystem>

#include "cath/score/score_type_aliases.hpp"
#include "cath/structure/structure_type_aliases.hpp"
#include "cath/structure/view_cache/filter/filter_vs_full_score.hpp"

#include <cstddef>

namespace cath { namespace index { namespace filter { class filter_vs_full_score; } } }
namespace cath { namespace score { class classn_stat; } }
namespace cath { namespace score { class true_false_pos_neg;   } }

namespace cath {
	namespace index {
		namespace filter {

			/// \brief A list of filter_vs_full_score objects
			///
			/// At present, this isn't far from a `const filter_vs_full_score_vec` that keeps all entries
			/// sorted by full_score.
			///
			/// Given a score that's slow to calculate for entries (eg SSAP score for previously unaligned
			/// pairs of structures), it can be helpful to use a fast filter to generate a score that
			/// identifies entries likely to get good scores. That are worth doing the slow calculations on.
			///
			/// Invariants:
			///  * the vector of filter_vs_full_scores is kept sorted on full_score (ascending)
			class filter_vs_full_score_list final {
			private:
				/// \brief The filter_vs_full_score objects, stored in a vector
				filter_vs_full_score_vec filter_vs_full_scores;

				void sort_filter_vs_full_scores();

			public:
				filter_vs_full_score_list();
				explicit filter_vs_full_score_list(filter_vs_full_score_vec);

				void add_filter_vs_full_score(const filter_vs_full_score &);
				size_t size() const;
				const filter_vs_full_score & operator[](const size_t &) const;

				/// \brief A const_iterator type alias as part of making filter_vs_full_score_list into a (const) range
				using const_iterator = filter_vs_full_score_vec::const_iterator;

				/// \brief A iterator type alias. Though filter_vs_full_score_list doesn't currently offer a mutable cbegin()/end(),
				///        this type alias is required by Boost's sub_range<>
				using iterator = filter_vs_full_score_vec::const_iterator;

				const_iterator begin() const;
				const_iterator end() const;
			};

			/// \brief Type alias for filter_vs_full_score_list's const_iterator
			using filter_vs_full_score_list_citr = filter_vs_full_score_list::const_iterator;

			double filter_score_full_score_with_sensitivity(const filter_vs_full_score_list &,
			                                                const double &,
			                                                const double &);

			filter_vs_full_score filter_attempt_full_score_with_sensitivity(const filter_vs_full_score_list &,
			                                                                const double &,
			                                                                const double &);

			score::true_false_pos_neg filter_result_full_score_with_sensitivity(const filter_vs_full_score_list &,
			                                                                    const double &,
			                                                                    const double &);

			filter_vs_full_score_vec filter_attempts_with_sensitivity(const filter_vs_full_score_list &,
			                                                          const double &);

			score::doub_true_false_pos_neg_pair_vec filter_results_with_sensitivity(const filter_vs_full_score_list &,
			                                                                        const double &);

			// double filter_score_with_sensitivity(const filter_vs_full_score_list &,
			//                                      const double &);
	
			// filter_vs_full_score filter_attempt_with_sensitivity(const filter_vs_full_score_list &,
			//                                                      const double &);

			score::true_false_pos_neg assess_results_on_filter_attempt(const filter_vs_full_score_list &,
			                                                           const filter_vs_full_score &);

			void gnuplot_data(const filter_vs_full_score_list &,
			                  const ::std::filesystem::path &,
			                  const filter_vs_full_score_list & = filter_vs_full_score_list());

			void gnuplot_classsn_stat_for_recall(const filter_vs_full_score_list &,
			                                     const ::std::filesystem::path &,
			                                     const score::classn_stat &,
			                                     const double &);

			void gnuplot_classsn_stat_for_recall(const score::doub_true_false_pos_neg_pair_vec &,
			                                     const ::std::filesystem::path &,
			                                     const score::classn_stat &);
		} // namespace filter
	} // namespace index
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_LIST_HPP
