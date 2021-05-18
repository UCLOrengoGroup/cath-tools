/// \file
/// \brief The ssaps_and_prcs_of_query class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAPS_AND_PRCS_OF_QUERY_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAPS_AND_PRCS_OF_QUERY_HPP

#include <iosfwd>
#include <vector>

#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/score/homcheck_tools/ssap_and_prc.hpp"
#include "cath/score/score_type_aliases.hpp"

namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class ssap_scores_entry; } }
namespace cath { namespace file { using prc_scores_entry_vec = std::vector<prc_scores_entry>; } }
namespace cath { namespace file { using ssap_scores_entry_vec = std::vector<ssap_scores_entry>; } }
namespace cath { namespace homcheck { class ssap_and_prc; } }
namespace cath { namespace homcheck { class superfamily_of_domain; } }
namespace cath { namespace score { class rbf_model; } }


namespace cath {
	namespace homcheck {

		/// \brief Represent the ssap_and_prc results associated with a single query
		///
		/// \invariant All ssap_and_prc_entries have the same query_id
		class ssaps_and_prcs_of_query final {
		private:
			/// \brief The SSAP and PRC entries
			ssap_and_prc_vec ssap_and_prc_entries;

			void sanity_check() const;

		public:
			using const_iterator = ssap_and_prc_vec::const_iterator;

			ssaps_and_prcs_of_query() = default;
			explicit ssaps_and_prcs_of_query(ssap_and_prc_vec);

			void calculate_all_svm_scores(const score::rbf_model &);

			[[nodiscard]] bool   empty() const;
			[[nodiscard]] size_t size() const;

			const ssap_and_prc & operator[](const size_t &) const;

			[[nodiscard]] const_iterator begin() const;
			[[nodiscard]] const_iterator end() const;
		};

		ssaps_and_prcs_of_query calculate_all_svm_scores_copy(ssaps_and_prcs_of_query,
		                                                      const score::rbf_model &);

		const std::string & get_query_id(const ssaps_and_prcs_of_query &);

		ssap_and_prc_cref_opt best_svm_assignable(const ssaps_and_prcs_of_query &,
		                                          const superfamily_of_domain &,
		                                          const double & = -0.1);

		ssap_and_prc_cref_opt best_magic_function_assignable(const ssaps_and_prcs_of_query &,
		                                                     const superfamily_of_domain &,
		                                                     const double & = -0.1);

		file::ssap_scores_entry_cref_opt best_fold_level_match(const file::ssap_scores_entry_vec &,
		                                                       const superfamily_of_domain &);

		ssaps_and_prcs_of_query make_ssaps_and_prcs_of_query(const file::ssap_scores_entry_vec &,
		                                                     const file::prc_scores_entry_vec &);

	} // namespace homcheck
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SSAPS_AND_PRCS_OF_QUERY_HPP
