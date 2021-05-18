/// \file
/// \brief The alignment_scan_comparison class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_ALIGNMENT_SCAN_COMPARISON_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_ALIGNMENT_SCAN_COMPARISON_HPP

#include "cath/scan/detail/check_scan/test_only/quad_and_rep_criteria_result.hpp"

#include <iosfwd>
#include <map>
#include <utility>
#include <vector>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			class alignment_scan_comparison final {
			public:
				/// \brief TODOCUMENT
				using crit_res_doub_map = std::map<quad_and_rep_criteria_result, double>;

				/// \brief TODOCUMENT
				using const_iterator = crit_res_doub_map::const_iterator;

			private:
				/// \brief TODOCUMENT
				crit_res_doub_map score_by_result;

			public:
				alignment_scan_comparison & operator+=(const std::pair<quad_and_rep_criteria_result, double> &);

				[[nodiscard]] bool   has_score_of_criteria_result( const quad_and_rep_criteria_result & ) const;
				[[nodiscard]] double get_score_of_criteria_result( const quad_and_rep_criteria_result & ) const;

				[[nodiscard]] const_iterator begin() const;
				[[nodiscard]] const_iterator end() const;
			};

			std::vector<quad_and_rep_criteria_result> get_criteria_results(const alignment_scan_comparison &);
			std::vector<quad_and_rep_criteria_result> get_criteria_results_sorted_by_score_desc(const alignment_scan_comparison &);

			std::ostream & operator<<(std::ostream &,
			                          const alignment_scan_comparison &);

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_ALIGNMENT_SCAN_COMPARISON_HPP
