/// \file
/// \brief The quad_and_rep_criteria_result class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_QUAD_AND_REP_CRITERIA_RESULT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_QUAD_AND_REP_CRITERIA_RESULT_HPP

#include "cath/scan/detail/check_scan/test_only/quad_criteria_result.hpp"

namespace cath::scan::detail {

	/// \brief TODOCUMENT
	class quad_and_rep_criteria_result final {
	private:
		/// \brief TODOCUMENT
		quad_criteria_result rep_status;

		/// \brief TODOCUMENT
		quad_criteria_result quad_status;

		void sanity_check() const;

	public:
		quad_and_rep_criteria_result(const quad_criteria_result &,
		                             const quad_criteria_result &);
		[[nodiscard]] const quad_criteria_result &get_rep_status() const;
		[[nodiscard]] const quad_criteria_result &get_quad_status() const;
	};

	std::string to_string(const quad_and_rep_criteria_result &);
	std::ostream & operator<<(std::ostream &,
	                          const quad_and_rep_criteria_result &);

	bool operator<(const quad_and_rep_criteria_result &,
	               const quad_and_rep_criteria_result &);

	bool fully_passes(const quad_and_rep_criteria_result &);

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_QUAD_AND_REP_CRITERIA_RESULT_HPP
