/// \file
/// \brief The check_scan_on_final_alignment class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_CHECK_SCAN_ON_FINAL_ALIGNMENT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_CHECK_SCAN_ON_FINAL_ALIGNMENT_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"

#include <cstddef>
#include <iosfwd>
#include <utility>

// clang-format off
namespace cath { class protein; }
namespace cath::align { class alignment; }
namespace cath::scan { class quad_criteria; }
namespace cath::scan { class scan_stride; }
namespace cath::scan::detail { class alignment_scan_comparison; }
namespace cath::scan::detail { class multi_struc_res_rep_pair; }
namespace cath::scan::detail { class quad_and_rep_criteria_result; }
namespace cath::scan::detail { class roled_scan_stride; }
namespace cath::scan::detail { class single_struc_res_pair; }
namespace cath::scan::detail { enum class quad_criteria_result : unsigned int; }
// clang-format on

namespace cath::scan::detail {

	/// \brief TODOCUMENT
	///
	/// A quad can fail because:
	///  * it doesn't pass the criteria itself
	///  * it *would* pass the criteria itself but its rep pair doesn't
	///
	/// \todo At present this doesn't actually check that the keyer parts
	///       correctly deliver all the matches that they should but
	///       this could be added if there's any indication that the
	///       predicted scores aren't being achieved.
	class check_scan_on_final_alignment final {
	private:
		static quad_criteria_result criteria_result_of(const quad_criteria &,
		                                               const multi_struc_res_rep_pair &,
		                                               const multi_struc_res_rep_pair &);

		static quad_criteria_result criteria_result_of(const quad_criteria &,
		                                               const single_struc_res_pair &,
		                                               const single_struc_res_pair &);

		static quad_criteria_result rep_quad_criteria_result_of(const protein &,
		                                                        const protein &,
		                                                        const quad_criteria &,
		                                                        const scan_stride &,
		                                                        const index_type &,
		                                                        const index_type &,
		                                                        const index_type &,
		                                                        const index_type &);

		static quad_criteria_result quad_criteria_result_of(const protein &,
		                                                    const protein &,
		                                                    const quad_criteria &,
		                                                    const index_type &,
		                                                    const index_type &,
		                                                    const index_type &,
		                                                    const index_type &);

		static quad_and_rep_criteria_result quad_and_rep_criteria_result_of(const protein &,
		                                                                    const protein &,
		                                                                    const quad_criteria &,
		                                                                    const scan_stride &,
		                                                                    const index_type &,
		                                                                    const index_type &,
		                                                                    const index_type &,
		                                                                    const index_type &);

	public:
		check_scan_on_final_alignment() = delete;

		/// \brief TODOCUMENT
		static constexpr size_t NUM_EXCLUDED_ON_SIDES = 5;

		[[nodiscard]] alignment_scan_comparison do_check( const align::alignment &,
		                                                  const protein &,
		                                                  const protein &,
		                                                  const quad_criteria &,
		                                                  const scan_stride & ) const;

		[[nodiscard]] std::pair<str_vec, str_vec> get_rep_name_lists( const protein &, const roled_scan_stride & ) const;
	};

	std::pair<index_vec, index_vec> get_rep_index_lists(const roled_scan_stride &,
	                                                    const index_type &);

	void print_highlight_rep_pymol_commands(std::ostream &,
	                                        const std::string &,
	                                        const std::pair<str_vec, str_vec> &);

	void print_highlight_rep_pymol_commands(std::ostream &,
	                                        const check_scan_on_final_alignment &,
	                                        const protein &,
	                                        const roled_scan_stride &);

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY_CHECK_SCAN_ON_FINAL_ALIGNMENT_HPP
