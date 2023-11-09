/// \file
/// \brief The path_step header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_PATH_STEP_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_PATH_STEP_HPP

#include <array>

#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
// clang-format on

namespace cath::align::detail {

	enum class path_step : char {
		ALIGN_PAIR,
		INSERT_INTO_FIRST,
		INSERT_INTO_SECOND
	};

	/// \brief TODOCUMENT
	using path_step_vec     = std::vector<path_step>;

	/// \brief TODOCUMENT
	using path_step_vec_vec = std::vector<path_step_vec>;

	/// \brief TODOCUMENT
	inline constexpr ::std::array ALL_PATH_STEPS = { path_step::ALIGN_PAIR,
		                                             path_step::INSERT_INTO_FIRST,
		                                             path_step::INSERT_INTO_SECOND };

	std::ostream & operator<<(std::ostream &,
	                          const path_step &);

	size_size_pair indices_of_path_step(const path_step &);
	size_size_pair indices_of_point_after_path_step(const path_step &,
	                                                const size_t &,
	                                                const size_t &);
	std::map<path_step, size_size_pair> indices_of_point_by_path_step(const size_t &,
				                                                      const size_t &);

	void append_path_step_to_pair_alignment_from_point(alignment &,
	                                                   const path_step &,
	                                                   const size_size_pair &);

} // namespace cath::align::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_PATH_STEP_HPP
