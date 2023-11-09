/// \file
/// \brief The overlap_type class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_TYPE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_TYPE_HPP

#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/score/length_getter/length_getter_enum.hpp"

#include <tuple>

namespace cath::score {

	/// \brief TODOCUMENT
	enum class overlap_type : char {
		SHORTER_OVER_LONGER,      ///< TODOCUMENT
		NUM_ALIGNED_OVER_SHORTER, ///< TODOCUMENT
		NUM_ALIGNED_OVER_LONGER   ///< TODOCUMENT
	};

	/// \brief TODOCUMENT
	static constexpr auto all_overlap_types = common::make_array(
		std::make_tuple( overlap_type::SHORTER_OVER_LONGER,      detail::length_getter_enum::SHORTER,     detail::length_getter_enum::LONGER  ),
		std::make_tuple( overlap_type::NUM_ALIGNED_OVER_SHORTER, detail::length_getter_enum::NUM_ALIGNED, detail::length_getter_enum::SHORTER ),
		std::make_tuple( overlap_type::NUM_ALIGNED_OVER_LONGER,  detail::length_getter_enum::NUM_ALIGNED, detail::length_getter_enum::LONGER  )
	);

	static_assert( cath::common::constexpr_is_uniq( all_overlap_types ), "all_overlap_types should not contain any repeated entries" );

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_TYPE_HPP
