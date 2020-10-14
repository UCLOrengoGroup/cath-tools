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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_TYPE_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_TYPE_HPP

#include "common/algorithm/constexpr_is_uniq.hpp"
#include "common/cpp20/make_array.hpp"
#include "score/length_getter/length_getter_enum.hpp"

#include <tuple>

namespace cath {
	namespace score {

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

	} // namespace score
} // namespace cath

#endif
