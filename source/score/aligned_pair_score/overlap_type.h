/// \file
/// \brief The overlap_type class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef OVERLAP_TYPE_H_INCLUDED
#define OVERLAP_TYPE_H_INCLUDED

#include "common/algorithm/constexpr_is_uniq.h"
#include "score/length_getter/length_getter_enum.h"

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		enum class overlap_type {
			SHORTER_OVER_LONGER,      ///< TODOCUMENT
			NUM_ALIGNED_OVER_SHORTER, ///< TODOCUMENT
			NUM_ALIGNED_OVER_LONGER   ///< TODOCUMENT
		};

		static constexpr std::array<std::tuple<overlap_type, detail::length_getter_enum, detail::length_getter_enum>, 3> all_overlap_types = { {
			std::make_tuple( overlap_type::SHORTER_OVER_LONGER,      detail::length_getter_enum::SHORTER,     detail::length_getter_enum::LONGER  ),
			std::make_tuple( overlap_type::NUM_ALIGNED_OVER_SHORTER, detail::length_getter_enum::NUM_ALIGNED, detail::length_getter_enum::SHORTER ),
			std::make_tuple( overlap_type::NUM_ALIGNED_OVER_LONGER,  detail::length_getter_enum::NUM_ALIGNED, detail::length_getter_enum::LONGER  )
		} };

		static_assert( cath::common::constexpr_is_uniq( all_overlap_types ), "all_overlap_types should not contain any repeated entries" );

	}
}

#endif
