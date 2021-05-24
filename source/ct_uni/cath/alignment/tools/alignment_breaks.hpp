/// \file
/// \brief The alignment_breaks header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TOOLS_ALIGNMENT_BREAKS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TOOLS_ALIGNMENT_BREAKS_HPP

#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
// clang-format on

namespace cath::align {
	namespace detail {

		/// \brief TODOCUMENT
		enum class break_pair_validity : bool {
			GOOD,
			BAD
		};

		/// \brief TODOCUMENT
		enum class break_pair_future : bool {
			NEVER_AGAIN,
			MAYBE_LATER
		};

		/// \brief TODOCUMENT
		using break_pair_validity_and_future = std::pair<break_pair_validity, break_pair_future>;

		break_pair_validity_and_future check_pair(const alignment &,
		                                          const size_t &,
		                                          const size_t &);
	} // namespace detail

	size_vec get_alignment_breaks(const alignment &);
	size_size_pair_vec get_alignment_break_pairs(const alignment &);

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TOOLS_ALIGNMENT_BREAKS_HPP
