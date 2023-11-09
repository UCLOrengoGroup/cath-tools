/// \file
/// \brief The alignment gap header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_ALIGNMENT_GAP_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_ALIGNMENT_GAP_HPP

#include <cstddef>

#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::align::gap { class gap_penalty; }
// clang-format on

namespace cath::align::gap {

	namespace detail {

		size_t get_naive_num_gaps_of_entry(const alignment &,
											const size_t &);

		size_size_pair gap_open_and_extend_counts_of_pair_in_alignment(const alignment &,
																		const size_t &,
																		const size_t &);

	} // namespace detail

	float_score_type gap_count_of_alignment(const alignment &);

	float_score_float_score_pair gap_open_and_extend_counts_of_alignment(const alignment &);

	float_score_type gap_penalty_value_of_alignment(const alignment &,
													const gap_penalty &);

	size_t get_naive_num_gaps(const alignment &);

} // namespace cath::align::gap

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_ALIGNMENT_GAP_HPP
