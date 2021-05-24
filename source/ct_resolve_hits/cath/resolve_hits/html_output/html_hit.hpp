/// \file
/// \brief The html_hit class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_HIT_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_HIT_HPP

#include "cath/resolve_hits/full_hit.hpp"

#include "cath/display_colour/display_colour.hpp"

#include <functional>

// clang-format off
namespace cath::rslv { class full_hit; }
// clang-format on

namespace cath::rslv {

	/// \brief Represent a hit to be rendered in HTML
	struct html_hit final {

		/// \brief A const-reference to the full_hit to be rendered
		std::reference_wrapper<const full_hit> hit_ref;

		/// \brief The index of the batch of data from which this hit came
		///        (where a batch is the bunch of hits relating to one query ID;
		///         the HTML can display multiple batches)
		///
		/// This is used to put a generate a unique ID for the hit in the HTML
		size_t batch_idx;

		/// \brief The index of this hit within its batch of data
		///
		/// This is used to put a generate a unique ID for the hit in the HTML
		size_t hit_idx;

		/// \brief The colour in which the hit should be rendered
		display_colour colour;

		/// \brief The resolved boundaries for the hit if this is being rendered as a result
		///        or nullopt otherwise
		seg_boundary_pair_vec_opt result_boundaries;

	};

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_HIT_HPP
