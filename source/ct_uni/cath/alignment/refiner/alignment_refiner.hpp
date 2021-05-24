/// \file
/// \brief The alignment_refiner class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_ALIGNMENT_REFINER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_ALIGNMENT_REFINER_HPP

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath { class protein_list; }
namespace cath::align { class alignment; }
namespace cath::align::detail { class alignment_split; }
namespace cath::align::detail { class alignment_split_list; }
namespace cath::align::gap { class gap_penalty; }
namespace cath::index { class view_cache_list; }
// clang-format on

namespace cath::align {

	/// \brief TODOCUMENT
	class alignment_refiner final {
	private:
		/// \brief TODOCUMENT
		float_score_vec_vec from_alignment_scores;

		/// \brief TODOCUMENT
		float_score_vec_vec to_alignment_scores;

		detail::bool_aln_pair iterate_step(const alignment &,
		                                   const protein_list &,
		                                   const index::view_cache_list &,
		                                   const gap::gap_penalty &);

		detail::bool_aln_pair iterate_step_for_alignment_split_list(const alignment &,
		                                                            const protein_list &,
		                                                            const index::view_cache_list &,
		                                                            const gap::gap_penalty &,
		                                                            const detail::alignment_split_list &);

		detail::bool_aln_pair iterate_step_for_alignment_split(const alignment &,
		                                                       const protein_list &,
		                                                       const index::view_cache_list &,
		                                                       const gap::gap_penalty &,
		                                                       const detail::alignment_split &);

	public:
		alignment iterate(const alignment &,
		                  const protein_list &,
		                  const gap::gap_penalty &);

		alignment iterate(const alignment &,
		                  const protein_list &,
		                  const index::view_cache_list &,
		                  const gap::gap_penalty &);

		alignment iterate_join(const alignment &,
		                       const protein_list &,
		                       const gap::gap_penalty &,
		                       const size_vec &);

		alignment iterate_join(const alignment &,
		                       const protein_list &,
		                       const index::view_cache_list &,
		                       const gap::gap_penalty &,
		                       const size_vec &);

	};

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_ALIGNMENT_REFINER_HPP
