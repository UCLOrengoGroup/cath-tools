/// \file
/// \brief The multi_align_group class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DETAIL_MULTI_ALIGN_GROUP_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DETAIL_MULTI_ALIGN_GROUP_HPP

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/refiner/alignment_refiner.hpp"

// clang-format off
namespace cath::file { class pdb_list; }
// clang-format on

namespace cath::align::detail {

	/// \brief Represent a group of entries and the alignment between them
	///
	/// This is an implementation class for build_alignment_from_parts()
	/// (via multi_align_builder)
	class multi_align_group final {
	private:
		/// \brief The alignment between the entries in the group
		///
		/// This is in the same order as the entries
		alignment the_alignment;

		/// \brief The list of entries represented in this group
		///
		/// These entries are in the same order as the entries in the alignment
		size_vec  entries;

	public:
		void refine_join(alignment_refiner &,
		                 const protein_list &,
		                 const gap::gap_penalty &,
		                 const size_vec &);

		explicit multi_align_group(const size_t &);
		multi_align_group(alignment,
		                  const size_t &,
		                  const size_t &);

		[[nodiscard]] const alignment &get_alignment() const;

		[[nodiscard]] const size_vec &get_entries() const;

		void add_alignment(const size_t &,
		                   const alignment &);

		void glue_in_copy_of_group(const multi_align_group &,
		                           const size_t &);

		void refine_alignment(alignment_refiner &,
		                      const protein_list &,
		                      const gap::gap_penalty &);
	};

	std::ostream & operator<<(std::ostream &,
	                          const multi_align_group &);

	size_t get_index_of_entry(const multi_align_group &,
	                          const size_t &);

	void glue_in_alignment(multi_align_group &,
	                       const alignment &,
	                       const size_t &,
	                       const size_t &);

} // namespace cath::align::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DETAIL_MULTI_ALIGN_GROUP_HPP
