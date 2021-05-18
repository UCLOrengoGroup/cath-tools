/// \file
/// \brief The alignment_split_list class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_LIST_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_LIST_HPP

#include "cath/alignment/refiner/detail/alignment_split.hpp"

#include <set>

namespace cath { namespace align { class alignment; } }

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			class alignment_split_list final {
			private:
				using alignment_split_set = std::set<alignment_split>;

				alignment_split_set splits;

			public:
				void insert(const alignment_split &);

				using iterator = alignment_split_set::iterator;
				using const_iterator = alignment_split_set::const_iterator;
				iterator begin();
				iterator end();
				[[nodiscard]] const_iterator begin() const;
				[[nodiscard]] const_iterator end() const;
			};

			alignment_split_list make_list_of_alignment_split(const alignment &,
			                                                  const size_vec &);
			alignment_split_list get_all_single_alignment_splits(const alignment &);
			alignment_split_list get_preexisting_alignment_splits(const alignment &);
			alignment_split_list get_standard_alignment_splits(const alignment &);
			void add_alignment_splits(alignment_split_list &,
			                          const alignment_split_list &);
			alignment_split_list add_alignment_splits_copy(alignment_split_list,
			                                               const alignment_split_list &);
		} // namespace detail
	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_LIST_HPP
