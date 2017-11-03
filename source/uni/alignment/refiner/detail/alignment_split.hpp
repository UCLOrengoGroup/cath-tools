/// \file
/// \brief The alignment_split class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_H

#include <boost/operators.hpp>

#include "alignment/refiner/detail/alignment_split_half.hpp"
#include "common/type_aliases.hpp"

#include <cstddef>

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			class alignment_split final : private boost::equivalent      < alignment_split,
			                                      boost::totally_ordered < alignment_split > > {
			private:
				/// \brief TODOCUMENT
				size_t num_entries;

				/// \brief TODOCUMENT
				size_set first_half_entries;

			public:
				explicit alignment_split(const size_t &);

				size_t get_num_entries() const;
				void add_first_half_entry(const size_t &);

				size_t get_num_first_half_entries() const;
				const size_set & get_first_half_entries() const;

				using const_iterator = size_set::const_iterator;
				const_iterator begin() const;
				const_iterator end() const;
			};

			bool operator<(const alignment_split &,
			               const alignment_split &);

			size_set entries_of_alignment_split_half(const alignment_split &,
			                                         const alignment_split_half &);

			alignment_split make_opposite_version(const alignment_split &);

			alignment_split get_least_version(const alignment_split &);

			bool is_valid_split(const alignment_split &);

			alignment_split make_single_alignment_split(const size_t &,
			                                            const size_t &);

			alignment_split make_alignment_split(const size_vec &,
			                                     const size_t &);

//			size_t get_split_length(const alignment &,
//			                        const protein_list &,
//			                        const alignment_split &,
//			                        const alignment_split_half &);

		} // namespace detail
	} // namespace align
} // namespace cath

#endif
