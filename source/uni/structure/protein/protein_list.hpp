/// \file
/// \brief The protein_list class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_LIST_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_LIST_H

#include "common/type_aliases.hpp"
#include "structure/structure_type_aliases.hpp"

#include <cstddef>

namespace cath {
	class protein;

	/// \brief TODOCUMENT
	class protein_list final {
	private:
		/// \brief TODOCUMENT
		protein_vec proteins;

	public:
		void push_back(const protein &);
		void reserve(const size_t &);

		size_t size() const noexcept;
		size_t max_size() const noexcept;
		bool empty() const;

		protein & operator[](const size_t &);
		const protein & operator[](const size_t &) const;

		using iterator        = protein_vec::iterator;
		using const_iterator  = protein_vec::const_iterator;
		using const_pointer   = protein_vec::const_pointer;
		using const_reference = protein_vec::const_reference;
		using difference_type = protein_vec::difference_type;
		using size_type       = protein_vec::size_type;
		using value_type      = protein_vec::value_type;

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;
	};

	protein_list make_protein_list(const protein_vec &);

	protein_list make_subset_protein_list(const protein_list &,
	                                      const size_vec &);

	amino_acid_vec_vec get_amino_acid_lists(const protein_list &);

	size_vec get_protein_lengths(const protein_list &);

	size_size_pair min_max_protein_length(const protein_list &);
} // namespace cath

#endif
