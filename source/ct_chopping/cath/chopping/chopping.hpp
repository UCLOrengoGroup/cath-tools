/// \file
/// \brief The chopping class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_HPP
#define CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_HPP

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"

#include <cstddef>
#include <vector>

namespace cath::chop {

	/// \brief TODOCUMENT
	class chopping final {
	private:
		/// \brief TODOCUMENT
		domain_vec domains;

		/// \brief TODOCUMENT
		region_vec fragments;

		void sanity_check() const;

	public:
		using iterator = domain_vec::iterator;
		using const_iterator = domain_vec::const_iterator;

		explicit chopping(domain_vec,
		                  region_vec = region_vec());

		[[nodiscard]] size_t num_domains() const;
		[[nodiscard]] size_t num_fragments() const;

		[[nodiscard]] const region &get_fragment_of_index( const size_t & ) const;

		const domain & operator[](const size_t &prm_index) const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

} // namespace cath::chop

#endif // CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_HPP
