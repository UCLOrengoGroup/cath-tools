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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_H

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/region/region.hpp"

#include <cstddef>
#include <vector>

namespace cath {
	namespace chop {

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

			explicit chopping(const domain_vec &,
			                  const region_vec & = region_vec() );

			size_t num_domains() const;
			size_t num_fragments() const;

			const region & get_fragment_of_index(const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

	} // namespace chop
} // namespace cath

#endif
