/// \file
/// \brief The residue_location class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATION_H
#define _CATH_TOOLS_SOURCE_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATION_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "biocore/residue_name.hpp"
#include "chopping/residue_location/residue_locating.hpp"
#include "common/type_aliases.hpp"
#include "structure/structure_type_aliases.hpp"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class residue_location final : private boost::totally_ordered<residue_location> {
			/// \brief TODOCUMENT
			residue_name_opt the_residue_name;

			/// \brief TODOCUMENT
			size_opt         residue_index;

		public:
			residue_location(const residue_name &);
			residue_location(const residue_name &,
			                 const size_t &);
			residue_location(const size_t &);

			const residue_name_opt & get_opt_residue_name() const;
			const size_opt & get_opt_residue_index() const;
		};

		bool has_residue_name(const residue_location &);
		bool has_residue_index(const residue_location &);

		residue_name get_residue_name(const residue_location &);
		size_t get_residue_index(const residue_location &);

		residue_locating get_residue_locating(const residue_location &);

		bool operator==(const residue_location &,
		                const residue_location &);
		bool operator<(const residue_location &,
		               const residue_location &);
	} // namespace chop
} // namespace cath

#endif
