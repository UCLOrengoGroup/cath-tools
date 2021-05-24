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

#ifndef _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATION_HPP
#define _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATION_HPP

#include <boost/operators.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/biocore/residue_name.hpp"
#include "cath/chopping/residue_location/residue_locating.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::chop {

	/// \brief TODOCUMENT
	class residue_location final : private boost::totally_ordered<residue_location> {
		/// \brief TODOCUMENT
		residue_name_opt the_residue_name;

		/// \brief TODOCUMENT
		size_opt         residue_index;

	public:
		explicit residue_location(const residue_name &);
		residue_location(const residue_name &,
		                 const size_t &);
		explicit residue_location(const size_t &);

		[[nodiscard]] const residue_name_opt &get_opt_residue_name() const;
		[[nodiscard]] const size_opt &        get_opt_residue_index() const;
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

} // namespace cath::chop

#endif // _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATION_HPP
