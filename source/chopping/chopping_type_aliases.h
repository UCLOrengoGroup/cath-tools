/// \file
/// \brief The chopping type_aliases header

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

#ifndef CHOPPING_TYPE_ALIASES_H_INCLUDED
#define CHOPPING_TYPE_ALIASES_H_INCLUDED

#include <boost/config.hpp> /// \todo Come a resolution for Boost Trac tickets 12142 & 12179, remove this #include
#include <boost/optional/optional_fwd.hpp>

#include "chopping/residue_location/residue_locating.h"

#include <vector>

namespace cath {
	namespace chop {
		class domain;
		class domain_definition;
		class region;

    	/// \brief TODOCUMENT
		using domain_vec = std::vector<domain>;

    	/// \brief TODOCUMENT
		using region_vec = std::vector<region>;


		/// \brief TODOCUMENT
		using domain_definition_vec = std::vector<domain_definition>;


		/// \brief TODOCUMENT
		using opt_residue_locating = boost::optional<residue_locating>;
	}
}

#endif
