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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_TYPE_ALIASES_H

#include <boost/config.hpp> /// \todo Come a resolution for Boost Trac tickets 12142 & 12179, remove this #include
#include <boost/optional/optional_fwd.hpp>

#include "chopping/residue_location/residue_locating.hpp"

#include <vector>

namespace cath {
	namespace chop {
		class domain;
		class domain_definition;
		class region;

		/// \brief Type alias for a vector of domains
		using domain_vec            = std::vector<domain>;

		/// \brief Type alias for a vector of regions
		using region_vec            = std::vector<region>;

		/// \brief Type alias for an optional region_vec
		using region_vec_opt        = boost::optional<region_vec>;

		/// \brief Type alias for a reference_wrapper of const region_vec
		using region_vec_cref       = std::reference_wrapper<const region_vec>;

		/// \brief Type alias for an optional region_vec_cref
		using region_vec_cref_opt   = boost::optional<region_vec_cref>;

		/// \brief Type alias for a vector of region_vecs
		using region_vec_vec        = std::vector<region_vec>;

		/// \brief Type alias for a vector of region_vec_opts
		using region_vec_opt_vec    = std::vector<region_vec_opt>;

		/// \brief Type alias for a vector of domain_definitions
		using domain_definition_vec = std::vector<domain_definition>;


		/// \brief Type alias for an optional residue_locating
		using residue_locating_opt  = boost::optional<residue_locating>;
	} // namespace chop
} // namespace cath

#endif
