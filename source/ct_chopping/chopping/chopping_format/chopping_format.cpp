/// \file
/// \brief The chopping_format class definitions

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

#include "chopping_format.hpp"

using namespace cath::common;
using namespace cath::chop;

using std::string;

/// \brief TODOCUMENT
///
/// \relates chopping_format
domain cath::chop::parse_domain(const chopping_format &prm_chopping_format,        ///< TODOCUMENT
                                const string          &prm_domain_chopping_string, ///< TODOCUMENT
                                const string          &prm_domain_id               ///< TODOCUMENT
                                ) {
	domain new_domain = prm_chopping_format.parse_domain( prm_domain_chopping_string );
	new_domain.set_opt_domain_id( prm_domain_id );
	return new_domain;
}
