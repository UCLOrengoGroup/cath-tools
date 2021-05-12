/// \file
/// \brief The domain_definition class definitions

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

#include "domain_definition.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::std::string;

/// \brief Ctor for domain_definition
domain_definition::domain_definition(domain prm_domain,  ///< TODOCUMENT
                                     string prm_pdb_name ///< TODOCUMENT
                                     ) : the_domain{ std::move( prm_domain   ) },
                                         pdb_name  { std::move( prm_pdb_name ) } {
}

/// \brief TODOCUMENT
const domain & domain_definition::get_domain() const {
	return the_domain;
}

/// \brief TODOCUMENT
const string & domain_definition::get_pdb_name() const {
	return pdb_name;
}
