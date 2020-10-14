/// \file
/// \brief The domall_chopping_format class definitions

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

#include "domall_chopping_format.hpp"

#include "chopping/domain/domain.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

#include <iostream> // ***** TEMPORARY *****

using namespace cath::chop;
using namespace cath::common;

using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<chopping_format> domall_chopping_format::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool domall_chopping_format::do_represents_fragments() const {
	return true;
}

domain domall_chopping_format::do_parse_domain(const string &prm_domain_chopping_string ///< TODOCUMENT
                                               ) const {
	std::cerr << "domain_chopping_string is " << prm_domain_chopping_string << "\n";

	BOOST_THROW_EXCEPTION(not_implemented_exception("domall_chopping_format::do_parse_domain()"));

	return domain( region_vec() );
}

/// \brief Concrete definition of this chopping_format writes a region to a string
string domall_chopping_format::do_write_region(const region &/*prm_region*/ ///< The region to write to a string
                                               ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("domall_chopping_format cannot currently write regions"));
}

/// \brief Concrete definition of this chopping_format writes a domain to a string
string domall_chopping_format::do_write_domain(const domain &/*prm_domain*/ ///< The domain to write to a string
                                               ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("domall_chopping_format cannot currently write domains"));
}
