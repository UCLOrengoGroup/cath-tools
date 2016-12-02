/// \file
/// \brief The jmol_selection_chopping_format class definitions

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

#include "jmol_selection_chopping_format.hpp"

#include "chopping/domain/domain.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

#include <iostream> // ***** TEMPORARY *****

using namespace cath::chop;
using namespace cath::common;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<chopping_format> jmol_selection_chopping_format::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool jmol_selection_chopping_format::do_represents_fragments() const {
	return false;
}

/// \brief TODOCUMENT
domain jmol_selection_chopping_format::do_parse_domain(const string &arg_domain_chopping_string ///< TODOCUMENT
                                                       ) const {
	cerr << "domain_chopping_string is " << arg_domain_chopping_string << endl;

	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_selection_chopping_format::do_parse_domain()"));

	return domain( region_vec() );
}
