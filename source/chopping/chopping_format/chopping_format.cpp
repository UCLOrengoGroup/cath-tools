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

#include "chopping/domain/domain.hpp"
#include "common/clone/check_uptr_clone_against_this.hpp"

//using namespace cath;
using namespace cath::common;
using namespace cath::chop;
using namespace std;

/// \brief TODOCUMENT
bool chopping_format::represents_fragments() const {
	return do_represents_fragments();
}

/// \brief TODOCUMENT
domain chopping_format::parse_domain(const string &arg_domain_chopping_string ///< TODOCUMENT
                                     ) const {
	return do_parse_domain( arg_domain_chopping_string );
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<chopping_format> chopping_format::clone() const {
 return check_uptr_clone_against_this( do_clone(), *this );
}


/// \brief TODOCUMENT
domain cath::chop::parse_domain(const chopping_format &arg_chopping_format,        ///< TODOCUMENT
                                const string          &arg_domain_chopping_string, ///< TODOCUMENT
                                const string          &arg_domain_id               ///< TODOCUMENT
                                ) {
	domain new_domain = arg_chopping_format.parse_domain( arg_domain_chopping_string );
	new_domain.set_opt_domain_id( arg_domain_id );
	return new_domain;
}



