/// \file
/// \brief The superposition_outputter class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "superposition_outputter.h"

#include "common/clone/check_uptr_clone_against_this.h"

#include <cassert>
#include <typeinfo>

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<superposition_outputter> superposition_outputter::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
void superposition_outputter::output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
                                                   ostream &arg_ostream ///< TODOCUMENT
                                                   ) const {
	do_output_superposition(
		arg_superposition_context,
		arg_ostream
	);
}

/// \brief TODOCUMENT
bool superposition_outputter::involves_display_spec() const {
	return do_involves_display_spec();
}
