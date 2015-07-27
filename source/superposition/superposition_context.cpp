/// \file
/// \brief The superposition_context class definitions

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

#include "superposition_context.h"

#include "alignment/alignment_context.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sup;
using namespace std;

/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,         ///< TODOCUMENT
                                             const str_vec       &arg_names,        ///< TODOCUMENT
                                             const superposition &arg_superposition ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition) {
}


/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                                             const str_vec       &arg_names,         ///< TODOCUMENT
                                             const superposition &arg_superposition, ///< TODOCUMENT
                                             const alignment     &arg_alignment      ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition),
                                                 any_alignment(1, arg_alignment) {

}

/// \brief TODOCUMENT
const pdb_list & superposition_context::get_pdbs_cref() const {
	return pdbs;
}

/// \brief TODOCUMENT
const str_vec & superposition_context::get_names_cref() const {
	return names;
}

/// \brief TODOCUMENT
const superposition & superposition_context::get_superposition_cref() const {
	return the_superposition;
}

/// \brief TODOCUMENT
bool superposition_context::has_alignment() const {
	return static_cast<bool>( any_alignment );
}

/// \brief TODOCUMENT
const alignment & superposition_context::get_alignment_cref() const {
	if ( ! has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get alignment from superposition_context that doesn't contain one"));
	}
	return *any_alignment;
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
///
/// \relatesalso alignment_context
alignment_context cath::sup::make_alignment_context(const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                    ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("WARNING: Unable to extract alignment from superposition_context with no alignment"));
	}
	return alignment_context(
		arg_superposition_context.get_pdbs_cref(),
		arg_superposition_context.get_names_cref(),
		arg_superposition_context.get_alignment_cref()
	);
}
