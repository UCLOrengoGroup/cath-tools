/// \file
/// \brief The alignment_context class definitions

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

#include "alignment_context.hpp"

#include "chopping/region/region.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/superposition_context.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::chop;
using namespace cath::file;
using namespace cath::sup;

/// \brief Ctor for alignment_context
alignment_context::alignment_context(const alignment          &arg_alignment, ///< TODOCUMENT
                                     const strucs_context     &arg_context    ///< TODOCUMENT
                                     ) : the_alignment { arg_alignment },
                                         context       { arg_context   } {
}

/// \brief Ctor for alignment_context
alignment_context::alignment_context(const alignment          &arg_alignment, ///< TODOCUMENT
                                     const pdb_list           &arg_pdbs,      ///< TODOCUMENT
                                     const str_vec            &arg_names,     ///< TODOCUMENT
                                     const region_vec_opt_vec &arg_regions    ///< The specification of the regions of the PDBs to which the alignment refers
                                     ) : the_alignment { arg_alignment                    },
                                         context       { arg_pdbs, arg_names, arg_regions } {
}

/// \brief TODOCUMENT
const alignment & alignment_context::get_alignment() const {
	return the_alignment;
}

/// \brief TODOCUMENT
const strucs_context & alignment_context::get_strucs_context() const {
	return context;
}

/// \brief TODOCUMENT
const pdb_list & cath::align::get_pdbs(const alignment_context &arg_align_context ///< TODOCUMENT
                                       ) {
	return arg_align_context.get_strucs_context().get_pdbs();
}

/// \brief TODOCUMENT
const str_vec & cath::align::get_names(const alignment_context &arg_align_context ///< TODOCUMENT
                                       ) {
	return arg_align_context.get_strucs_context().get_names();
}

/// \brief Getter for the specification of the regions of the PDBs to which the alignment refers
const region_vec_opt_vec & cath::align::get_regions(const alignment_context &arg_align_context ///< TODOCUMENT
                                                    ) {
	return arg_align_context.get_strucs_context().get_regions();
}

/// \brief TODOCUMENT
///
/// \relates alignment_context
///
/// \relatesalso superposition_context
superposition_context cath::align::make_superposition_context(const alignment_context &arg_alignment_context, ///< TODOCUMENT
                                                              const superposition     &arg_superposition      ///< TODOCUMENT
                                                              ) {
	return superposition_context(
		arg_superposition,
		arg_alignment_context.get_strucs_context(),
		arg_alignment_context.get_alignment()
	);
}
