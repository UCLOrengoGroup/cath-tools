/// \file
/// \brief The html_ostream_alignment_outputter class definitions

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


#include "html_ostream_alignment_outputter.h"

#include "alignment/alignment_context.h"
#include "alignment/io/alignment_io.h"
#include "alignment/io/outputter/html_align_outputter.h"
#include "common/clone/make_uptr_clone.h"
#include "display/display_colour/display_colour.h"
#include "display/display_spec/display_spec.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_outputter> html_ostream_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void html_ostream_alignment_outputter::do_output_alignment(const alignment_context &arg_alignment_context, ///< TODOCUMENT
                                                           ostream                 &arg_ostream            ///< TODOCUMENT
                                                           ) const {
	arg_ostream << html_align_outputter(
		arg_alignment_context.get_alignment(),
		arg_alignment_context.get_pdbs(),
		arg_alignment_context.get_names(),
		*colourer_ptr
	) << flush;
}

/// \brief TODOCUMENT
bool html_ostream_alignment_outputter::do_involves_display_spec() const {
	return true;
}

/// \brief Ctor for html_ostream_alignment_outputter
html_ostream_alignment_outputter::html_ostream_alignment_outputter(const display_colourer &arg_colourer ///< TODOCUMENT
                                                                   ) : colourer_ptr( arg_colourer.clone() ) {
}

/// \brief TODOCUMENT
html_ostream_alignment_outputter cath::opts::make_html_ostream_alignment_outputter(const display_spec &arg_display_spec ///< TODOCUMENT
                                                                                   ) {
	return html_ostream_alignment_outputter( *arg_display_spec.get_display_colourer( make_default_light_colour_gradient() ) );
}
