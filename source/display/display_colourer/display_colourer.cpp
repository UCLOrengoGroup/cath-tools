/// \file
/// \brief The display_colourer class definitions

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

#include "display_colourer.h"

#include "alignment/alignment_context.h"
#include "common/clone/check_uptr_clone_against_this.h"
#include "common/cpp14/make_unique.h"
#include "display/display_colourer/detail/score_colour_handler.h"
#include "display/display_colourer/display_colour_spec.h"
#include "display/display_colourer/display_colourer_alignment.h"
#include "display/display_colourer/display_colourer_consecutive.h"
#include "display/options/display_spec.h"
#include "display/viewer/viewer.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/superposition_context.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::file;
using namespace cath::sup;

using std::ostream;
using std::unique_ptr;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<display_colourer> display_colourer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}


/// \brief TODOCUMENT
display_colour_spec display_colourer::get_colour_spec(const alignment_context &arg_alignment_context ///< TODOCUMENT
                                                      ) const {
	return do_get_colour_spec( arg_alignment_context );
}

/// \brief TODOCUMENT
///
/// \relates display_colourer
display_colour_spec cath::get_colour_spec(const display_colourer &arg_colourer, ///< TODOCUMENT
                                          const pdb_list         &arg_pdbs,     ///< TODOCUMENT
                                          const str_vec          &arg_names,    ///< TODOCUMENT
                                          const alignment        &arg_alignment ///< TODOCUMENT
                                          ) {
	const alignment::size_type num_entries   = arg_alignment.num_entries();
	const alignment::size_type aln_length    = arg_alignment.length();

	if ( aln_length <= 0 || num_entries <= 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to colour the alignment_context because the alignment is empty"));
	}
	if ( num_entries != arg_pdbs.size()  ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to colour the alignment_context because the number of entries doesn't match the number of PDBs"));
	}
	if ( num_entries != arg_names.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to colour the alignment_context because the number of entries doesn't match the number of names"));
	}
	return arg_colourer.get_colour_spec( alignment_context(
		arg_pdbs,
		arg_names,
		arg_alignment
	) );
}

/// \brief TODOCUMENT
unique_ptr<const display_colourer> cath::get_display_colourer(const display_spec            &arg_display_spec,   ///< TODOCUMENT
                                                              const display_colour_gradient &arg_colour_gradient ///< TODOCUMENT
                                                              ) {
	const score_colour_handler colour_handler{
		arg_display_spec.get_show_scores_if_present(),
		arg_display_spec.get_scores_to_equivs(),
		arg_display_spec.get_normalise_scores()
	};
	if ( arg_display_spec.get_gradient_colour_alignment() ) {
		return { common::make_unique< display_colourer_alignment   >( colour_handler, arg_colour_gradient                 ) };
	}
	else {
		return { common::make_unique< display_colourer_consecutive >( colour_handler, get_colour_list( arg_display_spec ) ) };
	}
}

/// \brief TODOCUMENT
///
/// \relates display_colourer
void cath::colour_viewer(const display_colourer  &arg_colourer, ///< TODOCUMENT
                         ostream                 &arg_os,       ///< TODOCUMENT
                         const viewer            &arg_viewer,   ///< TODOCUMENT
                         const alignment_context &arg_aln_con   ///< TODOCUMENT
                         ) {
	const display_colour_spec the_spec = arg_colourer.get_colour_spec( arg_aln_con );
	colour_viewer_with_spec( the_spec, arg_viewer, arg_aln_con, arg_os );
}

///// \brief TODOCUMENT
/////
///// \relates display_colourer
//void cath::colour_alignment(const display_colourer      &arg_colourer, ///< TODOCUMENT
//                            ostream                     &arg_os,       ///< TODOCUMENT
//                            const superposition_context &arg_sup_con   ///< TODOCUMENT
//                            ) {
//	const display_colour_spec the_spec = arg_colourer.get_colour_spec( arg_sup_con );
//	colour_alignment_with_spec(
//		the_spec,
//		arg_sup_con.get_alignment_cref(),
//		arg_sup_con.get_pdbs_cref(),
//		clean_names_for_viewer( arg_sup_con ),
//		arg_os
//	);
//}
