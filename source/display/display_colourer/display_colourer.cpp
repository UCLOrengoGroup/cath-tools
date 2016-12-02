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

#include "display_colourer.hpp"

#include "alignment/alignment_context.hpp"
#include "common/clone/check_uptr_clone_against_this.hpp"
#include "common/cpp14/make_unique.hpp"
#include "display/display_colourer/detail/score_colour_handler.hpp"
#include "display/display_colour_spec/display_colour_spec.hpp"
#include "display/display_colourer/display_colourer_alignment.hpp"
#include "display/display_colourer/display_colourer_consecutive.hpp"
#include "display/options/display_spec.hpp"
#include "display/viewer/viewer.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/superposition_context.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::file;
using namespace cath::sup;

using std::ostream;
using std::unique_ptr;

/// \brief Ctor from specification for post-modifying the colouring based on scores
display_colourer::display_colourer(const score_colour_handler &arg_score_colour_handler ///< Specification for post-modifying the colouring based on scores
                                   ) : the_score_colour_handler{ arg_score_colour_handler } {
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<display_colourer> display_colourer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief Getter for the optional specification for post-modifying the colouring based on scores
const score_colour_handler_opt & display_colourer::get_score_colour_handler_opt() const {
	return the_score_colour_handler;
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
	auto &&result_spec = arg_colourer.get_colour_spec( alignment_context(
		arg_pdbs,
		arg_names,
		arg_alignment
	) );

	return has_score_colour_handler( arg_colourer )
		? adjust_display_colour_spec_copy(
			std::forward< decltype( result_spec ) >( result_spec ),
			get_score_colour_handler( arg_colourer ),
			arg_alignment
		)
		: result_spec;
}

/// \brief Return whether the specified display_colourer has a score_colour_handler
///
/// \relates display_colourer
bool cath::has_score_colour_handler(const display_colourer &arg_display_colourer ///< The display_colourer to query
                                    ) {
	return static_cast<bool>( arg_display_colourer.get_score_colour_handler_opt() );
}

/// \brief Get the specified display_colourer's score_colour_handler
///
/// \pre `has_score_colour_handler( arg_display_colourer )` else throws an invalid_argument_exception
///
/// \relates display_colourer
const score_colour_handler & cath::get_score_colour_handler(const display_colourer &arg_display_colourer ///< The display_colourer to query
                                                            ) {
	if ( ! has_score_colour_handler( arg_display_colourer ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot retrieve score_colour_handler from display_colourer which doesn't have one"));
	}
	return arg_display_colourer.get_score_colour_handler_opt().get();
}

/// \brief TODOCUMENT
///
/// \relates display_colourer
unique_ptr<const display_colourer> cath::get_display_colourer(const display_spec            &arg_display_spec,   ///< TODOCUMENT
                                                              const display_colour_gradient &arg_colour_gradient ///< TODOCUMENT
                                                              ) {
	const score_colour_handler colour_handler{
		arg_display_spec.get_show_scores_if_present(),
		arg_display_spec.get_scores_to_equivs(),
		arg_display_spec.get_normalise_scores()
	};
	if ( arg_display_spec.get_gradient_colour_alignment() ) {
		return { common::make_unique< display_colourer_alignment   >( arg_colour_gradient,                 colour_handler ) };
	}
	else {
		return { common::make_unique< display_colourer_consecutive >( get_colour_list( arg_display_spec ), colour_handler ) };
	}
}

/// \brief Write instructions for the specified viewer to the specified ostream
///        to represent the specified display_colourer in the context of the specified alignment_context
///
/// \relates display_colourer
void cath::colour_viewer(const display_colourer  &arg_colourer, ///< The display_colourer to write to the ostream
                         ostream                 &arg_os,       ///< The ostream to which the instructions should be written
                         const viewer            &arg_viewer,   ///< The viewer defining the instructions to be written
                         const alignment_context &arg_aln_con   ///< The alignment_context providing the context for the display_colourer
                         ) {
	colour_viewer_with_spec(
		arg_colourer.get_colour_spec( arg_aln_con ),
		arg_viewer,
		arg_aln_con,
		arg_os
	);
}

/// \brief Write instructions for the specified viewer to the specified ostream
///        to represent the specified alignment_free_display_colourer in the context of
///        the specified cleaned structure names
///
/// \relates display_colourer
void cath::colour_viewer(const alignment_free_display_colourer &arg_colourer,                ///< The alignment_free_display_colourer to write to the ostream
                         ostream                               &arg_os,                      ///< The ostream to which the instructions should be written
                         const viewer                          &arg_viewer,                  ///< The viewer defining the instructions to be written
                         const str_vec                         &arg_cleaned_names_for_viewer ///< The names of the structures (cleaned for the viewer)
                         ) {
	colour_viewer_with_spec(
		arg_colourer.get_colour_spec_from_num_entries( arg_cleaned_names_for_viewer.size() ),
		arg_viewer,
		arg_cleaned_names_for_viewer,
		arg_os
	);
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
