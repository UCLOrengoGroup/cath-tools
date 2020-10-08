/// \file
/// \brief The viewer class definitions

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

#include "viewer.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>

#include "alignment/alignment_context.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"
#include "display/display_colourer/display_colourer.hpp"
#include "display/display_colourer/display_colourer_consecutive.hpp"
#include "display/options/display_spec.hpp"
#include "display/viewer/pymol/pymol_tools.hpp"
#include "display_colour/display_colour.hpp"
#include "display_colour/display_colour_list.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/superposition_context.hpp"

#include <fstream>
#include <sstream>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::file;
using namespace cath::sup;

using boost::algorithm::is_space;
using boost::filesystem::path;
using boost::format;
using boost::is_alnum;
using boost::lexical_cast;
using boost::string_ref;
using std::flush;
using std::max;
using std::ofstream;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief Default is to no accept multiple colourings
bool viewer::do_accepts_multiple_colourings() const {
	return false;
}

/// \brief Default to writing no commands to the specified ostream when beginning a colouring with the specified colourer
void viewer::do_begin_colouring(ostream                &/*prm_ostream*/, ///< The ostream to which the PyMOL commands should be written
                                const display_colourer &/*prm_colourer*/ ///< The display_colourer to be used for the colouring that is beginning
                                ) {
}

/// \brief Default to writing no commands to the specified ostream when ending a colouring with the specified colourer
void viewer::do_end_colouring(ostream                &/*prm_ostream*/, ///< The ostream to which the PyMOL commands should be written
                              const display_colourer &/*prm_colourer*/ ///< The display_colourer to be used for the colouring that is ending
                              ) {
}

/// \brief TODOCUMENT
string viewer::default_executable() const {
	return do_default_executable();
}

/// \brief TODOCUMENT
string viewer::default_file_extension() const {
	return do_default_file_extension();
}

/// \brief TODOCUMENT
void viewer::write_start(ostream &prm_ostream ///< TODOCUMENT
                         ) const {
	do_write_start(prm_ostream);
}

/// \brief TODOCUMENT
void viewer::write_load_pdbs(ostream             &prm_ostream,       ///< TODOCUMENT
                             const superposition &prm_superposition, ///< TODOCUMENT
                             const pdb_list      &prm_pdbs,          ///< TODOCUMENT
                             const str_vec       &prm_names          ///< TODOCUMENT
                             ) const {
	return do_write_load_pdbs(
		prm_ostream,
		prm_superposition,
		prm_pdbs,
		prm_names
	);
}

/// \brief TODOCUMENT
void viewer::define_colour(ostream              &prm_ostream,    ///< TODOCUMENT
                           const display_colour &prm_colour,     ///< TODOCUMENT
                           const string         &prm_colour_name ///< TODOCUMENT
                           ) const {
	do_define_colour( prm_ostream, prm_colour, prm_colour_name );
}

/// \brief Get a string for colouring the base (ie everything) in the colour that has previously been defined with the specified name
string viewer::get_colour_base_str(const string &prm_colour_name ///< The previously-defined colour with which to colour the base
                                   ) const {
	return do_get_colour_base_str( prm_colour_name );
}

/// \brief TODOCUMENT
string viewer::get_colour_pdb_str(const string &prm_colour_name, ///< TODOCUMENT
                                  const string &prm_pdb_name     ///< TODOCUMENT
                                  ) const {
	return do_get_colour_pdb_str( prm_colour_name, prm_pdb_name );
}

/// \brief TODOCUMENT
string viewer::get_colour_pdb_residues_str(const string         &prm_colour_name, ///< TODOCUMENT
                                           const string         &prm_pdb_name,    ///< TODOCUMENT
                                           const residue_id_vec &prm_residues     ///< TODOCUMENT
                                           ) const {
	return do_get_colour_pdb_residues_str( prm_colour_name, prm_pdb_name, prm_residues );
}

/// \brief TODOCUMENT
void viewer::write_alignment_extras(ostream                     &prm_ostream,              ///< TODOCUMENT
                                    const superposition_context &prm_superposition_context ///< TODOCUMENT
                                    ) const {
	if ( ! prm_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write alignment extras for superposition_context that doesn't contain an alignment"));
	}
	return do_write_alignment_extras( prm_ostream, prm_superposition_context );
}

/// \brief TODOCUMENT
void viewer::write_end(ostream          &prm_os,  ///< TODOCUMENT
                       const string_ref &prm_name ///< TODOCUMENT
                       ) const {
	do_write_end( prm_os, prm_name );
}

/// \brief Strip
///
/// \relates viewer
///
/// Current heuristic:
///  - Keep alpha-numeric characters as they are
///  - Convert all space characters to underscores
///  - Drop all other characters
///
/// \todo Find a better way to do this. It doesn't seem immediately obvious that replace_all()
///       in Boost string algorithm will do the trick. Does it require regexps?
string cath::clean_name_for_viewer(const string &prm_name
                                   ) {
	string new_name;
	new_name.reserve(prm_name.size());
	for (const char &name_char : prm_name) {
		if (is_alnum()(name_char)) {
			new_name.push_back(name_char);
		}
		else if ( is_space()( name_char ) || name_char == '_' ) {
			new_name.push_back('_');
		}
	}
	return new_name;
}

/// \brief TODOCUMENT
///
/// \relates viewer
str_vec cath::clean_names_for_viewer(const str_vec &prm_names
                                     ) {
	str_vec new_names;
	new_names.reserve(prm_names.size());
	for (const string &name : prm_names) {
		new_names.push_back(clean_name_for_viewer(name));
	}
	return new_names;
}

/// \brief TODOCUMENT
///
/// \relates viewer
str_vec cath::clean_names_for_viewer(const superposition_context &prm_superposition_context ///< TODOCUMENT
                                     ) {
	return clean_names_for_viewer( get_viewer_names( get_name_sets( prm_superposition_context ) ) );
}

/// \brief TODOCUMENT
///
/// \relates viewer
str_vec cath::clean_names_for_viewer(const alignment_context &prm_alignment_context ///< TODOCUMENT
                                     ) {
	return clean_names_for_viewer( get_viewer_names( get_name_sets( prm_alignment_context ) ) );
}

/// \brief Output instructions from the specified viewer for the specified superposition_context to
///        the specified ostream, using the specified display_spec and only_warn flat
///
/// \relates viewer
void cath::output_superposition_to_viewer(ostream                          &prm_ostream,                 ///< The ostream to which the data should be written
                                          viewer                           &prm_viewer,                  ///< The viewer defining the instructions to be written
                                          const display_spec               &prm_display_spec,            ///< The specification for how to display the superposition
                                          const superposition_context      &prm_superposition_context,   ///< The superposition_context to output
                                          const superposition_content_spec &prm_content_spec,            ///< The specification of what should be included in the superposition
                                          const missing_aln_policy         &prm_missing_aln_policy,      ///< Whether to warn or throw if no alignment is present
                                          const string_ref                 &prm_name                     ///< A name for the superposition (to write out so users of the viewer know what they're looking at)
                                          ) {
	// Write the start of the viewer output
	prm_viewer.write_start(prm_ostream);

	// Write the text to load the PDBs
	prm_viewer.write_load_pdbs(
		prm_ostream,
		prm_superposition_context.get_superposition(),
		get_supn_content_pdbs( prm_superposition_context, prm_content_spec ),
		clean_names_for_viewer( prm_superposition_context )
	);

	const bool spec_is_consecutive            = is_consecutive( prm_display_spec );
	const bool spec_requires_alignment        = requires_alignment( prm_display_spec );
	const bool supn_has_alignment             = prm_superposition_context.has_alignment();
	const bool missing_wanted_alignment       = spec_requires_alignment && ! supn_has_alignment;
	const bool would_accept_extra_consecutive = prm_viewer.accepts_multiple_colourings() && ! spec_is_consecutive;

	if ( missing_wanted_alignment || would_accept_extra_consecutive) {
		if ( missing_wanted_alignment ) {
			const auto message = "Unable to apply an alignment-based colouring scheme to the superposition because it doesn't contain an alignment";
			if ( prm_missing_aln_policy == missing_aln_policy::WARN_AND_COLOUR_CONSECUTIVELY ) {
				BOOST_LOG_TRIVIAL( warning ) << message;
			}
			else {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(message));
			}
		}
		const display_colourer_consecutive the_colourer{ get_colour_list( prm_display_spec ) };

		colour_viewer(
			the_colourer,
			prm_ostream,
			prm_viewer,
			get_pdbs              ( prm_superposition_context ),
			clean_names_for_viewer( prm_superposition_context ),
			get_regions           ( prm_superposition_context )
		);
	}

	if ( prm_superposition_context.has_alignment() ) {
		// Apply the colour
		const unique_ptr<const display_colourer> display_colourer_ptr = get_display_colourer( prm_display_spec );

		colour_viewer(
			*display_colourer_ptr,
			prm_ostream,
			prm_viewer,
			make_restricted_alignment_context( prm_superposition_context )
		);

		// If there is an alignment then do magic with it
		// if ( prm_superposition_context.has_alignment() ) {
		prm_viewer.write_alignment_extras( prm_ostream, prm_superposition_context );
	}

	// Write the start of the viewer output
	prm_viewer.write_end( prm_ostream, prm_name );
}

/// \brief Output instructions from the specified viewer for the specified superposition_context to
///        the specified file, using the specified display_spec and only_warn flat
///
/// \relates viewer
void cath::output_superposition_to_viewer_file(const path                       &prm_out_file,                ///< The ostream to which the data should be written
                                               viewer                           &prm_viewer,                  ///< The viewer defining the instructions to be written
                                               const display_spec               &prm_display_spec,            ///< The specification for how to display the superposition
                                               const superposition_context      &prm_superposition_context,   ///< The superposition_context to output
                                               const superposition_content_spec &prm_content_spec,            ///< The specification of what should be included in the superposition
                                               const missing_aln_policy         &prm_missing_aln_policy,      ///< Whether to warn or throw if no alignment is present
                                               const string_ref                 &prm_name                     ///< A name for the superposition (to write out so users of the viewer know what they're looking at)
                                               ) {
	ofstream pymol_file_ostream;
	open_ofstream( pymol_file_ostream, prm_out_file );

	output_superposition_to_viewer(
		pymol_file_ostream,
		prm_viewer,
		prm_display_spec,
		prm_superposition_context,
		prm_content_spec,
		prm_missing_aln_policy,
		prm_name
	);
	pymol_file_ostream << flush;
	pymol_file_ostream.close();

}

/// \brief The name of the base-colour
string cath::base_colour_name() {
	return "base_colour";
}

/// \brief Generate a name to use in the viewer for the specified colour index
///        in the specified number of colours
string cath::generate_colour_name(const size_t          &prm_colour_index,   ///< The index of the colour to name
                                  const size_t          &prm_num_colours,    ///< The total number of colours
                                  const colour_category &prm_colour_category ///< The category of colouring (structure-only or structure-or-residue)
                                  ) {
	const size_t num_width  = lexical_cast<string>( max( 1_z, prm_num_colours ) - 1 ).length();
	const string format_str = R"(%0)" + ::std::to_string( num_width ) + "d";

	return "ct_"
		+ ( prm_colour_category == colour_category::STRUC_ONLY ? "strc"s : "strc_or_res"s )
		+ "_colr_"
		+ ( format( format_str ) % prm_colour_index ).str();
}

/// \brief Generate names to use in the viewer for the specified number of colours
str_vec cath::generate_colour_names(const size_t          &prm_num_colours,    ///< The total number of colours
                                    const colour_category &prm_colour_category ///< The category of colouring (structure-only or structure-or-residue)
                                    ) {
	return transform_build<str_vec>(
		indices( prm_num_colours ),
		[&] (const size_t &x) { return generate_colour_name(
			x,
			prm_num_colours,
			prm_colour_category
		); }
	);
}
