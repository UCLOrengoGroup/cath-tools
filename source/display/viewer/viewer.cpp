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

#include "viewer.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>

#include "alignment/alignment_context.h"
#include "display/display_colourer/display_colourer.h"
#include "display/options/display_spec.h"
#include "display/viewer/pymol/pymol_tools.h"
#include "display_colour/display_colour.h"
#include "display_colour/display_colour_list.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/superposition_context.h"

#include <sstream>

using namespace boost::algorithm;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sup;
using namespace std;

using boost::algorithm::is_space;
using boost::lexical_cast;

/// \brief TODOCUMENT
string viewer::default_executable() const {
	return do_default_executable();
}

/// \brief TODOCUMENT
string viewer::default_file_extension() const {
	return do_default_file_extension();
}

/// \brief TODOCUMENT
void viewer::write_start(ostream &arg_ostream ///< TODOCUMENT
                         ) const {
	do_write_start(arg_ostream);
}

/// \brief TODOCUMENT
void viewer::write_load_pdbs(ostream             &arg_ostream,       ///< TODOCUMENT
                             const superposition &arg_superposition, ///< TODOCUMENT
                             const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                             const str_vec       &arg_names          ///< TODOCUMENT
                             ) const {
	return do_write_load_pdbs(
		arg_ostream,
		arg_superposition,
		arg_pdbs,
		arg_names
	);
}

/// \brief TODOCUMENT
void viewer::define_colour(ostream             &arg_ostream,    ///< TODOCUMENT
                           const display_colour &arg_colour,     ///< TODOCUMENT
                           const string        &arg_colour_name ///< TODOCUMENT
                           ) const {
	do_define_colour(arg_ostream, arg_colour, arg_colour_name);
}

/// \brief TODOCUMENT
string viewer::get_colour_pdb_str(const string &arg_colour_name, ///< TODOCUMENT
                                  const string &arg_pdb_name     ///< TODOCUMENT
                                  ) const {
	return do_get_colour_pdb_str( arg_colour_name, arg_pdb_name );
}

/// \brief TODOCUMENT
string viewer::get_colour_pdb_residues_str(const string           &arg_colour_name, ///< TODOCUMENT
                                           const string           &arg_pdb_name,    ///< TODOCUMENT
                                           const residue_name_vec &arg_residues     ///< TODOCUMENT
                                           ) const {
	return do_get_colour_pdb_residues_str( arg_colour_name, arg_pdb_name, arg_residues );
}

/// \brief TODOCUMENT
void viewer::write_alignment_extras(ostream                     &arg_ostream,              ///< TODOCUMENT
                                    const superposition_context &arg_superposition_context ///< TODOCUMENT
                                    ) const {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot write alignment extras for superposition_context that doesn't contain an alignment"));
	}
	return do_write_alignment_extras( arg_ostream, arg_superposition_context );
}

/// \brief TODOCUMENT
void viewer::write_end(ostream &arg_ostream ///< TODOCUMENT
                       ) const {
	do_write_end(arg_ostream);
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
string cath::clean_name_for_viewer(const string &arg_name
                                   ) {
	string new_name;
	new_name.reserve(arg_name.size());
	for (const char &name_char : arg_name) {
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
str_vec cath::clean_names_for_viewer(const str_vec &arg_names
                                     ) {
	str_vec new_names;
	new_names.reserve(arg_names.size());
	for (const string &name : arg_names) {
		new_names.push_back(clean_name_for_viewer(name));
	}
	return new_names;
}

/// \brief TODOCUMENT
///
/// \relates viewer
str_vec cath::clean_names_for_viewer(const superposition_context &arg_superposition_context ///< TODOCUMENT
                                     ) {
	return clean_names_for_viewer( arg_superposition_context.get_names_cref() );
}

/// \brief TODOCUMENT
///
/// \relates viewer
str_vec cath::clean_names_for_viewer(const alignment_context &arg_alignment_context ///< TODOCUMENT
                                     ) {
	return clean_names_for_viewer( arg_alignment_context.get_names() );
}

/// \brief TODOCUMENT
///
/// \relates viewer
void cath::output_superposition_to_viewer(ostream                     &arg_ostream,              ///< TODOCUMENT
                                          const viewer                &arg_viewer,               ///< TODOCUMENT
                                          const display_spec          &arg_display_spec,         ///< TODOCUMENT
                                          const superposition_context &arg_superposition_context ///< TODOCUMENT
                                          ) {
	// Write the start of the viewer output
	arg_viewer.write_start(arg_ostream);

	// Write the text to load the PDBs
	arg_viewer.write_load_pdbs(
		arg_ostream,
		arg_superposition_context.get_superposition_cref(),
		arg_superposition_context.get_pdbs_cref(),
		clean_names_for_viewer( arg_superposition_context )
	);

	// Apply the colour
	const unique_ptr<const display_colourer> display_colourer_ptr = get_display_colourer( arg_display_spec );

//	colour_alignment(
//		*display_colourer_ptr,
//		cerr,
//		arg_superposition_context
//	);

	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("WARNING: Unable to colour the superposition because it does not contain an alignment"));
	}

	colour_viewer(
		*display_colourer_ptr,
		arg_ostream,
		arg_viewer,
		make_alignment_context( arg_superposition_context )
	);

	// If there is an alignment then do magic with it
	if ( arg_superposition_context.has_alignment() ) {
		arg_viewer.write_alignment_extras( arg_ostream, arg_superposition_context );
	}

	// Write the start of the viewer output
	arg_viewer.write_end(arg_ostream);
}

