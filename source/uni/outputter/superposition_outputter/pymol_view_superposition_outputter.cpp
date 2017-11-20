/// \file
/// \brief The pymol_view_superposition_outputter class definitions

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

#include "pymol_view_superposition_outputter.hpp"

#include <boost/log/trivial.hpp>

#include "common/clone/make_uptr_clone.hpp"
#include "common/command_executer.hpp"
#include "common/file/temp_file.hpp"
#include "outputter/superposition_outputter/pymol_file_superposition_outputter.hpp"

#include <iostream>

using namespace boost::log;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;

using boost::filesystem::path;
using boost::string_ref;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pymol_view_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pymol_view_superposition_outputter::do_output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
                                                                 ostream                     &arg_ostream,               ///< TODOCUMENT
                                                                 const string_ref            &arg_name                   ///< A name for the superposition (so users of the superposition know what it represents)
                                                                 ) const {
	const temp_file pymol_script_filename(".%%%%-%%%%-%%%%-%%%%.pml");
	const pymol_file_superposition_outputter pymol_file_outputter(
		get_filename( pymol_script_filename ),
		the_display_spec,
		content_spec
	);
	pymol_file_outputter.output_superposition( arg_superposition_context, arg_ostream, arg_name );

	const bool pymol_success = command_executer::execute(
		pymol_program,
		{ get_filename( pymol_script_filename ).string() }
	);
	if ( ! pymol_success ) {
		BOOST_LOG_TRIVIAL( warning ) << "PyMOL executable " + pymol_program.string() + " did not run/shutdown normally.";
	}
}

/// \brief TODOCUMENT
bool pymol_view_superposition_outputter::do_involves_display_spec() const {
	return true;
}

/// \brief Getter for the name of this superposition_outputter
string pymol_view_superposition_outputter::do_get_name() const {
	return "pymol_view_superposition_outputter";
}

/// \brief Ctor for pymol_view_superposition_outputter
pymol_view_superposition_outputter::pymol_view_superposition_outputter(const path                 &arg_pymol_program, ///< TODOCUMENT
                                                                       display_spec                arg_display_spec,  ///< TODOCUMENT
                                                                       superposition_content_spec  arg_content_spec   ///< The specification of what should be included in the superposition
                                                                       ) : pymol_program    { arg_pymol_program             },
                                                                           the_display_spec { std::move( arg_display_spec ) },
                                                                           content_spec     { std::move( arg_content_spec ) } {
}
