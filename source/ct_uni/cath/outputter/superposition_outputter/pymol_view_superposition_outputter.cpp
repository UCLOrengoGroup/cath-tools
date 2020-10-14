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

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/command_executer.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/outputter/run_pymol.hpp"
#include "cath/outputter/superposition_outputter/pymol_file_superposition_outputter.hpp"

#include <iostream>

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace cath::view;

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
void pymol_view_superposition_outputter::do_output_superposition(const superposition_context &prm_superposition_context, ///< TODOCUMENT
                                                                 ostream                     &prm_ostream                ///< TODOCUMENT
                                                                 ) const {
	const temp_file pymol_script_filename(".%%%%-%%%%-%%%%-%%%%.pml");
	const pymol_file_superposition_outputter pymol_file_outputter(
		get_filename( pymol_script_filename ),
		the_display_spec,
		content_spec
	);
	pymol_file_outputter.output_superposition( prm_superposition_context, prm_ostream );

	run_pymol(
		get_filename( pymol_script_filename ),
		pymol_program
	);
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
pymol_view_superposition_outputter::pymol_view_superposition_outputter(const path                 &prm_pymol_program, ///< TODOCUMENT
                                                                       display_spec                prm_display_spec,  ///< TODOCUMENT
                                                                       superposition_content_spec  prm_content_spec   ///< The specification of what should be included in the superposition
                                                                       ) : pymol_program    { prm_pymol_program             },
                                                                           the_display_spec { std::move( prm_display_spec ) },
                                                                           content_spec     { std::move( prm_content_spec ) } {
}
