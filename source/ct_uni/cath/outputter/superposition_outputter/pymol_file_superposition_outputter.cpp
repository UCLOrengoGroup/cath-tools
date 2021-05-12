/// \file
/// \brief The pymol_file_superposition_outputter class definitions

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

#include "pymol_file_superposition_outputter.hpp"

#include <filesystem>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/display/viewer/pymol_viewer.hpp"
#include "cath/superposition/superposition_context.hpp"

using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::string_ref;
using ::std::filesystem::path;
using ::std::ostream;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pymol_file_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pymol_file_superposition_outputter::do_output_superposition(const superposition_context &prm_superposition_context, ///< TODOCUMENT
                                                                 ostream                     &/*prm_ostream*/            ///< TODOCUMENT
                                                                 ) const {
	pymol_viewer the_pymol_viewer{};
	output_superposition_to_viewer_file(
		output_file,
		the_pymol_viewer,
		the_display_spec,
		prm_superposition_context,
		content_spec,
		missing_aln_policy::WARN_AND_COLOUR_CONSECUTIVELY
	);
}

/// \brief TODOCUMENT
bool pymol_file_superposition_outputter::do_involves_display_spec() const {
	return true;
}

/// \brief Getter for the name of this superposition_outputter
string pymol_file_superposition_outputter::do_get_name() const {
	return "pymol_file_superposition_outputter";
}

/// \brief Ctor for pymol_file_superposition_outputter
pymol_file_superposition_outputter::pymol_file_superposition_outputter(const path                 &prm_output_file,  ///< TODOCUMENT
                                                                       display_spec                prm_display_spec, ///< TODOCUMENT
                                                                       superposition_content_spec  prm_content_spec  ///< The specification of what should be included in the superposition
                                                                       ) : output_file      { prm_output_file               },
                                                                           the_display_spec { std::move( prm_display_spec ) },
                                                                           content_spec     { std::move( prm_content_spec ) } {
}

