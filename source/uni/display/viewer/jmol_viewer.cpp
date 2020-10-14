/// \file
/// \brief The jmol_viewer class definitions

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

#include "jmol_viewer.hpp"

#include "common/exception/not_implemented_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::sup;
using namespace std;

using boost::string_ref;

/// \brief TODOCUMENT
string jmol_viewer::do_default_executable() const {
	return "jmol";
}

/// \brief TODOCUMENT
void jmol_viewer::do_write_start(ostream &) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_write_start"));
}

/// \brief TODOCUMENT
void jmol_viewer::do_write_load_pdbs(ostream &,
                                     const sup::superposition &,
                                     const file::pdb_list &,
                                     const str_vec &) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_write_load_pdbs"));
}

/// \brief TODOCUMENT
void jmol_viewer::do_define_colour(ostream &,
                                   const display_colour &,
                                   const string &) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_define_colour"));
}

/// \brief Get a string for colouring the base (ie everything) in the colour that has previously been defined with the specified name
string jmol_viewer::do_get_colour_base_str(const string &
                                           ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_get_colour_base_str"));
}

/// \brief TODOCUMENT
string jmol_viewer::do_get_colour_pdb_str(const string &,
                                          const string &) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_colour_pdb"));
}

/// \brief TODOCUMENT
string jmol_viewer::do_get_colour_pdb_residues_str(const string &,
                                                   const string &,
                                                   const residue_id_vec &) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_colour_pdb_residues"));
}

/// \brief TODOCUMENT
void jmol_viewer::do_write_alignment_extras(ostream                     &/*prm_os*/,                   ///< TODOCUMENT
                                            const superposition_context &/*prm_superposition_context*/ ///< TODOCUMENT
                                            ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_write_alignment_extras"));
}

/// \brief TODOCUMENT
void jmol_viewer::do_write_end(ostream          &/* prm_os */,            ///< TODOCUMENT
                               const string_ref &/* prm_advert_message */ ///< TODOCUMENT
                               ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("jmol_viewer::do_write_end"));
}
