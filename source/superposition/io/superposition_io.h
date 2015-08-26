/// \file
/// \brief The superposition_io function headers

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

#ifndef SUPERPOSITION_IO_H_INCLUDED
#define SUPERPOSITION_IO_H_INCLUDED

#include <boost/filesystem.hpp>

#include "superposition/superposition.h"

#include <iosfwd>
#include <vector>

namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace sup {
		namespace detail {

			/// \brief Store constants to be used in superposition I/O code
			struct superposition_io_consts final {
				static const std::string ENTRIES_KEY;
				static const std::string NAME_KEY;
				static const std::string ROTATION_KEY;
				static const std::string TRANSFORMATION_KEY;
				static const std::string TRANSFORMATIONS_KEY;
				static const std::string TRANSLATION_KEY;
			};

		}

		void write_xml_sup(std::ostream &,
		                   const superposition &,
		                   const str_vec &);

		void write_xml_sup_filename(const superposition &,
		                            const boost::filesystem::path &,
		                            const str_vec &);

		std::ostream & write_superposed_pdb_to_ostream(std::ostream &,
		                                               const superposition &,
		                                               file::pdb,
		                                               const size_t &,
		                                               const bool &arg_relabel_chain = false);

		std::ostream & write_superposed_pdbs_to_ostream(std::ostream &,
		                                               const superposition &,
		                                               const file::pdb_list,
		                                               const bool &,
		                                               const bool &arg_relabel_chain = false);

		void write_superposed_pdb_to_file(const superposition &,
		                                  const boost::filesystem::path &,
		                                  const file::pdb &,
		                                  const size_t &,
		                                  const bool &arg_relabel_chain = false);

		void write_superposed_pdb_to_file(const superposition &,
		                                  const boost::filesystem::path &,
		                                  const file::pdb_list &,
		                                  const bool &,
		                                  const bool &arg_relabel_chain = false);

		void write_superposed_pdb_from_files(const superposition &,
		                                     const boost::filesystem::path &,
		                                     const path_vec &,
		                                     const bool &,
		                                     const bool &arg_relabel_chain = false);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const superposition &);

		boost::property_tree::ptree make_ptree_of(const superposition &);

		std::string to_json_string(const superposition &,
		                           const bool &arg_pretty_print = true);
	}
}

#endif
