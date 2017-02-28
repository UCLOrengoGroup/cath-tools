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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_IO_SUPERPOSITION_IO_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_IO_SUPERPOSITION_IO_H

#include <boost/filesystem.hpp>

#include "chopping/region/regions_limiter.hpp"
#include "common/path_type_aliases.hpp"
#include "file/pdb/pdb_write_mode.hpp"
#include "superposition/io/sup_pdbs_script_policy.hpp"
#include "superposition/superposition.hpp"

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

		} // namespace detail

		enum class chain_relabel_policy : bool {
			RELABEL,
			LEAVE
		};

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
		                                               const chain_relabel_policy & = chain_relabel_policy::LEAVE,
		                                               const chop::regions_limiter & = chop::regions_limiter{},
		                                               const file::pdb_write_mode & = file::pdb_write_mode::ONLY_OR_LAST_PDB);

		std::ostream & write_superposed_pdbs_to_ostream(std::ostream &,
		                                                const superposition &,
		                                                const file::pdb_list,
		                                                const sup_pdbs_script_policy &,
		                                                const chain_relabel_policy & = chain_relabel_policy::LEAVE,
		                                                const chop::regions_limiter & = chop::regions_limiter{});

		void write_superposed_pdb_to_file(const superposition &,
		                                  const boost::filesystem::path &,
		                                  const file::pdb &,
		                                  const size_t &,
		                                  const chain_relabel_policy & = chain_relabel_policy::LEAVE,
		                                  const chop::regions_limiter & = chop::regions_limiter{});

		void write_superposed_pdb_to_file(const superposition &,
		                                  const boost::filesystem::path &,
		                                  const file::pdb_list &,
		                                  const sup_pdbs_script_policy &,
		                                  const chain_relabel_policy & = chain_relabel_policy::LEAVE,
		                                  const chop::regions_limiter & = chop::regions_limiter{});

		void write_superposed_pdb_from_files(const superposition &,
		                                     const boost::filesystem::path &,
		                                     const path_vec &,
		                                     const sup_pdbs_script_policy &,
		                                     const chain_relabel_policy & = chain_relabel_policy::LEAVE,
		                                     const chop::regions_limiter & = chop::regions_limiter{});

		superposition superposition_from_ptree(const boost::property_tree::ptree &);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const superposition &);

		boost::property_tree::ptree make_ptree_of(const superposition &);

		superposition superposition_from_json_string(const std::string &);

		std::string to_json_string(const superposition &,
		                           const bool & = true);
	} // namespace sup
} // namespace cath

#endif
