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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_IO_SUPERPOSITION_IO_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_IO_SUPERPOSITION_IO_HPP

#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string_view>
#include <vector>

#include "cath/chopping/region/region.hpp"
#include "cath/chopping/region/regions_limiter.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/property_tree/read_from_ptree.hpp"
#include "cath/file/pdb/pdb_write_mode.hpp"
#include "cath/superposition/io/sup_pdbs_script_policy.hpp"
#include "cath/superposition/superposition.hpp"

// clang-format off
namespace cath::file { class pdb; }
namespace cath::file { class pdb_list; }
// clang-format on

namespace cath::sup {

	namespace detail {

		/// \brief Store constants to be used in superposition I/O code
		struct superposition_io_consts final {

			static constexpr ::std::string_view ENTRIES_KEY         = "entries";
			static constexpr ::std::string_view NAME_KEY            = "name";
			static constexpr ::std::string_view ROTATION_KEY        = "rotation";
			static constexpr ::std::string_view TRANSFORMATION_KEY  = "transformation";
			static constexpr ::std::string_view TRANSFORMATIONS_KEY = "transformations";
			static constexpr ::std::string_view TRANSLATION_KEY     = "translation";
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
	                            const ::std::filesystem::path &,
	                            const str_vec &);

	std::ostream & write_superposed_pdb_to_ostream(std::ostream &,
	                                               const superposition &,
	                                               file::pdb,
	                                               const size_t &,
	                                               const chain_relabel_policy & = chain_relabel_policy::LEAVE,
	                                               const chop::region_vec_opt & = ::std::nullopt,
	                                               const file::pdb_write_mode & = file::pdb_write_mode::ONLY_OR_LAST_PDB);

	std::ostream & write_superposed_pdbs_to_ostream(std::ostream &,
	                                                const superposition &,
	                                                const file::pdb_list &,
	                                                const sup_pdbs_script_policy &,
	                                                const chain_relabel_policy & = chain_relabel_policy::LEAVE,
	                                                const chop::region_vec_opt & = ::std::nullopt);

	void write_superposed_pdb_to_file(const superposition &,
	                                  const ::std::filesystem::path &,
	                                  const file::pdb &,
	                                  const size_t &,
	                                  const chain_relabel_policy & = chain_relabel_policy::LEAVE,
	                                  const chop::region_vec_opt & = ::std::nullopt);

	void write_superposed_pdb_to_file(const superposition &,
	                                  const ::std::filesystem::path &,
	                                  const file::pdb_list &,
	                                  const sup_pdbs_script_policy &,
	                                  const chain_relabel_policy & = chain_relabel_policy::LEAVE,
	                                  const chop::region_vec_opt & = ::std::nullopt);

	void write_superposed_pdb_from_files(const superposition &,
	                                     const ::std::filesystem::path &,
	                                     const path_vec &,
	                                     const sup_pdbs_script_policy &,
	                                     const chain_relabel_policy & = chain_relabel_policy::LEAVE,
	                                     const chop::region_vec_opt & = ::std::nullopt);

	superposition superposition_from_ptree(const boost::property_tree::ptree &);

	void save_to_ptree(boost::property_tree::ptree &,
	                   const superposition &);

} // namespace cath::sup

namespace cath::common {
	
	/// \brief Specialisation of cath::common::read_from_ptree for superposition
	template <>
	inline sup::superposition read_from_ptree<sup::superposition>(const boost::property_tree::ptree &prm_ptree ///< The ptree from which to read the superposition
	                                                              ) {
		return sup::superposition_from_ptree( prm_ptree );
	}
	
} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_IO_SUPERPOSITION_IO_HPP
