/// \file
/// \brief The domain_definition_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DOMAIN_DEFINITION_LIST_DOMAIN_DEFINITION_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DOMAIN_DEFINITION_LIST_DOMAIN_DEFINITION_LIST_HPP

#include <filesystem>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"

// clang-format off
namespace cath::opts { class data_dirs_spec; }
// clang-format on

namespace cath::file {

	/// \brief TODOCUMENT
	class domain_definition_list final {
	private:
		/// \brief TODOCUMENT
		chop::domain_definition_vec domain_definitions;

	public:
		using const_iterator = chop::domain_definition_vec::const_iterator;

		explicit domain_definition_list(chop::domain_definition_vec);

		[[nodiscard]] size_t size() const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	domain_definition_list parse_domain_definition_file(const ::std::filesystem::path &);
	domain_definition_list parse_domain_definition_file(std::istream &);

	pdb_list_name_set_list_pair read_domains_from_pdbs(const domain_definition_list &,
	                                                   const opts::data_dirs_spec &);

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DOMAIN_DEFINITION_LIST_DOMAIN_DEFINITION_LIST_HPP
