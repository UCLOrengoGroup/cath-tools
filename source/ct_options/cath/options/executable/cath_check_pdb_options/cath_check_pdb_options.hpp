/// \file
/// \brief The cath_check_pdb_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS_CATH_CHECK_PDB_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS_CATH_CHECK_PDB_OPTIONS_HPP

#include <filesystem>
#include <string_view>
#include <vector>

#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/check_pdb_options_block.hpp"

namespace cath::opts {

	/// \brief TODOCUMENT
	class cath_check_pdb_options final : public executable_options {
	private:
		using super = executable_options;

		/// \brief TODOCUMENT
		check_pdb_options_block the_check_pdb_options_block;

		[[nodiscard]] std::string_view                         do_get_program_name() const final;
		boost::program_options::positional_options_description get_positional_options() final;
		[[nodiscard]] str_opt                                  do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

	  public:
		cath_check_pdb_options();

		[[nodiscard]] ::std::filesystem::path get_pdb_file() const;
		[[nodiscard]] bool                    get_permit_no_atoms() const;

		/// \brief The name of the program that uses this executable_options
		static constexpr ::std::string_view PROGRAM_NAME{ "check-pdb" };
	};

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS_CATH_CHECK_PDB_OPTIONS_HPP
