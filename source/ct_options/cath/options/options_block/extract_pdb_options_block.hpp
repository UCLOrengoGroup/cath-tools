/// \file
/// \brief The extract_pdb_options_block class header

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

#ifndef CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_EXTRACT_PDB_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_EXTRACT_PDB_OPTIONS_BLOCK_HPP

#include <filesystem>
#include <string_view>

#include "cath/chopping/domain/domain.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath::opts {

	/// \brief Handle the options specific to extract-pdb
	class extract_pdb_options_block final : public options_block {
	private:
		using super = options_block;

		/// \brief The PDB file to extract
		::std::filesystem::path input_pdb_file;

		/// \brief TODOCUMENT
		path_opt output_pdb_file;

		/// \brief TODOCUMENT
		chop::domain_opt regions;

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		void do_add_hidden_options_to_description(boost::program_options::options_description &,
		                                          const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const ::std::filesystem::path &get_input_pdb_file() const;
		[[nodiscard]] const path_opt &               get_output_pdb_file() const;
		[[nodiscard]] const chop::domain_opt &       get_regions() const;

		static constexpr ::std::string_view PO_INPUT_PDB_FILE{ "input-pdb-file" };
		static constexpr ::std::string_view PO_OUTPUT_PDB_FILE{ "output-pdb-file" };
		static constexpr ::std::string_view PO_REGIONS{ "regions" };
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_EXTRACT_PDB_OPTIONS_BLOCK_HPP
