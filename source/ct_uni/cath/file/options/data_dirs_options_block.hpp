/// \file
/// \brief The data_dirs_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_OPTIONS_BLOCK_HPP

#include "cath/file/options/data_dirs_spec.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath {
	namespace opts {

		/// \brief Manage the program options to populate a data_dirs_spec
		///
		/// Note that there are two different uses of the word 'path' around this code:
		///  * the Boost filesystem type that is used to store the location of a file or directory
		///  * a prioritised list of directories to search for some sort of file
		class data_dirs_options_block final : public cath::opts::options_block {
		private:
			static const std::string DATA_OPTION_PATH_VARNAME;
			static const std::string DATA_OPTION_PREFIX_VARNAME;
			static const std::string DATA_OPTION_SUFFIX_VARNAME;

			static const data_option_str_map DATA_OPTION_SUFFIXES;
			static const data_option_str_map DATA_OPTION_VARNAME;
			static const data_option_str_map DATA_OPTION_DESCRIPTION_START;
			static const data_option_str_map DATA_OPTION_DESCRIPTION_END;

			[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
			[[nodiscard]] std::string                    do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			void do_add_hidden_options_to_description(boost::program_options::options_description &,
			                                          const size_t &) final;
			[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
			[[nodiscard]] str_vec do_get_all_options_names() const final;

			/// \brief The data_dirs_spec into which to parse options
			data_dirs_spec the_data_dirs_spec;

		public:
			data_dirs_options_block() = default;

			[[nodiscard]] const data_dirs_spec &get_data_dirs_spec() const;

			static const std::string PO_CATH_ROOT_DIR;
		};
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_OPTIONS_BLOCK_HPP
