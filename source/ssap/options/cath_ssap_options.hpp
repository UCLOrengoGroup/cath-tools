/// \file
/// \brief The cath_ssap_options class header

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

#ifndef _CATH_TOOLS_SOURCE_SSAP_OPTIONS_CATH_SSAP_OPTIONS_H
#define _CATH_TOOLS_SOURCE_SSAP_OPTIONS_CATH_SSAP_OPTIONS_H

#include "file/options/data_dirs_options_block.hpp"
#include "options/executable/executable_options.hpp"
#include "options/options_block/detail_help_options_block.hpp"
#include "ssap/options/old_ssap_options_block.hpp"

namespace cath {
	namespace opts {

		/// \brief Handle the options associated with the cath_ssap executable
		///
		/// This uses executable_options to achieve a lot of its work
		class cath_ssap_options final : public executable_options {
		private:
			using super = executable_options;

			// static const std::string STANDARD_USAGE_ERROR_STRING;
			static std::map<std::string, str_str_pair> detail_help_spec();

			static const std::string PO_CITATION_HELP;

			/// \brief TODOCUMENT
			old_ssap_options_block          the_ssap_options_block;

			/// \brief TODOCUMENT
			data_dirs_options_block         the_data_dirs_options_block;

			/// \brief TODOCUMENT
			detail_help_options_block       the_detail_help_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual boost::program_options::positional_options_description get_positional_options() override final;
			virtual str_opt do_get_error_or_help_string() const override final;

			virtual std::string do_get_help_prefix_string() const override final;
			virtual std::string do_get_help_suffix_string() const override final;
			virtual std::string do_get_overview_string() const override final;

			static std::string get_version_description_string();

			void check_ok_to_use() const;

		public:
			cath_ssap_options();

			const old_ssap_options_block & get_old_ssap_options() const;
			const data_dirs_spec &         get_data_dirs_spec() const;

			static const std::string PROGRAM_NAME;
		};

		std::string get_ssap_matches_format_help_string();

		std::string get_ssap_alignment_format_help_string();

		std::string get_ssap_citation_help_string();
		
	} // namespace opts
} // namespace cath

#endif
