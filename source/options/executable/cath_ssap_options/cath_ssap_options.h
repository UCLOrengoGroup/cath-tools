/// \file
/// \brief The cath_ssap_options class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef CATH_SSAP_OPTIONS_H_INCLUDED
#define CATH_SSAP_OPTIONS_H_INCLUDED

#include "options/executable/executable_options.h"
#include "options/options_block/data_dirs_options_block.h"
#include "options/options_block/detail_help_options_block.h"
#include "options/options_block/misc_help_version_options_block.h"
#include "options/options_block/old_ssap_options_block.h"

namespace cath {
	namespace opts {

		/// \brief Handle the options associated with the cath_ssap executable
		///
		/// This uses executable_options to achieve a lot of its work
		class cath_ssap_options final : public executable_options {
		private:
			using super = executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;
			static const std::map<std::string, str_str_pair> DETAIL_HELP_SPEC;

			/// \brief TODOCUMENT
			misc_help_version_options_block the_misc_options_block;

			/// \brief TODOCUMENT
			old_ssap_options_block          the_ssap_options_block;

			/// \brief TODOCUMENT
			data_dirs_options_block         the_data_dirs_options_block;

			/// \brief TODOCUMENT
			detail_help_options_block       the_detail_help_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual boost::program_options::positional_options_description get_positional_options() override final;
			virtual std::string do_update_error_or_help_string(const boost::program_options::options_description &) const override final;

			static std::string get_version_description_string();
			static std::string get_usage_string();

			void check_ok_to_use() const;

		public:
			cath_ssap_options();
			virtual ~cath_ssap_options() noexcept = default;

			const old_ssap_options_block  get_old_ssap_options()  const;
			const data_dirs_options_block get_data_dirs_options() const;

			static const std::string PROGRAM_NAME;
		};

		std::string get_ssap_matches_format_help_string();

		std::string get_ssap_alignment_format_help_string();
		
	}
}

#endif
