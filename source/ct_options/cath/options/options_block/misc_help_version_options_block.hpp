/// \file
/// \brief The misc_help_version_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_MISC_HELP_VERSION_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_MISC_HELP_VERSION_OPTIONS_BLOCK_HPP

#include "cath/options/options_block/options_block.hpp"

namespace cath {
	namespace opts {

		/// \brief A standard options block for the miscellaneous help/version information to be used in most programs
		///
		/// Note: I previously made this class generate the strings as part of its do_invalid_string()
		///       behaviour (which required altering all block_options classes to take an options_description argument
		///       to do_invalid_string() and storing all necessary extra information (help_message_prefix, help_message_suffix,
		///       program_name and version_program_description).
		///       I wasn't fond of the result, in large part because the help/version messages ended up getting appended with:
		///       > Try 'cath-superpose --help' for usage information.
		///       or similar.
		class misc_help_version_options_block final : public options_block {
		private:
			using super = options_block;

			/// \brief Whether hidden_help has been requested
			bool hidden_help = false;

			/// \brief Whether help has been requested
			bool help        = false;

			/// \brief Whether version information has been requested
			bool version     = false;

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
			void do_add_hidden_options_to_description(boost::program_options::options_description &,
			                                          const size_t &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;
			str_vec do_get_all_options_names() const final;

		public:
			const bool & get_hidden_help() const;
			const bool & get_help() const;
			const bool & get_version() const;

			static std::string get_help_string(const boost::program_options::options_description &,
			                                   const std::string &,
			                                   const std::string &);
			static std::string get_version_string(const std::string &,
			                                      const std::string &);

			static const std::string PO_HIDDEN_HELP;
			static const std::string PO_HELP;
			static const std::string PO_VERSION;

			static constexpr char PO_CHAR_HELP    = 'h';
			static constexpr char PO_CHAR_VERSION = 'v';
		};
	} // namespace opts
} // namespace cath

#endif
